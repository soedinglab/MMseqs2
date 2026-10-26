/**
 * @file gzip.c
 * @brief GZIP compression/decompression using zlib
 *
 * Parquet's GZIP codec is the RFC 1952 gzip format (zlib windowBits 15+16),
 * not raw DEFLATE.
 */

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* ============================================================================
 * Per-thread cached zlib streams
 * ============================================================================
 * inflateInit2 / deflateInit2 allocate and initialise the window and hash
 * tables (about 32 KB for inflate, ~256 KB for deflate at the default level)
 * on every call, which is a sizeable fixed cost per page. Streams are cached
 * per thread and reset between pages, like the ZSTD contexts in zstd.c:
 * POSIX via pthread_key with a destructor, Windows via the TLS API; any
 * failure to obtain a cached stream falls back to a per-call stream. */

typedef struct {
    z_stream inf;
    bool inf_ok;
    z_stream def;
    bool def_ok;
    int def_level;
} gzip_tls_t;

static void gzip_tls_free(void* p) {
    gzip_tls_t* t = (gzip_tls_t*)p;
    if (!t) return;
    if (t->inf_ok) inflateEnd(&t->inf);
    if (t->def_ok) deflateEnd(&t->def);
    free(t);
}

#if !defined(_WIN32)
#include <pthread.h>
static pthread_key_t gzip_tls_key;
static pthread_once_t gzip_tls_once = PTHREAD_ONCE_INIT;
static int gzip_tls_key_ok = 0;
static void gzip_tls_init_key(void) {
    gzip_tls_key_ok = (pthread_key_create(&gzip_tls_key, gzip_tls_free) == 0);
}
static gzip_tls_t* gzip_tls_get(void) {
    pthread_once(&gzip_tls_once, gzip_tls_init_key);
    if (!gzip_tls_key_ok) return NULL;
    gzip_tls_t* t = (gzip_tls_t*)pthread_getspecific(gzip_tls_key);
    if (!t) {
        t = (gzip_tls_t*)calloc(1, sizeof(*t));
        if (t && pthread_setspecific(gzip_tls_key, t) != 0) { free(t); t = NULL; }
    }
    return t;
}
void carquet_gzip_cleanup(void) {
    pthread_once(&gzip_tls_once, gzip_tls_init_key);
    if (!gzip_tls_key_ok) return;
    gzip_tls_t* t = (gzip_tls_t*)pthread_getspecific(gzip_tls_key);
    if (t) {
        gzip_tls_free(t);
        pthread_setspecific(gzip_tls_key, NULL);
    }
}
#else
#include <windows.h>
static DWORD gzip_tls_index = FLS_OUT_OF_INDEXES;
/* FLS destructor: frees the cached streams when a thread exits (see zstd.c). */
static void WINAPI gzip_fls_free(PVOID p) { gzip_tls_free(p); }
static INIT_ONCE gzip_tls_once = INIT_ONCE_STATIC_INIT;
static BOOL CALLBACK gzip_tls_init_index(PINIT_ONCE once, PVOID param, PVOID* ctx) {
    (void)once; (void)param; (void)ctx;
    gzip_tls_index = FlsAlloc(gzip_fls_free);
    return TRUE;
}
static gzip_tls_t* gzip_tls_get(void) {
    InitOnceExecuteOnce(&gzip_tls_once, gzip_tls_init_index, NULL, NULL);
    if (gzip_tls_index == FLS_OUT_OF_INDEXES) return NULL;
    gzip_tls_t* t = (gzip_tls_t*)FlsGetValue(gzip_tls_index);
    if (!t) {
        t = (gzip_tls_t*)calloc(1, sizeof(*t));
        if (t && !FlsSetValue(gzip_tls_index, t)) { free(t); t = NULL; }
    }
    return t;
}
/* Worker-thread streams leak at pool teardown (no TLS destructor); callers
 * arrange per-thread carquet_gzip_cleanup where it matters, as for ZSTD. */
void carquet_gzip_cleanup(void) {
    InitOnceExecuteOnce(&gzip_tls_once, gzip_tls_init_index, NULL, NULL);
    if (gzip_tls_index == FLS_OUT_OF_INDEXES) return;
    gzip_tls_t* t = (gzip_tls_t*)FlsGetValue(gzip_tls_index);
    if (t) {
        gzip_tls_free(t);
        FlsSetValue(gzip_tls_index, NULL);
    }
}
#endif

int carquet_gzip_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size > (size_t)UINT_MAX || dst_capacity > (size_t)UINT_MAX) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    z_stream local;
    z_stream* strm = NULL;
    bool cached = false;
    gzip_tls_t* t = gzip_tls_get();
    if (t) {
        if (t->inf_ok) {
            if (inflateReset(&t->inf) == Z_OK) { strm = &t->inf; cached = true; }
        } else {
            memset(&t->inf, 0, sizeof(t->inf));
            /* 15 + 16 = gzip format (RFC 1952) */
            if (inflateInit2(&t->inf, 15 + 16) == Z_OK) {
                t->inf_ok = true; strm = &t->inf; cached = true;
            }
        }
    }
    if (!strm) {
        memset(&local, 0, sizeof(local));
        if (inflateInit2(&local, 15 + 16) != Z_OK) {
            return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
        }
        strm = &local;
    }
    strm->next_in = (Bytef*)src;
    strm->avail_in = (uInt)src_size;
    strm->next_out = (Bytef*)dst;
    strm->avail_out = (uInt)dst_capacity;

    int ret = inflate(strm, Z_FINISH);
    size_t output_size = strm->total_out;
    if (!cached) inflateEnd(strm);

    if (ret != Z_STREAM_END) {
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }

    *dst_size = output_size;
    return CARQUET_OK;
}

int carquet_gzip_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size,
    int level) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size > (size_t)UINT_MAX || dst_capacity > (size_t)UINT_MAX) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (level < 1) level = 1;
    if (level > 9) level = 9;

    z_stream local;
    z_stream* strm = NULL;
    bool cached = false;
    gzip_tls_t* t = gzip_tls_get();
    if (t) {
        if (t->def_ok && t->def_level != level) {
            /* Level change: deflateParams would be possible but the memory
             * level / tables differ anyway; re-create. */
            deflateEnd(&t->def);
            t->def_ok = false;
        }
        if (t->def_ok) {
            if (deflateReset(&t->def) == Z_OK) { strm = &t->def; cached = true; }
        } else {
            memset(&t->def, 0, sizeof(t->def));
            /* 15 + 16 = gzip format (RFC 1952) */
            if (deflateInit2(&t->def, level, Z_DEFLATED, 15 + 16, 8,
                             Z_DEFAULT_STRATEGY) == Z_OK) {
                t->def_ok = true; t->def_level = level; strm = &t->def; cached = true;
            }
        }
    }
    if (!strm) {
        memset(&local, 0, sizeof(local));
        if (deflateInit2(&local, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            return CARQUET_ERROR_COMPRESSION;
        }
        strm = &local;
    }
    strm->next_in = (Bytef*)src;
    strm->avail_in = (uInt)src_size;
    strm->next_out = (Bytef*)dst;
    strm->avail_out = (uInt)dst_capacity;

    int ret = deflate(strm, Z_FINISH);
    size_t output_size = strm->total_out;
    if (!cached) deflateEnd(strm);

    if (ret != Z_STREAM_END) {
        return CARQUET_ERROR_COMPRESSION;
    }

    *dst_size = output_size;
    return CARQUET_OK;
}

size_t carquet_gzip_compress_bound(size_t src_size) {
    /* compressBound is for zlib format; gzip adds ~18 bytes header/trailer */
    return compressBound((uLong)src_size) + 18;
}

void carquet_gzip_init_tables(void) {
    /* No-op - zlib handles initialization internally */
}
