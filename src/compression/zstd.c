/**
 * @file zstd.c
 * @brief ZSTD compression/decompression using libzstd
 *
 * Uses streaming context for better performance on repeated decompressions.
 */

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <zstd.h>

/* ============================================================================
 * Thread-local ZSTD context management
 *
 * ZSTD contexts are expensive to create (~650KB each) so we cache them per
 * thread.  They MUST be per-thread: the batch reader decompresses pages on a
 * worker-thread pool, and a ZSTD_DCtx entered concurrently corrupts and
 * crashes.  The challenge is cleanup — TLS contexts have no destructor, so
 * contexts allocated by worker threads leak when the thread pool is torn down.
 *
 * Strategy:
 *   POSIX (any) -> pthread_key_create with destructors.  Works for both
 *                  OpenMP threads and worker pool pthreads.  The pthread
 *                  runtime calls the destructor when each thread exits.
 *   Windows     -> Win32 FLS API (FlsAlloc with destructor/FlsGetValue/FlsSetValue), which is
 *                  reliably per-thread for every thread including raw
 *                  CreateThread worker threads (native __declspec(thread) is
 *                  not), with explicit carquet_zstd_cleanup per thread.
 *
 * carquet_cleanup() (public API) calls carquet_zstd_cleanup() for the
 * calling thread.  On POSIX the worker-thread contexts are freed
 * automatically; on Windows callers must arrange per-thread cleanup.
 * ============================================================================ */

#if !defined(_WIN32)
/* ---- POSIX: pthread_key with destructors (works for OMP + worker pool) ---- */
#include <pthread.h>

static pthread_key_t tls_dctx_key;
static pthread_key_t tls_cctx_key;
static pthread_once_t tls_keys_once = PTHREAD_ONCE_INIT;

static void destroy_dctx(void* ctx) {
    if (ctx) ZSTD_freeDCtx((ZSTD_DCtx*)ctx);
}

static void destroy_cctx(void* ctx) {
    if (ctx) ZSTD_freeCCtx((ZSTD_CCtx*)ctx);
}

static void init_tls_keys(void) {
    pthread_key_create(&tls_dctx_key, destroy_dctx);
    pthread_key_create(&tls_cctx_key, destroy_cctx);
}

static ZSTD_DCtx* get_dctx(void) {
    pthread_once(&tls_keys_once, init_tls_keys);
    ZSTD_DCtx* dctx = (ZSTD_DCtx*)pthread_getspecific(tls_dctx_key);
    if (!dctx) {
        dctx = ZSTD_createDCtx();
        if (dctx) pthread_setspecific(tls_dctx_key, dctx);
    }
    return dctx;
}

static ZSTD_CCtx* get_cctx(void) {
    pthread_once(&tls_keys_once, init_tls_keys);
    ZSTD_CCtx* cctx = (ZSTD_CCtx*)pthread_getspecific(tls_cctx_key);
    if (!cctx) {
        cctx = ZSTD_createCCtx();
        if (cctx) pthread_setspecific(tls_cctx_key, cctx);
    }
    return cctx;
}

void carquet_zstd_cleanup(void) {
    pthread_once(&tls_keys_once, init_tls_keys);
    ZSTD_DCtx* dctx = (ZSTD_DCtx*)pthread_getspecific(tls_dctx_key);
    if (dctx) {
        ZSTD_freeDCtx(dctx);
        pthread_setspecific(tls_dctx_key, NULL);
    }
    ZSTD_CCtx* cctx = (ZSTD_CCtx*)pthread_getspecific(tls_cctx_key);
    if (cctx) {
        ZSTD_freeCCtx(cctx);
        pthread_setspecific(tls_cctx_key, NULL);
    }
}

#else
/* ---- Windows: Win32 TLS API (per-thread, works for every thread) ----
 *
 * Per-thread, NOT global: the batch reader runs its own worker pool
 * (carquet_worker_pool, plain Win32 threads) to decompress pages in
 * parallel whether or not OpenMP is present. A single shared ZSTD_DCtx
 * would then be entered concurrently by several threads — ZSTD_DCtx is
 * not thread-safe, so its internals corrupt and decode reads a wild
 * pointer (crash).
 *
 * We use TlsAlloc/TlsGetValue/TlsSetValue rather than __declspec(thread):
 * native TLS is NOT reliably allocated per-thread for threads created
 * with the raw CreateThread API under every loader/runtime (observed all
 * worker threads sharing one slot), whereas the explicit TLS API is
 * guaranteed per-thread for every thread. The FLS variant of that API adds
 * a per-thread destructor, so worker-thread contexts are freed when the
 * thread exits (v0.7.2; TlsAlloc leaked them at pool teardown). */
#include <windows.h>

/* Fiber-local storage (Vista+): same per-thread guarantee as TlsAlloc for
 * every thread, including raw CreateThread pool workers, plus a destructor
 * that the system runs when the thread exits, so worker contexts no longer
 * leak at pool teardown. FlsSetValue(idx, NULL) in the explicit cleanup
 * below does not invoke the callback, so there is no double free. */
static DWORD tls_dctx_index = FLS_OUT_OF_INDEXES;
static DWORD tls_cctx_index = FLS_OUT_OF_INDEXES;

static void WINAPI fls_free_dctx(PVOID p) { if (p) ZSTD_freeDCtx((ZSTD_DCtx*)p); }
static void WINAPI fls_free_cctx(PVOID p) { if (p) ZSTD_freeCCtx((ZSTD_CCtx*)p); }
static INIT_ONCE tls_index_once = INIT_ONCE_STATIC_INIT;

static BOOL CALLBACK init_tls_indices(PINIT_ONCE once, PVOID param, PVOID* ctx) {
    (void)once; (void)param; (void)ctx;
    tls_dctx_index = FlsAlloc(fls_free_dctx);
    tls_cctx_index = FlsAlloc(fls_free_cctx);
    return TRUE;
}

static void ensure_tls_indices(void) {
    InitOnceExecuteOnce(&tls_index_once, init_tls_indices, NULL, NULL);
}

static ZSTD_DCtx* get_dctx(void) {
    ensure_tls_indices();
    if (tls_dctx_index == FLS_OUT_OF_INDEXES) return NULL;
    ZSTD_DCtx* dctx = (ZSTD_DCtx*)FlsGetValue(tls_dctx_index);
    if (!dctx) {
        dctx = ZSTD_createDCtx();
        if (dctx) FlsSetValue(tls_dctx_index, dctx);
    }
    return dctx;
}

static ZSTD_CCtx* get_cctx(void) {
    ensure_tls_indices();
    if (tls_cctx_index == FLS_OUT_OF_INDEXES) return NULL;
    ZSTD_CCtx* cctx = (ZSTD_CCtx*)FlsGetValue(tls_cctx_index);
    if (!cctx) {
        cctx = ZSTD_createCCtx();
        if (cctx) FlsSetValue(tls_cctx_index, cctx);
    }
    return cctx;
}

void carquet_zstd_cleanup(void) {
    ensure_tls_indices();
    if (tls_dctx_index != FLS_OUT_OF_INDEXES) {
        ZSTD_DCtx* dctx = (ZSTD_DCtx*)FlsGetValue(tls_dctx_index);
        if (dctx) { ZSTD_freeDCtx(dctx); FlsSetValue(tls_dctx_index, NULL); }
    }
    if (tls_cctx_index != FLS_OUT_OF_INDEXES) {
        ZSTD_CCtx* cctx = (ZSTD_CCtx*)FlsGetValue(tls_cctx_index);
        if (cctx) { ZSTD_freeCCtx(cctx); FlsSetValue(tls_cctx_index, NULL); }
    }
}
#endif

int carquet_zstd_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Use streaming context for better buffer reuse */
    ZSTD_DCtx* dctx = get_dctx();
    if (!dctx) {
        /* Fallback to simple API */
        size_t result = ZSTD_decompress(dst, dst_capacity, src, src_size);
        if (ZSTD_isError(result)) {
            return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
        }
        *dst_size = result;
        return CARQUET_OK;
    }

    size_t result = ZSTD_decompressDCtx(dctx, dst, dst_capacity, src, src_size);
    if (ZSTD_isError(result)) {
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }

    *dst_size = result;
    return CARQUET_OK;
}

int carquet_zstd_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size,
    int level) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (level < 1) level = 1;
    if (level > ZSTD_maxCLevel()) level = ZSTD_maxCLevel();

    /* Use cached context for repeated compressions (e.g., per-page).
     * Only enable multi-threading for large inputs (>4MB) where the
     * parallelism overhead is worthwhile. */
    ZSTD_CCtx* cctx = get_cctx();
    if (cctx) {
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, level);

        /* Only use multi-threading for large inputs where parallelism
         * outweighs coordination overhead. For typical 1MB pages, single-
         * threaded with a cached context is faster. */
        if (src_size > 4 * 1024 * 1024) {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, 4);
        } else {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, 0);
        }

        size_t result = ZSTD_compress2(cctx, dst, dst_capacity, src, src_size);
        if (!ZSTD_isError(result)) {
            *dst_size = result;
            return CARQUET_OK;
        }
    }

    /* Fallback to simple API */
    size_t result = ZSTD_compress(dst, dst_capacity, src, src_size, level);
    if (ZSTD_isError(result)) {
        return CARQUET_ERROR_COMPRESSION;
    }

    *dst_size = result;
    return CARQUET_OK;
}

size_t carquet_zstd_compress_bound(size_t src_size) {
    return ZSTD_compressBound(src_size);
}

void carquet_zstd_init_tables(void) {
    /* No-op - libzstd handles initialization internally */
}
