/**
 * @file custom.c
 * @brief Pluggable codec registration table.
 *
 * Stores at most one user-supplied implementation per `carquet_compression_t`
 * value. The reader and writer check this table before falling through to the
 * built-in implementation, so registering a codec overrides built-ins as well
 * as filling slots that have no built-in (LZO, BROTLI).
 *
 * Registration mutates a process-wide table and is not safe to interleave
 * with concurrent compress/decompress calls; the public API documents this.
 */
#include "custom.h"
#include <string.h>

/* Number of codec slots in carquet_compression_t. Keep in sync with the enum
 * in include/carquet/types.h; the largest value today is LZ4_RAW = 7. */
#define CARQUET_CODEC_SLOTS 8

static carquet_custom_codec_t g_custom_codecs[CARQUET_CODEC_SLOTS];
static bool g_custom_codec_set[CARQUET_CODEC_SLOTS];

static bool slot_in_range(carquet_compression_t codec) {
    return (int)codec >= 0 && (int)codec < CARQUET_CODEC_SLOTS;
}

carquet_status_t carquet_register_codec(
    carquet_compression_t codec,
    const carquet_custom_codec_t* impl) {
    if (!slot_in_range(codec)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    /* UNCOMPRESSED is a no-op path with a no-copy fast lane in the writer;
     * overriding it would only confuse things. */
    if (codec == CARQUET_COMPRESSION_UNCOMPRESSED) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (impl == NULL) {
        g_custom_codec_set[codec] = false;
        memset(&g_custom_codecs[codec], 0, sizeof(g_custom_codecs[codec]));
        return CARQUET_OK;
    }
    /* compress, decompress, and compress_bound are all required for a usable
     * codec — partial registration would surface as a NULL deref later. */
    if (!impl->compress || !impl->decompress || !impl->compress_bound) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    g_custom_codecs[codec] = *impl;
    g_custom_codec_set[codec] = true;
    return CARQUET_OK;
}

bool carquet_custom_codec_lookup(carquet_compression_t codec,
                                 carquet_custom_codec_t* out) {
    if (!slot_in_range(codec) || !g_custom_codec_set[codec]) {
        return false;
    }
    if (out) *out = g_custom_codecs[codec];
    return true;
}
