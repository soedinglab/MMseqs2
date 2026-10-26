/**
 * @file brotli.c
 * @brief BROTLI (Parquet codec 4) via libbrotli.
 *
 * Parquet stores a raw Brotli stream per page (no framing). The library is
 * optional: when carquet is built without it (CARQUET_HAVE_BROTLI undefined,
 * i.e. -DCARQUET_ENABLE_BROTLI=OFF or no brotli found), both directions report
 * CARQUET_ERROR_UNSUPPORTED_CODEC and a user codec registered with
 * carquet_register_codec() can still fill the slot.
 *
 * Quality: Parquet/Arrow map the writer's compression level straight onto
 * Brotli's quality (0..11); Arrow's default is 8, which is also what a level
 * of 0 ("default") selects here. The window stays at BROTLI_DEFAULT_WINDOW
 * (22, 4 MB), matching Arrow, so files are byte-compatible with pyarrow's.
 *
 * libbrotli's one-shot encoder/decoder keep no state between calls, so no
 * per-thread context is needed (unlike zstd / gzip).
 */

#include <carquet/error.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#if defined(CARQUET_HAVE_BROTLI)
#include <brotli/decode.h>
#include <brotli/encode.h>
#endif

int carquet_brotli_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
#if defined(CARQUET_HAVE_BROTLI)
    size_t out = dst_capacity;
    BrotliDecoderResult r = BrotliDecoderDecompress(src_size, src, &out, dst);
    if (r != BROTLI_DECODER_RESULT_SUCCESS) {
        /* NEEDS_MORE_OUTPUT (page header lied about the uncompressed size),
         * NEEDS_MORE_INPUT (truncated) and ERROR all mean the page is bad. */
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }
    *dst_size = out;
    return CARQUET_OK;
#else
    (void)src_size; (void)dst_capacity;
    return CARQUET_ERROR_UNSUPPORTED_CODEC;
#endif
}

int carquet_brotli_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size,
    int level) {

    if (!src || !dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
#if defined(CARQUET_HAVE_BROTLI)
    int quality = level;
    if (quality <= 0) quality = 8;                      /* Arrow's default */
    if (quality > BROTLI_MAX_QUALITY) quality = BROTLI_MAX_QUALITY;
    size_t out = dst_capacity;
    if (!BrotliEncoderCompress(quality, BROTLI_DEFAULT_WINDOW, BROTLI_MODE_GENERIC,
                               src_size, src, &out, dst)) {
        return CARQUET_ERROR_COMPRESSION;
    }
    *dst_size = out;
    return CARQUET_OK;
#else
    (void)src_size; (void)dst_capacity; (void)level;
    return CARQUET_ERROR_UNSUPPORTED_CODEC;
#endif
}

size_t carquet_brotli_compress_bound(size_t src_size) {
#if defined(CARQUET_HAVE_BROTLI)
    /* BrotliEncoderMaxCompressedSize returns 0 when the input is too large to
     * bound; an uncompressed meta-block costs at most src_size + a few bytes
     * per 16 MiB block, so fall back to a generous linear bound. */
    size_t b = BrotliEncoderMaxCompressedSize(src_size);
    if (b == 0) b = src_size + src_size / 1024 + 1024;
    return b;
#else
    (void)src_size;
    return 0;
#endif
}

bool carquet_brotli_available(void) {
#if defined(CARQUET_HAVE_BROTLI)
    return true;
#else
    return false;
#endif
}
