/**
 * @file lz4.c
 * @brief LZ4 compression/decompression wrapper using the official lz4 library
 *
 * Implements LZ4 block format (LZ4_RAW) as used by Apache Parquet.
 */

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <lz4.h>

/* ============================================================================
 * LZ4 Decompression
 * ============================================================================
 */

carquet_status_t carquet_lz4_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size == 0) {
        *dst_size = 0;
        return CARQUET_OK;
    }

    if (!src) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size > (size_t)INT_MAX || dst_capacity > (size_t)INT_MAX) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int result = LZ4_decompress_safe(
        (const char*)src, (char*)dst,
        (int)src_size, (int)dst_capacity);

    if (result < 0) {
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }

    *dst_size = (size_t)result;
    return CARQUET_OK;
}

/* ============================================================================
 * LZ4 Compression
 * ============================================================================
 */

carquet_status_t carquet_lz4_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!dst || !dst_size) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size == 0) {
        *dst_size = 0;
        return CARQUET_OK;
    }

    if (!src) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (src_size > (size_t)INT_MAX || dst_capacity > (size_t)INT_MAX) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int result = LZ4_compress_default(
        (const char*)src, (char*)dst,
        (int)src_size, (int)dst_capacity);

    if (result <= 0) {
        return CARQUET_ERROR_COMPRESSION;
    }

    *dst_size = (size_t)result;
    return CARQUET_OK;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

size_t carquet_lz4_compress_bound(size_t src_size) {
    if (src_size > (size_t)INT_MAX) return 0;
    return (size_t)LZ4_compressBound((int)src_size);
}

/* ============================================================================
 * Hadoop-framed LZ4 (Parquet codec 5, the deprecated "LZ4")
 * ============================================================================
 *
 * The frame is a sequence of outer blocks, each:
 *   uint32 big-endian  total decompressed length of the outer block
 *   one or more inner blocks:
 *     uint32 big-endian  compressed length
 *     <compressed length> bytes of a raw LZ4 block
 * concatenated until the outer block's decompressed length is reached.
 * We emit the minimal conformant shape (one outer block, one inner block);
 * the decoder handles the fully general multi-block layout that legacy
 * Hadoop/Spark writers produce.
 */

static void put_be32(uint8_t* p, uint32_t v) {
    p[0] = (uint8_t)(v >> 24);
    p[1] = (uint8_t)(v >> 16);
    p[2] = (uint8_t)(v >> 8);
    p[3] = (uint8_t)v;
}

static uint32_t get_be32(const uint8_t* p) {
    return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) |
           ((uint32_t)p[2] << 8) | (uint32_t)p[3];
}

size_t carquet_lz4_hadoop_compress_bound(size_t src_size) {
    size_t inner = carquet_lz4_compress_bound(src_size);
    if (inner == 0 && src_size != 0) return 0;
    return 8 + inner;  /* outer length + inner length prefixes */
}

carquet_status_t carquet_lz4_hadoop_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!dst || !dst_size) return CARQUET_ERROR_INVALID_ARGUMENT;
    if (src_size == 0) { *dst_size = 0; return CARQUET_OK; }
    if (!src) return CARQUET_ERROR_INVALID_ARGUMENT;
    if (src_size > (size_t)INT_MAX || dst_capacity < 8) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    size_t body_cap = dst_capacity - 8;
    size_t body_size = 0;
    carquet_status_t s = carquet_lz4_compress(src, src_size, dst + 8,
                                              body_cap, &body_size);
    if (s != CARQUET_OK) return s;

    put_be32(dst, (uint32_t)src_size);         /* outer decompressed length */
    put_be32(dst + 4, (uint32_t)body_size);    /* inner compressed length */
    *dst_size = 8 + body_size;
    return CARQUET_OK;
}

carquet_status_t carquet_lz4_hadoop_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!dst || !dst_size) return CARQUET_ERROR_INVALID_ARGUMENT;
    if (src_size == 0) { *dst_size = 0; return CARQUET_OK; }
    if (!src) return CARQUET_ERROR_INVALID_ARGUMENT;

    size_t in_off = 0;
    size_t out_off = 0;

    while (in_off < src_size) {
        if (in_off + 4 > src_size) return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
        uint32_t outer_len = get_be32(src + in_off);
        in_off += 4;
        if (outer_len > dst_capacity - out_off) {
            return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
        }

        size_t outer_produced = 0;
        while (outer_produced < outer_len) {
            if (in_off + 4 > src_size) return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
            uint32_t comp_len = get_be32(src + in_off);
            in_off += 4;
            if (comp_len == 0 || in_off + comp_len > src_size ||
                comp_len > (size_t)INT_MAX) {
                return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
            }
            size_t remaining_out = dst_capacity - out_off;
            if (remaining_out > (size_t)INT_MAX) remaining_out = (size_t)INT_MAX;
            int r = LZ4_decompress_safe((const char*)(src + in_off),
                                        (char*)(dst + out_off),
                                        (int)comp_len, (int)remaining_out);
            if (r < 0) return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
            in_off += comp_len;
            out_off += (size_t)r;
            outer_produced += (size_t)r;
        }
        if (outer_produced != outer_len) {
            return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
        }
    }

    *dst_size = out_off;
    return CARQUET_OK;
}
