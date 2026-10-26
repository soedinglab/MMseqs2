/**
 * @file snappy.c
 * @brief C Snappy compression/decompression
 *
 * Based on Google's Snappy (BSD-3-Clause license).
 * C implementation of the Snappy format with NEON/SSSE3 SIMD support for
 * pattern extension in overlapping copies. Uses a fixed 16K hash table
 * (upstream uses adaptive 16K-32K); compressed output may differ byte-for-byte
 * from upstream but always decompresses to the same result.
 *
 * Reference: https://github.com/google/snappy/blob/main/format_description.txt
 *
 * Copyright 2005 Google Inc. All Rights Reserved.
 * Copyright 2024 carquet contributors.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright notice.
 *   * Redistributions in binary form must reproduce the above copyright notice
 *     in the documentation and/or other materials provided with the distribution.
 *   * Neither the name of Google Inc. nor the names of its contributors may be
 *     used to endorse or promote products derived from this software without
 *     specific prior written permission.
 */

#include <carquet/error.h>
#include "snappy.h"
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define SNAPPY_HAVE_NEON 1
#else
#define SNAPPY_HAVE_NEON 0
#endif

#if defined(__SSE2__)
#include <emmintrin.h>
#define SNAPPY_HAVE_SSE2 1
#else
#define SNAPPY_HAVE_SSE2 0
#endif

#if defined(__SSSE3__)
#include <tmmintrin.h>
#define SNAPPY_HAVE_SSSE3 1
#else
#define SNAPPY_HAVE_SSSE3 0
#endif

#if SNAPPY_HAVE_SSSE3 || SNAPPY_HAVE_NEON
#define SNAPPY_HAVE_VECTOR_SHUFFLE 1
#else
#define SNAPPY_HAVE_VECTOR_SHUFFLE 0
#endif

#if defined(__GNUC__) || defined(__clang__)
#define SNAPPY_PREDICT_TRUE(x)  __builtin_expect(!!(x), 1)
#define SNAPPY_PREDICT_FALSE(x) __builtin_expect(!!(x), 0)
#define SNAPPY_PREFETCH(addr)   __builtin_prefetch((addr), 0, 1)
#define SNAPPY_CTZ64(x)         __builtin_ctzll(x)
#define SNAPPY_CLZ32(x)         __builtin_clz(x)
#elif defined(_MSC_VER)
#include <intrin.h>
#define SNAPPY_PREDICT_TRUE(x)  (x)
#define SNAPPY_PREDICT_FALSE(x) (x)
#define SNAPPY_PREFETCH(addr)   ((void)0)
static inline int snappy_ctz64(uint64_t x) {
    unsigned long idx; _BitScanForward64(&idx, x); return (int)idx;
}
#define SNAPPY_CTZ64(x) snappy_ctz64(x)
static inline int snappy_clz32(uint32_t x) {
    unsigned long idx; _BitScanReverse(&idx, x); return 31 - (int)idx;
}
#define SNAPPY_CLZ32(x) snappy_clz32(x)
#else
#define SNAPPY_PREDICT_TRUE(x)  (x)
#define SNAPPY_PREDICT_FALSE(x) (x)
#define SNAPPY_PREFETCH(addr)   ((void)0)
#define SNAPPY_CTZ64(x)         snappy_ctz64_fallback(x)
#define SNAPPY_CLZ32(x)         snappy_clz32_fallback(x)
static inline int snappy_ctz64_fallback(uint64_t x) {
    int n = 0; while (!(x & 1)) { x >>= 1; n++; } return n;
}
static inline int snappy_clz32_fallback(uint32_t x) {
    int n = 0; while (!(x & 0x80000000u)) { x <<= 1; n++; } return n;
}
#endif

/* Tag types */
#define SNAPPY_LITERAL   0
#define SNAPPY_COPY_1    1
#define SNAPPY_COPY_2    2
#define SNAPPY_COPY_4    3

/* Compression constants */
#define SNAPPY_HASH_LOG  14
#define SNAPPY_HASH_SIZE (1 << SNAPPY_HASH_LOG)
#define SNAPPY_MAX_OFFSET 65535
#define SNAPPY_BLOCK_SIZE (1 << 16)

/* Slop bytes for unconditional copies in decompression */
#define SNAPPY_SLOP_BYTES 64

/* Forward declaration */
size_t carquet_snappy_compress_bound(size_t src_size);

/* ============================================================================
 * Unaligned load/store helpers
 * ============================================================================ */

static inline uint32_t load32(const void* p) {
    uint32_t v; memcpy(&v, p, 4); return v;
}

static inline uint64_t load64(const void* p) {
    uint64_t v; memcpy(&v, p, 8); return v;
}

static inline void store32(void* p, uint32_t v) {
    memcpy(p, &v, 4);
}

static inline void copy64(const void* src, void* dst) {
    uint64_t v; memcpy(&v, src, 8); memcpy(dst, &v, 8);
}

static inline void copy128(const void* src, void* dst) {
    uint64_t lo, hi;
    memcpy(&lo, src, 8);
    memcpy(&hi, (const char*)src + 8, 8);
    memcpy(dst, &lo, 8);
    memcpy((char*)dst + 8, &hi, 8);
}

/* ============================================================================
 * kLengthMinusOffset — Tag decode lookup table
 * Encodes length - (offset << 8) for copy-1/copy-2 length extraction.
 * From Google Snappy (BSD-3-Clause). Low byte = copy length.
 * ============================================================================ */

static const int16_t kLengthMinusOffset[256] = {
    /* Generated from: LengthMinusOffset(tag>>2, tag&3) for tag 0..255
     * Low byte = copy length. Used for fast copy-1/copy-2 length decode. */
     -255,     4,     1,   255,  -254,     5,     2,   255,
     -253,     6,     3,   255,  -252,     7,     4,   255,
     -251,     8,     5,   255,  -250,     9,     6,   255,
     -249,    10,     7,   255,  -248,    11,     8,   255,
     -247,  -252,     9,   255,  -246,  -251,    10,   255,
     -245,  -250,    11,   255,  -244,  -249,    12,   255,
     -243,  -248,    13,   255,  -242,  -247,    14,   255,
     -241,  -246,    15,   255,  -240,  -245,    16,   255,
     -239,  -508,    17,   255,  -238,  -507,    18,   255,
     -237,  -506,    19,   255,  -236,  -505,    20,   255,
     -235,  -504,    21,   255,  -234,  -503,    22,   255,
     -233,  -502,    23,   255,  -232,  -501,    24,   255,
     -231,  -764,    25,   255,  -230,  -763,    26,   255,
     -229,  -762,    27,   255,  -228,  -761,    28,   255,
     -227,  -760,    29,   255,  -226,  -759,    30,   255,
     -225,  -758,    31,   255,  -224,  -757,    32,   255,
     -223, -1020,    33,   255,  -222, -1019,    34,   255,
     -221, -1018,    35,   255,  -220, -1017,    36,   255,
     -219, -1016,    37,   255,  -218, -1015,    38,   255,
     -217, -1014,    39,   255,  -216, -1013,    40,   255,
     -215, -1276,    41,   255,  -214, -1275,    42,   255,
     -213, -1274,    43,   255,  -212, -1273,    44,   255,
     -211, -1272,    45,   255,  -210, -1271,    46,   255,
     -209, -1270,    47,   255,  -208, -1269,    48,   255,
     -207, -1532,    49,   255,  -206, -1531,    50,   255,
     -205, -1530,    51,   255,  -204, -1529,    52,   255,
     -203, -1528,    53,   255,  -202, -1527,    54,   255,
     -201, -1526,    55,   255,  -200, -1525,    56,   255,
     -199, -1788,    57,   255,  -198, -1787,    58,   255,
     -197, -1786,    59,   255,  -196, -1785,    60,   255,
      255, -1784,    61,   255,   255, -1783,    62,   255,
      255, -1782,    63,   255,   255, -1781,    64,   255,
};

/* ============================================================================
 * Varint Encoding/Decoding
 * ============================================================================ */

/**
 * Read a varint-encoded uint32. Matches upstream Google Snappy's
 * Varint::Parse32WithLimit: at most 5 bytes, and the 5th byte must
 * have value < 16 (i.e. contributes at most 4 bits at position 28,
 * keeping the result within uint32 range).
 */
static size_t snappy_read_varint(const uint8_t* p, const uint8_t* end, uint32_t* value) {
    uint32_t result = 0;
    const uint8_t* start = p;

    if (p >= end) return 0;
    uint8_t b = *p++; result = b & 0x7F;           if (b < 128) goto done;
    if (p >= end) return 0;
    b = *p++; result |= (uint32_t)(b & 0x7F) << 7; if (b < 128) goto done;
    if (p >= end) return 0;
    b = *p++; result |= (uint32_t)(b & 0x7F) << 14; if (b < 128) goto done;
    if (p >= end) return 0;
    b = *p++; result |= (uint32_t)(b & 0x7F) << 21; if (b < 128) goto done;
    if (p >= end) return 0;
    b = *p++; result |= (uint32_t)(b & 0x7F) << 28; if (b < 16) goto done;
    return 0; /* Overflow: 5th byte >= 16 would exceed uint32 */

done:
    *value = result;
    return (size_t)(p - start);
}

static size_t snappy_write_varint(uint8_t* p, uint32_t value) {
    uint8_t* start = p;
    while (value >= 0x80) {
        *p++ = (uint8_t)(value | 0x80);
        value >>= 7;
    }
    *p++ = (uint8_t)value;
    return (size_t)(p - start);
}

/* ============================================================================
 * SIMD Pattern Extension for Decompression
 *
 * Precomputed shuffle masks eliminate runtime modulo operations.
 * pattern_size ranges from 1..15 (the < 16 branch of incremental_copy).
 * Two tables: offset-0 masks (for initial load) and offset-16 masks (reshuffle).
 * ============================================================================ */

/* masks_offset0[ps][i] = i % ps, for ps = 1..15 */
static const uint8_t snappy_masks_offset0[16][16] = {
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ps=0 (unused) */
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ps=1 */
    {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}, /* ps=2 */
    {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0}, /* ps=3 */
    {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3}, /* ps=4 */
    {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0}, /* ps=5 */
    {0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3}, /* ps=6 */
    {0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1}, /* ps=7 */
    {0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7}, /* ps=8 */
    {0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6}, /* ps=9 */
    {0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5}, /* ps=10 */
    {0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4}, /* ps=11 */
    {0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3}, /* ps=12 */
    {0,1,2,3,4,5,6,7,8,9,10,11,12,0,1,2}, /* ps=13 */
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,1}, /* ps=14 */
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,0}, /* ps=15 */
};

/* masks_offset16[ps][i] = (16 + i) % ps, for ps = 1..15 */
static const uint8_t snappy_masks_offset16[16][16] = {
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ps=0 (unused) */
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* ps=1 */
    {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}, /* ps=2 */
    {1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1}, /* ps=3 */
    {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3}, /* ps=4 */
    {1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1}, /* ps=5 */
    {4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1}, /* ps=6 */
    {2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3}, /* ps=7 */
    {0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7}, /* ps=8 */
    {7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4}, /* ps=9 */
    {6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1}, /* ps=10 */
    {5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9}, /* ps=11 */
    {4,5,6,7,8,9,10,11,0,1,2,3,4,5,6,7}, /* ps=12 */
    {3,4,5,6,7,8,9,10,11,12,0,1,2,3,4,5}, /* ps=13 */
    {2,3,4,5,6,7,8,9,10,11,12,13,0,1,2,3}, /* ps=14 */
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,0,1}, /* ps=15 */
};

#if SNAPPY_HAVE_NEON

static inline uint8x16_t neon_load_pattern(const uint8_t* src, int pattern_size) {
    uint8x16_t gen_mask = vld1q_u8(snappy_masks_offset0[pattern_size]);
    uint8x16_t raw = vld1q_u8(src);
    return vqtbl1q_u8(raw, gen_mask);
}

static inline uint8x16_t neon_reshuffle_mask(int pattern_size) {
    return vld1q_u8(snappy_masks_offset16[pattern_size]);
}

#elif SNAPPY_HAVE_SSSE3

static inline __m128i ssse3_load_pattern(const uint8_t* src, int pattern_size) {
    __m128i gen_mask = _mm_loadu_si128((const __m128i*)snappy_masks_offset0[pattern_size]);
    __m128i raw = _mm_loadu_si128((const __m128i*)src);
    return _mm_shuffle_epi8(raw, gen_mask);
}

static inline __m128i ssse3_reshuffle_mask(int pattern_size) {
    return _mm_loadu_si128((const __m128i*)snappy_masks_offset16[pattern_size]);
}

#endif

/* ============================================================================
 * IncrementalCopy — Overlapping copy for match expansion
 * ============================================================================ */

static inline uint8_t* incremental_copy_slow(const uint8_t* src, uint8_t* op,
                                              uint8_t* const op_limit) {
    while (op < op_limit) *op++ = *src++;
    return op_limit;
}

static inline uint8_t* incremental_copy(const uint8_t* src, uint8_t* op,
                                         uint8_t* const op_limit,
                                         uint8_t* const buf_limit) {
    size_t pattern_size = (size_t)(op - src);

    /* The pattern_size >= big_pattern block below extends the match with 16-byte
     * copy128 reads from `src`, so it is only correct once at least 16 valid
     * pattern bytes precede `op` — i.e. big_pattern must be 16 on every build.
     * With big_pattern == 8, pattern_size in [8,15] reaches copy128 and reads
     * 16 - pattern_size bytes past `op` (uninitialised output), corrupting the
     * result. Those sizes are instead served safely by the scalar doubling/
     * 8-byte-copy path, which the SIMD build already exercises via its own 16. */
    const int big_pattern = 16;

    if (pattern_size < (size_t)big_pattern) {
#if SNAPPY_HAVE_VECTOR_SHUFFLE
        if (SNAPPY_PREDICT_TRUE(op_limit <= buf_limit - 15)) {
#if SNAPPY_HAVE_NEON
            uint8x16_t pattern = neon_load_pattern(src, (int)pattern_size);
            uint8x16_t reshuffle = neon_reshuffle_mask((int)pattern_size);
            vst1q_u8(op, pattern);
            if (op + 16 < op_limit) {
                pattern = vqtbl1q_u8(pattern, reshuffle);
                vst1q_u8(op + 16, pattern);
            }
            if (op + 32 < op_limit) {
                pattern = vqtbl1q_u8(pattern, reshuffle);
                vst1q_u8(op + 32, pattern);
            }
            if (op + 48 < op_limit) {
                pattern = vqtbl1q_u8(pattern, reshuffle);
                vst1q_u8(op + 48, pattern);
            }
#else
            __m128i pattern = ssse3_load_pattern(src, (int)pattern_size);
            __m128i reshuffle = ssse3_reshuffle_mask((int)pattern_size);
            _mm_storeu_si128((__m128i*)op, pattern);
            if (op + 16 < op_limit) {
                pattern = _mm_shuffle_epi8(pattern, reshuffle);
                _mm_storeu_si128((__m128i*)(op + 16), pattern);
            }
            if (op + 32 < op_limit) {
                pattern = _mm_shuffle_epi8(pattern, reshuffle);
                _mm_storeu_si128((__m128i*)(op + 32), pattern);
            }
            if (op + 48 < op_limit) {
                pattern = _mm_shuffle_epi8(pattern, reshuffle);
                _mm_storeu_si128((__m128i*)(op + 48), pattern);
            }
#endif
            return op_limit;
        }
        return incremental_copy_slow(src, op, op_limit);
#else  /* !SNAPPY_HAVE_VECTOR_SHUFFLE */
        /* Non-SIMD: expand pattern to at least 8 bytes by doubling */
        if (SNAPPY_PREDICT_TRUE(op <= buf_limit - 11)) {
            while (pattern_size < 8) {
                copy64(src, op);
                op += pattern_size;
                pattern_size *= 2;
            }
            if (SNAPPY_PREDICT_TRUE(op >= op_limit)) return op_limit;
            /* Pattern is now 8 bytes wide — use 8-byte block copies.
               We must NOT fall through to copy128 since only 8 bytes of the
               pattern are valid; a 16-byte read would pick up garbage. */
            src = op - pattern_size;
            while (op + 8 <= op_limit && op + 8 <= buf_limit) {
                copy64(src, op);
                src += 8;
                op += 8;
            }
            if (op >= op_limit) return op_limit;
            return incremental_copy_slow(src, op, op_limit);
        } else {
            return incremental_copy_slow(src, op, op_limit);
        }
#endif
    }

    /* pattern_size >= big_pattern (>= 16 with SIMD): simple block copies */
    if (SNAPPY_PREDICT_TRUE(op_limit <= buf_limit - 15)) {
        copy128(src, op);
        if (op + 16 < op_limit) copy128(src + 16, op + 16);
        if (op + 32 < op_limit) copy128(src + 32, op + 32);
        if (op + 48 < op_limit) copy128(src + 48, op + 48);
        return op_limit;
    }

    /* Near end of buffer: 16-byte copies until we run out of slop */
    {
        uint8_t* op_end = buf_limit - 16;
        while (op < op_end) {
            copy128(src, op);
            op += 16;
            src += 16;
        }
        if (op >= op_limit) return op_limit;
    }

    if (SNAPPY_PREDICT_FALSE(op <= buf_limit - 8)) {
        copy64(src, op);
        src += 8;
        op += 8;
    }
    return incremental_copy_slow(src, op, op_limit);
}

/* ============================================================================
 * Snappy Decompression
 * ============================================================================ */

const char* carquet_snappy_reason_string(carquet_snappy_reason_t reason) {
    switch (reason) {
        case CARQUET_SNAPPY_OK:
            return "ok";
        case CARQUET_SNAPPY_ERR_LENGTH_VARINT:
            return "malformed uncompressed-length varint";
        case CARQUET_SNAPPY_ERR_OUTPUT_TOO_SMALL:
            return "declared uncompressed length exceeds output capacity";
        case CARQUET_SNAPPY_ERR_LITERAL_LEN_TRUNC:
            return "long-literal length bytes run past end of input";
        case CARQUET_SNAPPY_ERR_LITERAL_INPUT:
            return "literal runs past end of input";
        case CARQUET_SNAPPY_ERR_LITERAL_OUTPUT:
            return "literal runs past end of output";
        case CARQUET_SNAPPY_ERR_COPY1_TRUNC:
            return "COPY_1 trailer byte missing";
        case CARQUET_SNAPPY_ERR_COPY2_TRUNC:
            return "COPY_2 trailer bytes missing";
        case CARQUET_SNAPPY_ERR_COPY4_TRUNC:
            return "COPY_4 trailer bytes missing";
        case CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO:
            return "copy offset is zero";
        case CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE:
            return "copy offset points before start of output";
        case CARQUET_SNAPPY_ERR_COPY_OUTPUT:
            return "copy runs past end of output";
        case CARQUET_SNAPPY_ERR_SHORT_OUTPUT:
            return "input exhausted before the declared uncompressed length";
        case CARQUET_SNAPPY_ERR_TRAILING_INPUT:
            return "output complete but input bytes remain";
        default:
            return "unknown";
    }
}

const char* carquet_snappy_diag_format(const carquet_snappy_diag_t* diag,
                                       char* buf, size_t buf_size) {
    if (!buf || buf_size == 0) return buf;
    if (!diag) {
        snprintf(buf, buf_size, "no snappy diagnostic");
        return buf;
    }
    /* Kept compact on purpose: callers embed this in a carquet_error_t message,
     * which is capped at CARQUET_ERROR_MESSAGE_MAX (256) bytes including their
     * own prefix. %llu, not %zu — the MinGW/msvcrt printf has no z modifier. */
    snprintf(buf, buf_size,
             "snappy[%d]: %s [src=%llu dst=%llu declared=%llu tag=%d op=%llu]",
             (int)diag->reason, carquet_snappy_reason_string(diag->reason),
             (unsigned long long)diag->input_pos,
             (unsigned long long)diag->output_pos,
             (unsigned long long)diag->declared_length,
             diag->tag, (unsigned long long)diag->operand);
    return buf;
}

/* Record the failing check and bail out of the decode loop. */
#define SNAPPY_FAIL(reason_code, operand_value)                       \
    do {                                                              \
        if (diag) {                                                   \
            diag->reason = (reason_code);                             \
            diag->input_pos = (size_t)(tag_start - src);              \
            diag->output_pos = (size_t)(op - dst);                    \
            diag->declared_length = uncompressed_len;                 \
            diag->tag = (int)tag;                                     \
            diag->operand = (uint64_t)(operand_value);                \
        }                                                             \
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;                 \
    } while (0)

carquet_status_t carquet_snappy_decompress_diag(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size,
    carquet_snappy_diag_t* diag) {

    if (diag) {
        memset(diag, 0, sizeof(*diag));
        diag->tag = -1;
    }

    if (!src || !dst || !dst_size)
        return CARQUET_ERROR_INVALID_ARGUMENT;

    if (src_size == 0) {
        *dst_size = 0;
        return CARQUET_OK;
    }

    const uint8_t* ip = src;
    const uint8_t* const ip_end = src + src_size;

    /* Read uncompressed length */
    uint32_t uncompressed_len;
    size_t varint_len = snappy_read_varint(ip, ip_end, &uncompressed_len);
    if (varint_len == 0) {
        if (diag) {
            diag->reason = CARQUET_SNAPPY_ERR_LENGTH_VARINT;
            diag->tag = -1;
            diag->operand = src_size;
        }
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }
    ip += varint_len;

    if (uncompressed_len > dst_capacity) {
        if (diag) {
            diag->reason = CARQUET_SNAPPY_ERR_OUTPUT_TOO_SMALL;
            diag->input_pos = varint_len;
            diag->declared_length = uncompressed_len;
            diag->tag = -1;
            diag->operand = dst_capacity;
        }
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }

    uint8_t* op = dst;
    uint8_t* const op_end = dst + uncompressed_len;
    /* For safe SIMD writes, we need slop at the end */
    uint8_t* const op_limit_min_slop = (uncompressed_len >= SNAPPY_SLOP_BYTES)
        ? (op_end - SNAPPY_SLOP_BYTES + 1) : dst;

    while (ip < ip_end && op < op_end) {
        const uint8_t* const tag_start = ip;
        const uint8_t tag = *ip++;
        const uint8_t type = tag & 0x03;

        if (type == SNAPPY_LITERAL) {
            size_t literal_len = (tag >> 2) + 1;
            if (SNAPPY_PREDICT_FALSE(literal_len >= 61)) {
                /* Long literal: length is encoded in 1-4 following bytes */
                size_t extra_bytes = literal_len - 60;
                if (SNAPPY_PREDICT_FALSE(extra_bytes > (size_t)(ip_end - ip)))
                    SNAPPY_FAIL(CARQUET_SNAPPY_ERR_LITERAL_LEN_TRUNC, extra_bytes);
                /* Use a 32-bit load and mask (like Google Snappy) */
                uint32_t raw = 0;
                memcpy(&raw, ip, extra_bytes <= 4 ? extra_bytes : 4);
                /* Mask to the relevant bytes */
                uint64_t mask64 = 0xFFFFFFFF;
                literal_len = (raw & (uint32_t)~(mask64 << (8 * extra_bytes))) + 1;
                ip += extra_bytes;
            }

            /* Fast path for short literals with enough room */
            if (SNAPPY_PREDICT_TRUE(literal_len <= 16 &&
                                    ip + 16 <= ip_end &&
                                    op + 16 <= op_end)) {
                copy128(ip, op);
                ip += literal_len;
                op += literal_len;
                continue;
            }

            if (SNAPPY_PREDICT_FALSE(literal_len > (size_t)(ip_end - ip)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_LITERAL_INPUT, literal_len);
            if (SNAPPY_PREDICT_FALSE(literal_len > (size_t)(op_end - op)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_LITERAL_OUTPUT, literal_len);
            memcpy(op, ip, literal_len);
            ip += literal_len;
            op += literal_len;

        } else if (SNAPPY_PREDICT_TRUE(type != SNAPPY_COPY_4)) {
            /* COPY_1 or COPY_2 — use kLengthMinusOffset for branchless decode */
            int16_t entry = kLengthMinusOffset[tag];
            uint32_t trailer;
            size_t length;
            size_t copy_offset;

            if (type == SNAPPY_COPY_1) {
                if (SNAPPY_PREDICT_FALSE(ip >= ip_end))
                    SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY1_TRUNC, 1);
                trailer = ((uint32_t)(tag & 0xE0) << 3) | *ip++;
                length = (size_t)(entry & 0xFF);
                copy_offset = trailer;
            } else { /* SNAPPY_COPY_2 */
                if (SNAPPY_PREDICT_FALSE((size_t)(ip_end - ip) < 2))
                    SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY2_TRUNC, 2);
                trailer = (uint32_t)ip[0] | ((uint32_t)ip[1] << 8);
                ip += 2;
                length = (size_t)(entry & 0xFF);
                copy_offset = trailer;
            }

            if (SNAPPY_PREDICT_FALSE(copy_offset == 0))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO, 0);
            if (SNAPPY_PREDICT_FALSE(copy_offset > (size_t)(op - dst)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE, copy_offset);
            if (SNAPPY_PREDICT_FALSE(length > (size_t)(op_end - op)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OUTPUT, length);

            const uint8_t* match_src = op - copy_offset;

            /* Fast path: offset >= 16, just copy non-overlapping blocks */
            if (SNAPPY_PREDICT_TRUE(copy_offset >= 16 && op + length <= op_limit_min_slop)) {
                copy128(match_src, op);
                if (length > 16) copy128(match_src + 16, op + 16);
                if (length > 32) copy128(match_src + 32, op + 32);
                if (length > 48) copy128(match_src + 48, op + 48);
                op += length;
            } else if (SNAPPY_PREDICT_TRUE(length <= SNAPPY_SLOP_BYTES &&
                                           op + length <= op_limit_min_slop &&
                                           copy_offset >= length)) {
                /* Non-overlapping but small offset: use memmove */
                memmove(op, match_src, SNAPPY_SLOP_BYTES);
                op += length;
            } else {
                (void)incremental_copy(match_src, op, op + length, op_end);
                op += length;
            }

        } else {
            /* COPY_4: 4-byte offset (rare) */
            if (SNAPPY_PREDICT_FALSE((size_t)(ip_end - ip) < 4))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY4_TRUNC, 4);
            size_t length = ((tag >> 2) & 0x3F) + 1;
            size_t copy_offset = (size_t)load32(ip);
            ip += 4;
            if (SNAPPY_PREDICT_FALSE(copy_offset == 0))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO, 0);
            if (SNAPPY_PREDICT_FALSE(copy_offset > (size_t)(op - dst)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE, copy_offset);
            if (SNAPPY_PREDICT_FALSE(length > (size_t)(op_end - op)))
                SNAPPY_FAIL(CARQUET_SNAPPY_ERR_COPY_OUTPUT, length);
            const uint8_t* match_src = op - copy_offset;
            (void)incremental_copy(match_src, op, op + length, op_end);
            op += length;
        }
    }

    /* Upstream requires BOTH output length match AND full input consumption.
     * Without the ip check, trailing garbage after valid data is accepted. */
    if ((size_t)(op - dst) != uncompressed_len) {
        if (diag) {
            diag->reason = CARQUET_SNAPPY_ERR_SHORT_OUTPUT;
            diag->input_pos = (size_t)(ip - src);
            diag->output_pos = (size_t)(op - dst);
            diag->declared_length = uncompressed_len;
            diag->tag = -1;
            diag->operand = (uint64_t)(ip_end - ip);
        }
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }
    if (ip != ip_end) {
        if (diag) {
            diag->reason = CARQUET_SNAPPY_ERR_TRAILING_INPUT;
            diag->input_pos = (size_t)(ip - src);
            diag->output_pos = (size_t)(op - dst);
            diag->declared_length = uncompressed_len;
            diag->tag = -1;
            diag->operand = (uint64_t)(ip_end - ip);
        }
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;
    }

    *dst_size = uncompressed_len;
    return CARQUET_OK;
}

#undef SNAPPY_FAIL

carquet_status_t carquet_snappy_decompress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {
    return carquet_snappy_decompress_diag(src, src_size, dst, dst_capacity,
                                          dst_size, NULL);
}

/* ============================================================================
 * Snappy Compression
 * ============================================================================ */

static inline uint32_t snappy_hash(uint32_t val) {
    return (val * 0x1e35a7bd) >> (32 - SNAPPY_HASH_LOG);
}

/* Fast match length using 64-bit XOR comparison */
static inline size_t fast_match_length(const uint8_t* p, const uint8_t* match,
                                        const uint8_t* limit) {
    const uint8_t* start = p;
    while (p + 8 <= limit) {
        uint64_t a = load64(p);
        uint64_t b = load64(match);
        uint64_t diff = a ^ b;
        if (diff)
            return (size_t)(p - start) + ((size_t)SNAPPY_CTZ64(diff) >> 3);
        p += 8;
        match += 8;
    }
    while (p < limit && *p == *match) { p++; match++; }
    return (size_t)(p - start);
}

static uint8_t* snappy_emit_literal(uint8_t* op, const uint8_t* literal, size_t len) {
    size_t n = len - 1;
    if (n < 60) {
        *op++ = (uint8_t)(n << 2);
    } else {
        /* Encode length in 1-4 extra bytes, like Google Snappy */
        int count = (31 - SNAPPY_CLZ32((uint32_t)n)) / 8 + 1;
        *op++ = (uint8_t)((59 + count) << 2);
        store32(op, (uint32_t)n);
        op += count;
    }
    memcpy(op, literal, len);
    return op + len;
}

static inline uint8_t* snappy_emit_copy(uint8_t* op, size_t offset, size_t len) {
    /* Emit 64-byte chunks */
    while (SNAPPY_PREDICT_FALSE(len >= 68)) {
        *op++ = (uint8_t)((63 << 2) | SNAPPY_COPY_2);
        *op++ = (uint8_t)(offset & 0xFF);
        *op++ = (uint8_t)(offset >> 8);
        len -= 64;
    }

    if (len > 64) {
        *op++ = (uint8_t)((59 << 2) | SNAPPY_COPY_2);
        *op++ = (uint8_t)(offset & 0xFF);
        *op++ = (uint8_t)(offset >> 8);
        len -= 60;
    }

    /* Branchless offset type selection (like Google Snappy) */
    if (len < 12 && offset < 2048) {
        /* 1-byte offset copy */
        uint32_t u = ((uint32_t)len << 2) + ((uint32_t)offset << 8);
        uint32_t copy1 = SNAPPY_COPY_1 - (4 << 2) + (((uint32_t)offset >> 3) & 0xe0);
        u += copy1;
        store32(op, u);
        op += 2;
    } else if (len < 12) {
        /* 2-byte offset copy for small length, large offset */
        uint32_t u = SNAPPY_COPY_2 + (((uint32_t)len - 1) << 2) + ((uint32_t)offset << 8);
        store32(op, u);
        op += 3;
    } else {
        /* 2-byte offset copy */
        *op++ = (uint8_t)(((len - 1) << 2) | SNAPPY_COPY_2);
        *op++ = (uint8_t)(offset & 0xFF);
        *op++ = (uint8_t)(offset >> 8);
    }

    return op;
}

static uint8_t* snappy_compress_block(
    const uint8_t* src,
    size_t src_size,
    uint8_t* op) {

    if (src_size == 0) return op;
    if (src_size < 15) return snappy_emit_literal(op, src, src_size);

    uint16_t hash_table[SNAPPY_HASH_SIZE];
    memset(hash_table, 0, sizeof(hash_table));

    const uint8_t* const iend = src + src_size;
    const uint8_t* const ilimit = iend - 15;
    const uint8_t* ip = src + 1;
    const uint8_t* anchor = src;
    const uint8_t* candidate;

    /* Pre-seed position 0 */
    hash_table[snappy_hash(load32(src))] = 0;

    /* Try to match 16 bytes at positions 0..15 for fast startup */
    if (ilimit - ip >= 16) {
        uint64_t data = load64(ip);
        ptrdiff_t delta = (ptrdiff_t)(ip - src);
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                int i = 4 * j + k;
                uint32_t dword = (i == 0) ? load32(ip) : (uint32_t)data;
                uint32_t h = snappy_hash(dword);
                candidate = src + hash_table[h];
                hash_table[h] = (uint16_t)(delta + i);
                if (SNAPPY_PREDICT_FALSE(load32(candidate) == dword)) {
                    *op = (uint8_t)(SNAPPY_LITERAL | (i << 2));
                    copy128(anchor, op + 1);
                    ip += i;
                    op = op + i + 2;
                    goto emit_match;
                }
                data >>= 8;
            }
            data = load64(ip + 4 * j + 4);
        }
        ip += 16;
    }

    {
        uint32_t skip = 32;
        for (;;) {
            uint32_t h = snappy_hash(load32(ip));
            uint32_t bytes_between = skip >> 5;
            skip += bytes_between;
            const uint8_t* next_ip = ip + bytes_between;

            if (SNAPPY_PREDICT_FALSE(next_ip > ilimit)) {
                ip = anchor;
                goto emit_remainder;
            }

            candidate = src + hash_table[h];
            hash_table[h] = (uint16_t)(ip - src);

            if (SNAPPY_PREDICT_FALSE(load32(ip) == load32(candidate)))
                break;

            ip = next_ip;
        }
    }

    /* Emit pending literal */
    if (ip > anchor)
        op = snappy_emit_literal(op, anchor, (size_t)(ip - anchor));

emit_match:
    do {
        size_t match_len = 4 + fast_match_length(ip + 4, candidate + 4, iend);
        size_t offset = (size_t)(ip - candidate);
        ip += match_len;
        op = snappy_emit_copy(op, offset, match_len);

        if (SNAPPY_PREDICT_FALSE(ip >= ilimit)) {
            anchor = ip;
            goto emit_remainder;
        }

        /* Insert hash entries near match end */
        hash_table[snappy_hash(load32(ip - 1))] = (uint16_t)(ip - 1 - src);
        uint32_t h = snappy_hash(load32(ip));
        candidate = src + hash_table[h];
        hash_table[h] = (uint16_t)(ip - src);
    } while (load32(ip) == load32(candidate) &&
             (size_t)(ip - candidate) <= SNAPPY_MAX_OFFSET);

    anchor = ip++;
    {
        uint32_t skip = 32;
        for (;;) {
            uint32_t h = snappy_hash(load32(ip));
            uint32_t bytes_between = skip >> 5;
            skip += bytes_between;
            const uint8_t* next_ip = ip + bytes_between;

            if (SNAPPY_PREDICT_FALSE(next_ip > ilimit))
                goto emit_remainder;

            candidate = src + hash_table[h];
            hash_table[h] = (uint16_t)(ip - src);

            if (SNAPPY_PREDICT_FALSE(load32(ip) == load32(candidate))) {
                if (ip > anchor)
                    op = snappy_emit_literal(op, anchor, (size_t)(ip - anchor));
                goto emit_match;
            }

            ip = next_ip;
        }
    }

emit_remainder:
    if (anchor < iend)
        op = snappy_emit_literal(op, anchor, (size_t)(iend - anchor));
    return op;
}

carquet_status_t carquet_snappy_compress(
    const uint8_t* src,
    size_t src_size,
    uint8_t* dst,
    size_t dst_capacity,
    size_t* dst_size) {

    if (!dst || !dst_size)
        return CARQUET_ERROR_INVALID_ARGUMENT;

    size_t max_output = carquet_snappy_compress_bound(src_size);
    if (dst_capacity < max_output)
        return CARQUET_ERROR_COMPRESSION;

    uint8_t* op = dst;
    op += snappy_write_varint(op, (uint32_t)src_size);

    if (src_size == 0) {
        *dst_size = (size_t)(op - dst);
        return CARQUET_OK;
    }

    if (!src)
        return CARQUET_ERROR_INVALID_ARGUMENT;

    size_t pos = 0;
    while (pos < src_size) {
        size_t block_size = src_size - pos;
        if (block_size > SNAPPY_BLOCK_SIZE)
            block_size = SNAPPY_BLOCK_SIZE;
        op = snappy_compress_block(src + pos, block_size, op);
        pos += block_size;
    }

    *dst_size = (size_t)(op - dst);
    return CARQUET_OK;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

size_t carquet_snappy_compress_bound(size_t src_size) {
    return 32 + src_size + src_size / 6;
}

carquet_status_t carquet_snappy_get_uncompressed_length(
    const uint8_t* src,
    size_t src_size,
    size_t* length) {

    if (!src || !length)
        return CARQUET_ERROR_INVALID_ARGUMENT;

    uint32_t len;
    size_t varint_len = snappy_read_varint(src, src + src_size, &len);
    if (varint_len == 0)
        return CARQUET_ERROR_INVALID_COMPRESSED_DATA;

    *length = len;
    return CARQUET_OK;
}
