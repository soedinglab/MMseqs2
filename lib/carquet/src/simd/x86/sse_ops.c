/**
 * @file sse_ops.c
 * @brief SSE4.2 optimized operations for x86 processors
 *
 * Provides SIMD-accelerated implementations of:
 * - Bit unpacking for common bit widths
 * - Byte stream split/merge (for BYTE_STREAM_SPLIT encoding)
 * - Delta decoding (prefix sums)
 * - Dictionary gather operations
 */

#include <carquet/error.h>
#include "simd/simd_unaligned.h"
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
/* SSE4.2 is always available on x64 MSVC, check __SSE4_2__ for GCC/Clang */
#if defined(__SSE4_2__) || defined(_M_X64) || defined(_M_IX86)

#ifdef _MSC_VER
#include <intrin.h>

/* MSVC doesn't have __builtin_prefetch, use _mm_prefetch instead */
#define __builtin_prefetch(addr, rw, locality) \
    _mm_prefetch((const char*)(addr), _MM_HINT_T0)

/* MSVC doesn't have __builtin_ctz (count trailing zeros) */
static inline int msvc_ctz(unsigned int x) {
    unsigned long index;
    _BitScanForward(&index, x);
    return (int)index;
}
#define __builtin_ctz(x) msvc_ctz(x)

/* MSVC doesn't have __builtin_popcount */
#define __builtin_popcount(x) __popcnt(x)

#endif
#include <smmintrin.h>
#include <nmmintrin.h>

static inline uint16_t sse_read_le16(const uint8_t* p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t sse_read_le24(const uint8_t* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16);
}

static inline uint64_t sse_read_le40(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32);
}

static inline uint64_t sse_read_le48(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40);
}

static inline uint64_t sse_read_le56(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48);
}

/* ============================================================================
 * Bit Unpacking - SSE Optimized
 * ============================================================================
 */

/**
 * Unpack 8 1-bit values using SSE.
 */
void carquet_sse_bitunpack8_1bit(const uint8_t* input, uint32_t* values) {
    __m128i bytes = _mm_set1_epi8((char)input[0]);
    const __m128i bit_mask = _mm_setr_epi8(
        0x01, 0x02, 0x04, 0x08,
        0x10, 0x20, 0x40, (char)0x80,
        0, 0, 0, 0, 0, 0, 0, 0
    );
    __m128i masked = _mm_and_si128(bytes, bit_mask);
    __m128i cmp = _mm_cmpeq_epi8(masked, bit_mask);
    __m128i ones = _mm_set1_epi8(1);
    __m128i result8 = _mm_and_si128(cmp, ones);
    __m128i zero = _mm_setzero_si128();
    __m128i words = _mm_unpacklo_epi8(result8, zero);
    __m128i v0 = _mm_unpacklo_epi16(words, zero);
    __m128i v1 = _mm_unpackhi_epi16(words, zero);

    _mm_storeu_si128((__m128i*)(values + 0), v0);
    _mm_storeu_si128((__m128i*)(values + 4), v1);
}

void carquet_sse_bitunpack8_2bit(const uint8_t* input, uint32_t* values) {
    uint16_t v = sse_read_le16(input);
    __m128i lo = _mm_setr_epi32(
        (int)((v >> 0) & 0x3), (int)((v >> 2) & 0x3),
        (int)((v >> 4) & 0x3), (int)((v >> 6) & 0x3));
    __m128i hi = _mm_setr_epi32(
        (int)((v >> 8) & 0x3), (int)((v >> 10) & 0x3),
        (int)((v >> 12) & 0x3), (int)((v >> 14) & 0x3));
    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

void carquet_sse_bitunpack8_3bit(const uint8_t* input, uint32_t* values) {
    uint32_t v = sse_read_le24(input);
    __m128i lo = _mm_setr_epi32(
        (int)((v >> 0) & 0x7), (int)((v >> 3) & 0x7),
        (int)((v >> 6) & 0x7), (int)((v >> 9) & 0x7));
    __m128i hi = _mm_setr_epi32(
        (int)((v >> 12) & 0x7), (int)((v >> 15) & 0x7),
        (int)((v >> 18) & 0x7), (int)((v >> 21) & 0x7));
    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

/**
 * Unpack 32 1-bit values using SSE.
 * Input: 4 bytes, Output: 32 x uint32_t
 */
void carquet_sse_bitunpack32_1bit(const uint8_t* input, uint32_t* values) {
    /* Load 4 bytes and expand */
    __m128i bytes = _mm_cvtsi32_si128(*(const int32_t*)input);

    /* Shuffle to repeat each byte 8 times for masking */
    static const int8_t shuffle_mask[16] = {
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1
    };
    __m128i shuf = _mm_loadu_si128((const __m128i*)shuffle_mask);

    /* Process first 16 bits */
    __m128i expanded = _mm_shuffle_epi8(bytes, shuf);
    __m128i bit_mask = _mm_set_epi8(
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01,
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01
    );

    __m128i masked = _mm_and_si128(expanded, bit_mask);
    __m128i result = _mm_min_epu8(masked, _mm_set1_epi8(1));

    /* Unpack to 32-bit */
    __m128i zero = _mm_setzero_si128();
    __m128i lo8 = _mm_unpacklo_epi8(result, zero);
    __m128i hi8 = _mm_unpackhi_epi8(result, zero);

    __m128i v0 = _mm_unpacklo_epi16(lo8, zero);
    __m128i v1 = _mm_unpackhi_epi16(lo8, zero);
    __m128i v2 = _mm_unpacklo_epi16(hi8, zero);
    __m128i v3 = _mm_unpackhi_epi16(hi8, zero);

    _mm_storeu_si128((__m128i*)(values + 0), v0);
    _mm_storeu_si128((__m128i*)(values + 4), v1);
    _mm_storeu_si128((__m128i*)(values + 8), v2);
    _mm_storeu_si128((__m128i*)(values + 12), v3);

    /* Process bytes 2-3 for values 16-31 */
    static const int8_t shuffle_mask2[16] = {
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3
    };
    shuf = _mm_loadu_si128((const __m128i*)shuffle_mask2);
    expanded = _mm_shuffle_epi8(bytes, shuf);
    masked = _mm_and_si128(expanded, bit_mask);
    result = _mm_min_epu8(masked, _mm_set1_epi8(1));

    lo8 = _mm_unpacklo_epi8(result, zero);
    hi8 = _mm_unpackhi_epi8(result, zero);

    v0 = _mm_unpacklo_epi16(lo8, zero);
    v1 = _mm_unpackhi_epi16(lo8, zero);
    v2 = _mm_unpacklo_epi16(hi8, zero);
    v3 = _mm_unpackhi_epi16(hi8, zero);

    _mm_storeu_si128((__m128i*)(values + 16), v0);
    _mm_storeu_si128((__m128i*)(values + 20), v1);
    _mm_storeu_si128((__m128i*)(values + 24), v2);
    _mm_storeu_si128((__m128i*)(values + 28), v3);
}

/**
 * Unpack 8 4-bit values using SSE.
 */
void carquet_sse_bitunpack8_4bit(const uint8_t* input, uint32_t* values) {
    /* Load 4 bytes containing 8 x 4-bit values */
    __m128i bytes = _mm_cvtsi32_si128(*(const int32_t*)input);

    /* Split into low and high nibbles */
    __m128i lo_nibbles = _mm_and_si128(bytes, _mm_set1_epi8(0x0F));
    __m128i hi_nibbles = _mm_srli_epi16(bytes, 4);
    hi_nibbles = _mm_and_si128(hi_nibbles, _mm_set1_epi8(0x0F));

    /* Interleave */
    __m128i interleaved = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);

    /* Expand to 32-bit */
    __m128i zero = _mm_setzero_si128();
    __m128i words = _mm_unpacklo_epi8(interleaved, zero);

    __m128i v0 = _mm_unpacklo_epi16(words, zero);
    __m128i v1 = _mm_unpackhi_epi16(words, zero);

    _mm_storeu_si128((__m128i*)(values + 0), v0);
    _mm_storeu_si128((__m128i*)(values + 4), v1);
}

void carquet_sse_bitunpack8_5bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = sse_read_le40(input);
    __m128i lo = _mm_setr_epi32(
        (int)((v >> 0) & 0x1F), (int)((v >> 5) & 0x1F),
        (int)((v >> 10) & 0x1F), (int)((v >> 15) & 0x1F));
    __m128i hi = _mm_setr_epi32(
        (int)((v >> 20) & 0x1F), (int)((v >> 25) & 0x1F),
        (int)((v >> 30) & 0x1F), (int)((v >> 35) & 0x1F));
    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

void carquet_sse_bitunpack8_6bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = sse_read_le48(input);
    __m128i lo = _mm_setr_epi32(
        (int)((v >> 0) & 0x3F), (int)((v >> 6) & 0x3F),
        (int)((v >> 12) & 0x3F), (int)((v >> 18) & 0x3F));
    __m128i hi = _mm_setr_epi32(
        (int)((v >> 24) & 0x3F), (int)((v >> 30) & 0x3F),
        (int)((v >> 36) & 0x3F), (int)((v >> 42) & 0x3F));
    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

void carquet_sse_bitunpack8_7bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = sse_read_le56(input);
    __m128i lo = _mm_setr_epi32(
        (int)((v >> 0) & 0x7F), (int)((v >> 7) & 0x7F),
        (int)((v >> 14) & 0x7F), (int)((v >> 21) & 0x7F));
    __m128i hi = _mm_setr_epi32(
        (int)((v >> 28) & 0x7F), (int)((v >> 35) & 0x7F),
        (int)((v >> 42) & 0x7F), (int)((v >> 49) & 0x7F));
    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

/**
 * Unpack 8 8-bit values using SSE (widen u8 to u32).
 */
void carquet_sse_bitunpack8_8bit(const uint8_t* input, uint32_t* values) {
    /* Load 8 bytes */
    __m128i bytes = _mm_loadl_epi64((const __m128i*)input);

    /* Expand to 32-bit */
    __m128i zero = _mm_setzero_si128();
    __m128i words = _mm_unpacklo_epi8(bytes, zero);

    __m128i v0 = _mm_unpacklo_epi16(words, zero);
    __m128i v1 = _mm_unpackhi_epi16(words, zero);

    _mm_storeu_si128((__m128i*)(values + 0), v0);
    _mm_storeu_si128((__m128i*)(values + 4), v1);
}

/**
 * Unpack 8 16-bit values using SSE.
 */
void carquet_sse_bitunpack8_16bit(const uint8_t* input, uint32_t* values) {
    __m128i words = _mm_loadu_si128((const __m128i*)input);
    __m128i lo = _mm_cvtepu16_epi32(words);
    __m128i hi = _mm_cvtepu16_epi32(_mm_srli_si128(words, 8));

    _mm_storeu_si128((__m128i*)(values + 0), lo);
    _mm_storeu_si128((__m128i*)(values + 4), hi);
}

/* ============================================================================
 * Byte Stream Split - SSE Optimized
 * ============================================================================
 */

/**
 * Encode floats using byte stream split with SSE.
 * Transposes: puts all byte 0s together, then all byte 1s, etc.
 */
void carquet_sse_byte_stream_split_encode_float(
    const float* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;
    const __m128i s0 = _mm_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s1 = _mm_setr_epi8(1, 5, 9, 13, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s2 = _mm_setr_epi8(2, 6, 10, 14, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s3 = _mm_setr_epi8(3, 7, 11, 15, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);

    /* Process 4 floats (16 bytes) at a time */
    for (; i + 4 <= count; i += 4) {
        __m128i v = _mm_loadu_si128((const __m128i*)(src + i * 4));

        __m128i out0 = _mm_shuffle_epi8(v, s0);
        __m128i out1 = _mm_shuffle_epi8(v, s1);
        __m128i out2 = _mm_shuffle_epi8(v, s2);
        __m128i out3 = _mm_shuffle_epi8(v, s3);

        /* Store to transposed positions (use memcpy for unaligned access) */
        uint32_t t0 = (uint32_t)_mm_cvtsi128_si32(out0);
        uint32_t t1 = (uint32_t)_mm_cvtsi128_si32(out1);
        uint32_t t2 = (uint32_t)_mm_cvtsi128_si32(out2);
        uint32_t t3 = (uint32_t)_mm_cvtsi128_si32(out3);
        memcpy(output + 0 * count + i, &t0, sizeof(uint32_t));
        memcpy(output + 1 * count + i, &t1, sizeof(uint32_t));
        memcpy(output + 2 * count + i, &t2, sizeof(uint32_t));
        memcpy(output + 3 * count + i, &t3, sizeof(uint32_t));
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            output[b * count + i] = src[i * 4 + b];
        }
    }
}

/**
 * Decode byte stream split floats using SSE.
 */
void carquet_sse_byte_stream_split_decode_float(
    const uint8_t* data,
    int64_t count,
    float* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    /* Process 4 floats at a time */
    for (; i + 4 <= count; i += 4) {
        /* Load 4 bytes from each stream (use memcpy for unaligned access) */
        uint32_t b0, b1, b2, b3;
        memcpy(&b0, data + 0 * count + i, sizeof(uint32_t));
        memcpy(&b1, data + 1 * count + i, sizeof(uint32_t));
        memcpy(&b2, data + 2 * count + i, sizeof(uint32_t));
        memcpy(&b3, data + 3 * count + i, sizeof(uint32_t));

        __m128i v0 = _mm_cvtsi32_si128((int)b0);
        __m128i v1 = _mm_cvtsi32_si128((int)b1);
        __m128i v2 = _mm_cvtsi32_si128((int)b2);
        __m128i v3 = _mm_cvtsi32_si128((int)b3);

        /* Interleave bytes back into floats */
        __m128i lo01 = _mm_unpacklo_epi8(v0, v1);  /* a0b0 a1b1 a2b2 a3b3 ... */
        __m128i lo23 = _mm_unpacklo_epi8(v2, v3);  /* c0d0 c1d1 c2d2 c3d3 ... */
        __m128i result = _mm_unpacklo_epi16(lo01, lo23);  /* a0b0c0d0 a1b1c1d1 ... */

        _mm_storeu_si128((__m128i*)(dst + i * 4), result);
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            dst[i * 4 + b] = data[b * count + i];
        }
    }
}

/**
 * Encode doubles using byte stream split with SSE.
 */
void carquet_sse_byte_stream_split_encode_double(
    const double* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;
    const __m128i s0 = _mm_setr_epi8(0, 8, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s1 = _mm_setr_epi8(1, 9, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s2 = _mm_setr_epi8(2, 10, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s3 = _mm_setr_epi8(3, 11, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s4 = _mm_setr_epi8(4, 12, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s5 = _mm_setr_epi8(5, 13, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s6 = _mm_setr_epi8(6, 14, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i s7 = _mm_setr_epi8(7, 15, -1, -1, -1, -1, -1, -1,
                                     -1, -1, -1, -1, -1, -1, -1, -1);

    /* Process 2 doubles (16 bytes) at a time */
    for (; i + 2 <= count; i += 2) {
        __m128i v = _mm_loadu_si128((const __m128i*)(src + i * 8));

        uint16_t t0 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s0), 0);
        uint16_t t1 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s1), 0);
        uint16_t t2 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s2), 0);
        uint16_t t3 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s3), 0);
        uint16_t t4 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s4), 0);
        uint16_t t5 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s5), 0);
        uint16_t t6 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s6), 0);
        uint16_t t7 = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(v, s7), 0);

        memcpy(output + 0 * count + i, &t0, sizeof(t0));
        memcpy(output + 1 * count + i, &t1, sizeof(t1));
        memcpy(output + 2 * count + i, &t2, sizeof(t2));
        memcpy(output + 3 * count + i, &t3, sizeof(t3));
        memcpy(output + 4 * count + i, &t4, sizeof(t4));
        memcpy(output + 5 * count + i, &t5, sizeof(t5));
        memcpy(output + 6 * count + i, &t6, sizeof(t6));
        memcpy(output + 7 * count + i, &t7, sizeof(t7));
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            output[b * count + i] = src[i * 8 + b];
        }
    }
}

/**
 * Decode byte stream split doubles using SSE.
 */
void carquet_sse_byte_stream_split_decode_double(
    const uint8_t* data,
    int64_t count,
    double* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    for (; i + 2 <= count; i += 2) {
        uint16_t b0, b1, b2, b3, b4, b5, b6, b7;
        memcpy(&b0, data + 0 * count + i, sizeof(b0));
        memcpy(&b1, data + 1 * count + i, sizeof(b1));
        memcpy(&b2, data + 2 * count + i, sizeof(b2));
        memcpy(&b3, data + 3 * count + i, sizeof(b3));
        memcpy(&b4, data + 4 * count + i, sizeof(b4));
        memcpy(&b5, data + 5 * count + i, sizeof(b5));
        memcpy(&b6, data + 6 * count + i, sizeof(b6));
        memcpy(&b7, data + 7 * count + i, sizeof(b7));

        __m128i s0 = _mm_cvtsi32_si128((int)b0);
        __m128i s1 = _mm_cvtsi32_si128((int)b1);
        __m128i s2 = _mm_cvtsi32_si128((int)b2);
        __m128i s3 = _mm_cvtsi32_si128((int)b3);
        __m128i s4 = _mm_cvtsi32_si128((int)b4);
        __m128i s5 = _mm_cvtsi32_si128((int)b5);
        __m128i s6 = _mm_cvtsi32_si128((int)b6);
        __m128i s7 = _mm_cvtsi32_si128((int)b7);

        __m128i u01 = _mm_unpacklo_epi8(s0, s1);
        __m128i u23 = _mm_unpacklo_epi8(s2, s3);
        __m128i u45 = _mm_unpacklo_epi8(s4, s5);
        __m128i u67 = _mm_unpacklo_epi8(s6, s7);
        __m128i v0 = _mm_unpacklo_epi16(u01, u23);
        __m128i v1 = _mm_unpacklo_epi16(u45, u67);
        __m128i result = _mm_unpacklo_epi32(v0, v1);

        _mm_storeu_si128((__m128i*)(dst + i * 8), result);
    }

    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

/* ============================================================================
 * Delta Decoding - SSE Optimized (Prefix Sum)
 * ============================================================================
 */

/**
 * Apply prefix sum (cumulative sum) to int32 array using SSE.
 */
void carquet_sse_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB.
     * Delta encoding relies on modular arithmetic — _mm_add_epi32 is
     * already modular, so only the scalar accumulator needs fixing. */
    uint32_t sum = (uint32_t)initial;
    int64_t i = 0;

    /* SSE prefix sum for 4 elements at a time */
    for (; i + 4 <= count; i += 4) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));

        /* Partial prefix sums within the vector */
        /* v = [a, b, c, d] */
        /* After step 1: [a, a+b, c, c+d] */
        __m128i shifted1 = _mm_slli_si128(v, 4);
        v = _mm_add_epi32(v, shifted1);

        /* After step 2: [a, a+b, a+c, a+b+c+d] */
        __m128i shifted2 = _mm_slli_si128(v, 8);
        v = _mm_add_epi32(v, shifted2);

        /* Add running sum */
        __m128i sums = _mm_set1_epi32((int32_t)sum);
        v = _mm_add_epi32(v, sums);
        _mm_storeu_si128((__m128i*)(values + i), v);

        /* Update running sum to last element */
        sum = (uint32_t)_mm_extract_epi32(v, 3);
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint32_t)values[i];
        values[i] = (int32_t)sum;
    }
}

/**
 * Apply prefix sum to int64 array using SSE.
 */
void carquet_sse_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB. */
    uint64_t sum = (uint64_t)initial;
    int64_t i = 0;

    /* SSE prefix sum for 2 elements at a time */
    for (; i + 2 <= count; i += 2) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));

        /* v = [a, b] -> [a, a+b] */
        __m128i shifted = _mm_slli_si128(v, 8);
        v = _mm_add_epi64(v, shifted);

        /* Add running sum */
        __m128i sums = _mm_set1_epi64x((int64_t)sum);
        v = _mm_add_epi64(v, sums);
        _mm_storeu_si128((__m128i*)(values + i), v);

        /* Update running sum without spilling the whole vector */
        sum = (uint64_t)_mm_extract_epi64(v, 1);
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint64_t)values[i];
        values[i] = (int64_t)sum;
    }
}

/* ============================================================================
 * Dictionary Gather - SSE Optimized
 * ============================================================================
 */

/**
 * Gather int32 values from dictionary using indices (SSE).
 * Uses prefetching for better memory access patterns (matching NEON implementation).
 */
void carquet_sse_gather_i32(const int32_t* dict, const uint32_t* indices,
                             int64_t count, int32_t* output) {
    int64_t i = 0;

    /* Process 8 at a time with prefetching (like NEON) */
    for (; i + 8 <= count; i += 8) {
        /* Prefetch future indices */
        __builtin_prefetch(indices + i + 16, 0, 1);

        /* Prefetch dictionary entries for current batch */
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);
        __builtin_prefetch(dict + indices[i + 4], 0, 0);
        __builtin_prefetch(dict + indices[i + 6], 0, 0);

        /* First 4 values */
        int32_t v0 = cq_loadu(dict + (indices[i + 0]));
        int32_t v1 = cq_loadu(dict + (indices[i + 1]));
        int32_t v2 = cq_loadu(dict + (indices[i + 2]));
        int32_t v3 = cq_loadu(dict + (indices[i + 3]));
        __m128i result0 = _mm_set_epi32(v3, v2, v1, v0);
        _mm_storeu_si128((__m128i*)(output + i), result0);

        /* Second 4 values */
        int32_t v4 = cq_loadu(dict + (indices[i + 4]));
        int32_t v5 = cq_loadu(dict + (indices[i + 5]));
        int32_t v6 = cq_loadu(dict + (indices[i + 6]));
        int32_t v7 = cq_loadu(dict + (indices[i + 7]));
        __m128i result1 = _mm_set_epi32(v7, v6, v5, v4);
        _mm_storeu_si128((__m128i*)(output + i + 4), result1);
    }

    /* Process remaining 4 at a time */
    for (; i + 4 <= count; i += 4) {
        int32_t v0 = cq_loadu(dict + (indices[i + 0]));
        int32_t v1 = cq_loadu(dict + (indices[i + 1]));
        int32_t v2 = cq_loadu(dict + (indices[i + 2]));
        int32_t v3 = cq_loadu(dict + (indices[i + 3]));

        __m128i result = _mm_set_epi32(v3, v2, v1, v0);
        _mm_storeu_si128((__m128i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather float values from dictionary using indices (SSE).
 * Uses prefetching for better memory access patterns.
 */
void carquet_sse_gather_float(const float* dict, const uint32_t* indices,
                               int64_t count, float* output) {
    int64_t i = 0;

    /* Process 8 at a time with prefetching */
    for (; i + 8 <= count; i += 8) {
        /* Prefetch future indices */
        __builtin_prefetch(indices + i + 16, 0, 1);

        /* Prefetch dictionary entries */
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);
        __builtin_prefetch(dict + indices[i + 4], 0, 0);
        __builtin_prefetch(dict + indices[i + 6], 0, 0);

        /* First 4 values */
        float v0 = cq_loadu(dict + (indices[i + 0]));
        float v1 = cq_loadu(dict + (indices[i + 1]));
        float v2 = cq_loadu(dict + (indices[i + 2]));
        float v3 = cq_loadu(dict + (indices[i + 3]));
        __m128 result0 = _mm_set_ps(v3, v2, v1, v0);
        _mm_storeu_ps(output + i, result0);

        /* Second 4 values */
        float v4 = cq_loadu(dict + (indices[i + 4]));
        float v5 = cq_loadu(dict + (indices[i + 5]));
        float v6 = cq_loadu(dict + (indices[i + 6]));
        float v7 = cq_loadu(dict + (indices[i + 7]));
        __m128 result1 = _mm_set_ps(v7, v6, v5, v4);
        _mm_storeu_ps(output + i + 4, result1);
    }

    /* Process remaining 4 at a time */
    for (; i + 4 <= count; i += 4) {
        float v0 = cq_loadu(dict + (indices[i + 0]));
        float v1 = cq_loadu(dict + (indices[i + 1]));
        float v2 = cq_loadu(dict + (indices[i + 2]));
        float v3 = cq_loadu(dict + (indices[i + 3]));

        __m128 result = _mm_set_ps(v3, v2, v1, v0);
        _mm_storeu_ps(output + i, result);
    }

    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather int64 values from dictionary using indices (SSE).
 * Uses prefetching for better memory access patterns.
 */
void carquet_sse_gather_i64(const int64_t* dict, const uint32_t* indices,
                             int64_t count, int64_t* output) {
    int64_t i = 0;

    /* Process 4 at a time with prefetching */
    for (; i + 4 <= count; i += 4) {
        /* Prefetch future indices */
        __builtin_prefetch(indices + i + 8, 0, 1);

        /* Prefetch dictionary entries */
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);

        int64_t v0 = cq_loadu(dict + (indices[i + 0]));
        int64_t v1 = cq_loadu(dict + (indices[i + 1]));
        int64_t v2 = cq_loadu(dict + (indices[i + 2]));
        int64_t v3 = cq_loadu(dict + (indices[i + 3]));

        __m128i result0 = _mm_set_epi64x(v1, v0);
        __m128i result1 = _mm_set_epi64x(v3, v2);
        _mm_storeu_si128((__m128i*)(output + i), result0);
        _mm_storeu_si128((__m128i*)(output + i + 2), result1);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather double values from dictionary using indices (SSE).
 * Uses prefetching for better memory access patterns.
 */
void carquet_sse_gather_double(const double* dict, const uint32_t* indices,
                                int64_t count, double* output) {
    int64_t i = 0;

    /* Process 4 at a time with prefetching */
    for (; i + 4 <= count; i += 4) {
        /* Prefetch future indices */
        __builtin_prefetch(indices + i + 8, 0, 1);

        /* Prefetch dictionary entries */
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);

        double v0 = cq_loadu(dict + (indices[i + 0]));
        double v1 = cq_loadu(dict + (indices[i + 1]));
        double v2 = cq_loadu(dict + (indices[i + 2]));
        double v3 = cq_loadu(dict + (indices[i + 3]));

        __m128d result0 = _mm_set_pd(v1, v0);
        __m128d result1 = _mm_set_pd(v3, v2);
        _mm_storeu_pd(output + i, result0);
        _mm_storeu_pd(output + i + 2, result1);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

static inline int sse_indices_in_bounds_4(const uint32_t* indices, uint32_t limit) {
    __m128i idx = _mm_loadu_si128((const __m128i*)indices);
    __m128i bias = _mm_set1_epi32((int)0x80000000u);
    __m128i idx_biased = _mm_xor_si128(idx, bias);
    __m128i limit_biased = _mm_set1_epi32((int)(limit ^ 0x80000000u));
    __m128i cmp = _mm_cmplt_epi32(idx_biased, limit_biased);
    return _mm_movemask_epi8(cmp) == 0xFFFF;
}

bool carquet_sse_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                     const uint32_t* indices, int64_t count,
                                     int32_t* output) {
    int64_t i = 0;
    uint32_t limit = (uint32_t)dict_count;

    for (; i + 4 <= count; i += 4) {
        if (!sse_indices_in_bounds_4(indices + i, limit)) {
            return false;
        }

        __builtin_prefetch(indices + i + 8, 0, 1);
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);

        __m128i result = _mm_set_epi32(
            cq_loadu(dict + (indices[i + 3])),
            cq_loadu(dict + (indices[i + 2])),
            cq_loadu(dict + (indices[i + 1])),
            cq_loadu(dict + (indices[i + 0])));
        _mm_storeu_si128((__m128i*)(output + i), result);
    }

    for (; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= limit) {
            return false;
        }
        output[i] = cq_loadu(dict + (idx));
    }

    return true;
}

bool carquet_sse_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                     const uint32_t* indices, int64_t count,
                                     int64_t* output) {
    int64_t i = 0;
    uint32_t limit = (uint32_t)dict_count;

    for (; i + 4 <= count; i += 4) {
        if (!sse_indices_in_bounds_4(indices + i, limit)) {
            return false;
        }

        __builtin_prefetch(indices + i + 8, 0, 1);
        __builtin_prefetch(dict + indices[i], 0, 0);
        __builtin_prefetch(dict + indices[i + 2], 0, 0);

        __m128i result0 = _mm_set_epi64x(cq_loadu(dict + (indices[i + 1])), cq_loadu(dict + (indices[i + 0])));
        __m128i result1 = _mm_set_epi64x(cq_loadu(dict + (indices[i + 3])), cq_loadu(dict + (indices[i + 2])));
        _mm_storeu_si128((__m128i*)(output + i), result0);
        _mm_storeu_si128((__m128i*)(output + i + 2), result1);
    }

    for (; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= limit) {
            return false;
        }
        output[i] = cq_loadu(dict + (idx));
    }

    return true;
}

bool carquet_sse_checked_gather_float(const float* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       float* output) {
    return carquet_sse_checked_gather_i32(
        (const int32_t*)dict, dict_count, indices, count, (int32_t*)output);
}

bool carquet_sse_checked_gather_double(const double* dict, int32_t dict_count,
                                        const uint32_t* indices, int64_t count,
                                        double* output) {
    return carquet_sse_checked_gather_i64(
        (const int64_t*)dict, dict_count, indices, count, (int64_t*)output);
}

/* ============================================================================
 * Memcpy/Memset - SSE Optimized
 * ============================================================================
 */

/**
 * Fast memset for small-medium buffers using SSE.
 */
void carquet_sse_memset_small(void* dest, uint8_t value, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    __m128i v = _mm_set1_epi8((char)value);

    while (n >= 64) {
        _mm_storeu_si128((__m128i*)(d + 0), v);
        _mm_storeu_si128((__m128i*)(d + 16), v);
        _mm_storeu_si128((__m128i*)(d + 32), v);
        _mm_storeu_si128((__m128i*)(d + 48), v);
        d += 64;
        n -= 64;
    }

    while (n >= 16) {
        _mm_storeu_si128((__m128i*)d, v);
        d += 16;
        n -= 16;
    }

    while (n > 0) {
        *d++ = value;
        n--;
    }
}

/**
 * Fast memcpy for small-medium buffers using SSE.
 */
void carquet_sse_memcpy_small(void* dest, const void* src, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    const uint8_t* s = (const uint8_t*)src;

    while (n >= 64) {
        __m128i v0 = _mm_loadu_si128((const __m128i*)(s + 0));
        __m128i v1 = _mm_loadu_si128((const __m128i*)(s + 16));
        __m128i v2 = _mm_loadu_si128((const __m128i*)(s + 32));
        __m128i v3 = _mm_loadu_si128((const __m128i*)(s + 48));
        _mm_storeu_si128((__m128i*)(d + 0), v0);
        _mm_storeu_si128((__m128i*)(d + 16), v1);
        _mm_storeu_si128((__m128i*)(d + 32), v2);
        _mm_storeu_si128((__m128i*)(d + 48), v3);
        d += 64;
        s += 64;
        n -= 64;
    }

    while (n >= 16) {
        _mm_storeu_si128((__m128i*)d, _mm_loadu_si128((const __m128i*)s));
        d += 16;
        s += 16;
        n -= 16;
    }

    while (n > 0) {
        *d++ = *s++;
        n--;
    }
}

/* ============================================================================
 * Boolean Unpacking - SSE Optimized
 * ============================================================================
 */

/**
 * Unpack boolean values from packed bits to byte array.
 * Each output byte is 0 or 1.
 */
void carquet_sse_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;
    const __m128i mask = _mm_set_epi8(
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01,
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01
    );
    const __m128i shuf = _mm_setr_epi8(
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1
    );

    /* Process 16 bools (2 bytes) at a time */
    for (; i + 16 <= count; i += 16) {
        int byte_idx = (int)(i / 8);
        uint16_t packed;
        memcpy(&packed, input + byte_idx, 2);

        __m128i bits = _mm_set1_epi16(packed);

        /* Expand each byte for its corresponding bits */
        __m128i shuffled = _mm_shuffle_epi8(bits, shuf);

        /* AND with mask and normalize to 0/1 */
        __m128i masked = _mm_and_si128(shuffled, mask);
        __m128i result = _mm_min_epu8(masked, _mm_set1_epi8(1));

        _mm_storeu_si128((__m128i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        int byte_idx = (int)(i / 8);
        int bit_idx = (int)(i % 8);
        output[i] = (input[byte_idx] >> bit_idx) & 1;
    }
}

/**
 * Pack boolean values from byte array to packed bits.
 * Input bytes should be 0 or 1.
 */
void carquet_sse_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    /* Process 8 bools (1 output byte) at a time using the compare + movemask
     * trick: cmpeq(v, 0) is 0xFF for zero bytes, movemask gathers bit 7 of
     * each byte (byte k -> mask bit k), so ~mask has bit k set iff byte k is
     * non-zero — the scalar reference rule, for any non-zero byte. The earlier
     * `_mm_slli_epi32(bools, 7)` form was only correct for 0/1 inputs: a byte
     * such as 0x02 shifted its set bit into the neighbouring lane's bit 7 and
     * corrupted the *next* bool's bit (verified against the scalar reference
     * under Rosetta). loadl_epi64 zero-fills the upper 8 lanes; their mask
     * bits land above bit 7 and are dropped by the uint8_t cast. */
    const __m128i zero = _mm_setzero_si128();
    for (; i + 8 <= count; i += 8) {
        __m128i bools = _mm_loadl_epi64((const __m128i*)(input + i));
        int eq = _mm_movemask_epi8(_mm_cmpeq_epi8(bools, zero));
        output[i / 8] = (uint8_t)~eq;
    }

    /* Handle remaining */
    if (i < count) {
        uint8_t byte = 0;
        for (int64_t j = 0; j < count - i && j < 8; j++) {
            if (input[i + j]) {
                byte |= (1 << j);
            }
        }
        output[i / 8] = byte;
    }
}

/* ============================================================================
 * Compression Helpers
 * ============================================================================
 */

/**
 * Fast match copy for LZ4/Snappy decompression.
 * Handles overlapping copies correctly using SSE optimizations.
 */
void carquet_sse_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset) {
    if (offset >= 16) {
        /* Non-overlapping: use full SSE copies */
        while (len >= 16) {
            _mm_storeu_si128((__m128i*)dst, _mm_loadu_si128((const __m128i*)src));
            dst += 16;
            src += 16;
            len -= 16;
        }

        if (len >= 8) {
            _mm_storel_epi64((__m128i*)dst, _mm_loadl_epi64((const __m128i*)src));
            dst += 8;
            src += 8;
            len -= 8;
        }

        while (len > 0) {
            *dst++ = *src++;
            len--;
        }
    } else if (offset == 1) {
        /* Common pattern: fill with single byte */
        uint8_t val = *src;
        __m128i v = _mm_set1_epi8((char)val);

        while (len >= 16) {
            _mm_storeu_si128((__m128i*)dst, v);
            dst += 16;
            len -= 16;
        }

        while (len > 0) {
            *dst++ = val;
            len--;
        }
    } else if (offset == 2) {
        /* Fill with 2-byte pattern */
        uint8_t v0 = src[0], v1 = src[1];
        while (len >= 2) {
            *dst++ = v0;
            *dst++ = v1;
            len -= 2;
        }
        if (len) *dst = v0;
    } else if (offset == 4) {
        /* Fill with 4-byte pattern */
        uint32_t pattern;
        memcpy(&pattern, src, 4);
        __m128i v = _mm_set1_epi32((int32_t)pattern);

        while (len >= 16) {
            _mm_storeu_si128((__m128i*)dst, v);
            dst += 16;
            len -= 16;
        }

        while (len >= 4) {
            memcpy(dst, &pattern, 4);
            dst += 4;
            len -= 4;
        }

        for (size_t i = 0; i < len; i++) {
            dst[i] = src[i];
        }
    } else if (offset >= 8) {
        /* Offset 8-15: copy 8 bytes at a time (no SSE, but 64-bit safe) */
        while (len >= 8) {
            /* Must use memmove-style since src+8 may overlap dst */
            uint64_t v;
            memcpy(&v, src, 8);
            memcpy(dst, &v, 8);
            dst += 8;
            src += 8;
            len -= 8;
        }
        while (len > 0) {
            *dst++ = *src++;
            len--;
        }
    } else {
        /* Offset 3, 5, 6, 7: the period does not divide 16, so a tiled
         * 16-byte vector restarts the pattern at the wrong phase on every
         * store after the first. Copy byte by byte (src trails dst by
         * `offset`, so the overlap reproduces the pattern exactly). */
        while (len > 0) {
            *dst++ = *src++;
            len--;
        }
    }
}

/**
 * Count matching bytes between two buffers using SSE.
 * Returns the number of matching bytes from the start.
 */
size_t carquet_sse_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit) {
    const uint8_t* start = p;

    /* Fast path: compare 16 bytes at a time */
    while (p + 16 <= limit) {
        __m128i a = _mm_loadu_si128((const __m128i*)p);
        __m128i b = _mm_loadu_si128((const __m128i*)match);
        __m128i cmp = _mm_cmpeq_epi8(a, b);
        int mask = _mm_movemask_epi8(cmp);

        if (mask != 0xFFFF) {
            /* Find first differing byte */
            int first_diff = __builtin_ctz(~mask);
            return (size_t)(p - start) + (size_t)first_diff;
        }

        p += 16;
        match += 16;
    }

    /* Byte-by-byte for remaining */
    while (p < limit && *p == *match) {
        p++;
        match++;
    }

    return (size_t)(p - start);
}

/* ============================================================================
 * Definition Level Processing (Critical for Read Performance)
 * ============================================================================
 */

/**
 * Count non-null values using SIMD.
 * Counts how many def_levels[i] == max_def_level.
 */
int64_t carquet_sse_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    int64_t non_null_count = 0;
    int64_t i = 0;

    __m128i max_vec = _mm_set1_epi16(max_def_level);

    /* Process 8 int16_t values at a time */
    for (; i + 8 <= count; i += 8) {
        __m128i levels = _mm_loadu_si128((const __m128i*)(def_levels + i));
        __m128i cmp = _mm_cmpeq_epi16(levels, max_vec);
        int mask = _mm_movemask_epi8(cmp);
        /* Each matching int16 produces 2 bits set, so count and divide by 2 */
        non_null_count += __builtin_popcount(mask) >> 1;
    }

    /* Handle remaining */
    for (; i < count; i++) {
        if (def_levels[i] == max_def_level) {
            non_null_count++;
        }
    }

    return non_null_count;
}

/**
 * Build null bitmap from definition levels using SIMD.
 * Sets bit to 1 if def_levels[i] == max_def_level (present).
 */
void carquet_sse_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                    int16_t max_def_level, uint8_t* null_bitmap) {
    int64_t i = 0;

    __m128i max_vec = _mm_set1_epi16(max_def_level);
    __m128i zero = _mm_setzero_si128();

    /* Process 8 int16_t values -> 1 byte of bitmap */
    int64_t full_bytes = count / 8;
    for (int64_t b = 0; b < full_bytes; b++) {
        __m128i levels = _mm_loadu_si128((const __m128i*)(def_levels + b * 8));
        /* levels == max_def means present: result is 0x0000 or 0xFFFF per lane */
        __m128i cmp = _mm_cmpeq_epi16(levels, max_vec);
        /* Pack 8 int16 results (0x0000 or 0xFFFF) to 8 int8 (0x00 or 0xFF) */
        __m128i packed = _mm_packs_epi16(cmp, zero);
        /* movemask extracts bit 7 from each byte -> 8-bit result in low byte */
        null_bitmap[b] = (uint8_t)_mm_movemask_epi8(packed);
        i += 8;
    }

    /* Handle remaining bits */
    if (i < count) {
        uint8_t present_bits = 0;
        for (int64_t j = 0; i + j < count && j < 8; j++) {
            if (def_levels[i + j] == max_def_level) {
                present_bits |= (1 << j);
            }
        }
        null_bitmap[full_bytes] = present_bits;
    }
}

/**
 * Find the length of a run of identical int32 values starting from values[0].
 * Returns the number of consecutive values equal to values[0].
 * Uses SSE4.2 to compare 4 int32 values at a time.
 */
int64_t carquet_sse_find_run_length_i32(const int32_t* values, int64_t count) {
    if (count == 0) return 0;

    int32_t first = values[0];
    __m128i target = _mm_set1_epi32(first);
    int64_t i = 0;

    /* Check 4 at a time */
    for (; i + 4 <= count; i += 4) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));
        __m128i cmp = _mm_cmpeq_epi32(v, target);
        int mask = _mm_movemask_epi8(cmp);

        if (mask != 0xFFFF) {
            /* Not all equal - find first mismatch.
             * Each int32 produces 4 bits in the mask (all-1s if equal, all-0s if not).
             * Find first zero bit and divide by 4 to get the lane index. */
            int first_zero = __builtin_ctz(~mask);
            return i + (first_zero >> 2);
        }
    }

    /* Handle remaining */
    for (; i < count; i++) {
        if (values[i] != first) {
            return i;
        }
    }

    return count;
}

/**
 * Fill definition levels with a constant value using SIMD.
 */
void carquet_sse_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    int64_t i = 0;
    __m128i val_vec = _mm_set1_epi16(value);

    /* Process 8 int16_t values at a time */
    for (; i + 8 <= count; i += 8) {
        _mm_storeu_si128((__m128i*)(def_levels + i), val_vec);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        def_levels[i] = value;
    }
}

void carquet_sse_minmax_i32(const int32_t* values, int64_t count,
                             int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    __m128i min_vec = _mm_set1_epi32(min_v);
    __m128i max_vec = _mm_set1_epi32(max_v);
    int64_t i = 1;

    for (; i + 4 <= count; i += 4) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));
        min_vec = _mm_min_epi32(min_vec, v);
        max_vec = _mm_max_epi32(max_vec, v);
    }

    int32_t tmp_min[4];
    int32_t tmp_max[4];
    _mm_storeu_si128((__m128i*)tmp_min, min_vec);
    _mm_storeu_si128((__m128i*)tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_minmax_i64(const int64_t* values, int64_t count,
                             int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    __m128i min_vec = _mm_set1_epi64x(min_v);
    __m128i max_vec = _mm_set1_epi64x(max_v);
    int64_t i = 1;

    for (; i + 2 <= count; i += 2) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));
        __m128i lt = _mm_cmpgt_epi64(min_vec, v);
        __m128i gt = _mm_cmpgt_epi64(v, max_vec);
        min_vec = _mm_blendv_epi8(min_vec, v, lt);
        max_vec = _mm_blendv_epi8(max_vec, v, gt);
    }

    int64_t tmp_min[2];
    int64_t tmp_max[2];
    _mm_storeu_si128((__m128i*)tmp_min, min_vec);
    _mm_storeu_si128((__m128i*)tmp_max, max_vec);
    for (int j = 0; j < 2; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_minmax_float(const float* values, int64_t count,
                               float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m128 min_vec = _mm_set1_ps(min_v);
    __m128 max_vec = _mm_set1_ps(max_v);
    int64_t i = 1;

    for (; i + 4 <= count; i += 4) {
        __m128 v = _mm_loadu_ps(values + i);
        __m128 lt = _mm_cmplt_ps(v, min_vec);
        __m128 gt = _mm_cmpgt_ps(v, max_vec);
        min_vec = _mm_blendv_ps(min_vec, v, lt);
        max_vec = _mm_blendv_ps(max_vec, v, gt);
    }

    float tmp_min[4];
    float tmp_max[4];
    _mm_storeu_ps(tmp_min, min_vec);
    _mm_storeu_ps(tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_minmax_double(const double* values, int64_t count,
                                double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m128d min_vec = _mm_set1_pd(min_v);
    __m128d max_vec = _mm_set1_pd(max_v);
    int64_t i = 1;

    for (; i + 2 <= count; i += 2) {
        __m128d v = _mm_loadu_pd(values + i);
        __m128d lt = _mm_cmplt_pd(v, min_vec);
        __m128d gt = _mm_cmpgt_pd(v, max_vec);
        min_vec = _mm_blendv_pd(min_vec, v, lt);
        max_vec = _mm_blendv_pd(max_vec, v, gt);
    }

    double tmp_min[2];
    double tmp_max[2];
    _mm_storeu_pd(tmp_min, min_vec);
    _mm_storeu_pd(tmp_max, max_vec);
    for (int j = 0; j < 2; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                  int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    __m128i min_vec = _mm_set1_epi32(min_v);
    __m128i max_vec = _mm_set1_epi32(max_v);
    int64_t i = 0;

    for (; i + 4 <= count; i += 4) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));
        _mm_storeu_si128((__m128i*)(output + i), v);
        min_vec = _mm_min_epi32(min_vec, v);
        max_vec = _mm_max_epi32(max_vec, v);
    }

    int32_t tmp_min[4];
    int32_t tmp_max[4];
    _mm_storeu_si128((__m128i*)tmp_min, min_vec);
    _mm_storeu_si128((__m128i*)tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        int32_t v = values[i];
        output[i] = v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                  int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    __m128i min_vec = _mm_set1_epi64x(min_v);
    __m128i max_vec = _mm_set1_epi64x(max_v);
    int64_t i = 0;

    for (; i + 2 <= count; i += 2) {
        __m128i v = _mm_loadu_si128((const __m128i*)(values + i));
        _mm_storeu_si128((__m128i*)(output + i), v);
        __m128i lt = _mm_cmpgt_epi64(min_vec, v);
        __m128i gt = _mm_cmpgt_epi64(v, max_vec);
        min_vec = _mm_blendv_epi8(min_vec, v, lt);
        max_vec = _mm_blendv_epi8(max_vec, v, gt);
    }

    int64_t tmp_min[2];
    int64_t tmp_max[2];
    _mm_storeu_si128((__m128i*)tmp_min, min_vec);
    _mm_storeu_si128((__m128i*)tmp_max, max_vec);
    for (int j = 0; j < 2; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        int64_t v = values[i];
        output[i] = v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_copy_minmax_float(const float* values, int64_t count, float* output,
                                    float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m128 min_vec = _mm_set1_ps(min_v);
    __m128 max_vec = _mm_set1_ps(max_v);
    int64_t i = 0;

    for (; i + 4 <= count; i += 4) {
        __m128 v = _mm_loadu_ps(values + i);
        _mm_storeu_ps(output + i, v);
        __m128 lt = _mm_cmplt_ps(v, min_vec);
        __m128 gt = _mm_cmpgt_ps(v, max_vec);
        min_vec = _mm_blendv_ps(min_vec, v, lt);
        max_vec = _mm_blendv_ps(max_vec, v, gt);
    }

    float tmp_min[4];
    float tmp_max[4];
    _mm_storeu_ps(tmp_min, min_vec);
    _mm_storeu_ps(tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        float v = values[i];
        output[i] = v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sse_copy_minmax_double(const double* values, int64_t count, double* output,
                                     double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m128d min_vec = _mm_set1_pd(min_v);
    __m128d max_vec = _mm_set1_pd(max_v);
    int64_t i = 0;

    for (; i + 2 <= count; i += 2) {
        __m128d v = _mm_loadu_pd(values + i);
        _mm_storeu_pd(output + i, v);
        __m128d lt = _mm_cmplt_pd(v, min_vec);
        __m128d gt = _mm_cmpgt_pd(v, max_vec);
        min_vec = _mm_blendv_pd(min_vec, v, lt);
        max_vec = _mm_blendv_pd(max_vec, v, gt);
    }

    double tmp_min[2];
    double tmp_max[2];
    _mm_storeu_pd(tmp_min, min_vec);
    _mm_storeu_pd(tmp_max, max_vec);
    for (int j = 0; j < 2; j++) {
        if (tmp_min[j] < min_v) min_v = tmp_min[j];
        if (tmp_max[j] > max_v) max_v = tmp_max[j];
    }
    for (; i < count; i++) {
        double v = values[i];
        output[i] = v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    *min_value = min_v;
    *max_value = max_v;
}

#endif /* __SSE4_2__ */
#endif /* x86 */
