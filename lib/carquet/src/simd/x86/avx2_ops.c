/**
 * @file avx2_ops.c
 * @brief AVX2 optimized operations for x86-64 processors
 *
 * Provides SIMD-accelerated implementations using 256-bit vectors:
 * - Bit unpacking for common bit widths
 * - Byte stream split/merge (for BYTE_STREAM_SPLIT encoding)
 * - Delta decoding (prefix sums)
 * - Dictionary gather operations (using AVX2 gather instructions)
 * - Boolean packing/unpacking
 */

#include <carquet/error.h>
#include "simd/simd_unaligned.h"
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
/* Check for AVX2 support - MSVC defines __AVX2__ when /arch:AVX2 is used */
#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))

#ifdef _MSC_VER
#include <intrin.h>

static inline int msvc_ctz(unsigned int x) {
    unsigned long index;
    _BitScanForward(&index, x);
    return (int)index;
}
#define __builtin_ctz(x) msvc_ctz(x)
#define __builtin_popcount(x) __popcnt(x)
#endif
#include <immintrin.h>

static inline uint16_t avx2_read_le16(const uint8_t* p) {
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}

static inline uint32_t avx2_read_le24(const uint8_t* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16);
}

static inline uint64_t avx2_read_le40(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32);
}

static inline uint64_t avx2_read_le48(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40);
}

static inline uint64_t avx2_read_le56(const uint8_t* p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48);
}

/* ============================================================================
 * Bit Unpacking - AVX2 Optimized
 * ============================================================================
 */

/**
 * Unpack 8 1-bit values using AVX2.
 */
void carquet_avx2_bitunpack8_1bit(const uint8_t* input, uint32_t* values) {
    __m128i bytes = _mm_set1_epi8((char)input[0]);
    const __m128i bit_mask = _mm_setr_epi8(
        0x01, 0x02, 0x04, 0x08,
        0x10, 0x20, 0x40, (char)0x80,
        0, 0, 0, 0, 0, 0, 0, 0
    );
    __m128i masked = _mm_and_si128(bytes, bit_mask);
    __m128i cmp = _mm_cmpeq_epi8(masked, bit_mask);
    __m128i result8 = _mm_and_si128(cmp, _mm_set1_epi8(1));
    __m256i result = _mm256_cvtepu8_epi32(result8);
    _mm256_storeu_si256((__m256i*)values, result);
}

void carquet_avx2_bitunpack8_2bit(const uint8_t* input, uint32_t* values) {
    uint16_t v = avx2_read_le16(input);
    __m256i result = _mm256_setr_epi32(
        (int)((v >> 0) & 0x3), (int)((v >> 2) & 0x3),
        (int)((v >> 4) & 0x3), (int)((v >> 6) & 0x3),
        (int)((v >> 8) & 0x3), (int)((v >> 10) & 0x3),
        (int)((v >> 12) & 0x3), (int)((v >> 14) & 0x3));
    _mm256_storeu_si256((__m256i*)values, result);
}

void carquet_avx2_bitunpack8_3bit(const uint8_t* input, uint32_t* values) {
    uint32_t v = avx2_read_le24(input);
    __m256i result = _mm256_setr_epi32(
        (int)((v >> 0) & 0x7), (int)((v >> 3) & 0x7),
        (int)((v >> 6) & 0x7), (int)((v >> 9) & 0x7),
        (int)((v >> 12) & 0x7), (int)((v >> 15) & 0x7),
        (int)((v >> 18) & 0x7), (int)((v >> 21) & 0x7));
    _mm256_storeu_si256((__m256i*)values, result);
}


/**
 * Unpack 8 4-bit values using AVX2.
 */
void carquet_avx2_bitunpack8_4bit(const uint8_t* input, uint32_t* values) {
    __m128i bytes = _mm_cvtsi32_si128(*(const int32_t*)input);
    __m128i lo_nibbles = _mm_and_si128(bytes, _mm_set1_epi8(0x0F));
    __m128i hi_nibbles = _mm_and_si128(_mm_srli_epi16(bytes, 4), _mm_set1_epi8(0x0F));
    __m128i interleaved = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
    __m256i result = _mm256_cvtepu8_epi32(interleaved);
    _mm256_storeu_si256((__m256i*)values, result);
}

void carquet_avx2_bitunpack8_5bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = avx2_read_le40(input);
    __m256i result = _mm256_setr_epi32(
        (int)((v >> 0) & 0x1F), (int)((v >> 5) & 0x1F),
        (int)((v >> 10) & 0x1F), (int)((v >> 15) & 0x1F),
        (int)((v >> 20) & 0x1F), (int)((v >> 25) & 0x1F),
        (int)((v >> 30) & 0x1F), (int)((v >> 35) & 0x1F));
    _mm256_storeu_si256((__m256i*)values, result);
}

void carquet_avx2_bitunpack8_6bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = avx2_read_le48(input);
    __m256i result = _mm256_setr_epi32(
        (int)((v >> 0) & 0x3F), (int)((v >> 6) & 0x3F),
        (int)((v >> 12) & 0x3F), (int)((v >> 18) & 0x3F),
        (int)((v >> 24) & 0x3F), (int)((v >> 30) & 0x3F),
        (int)((v >> 36) & 0x3F), (int)((v >> 42) & 0x3F));
    _mm256_storeu_si256((__m256i*)values, result);
}

void carquet_avx2_bitunpack8_7bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = avx2_read_le56(input);
    __m256i result = _mm256_setr_epi32(
        (int)((v >> 0) & 0x7F), (int)((v >> 7) & 0x7F),
        (int)((v >> 14) & 0x7F), (int)((v >> 21) & 0x7F),
        (int)((v >> 28) & 0x7F), (int)((v >> 35) & 0x7F),
        (int)((v >> 42) & 0x7F), (int)((v >> 49) & 0x7F));
    _mm256_storeu_si256((__m256i*)values, result);
}

/**
 * Unpack 16 4-bit values using AVX2.
 */
void carquet_avx2_bitunpack16_4bit(const uint8_t* input, uint32_t* values) {
    /* Load 8 bytes containing 16 x 4-bit values */
    __m128i bytes = _mm_loadl_epi64((const __m128i*)input);

    /* Split nibbles */
    __m128i lo_nibbles = _mm_and_si128(bytes, _mm_set1_epi8(0x0F));
    __m128i hi_nibbles = _mm_srli_epi16(bytes, 4);
    hi_nibbles = _mm_and_si128(hi_nibbles, _mm_set1_epi8(0x0F));

    /* Interleave */
    __m128i interleaved = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);

    /* Expand to 32-bit using AVX2 */
    __m256i result = _mm256_cvtepu8_epi32(interleaved);
    _mm256_storeu_si256((__m256i*)values, result);

    /* Process second half */
    __m128i second_half = _mm_unpackhi_epi64(interleaved, interleaved);
    result = _mm256_cvtepu8_epi32(second_half);
    _mm256_storeu_si256((__m256i*)(values + 8), result);
}

/**
 * Unpack 8 8-bit values using AVX2.
 */
void carquet_avx2_bitunpack8_8bit(const uint8_t* input, uint32_t* values) {
    __m128i bytes = _mm_loadl_epi64((const __m128i*)input);
    __m256i result = _mm256_cvtepu8_epi32(bytes);
    _mm256_storeu_si256((__m256i*)values, result);
}

/**
 * Unpack 16 8-bit values using AVX2 (widen u8 to u32).
 */
void carquet_avx2_bitunpack16_8bit(const uint8_t* input, uint32_t* values) {
    /* Load 16 bytes */
    __m128i bytes = _mm_loadu_si128((const __m128i*)input);

    /* Expand low 8 bytes to 8 x 32-bit */
    __m256i lo = _mm256_cvtepu8_epi32(bytes);
    _mm256_storeu_si256((__m256i*)values, lo);

    /* Expand high 8 bytes to 8 x 32-bit */
    __m128i hi_bytes = _mm_srli_si128(bytes, 8);
    __m256i hi = _mm256_cvtepu8_epi32(hi_bytes);
    _mm256_storeu_si256((__m256i*)(values + 8), hi);
}

/**
 * Unpack 8 16-bit values to 32-bit using AVX2.
 */
void carquet_avx2_bitunpack8_16bit(const uint8_t* input, uint32_t* values) {
    __m128i words = _mm_loadu_si128((const __m128i*)input);
    __m256i result = _mm256_cvtepu16_epi32(words);
    _mm256_storeu_si256((__m256i*)values, result);
}

/* ============================================================================
 * Byte Stream Split - AVX2 Optimized
 * ============================================================================
 */

/**
 * Encode floats using byte stream split with AVX2.
 * Processes 8 floats (32 bytes) at a time.
 */
void carquet_avx2_byte_stream_split_encode_float(
    const float* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;
    const __m256i s0 = _mm256_setr_epi8(
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s1 = _mm256_setr_epi8(
        1, 5, 9, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        1, 5, 9, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s2 = _mm256_setr_epi8(
        2, 6, 10, 14, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        2, 6, 10, 14, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s3 = _mm256_setr_epi8(
        3, 7, 11, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        3, 7, 11, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);

    /* Process 8 floats (32 bytes) at a time */
    for (; i + 8 <= count; i += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(src + i * 4));

        __m256i out0 = _mm256_shuffle_epi8(v, s0);
        __m256i out1 = _mm256_shuffle_epi8(v, s1);
        __m256i out2 = _mm256_shuffle_epi8(v, s2);
        __m256i out3 = _mm256_shuffle_epi8(v, s3);

        /* Extract and combine low and high 128-bit lanes */
        uint32_t b0_lo = _mm256_extract_epi32(out0, 0);
        uint32_t b0_hi = _mm256_extract_epi32(out0, 4);
        uint32_t b1_lo = _mm256_extract_epi32(out1, 0);
        uint32_t b1_hi = _mm256_extract_epi32(out1, 4);
        uint32_t b2_lo = _mm256_extract_epi32(out2, 0);
        uint32_t b2_hi = _mm256_extract_epi32(out2, 4);
        uint32_t b3_lo = _mm256_extract_epi32(out3, 0);
        uint32_t b3_hi = _mm256_extract_epi32(out3, 4);

        /* Store to transposed positions (use memcpy for unaligned access) */
        memcpy(output + 0 * count + i, &b0_lo, sizeof(uint32_t));
        memcpy(output + 0 * count + i + 4, &b0_hi, sizeof(uint32_t));
        memcpy(output + 1 * count + i, &b1_lo, sizeof(uint32_t));
        memcpy(output + 1 * count + i + 4, &b1_hi, sizeof(uint32_t));
        memcpy(output + 2 * count + i, &b2_lo, sizeof(uint32_t));
        memcpy(output + 2 * count + i + 4, &b2_hi, sizeof(uint32_t));
        memcpy(output + 3 * count + i, &b3_lo, sizeof(uint32_t));
        memcpy(output + 3 * count + i + 4, &b3_hi, sizeof(uint32_t));
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            output[b * count + i] = src[i * 4 + b];
        }
    }
}

/**
 * Decode byte stream split floats using full-width AVX2.
 * Processes 16 floats per iteration using 256-bit loads and unpack cascade.
 */
void carquet_avx2_byte_stream_split_decode_float(
    const uint8_t* data,
    int64_t count,
    float* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    /* Process 16 floats at a time with full 256-bit AVX2 operations.
     * Load 16 bytes from each of 4 streams into __m256i (via two 8-byte halves),
     * then use 256-bit unpack cascade + lane-fix permute. */
    for (; i + 16 <= count; i += 16) {
        __m128i s0_lo = _mm_loadl_epi64((const __m128i*)(data + 0 * count + i));
        __m128i s0_hi = _mm_loadl_epi64((const __m128i*)(data + 0 * count + i + 8));
        __m128i s1_lo = _mm_loadl_epi64((const __m128i*)(data + 1 * count + i));
        __m128i s1_hi = _mm_loadl_epi64((const __m128i*)(data + 1 * count + i + 8));
        __m128i s2_lo = _mm_loadl_epi64((const __m128i*)(data + 2 * count + i));
        __m128i s2_hi = _mm_loadl_epi64((const __m128i*)(data + 2 * count + i + 8));
        __m128i s3_lo = _mm_loadl_epi64((const __m128i*)(data + 3 * count + i));
        __m128i s3_hi = _mm_loadl_epi64((const __m128i*)(data + 3 * count + i + 8));

        __m256i b0 = _mm256_inserti128_si256(_mm256_castsi128_si256(s0_lo), s0_hi, 1);
        __m256i b1 = _mm256_inserti128_si256(_mm256_castsi128_si256(s1_lo), s1_hi, 1);
        __m256i b2 = _mm256_inserti128_si256(_mm256_castsi128_si256(s2_lo), s2_hi, 1);
        __m256i b3 = _mm256_inserti128_si256(_mm256_castsi128_si256(s3_lo), s3_hi, 1);

        /* AVX2 unpack operates per-lane */
        __m256i lo01 = _mm256_unpacklo_epi8(b0, b1);
        __m256i lo23 = _mm256_unpacklo_epi8(b2, b3);

        __m256i r0 = _mm256_unpacklo_epi16(lo01, lo23);
        __m256i r1 = _mm256_unpackhi_epi16(lo01, lo23);

        /* Fix lane ordering for sequential output */
        __m256i out0 = _mm256_permute2x128_si256(r0, r1, 0x20);
        __m256i out1 = _mm256_permute2x128_si256(r0, r1, 0x31);

        _mm256_storeu_si256((__m256i*)(dst + i * 4), out0);
        _mm256_storeu_si256((__m256i*)(dst + i * 4 + 32), out1);
    }

    /* 8-float fallback using 128-bit ops */
    for (; i + 8 <= count; i += 8) {
        uint64_t t0, t1, t2, t3;
        memcpy(&t0, data + 0 * count + i, sizeof(uint64_t));
        memcpy(&t1, data + 1 * count + i, sizeof(uint64_t));
        memcpy(&t2, data + 2 * count + i, sizeof(uint64_t));
        memcpy(&t3, data + 3 * count + i, sizeof(uint64_t));
        __m128i b0 = _mm_cvtsi64_si128((long long)t0);
        __m128i b1 = _mm_cvtsi64_si128((long long)t1);
        __m128i b2 = _mm_cvtsi64_si128((long long)t2);
        __m128i b3 = _mm_cvtsi64_si128((long long)t3);

        __m128i lo01 = _mm_unpacklo_epi8(b0, b1);
        __m128i lo23 = _mm_unpacklo_epi8(b2, b3);
        __m128i result_lo = _mm_unpacklo_epi16(lo01, lo23);
        __m128i result_hi = _mm_unpackhi_epi16(lo01, lo23);

        _mm_storeu_si128((__m128i*)(dst + i * 4), result_lo);
        _mm_storeu_si128((__m128i*)(dst + i * 4 + 16), result_hi);
    }

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            dst[i * 4 + b] = data[b * count + i];
        }
    }
}

/**
 * Encode doubles using byte stream split with AVX2.
 */
void carquet_avx2_byte_stream_split_encode_double(
    const double* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;
    const __m256i s0 = _mm256_setr_epi8(
        0, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s1 = _mm256_setr_epi8(
        1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s2 = _mm256_setr_epi8(
        2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s3 = _mm256_setr_epi8(
        3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s4 = _mm256_setr_epi8(
        4, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        4, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s5 = _mm256_setr_epi8(
        5, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        5, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s6 = _mm256_setr_epi8(
        6, 14, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        6, 14, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i s7 = _mm256_setr_epi8(
        7, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        7, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);

    /* Process 4 doubles (32 bytes) at a time */
    for (; i + 4 <= count; i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(src + i * 8));
        __m256i out0 = _mm256_shuffle_epi8(v, s0);
        __m256i out1 = _mm256_shuffle_epi8(v, s1);
        __m256i out2 = _mm256_shuffle_epi8(v, s2);
        __m256i out3 = _mm256_shuffle_epi8(v, s3);
        __m256i out4 = _mm256_shuffle_epi8(v, s4);
        __m256i out5 = _mm256_shuffle_epi8(v, s5);
        __m256i out6 = _mm256_shuffle_epi8(v, s6);
        __m256i out7 = _mm256_shuffle_epi8(v, s7);

        __m128i lo0 = _mm256_castsi256_si128(out0);
        __m128i lo1 = _mm256_castsi256_si128(out1);
        __m128i lo2 = _mm256_castsi256_si128(out2);
        __m128i lo3 = _mm256_castsi256_si128(out3);
        __m128i lo4 = _mm256_castsi256_si128(out4);
        __m128i lo5 = _mm256_castsi256_si128(out5);
        __m128i lo6 = _mm256_castsi256_si128(out6);
        __m128i lo7 = _mm256_castsi256_si128(out7);
        __m128i hi0 = _mm256_extracti128_si256(out0, 1);
        __m128i hi1 = _mm256_extracti128_si256(out1, 1);
        __m128i hi2 = _mm256_extracti128_si256(out2, 1);
        __m128i hi3 = _mm256_extracti128_si256(out3, 1);
        __m128i hi4 = _mm256_extracti128_si256(out4, 1);
        __m128i hi5 = _mm256_extracti128_si256(out5, 1);
        __m128i hi6 = _mm256_extracti128_si256(out6, 1);
        __m128i hi7 = _mm256_extracti128_si256(out7, 1);

        uint32_t t0 = (uint32_t)_mm_extract_epi16(lo0, 0) |
                      ((uint32_t)_mm_extract_epi16(hi0, 0) << 16);
        uint32_t t1 = (uint32_t)_mm_extract_epi16(lo1, 0) |
                      ((uint32_t)_mm_extract_epi16(hi1, 0) << 16);
        uint32_t t2 = (uint32_t)_mm_extract_epi16(lo2, 0) |
                      ((uint32_t)_mm_extract_epi16(hi2, 0) << 16);
        uint32_t t3 = (uint32_t)_mm_extract_epi16(lo3, 0) |
                      ((uint32_t)_mm_extract_epi16(hi3, 0) << 16);
        uint32_t t4 = (uint32_t)_mm_extract_epi16(lo4, 0) |
                      ((uint32_t)_mm_extract_epi16(hi4, 0) << 16);
        uint32_t t5 = (uint32_t)_mm_extract_epi16(lo5, 0) |
                      ((uint32_t)_mm_extract_epi16(hi5, 0) << 16);
        uint32_t t6 = (uint32_t)_mm_extract_epi16(lo6, 0) |
                      ((uint32_t)_mm_extract_epi16(hi6, 0) << 16);
        uint32_t t7 = (uint32_t)_mm_extract_epi16(lo7, 0) |
                      ((uint32_t)_mm_extract_epi16(hi7, 0) << 16);

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
 * Decode byte stream split doubles using full-width AVX2.
 * Processes 32 doubles per iteration using 256-bit loads, a 3-stage unpack
 * cascade, and lane-fix permutations — 8x throughput vs the old 4-double path.
 */
void carquet_avx2_byte_stream_split_decode_double(
    const uint8_t* data,
    int64_t count,
    double* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    /* Process 32 doubles at a time with full-width AVX2.
     * Load 32 bytes from each of 8 streams, do 256-bit unpack cascade,
     * then fix lane ordering with permute2x128. */
    for (; i + 32 <= count; i += 32) {
        __m256i s0 = _mm256_loadu_si256((const __m256i*)(data + 0 * count + i));
        __m256i s1 = _mm256_loadu_si256((const __m256i*)(data + 1 * count + i));
        __m256i s2 = _mm256_loadu_si256((const __m256i*)(data + 2 * count + i));
        __m256i s3 = _mm256_loadu_si256((const __m256i*)(data + 3 * count + i));
        __m256i s4 = _mm256_loadu_si256((const __m256i*)(data + 4 * count + i));
        __m256i s5 = _mm256_loadu_si256((const __m256i*)(data + 5 * count + i));
        __m256i s6 = _mm256_loadu_si256((const __m256i*)(data + 6 * count + i));
        __m256i s7 = _mm256_loadu_si256((const __m256i*)(data + 7 * count + i));

        /* Stage 1: byte-level interleave of stream pairs (per 128-bit lane) */
        __m256i a01_lo = _mm256_unpacklo_epi8(s0, s1);
        __m256i a01_hi = _mm256_unpackhi_epi8(s0, s1);
        __m256i a23_lo = _mm256_unpacklo_epi8(s2, s3);
        __m256i a23_hi = _mm256_unpackhi_epi8(s2, s3);
        __m256i a45_lo = _mm256_unpacklo_epi8(s4, s5);
        __m256i a45_hi = _mm256_unpackhi_epi8(s4, s5);
        __m256i a67_lo = _mm256_unpacklo_epi8(s6, s7);
        __m256i a67_hi = _mm256_unpackhi_epi8(s6, s7);

        /* Stage 2: 16-bit interleave of quad groups */
        __m256i b0 = _mm256_unpacklo_epi16(a01_lo, a23_lo);
        __m256i b1 = _mm256_unpackhi_epi16(a01_lo, a23_lo);
        __m256i b2 = _mm256_unpacklo_epi16(a01_hi, a23_hi);
        __m256i b3 = _mm256_unpackhi_epi16(a01_hi, a23_hi);
        __m256i b4 = _mm256_unpacklo_epi16(a45_lo, a67_lo);
        __m256i b5 = _mm256_unpackhi_epi16(a45_lo, a67_lo);
        __m256i b6 = _mm256_unpacklo_epi16(a45_hi, a67_hi);
        __m256i b7 = _mm256_unpackhi_epi16(a45_hi, a67_hi);

        /* Stage 3: 32-bit interleave assembles full 8-byte doubles */
        __m256i c0 = _mm256_unpacklo_epi32(b0, b4);
        __m256i c1 = _mm256_unpackhi_epi32(b0, b4);
        __m256i c2 = _mm256_unpacklo_epi32(b1, b5);
        __m256i c3 = _mm256_unpackhi_epi32(b1, b5);
        __m256i c4 = _mm256_unpacklo_epi32(b2, b6);
        __m256i c5 = _mm256_unpackhi_epi32(b2, b6);
        __m256i c6 = _mm256_unpacklo_epi32(b3, b7);
        __m256i c7 = _mm256_unpackhi_epi32(b3, b7);

        /* Fix lane ordering and store 32 doubles sequentially */
        _mm256_storeu_si256((__m256i*)(dst + i * 8),       _mm256_permute2x128_si256(c0, c1, 0x20));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 32),  _mm256_permute2x128_si256(c2, c3, 0x20));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 64),  _mm256_permute2x128_si256(c4, c5, 0x20));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 96),  _mm256_permute2x128_si256(c6, c7, 0x20));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 128), _mm256_permute2x128_si256(c0, c1, 0x31));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 160), _mm256_permute2x128_si256(c2, c3, 0x31));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 192), _mm256_permute2x128_si256(c4, c5, 0x31));
        _mm256_storeu_si256((__m256i*)(dst + i * 8 + 224), _mm256_permute2x128_si256(c6, c7, 0x31));
    }

    /* 4-double fallback using 128-bit ops */
    for (; i + 4 <= count; i += 4) {
        uint32_t b0, b1, b2, b3, b4, b5, b6, b7;
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
        __m128i lo_ab = _mm_unpacklo_epi32(v0, v1);
        __m128i hi_cd = _mm_unpackhi_epi32(v0, v1);

        _mm_storeu_si128((__m128i*)(dst + i * 8), lo_ab);
        _mm_storeu_si128((__m128i*)(dst + i * 8 + 16), hi_cd);
    }

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

/* ============================================================================
 * Delta Decoding - AVX2 Optimized (Prefix Sum)
 * ============================================================================
 */

/**
 * Apply prefix sum (cumulative sum) to int32 array using AVX2.
 */
void carquet_avx2_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB.
     * Delta encoding relies on modular arithmetic — _mm256_add_epi32 is
     * already modular, so only the scalar accumulator needs fixing. */
    uint32_t sum = (uint32_t)initial;
    int64_t i = 0;

    /* AVX2 prefix sum for 8 elements at a time */
    for (; i + 8 <= count; i += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));

        /* Partial prefix sums within the vector */
        /* Step 1: Add adjacent pairs */
        __m256i shifted1 = _mm256_slli_si256(v, 4);
        v = _mm256_add_epi32(v, shifted1);

        /* Step 2: Add pairs that are 2 apart */
        __m256i shifted2 = _mm256_slli_si256(v, 8);
        v = _mm256_add_epi32(v, shifted2);

        /* Step 3: Handle cross-lane (bit tricky with AVX2) */
        /* Extract lane 0's last value and add to all of lane 1 */
        __m128i lo = _mm256_extracti128_si256(v, 0);
        __m128i hi = _mm256_extracti128_si256(v, 1);

        int32_t lane0_sum = _mm_extract_epi32(lo, 3);
        __m128i lane0_broadcast = _mm_set1_epi32(lane0_sum);
        hi = _mm_add_epi32(hi, lane0_broadcast);

        v = _mm256_inserti128_si256(v, hi, 1);

        /* Add running sum */
        __m256i sums = _mm256_set1_epi32((int32_t)sum);
        v = _mm256_add_epi32(v, sums);
        _mm256_storeu_si256((__m256i*)(values + i), v);

        /* Update running sum to last element */
        sum = (uint32_t)_mm256_extract_epi32(v, 7);
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint32_t)values[i];
        values[i] = (int32_t)sum;
    }
}

/**
 * Apply prefix sum to int64 array using AVX2.
 */
void carquet_avx2_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB. */
    uint64_t sum = (uint64_t)initial;
    int64_t i = 0;

    /* AVX2 prefix sum for 4 elements at a time */
    for (; i + 4 <= count; i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));

        /* Partial prefix sums */
        __m256i shifted1 = _mm256_slli_si256(v, 8);
        v = _mm256_add_epi64(v, shifted1);

        /* Cross-lane fixup */
        __m128i lo = _mm256_extracti128_si256(v, 0);
        __m128i hi = _mm256_extracti128_si256(v, 1);

        int64_t lane0_last;
        _mm_storel_epi64((__m128i*)&lane0_last, _mm_srli_si128(lo, 8));
        __m128i lane0_broadcast = _mm_set1_epi64x(lane0_last);
        hi = _mm_add_epi64(hi, lane0_broadcast);

        v = _mm256_inserti128_si256(v, hi, 1);

        /* Add running sum */
        __m256i sums = _mm256_set1_epi64x((int64_t)sum);
        v = _mm256_add_epi64(v, sums);
        _mm256_storeu_si256((__m256i*)(values + i), v);

        /* Update running sum */
        sum = (uint64_t)_mm256_extract_epi64(v, 3);
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint64_t)values[i];
        values[i] = (int64_t)sum;
    }
}

/* ============================================================================
 * Dictionary Gather - AVX2 Optimized (True Hardware Gather)
 * ============================================================================
 */

/**
 * Gather int32 values from dictionary using AVX2 gather instructions.
 */
void carquet_avx2_gather_i32(const int32_t* dict, const uint32_t* indices,
                              int64_t count, int32_t* output) {
    int64_t i = 0;

    /* Process 8 at a time using AVX2 gather */
    for (; i + 8 <= count; i += 8) {
        __m256i idx = _mm256_loadu_si256((const __m256i*)(indices + i));
        __m256i result = _mm256_i32gather_epi32(dict, idx, 4);  /* Scale = 4 bytes per int32 */
        _mm256_storeu_si256((__m256i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather int64 values from dictionary using AVX2 gather instructions.
 */
void carquet_avx2_gather_i64(const int64_t* dict, const uint32_t* indices,
                              int64_t count, int64_t* output) {
    int64_t i = 0;

    /* Process 4 at a time using AVX2 gather */
    for (; i + 4 <= count; i += 4) {
        __m128i idx = _mm_loadu_si128((const __m128i*)(indices + i));
        __m256i result = _mm256_i32gather_epi64((const long long*)dict, idx, 8);
        _mm256_storeu_si256((__m256i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather float values from dictionary using AVX2 gather instructions.
 * Note: float and int32 are both 4 bytes, so we reuse gather_i32 via cast.
 */
void carquet_avx2_gather_float(const float* dict, const uint32_t* indices,
                                int64_t count, float* output) {
    /* Data movement doesn't care about type - reuse int32 implementation */
    carquet_avx2_gather_i32((const int32_t*)dict, indices, count, (int32_t*)output);
}

/**
 * Gather double values from dictionary using AVX2 gather instructions.
 * Note: double and int64 are both 8 bytes, so we reuse gather_i64 via cast.
 */
void carquet_avx2_gather_double(const double* dict, const uint32_t* indices,
                                 int64_t count, double* output) {
    /* Data movement doesn't care about type - reuse int64 implementation */
    carquet_avx2_gather_i64((const int64_t*)dict, indices, count, (int64_t*)output);
}

static inline int avx2_indices_in_bounds_8(const uint32_t* indices, uint32_t limit) {
    __m256i idx = _mm256_loadu_si256((const __m256i*)indices);
    __m256i bias = _mm256_set1_epi32((int)0x80000000u);
    __m256i idx_biased = _mm256_xor_si256(idx, bias);
    __m256i limit_biased = _mm256_set1_epi32((int)(limit ^ 0x80000000u));
    __m256i cmp = _mm256_cmpgt_epi32(limit_biased, idx_biased);
    return _mm256_movemask_epi8(cmp) == -1;
}

bool carquet_avx2_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                      const uint32_t* indices, int64_t count,
                                      int32_t* output) {
    int64_t i = 0;
    uint32_t limit = (uint32_t)dict_count;

    for (; i + 8 <= count; i += 8) {
        if (!avx2_indices_in_bounds_8(indices + i, limit)) {
            return false;
        }

        __m256i idx = _mm256_loadu_si256((const __m256i*)(indices + i));
        __m256i result = _mm256_i32gather_epi32(dict, idx, 4);
        _mm256_storeu_si256((__m256i*)(output + i), result);
    }

    for (; i + 4 <= count; i += 4) {
        uint32_t a = indices[i + 0];
        uint32_t b = indices[i + 1];
        uint32_t c = indices[i + 2];
        uint32_t d = indices[i + 3];
        if (a >= limit || b >= limit || c >= limit || d >= limit) {
            return false;
        }
        __m128i result = _mm_set_epi32(cq_loadu(dict + (d)), cq_loadu(dict + (c)), cq_loadu(dict + (b)), cq_loadu(dict + (a)));
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

bool carquet_avx2_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                      const uint32_t* indices, int64_t count,
                                      int64_t* output) {
    int64_t i = 0;
    uint32_t limit = (uint32_t)dict_count;

    for (; i + 8 <= count; i += 8) {
        if (!avx2_indices_in_bounds_8(indices + i, limit)) {
            return false;
        }

        __m128i idx0 = _mm_loadu_si128((const __m128i*)(indices + i));
        __m128i idx1 = _mm_loadu_si128((const __m128i*)(indices + i + 4));
        __m256i result0 = _mm256_i32gather_epi64((const long long*)dict, idx0, 8);
        __m256i result1 = _mm256_i32gather_epi64((const long long*)dict, idx1, 8);
        _mm256_storeu_si256((__m256i*)(output + i), result0);
        _mm256_storeu_si256((__m256i*)(output + i + 4), result1);
    }

    for (; i + 4 <= count; i += 4) {
        uint32_t a = indices[i + 0];
        uint32_t b = indices[i + 1];
        uint32_t c = indices[i + 2];
        uint32_t d = indices[i + 3];
        if (a >= limit || b >= limit || c >= limit || d >= limit) {
            return false;
        }
        __m256i result = _mm256_set_epi64x(cq_loadu(dict + (d)), cq_loadu(dict + (c)), cq_loadu(dict + (b)), cq_loadu(dict + (a)));
        _mm256_storeu_si256((__m256i*)(output + i), result);
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

bool carquet_avx2_checked_gather_float(const float* dict, int32_t dict_count,
                                        const uint32_t* indices, int64_t count,
                                        float* output) {
    return carquet_avx2_checked_gather_i32(
        (const int32_t*)dict, dict_count, indices, count, (int32_t*)output);
}

bool carquet_avx2_checked_gather_double(const double* dict, int32_t dict_count,
                                         const uint32_t* indices, int64_t count,
                                         double* output) {
    return carquet_avx2_checked_gather_i64(
        (const int64_t*)dict, dict_count, indices, count, (int64_t*)output);
}

void carquet_avx2_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset) {
    if (offset >= 32) {
        while (len >= 32) {
            _mm256_storeu_si256((__m256i*)dst, _mm256_loadu_si256((const __m256i*)src));
            dst += 32;
            src += 32;
            len -= 32;
        }
        while (len >= 16) {
            _mm_storeu_si128((__m128i*)dst, _mm_loadu_si128((const __m128i*)src));
            dst += 16;
            src += 16;
            len -= 16;
        }
    } else if (offset == 1) {
        __m256i v = _mm256_set1_epi8((char)*src);
        while (len >= 32) {
            _mm256_storeu_si256((__m256i*)dst, v);
            dst += 32;
            len -= 32;
        }
    } else if (offset == 2) {
        uint16_t pattern;
        memcpy(&pattern, src, sizeof(pattern));
        while (len >= 2) {
            memcpy(dst, &pattern, sizeof(pattern));
            dst += 2;
            len -= 2;
        }
        if (len) {
            *dst = *(const uint8_t*)&pattern;
            return;
        }
        return;
    } else if (offset == 4) {
        uint32_t pattern;
        memcpy(&pattern, src, sizeof(pattern));
        __m256i v = _mm256_set1_epi32((int32_t)pattern);
        while (len >= 32) {
            _mm256_storeu_si256((__m256i*)dst, v);
            dst += 32;
            len -= 32;
        }
    } else if (offset == 8) {
        uint64_t pattern;
        memcpy(&pattern, src, sizeof(pattern));
        __m256i v = _mm256_set1_epi64x((long long)pattern);
        while (len >= 32) {
            _mm256_storeu_si256((__m256i*)dst, v);
            dst += 32;
            len -= 32;
        }
    }

    while (len > 0) {
        *dst++ = *src++;
        len--;
    }
}

size_t carquet_avx2_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit) {
    const uint8_t* start = p;

    while (p + 32 <= limit) {
        __m256i a = _mm256_loadu_si256((const __m256i*)p);
        __m256i b = _mm256_loadu_si256((const __m256i*)match);
        __m256i cmp = _mm256_cmpeq_epi8(a, b);
        uint32_t mask = (uint32_t)_mm256_movemask_epi8(cmp);

        if (mask != 0xFFFFFFFFu) {
            return (size_t)(p - start) + (size_t)__builtin_ctz(~mask);
        }

        p += 32;
        match += 32;
    }

    while (p < limit && *p == *match) {
        p++;
        match++;
    }

    return (size_t)(p - start);
}

/* ============================================================================
 * Memcpy/Memset - AVX2 Optimized
 * ============================================================================
 */

/**
 * Fast memset for buffers using AVX2.
 */
void carquet_avx2_memset(void* dest, uint8_t value, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    __m256i v = _mm256_set1_epi8((char)value);

    while (n >= 128) {
        _mm256_storeu_si256((__m256i*)(d + 0), v);
        _mm256_storeu_si256((__m256i*)(d + 32), v);
        _mm256_storeu_si256((__m256i*)(d + 64), v);
        _mm256_storeu_si256((__m256i*)(d + 96), v);
        d += 128;
        n -= 128;
    }

    while (n >= 32) {
        _mm256_storeu_si256((__m256i*)d, v);
        d += 32;
        n -= 32;
    }

    /* Handle tail with SSE */
    __m128i v128 = _mm_set1_epi8((char)value);
    while (n >= 16) {
        _mm_storeu_si128((__m128i*)d, v128);
        d += 16;
        n -= 16;
    }

    while (n > 0) {
        *d++ = value;
        n--;
    }
}

/**
 * Fast memcpy for buffers using AVX2.
 */
void carquet_avx2_memcpy(void* dest, const void* src, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    const uint8_t* s = (const uint8_t*)src;

    while (n >= 128) {
        __m256i v0 = _mm256_loadu_si256((const __m256i*)(s + 0));
        __m256i v1 = _mm256_loadu_si256((const __m256i*)(s + 32));
        __m256i v2 = _mm256_loadu_si256((const __m256i*)(s + 64));
        __m256i v3 = _mm256_loadu_si256((const __m256i*)(s + 96));
        _mm256_storeu_si256((__m256i*)(d + 0), v0);
        _mm256_storeu_si256((__m256i*)(d + 32), v1);
        _mm256_storeu_si256((__m256i*)(d + 64), v2);
        _mm256_storeu_si256((__m256i*)(d + 96), v3);
        d += 128;
        s += 128;
        n -= 128;
    }

    while (n >= 32) {
        _mm256_storeu_si256((__m256i*)d, _mm256_loadu_si256((const __m256i*)s));
        d += 32;
        s += 32;
        n -= 32;
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
 * Boolean Unpacking - AVX2 Optimized
 * ============================================================================
 */

/**
 * Unpack boolean values from packed bits to byte array using AVX2.
 * Each output byte is 0 or 1.
 */
void carquet_avx2_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;
    const __m256i mask = _mm256_set_epi8(
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01,
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01,
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01,
        (char)0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01
    );
    const __m256i shuf = _mm256_setr_epi8(
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3
    );

    /* Process 32 bools (4 bytes) at a time */
    for (; i + 32 <= count; i += 32) {
        int byte_idx = (int)(i / 8);
        uint32_t packed;
        memcpy(&packed, input + byte_idx, 4);

        __m256i bits = _mm256_set1_epi32(packed);

        /* Create masks for each bit position */
        __m256i shuffled = _mm256_shuffle_epi8(bits, shuf);

        /* AND with mask and normalize to 0/1 */
        __m256i masked = _mm256_and_si256(shuffled, mask);
        __m256i result = _mm256_min_epu8(masked, _mm256_set1_epi8(1));

        _mm256_storeu_si256((__m256i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        int byte_idx = (int)(i / 8);
        int bit_idx = (int)(i % 8);
        output[i] = (input[byte_idx] >> bit_idx) & 1;
    }
}

/**
 * Pack boolean values from byte array to packed bits using AVX2.
 */
void carquet_avx2_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    /* 32 bools -> 4 output bytes per iteration via compare + movemask.
     * cmpeq(v, 0) yields 0xFF for every zero byte; movemask gathers the top
     * bit of each byte (byte k -> mask bit k), so the inverted mask has bit k
     * set iff byte k is non-zero. That is exactly the scalar reference rule
     * (any non-zero byte is true, not just 0x01) and movemask's byte order is
     * the LSB-first layout the PLAIN BOOLEAN encoding needs. This replaces a
     * 128-bit mullo/horizontal-add reduction that handled only 8 bools per
     * iteration. */
    const __m256i zero256 = _mm256_setzero_si256();
    for (; i + 32 <= count; i += 32) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(input + i));
        uint32_t nz = ~(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v, zero256));
        uint8_t* o = output + i / 8;
        o[0] = (uint8_t)nz;
        o[1] = (uint8_t)(nz >> 8);
        o[2] = (uint8_t)(nz >> 16);
        o[3] = (uint8_t)(nz >> 24);
    }

    /* 8 bools -> 1 byte with the 128-bit half of the same trick. loadl_epi64
     * zero-fills the upper 8 lanes; their mask bits land above bit 7 and are
     * dropped by the uint8_t cast. */
    const __m128i zero128 = _mm_setzero_si128();
    for (; i + 8 <= count; i += 8) {
        __m128i v = _mm_loadl_epi64((const __m128i*)(input + i));
        int eq = _mm_movemask_epi8(_mm_cmpeq_epi8(v, zero128));
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
 * RLE Run Detection - AVX2 Optimized
 * ============================================================================
 */

/**
 * Find the length of a run of repeated values.
 * Returns the number of consecutive identical values starting at the given position.
 */
int64_t carquet_avx2_find_run_length_i32(const int32_t* values, int64_t count) {
    if (count == 0) return 0;

    int32_t first = values[0];
    __m256i target = _mm256_set1_epi32(first);
    int64_t i = 0;

    /* Check 8 at a time */
    for (; i + 8 <= count; i += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));
        __m256i cmp = _mm256_cmpeq_epi32(v, target);
        uint32_t mask = (uint32_t)_mm256_movemask_epi8(cmp);

        if (mask != 0xFFFFFFFFu) {
            return i + (__builtin_ctz(~mask) >> 2);
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

int64_t carquet_avx2_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    int64_t non_null_count = 0;
    int64_t i = 0;
    __m256i max_vec = _mm256_set1_epi16(max_def_level);

    for (; i + 16 <= count; i += 16) {
        __m256i levels = _mm256_loadu_si256((const __m256i*)(def_levels + i));
        __m256i cmp = _mm256_cmpeq_epi16(levels, max_vec);
        uint32_t mask = (uint32_t)_mm256_movemask_epi8(cmp);
        non_null_count += __builtin_popcount(mask) >> 1;
    }

    for (; i < count; i++) {
        if (def_levels[i] == max_def_level) {
            non_null_count++;
        }
    }

    return non_null_count;
}

void carquet_avx2_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                     int16_t max_def_level, uint8_t* null_bitmap) {
    int64_t i = 0;
    int64_t full_bytes = count / 8;
    __m256i max_vec = _mm256_set1_epi16(max_def_level);
    __m128i zero = _mm_setzero_si128();

    for (int64_t b = 0; b + 1 < full_bytes; b += 2) {
        __m256i levels = _mm256_loadu_si256((const __m256i*)(def_levels + i));
        __m256i cmp = _mm256_cmpeq_epi16(levels, max_vec);
        __m128i lo = _mm256_castsi256_si128(cmp);
        __m128i hi = _mm256_extracti128_si256(cmp, 1);
        __m128i packed = _mm_packs_epi16(lo, hi);
        int mask = _mm_movemask_epi8(packed);
        null_bitmap[b] = (uint8_t)(mask & 0xFF);
        null_bitmap[b + 1] = (uint8_t)((mask >> 8) & 0xFF);
        i += 16;
    }

    for (int64_t b = (full_bytes & ~1LL); b < full_bytes; b++) {
        __m128i levels = _mm_loadu_si128((const __m128i*)(def_levels + i));
        __m128i max128 = _mm256_castsi256_si128(max_vec);
        __m128i cmp = _mm_cmpeq_epi16(levels, max128);
        __m128i packed = _mm_packs_epi16(cmp, zero);
        null_bitmap[b] = (uint8_t)_mm_movemask_epi8(packed);
        i += 8;
    }

    if (i < count) {
        uint8_t present_bits = 0;
        for (int64_t j = 0; i + j < count && j < 8; j++) {
            if (def_levels[i + j] == max_def_level) {
                present_bits |= (uint8_t)(1u << j);
            }
        }
        null_bitmap[full_bytes] = present_bits;
    }
}

void carquet_avx2_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    int64_t i = 0;
    __m256i val_vec = _mm256_set1_epi16(value);

    for (; i + 16 <= count; i += 16) {
        _mm256_storeu_si256((__m256i*)(def_levels + i), val_vec);
    }
    for (; i + 8 <= count; i += 8) {
        _mm_storeu_si128((__m128i*)(def_levels + i), _mm256_castsi256_si128(val_vec));
    }
    for (; i < count; i++) {
        def_levels[i] = value;
    }
}

void carquet_avx2_minmax_i32(const int32_t* values, int64_t count,
                              int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    __m256i min_vec = _mm256_set1_epi32(min_v);
    __m256i max_vec = _mm256_set1_epi32(max_v);
    int64_t i = 1;

    for (; i + 8 <= count; i += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));
        min_vec = _mm256_min_epi32(min_vec, v);
        max_vec = _mm256_max_epi32(max_vec, v);
    }

    int32_t tmp_min[8];
    int32_t tmp_max[8];
    _mm256_storeu_si256((__m256i*)tmp_min, min_vec);
    _mm256_storeu_si256((__m256i*)tmp_max, max_vec);
    for (int j = 0; j < 8; j++) {
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

void carquet_avx2_minmax_i64(const int64_t* values, int64_t count,
                              int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    __m256i min_vec = _mm256_set1_epi64x(min_v);
    __m256i max_vec = _mm256_set1_epi64x(max_v);
    int64_t i = 1;

    for (; i + 4 <= count; i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));
        __m256i lt = _mm256_cmpgt_epi64(min_vec, v);
        __m256i gt = _mm256_cmpgt_epi64(v, max_vec);
        min_vec = _mm256_blendv_epi8(min_vec, v, lt);
        max_vec = _mm256_blendv_epi8(max_vec, v, gt);
    }

    int64_t tmp_min[4];
    int64_t tmp_max[4];
    _mm256_storeu_si256((__m256i*)tmp_min, min_vec);
    _mm256_storeu_si256((__m256i*)tmp_max, max_vec);
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

void carquet_avx2_minmax_float(const float* values, int64_t count,
                                float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m256 min_vec = _mm256_set1_ps(min_v);
    __m256 max_vec = _mm256_set1_ps(max_v);
    int64_t i = 1;

    for (; i + 8 <= count; i += 8) {
        __m256 v = _mm256_loadu_ps(values + i);
        __m256 lt = _mm256_cmp_ps(v, min_vec, _CMP_LT_OQ);
        __m256 gt = _mm256_cmp_ps(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm256_blendv_ps(min_vec, v, lt);
        max_vec = _mm256_blendv_ps(max_vec, v, gt);
    }

    float tmp_min[8];
    float tmp_max[8];
    _mm256_storeu_ps(tmp_min, min_vec);
    _mm256_storeu_ps(tmp_max, max_vec);
    for (int j = 0; j < 8; j++) {
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

void carquet_avx2_minmax_double(const double* values, int64_t count,
                                 double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m256d min_vec = _mm256_set1_pd(min_v);
    __m256d max_vec = _mm256_set1_pd(max_v);
    int64_t i = 1;

    for (; i + 4 <= count; i += 4) {
        __m256d v = _mm256_loadu_pd(values + i);
        __m256d lt = _mm256_cmp_pd(v, min_vec, _CMP_LT_OQ);
        __m256d gt = _mm256_cmp_pd(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm256_blendv_pd(min_vec, v, lt);
        max_vec = _mm256_blendv_pd(max_vec, v, gt);
    }

    double tmp_min[4];
    double tmp_max[4];
    _mm256_storeu_pd(tmp_min, min_vec);
    _mm256_storeu_pd(tmp_max, max_vec);
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

void carquet_avx2_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                   int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    __m256i min_vec = _mm256_set1_epi32(min_v);
    __m256i max_vec = _mm256_set1_epi32(max_v);
    int64_t i = 0;

    for (; i + 8 <= count; i += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));
        _mm256_storeu_si256((__m256i*)(output + i), v);
        min_vec = _mm256_min_epi32(min_vec, v);
        max_vec = _mm256_max_epi32(max_vec, v);
    }

    int32_t tmp_min[8];
    int32_t tmp_max[8];
    _mm256_storeu_si256((__m256i*)tmp_min, min_vec);
    _mm256_storeu_si256((__m256i*)tmp_max, max_vec);
    for (int j = 0; j < 8; j++) {
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

void carquet_avx2_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                   int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    __m256i min_vec = _mm256_set1_epi64x(min_v);
    __m256i max_vec = _mm256_set1_epi64x(max_v);
    int64_t i = 0;

    for (; i + 4 <= count; i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(values + i));
        _mm256_storeu_si256((__m256i*)(output + i), v);
        __m256i lt = _mm256_cmpgt_epi64(min_vec, v);
        __m256i gt = _mm256_cmpgt_epi64(v, max_vec);
        min_vec = _mm256_blendv_epi8(min_vec, v, lt);
        max_vec = _mm256_blendv_epi8(max_vec, v, gt);
    }

    int64_t tmp_min[4];
    int64_t tmp_max[4];
    _mm256_storeu_si256((__m256i*)tmp_min, min_vec);
    _mm256_storeu_si256((__m256i*)tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
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

void carquet_avx2_copy_minmax_float(const float* values, int64_t count, float* output,
                                     float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m256 min_vec = _mm256_set1_ps(min_v);
    __m256 max_vec = _mm256_set1_ps(max_v);
    int64_t i = 0;

    for (; i + 8 <= count; i += 8) {
        __m256 v = _mm256_loadu_ps(values + i);
        _mm256_storeu_ps(output + i, v);
        __m256 lt = _mm256_cmp_ps(v, min_vec, _CMP_LT_OQ);
        __m256 gt = _mm256_cmp_ps(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm256_blendv_ps(min_vec, v, lt);
        max_vec = _mm256_blendv_ps(max_vec, v, gt);
    }

    float tmp_min[8];
    float tmp_max[8];
    _mm256_storeu_ps(tmp_min, min_vec);
    _mm256_storeu_ps(tmp_max, max_vec);
    for (int j = 0; j < 8; j++) {
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

void carquet_avx2_copy_minmax_double(const double* values, int64_t count, double* output,
                                      double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m256d min_vec = _mm256_set1_pd(min_v);
    __m256d max_vec = _mm256_set1_pd(max_v);
    int64_t i = 0;

    for (; i + 4 <= count; i += 4) {
        __m256d v = _mm256_loadu_pd(values + i);
        _mm256_storeu_pd(output + i, v);
        __m256d lt = _mm256_cmp_pd(v, min_vec, _CMP_LT_OQ);
        __m256d gt = _mm256_cmp_pd(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm256_blendv_pd(min_vec, v, lt);
        max_vec = _mm256_blendv_pd(max_vec, v, gt);
    }

    double tmp_min[4];
    double tmp_max[4];
    _mm256_storeu_pd(tmp_min, min_vec);
    _mm256_storeu_pd(tmp_max, max_vec);
    for (int j = 0; j < 4; j++) {
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

#endif /* __AVX2__ */
#endif /* x86 */
