/**
 * @file avx512_ops.c
 * @brief AVX-512 optimized operations for x86-64 processors
 *
 * Provides SIMD-accelerated implementations using 512-bit vectors:
 * - Bit unpacking for various bit widths
 * - Byte stream split/merge (for BYTE_STREAM_SPLIT encoding)
 * - Delta decoding (prefix sums)
 * - Dictionary gather operations (using AVX-512 scatter/gather)
 * - Boolean packing/unpacking
 * - Masked operations for predicated processing
 */

#include <carquet/error.h>
#include "simd/simd_unaligned.h"
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__x86_64__) || defined(_M_X64)
/* Check for AVX-512 support */
#if defined(__AVX512F__) || (defined(_MSC_VER) && defined(__AVX512F__))

#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <immintrin.h>

/* Portable count trailing zeros */
static inline int portable_ctz(unsigned int v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(v);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward(&index, v);
    return (int)index;
#else
    int n = 0;
    if (!(v & 0xFFFF)) { n += 16; v >>= 16; }
    if (!(v & 0xFF)) { n += 8; v >>= 8; }
    if (!(v & 0xF)) { n += 4; v >>= 4; }
    if (!(v & 0x3)) { n += 2; v >>= 2; }
    if (!(v & 0x1)) { n += 1; }
    return n;
#endif
}

static inline int portable_popcount(unsigned int v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(v);
#elif defined(_MSC_VER)
    return (int)__popcnt(v);
#else
    int count = 0;
    while (v) {
        v &= v - 1;
        count++;
    }
    return count;
#endif
}

/* Portable 64-bit count trailing zeros (needed for 64-byte mask operations) */
static inline int portable_ctz64(uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(v);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward64(&index, (unsigned __int64)v);
    return (int)index;
#else
    if ((uint32_t)v) return portable_ctz((unsigned int)v);
    return 32 + portable_ctz((unsigned int)(v >> 32));
#endif
}

/* ============================================================================
 * Bit Unpacking - AVX-512 Optimized
 * ============================================================================
 */

void carquet_avx512_bitunpack32_8bit(const uint8_t* input, uint32_t* values);

/**
 * Unpack 8 8-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack8_8bit(const uint8_t* input, uint32_t* values) {
    __m128i bytes = _mm_loadl_epi64((const __m128i*)input);
    __m512i expanded = _mm512_cvtepu8_epi32(bytes);
    __m256i result = _mm512_castsi512_si256(expanded);
    _mm256_storeu_si256((__m256i*)values, result);
}

/**
 * Unpack 8 16-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack8_16bit(const uint8_t* input, uint32_t* values) {
    __m128i words = _mm_loadu_si128((const __m128i*)input);
    __m256i result = _mm256_cvtepu16_epi32(words);
    _mm256_storeu_si256((__m256i*)values, result);
}

/**
 * Unpack 8 4-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack8_4bit(const uint8_t* input, uint32_t* values) {
    /* 8 x 4-bit values = 4 input bytes, expand nibbles to bytes then widen */
    uint8_t expanded[8];
    for (int i = 0; i < 4; i++) {
        uint8_t byte = input[i];
        expanded[i * 2] = (uint8_t)(byte & 0x0F);
        expanded[i * 2 + 1] = (uint8_t)(byte >> 4);
    }
    __m128i bytes = _mm_loadl_epi64((const __m128i*)expanded);
    __m256i result = _mm256_cvtepu8_epi32(bytes);
    _mm256_storeu_si256((__m256i*)values, result);
}

/**
 * Unpack 32 8-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack32_8bit(const uint8_t* input, uint32_t* values) {
    /* Load 32 bytes as two 128-bit halves */
    __m128i bytes_lo = _mm_loadu_si128((const __m128i*)input);
    __m128i bytes_hi = _mm_loadu_si128((const __m128i*)(input + 16));

    /* Expand each half to 32-bit using AVX-512 (16 x 8-bit -> 16 x 32-bit) */
    __m512i result_lo = _mm512_cvtepu8_epi32(bytes_lo);
    __m512i result_hi = _mm512_cvtepu8_epi32(bytes_hi);

    _mm512_storeu_si512((__m512i*)values, result_lo);
    _mm512_storeu_si512((__m512i*)(values + 16), result_hi);
}

/**
 * Unpack 16 16-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack16_16bit(const uint8_t* input, uint32_t* values) {
    __m256i words = _mm256_loadu_si256((const __m256i*)input);
    __m512i result = _mm512_cvtepu16_epi32(words);
    _mm512_storeu_si512((__m512i*)values, result);
}

/**
 * Unpack 32 4-bit values to 32-bit using AVX-512.
 */
void carquet_avx512_bitunpack32_4bit(const uint8_t* input, uint32_t* values) {
    /* Load 16 bytes containing 32 x 4-bit values */
    __m128i bytes = _mm_loadu_si128((const __m128i*)input);

    /* Split nibbles */
    __m128i lo_nibbles = _mm_and_si128(bytes, _mm_set1_epi8(0x0F));
    __m128i hi_nibbles = _mm_srli_epi16(bytes, 4);
    hi_nibbles = _mm_and_si128(hi_nibbles, _mm_set1_epi8(0x0F));

    /* Interleave to get correct order - produces two 128-bit results */
    __m128i interleaved_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
    __m128i interleaved_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);

    /* Expand each half to 32-bit using AVX-512 (16 x 8-bit -> 16 x 32-bit) */
    __m512i result_lo = _mm512_cvtepu8_epi32(interleaved_lo);
    __m512i result_hi = _mm512_cvtepu8_epi32(interleaved_hi);

    _mm512_storeu_si512((__m512i*)values, result_lo);
    _mm512_storeu_si512((__m512i*)(values + 16), result_hi);
}

/* ============================================================================
 * Byte Stream Split - AVX-512 Optimized
 * ============================================================================
 */

/**
 * Encode floats using byte stream split with AVX-512.
 * Processes 16 floats (64 bytes) at a time using VBMI byte permutation.
 */
void carquet_avx512_byte_stream_split_encode_float(
    const float* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;

#ifdef __AVX512VBMI__
    /* Single permutation that places all 4 byte streams in the 4 128-bit lanes:
     * Lane 0 (bits 0-127):   byte 0 from each of 16 floats
     * Lane 1 (bits 128-255): byte 1 from each of 16 floats
     * Lane 2 (bits 256-383): byte 2 from each of 16 floats
     * Lane 3 (bits 384-511): byte 3 from each of 16 floats
     */
    /* Use _mm512_set_epi32 instead of _mm512_set_epi8 for GCC 8 compatibility */
    const __m512i perm_all = _mm512_set_epi32(
        0x3F3B3733, 0x2F2B2723, 0x1F1B1713, 0x0F0B0703,  /* byte 3s */
        0x3E3A3632, 0x2E2A2622, 0x1E1A1612, 0x0E0A0602,  /* byte 2s */
        0x3D393531, 0x2D292521, 0x1D191511, 0x0D090501,   /* byte 1s */
        0x3C383430, 0x2C282420, 0x1C181410, 0x0C080400);  /* byte 0s */

    for (; i + 16 <= count; i += 16) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(src + i * 4));

        /* Single permutation gathers all 4 streams */
        __m512i transposed = _mm512_permutexvar_epi8(perm_all, v);

        /* Extract and store each 128-bit lane to its stream */
        _mm_storeu_si128((__m128i*)(output + 0 * count + i), _mm512_castsi512_si128(transposed));
        _mm_storeu_si128((__m128i*)(output + 1 * count + i), _mm512_extracti32x4_epi32(transposed, 1));
        _mm_storeu_si128((__m128i*)(output + 2 * count + i), _mm512_extracti32x4_epi32(transposed, 2));
        _mm_storeu_si128((__m128i*)(output + 3 * count + i), _mm512_extracti32x4_epi32(transposed, 3));
    }
#else
    /* Fallback without VBMI: use shuffle + permutexvar approach
     * Step 1: shuffle_epi8 transposes within each 128-bit lane (4 floats -> 4 bytes per stream)
     * Step 2: permutexvar_epi32 rearranges dwords to group all byte 0s, byte 1s, etc.
     */
    /* Use _mm512_set_epi32 instead of _mm512_set_epi8 for GCC 8 compatibility */
    const __m512i intra_lane_shuf = _mm512_set_epi32(
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400);
    const __m512i cross_lane_perm = _mm512_set_epi32(
        15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0);

    for (; i + 16 <= count; i += 16) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(src + i * 4));

        /* Transpose within each 128-bit lane */
        __m512i shuffled = _mm512_shuffle_epi8(v, intra_lane_shuf);

        /* Rearrange dwords across lanes to group streams */
        __m512i transposed = _mm512_permutexvar_epi32(cross_lane_perm, shuffled);

        /* Extract and store each 128-bit lane to its stream */
        _mm_storeu_si128((__m128i*)(output + 0 * count + i), _mm512_castsi512_si128(transposed));
        _mm_storeu_si128((__m128i*)(output + 1 * count + i), _mm512_extracti32x4_epi32(transposed, 1));
        _mm_storeu_si128((__m128i*)(output + 2 * count + i), _mm512_extracti32x4_epi32(transposed, 2));
        _mm_storeu_si128((__m128i*)(output + 3 * count + i), _mm512_extracti32x4_epi32(transposed, 3));
    }
#endif

    /* Handle remaining values */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            output[b * count + i] = src[i * 4 + b];
        }
    }
}

/**
 * Decode byte stream split floats using AVX-512.
 * Processes 16 floats (64 bytes) at a time using 512-bit operations.
 */
void carquet_avx512_byte_stream_split_decode_float(
    const uint8_t* data,
    int64_t count,
    float* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

#ifdef __AVX512VBMI__
    /* VBMI path: single vpermb for 16 floats.
     * Input in __m512i:
     *   bytes  0-15: stream 0 (byte 0 of each float)
     *   bytes 16-31: stream 1 (byte 1 of each float)
     *   bytes 32-47: stream 2 (byte 2 of each float)
     *   bytes 48-63: stream 3 (byte 3 of each float)
     * Output: float[k] = {byte0[k], byte1[k], byte2[k], byte3[k]}
     *   output byte 4*k+0 = input byte k
     *   output byte 4*k+1 = input byte 16+k
     *   output byte 4*k+2 = input byte 32+k
     *   output byte 4*k+3 = input byte 48+k */
    const __m512i perm = _mm512_set_epi32(
        0x3F2F1F0F, 0x3E2E1E0E, 0x3D2D1D0D, 0x3C2C1C0C,
        0x3B2B1B0B, 0x3A2A1A0A, 0x39291909, 0x38281808,
        0x37271707, 0x36261606, 0x35251505, 0x34241404,
        0x33231303, 0x32221202, 0x31211101, 0x30201000);

    for (; i + 16 <= count; i += 16) {
        __m128i s0 = _mm_loadu_si128((const __m128i*)(data + 0 * count + i));
        __m128i s1 = _mm_loadu_si128((const __m128i*)(data + 1 * count + i));
        __m128i s2 = _mm_loadu_si128((const __m128i*)(data + 2 * count + i));
        __m128i s3 = _mm_loadu_si128((const __m128i*)(data + 3 * count + i));

        __m512i combined = _mm512_castsi128_si512(s0);
        combined = _mm512_inserti32x4(combined, s1, 1);
        combined = _mm512_inserti32x4(combined, s2, 2);
        combined = _mm512_inserti32x4(combined, s3, 3);

        __m512i result = _mm512_permutexvar_epi8(perm, combined);
        _mm512_storeu_si512((__m512i*)(dst + i * 4), result);
    }
#else
    /* Non-VBMI fallback: use the inverse of the encode approach.
     * The encode uses: (1) intra-lane shuffle to group bytes, (2) cross-lane dword permute.
     * For decode, reverse the process:
     * (1) Load streams into 512-bit lanes, (2) cross-lane permute, (3) intra-lane shuffle. */

    /* Step 1 cross-lane: inverse of the encode's cross_lane_perm.
     * Encode permutes dwords as {12,8,4,0, 13,9,5,1, 14,10,6,2, 15,11,7,3}.
     * The inverse moves dword K to position inverse[K]. */
    const __m512i cross_lane_inv = _mm512_set_epi32(
        15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0);

    /* Step 2 intra-lane: inverse of the encode's intra_lane_shuf.
     * Encode: byte[4k+j] -> position[j*4+k] (for k=0..3, j=0..3)
     * Decode: position[j*4+k] -> byte[4k+j], i.e. byte[p] -> pos[((p%4)*4 + p/4)] */
    const __m512i intra_lane_inv = _mm512_set_epi32(
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400,
        0x0F0B0703, 0x0E0A0602, 0x0D090501, 0x0C080400);

    for (; i + 16 <= count; i += 16) {
        __m128i s0 = _mm_loadu_si128((const __m128i*)(data + 0 * count + i));
        __m128i s1 = _mm_loadu_si128((const __m128i*)(data + 1 * count + i));
        __m128i s2 = _mm_loadu_si128((const __m128i*)(data + 2 * count + i));
        __m128i s3 = _mm_loadu_si128((const __m128i*)(data + 3 * count + i));

        __m512i combined = _mm512_castsi128_si512(s0);
        combined = _mm512_inserti32x4(combined, s1, 1);
        combined = _mm512_inserti32x4(combined, s2, 2);
        combined = _mm512_inserti32x4(combined, s3, 3);

        /* Rearrange dwords across lanes */
        __m512i permuted = _mm512_permutexvar_epi32(cross_lane_inv, combined);
        /* Shuffle bytes within each lane to reconstruct floats */
        __m512i result = _mm512_shuffle_epi8(permuted, intra_lane_inv);

        _mm512_storeu_si512((__m512i*)(dst + i * 4), result);
    }
#endif

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            dst[i * 4 + b] = data[b * count + i];
        }
    }
}

/**
 * Encode doubles using byte stream split with AVX-512.
 * Processes 8 doubles (64 bytes) at a time using 512-bit operations.
 */
void carquet_avx512_byte_stream_split_encode_double(
    const double* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;

#ifdef __AVX512VBMI__
    /* Single vpermb transposes 8 doubles (64 bytes) into 8 byte streams.
     * Output layout: 8 lanes of 8 bytes, each lane is one byte stream. */
    const __m512i perm = _mm512_set_epi32(
        0x3F372F27, 0x1F170F07,  /* stream 7 */
        0x3E362E26, 0x1E160E06,  /* stream 6 */
        0x3D352D25, 0x1D150D05,  /* stream 5 */
        0x3C342C24, 0x1C140C04,  /* stream 4 */
        0x3B332B23, 0x1B130B03,  /* stream 3 */
        0x3A322A22, 0x1A120A02,  /* stream 2 */
        0x39312921, 0x19110901,  /* stream 1 */
        0x38302820, 0x18100800); /* stream 0 */

    for (; i + 8 <= count; i += 8) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(src + i * 8));
        __m512i transposed = _mm512_permutexvar_epi8(perm, v);

        /* Extract 8 bytes per stream from the 512-bit result.
         * After permutation, the result is organized as:
         * bytes  0-7:  stream 0 (byte 0 from each of 8 doubles)
         * bytes  8-15: stream 1
         * ... etc */
        _mm_storel_epi64((__m128i*)(output + 0 * count + i),
                         _mm512_castsi512_si128(transposed));
        _mm_storel_epi64((__m128i*)(output + 1 * count + i),
                         _mm_srli_si128(_mm512_castsi512_si128(transposed), 8));
        __m128i lane1 = _mm512_extracti32x4_epi32(transposed, 1);
        _mm_storel_epi64((__m128i*)(output + 2 * count + i), lane1);
        _mm_storel_epi64((__m128i*)(output + 3 * count + i),
                         _mm_srli_si128(lane1, 8));
        __m128i lane2 = _mm512_extracti32x4_epi32(transposed, 2);
        _mm_storel_epi64((__m128i*)(output + 4 * count + i), lane2);
        _mm_storel_epi64((__m128i*)(output + 5 * count + i),
                         _mm_srli_si128(lane2, 8));
        __m128i lane3 = _mm512_extracti32x4_epi32(transposed, 3);
        _mm_storel_epi64((__m128i*)(output + 6 * count + i), lane3);
        _mm_storel_epi64((__m128i*)(output + 7 * count + i),
                         _mm_srli_si128(lane3, 8));
    }
#else
    /* Non-VBMI: use AVX2-style shuffle approach with 256-bit halves */
    const __m256i sh0 = _mm256_setr_epi8(0,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh1 = _mm256_setr_epi8(1,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh2 = _mm256_setr_epi8(2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh3 = _mm256_setr_epi8(3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh4 = _mm256_setr_epi8(4,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh5 = _mm256_setr_epi8(5,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh6 = _mm256_setr_epi8(6,14,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,14,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const __m256i sh7 = _mm256_setr_epi8(7,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 7,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);

    for (; i + 8 <= count; i += 8) {
        __m256i lo = _mm256_loadu_si256((const __m256i*)(src + i * 8));
        __m256i hi = _mm256_loadu_si256((const __m256i*)(src + i * 8 + 32));

        /* For each byte position, shuffle both halves and combine to 8 bytes */
        #define DO_STREAM(B, SH) do { \
            __m256i slo = _mm256_shuffle_epi8(lo, SH); \
            __m256i shi = _mm256_shuffle_epi8(hi, SH); \
            uint32_t wlo = (uint32_t)_mm_cvtsi128_si32(_mm256_castsi256_si128(slo)); \
            uint32_t whi_lane = (uint32_t)_mm_cvtsi128_si32(_mm256_extracti128_si256(slo, 1)); \
            uint32_t who = (uint32_t)_mm_cvtsi128_si32(_mm256_castsi256_si128(shi)); \
            uint32_t who_lane = (uint32_t)_mm_cvtsi128_si32(_mm256_extracti128_si256(shi, 1)); \
            uint64_t val = (uint64_t)(uint16_t)wlo | ((uint64_t)(uint16_t)whi_lane << 16) | \
                          ((uint64_t)(uint16_t)who << 32) | ((uint64_t)(uint16_t)who_lane << 48); \
            memcpy(output + (B) * count + i, &val, sizeof(uint64_t)); \
        } while(0)

        DO_STREAM(0, sh0);
        DO_STREAM(1, sh1);
        DO_STREAM(2, sh2);
        DO_STREAM(3, sh3);
        DO_STREAM(4, sh4);
        DO_STREAM(5, sh5);
        DO_STREAM(6, sh6);
        DO_STREAM(7, sh7);

        #undef DO_STREAM
    }
#endif

    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            output[b * count + i] = src[i * 8 + b];
        }
    }
}

/**
 * Decode byte stream split doubles using AVX-512.
 * Processes 8 doubles (64 bytes) at a time using 512-bit operations.
 */
void carquet_avx512_byte_stream_split_decode_double(
    const uint8_t* data,
    int64_t count,
    double* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

#ifdef __AVX512VBMI__
    /* VBMI path: load 8 bytes from each of 8 streams into __m512i, single vpermb.
     * Input layout:
     *   bytes  0-7:  stream 0
     *   bytes  8-15: stream 1
     *   bytes 16-23: stream 2
     *   bytes 24-31: stream 3
     *   bytes 32-39: stream 4
     *   bytes 40-47: stream 5
     *   bytes 48-55: stream 6
     *   bytes 56-63: stream 7
     * Output: double[k] = {s0[k], s1[k], s2[k], s3[k], s4[k], s5[k], s6[k], s7[k]}
     *   output byte 8*k+j = input byte j*8+k */
    const __m512i perm = _mm512_set_epi32(
        0x3F372F27, 0x1F170F07,  /* double 7 */
        0x3E362E26, 0x1E160E06,  /* double 6 */
        0x3D352D25, 0x1D150D05,  /* double 5 */
        0x3C342C24, 0x1C140C04,  /* double 4 */
        0x3B332B23, 0x1B130B03,  /* double 3 */
        0x3A322A22, 0x1A120A02,  /* double 2 */
        0x39312921, 0x19110901,  /* double 1 */
        0x38302820, 0x18100800); /* double 0 */

    for (; i + 8 <= count; i += 8) {
        /* Pack 8 streams into one __m512i using 64-bit lane loads */
        __m128i p01 = _mm_unpacklo_epi64(
            _mm_loadl_epi64((const __m128i*)(data + 0 * count + i)),
            _mm_loadl_epi64((const __m128i*)(data + 1 * count + i)));
        __m128i p23 = _mm_unpacklo_epi64(
            _mm_loadl_epi64((const __m128i*)(data + 2 * count + i)),
            _mm_loadl_epi64((const __m128i*)(data + 3 * count + i)));
        __m128i p45 = _mm_unpacklo_epi64(
            _mm_loadl_epi64((const __m128i*)(data + 4 * count + i)),
            _mm_loadl_epi64((const __m128i*)(data + 5 * count + i)));
        __m128i p67 = _mm_unpacklo_epi64(
            _mm_loadl_epi64((const __m128i*)(data + 6 * count + i)),
            _mm_loadl_epi64((const __m128i*)(data + 7 * count + i)));

        __m512i combined = _mm512_castsi128_si512(p01);
        combined = _mm512_inserti32x4(combined, p23, 1);
        combined = _mm512_inserti32x4(combined, p45, 2);
        combined = _mm512_inserti32x4(combined, p67, 3);

        __m512i result = _mm512_permutexvar_epi8(perm, combined);
        _mm512_storeu_si512((__m512i*)(dst + i * 8), result);
    }
#else
    /* Non-VBMI: use 128-bit unpack cascade for 8 doubles, plus 512-bit stores */
    for (; i + 8 <= count; i += 8) {
        __m128i s0 = _mm_loadl_epi64((const __m128i*)(data + 0 * count + i));
        __m128i s1 = _mm_loadl_epi64((const __m128i*)(data + 1 * count + i));
        __m128i s2 = _mm_loadl_epi64((const __m128i*)(data + 2 * count + i));
        __m128i s3 = _mm_loadl_epi64((const __m128i*)(data + 3 * count + i));
        __m128i s4 = _mm_loadl_epi64((const __m128i*)(data + 4 * count + i));
        __m128i s5 = _mm_loadl_epi64((const __m128i*)(data + 5 * count + i));
        __m128i s6 = _mm_loadl_epi64((const __m128i*)(data + 6 * count + i));
        __m128i s7 = _mm_loadl_epi64((const __m128i*)(data + 7 * count + i));

        /* Stage 1: interleave bytes */
        __m128i u01 = _mm_unpacklo_epi8(s0, s1);
        __m128i u23 = _mm_unpacklo_epi8(s2, s3);
        __m128i u45 = _mm_unpacklo_epi8(s4, s5);
        __m128i u67 = _mm_unpacklo_epi8(s6, s7);

        /* Stage 2: interleave 16-bit words */
        __m128i v0 = _mm_unpacklo_epi16(u01, u23);
        __m128i v1 = _mm_unpackhi_epi16(u01, u23);
        __m128i v2 = _mm_unpacklo_epi16(u45, u67);
        __m128i v3 = _mm_unpackhi_epi16(u45, u67);

        /* Stage 3: interleave 32-bit dwords */
        __m128i d0 = _mm_unpacklo_epi32(v0, v2);
        __m128i d1 = _mm_unpackhi_epi32(v0, v2);
        __m128i d2 = _mm_unpacklo_epi32(v1, v3);
        __m128i d3 = _mm_unpackhi_epi32(v1, v3);

        /* Use 512-bit store: combine 4 x 128-bit results into one 512-bit write */
        __m512i out = _mm512_castsi128_si512(d0);
        out = _mm512_inserti32x4(out, d1, 1);
        out = _mm512_inserti32x4(out, d2, 2);
        out = _mm512_inserti32x4(out, d3, 3);
        _mm512_storeu_si512((__m512i*)(dst + i * 8), out);
    }
#endif

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

/* ============================================================================
 * Delta Decoding - AVX-512 Optimized (Prefix Sum)
 * ============================================================================
 */

/**
 * Apply prefix sum (cumulative sum) to int32 array using AVX-512.
 */
void carquet_avx512_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB.
     * Delta encoding relies on modular arithmetic — _mm512_add_epi32 is
     * already modular, so only the scalar accumulator needs fixing. */
    uint32_t sum = (uint32_t)initial;
    int64_t i = 0;

    /* AVX-512 prefix sum for 16 elements at a time */
    for (; i + 16 <= count; i += 16) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(values + i));

        /* Multi-step prefix sum within vector */
        /* Step 1: Add adjacent pairs */
        __m512i shifted1 = _mm512_maskz_alignr_epi32(0xFFFE, v, _mm512_setzero_si512(), 15);
        v = _mm512_add_epi32(v, shifted1);

        /* Step 2: Add elements 2 apart */
        __m512i shifted2 = _mm512_maskz_alignr_epi32(0xFFFC, v, _mm512_setzero_si512(), 14);
        v = _mm512_add_epi32(v, shifted2);

        /* Step 3: Add elements 4 apart */
        __m512i shifted4 = _mm512_maskz_alignr_epi32(0xFFF0, v, _mm512_setzero_si512(), 12);
        v = _mm512_add_epi32(v, shifted4);

        /* Step 4: Add elements 8 apart */
        __m512i shifted8 = _mm512_maskz_alignr_epi32(0xFF00, v, _mm512_setzero_si512(), 8);
        v = _mm512_add_epi32(v, shifted8);

        /* Add running sum */
        __m512i sums = _mm512_set1_epi32((int32_t)sum);
        v = _mm512_add_epi32(v, sums);
        _mm512_storeu_si512((__m512i*)(values + i), v);

        /* Update running sum to last element */
        sum = (uint32_t)values[i + 15];
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint32_t)values[i];
        values[i] = (int32_t)sum;
    }
}

/**
 * Apply prefix sum to int64 array using AVX-512.
 */
void carquet_avx512_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB. */
    uint64_t sum = (uint64_t)initial;
    int64_t i = 0;

    /* AVX-512 prefix sum for 8 elements at a time */
    for (; i + 8 <= count; i += 8) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(values + i));

        /* Multi-step prefix sum */
        __m512i shifted1 = _mm512_maskz_alignr_epi64(0xFE, v, _mm512_setzero_si512(), 7);
        v = _mm512_add_epi64(v, shifted1);

        __m512i shifted2 = _mm512_maskz_alignr_epi64(0xFC, v, _mm512_setzero_si512(), 6);
        v = _mm512_add_epi64(v, shifted2);

        __m512i shifted4 = _mm512_maskz_alignr_epi64(0xF0, v, _mm512_setzero_si512(), 4);
        v = _mm512_add_epi64(v, shifted4);

        /* Add running sum */
        __m512i sums = _mm512_set1_epi64((int64_t)sum);
        v = _mm512_add_epi64(v, sums);
        _mm512_storeu_si512((__m512i*)(values + i), v);

        /* Update running sum */
        sum = (uint64_t)values[i + 7];
    }

    /* Handle remaining values */
    for (; i < count; i++) {
        sum += (uint64_t)values[i];
        values[i] = (int64_t)sum;
    }
}

/* ============================================================================
 * Dictionary Gather - AVX-512 Optimized
 * ============================================================================
 */

/**
 * Gather int32 values from dictionary using AVX-512 gather instructions.
 */
void carquet_avx512_gather_i32(const int32_t* dict, const uint32_t* indices,
                               int64_t count, int32_t* output) {
    int64_t i = 0;

    /* Process 16 at a time using AVX-512 gather */
    for (; i + 16 <= count; i += 16) {
        __m512i idx = _mm512_loadu_si512((const __m512i*)(indices + i));
        __m512i result = _mm512_i32gather_epi32(idx, dict, 4);
        _mm512_storeu_si512((__m512i*)(output + i), result);
    }

    /* Handle remaining with AVX2 */
    for (; i + 8 <= count; i += 8) {
        __m256i idx = _mm256_loadu_si256((const __m256i*)(indices + i));
        __m256i result = _mm256_i32gather_epi32(dict, idx, 4);
        _mm256_storeu_si256((__m256i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather int64 values from dictionary using AVX-512 gather instructions.
 */
void carquet_avx512_gather_i64(const int64_t* dict, const uint32_t* indices,
                               int64_t count, int64_t* output) {
    int64_t i = 0;

    /* Process 8 at a time using AVX-512 gather */
    for (; i + 8 <= count; i += 8) {
        __m256i idx = _mm256_loadu_si256((const __m256i*)(indices + i));
        __m512i result = _mm512_i32gather_epi64(idx, dict, 8);
        _mm512_storeu_si512((__m512i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        output[i] = cq_loadu(dict + (indices[i]));
    }
}

/**
 * Gather float values from dictionary using AVX-512 gather instructions.
 * Note: float and int32 are both 4 bytes, so we reuse gather_i32 via cast.
 */
void carquet_avx512_gather_float(const float* dict, const uint32_t* indices,
                                  int64_t count, float* output) {
    /* Data movement doesn't care about type - reuse int32 implementation */
    carquet_avx512_gather_i32((const int32_t*)dict, indices, count, (int32_t*)output);
}

/**
 * Gather double values from dictionary using AVX-512 gather instructions.
 * Note: double and int64 are both 8 bytes, so we reuse gather_i64 via cast.
 */
void carquet_avx512_gather_double(const double* dict, const uint32_t* indices,
                                   int64_t count, double* output) {
    /* Data movement doesn't care about type - reuse int64 implementation */
    carquet_avx512_gather_i64((const int64_t*)dict, indices, count, (int64_t*)output);
}

bool carquet_avx512_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                        const uint32_t* indices, int64_t count,
                                        int32_t* output) {
    int64_t i = 0;
    __m512i limit = _mm512_set1_epi32(dict_count);

    for (; i + 16 <= count; i += 16) {
        __m512i idx = _mm512_loadu_si512((const void*)(indices + i));
        __mmask16 valid = _mm512_cmp_epu32_mask(idx, limit, _MM_CMPINT_LT);
        if (valid != 0xFFFFu) {
            return false;
        }
        __m512i result = _mm512_i32gather_epi32(idx, dict, 4);
        _mm512_storeu_si512((void*)(output + i), result);
    }

    for (; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            return false;
        }
        output[i] = cq_loadu(dict + (idx));
    }

    return true;
}

bool carquet_avx512_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                        const uint32_t* indices, int64_t count,
                                        int64_t* output) {
    int64_t i = 0;
    __m256i limit = _mm256_set1_epi32(dict_count);

    for (; i + 8 <= count; i += 8) {
        __m256i idx = _mm256_loadu_si256((const __m256i*)(indices + i));
        __mmask8 valid = _mm256_cmp_epu32_mask(idx, limit, _MM_CMPINT_LT);
        if (valid != 0xFFu) {
            return false;
        }
        __m512i result = _mm512_i32gather_epi64(idx, dict, 8);
        _mm512_storeu_si512((void*)(output + i), result);
    }

    for (; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            return false;
        }
        output[i] = cq_loadu(dict + (idx));
    }

    return true;
}

bool carquet_avx512_checked_gather_float(const float* dict, int32_t dict_count,
                                          const uint32_t* indices, int64_t count,
                                          float* output) {
    return carquet_avx512_checked_gather_i32(
        (const int32_t*)dict, dict_count, indices, count, (int32_t*)output);
}

bool carquet_avx512_checked_gather_double(const double* dict, int32_t dict_count,
                                           const uint32_t* indices, int64_t count,
                                           double* output) {
    return carquet_avx512_checked_gather_i64(
        (const int64_t*)dict, dict_count, indices, count, (int64_t*)output);
}

/* ============================================================================
 * Memcpy/Memset - AVX-512 Optimized
 * ============================================================================
 */

/**
 * Fast memset for large buffers using AVX-512.
 */
void carquet_avx512_memset(void* dest, uint8_t value, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    __m512i v = _mm512_set1_epi8((char)value);

    while (n >= 256) {
        _mm512_storeu_si512((__m512i*)(d + 0), v);
        _mm512_storeu_si512((__m512i*)(d + 64), v);
        _mm512_storeu_si512((__m512i*)(d + 128), v);
        _mm512_storeu_si512((__m512i*)(d + 192), v);
        d += 256;
        n -= 256;
    }

    while (n >= 64) {
        _mm512_storeu_si512((__m512i*)d, v);
        d += 64;
        n -= 64;
    }

    /* Handle tail with AVX2/SSE */
    __m256i v256 = _mm256_set1_epi8((char)value);
    while (n >= 32) {
        _mm256_storeu_si256((__m256i*)d, v256);
        d += 32;
        n -= 32;
    }

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
 * Fast memcpy for large buffers using AVX-512.
 */
void carquet_avx512_memcpy(void* dest, const void* src, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    const uint8_t* s = (const uint8_t*)src;

    while (n >= 256) {
        __m512i v0 = _mm512_loadu_si512((const __m512i*)(s + 0));
        __m512i v1 = _mm512_loadu_si512((const __m512i*)(s + 64));
        __m512i v2 = _mm512_loadu_si512((const __m512i*)(s + 128));
        __m512i v3 = _mm512_loadu_si512((const __m512i*)(s + 192));
        _mm512_storeu_si512((__m512i*)(d + 0), v0);
        _mm512_storeu_si512((__m512i*)(d + 64), v1);
        _mm512_storeu_si512((__m512i*)(d + 128), v2);
        _mm512_storeu_si512((__m512i*)(d + 192), v3);
        d += 256;
        s += 256;
        n -= 256;
    }

    while (n >= 64) {
        _mm512_storeu_si512((__m512i*)d, _mm512_loadu_si512((const __m512i*)s));
        d += 64;
        s += 64;
        n -= 64;
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

void carquet_avx512_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset) {
    if (offset >= 64) {
        while (len >= 64) {
            _mm512_storeu_si512((void*)dst, _mm512_loadu_si512((const void*)src));
            dst += 64;
            src += 64;
            len -= 64;
        }
    } else if (offset == 1) {
        __m512i v = _mm512_set1_epi8((char)*src);
        while (len >= 64) {
            _mm512_storeu_si512((void*)dst, v);
            dst += 64;
            len -= 64;
        }
    } else if (offset == 4) {
        uint32_t pattern;
        memcpy(&pattern, src, sizeof(pattern));
        __m512i v = _mm512_set1_epi32((int32_t)pattern);
        while (len >= 64) {
            _mm512_storeu_si512((void*)dst, v);
            dst += 64;
            len -= 64;
        }
    } else if (offset == 8) {
        uint64_t pattern;
        memcpy(&pattern, src, sizeof(pattern));
        __m512i v = _mm512_set1_epi64((long long)pattern);
        while (len >= 64) {
            _mm512_storeu_si512((void*)dst, v);
            dst += 64;
            len -= 64;
        }
    }

    while (len > 0) {
        *dst++ = *src++;
        len--;
    }
}

size_t carquet_avx512_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit) {
    const uint8_t* start = p;

    while (p + 64 <= limit) {
        __m512i a = _mm512_loadu_si512((const void*)p);
        __m512i b = _mm512_loadu_si512((const void*)match);
        __mmask64 mask = _mm512_cmpeq_epi8_mask(a, b);
        if (mask != ~0ULL) {
            return (size_t)(p - start) + (size_t)portable_ctz64(~mask);
        }
        p += 64;
        match += 64;
    }

    while (p < limit && *p == *match) {
        p++;
        match++;
    }

    return (size_t)(p - start);
}

/* ============================================================================
 * Boolean Operations - AVX-512 Optimized
 * ============================================================================
 */

/**
 * Unpack boolean values from packed bits to byte array using AVX-512.
 */
void carquet_avx512_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    /* Process 64 bools (8 bytes) at a time using AVX-512 mask */
    for (; i + 64 <= count; i += 64) {
        int byte_idx = (int)(i / 8);
        uint64_t packed;
        memcpy(&packed, input + byte_idx, 8);

        /* Convert to mask and create result with maskz_set1 (1 where set, 0 otherwise) */
        __m512i result = _mm512_maskz_set1_epi8((__mmask64)packed, 1);

        _mm512_storeu_si512((__m512i*)(output + i), result);
    }

    /* Handle remaining */
    for (; i < count; i++) {
        int byte_idx = (int)(i / 8);
        int bit_idx = (int)(i % 8);
        output[i] = (input[byte_idx] >> bit_idx) & 1;
    }
}

/**
 * Pack boolean values from byte array to packed bits using AVX-512.
 */
void carquet_avx512_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    /* Process 64 bools at a time */
    for (; i + 64 <= count; i += 64) {
        __m512i bools = _mm512_loadu_si512((const __m512i*)(input + i));

        /* Use test_epi8_mask: bit is set if (a & b) != 0, i.e., if bool is non-zero */
        __mmask64 mask = _mm512_test_epi8_mask(bools, bools);

        /* Store mask as 8 bytes */
        uint64_t packed = (uint64_t)mask;
        memcpy(output + i / 8, &packed, 8);
    }

    /* Handle remaining elements with masked load */
    if (i < count) {
        int64_t remaining = count - i;
        /* Create mask for remaining elements: set bits 0..(remaining-1) */
        __mmask64 load_mask = (remaining >= 64) ? ~0ULL : ((1ULL << remaining) - 1);

        /* Masked load zeros out elements beyond the mask */
        __m512i bools = _mm512_maskz_loadu_epi8(load_mask, input + i);

        /* Test for non-zero values */
        __mmask64 result_mask = _mm512_test_epi8_mask(bools, bools);

        /* Write only the bytes we need */
        int64_t bytes_to_write = (remaining + 7) / 8;
        uint64_t packed = (uint64_t)result_mask;
        memcpy(output + i / 8, &packed, (size_t)bytes_to_write);
    }
}

/* ============================================================================
 * Run Detection - AVX-512 Optimized
 * ============================================================================
 */

/**
 * Find the length of a run of repeated int32 values.
 */
int64_t carquet_avx512_find_run_length_i32(const int32_t* values, int64_t count) {
    if (count == 0) return 0;

    int32_t first = values[0];
    __m512i target = _mm512_set1_epi32(first);
    int64_t i = 0;

    /* Check 16 at a time */
    for (; i + 16 <= count; i += 16) {
        __m512i v = _mm512_loadu_si512((const __m512i*)(values + i));
        __mmask16 cmp = _mm512_cmpeq_epi32_mask(v, target);

        if (cmp != 0xFFFF) {  /* Not all equal */
            /* Find first mismatch using trailing zeros */
            int tz = portable_ctz(~cmp);
            return i + tz;
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

int64_t carquet_avx512_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    int64_t non_null_count = 0;
    int64_t i = 0;
    __m512i max_vec = _mm512_set1_epi16(max_def_level);

    for (; i + 32 <= count; i += 32) {
        __m512i levels = _mm512_loadu_si512((const void*)(def_levels + i));
        __mmask32 mask = _mm512_cmpeq_epi16_mask(levels, max_vec);
        non_null_count += portable_popcount((unsigned int)mask);
    }

    for (; i < count; i++) {
        if (def_levels[i] == max_def_level) {
            non_null_count++;
        }
    }

    return non_null_count;
}

void carquet_avx512_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                       int16_t max_def_level, uint8_t* null_bitmap) {
    int64_t i = 0;
    int64_t byte_index = 0;
    __m512i max_vec = _mm512_set1_epi16(max_def_level);

    for (; i + 32 <= count; i += 32, byte_index += 4) {
        __m512i levels = _mm512_loadu_si512((const void*)(def_levels + i));
        __mmask32 mask = _mm512_cmp_epi16_mask(levels, max_vec, _MM_CMPINT_EQ);
        uint32_t bits = (uint32_t)mask;
        memcpy(null_bitmap + byte_index, &bits, sizeof(bits));
    }

    for (; i < count; byte_index++) {
        uint8_t bits = 0;
        for (int j = 0; j < 8 && i < count; j++, i++) {
            if (def_levels[i] == max_def_level) {
                bits |= (uint8_t)(1u << j);
            }
        }
        null_bitmap[byte_index] = bits;
    }
}

void carquet_avx512_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    int64_t i = 0;
    __m512i val_vec = _mm512_set1_epi16(value);

    for (; i + 32 <= count; i += 32) {
        _mm512_storeu_si512((void*)(def_levels + i), val_vec);
    }
    for (; i < count; i++) {
        def_levels[i] = value;
    }
}

void carquet_avx512_minmax_i32(const int32_t* values, int64_t count,
                                int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    __m512i min_vec = _mm512_set1_epi32(min_v);
    __m512i max_vec = _mm512_set1_epi32(max_v);
    int64_t i = 1;

    for (; i + 16 <= count; i += 16) {
        __m512i v = _mm512_loadu_si512((const void*)(values + i));
        min_vec = _mm512_min_epi32(min_vec, v);
        max_vec = _mm512_max_epi32(max_vec, v);
    }

    int32_t tmp_min[16];
    int32_t tmp_max[16];
    _mm512_storeu_si512((void*)tmp_min, min_vec);
    _mm512_storeu_si512((void*)tmp_max, max_vec);
    for (int j = 0; j < 16; j++) {
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

void carquet_avx512_minmax_i64(const int64_t* values, int64_t count,
                                int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    __m512i min_vec = _mm512_set1_epi64(min_v);
    __m512i max_vec = _mm512_set1_epi64(max_v);
    int64_t i = 1;

    for (; i + 8 <= count; i += 8) {
        __m512i v = _mm512_loadu_si512((const void*)(values + i));
        __mmask8 lt = _mm512_cmpgt_epi64_mask(min_vec, v);
        __mmask8 gt = _mm512_cmpgt_epi64_mask(v, max_vec);
        min_vec = _mm512_mask_mov_epi64(min_vec, lt, v);
        max_vec = _mm512_mask_mov_epi64(max_vec, gt, v);
    }

    int64_t tmp_min[8];
    int64_t tmp_max[8];
    _mm512_storeu_si512((void*)tmp_min, min_vec);
    _mm512_storeu_si512((void*)tmp_max, max_vec);
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

void carquet_avx512_minmax_float(const float* values, int64_t count,
                                  float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m512 min_vec = _mm512_set1_ps(min_v);
    __m512 max_vec = _mm512_set1_ps(max_v);
    int64_t i = 1;

    for (; i + 16 <= count; i += 16) {
        __m512 v = _mm512_loadu_ps(values + i);
        __mmask16 lt = _mm512_cmplt_ps_mask(v, min_vec);
        __mmask16 gt = _mm512_cmp_ps_mask(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm512_mask_mov_ps(min_vec, lt, v);
        max_vec = _mm512_mask_mov_ps(max_vec, gt, v);
    }

    float tmp_min[16];
    float tmp_max[16];
    _mm512_storeu_ps(tmp_min, min_vec);
    _mm512_storeu_ps(tmp_max, max_vec);
    for (int j = 0; j < 16; j++) {
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

void carquet_avx512_minmax_double(const double* values, int64_t count,
                                   double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m512d min_vec = _mm512_set1_pd(min_v);
    __m512d max_vec = _mm512_set1_pd(max_v);
    int64_t i = 1;

    for (; i + 8 <= count; i += 8) {
        __m512d v = _mm512_loadu_pd(values + i);
        __mmask8 lt = _mm512_cmplt_pd_mask(v, min_vec);
        __mmask8 gt = _mm512_cmp_pd_mask(v, max_vec, _CMP_GT_OQ);
        min_vec = _mm512_mask_mov_pd(min_vec, lt, v);
        max_vec = _mm512_mask_mov_pd(max_vec, gt, v);
    }

    double tmp_min[8];
    double tmp_max[8];
    _mm512_storeu_pd(tmp_min, min_vec);
    _mm512_storeu_pd(tmp_max, max_vec);
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

/* ============================================================================
 * Conflict Detection - AVX-512 Specific
 * ============================================================================
 */

#ifdef __AVX512CD__

/**
 * Detect conflicts in indices for scatter operations.
 * Returns a mask where bit i is set if indices[i] conflicts with any earlier index.
 */
__mmask16 carquet_avx512_detect_conflicts_i32(const uint32_t* indices) {
    __m512i idx = _mm512_loadu_si512((const __m512i*)indices);
    __m512i conflicts = _mm512_conflict_epi32(idx);

    /* Non-zero conflict value means there's a conflict */
    return _mm512_cmpneq_epi32_mask(conflicts, _mm512_setzero_si512());
}

#endif /* __AVX512CD__ */

#endif /* __AVX512F__ */
#endif /* x86_64 */
