/**
 * @file avx_ops.c
 * @brief AVX optimized operations for x86-64 processors
 *
 * Provides a dedicated AVX tier for hosts that support AVX but not AVX2.
 * The main win here is wider byte-stream-split processing for floating-point
 * columns while reusing 128-bit shuffle operations within the AVX lanes.
 */

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#if defined(__AVX__) || (defined(_MSC_VER) && defined(__AVX__))

#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <immintrin.h>

/* min/max accumulators below put the loaded vector FIRST: MINPS/MAXPS return
 * the second operand when either is NaN, so a NaN value leaves the
 * accumulator untouched (NaN-skipping, like the AVX2 compare+blend path).
 * With (acc, v) a NaN would replace the accumulator and silently discard
 * every value seen so far in that lane. */

void carquet_avx_minmax_float(const float* values, int64_t count,
                               float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m256 min_vec = _mm256_set1_ps(min_v);
    __m256 max_vec = _mm256_set1_ps(max_v);
    int64_t i = 1;

    for (; i + 8 <= count; i += 8) {
        __m256 v = _mm256_loadu_ps(values + i);
        min_vec = _mm256_min_ps(v, min_vec);
        max_vec = _mm256_max_ps(v, max_vec);
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

void carquet_avx_minmax_double(const double* values, int64_t count,
                                double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m256d min_vec = _mm256_set1_pd(min_v);
    __m256d max_vec = _mm256_set1_pd(max_v);
    int64_t i = 1;

    for (; i + 4 <= count; i += 4) {
        __m256d v = _mm256_loadu_pd(values + i);
        min_vec = _mm256_min_pd(v, min_vec);
        max_vec = _mm256_max_pd(v, max_vec);
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

void carquet_avx_copy_minmax_float(const float* values, int64_t count, float* output,
                                    float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    __m256 min_vec = _mm256_set1_ps(min_v);
    __m256 max_vec = _mm256_set1_ps(max_v);
    int64_t i = 0;

    for (; i + 8 <= count; i += 8) {
        __m256 v = _mm256_loadu_ps(values + i);
        _mm256_storeu_ps(output + i, v);
        min_vec = _mm256_min_ps(v, min_vec);
        max_vec = _mm256_max_ps(v, max_vec);
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

void carquet_avx_copy_minmax_double(const double* values, int64_t count, double* output,
                                     double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    __m256d min_vec = _mm256_set1_pd(min_v);
    __m256d max_vec = _mm256_set1_pd(max_v);
    int64_t i = 0;

    for (; i + 4 <= count; i += 4) {
        __m256d v = _mm256_loadu_pd(values + i);
        _mm256_storeu_pd(output + i, v);
        min_vec = _mm256_min_pd(v, min_vec);
        max_vec = _mm256_max_pd(v, max_vec);
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

void carquet_avx_byte_stream_split_encode_float(
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

    for (; i + 8 <= count; i += 8) {
        __m256 v = _mm256_loadu_ps(values + i);
        __m128 lo_ps = _mm256_castps256_ps128(v);
        __m128 hi_ps = _mm256_extractf128_ps(v, 1);
        __m128i lo = _mm_castps_si128(lo_ps);
        __m128i hi = _mm_castps_si128(hi_ps);

        uint32_t t0_lo = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(lo, s0), 0);
        uint32_t t1_lo = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(lo, s1), 0);
        uint32_t t2_lo = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(lo, s2), 0);
        uint32_t t3_lo = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(lo, s3), 0);
        uint32_t t0_hi = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(hi, s0), 0);
        uint32_t t1_hi = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(hi, s1), 0);
        uint32_t t2_hi = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(hi, s2), 0);
        uint32_t t3_hi = (uint32_t)_mm_extract_epi32(_mm_shuffle_epi8(hi, s3), 0);

        memcpy(output + 0 * count + i, &t0_lo, sizeof(t0_lo));
        memcpy(output + 0 * count + i + 4, &t0_hi, sizeof(t0_hi));
        memcpy(output + 1 * count + i, &t1_lo, sizeof(t1_lo));
        memcpy(output + 1 * count + i + 4, &t1_hi, sizeof(t1_hi));
        memcpy(output + 2 * count + i, &t2_lo, sizeof(t2_lo));
        memcpy(output + 2 * count + i + 4, &t2_hi, sizeof(t2_hi));
        memcpy(output + 3 * count + i, &t3_lo, sizeof(t3_lo));
        memcpy(output + 3 * count + i + 4, &t3_hi, sizeof(t3_hi));
    }

    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            output[b * count + i] = src[i * 4 + b];
        }
    }
}

void carquet_avx_byte_stream_split_decode_float(
    const uint8_t* data,
    int64_t count,
    float* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    for (; i + 8 <= count; i += 8) {
        uint32_t b0_lo, b0_hi, b1_lo, b1_hi, b2_lo, b2_hi, b3_lo, b3_hi;
        memcpy(&b0_lo, data + 0 * count + i, sizeof(b0_lo));
        memcpy(&b0_hi, data + 0 * count + i + 4, sizeof(b0_hi));
        memcpy(&b1_lo, data + 1 * count + i, sizeof(b1_lo));
        memcpy(&b1_hi, data + 1 * count + i + 4, sizeof(b1_hi));
        memcpy(&b2_lo, data + 2 * count + i, sizeof(b2_lo));
        memcpy(&b2_hi, data + 2 * count + i + 4, sizeof(b2_hi));
        memcpy(&b3_lo, data + 3 * count + i, sizeof(b3_lo));
        memcpy(&b3_hi, data + 3 * count + i + 4, sizeof(b3_hi));

        __m128i b0l = _mm_cvtsi32_si128((int)b0_lo);
        __m128i b1l = _mm_cvtsi32_si128((int)b1_lo);
        __m128i b2l = _mm_cvtsi32_si128((int)b2_lo);
        __m128i b3l = _mm_cvtsi32_si128((int)b3_lo);
        __m128i b0h = _mm_cvtsi32_si128((int)b0_hi);
        __m128i b1h = _mm_cvtsi32_si128((int)b1_hi);
        __m128i b2h = _mm_cvtsi32_si128((int)b2_hi);
        __m128i b3h = _mm_cvtsi32_si128((int)b3_hi);

        __m128i lo01 = _mm_unpacklo_epi8(b0l, b1l);
        __m128i lo23 = _mm_unpacklo_epi8(b2l, b3l);
        __m128i hi01 = _mm_unpacklo_epi8(b0h, b1h);
        __m128i hi23 = _mm_unpacklo_epi8(b2h, b3h);
        __m128i result_lo = _mm_unpacklo_epi16(lo01, lo23);
        __m128i result_hi = _mm_unpacklo_epi16(hi01, hi23);

        _mm_storeu_si128((__m128i*)(dst + i * 4), result_lo);
        _mm_storeu_si128((__m128i*)(dst + i * 4 + 16), result_hi);
    }

    for (; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            dst[i * 4 + b] = data[b * count + i];
        }
    }
}

void carquet_avx_byte_stream_split_encode_double(
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

    for (; i + 4 <= count; i += 4) {
        __m256d v = _mm256_loadu_pd(values + i);
        __m128i lo = _mm_castpd_si128(_mm256_castpd256_pd128(v));
        __m128i hi = _mm_castpd_si128(_mm256_extractf128_pd(v, 1));

        uint16_t t0_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s0), 0);
        uint16_t t1_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s1), 0);
        uint16_t t2_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s2), 0);
        uint16_t t3_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s3), 0);
        uint16_t t4_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s4), 0);
        uint16_t t5_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s5), 0);
        uint16_t t6_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s6), 0);
        uint16_t t7_lo = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(lo, s7), 0);
        uint16_t t0_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s0), 0);
        uint16_t t1_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s1), 0);
        uint16_t t2_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s2), 0);
        uint16_t t3_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s3), 0);
        uint16_t t4_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s4), 0);
        uint16_t t5_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s5), 0);
        uint16_t t6_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s6), 0);
        uint16_t t7_hi = (uint16_t)_mm_extract_epi16(_mm_shuffle_epi8(hi, s7), 0);

        uint32_t t0 = (uint32_t)t0_lo | ((uint32_t)t0_hi << 16);
        uint32_t t1 = (uint32_t)t1_lo | ((uint32_t)t1_hi << 16);
        uint32_t t2 = (uint32_t)t2_lo | ((uint32_t)t2_hi << 16);
        uint32_t t3 = (uint32_t)t3_lo | ((uint32_t)t3_hi << 16);
        uint32_t t4 = (uint32_t)t4_lo | ((uint32_t)t4_hi << 16);
        uint32_t t5 = (uint32_t)t5_lo | ((uint32_t)t5_hi << 16);
        uint32_t t6 = (uint32_t)t6_lo | ((uint32_t)t6_hi << 16);
        uint32_t t7 = (uint32_t)t7_lo | ((uint32_t)t7_hi << 16);

        memcpy(output + 0 * count + i, &t0, sizeof(t0));
        memcpy(output + 1 * count + i, &t1, sizeof(t1));
        memcpy(output + 2 * count + i, &t2, sizeof(t2));
        memcpy(output + 3 * count + i, &t3, sizeof(t3));
        memcpy(output + 4 * count + i, &t4, sizeof(t4));
        memcpy(output + 5 * count + i, &t5, sizeof(t5));
        memcpy(output + 6 * count + i, &t6, sizeof(t6));
        memcpy(output + 7 * count + i, &t7, sizeof(t7));
    }

    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            output[b * count + i] = src[i * 8 + b];
        }
    }
}

void carquet_avx_byte_stream_split_decode_double(
    const uint8_t* data,
    int64_t count,
    double* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

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
        __m128i lo = _mm_unpacklo_epi32(v0, v1);
        __m128i hi = _mm_unpackhi_epi32(v0, v1);

        _mm_storeu_si128((__m128i*)(dst + i * 8), lo);
        _mm_storeu_si128((__m128i*)(dst + i * 8 + 16), hi);
    }

    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

#endif /* __AVX__ */
#endif /* x86 */
