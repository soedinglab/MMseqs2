/**
 * @file dispatch.c
 * @brief SIMD function dispatch
 *
 * This file provides runtime dispatch for SIMD-optimized functions based on
 * detected CPU features. Functions are selected at initialization time and
 * stored in function pointer tables for efficient runtime access.
 */

#include <carquet/carquet.h>
#include "core/bitpack.h"
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* ============================================================================
 * Function Pointer Types
 * ============================================================================
 */

typedef void (*prefix_sum_i32_fn)(int32_t* values, int64_t count, int32_t initial);
typedef void (*prefix_sum_i64_fn)(int64_t* values, int64_t count, int64_t initial);

typedef void (*gather_i32_fn)(const int32_t* dict, const uint32_t* indices,
                               int64_t count, int32_t* output);
typedef void (*gather_i64_fn)(const int64_t* dict, const uint32_t* indices,
                               int64_t count, int64_t* output);
typedef void (*gather_float_fn)(const float* dict, const uint32_t* indices,
                                 int64_t count, float* output);
typedef void (*gather_double_fn)(const double* dict, const uint32_t* indices,
                                  int64_t count, double* output);
typedef bool (*checked_gather_i32_fn)(const int32_t* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       int32_t* output);
typedef bool (*checked_gather_i64_fn)(const int64_t* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       int64_t* output);
typedef bool (*checked_gather_float_fn)(const float* dict, int32_t dict_count,
                                         const uint32_t* indices, int64_t count,
                                         float* output);
typedef bool (*checked_gather_double_fn)(const double* dict, int32_t dict_count,
                                          const uint32_t* indices, int64_t count,
                                          double* output);

typedef void (*byte_split_encode_float_fn)(const float* values, int64_t count,
                                            uint8_t* output);
typedef void (*byte_split_decode_float_fn)(const uint8_t* data, int64_t count,
                                            float* values);
typedef void (*byte_split_encode_double_fn)(const double* values, int64_t count,
                                             uint8_t* output);
typedef void (*byte_split_decode_double_fn)(const uint8_t* data, int64_t count,
                                             double* values);

typedef void (*memset_fn)(void* dest, uint8_t value, size_t n);
typedef void (*memcpy_fn)(void* dest, const void* src, size_t n);

typedef void (*unpack_bools_fn)(const uint8_t* input, uint8_t* output, int64_t count);
typedef void (*pack_bools_fn)(const uint8_t* input, uint8_t* output, int64_t count);
typedef void (*bitunpack8_u32_fn)(const uint8_t* input, uint32_t* values);

typedef int64_t (*find_run_length_i32_fn)(const int32_t* values, int64_t count);

typedef void (*match_copy_fn)(uint8_t* dst, const uint8_t* src, size_t len, size_t offset);
typedef size_t (*match_length_fn)(const uint8_t* p, const uint8_t* match, const uint8_t* limit);

typedef int64_t (*count_non_nulls_fn)(const int16_t* def_levels, int64_t count, int16_t max_def_level);
typedef void (*build_null_bitmap_fn)(const int16_t* def_levels, int64_t count,
                                      int16_t max_def_level, uint8_t* null_bitmap);
typedef void (*fill_def_levels_fn)(int16_t* def_levels, int64_t count, int16_t value);
typedef void (*minmax_i32_fn)(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
typedef void (*minmax_i64_fn)(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
typedef void (*minmax_float_fn)(const float* values, int64_t count, float* min_value, float* max_value);
typedef void (*minmax_double_fn)(const double* values, int64_t count, double* min_value, double* max_value);
typedef void (*copy_minmax_i32_fn)(const int32_t* values, int64_t count, int32_t* output,
                                    int32_t* min_value, int32_t* max_value);
typedef void (*copy_minmax_i64_fn)(const int64_t* values, int64_t count, int64_t* output,
                                    int64_t* min_value, int64_t* max_value);
typedef void (*copy_minmax_float_fn)(const float* values, int64_t count, float* output,
                                      float* min_value, float* max_value);
typedef void (*copy_minmax_double_fn)(const double* values, int64_t count, double* output,
                                       double* min_value, double* max_value);

/* ============================================================================
 * Scalar Fallback Implementations
 * ============================================================================
 */

/* Portable software prefetch */
#if defined(_MSC_VER)
#include <intrin.h>
#define CARQUET_PREFETCH(addr) _mm_prefetch((const char*)(addr), _MM_HINT_T1)
#elif defined(__GNUC__) || defined(__clang__)
#define CARQUET_PREFETCH(addr) __builtin_prefetch((addr), 0, 1)
#else
#define CARQUET_PREFETCH(addr) ((void)0)
#endif

static void scalar_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    uint32_t sum = (uint32_t)initial;
    for (int64_t i = 0; i < count; i++) {
        sum += (uint32_t)values[i];
        values[i] = (int32_t)sum;
    }
}

static void scalar_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    uint64_t sum = (uint64_t)initial;
    for (int64_t i = 0; i < count; i++) {
        sum += (uint64_t)values[i];
        values[i] = (int64_t)sum;
    }
}

/* Unaligned dictionary loads (portable; see header). */
#include "simd/simd_unaligned.h"

static void scalar_gather_i32(const int32_t* dict, const uint32_t* indices,
                               int64_t count, int32_t* output) {
    const int64_t prefetch_dist = 8;
    for (int64_t i = 0; i < count; i++) {
        if (i + prefetch_dist < count) {
            CARQUET_PREFETCH(&dict[indices[i + prefetch_dist]]);
        }
        output[i] = cq_load_i32u(dict + indices[i]);
    }
}

static void scalar_gather_i64(const int64_t* dict, const uint32_t* indices,
                               int64_t count, int64_t* output) {
    const int64_t prefetch_dist = 8;
    for (int64_t i = 0; i < count; i++) {
        if (i + prefetch_dist < count) {
            CARQUET_PREFETCH(&dict[indices[i + prefetch_dist]]);
        }
        output[i] = cq_load_i64u(dict + indices[i]);
    }
}

static void scalar_gather_float(const float* dict, const uint32_t* indices,
                                 int64_t count, float* output) {
    const int64_t prefetch_dist = 8;
    for (int64_t i = 0; i < count; i++) {
        if (i + prefetch_dist < count) {
            CARQUET_PREFETCH(&dict[indices[i + prefetch_dist]]);
        }
        output[i] = cq_load_f32u(dict + indices[i]);
    }
}

static void scalar_gather_double(const double* dict, const uint32_t* indices,
                                  int64_t count, double* output) {
    const int64_t prefetch_dist = 8;
    for (int64_t i = 0; i < count; i++) {
        if (i + prefetch_dist < count) {
            CARQUET_PREFETCH(&dict[indices[i + prefetch_dist]]);
        }
        output[i] = cq_load_f64u(dict + indices[i]);
    }
}

static bool scalar_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       int32_t* output) {
    for (int64_t i = 0; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            return false;
        }
        output[i] = cq_load_i32u(dict + idx);
    }
    return true;
}

static bool scalar_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       int64_t* output) {
    for (int64_t i = 0; i < count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            return false;
        }
        output[i] = cq_load_i64u(dict + idx);
    }
    return true;
}

static bool scalar_checked_gather_float(const float* dict, int32_t dict_count,
                                         const uint32_t* indices, int64_t count,
                                         float* output) {
    return scalar_checked_gather_i32((const int32_t*)dict, dict_count, indices,
                                     count, (int32_t*)output);
}

static bool scalar_checked_gather_double(const double* dict, int32_t dict_count,
                                          const uint32_t* indices, int64_t count,
                                          double* output) {
    return scalar_checked_gather_i64((const int64_t*)dict, dict_count, indices,
                                     count, (int64_t*)output);
}

static bool validate_gather_indices(const uint32_t* indices, int64_t count, int32_t dict_count) {
    uint32_t limit = (uint32_t)dict_count;
    for (int64_t i = 0; i < count; i++) {
        if (indices[i] >= limit) {
            return false;
        }
    }
    return true;
}

static void scalar_byte_split_encode_float(const float* values, int64_t count,
                                            uint8_t* output) {
    const uint8_t* src = (const uint8_t*)values;
    for (int64_t i = 0; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            output[b * count + i] = src[i * 4 + b];
        }
    }
}

static void scalar_byte_split_decode_float(const uint8_t* data, int64_t count,
                                            float* values) {
    uint8_t* dst = (uint8_t*)values;
    for (int64_t i = 0; i < count; i++) {
        for (int b = 0; b < 4; b++) {
            dst[i * 4 + b] = data[b * count + i];
        }
    }
}

static void scalar_byte_split_encode_double(const double* values, int64_t count,
                                             uint8_t* output) {
    const uint8_t* src = (const uint8_t*)values;
    for (int64_t i = 0; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            output[b * count + i] = src[i * 8 + b];
        }
    }
}

static void scalar_byte_split_decode_double(const uint8_t* data, int64_t count,
                                             double* values) {
    uint8_t* dst = (uint8_t*)values;
    for (int64_t i = 0; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

static void scalar_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        int byte_idx = (int)(i / 8);
        int bit_idx = (int)(i % 8);
        output[i] = (input[byte_idx] >> bit_idx) & 1;
    }
}

static void scalar_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    for (int64_t i = 0; i < count; i += 8) {
        uint8_t byte = 0;
        for (int64_t j = 0; j < 8 && i + j < count; j++) {
            if (input[i + j]) {
                byte |= (1 << j);
            }
        }
        output[i / 8] = byte;
    }
}

static int64_t scalar_find_run_length_i32(const int32_t* values, int64_t count) {
    if (count == 0) return 0;
    int32_t first = values[0];
    for (int64_t i = 1; i < count; i++) {
        if (values[i] != first) return i;
    }
    return count;
}

static void scalar_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset) {
    if (offset >= 8) {
        /* Non-overlapping: copy 8 bytes at a time */
        while (len >= 8) {
            memcpy(dst, src, 8);
            dst += 8;
            src += 8;
            len -= 8;
        }
        while (len > 0) {
            *dst++ = *src++;
            len--;
        }
    } else {
        /* Overlapping: byte by byte */
        while (len > 0) {
            *dst++ = *src++;
            len--;
        }
    }
}

static size_t scalar_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit) {
    const uint8_t* start = p;
    while (p < limit && *p == *match) {
        p++;
        match++;
    }
    return (size_t)(p - start);
}

static int64_t scalar_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    int64_t non_null_count = 0;
    for (int64_t i = 0; i < count; i++) {
        if (def_levels[i] == max_def_level) {
            non_null_count++;
        }
    }
    return non_null_count;
}

static void scalar_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                      int16_t max_def_level, uint8_t* null_bitmap) {
    int64_t full_bytes = count / 8;
    for (int64_t b = 0; b < full_bytes; b++) {
        uint8_t present_bits = 0;
        int64_t base = b * 8;
        if (def_levels[base + 0] == max_def_level) present_bits |= 0x01;
        if (def_levels[base + 1] == max_def_level) present_bits |= 0x02;
        if (def_levels[base + 2] == max_def_level) present_bits |= 0x04;
        if (def_levels[base + 3] == max_def_level) present_bits |= 0x08;
        if (def_levels[base + 4] == max_def_level) present_bits |= 0x10;
        if (def_levels[base + 5] == max_def_level) present_bits |= 0x20;
        if (def_levels[base + 6] == max_def_level) present_bits |= 0x40;
        if (def_levels[base + 7] == max_def_level) present_bits |= 0x80;
        null_bitmap[b] = present_bits;
    }
    for (int64_t j = full_bytes * 8; j < count; j++) {
        if (def_levels[j] == max_def_level) {
            null_bitmap[j / 8] |= (1 << (j % 8));
        }
    }
}

static void scalar_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    for (int64_t i = 0; i < count; i++) {
        def_levels[i] = value;
    }
}

static void scalar_minmax_i32(const int32_t* values, int64_t count,
                               int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    for (int64_t i = 1; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }
    *min_value = min_v;
    *max_value = max_v;
}

static void scalar_minmax_i64(const int64_t* values, int64_t count,
                               int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    for (int64_t i = 1; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }
    *min_value = min_v;
    *max_value = max_v;
}

static void scalar_minmax_float(const float* values, int64_t count,
                                 float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    for (int64_t i = 1; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }
    *min_value = min_v;
    *max_value = max_v;
}

static void scalar_minmax_double(const double* values, int64_t count,
                                  double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    for (int64_t i = 1; i < count; i++) {
        if (values[i] < min_v) min_v = values[i];
        if (values[i] > max_v) max_v = values[i];
    }
    *min_value = min_v;
    *max_value = max_v;
}

static void scalar_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                    int32_t* min_value, int32_t* max_value) {
    memcpy(output, values, (size_t)count * sizeof(int32_t));
    scalar_minmax_i32(values, count, min_value, max_value);
}

static void scalar_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                    int64_t* min_value, int64_t* max_value) {
    memcpy(output, values, (size_t)count * sizeof(int64_t));
    scalar_minmax_i64(values, count, min_value, max_value);
}

static void scalar_copy_minmax_float(const float* values, int64_t count, float* output,
                                      float* min_value, float* max_value) {
    memcpy(output, values, (size_t)count * sizeof(float));
    scalar_minmax_float(values, count, min_value, max_value);
}

static void scalar_copy_minmax_double(const double* values, int64_t count, double* output,
                                       double* min_value, double* max_value) {
    memcpy(output, values, (size_t)count * sizeof(double));
    scalar_minmax_double(values, count, min_value, max_value);
}

/* ============================================================================
 * External SIMD Function Declarations
 * ============================================================================
 */

/* Use CMake defines instead of compiler intrinsic macros, since dispatch.c
 * is not compiled with -msse4.2/-mavx2/-mavx512f flags */
#if defined(CARQUET_ARCH_X86)

#ifdef CARQUET_ENABLE_SSE
extern void carquet_sse_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial);
extern void carquet_sse_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial);
extern void carquet_sse_gather_i32(const int32_t* dict, const uint32_t* indices,
                                    int64_t count, int32_t* output);
extern void carquet_sse_gather_i64(const int64_t* dict, const uint32_t* indices,
                                    int64_t count, int64_t* output);
extern void carquet_sse_gather_float(const float* dict, const uint32_t* indices,
                                      int64_t count, float* output);
extern void carquet_sse_gather_double(const double* dict, const uint32_t* indices,
                                       int64_t count, double* output);
extern bool carquet_sse_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                            const uint32_t* indices, int64_t count,
                                            int32_t* output);
extern bool carquet_sse_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                            const uint32_t* indices, int64_t count,
                                            int64_t* output);
extern bool carquet_sse_checked_gather_float(const float* dict, int32_t dict_count,
                                              const uint32_t* indices, int64_t count,
                                              float* output);
extern bool carquet_sse_checked_gather_double(const double* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               double* output);
extern void carquet_sse_byte_stream_split_encode_float(const float* values, int64_t count,
                                                        uint8_t* output);
extern void carquet_sse_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                        float* values);
extern void carquet_sse_byte_stream_split_encode_double(const double* values, int64_t count,
                                                         uint8_t* output);
extern void carquet_sse_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                         double* values);
extern void carquet_sse_bitunpack8_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_2bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_3bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_5bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_6bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_7bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack8_16bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_bitunpack32_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_sse_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_sse_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_sse_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset);
extern size_t carquet_sse_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit);
extern int64_t carquet_sse_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level);
extern void carquet_sse_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                           int16_t max_def_level, uint8_t* null_bitmap);
extern void carquet_sse_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
extern void carquet_sse_minmax_i32(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
extern void carquet_sse_minmax_i64(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
extern void carquet_sse_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_sse_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern void carquet_sse_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                         int32_t* min_value, int32_t* max_value);
extern void carquet_sse_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                         int64_t* min_value, int64_t* max_value);
extern void carquet_sse_copy_minmax_float(const float* values, int64_t count, float* output,
                                           float* min_value, float* max_value);
extern void carquet_sse_copy_minmax_double(const double* values, int64_t count, double* output,
                                            double* min_value, double* max_value);
extern int64_t carquet_sse_find_run_length_i32(const int32_t* values, int64_t count);
#endif

#ifdef CARQUET_ENABLE_AVX
extern void carquet_avx_byte_stream_split_encode_float(const float* values, int64_t count,
                                                        uint8_t* output);
extern void carquet_avx_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                        float* values);
extern void carquet_avx_byte_stream_split_encode_double(const double* values, int64_t count,
                                                         uint8_t* output);
extern void carquet_avx_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                         double* values);
extern void carquet_avx_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_avx_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern void carquet_avx_copy_minmax_float(const float* values, int64_t count, float* output,
                                           float* min_value, float* max_value);
extern void carquet_avx_copy_minmax_double(const double* values, int64_t count, double* output,
                                            double* min_value, double* max_value);
#endif

#ifdef CARQUET_ENABLE_AVX2
extern void carquet_avx2_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial);
extern void carquet_avx2_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial);
extern void carquet_avx2_gather_i32(const int32_t* dict, const uint32_t* indices,
                                     int64_t count, int32_t* output);
extern void carquet_avx2_gather_i64(const int64_t* dict, const uint32_t* indices,
                                     int64_t count, int64_t* output);
extern void carquet_avx2_gather_float(const float* dict, const uint32_t* indices,
                                       int64_t count, float* output);
extern void carquet_avx2_gather_double(const double* dict, const uint32_t* indices,
                                        int64_t count, double* output);
extern bool carquet_avx2_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                             const uint32_t* indices, int64_t count,
                                             int32_t* output);
extern bool carquet_avx2_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                             const uint32_t* indices, int64_t count,
                                             int64_t* output);
extern bool carquet_avx2_checked_gather_float(const float* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               float* output);
extern bool carquet_avx2_checked_gather_double(const double* dict, int32_t dict_count,
                                                const uint32_t* indices, int64_t count,
                                                double* output);
extern void carquet_avx2_byte_stream_split_encode_float(const float* values, int64_t count,
                                                         uint8_t* output);
extern void carquet_avx2_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                         float* values);
extern void carquet_avx2_byte_stream_split_encode_double(const double* values, int64_t count,
                                                          uint8_t* output);
extern void carquet_avx2_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                          double* values);
extern void carquet_avx2_bitunpack8_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_2bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_3bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_5bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_6bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_7bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack8_16bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack16_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_bitunpack16_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx2_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_avx2_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_avx2_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset);
extern size_t carquet_avx2_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit);
extern int64_t carquet_avx2_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level);
extern void carquet_avx2_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                            int16_t max_def_level, uint8_t* null_bitmap);
extern void carquet_avx2_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
extern void carquet_avx2_minmax_i32(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
extern void carquet_avx2_minmax_i64(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
extern void carquet_avx2_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_avx2_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern void carquet_avx2_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                          int32_t* min_value, int32_t* max_value);
extern void carquet_avx2_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                          int64_t* min_value, int64_t* max_value);
extern void carquet_avx2_copy_minmax_float(const float* values, int64_t count, float* output,
                                            float* min_value, float* max_value);
extern void carquet_avx2_copy_minmax_double(const double* values, int64_t count, double* output,
                                             double* min_value, double* max_value);
extern int64_t carquet_avx2_find_run_length_i32(const int32_t* values, int64_t count);
#endif

#ifdef CARQUET_ENABLE_AVX512
extern void carquet_avx512_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial);
extern void carquet_avx512_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial);
extern void carquet_avx512_gather_i32(const int32_t* dict, const uint32_t* indices,
                                       int64_t count, int32_t* output);
extern void carquet_avx512_gather_i64(const int64_t* dict, const uint32_t* indices,
                                       int64_t count, int64_t* output);
extern void carquet_avx512_gather_float(const float* dict, const uint32_t* indices,
                                         int64_t count, float* output);
extern void carquet_avx512_gather_double(const double* dict, const uint32_t* indices,
                                          int64_t count, double* output);
extern bool carquet_avx512_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               int32_t* output);
extern bool carquet_avx512_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               int64_t* output);
extern bool carquet_avx512_checked_gather_float(const float* dict, int32_t dict_count,
                                                 const uint32_t* indices, int64_t count,
                                                 float* output);
extern bool carquet_avx512_checked_gather_double(const double* dict, int32_t dict_count,
                                                  const uint32_t* indices, int64_t count,
                                                  double* output);
extern void carquet_avx512_byte_stream_split_encode_float(const float* values, int64_t count,
                                                           uint8_t* output);
extern void carquet_avx512_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                           float* values);
extern void carquet_avx512_byte_stream_split_encode_double(const double* values, int64_t count,
                                                            uint8_t* output);
extern void carquet_avx512_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                            double* values);
extern void carquet_avx512_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_bitunpack8_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_bitunpack8_16bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_bitunpack32_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_bitunpack32_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_bitunpack16_16bit(const uint8_t* input, uint32_t* values);
extern void carquet_avx512_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_avx512_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_avx512_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset);
extern size_t carquet_avx512_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit);
extern int64_t carquet_avx512_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level);
extern void carquet_avx512_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                              int16_t max_def_level, uint8_t* null_bitmap);
extern void carquet_avx512_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
extern void carquet_avx512_minmax_i32(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
extern void carquet_avx512_minmax_i64(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
extern void carquet_avx512_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_avx512_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern int64_t carquet_avx512_find_run_length_i32(const int32_t* values, int64_t count);
#endif

#endif /* CARQUET_ARCH_X86 */

#if defined(CARQUET_ARCH_ARM)

/* NEON declarations - compiled when the ARM NEON backend is enabled. */
#if defined(CARQUET_ENABLE_NEON) && (defined(__ARM_NEON) || defined(__ARM_NEON__))
extern void carquet_neon_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial);
extern void carquet_neon_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial);
extern void carquet_neon_gather_i32(const int32_t* dict, const uint32_t* indices,
                                     int64_t count, int32_t* output);
extern void carquet_neon_gather_i64(const int64_t* dict, const uint32_t* indices,
                                     int64_t count, int64_t* output);
extern void carquet_neon_gather_float(const float* dict, const uint32_t* indices,
                                       int64_t count, float* output);
extern void carquet_neon_gather_double(const double* dict, const uint32_t* indices,
                                        int64_t count, double* output);
extern bool carquet_neon_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                             const uint32_t* indices, int64_t count,
                                             int32_t* output);
extern bool carquet_neon_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                             const uint32_t* indices, int64_t count,
                                             int64_t* output);
extern bool carquet_neon_checked_gather_float(const float* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               float* output);
extern bool carquet_neon_checked_gather_double(const double* dict, int32_t dict_count,
                                                const uint32_t* indices, int64_t count,
                                                double* output);
extern void carquet_neon_byte_stream_split_encode_float(const float* values, int64_t count,
                                                         uint8_t* output);
extern void carquet_neon_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                         float* values);
extern void carquet_neon_byte_stream_split_encode_double(const double* values, int64_t count,
                                                          uint8_t* output);
extern void carquet_neon_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                          double* values);
extern void carquet_neon_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_neon_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern int64_t carquet_neon_find_run_length_i32(const int32_t* values, int64_t count);
extern void carquet_neon_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset);
extern size_t carquet_neon_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit);
extern int64_t carquet_neon_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level);
extern void carquet_neon_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                            int16_t max_def_level, uint8_t* null_bitmap);
extern void carquet_neon_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
extern void carquet_neon_minmax_i32(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
extern void carquet_neon_minmax_i64(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
extern void carquet_neon_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_neon_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern void carquet_neon_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                          int32_t* min_value, int32_t* max_value);
extern void carquet_neon_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                          int64_t* min_value, int64_t* max_value);
extern void carquet_neon_copy_minmax_float(const float* values, int64_t count, float* output,
                                            float* min_value, float* max_value);
extern void carquet_neon_copy_minmax_double(const double* values, int64_t count, double* output,
                                             double* min_value, double* max_value);
extern void carquet_neon_bitunpack8_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_2bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_3bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_5bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_6bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_7bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack8_16bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack32_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack32_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack16_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_neon_bitunpack16_16bit(const uint8_t* input, uint32_t* values);
#endif

#if defined(CARQUET_ENABLE_SVE) && defined(__ARM_FEATURE_SVE)
extern void carquet_sve_gather_i32(const int32_t* dict, const uint32_t* indices,
                                    int64_t count, int32_t* output);
extern void carquet_sve_gather_i64(const int64_t* dict, const uint32_t* indices,
                                    int64_t count, int64_t* output);
extern void carquet_sve_gather_float(const float* dict, const uint32_t* indices,
                                      int64_t count, float* output);
extern void carquet_sve_gather_double(const double* dict, const uint32_t* indices,
                                       int64_t count, double* output);
extern bool carquet_sve_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                            const uint32_t* indices, int64_t count,
                                            int32_t* output);
extern bool carquet_sve_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                            const uint32_t* indices, int64_t count,
                                            int64_t* output);
extern bool carquet_sve_checked_gather_float(const float* dict, int32_t dict_count,
                                              const uint32_t* indices, int64_t count,
                                              float* output);
extern bool carquet_sve_checked_gather_double(const double* dict, int32_t dict_count,
                                               const uint32_t* indices, int64_t count,
                                               double* output);
extern void carquet_sve_byte_stream_split_encode_float(const float* values, int64_t count,
                                                        uint8_t* output);
extern void carquet_sve_byte_stream_split_decode_float(const uint8_t* data, int64_t count,
                                                        float* values);
extern void carquet_sve_byte_stream_split_encode_double(const double* values, int64_t count,
                                                         uint8_t* output);
extern void carquet_sve_byte_stream_split_decode_double(const uint8_t* data, int64_t count,
                                                         double* values);
extern void carquet_sve_bitunpack8_1bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_2bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_3bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_5bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_6bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_7bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_8bit(const uint8_t* input, uint32_t* values);
extern void carquet_sve_bitunpack8_16bit(const uint8_t* input, uint32_t* values);
extern int64_t carquet_sve_find_run_length_i32(const int32_t* values, int64_t count);
extern int64_t carquet_sve_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level);
extern void carquet_sve_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
extern void carquet_sve_minmax_i32(const int32_t* values, int64_t count, int32_t* min_value, int32_t* max_value);
extern void carquet_sve_minmax_i64(const int64_t* values, int64_t count, int64_t* min_value, int64_t* max_value);
extern void carquet_sve_minmax_float(const float* values, int64_t count, float* min_value, float* max_value);
extern void carquet_sve_minmax_double(const double* values, int64_t count, double* min_value, double* max_value);
extern void carquet_sve_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                         int32_t* min_value, int32_t* max_value);
extern void carquet_sve_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                         int64_t* min_value, int64_t* max_value);
extern void carquet_sve_copy_minmax_float(const float* values, int64_t count, float* output,
                                           float* min_value, float* max_value);
extern void carquet_sve_copy_minmax_double(const double* values, int64_t count, double* output,
                                            double* min_value, double* max_value);
#endif

#endif /* AArch64 */

/* ============================================================================
 * Dispatch Table
 * ============================================================================
 */

typedef struct {
    prefix_sum_i32_fn prefix_sum_i32;
    prefix_sum_i64_fn prefix_sum_i64;
    gather_i32_fn gather_i32;
    gather_i64_fn gather_i64;
    gather_float_fn gather_float;
    gather_double_fn gather_double;
    checked_gather_i32_fn checked_gather_i32;
    checked_gather_i64_fn checked_gather_i64;
    checked_gather_float_fn checked_gather_float;
    checked_gather_double_fn checked_gather_double;
    byte_split_encode_float_fn byte_split_encode_float;
    byte_split_decode_float_fn byte_split_decode_float;
    byte_split_encode_double_fn byte_split_encode_double;
    byte_split_decode_double_fn byte_split_decode_double;
    unpack_bools_fn unpack_bools;
    pack_bools_fn pack_bools;
    bitunpack8_u32_fn bitunpack8_u32[33];
    /* Optional wider kernels: produce bitunpack_wide_vals[bw] values (a
     * multiple of 8, identical to that many / 8 calls of bitunpack8_u32[bw])
     * per call. fn == NULL / vals == 0 means "no wide kernel for this width".
     * Only verified-correct, genuinely-SIMD kernels are installed here. */
    bitunpack8_u32_fn bitunpack_wide_fn[33];
    uint16_t bitunpack_wide_vals[33];
    find_run_length_i32_fn find_run_length_i32;
    match_copy_fn match_copy;
    match_length_fn match_length;
    count_non_nulls_fn count_non_nulls;
    build_null_bitmap_fn build_null_bitmap;
    fill_def_levels_fn fill_def_levels;
    minmax_i32_fn minmax_i32;
    minmax_i64_fn minmax_i64;
    minmax_float_fn minmax_float;
    minmax_double_fn minmax_double;
    copy_minmax_i32_fn copy_minmax_i32;
    copy_minmax_i64_fn copy_minmax_i64;
    copy_minmax_float_fn copy_minmax_float;
    copy_minmax_double_fn copy_minmax_double;
} carquet_simd_dispatch_t;

static carquet_simd_dispatch_t g_dispatch = {0};
static int g_dispatch_initialized = 0;

/* Acquire/release accessors for the init flag. The release store in
 * carquet_simd_dispatch_init() publishes all g_dispatch writes; the acquire
 * load in DISPATCH_ENSURE_INIT() ensures a thread that observes the flag set
 * also observes the fully-populated dispatch table. Concurrent first-use is
 * safe: init is idempotent (it always writes the same function pointers). */
static int dispatch_is_initialized(void) {
#if defined(__GNUC__) || defined(__clang__)
    return __atomic_load_n(&g_dispatch_initialized, __ATOMIC_ACQUIRE);
#elif defined(_MSC_VER)
    return _InterlockedCompareExchange((volatile long*)&g_dispatch_initialized, 1, 1);
#else
    return g_dispatch_initialized;
#endif
}

static void dispatch_set_initialized(void) {
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&g_dispatch_initialized, 1, __ATOMIC_RELEASE);
#elif defined(_MSC_VER)
    _InterlockedExchange((volatile long*)&g_dispatch_initialized, 1);
#else
    g_dispatch_initialized = 1;
#endif
}

/* Serialize the table-population critical section, mirroring carquet_init()'s
 * lock in detect.c. The acquire/release flag alone only gates the fast path;
 * without this lock two threads that both observe the flag unset (e.g. on first
 * SIMD use inside an OpenMP parallel column loop) would race writing the ~90
 * g_dispatch function pointers. */
static volatile int g_dispatch_lock = 0;

static void dispatch_lock_acquire(void) {
#if defined(__GNUC__) || defined(__clang__)
    while (__atomic_exchange_n(&g_dispatch_lock, 1, __ATOMIC_ACQUIRE)) { /* spin */ }
#elif defined(_MSC_VER)
    while (_InterlockedExchange((volatile long*)&g_dispatch_lock, 1)) { /* spin */ }
#endif
}

static void dispatch_lock_release(void) {
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&g_dispatch_lock, 0, __ATOMIC_RELEASE);
#elif defined(_MSC_VER)
    _InterlockedExchange((volatile long*)&g_dispatch_lock, 0);
#endif
}

/* ============================================================================
 * Dispatch Initialization
 * ============================================================================
 */

void carquet_simd_dispatch_init(void) {
    /* Fast path: already initialized */
    if (dispatch_is_initialized()) {
        return;
    }

    dispatch_lock_acquire();
    /* Re-check under the lock: another thread may have finished while we spun. */
    if (dispatch_is_initialized()) {
        dispatch_lock_release();
        return;
    }

    const carquet_cpu_info_t* cpu = carquet_get_cpu_info();
    (void)cpu;  /* May be unused on some platforms */

    /* Start with scalar fallbacks */
    g_dispatch.prefix_sum_i32 = scalar_prefix_sum_i32;
    g_dispatch.prefix_sum_i64 = scalar_prefix_sum_i64;
    g_dispatch.gather_i32 = scalar_gather_i32;
    g_dispatch.gather_i64 = scalar_gather_i64;
    g_dispatch.gather_float = scalar_gather_float;
    g_dispatch.gather_double = scalar_gather_double;
    g_dispatch.checked_gather_i32 = scalar_checked_gather_i32;
    g_dispatch.checked_gather_i64 = scalar_checked_gather_i64;
    g_dispatch.checked_gather_float = scalar_checked_gather_float;
    g_dispatch.checked_gather_double = scalar_checked_gather_double;
    g_dispatch.byte_split_encode_float = scalar_byte_split_encode_float;
    g_dispatch.byte_split_decode_float = scalar_byte_split_decode_float;
    g_dispatch.byte_split_encode_double = scalar_byte_split_encode_double;
    g_dispatch.byte_split_decode_double = scalar_byte_split_decode_double;
    g_dispatch.unpack_bools = scalar_unpack_bools;
    g_dispatch.pack_bools = scalar_pack_bools;
    g_dispatch.find_run_length_i32 = scalar_find_run_length_i32;
    g_dispatch.match_copy = scalar_match_copy;
    g_dispatch.match_length = scalar_match_length;
    g_dispatch.count_non_nulls = scalar_count_non_nulls;
    g_dispatch.build_null_bitmap = scalar_build_null_bitmap;
    g_dispatch.fill_def_levels = scalar_fill_def_levels;
    g_dispatch.minmax_i32 = scalar_minmax_i32;
    g_dispatch.minmax_i64 = scalar_minmax_i64;
    g_dispatch.minmax_float = scalar_minmax_float;
    g_dispatch.minmax_double = scalar_minmax_double;
    g_dispatch.copy_minmax_i32 = scalar_copy_minmax_i32;
    g_dispatch.copy_minmax_i64 = scalar_copy_minmax_i64;
    g_dispatch.copy_minmax_float = scalar_copy_minmax_float;
    g_dispatch.copy_minmax_double = scalar_copy_minmax_double;

#if defined(CARQUET_ARCH_X86)

#ifdef CARQUET_ENABLE_SSE
    if (cpu->has_sse42) {
        g_dispatch.prefix_sum_i32 = carquet_sse_prefix_sum_i32;
        g_dispatch.prefix_sum_i64 = carquet_sse_prefix_sum_i64;
        g_dispatch.gather_i32 = carquet_sse_gather_i32;
        g_dispatch.gather_i64 = carquet_sse_gather_i64;
        g_dispatch.gather_float = carquet_sse_gather_float;
        g_dispatch.gather_double = carquet_sse_gather_double;
        g_dispatch.checked_gather_i32 = carquet_sse_checked_gather_i32;
        g_dispatch.checked_gather_i64 = carquet_sse_checked_gather_i64;
        g_dispatch.checked_gather_float = carquet_sse_checked_gather_float;
        g_dispatch.checked_gather_double = carquet_sse_checked_gather_double;
        g_dispatch.byte_split_encode_float = carquet_sse_byte_stream_split_encode_float;
        g_dispatch.byte_split_decode_float = carquet_sse_byte_stream_split_decode_float;
        g_dispatch.byte_split_encode_double = carquet_sse_byte_stream_split_encode_double;
        g_dispatch.byte_split_decode_double = carquet_sse_byte_stream_split_decode_double;
        g_dispatch.unpack_bools = carquet_sse_unpack_bools;
        g_dispatch.pack_bools = carquet_sse_pack_bools;
        g_dispatch.bitunpack8_u32[1] = carquet_sse_bitunpack8_1bit;
        g_dispatch.bitunpack8_u32[2] = carquet_sse_bitunpack8_2bit;
        g_dispatch.bitunpack8_u32[3] = carquet_sse_bitunpack8_3bit;
        g_dispatch.bitunpack8_u32[4] = carquet_sse_bitunpack8_4bit;
        g_dispatch.bitunpack8_u32[5] = carquet_sse_bitunpack8_5bit;
        g_dispatch.bitunpack8_u32[6] = carquet_sse_bitunpack8_6bit;
        g_dispatch.bitunpack8_u32[7] = carquet_sse_bitunpack8_7bit;
        g_dispatch.bitunpack8_u32[8] = carquet_sse_bitunpack8_8bit;
        g_dispatch.bitunpack8_u32[16] = carquet_sse_bitunpack8_16bit;
        /* Wide: 32 x 1-bit per call (4 input bytes). Verified == 4 x the
         * scalar 1-bit unpacker by test_bitunpack_wide. */
        g_dispatch.bitunpack_wide_fn[1] = carquet_sse_bitunpack32_1bit;
        g_dispatch.bitunpack_wide_vals[1] = 32;
        g_dispatch.match_copy = carquet_sse_match_copy;
        g_dispatch.match_length = carquet_sse_match_length;
        g_dispatch.count_non_nulls = carquet_sse_count_non_nulls;
        g_dispatch.build_null_bitmap = carquet_sse_build_null_bitmap;
        g_dispatch.fill_def_levels = carquet_sse_fill_def_levels;
        g_dispatch.minmax_i32 = carquet_sse_minmax_i32;
        g_dispatch.minmax_i64 = carquet_sse_minmax_i64;
        g_dispatch.minmax_float = carquet_sse_minmax_float;
        g_dispatch.minmax_double = carquet_sse_minmax_double;
        g_dispatch.copy_minmax_i32 = carquet_sse_copy_minmax_i32;
        g_dispatch.copy_minmax_i64 = carquet_sse_copy_minmax_i64;
        g_dispatch.copy_minmax_float = carquet_sse_copy_minmax_float;
        g_dispatch.copy_minmax_double = carquet_sse_copy_minmax_double;
        g_dispatch.find_run_length_i32 = carquet_sse_find_run_length_i32;
    }
#endif

#ifdef CARQUET_ENABLE_AVX
    if (cpu->has_avx) {
        g_dispatch.byte_split_encode_float = carquet_avx_byte_stream_split_encode_float;
        g_dispatch.byte_split_decode_float = carquet_avx_byte_stream_split_decode_float;
        g_dispatch.byte_split_encode_double = carquet_avx_byte_stream_split_encode_double;
        g_dispatch.byte_split_decode_double = carquet_avx_byte_stream_split_decode_double;
        g_dispatch.minmax_float = carquet_avx_minmax_float;
        g_dispatch.minmax_double = carquet_avx_minmax_double;
        g_dispatch.copy_minmax_float = carquet_avx_copy_minmax_float;
        g_dispatch.copy_minmax_double = carquet_avx_copy_minmax_double;
    }
#endif

#ifdef CARQUET_ENABLE_AVX2
    if (cpu->has_avx2) {
        g_dispatch.prefix_sum_i32 = carquet_avx2_prefix_sum_i32;
        g_dispatch.prefix_sum_i64 = carquet_avx2_prefix_sum_i64;
        g_dispatch.gather_i32 = carquet_avx2_gather_i32;
        g_dispatch.gather_i64 = carquet_avx2_gather_i64;
        g_dispatch.gather_float = carquet_avx2_gather_float;
        g_dispatch.gather_double = carquet_avx2_gather_double;
        g_dispatch.checked_gather_i32 = carquet_avx2_checked_gather_i32;
        g_dispatch.checked_gather_i64 = carquet_avx2_checked_gather_i64;
        g_dispatch.checked_gather_float = carquet_avx2_checked_gather_float;
        g_dispatch.checked_gather_double = carquet_avx2_checked_gather_double;
        g_dispatch.byte_split_encode_float = carquet_avx2_byte_stream_split_encode_float;
        g_dispatch.byte_split_decode_float = carquet_avx2_byte_stream_split_decode_float;
        g_dispatch.byte_split_encode_double = carquet_avx2_byte_stream_split_encode_double;
        g_dispatch.byte_split_decode_double = carquet_avx2_byte_stream_split_decode_double;
        g_dispatch.unpack_bools = carquet_avx2_unpack_bools;
        g_dispatch.pack_bools = carquet_avx2_pack_bools;
        g_dispatch.match_copy = carquet_avx2_match_copy;
        g_dispatch.match_length = carquet_avx2_match_length;
        g_dispatch.count_non_nulls = carquet_avx2_count_non_nulls;
        g_dispatch.build_null_bitmap = carquet_avx2_build_null_bitmap;
        g_dispatch.fill_def_levels = carquet_avx2_fill_def_levels;
        g_dispatch.minmax_i32 = carquet_avx2_minmax_i32;
        g_dispatch.minmax_i64 = carquet_avx2_minmax_i64;
        g_dispatch.minmax_float = carquet_avx2_minmax_float;
        g_dispatch.minmax_double = carquet_avx2_minmax_double;
        g_dispatch.copy_minmax_i32 = carquet_avx2_copy_minmax_i32;
        g_dispatch.copy_minmax_i64 = carquet_avx2_copy_minmax_i64;
        g_dispatch.copy_minmax_float = carquet_avx2_copy_minmax_float;
        g_dispatch.copy_minmax_double = carquet_avx2_copy_minmax_double;
        g_dispatch.bitunpack8_u32[1] = carquet_avx2_bitunpack8_1bit;
        g_dispatch.bitunpack8_u32[2] = carquet_avx2_bitunpack8_2bit;
        g_dispatch.bitunpack8_u32[3] = carquet_avx2_bitunpack8_3bit;
        g_dispatch.bitunpack8_u32[4] = carquet_avx2_bitunpack8_4bit;
        g_dispatch.bitunpack8_u32[5] = carquet_avx2_bitunpack8_5bit;
        g_dispatch.bitunpack8_u32[6] = carquet_avx2_bitunpack8_6bit;
        g_dispatch.bitunpack8_u32[7] = carquet_avx2_bitunpack8_7bit;
        g_dispatch.bitunpack8_u32[8] = carquet_avx2_bitunpack8_8bit;
        g_dispatch.bitunpack8_u32[16] = carquet_avx2_bitunpack8_16bit;
        /* Wide: 16 values per call. Verified == 2 x the scalar unpacker
         * for these widths by test_bitunpack_wide. (1-bit stays SSE-32.) */
        g_dispatch.bitunpack_wide_fn[4] = carquet_avx2_bitunpack16_4bit;
        g_dispatch.bitunpack_wide_vals[4] = 16;
        g_dispatch.bitunpack_wide_fn[8] = carquet_avx2_bitunpack16_8bit;
        g_dispatch.bitunpack_wide_vals[8] = 16;
        g_dispatch.find_run_length_i32 = carquet_avx2_find_run_length_i32;
    }
#endif

#ifdef CARQUET_ENABLE_AVX512
    /* The AVX-512 objects are compiled with -mavx512bw/-mavx512vl, so all
     * three feature bits must be present (and OS-enabled, see detect.c). */
    if (cpu->has_avx512f && cpu->has_avx512bw && cpu->has_avx512vl) {
        g_dispatch.prefix_sum_i32 = carquet_avx512_prefix_sum_i32;
        g_dispatch.prefix_sum_i64 = carquet_avx512_prefix_sum_i64;
        g_dispatch.gather_i32 = carquet_avx512_gather_i32;
        g_dispatch.gather_i64 = carquet_avx512_gather_i64;
        g_dispatch.gather_float = carquet_avx512_gather_float;
        g_dispatch.gather_double = carquet_avx512_gather_double;
        g_dispatch.checked_gather_i32 = carquet_avx512_checked_gather_i32;
        g_dispatch.checked_gather_i64 = carquet_avx512_checked_gather_i64;
        g_dispatch.checked_gather_float = carquet_avx512_checked_gather_float;
        g_dispatch.checked_gather_double = carquet_avx512_checked_gather_double;
        g_dispatch.byte_split_encode_float = carquet_avx512_byte_stream_split_encode_float;
        g_dispatch.byte_split_decode_float = carquet_avx512_byte_stream_split_decode_float;
        g_dispatch.byte_split_encode_double = carquet_avx512_byte_stream_split_encode_double;
        g_dispatch.byte_split_decode_double = carquet_avx512_byte_stream_split_decode_double;
        g_dispatch.bitunpack8_u32[4] = carquet_avx512_bitunpack8_4bit;
        g_dispatch.bitunpack8_u32[8] = carquet_avx512_bitunpack8_8bit;
        g_dispatch.bitunpack8_u32[16] = carquet_avx512_bitunpack8_16bit;
        /* Wide: 32 (4/8-bit) or 16 (16-bit) values per call. Verified
         * against the scalar unpacker by test_bitunpack_wide. */
        g_dispatch.bitunpack_wide_fn[4] = carquet_avx512_bitunpack32_4bit;
        g_dispatch.bitunpack_wide_vals[4] = 32;
        g_dispatch.bitunpack_wide_fn[8] = carquet_avx512_bitunpack32_8bit;
        g_dispatch.bitunpack_wide_vals[8] = 32;
        g_dispatch.bitunpack_wide_fn[16] = carquet_avx512_bitunpack16_16bit;
        g_dispatch.bitunpack_wide_vals[16] = 16;
        g_dispatch.unpack_bools = carquet_avx512_unpack_bools;
        g_dispatch.pack_bools = carquet_avx512_pack_bools;
        g_dispatch.match_copy = carquet_avx512_match_copy;
        g_dispatch.match_length = carquet_avx512_match_length;
        g_dispatch.count_non_nulls = carquet_avx512_count_non_nulls;
        g_dispatch.build_null_bitmap = carquet_avx512_build_null_bitmap;
        g_dispatch.fill_def_levels = carquet_avx512_fill_def_levels;
        g_dispatch.minmax_i32 = carquet_avx512_minmax_i32;
        g_dispatch.minmax_i64 = carquet_avx512_minmax_i64;
        g_dispatch.minmax_float = carquet_avx512_minmax_float;
        g_dispatch.minmax_double = carquet_avx512_minmax_double;
        g_dispatch.find_run_length_i32 = carquet_avx512_find_run_length_i32;
    }
#endif

#endif /* CARQUET_ARCH_X86 */

#if defined(CARQUET_ARCH_ARM)

    /* Register NEON functions when the compiler can emit them and the CPU has NEON. */
#if defined(CARQUET_ENABLE_NEON) && (defined(__ARM_NEON) || defined(__ARM_NEON__))
    if (cpu->has_neon) {
    /* Prefix sums are loop-carried dependency chains. On Apple Silicon the
     * scalar compiler-generated loop is faster than the NEON shuffle-based
     * version, so keep the scalar fallback here while installing NEON where
     * it provides real throughput wins. */
    g_dispatch.gather_i32 = carquet_neon_gather_i32;
    g_dispatch.gather_i64 = carquet_neon_gather_i64;
    g_dispatch.gather_float = carquet_neon_gather_float;
    g_dispatch.gather_double = carquet_neon_gather_double;
    g_dispatch.checked_gather_i32 = carquet_neon_checked_gather_i32;
    g_dispatch.checked_gather_i64 = carquet_neon_checked_gather_i64;
    g_dispatch.checked_gather_float = carquet_neon_checked_gather_float;
    g_dispatch.checked_gather_double = carquet_neon_checked_gather_double;
    /* Byte-stream-split: keep the scalar/compiler path for float (the
     * auto-vectorized 4-byte transpose matches or beats hand-written NEON on
     * Apple Silicon, and vld4/vqtbl both regress small in-cache float decode),
     * but use NEON for double via vld4q_u16/vst4q_u16 structure load-stores.
     * Measured on M3 vs the prior vqtbl path: double encode +60-100%, decode
     * +47-69%; vs scalar, decode is +47-86% across cache regimes. Byte-exact. */
    g_dispatch.byte_split_encode_double = carquet_neon_byte_stream_split_encode_double;
    g_dispatch.byte_split_decode_double = carquet_neon_byte_stream_split_decode_double;
    g_dispatch.unpack_bools = carquet_neon_unpack_bools;
    g_dispatch.pack_bools = carquet_neon_pack_bools;
    g_dispatch.find_run_length_i32 = carquet_neon_find_run_length_i32;
    g_dispatch.match_copy = carquet_neon_match_copy;
    g_dispatch.match_length = carquet_neon_match_length;
    g_dispatch.count_non_nulls = carquet_neon_count_non_nulls;
    g_dispatch.build_null_bitmap = carquet_neon_build_null_bitmap;
    g_dispatch.fill_def_levels = carquet_neon_fill_def_levels;
    g_dispatch.minmax_i32 = carquet_neon_minmax_i32;
    /* i64 min/max is a short loop with scalar compares on NEON, and measured
     * slower than the compiler-generated scalar path on Apple Silicon. */
    g_dispatch.minmax_float = carquet_neon_minmax_float;
    g_dispatch.minmax_double = carquet_neon_minmax_double;
    g_dispatch.copy_minmax_i32 = carquet_neon_copy_minmax_i32;
    /* Same for i64 copy+minmax: keep scalar copy/min/max, which measured
     * substantially faster than the NEON implementation. */
    g_dispatch.copy_minmax_float = carquet_neon_copy_minmax_float;
    g_dispatch.copy_minmax_double = carquet_neon_copy_minmax_double;
    g_dispatch.bitunpack8_u32[1] = carquet_neon_bitunpack8_1bit;
    g_dispatch.bitunpack8_u32[2] = carquet_neon_bitunpack8_2bit;
    g_dispatch.bitunpack8_u32[3] = carquet_neon_bitunpack8_3bit;
    g_dispatch.bitunpack8_u32[4] = carquet_neon_bitunpack8_4bit;
    g_dispatch.bitunpack8_u32[5] = carquet_neon_bitunpack8_5bit;
    g_dispatch.bitunpack8_u32[6] = carquet_neon_bitunpack8_6bit;
    g_dispatch.bitunpack8_u32[7] = carquet_neon_bitunpack8_7bit;
    g_dispatch.bitunpack8_u32[8] = carquet_neon_bitunpack8_8bit;
    g_dispatch.bitunpack8_u32[16] = carquet_neon_bitunpack8_16bit;
    /* Wide: 32 x 1-bit per call (4 input bytes), == 4 calls of the
     * scalar 1-bit unpacker; verified by test_bitunpack_wide. */
    g_dispatch.bitunpack_wide_fn[1] = carquet_neon_bitunpack32_1bit;
    g_dispatch.bitunpack_wide_vals[1] = 32;
    g_dispatch.bitunpack_wide_fn[4] = carquet_neon_bitunpack32_4bit;
    g_dispatch.bitunpack_wide_vals[4] = 32;
    g_dispatch.bitunpack_wide_fn[8] = carquet_neon_bitunpack16_8bit;
    g_dispatch.bitunpack_wide_vals[8] = 16;
    g_dispatch.bitunpack_wide_fn[16] = carquet_neon_bitunpack16_16bit;
    g_dispatch.bitunpack_wide_vals[16] = 16;
    }
#endif

    /* SVE overrides NEON where SVE is genuinely better.
     * prefix_sum, unpack/pack_bools, build_null_bitmap are left as NEON
     * because their SVE implementations were pure scalar (no real benefit).
     * match_copy, match_length inherit from NEON. */
#if defined(CARQUET_ENABLE_SVE) && defined(__ARM_FEATURE_SVE)
    if (cpu->has_sve) {
        /* Gather: SVE has true hardware gather instructions */
        g_dispatch.gather_i32 = carquet_sve_gather_i32;
        g_dispatch.gather_i64 = carquet_sve_gather_i64;
        g_dispatch.gather_float = carquet_sve_gather_float;
        g_dispatch.gather_double = carquet_sve_gather_double;
        g_dispatch.checked_gather_i32 = carquet_sve_checked_gather_i32;
        g_dispatch.checked_gather_i64 = carquet_sve_checked_gather_i64;
        g_dispatch.checked_gather_float = carquet_sve_checked_gather_float;
        g_dispatch.checked_gather_double = carquet_sve_checked_gather_double;

        /* Byte stream split: SVE structure load/store (svld4/svst4) */
        g_dispatch.byte_split_encode_float = carquet_sve_byte_stream_split_encode_float;
        g_dispatch.byte_split_decode_float = carquet_sve_byte_stream_split_decode_float;
        g_dispatch.byte_split_encode_double = carquet_sve_byte_stream_split_encode_double;
        g_dispatch.byte_split_decode_double = carquet_sve_byte_stream_split_decode_double;

        /* Bit unpacking: all widths */
        g_dispatch.bitunpack8_u32[1] = carquet_sve_bitunpack8_1bit;
        g_dispatch.bitunpack8_u32[2] = carquet_sve_bitunpack8_2bit;
        g_dispatch.bitunpack8_u32[3] = carquet_sve_bitunpack8_3bit;
        g_dispatch.bitunpack8_u32[4] = carquet_sve_bitunpack8_4bit;
        g_dispatch.bitunpack8_u32[5] = carquet_sve_bitunpack8_5bit;
        g_dispatch.bitunpack8_u32[6] = carquet_sve_bitunpack8_6bit;
        g_dispatch.bitunpack8_u32[7] = carquet_sve_bitunpack8_7bit;
        g_dispatch.bitunpack8_u32[8] = carquet_sve_bitunpack8_8bit;
        g_dispatch.bitunpack8_u32[16] = carquet_sve_bitunpack8_16bit;

        /* Run detection: SVE comparison + first-fault */
        g_dispatch.find_run_length_i32 = carquet_sve_find_run_length_i32;

        /* Def levels: SVE vectorized comparison and fill */
        g_dispatch.count_non_nulls = carquet_sve_count_non_nulls;
        g_dispatch.fill_def_levels = carquet_sve_fill_def_levels;

        /* Min/max: SVE horizontal reduction */
        g_dispatch.minmax_i32 = carquet_sve_minmax_i32;
        g_dispatch.minmax_i64 = carquet_sve_minmax_i64;
        g_dispatch.minmax_float = carquet_sve_minmax_float;
        g_dispatch.minmax_double = carquet_sve_minmax_double;
        g_dispatch.copy_minmax_i32 = carquet_sve_copy_minmax_i32;
        g_dispatch.copy_minmax_i64 = carquet_sve_copy_minmax_i64;
        g_dispatch.copy_minmax_float = carquet_sve_copy_minmax_float;
        g_dispatch.copy_minmax_double = carquet_sve_copy_minmax_double;
    }
#endif

#endif /* ARM */

    dispatch_set_initialized();
    dispatch_lock_release();
}

/* ============================================================================
 * Public Dispatch Functions
 * ============================================================================
 */

/* Ensure dispatch is initialized. Uses __builtin_expect to hint that the
 * fast path (already initialized) is taken >99.99% of the time, eliminating
 * branch misprediction overhead on every dispatch call. */
#if defined(__GNUC__) || defined(__clang__)
#define DISPATCH_ENSURE_INIT() \
    do { if (__builtin_expect(!dispatch_is_initialized(), 0)) carquet_simd_dispatch_init(); } while(0)
#else
#define DISPATCH_ENSURE_INIT() \
    do { if (!dispatch_is_initialized()) carquet_simd_dispatch_init(); } while(0)
#endif

void carquet_dispatch_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.prefix_sum_i32(values, count, initial);
}

void carquet_dispatch_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.prefix_sum_i64(values, count, initial);
}

void carquet_dispatch_gather_i32(const int32_t* dict, const uint32_t* indices,
                                  int64_t count, int32_t* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.gather_i32(dict, indices, count, output);
}

void carquet_dispatch_gather_i64(const int64_t* dict, const uint32_t* indices,
                                  int64_t count, int64_t* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.gather_i64(dict, indices, count, output);
}

void carquet_dispatch_gather_float(const float* dict, const uint32_t* indices,
                                    int64_t count, float* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.gather_float(dict, indices, count, output);
}

void carquet_dispatch_gather_double(const double* dict, const uint32_t* indices,
                                     int64_t count, double* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.gather_double(dict, indices, count, output);
}

bool carquet_dispatch_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                          const uint32_t* indices, int64_t count,
                                          int32_t* output) {
    DISPATCH_ENSURE_INIT();
    /* No dictionary entries ⇒ every index is out of bounds. Guard here so the
     * SIMD fast paths never run unchecked: they compute the upper bound as
     * (uint32_t)dict_count - 1, which underflows to UINT32_MAX when
     * dict_count == 0, defeating the in-bound check and gathering from
     * arbitrary offsets. */
    if (dict_count <= 0) {
        return count == 0;
    }
#if defined(CARQUET_ARCH_ARM)
    if (g_dispatch.checked_gather_i32 &&
        g_dispatch.checked_gather_i32 != scalar_checked_gather_i32) {
        return g_dispatch.checked_gather_i32(dict, dict_count, indices, count, output);
    }
#endif
    if (!validate_gather_indices(indices, count, dict_count)) {
        return false;
    }
    g_dispatch.gather_i32(dict, indices, count, output);
    return true;
}

bool carquet_dispatch_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                          const uint32_t* indices, int64_t count,
                                          int64_t* output) {
    DISPATCH_ENSURE_INIT();
    if (dict_count <= 0) {  /* see carquet_dispatch_checked_gather_i32 */
        return count == 0;
    }
#if defined(CARQUET_ARCH_ARM)
    if (g_dispatch.checked_gather_i64 &&
        g_dispatch.checked_gather_i64 != scalar_checked_gather_i64) {
        return g_dispatch.checked_gather_i64(dict, dict_count, indices, count, output);
    }
#endif
    if (!validate_gather_indices(indices, count, dict_count)) {
        return false;
    }
    g_dispatch.gather_i64(dict, indices, count, output);
    return true;
}

bool carquet_dispatch_checked_gather_float(const float* dict, int32_t dict_count,
                                            const uint32_t* indices, int64_t count,
                                            float* output) {
    DISPATCH_ENSURE_INIT();
    if (dict_count <= 0) {  /* see carquet_dispatch_checked_gather_i32 */
        return count == 0;
    }
#if defined(CARQUET_ARCH_ARM)
    if (g_dispatch.checked_gather_float &&
        g_dispatch.checked_gather_float != scalar_checked_gather_float) {
        return g_dispatch.checked_gather_float(dict, dict_count, indices, count, output);
    }
#endif
    if (!validate_gather_indices(indices, count, dict_count)) {
        return false;
    }
    g_dispatch.gather_float(dict, indices, count, output);
    return true;
}

bool carquet_dispatch_checked_gather_double(const double* dict, int32_t dict_count,
                                             const uint32_t* indices, int64_t count,
                                             double* output) {
    DISPATCH_ENSURE_INIT();
    if (dict_count <= 0) {  /* see carquet_dispatch_checked_gather_i32 */
        return count == 0;
    }
#if defined(CARQUET_ARCH_ARM)
    if (g_dispatch.checked_gather_double &&
        g_dispatch.checked_gather_double != scalar_checked_gather_double) {
        return g_dispatch.checked_gather_double(dict, dict_count, indices, count, output);
    }
#endif
    if (!validate_gather_indices(indices, count, dict_count)) {
        return false;
    }
    g_dispatch.gather_double(dict, indices, count, output);
    return true;
}

void carquet_dispatch_byte_split_encode_float(const float* values, int64_t count,
                                               uint8_t* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.byte_split_encode_float(values, count, output);
}

void carquet_dispatch_byte_split_decode_float(const uint8_t* data, int64_t count,
                                               float* values) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.byte_split_decode_float(data, count, values);
}

void carquet_dispatch_byte_split_encode_double(const double* values, int64_t count,
                                                uint8_t* output) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.byte_split_encode_double(values, count, output);
}

void carquet_dispatch_byte_split_decode_double(const uint8_t* data, int64_t count,
                                                double* values) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.byte_split_decode_double(data, count, values);
}

void carquet_dispatch_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.unpack_bools(input, output, count);
}

void carquet_dispatch_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.pack_bools(input, output, count);
}

int64_t carquet_dispatch_find_run_length_i32(const int32_t* values, int64_t count) {
    DISPATCH_ENSURE_INIT();
    return g_dispatch.find_run_length_i32(values, count);
}

carquet_bitunpack8_fn carquet_dispatch_get_bitunpack8_fn(int bit_width) {
    DISPATCH_ENSURE_INIT();
    if (bit_width < 0 || bit_width > 32) {
        return NULL;
    }
    return g_dispatch.bitunpack8_u32[bit_width];
}

/* Wide bit-unpack accessor. Returns the number of values the wide kernel
 * for @p bit_width produces per call (a multiple of 8, identical to that
 * many / 8 scalar unpacks) and stores the kernel in *fn, or returns 0 and
 * leaves *fn untouched when there is no wide kernel for this width/ISA. */
int carquet_dispatch_get_bitunpack_wide(int bit_width, carquet_bitunpack8_fn* fn) {
    DISPATCH_ENSURE_INIT();
    if (bit_width < 1 || bit_width > 32) {
        return 0;
    }
    if (g_dispatch.bitunpack_wide_fn[bit_width] == NULL) {
        return 0;
    }
    *fn = g_dispatch.bitunpack_wide_fn[bit_width];
    return (int)g_dispatch.bitunpack_wide_vals[bit_width];
}

void carquet_dispatch_match_copy(uint8_t* dst, const uint8_t* src, size_t len, size_t offset) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.match_copy(dst, src, len, offset);
}

size_t carquet_dispatch_match_length(const uint8_t* p, const uint8_t* match, const uint8_t* limit) {
    DISPATCH_ENSURE_INIT();
    return g_dispatch.match_length(p, match, limit);
}

int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    DISPATCH_ENSURE_INIT();
    return g_dispatch.count_non_nulls(def_levels, count, max_def_level);
}

void carquet_dispatch_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                         int16_t max_def_level, uint8_t* null_bitmap) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.build_null_bitmap(def_levels, count, max_def_level, null_bitmap);
}

void carquet_dispatch_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.fill_def_levels(def_levels, count, value);
}

void carquet_dispatch_minmax_i32(const int32_t* values, int64_t count,
                                  int32_t* min_value, int32_t* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.minmax_i32(values, count, min_value, max_value);
}

void carquet_dispatch_minmax_i64(const int64_t* values, int64_t count,
                                  int64_t* min_value, int64_t* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.minmax_i64(values, count, min_value, max_value);
}

void carquet_dispatch_minmax_float(const float* values, int64_t count,
                                    float* min_value, float* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.minmax_float(values, count, min_value, max_value);
}

void carquet_dispatch_minmax_double(const double* values, int64_t count,
                                     double* min_value, double* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.minmax_double(values, count, min_value, max_value);
}

void carquet_dispatch_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                       int32_t* min_value, int32_t* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.copy_minmax_i32(values, count, output, min_value, max_value);
}

void carquet_dispatch_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                       int64_t* min_value, int64_t* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.copy_minmax_i64(values, count, output, min_value, max_value);
}

void carquet_dispatch_copy_minmax_float(const float* values, int64_t count, float* output,
                                         float* min_value, float* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.copy_minmax_float(values, count, output, min_value, max_value);
}

void carquet_dispatch_copy_minmax_double(const double* values, int64_t count, double* output,
                                          double* min_value, double* max_value) {
    DISPATCH_ENSURE_INIT();
    g_dispatch.copy_minmax_double(values, count, output, min_value, max_value);
}
