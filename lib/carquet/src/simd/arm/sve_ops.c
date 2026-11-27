/**
 * @file sve_ops.c
 * @brief SVE (Scalable Vector Extension) optimized operations for AArch64
 *
 * SVE provides scalable vectors that can be 128-2048 bits. These implementations
 * are vector-length agnostic and will automatically use the full vector width
 * available on the hardware.
 *
 * Provides SIMD-accelerated implementations of:
 * - Bit unpacking for common bit widths
 * - Byte stream split/merge (for BYTE_STREAM_SPLIT encoding)
 * - Delta decoding (prefix sums)
 * - Dictionary gather operations
 * - Boolean packing/unpacking
 */

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__aarch64__)
#ifdef __ARM_FEATURE_SVE

#include <arm_sve.h>

/* ============================================================================
 * Bit Unpacking - SVE Optimized
 * ============================================================================
 */

/**
 * Unpack 8-bit values to 32-bit using SVE.
 * Processes svcntw() elements per iteration (vector length dependent).
 */
void carquet_sve_bitunpack_8to32(const uint8_t* input, uint32_t* output, int64_t count) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* Load 8-bit values */
        svuint8_t bytes = svld1_u8(svwhilelt_b8(i, count), input + i);

        /* Widen to 16-bit, then to 32-bit */
        svuint16_t words = svunpklo_u16(bytes);
        svuint32_t dwords = svunpklo_u32(words);

        /* Store 32-bit values */
        svst1_u32(pg, output + i, dwords);

        i += svcntw();
    }
}

/**
 * Unpack 16-bit values to 32-bit using SVE.
 */
void carquet_sve_bitunpack_16to32(const uint16_t* input, uint32_t* output, int64_t count) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* Load 16-bit values */
        svuint16_t words = svld1_u16(svwhilelt_b16(i, count), input + i);

        /* Widen to 32-bit */
        svuint32_t dwords = svunpklo_u32(words);

        /* Store 32-bit values */
        svst1_u32(pg, output + i, dwords);

        i += svcntw();
    }
}

/* ============================================================================
 * Byte Stream Split - SVE Optimized
 * ============================================================================
 */

/**
 * Encode floats using byte stream split with SVE.
 * Uses SVE 4-way structure load to deinterleave bytes efficiently.
 * svld4_u8 naturally separates the 4 byte lanes of each float.
 */
void carquet_sve_byte_stream_split_encode_float(
    const float* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;

    /* svcntb() floats per iteration: svld4 produces svcntb()-element vectors */
    uint64_t vl = svcntb();

    while (i + (int64_t)vl <= count) {
        svbool_t pg = svptrue_b8();
        svuint8x4_t loaded = svld4_u8(pg, src + i * 4);
        svst1_u8(pg, output + 0 * count + i, svget4_u8(loaded, 0));
        svst1_u8(pg, output + 1 * count + i, svget4_u8(loaded, 1));
        svst1_u8(pg, output + 2 * count + i, svget4_u8(loaded, 2));
        svst1_u8(pg, output + 3 * count + i, svget4_u8(loaded, 3));
        i += (int64_t)vl;
    }

    /* Predicated tail */
    if (i < count) {
        svbool_t pg = svwhilelt_b8(i, count);
        svuint8x4_t loaded = svld4_u8(pg, src + i * 4);
        svst1_u8(pg, output + 0 * count + i, svget4_u8(loaded, 0));
        svst1_u8(pg, output + 1 * count + i, svget4_u8(loaded, 1));
        svst1_u8(pg, output + 2 * count + i, svget4_u8(loaded, 2));
        svst1_u8(pg, output + 3 * count + i, svget4_u8(loaded, 3));
    }
}

/**
 * Decode byte stream split floats using SVE.
 * Uses SVE 4-way structure store to interleave bytes efficiently.
 */
void carquet_sve_byte_stream_split_decode_float(
    const uint8_t* data,
    int64_t count,
    float* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;

    uint64_t vl = svcntb();

    while (i + (int64_t)vl <= count) {
        svbool_t pg = svptrue_b8();
        svuint8_t s0 = svld1_u8(pg, data + 0 * count + i);
        svuint8_t s1 = svld1_u8(pg, data + 1 * count + i);
        svuint8_t s2 = svld1_u8(pg, data + 2 * count + i);
        svuint8_t s3 = svld1_u8(pg, data + 3 * count + i);
        svuint8x4_t tuple = svcreate4_u8(s0, s1, s2, s3);
        svst4_u8(pg, dst + i * 4, tuple);
        i += (int64_t)vl;
    }

    if (i < count) {
        svbool_t pg = svwhilelt_b8(i, count);
        svuint8_t s0 = svld1_u8(pg, data + 0 * count + i);
        svuint8_t s1 = svld1_u8(pg, data + 1 * count + i);
        svuint8_t s2 = svld1_u8(pg, data + 2 * count + i);
        svuint8_t s3 = svld1_u8(pg, data + 3 * count + i);
        svuint8x4_t tuple = svcreate4_u8(s0, s1, s2, s3);
        svst4_u8(pg, dst + i * 4, tuple);
    }
}

/**
 * Encode doubles using byte stream split with SVE.
 * Uses svld4_u16 to deinterleave at 16-bit level (4 words per double),
 * then svuzp1/svuzp2 to split each uint16 into its two byte streams.
 * Processes svcnth() doubles per iteration.
 */
void carquet_sve_byte_stream_split_encode_double(
    const double* values,
    int64_t count,
    uint8_t* output) {

    const uint8_t* src = (const uint8_t*)values;
    int64_t i = 0;
    uint64_t vl16 = svcnth(); /* doubles per iteration */

    while (i + (int64_t)vl16 <= count) {
        svbool_t pg16 = svptrue_b16();
        svuint16x4_t loaded = svld4_u16(pg16, (const uint16_t*)(src + i * 8));

        /* Each u16 vector holds one word-position from each double.
         * Reinterpret as bytes and deinterleave low/high bytes with uzp. */
        svuint8_t v0 = svreinterpret_u8_u16(svget4_u16(loaded, 0));
        svuint8_t v1 = svreinterpret_u8_u16(svget4_u16(loaded, 1));
        svuint8_t v2 = svreinterpret_u8_u16(svget4_u16(loaded, 2));
        svuint8_t v3 = svreinterpret_u8_u16(svget4_u16(loaded, 3));

        svuint8_t zeros = svdup_n_u8(0);
        svbool_t pg_half = svwhilelt_b8((int64_t)0, (int64_t)vl16);

        svst1_u8(pg_half, output + 0 * count + i, svuzp1_u8(v0, zeros));
        svst1_u8(pg_half, output + 1 * count + i, svuzp2_u8(v0, zeros));
        svst1_u8(pg_half, output + 2 * count + i, svuzp1_u8(v1, zeros));
        svst1_u8(pg_half, output + 3 * count + i, svuzp2_u8(v1, zeros));
        svst1_u8(pg_half, output + 4 * count + i, svuzp1_u8(v2, zeros));
        svst1_u8(pg_half, output + 5 * count + i, svuzp2_u8(v2, zeros));
        svst1_u8(pg_half, output + 6 * count + i, svuzp1_u8(v3, zeros));
        svst1_u8(pg_half, output + 7 * count + i, svuzp2_u8(v3, zeros));

        i += (int64_t)vl16;
    }

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            output[b * count + i] = src[i * 8 + b];
        }
    }
}

/**
 * Decode byte stream split doubles using SVE.
 * Reverses the encode: loads from 8 byte streams, zips pairs into uint16
 * vectors, then uses svst4_u16 to interleave back into doubles.
 */
void carquet_sve_byte_stream_split_decode_double(
    const uint8_t* data,
    int64_t count,
    double* values) {

    uint8_t* dst = (uint8_t*)values;
    int64_t i = 0;
    uint64_t vl16 = svcnth();

    while (i + (int64_t)vl16 <= count) {
        svbool_t pg_half = svwhilelt_b8((int64_t)0, (int64_t)vl16);

        /* Load from 8 byte streams */
        svuint8_t s0 = svld1_u8(pg_half, data + 0 * count + i);
        svuint8_t s1 = svld1_u8(pg_half, data + 1 * count + i);
        svuint8_t s2 = svld1_u8(pg_half, data + 2 * count + i);
        svuint8_t s3 = svld1_u8(pg_half, data + 3 * count + i);
        svuint8_t s4 = svld1_u8(pg_half, data + 4 * count + i);
        svuint8_t s5 = svld1_u8(pg_half, data + 5 * count + i);
        svuint8_t s6 = svld1_u8(pg_half, data + 6 * count + i);
        svuint8_t s7 = svld1_u8(pg_half, data + 7 * count + i);

        /* Zip pairs of byte streams into uint16 vectors */
        svuint16_t w0 = svreinterpret_u16_u8(svzip1_u8(s0, s1));
        svuint16_t w1 = svreinterpret_u16_u8(svzip1_u8(s2, s3));
        svuint16_t w2 = svreinterpret_u16_u8(svzip1_u8(s4, s5));
        svuint16_t w3 = svreinterpret_u16_u8(svzip1_u8(s6, s7));

        /* Interleave 4 uint16 vectors back into doubles */
        svbool_t pg16 = svptrue_b16();
        svuint16x4_t tuple = svcreate4_u16(w0, w1, w2, w3);
        svst4_u16(pg16, (uint16_t*)(dst + i * 8), tuple);

        i += (int64_t)vl16;
    }

    /* Scalar tail */
    for (; i < count; i++) {
        for (int b = 0; b < 8; b++) {
            dst[i * 8 + b] = data[b * count + i];
        }
    }
}

/* ============================================================================
 * Delta Decoding - SVE Optimized (Prefix Sum)
 * ============================================================================
 */

/**
 * Apply prefix sum (cumulative sum) to int32 array using SVE.
 */
void carquet_sve_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB.
     * Delta encoding relies on modular arithmetic. */
    uint32_t sum = (uint32_t)initial;
    int64_t i = 0;

    /* SVE prefix sum using vector-length chunks */
    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* For correctness, we need to compute element-wise prefix */
        int64_t active = svcntp_b32(pg, pg);
        for (int64_t j = 0; j < active; j++) {
            sum += (uint32_t)values[i + j];
            values[i + j] = (int32_t)sum;
        }

        i += svcntw();
    }
}

/**
 * Apply prefix sum to int64 array using SVE.
 */
void carquet_sve_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial) {
    /* Use unsigned arithmetic to avoid signed overflow UB. */
    uint64_t sum = (uint64_t)initial;
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);

        int64_t active = svcntp_b64(pg, pg);
        for (int64_t j = 0; j < active; j++) {
            sum += (uint64_t)values[i + j];
            values[i + j] = (int64_t)sum;
        }

        i += svcntd();
    }
}

/* ============================================================================
 * Dictionary Gather - SVE Optimized
 * ============================================================================
 */

/**
 * Gather int32 values from dictionary using SVE gather instructions.
 */
void carquet_sve_gather_i32(const int32_t* dict, const uint32_t* indices,
                             int64_t count, int32_t* output) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* Load indices */
        svuint32_t idx = svld1_u32(pg, indices + i);

        /* Scale indices by 4 (sizeof(int32_t)) */
        svuint32_t offsets = svlsl_n_u32_x(pg, idx, 2);

        /* Gather values */
        svint32_t result = svld1_gather_u32offset_s32(pg, dict, offsets);

        /* Store results */
        svst1_s32(pg, output + i, result);

        i += svcntw();
    }
}

/**
 * Gather int64 values from dictionary using SVE gather instructions.
 */
void carquet_sve_gather_i64(const int64_t* dict, const uint32_t* indices,
                             int64_t count, int64_t* output) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);

        /* Load indices and extend to 64-bit */
        svuint32_t idx32 = svld1_u32(svwhilelt_b32(i, count), indices + i);
        svuint64_t idx = svunpklo_u64(idx32);

        /* Scale indices by 8 (sizeof(int64_t)) */
        svuint64_t offsets = svlsl_n_u64_x(pg, idx, 3);

        /* Gather values */
        svint64_t result = svld1_gather_u64offset_s64(pg, dict, offsets);

        /* Store results */
        svst1_s64(pg, output + i, result);

        i += svcntd();
    }
}

/**
 * Gather float values from dictionary using SVE gather instructions.
 */
void carquet_sve_gather_float(const float* dict, const uint32_t* indices,
                               int64_t count, float* output) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* Load indices */
        svuint32_t idx = svld1_u32(pg, indices + i);

        /* Scale indices by 4 (sizeof(float)) */
        svuint32_t offsets = svlsl_n_u32_x(pg, idx, 2);

        /* Gather values */
        svfloat32_t result = svld1_gather_u32offset_f32(pg, dict, offsets);

        /* Store results */
        svst1_f32(pg, output + i, result);

        i += svcntw();
    }
}

/**
 * Gather double values from dictionary using SVE gather instructions.
 */
void carquet_sve_gather_double(const double* dict, const uint32_t* indices,
                                int64_t count, double* output) {
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);

        /* Load indices and extend to 64-bit */
        svuint32_t idx32 = svld1_u32(svwhilelt_b32(i, count), indices + i);
        svuint64_t idx = svunpklo_u64(idx32);

        /* Scale indices by 8 (sizeof(double)) */
        svuint64_t offsets = svlsl_n_u64_x(pg, idx, 3);

        /* Gather values */
        svfloat64_t result = svld1_gather_u64offset_f64(pg, dict, offsets);

        /* Store results */
        svst1_f64(pg, output + i, result);

        i += svcntd();
    }
}

bool carquet_sve_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                     const uint32_t* indices, int64_t count,
                                     int32_t* output) {
    svuint32_t limit = svdup_n_u32((uint32_t)dict_count);
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);
        svuint32_t idx = svld1_u32(pg, indices + i);

        /* Bounds check: any index >= dict_count? */
        svbool_t bad = svcmpge_u32(pg, idx, limit);
        if (svptest_any(pg, bad)) {
            return false;
        }

        /* Gather in the same pass */
        svuint32_t offsets = svlsl_n_u32_x(pg, idx, 2);
        svint32_t result = svld1_gather_u32offset_s32(pg, dict, offsets);
        svst1_s32(pg, output + i, result);

        i += svcntw();
    }
    return true;
}

bool carquet_sve_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                     const uint32_t* indices, int64_t count,
                                     int64_t* output) {
    svuint32_t limit = svdup_n_u32((uint32_t)dict_count);
    int64_t i = 0;

    while (i < count) {
        svbool_t pg64 = svwhilelt_b64(i, count);
        int64_t active = svcntp_b64(pg64, pg64);
        svbool_t pg32 = svwhilelt_b32((int64_t)0, active);

        svuint32_t idx32 = svld1_u32(pg32, indices + i);

        /* Bounds check only the lanes we will actually use */
        svbool_t bad = svcmpge_u32(pg32, idx32, limit);
        if (svptest_any(pg32, bad)) {
            return false;
        }

        svuint64_t idx = svunpklo_u64(idx32);
        svuint64_t offsets = svlsl_n_u64_x(pg64, idx, 3);
        svint64_t result = svld1_gather_u64offset_s64(pg64, dict, offsets);
        svst1_s64(pg64, output + i, result);

        i += svcntd();
    }
    return true;
}

bool carquet_sve_checked_gather_float(const float* dict, int32_t dict_count,
                                       const uint32_t* indices, int64_t count,
                                       float* output) {
    return carquet_sve_checked_gather_i32(
        (const int32_t*)dict, dict_count, indices, count, (int32_t*)output);
}

bool carquet_sve_checked_gather_double(const double* dict, int32_t dict_count,
                                        const uint32_t* indices, int64_t count,
                                        double* output) {
    return carquet_sve_checked_gather_i64(
        (const int64_t*)dict, dict_count, indices, count, (int64_t*)output);
}

/* ============================================================================
 * Memcpy/Memset - SVE Optimized
 * ============================================================================
 */

/**
 * Fast memset using SVE.
 */
void carquet_sve_memset(void* dest, uint8_t value, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    svuint8_t v = svdup_n_u8(value);
    uint64_t i = 0;
    uint64_t len = (uint64_t)n;

    while (i < len) {
        svbool_t pg = svwhilelt_b8(i, len);
        svst1_u8(pg, d + i, v);
        i += svcntb();
    }
}

/**
 * Fast memcpy using SVE.
 */
void carquet_sve_memcpy(void* dest, const void* src, size_t n) {
    uint8_t* d = (uint8_t*)dest;
    const uint8_t* s = (const uint8_t*)src;
    uint64_t i = 0;
    uint64_t len = (uint64_t)n;

    while (i < len) {
        svbool_t pg = svwhilelt_b8(i, len);
        svuint8_t v = svld1_u8(pg, s + i);
        svst1_u8(pg, d + i, v);
        i += svcntb();
    }
}

/* ============================================================================
 * Boolean Operations - SVE Optimized
 * ============================================================================
 */

/**
 * Unpack boolean values from packed bits to byte array using SVE.
 */
void carquet_sve_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    /* Process one byte at a time, unpack to 8 output bytes */
    while (i < count) {
        int byte_idx = (int)(i / 8);
        uint8_t packed = input[byte_idx];

        /* Unpack 8 bits */
        int64_t remaining = count - i;
        int64_t bits_to_unpack = remaining < 8 ? remaining : 8;

        for (int64_t j = 0; j < bits_to_unpack; j++) {
            output[i + j] = (packed >> j) & 1;
        }

        i += 8;
    }
}

/**
 * Pack boolean values from byte array to packed bits using SVE.
 */
void carquet_sve_pack_bools(const uint8_t* input, uint8_t* output, int64_t count) {
    int64_t i = 0;

    while (i < count) {
        uint8_t byte = 0;
        int64_t remaining = count - i;
        int64_t bits_to_pack = remaining < 8 ? remaining : 8;

        for (int64_t j = 0; j < bits_to_pack; j++) {
            if (input[i + j]) {
                byte |= (1 << j);
            }
        }

        output[i / 8] = byte;
        i += 8;
    }
}

int64_t carquet_sve_count_non_nulls(const int16_t* def_levels, int64_t count, int16_t max_def_level) {
    int64_t non_null_count = 0;
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b16(i, count);
        svint16_t levels = svld1_s16(pg, def_levels + i);
        svbool_t matches = svcmpeq_n_s16(pg, levels, max_def_level);
        non_null_count += (int64_t)svcntp_b16(pg, matches);
        i += svcnth();
    }

    return non_null_count;
}

void carquet_sve_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                    int16_t max_def_level, uint8_t* null_bitmap) {
    int64_t i = 0;
    int64_t byte_index = 0;

    while (i < count) {
        uint8_t bits = 0;
        for (int j = 0; j < 8 && i < count; j++, i++) {
            if (def_levels[i] == max_def_level) {
                bits |= (uint8_t)(1u << j);
            }
        }
        null_bitmap[byte_index++] = bits;
    }
}

void carquet_sve_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value) {
    int64_t i = 0;
    svint16_t val = svdup_n_s16(value);

    while (i < count) {
        svbool_t pg = svwhilelt_b16(i, count);
        svst1_s16(pg, def_levels + i, val);
        i += svcnth();
    }
}

void carquet_sve_minmax_i32(const int32_t* values, int64_t count,
                             int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    int64_t i = 1;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);
        svint32_t v = svld1_s32(pg, values + i);
        int32_t chunk_min = svminv_s32(pg, v);
        int32_t chunk_max = svmaxv_s32(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntw();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_minmax_i64(const int64_t* values, int64_t count,
                             int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    int64_t i = 1;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);
        svint64_t v = svld1_s64(pg, values + i);
        int64_t chunk_min = svminv_s64(pg, v);
        int64_t chunk_max = svmaxv_s64(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntd();
    }

    *min_value = min_v;
    *max_value = max_v;
}

/* Float chunks reduce with FMINNMV/FMAXNMV, which ignore NaN lanes (NaN only
 * if every active lane is NaN, then the compare below skips it). FMINV/FMAXV
 * would return NaN for the whole chunk and the compare would drop all of its
 * values. A NaN values[0] still seeds NaN, which the writer rescans. */
void carquet_sve_minmax_float(const float* values, int64_t count,
                               float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    int64_t i = 1;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);
        svfloat32_t v = svld1_f32(pg, values + i);
        float chunk_min = svminnmv_f32(pg, v);
        float chunk_max = svmaxnmv_f32(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntw();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_minmax_double(const double* values, int64_t count,
                                double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    int64_t i = 1;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);
        svfloat64_t v = svld1_f64(pg, values + i);
        double chunk_min = svminnmv_f64(pg, v);
        double chunk_max = svmaxnmv_f64(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntd();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                  int32_t* min_value, int32_t* max_value) {
    int32_t min_v = values[0];
    int32_t max_v = values[0];
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);
        svint32_t v = svld1_s32(pg, values + i);
        svst1_s32(pg, output + i, v);
        int32_t chunk_min = svminv_s32(pg, v);
        int32_t chunk_max = svmaxv_s32(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntw();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_copy_minmax_i64(const int64_t* values, int64_t count, int64_t* output,
                                  int64_t* min_value, int64_t* max_value) {
    int64_t min_v = values[0];
    int64_t max_v = values[0];
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);
        svint64_t v = svld1_s64(pg, values + i);
        svst1_s64(pg, output + i, v);
        int64_t chunk_min = svminv_s64(pg, v);
        int64_t chunk_max = svmaxv_s64(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntd();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_copy_minmax_float(const float* values, int64_t count, float* output,
                                    float* min_value, float* max_value) {
    float min_v = values[0];
    float max_v = values[0];
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);
        svfloat32_t v = svld1_f32(pg, values + i);
        svst1_f32(pg, output + i, v);
        float chunk_min = svminnmv_f32(pg, v);
        float chunk_max = svmaxnmv_f32(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntw();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_copy_minmax_double(const double* values, int64_t count, double* output,
                                     double* min_value, double* max_value) {
    double min_v = values[0];
    double max_v = values[0];
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b64(i, count);
        svfloat64_t v = svld1_f64(pg, values + i);
        svst1_f64(pg, output + i, v);
        double chunk_min = svminnmv_f64(pg, v);
        double chunk_max = svmaxnmv_f64(pg, v);
        if (chunk_min < min_v) min_v = chunk_min;
        if (chunk_max > max_v) max_v = chunk_max;
        i += svcntd();
    }

    *min_value = min_v;
    *max_value = max_v;
}

void carquet_sve_bitunpack8_1bit(const uint8_t* input, uint32_t* values) {
    uint8_t byte_val = input[0];
    for (int i = 0; i < 8; i++) {
        values[i] = (byte_val >> i) & 1;
    }
}

void carquet_sve_bitunpack8_2bit(const uint8_t* input, uint32_t* values) {
    uint16_t v;
    memcpy(&v, input, 2);
    for (int i = 0; i < 8; i++) {
        values[i] = (v >> (i * 2)) & 0x3;
    }
}

void carquet_sve_bitunpack8_3bit(const uint8_t* input, uint32_t* values) {
    uint32_t v = 0;
    memcpy(&v, input, 3);
    for (int i = 0; i < 8; i++) {
        values[i] = (v >> (i * 3)) & 0x7;
    }
}

void carquet_sve_bitunpack8_4bit(const uint8_t* input, uint32_t* values) {
    uint32_t v = (uint32_t)input[0] | ((uint32_t)input[1] << 8) |
                 ((uint32_t)input[2] << 16) | ((uint32_t)input[3] << 24);
    for (int i = 0; i < 8; i++) {
        values[i] = (v >> (i * 4)) & 0xF;
    }
}

void carquet_sve_bitunpack8_5bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = 0;
    memcpy(&v, input, 5);
    for (int i = 0; i < 8; i++) {
        values[i] = (uint32_t)((v >> (i * 5)) & 0x1F);
    }
}

void carquet_sve_bitunpack8_6bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = 0;
    memcpy(&v, input, 6);
    for (int i = 0; i < 8; i++) {
        values[i] = (uint32_t)((v >> (i * 6)) & 0x3F);
    }
}

void carquet_sve_bitunpack8_7bit(const uint8_t* input, uint32_t* values) {
    uint64_t v = 0;
    memcpy(&v, input, 7);
    for (int i = 0; i < 8; i++) {
        values[i] = (uint32_t)((v >> (i * 7)) & 0x7F);
    }
}

void carquet_sve_bitunpack8_8bit(const uint8_t* input, uint32_t* values) {
    carquet_sve_bitunpack_8to32(input, values, 8);
}

void carquet_sve_bitunpack8_16bit(const uint8_t* input, uint32_t* values) {
    carquet_sve_bitunpack_16to32((const uint16_t*)input, values, 8);
}

/* ============================================================================
 * Run Detection - SVE Optimized
 * ============================================================================
 */

/**
 * Find the length of a run of repeated int32 values.
 */
int64_t carquet_sve_find_run_length_i32(const int32_t* values, int64_t count) {
    if (count == 0) return 0;

    int32_t first = values[0];
    svint32_t target = svdup_n_s32(first);
    int64_t i = 0;

    while (i < count) {
        svbool_t pg = svwhilelt_b32(i, count);

        /* Load values */
        svint32_t v = svld1_s32(pg, values + i);

        /* Compare with target */
        svbool_t cmp = svcmpeq_s32(pg, v, target);

        /* Check if all active elements match */
        if (!svptest_first(pg, svnot_b_z(pg, cmp))) {
            /* All match, continue */
            i += svcntw();
        } else {
            /* Found mismatch, find exact position */
            for (int64_t j = i; j < count && j < i + (int64_t)svcntw(); j++) {
                if (values[j] != first) {
                    return j;
                }
            }
            break;
        }
    }

    return count;
}

/* ============================================================================
 * Vector-Length Query
 * ============================================================================
 */

/**
 * Get the SVE vector length in bytes.
 */
size_t carquet_sve_get_vector_length_bytes(void) {
    return svcntb();
}

/**
 * Get the SVE vector length in 32-bit elements.
 */
size_t carquet_sve_get_vector_length_32(void) {
    return svcntw();
}

/**
 * Get the SVE vector length in 64-bit elements.
 */
size_t carquet_sve_get_vector_length_64(void) {
    return svcntd();
}

#endif /* __ARM_FEATURE_SVE */
#endif /* AArch64 */
