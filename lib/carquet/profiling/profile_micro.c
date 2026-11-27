/**
 * @file profile_micro.c
 * @brief Micro-benchmarks for isolated component profiling
 *
 * Isolates specific bottlenecks for detailed perf/sample analysis:
 * - RLE encoding & decoding (levels)
 * - Dictionary gather operations
 * - Null bitmap construction
 * - Compression/decompression (LZ4, Snappy, ZSTD, GZIP)
 * - SIMD dispatch overhead
 * - MinMax computation (statistics)
 * - Byte-stream-split encode/decode
 * - Prefix sum (delta decoding)
 * - Bloom filter insert/check
 * - CRC32 computation
 * - Encoding benchmarks (plain, delta, dictionary, BSS)
 *
 * Usage:
 *   ./profile_micro --component rle --iterations 1000000
 *   ./profile_micro --component all
 *   sample ./profile_micro --component encoding   # macOS
 *   perf record -g ./profile_micro --component all # Linux
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

/* Include internal headers for direct component access */
#include "carquet/carquet.h"
#include "encoding/rle.h"
#include "encoding/plain.h"
#include "core/buffer.h"

/* Forward declarations for internal SIMD dispatch */
extern void carquet_dispatch_gather_i32(const int32_t* dict, const uint32_t* indices,
                                         int64_t count, int32_t* output);
extern void carquet_dispatch_gather_i64(const int64_t* dict, const uint32_t* indices,
                                         int64_t count, int64_t* output);
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                 int16_t max_def);
extern void carquet_dispatch_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                                int16_t max_def, uint8_t* bitmap);

/* MinMax dispatch */
extern void carquet_dispatch_minmax_i32(const int32_t* values, int64_t count,
                                         int32_t* min_value, int32_t* max_value);
extern void carquet_dispatch_minmax_i64(const int64_t* values, int64_t count,
                                         int64_t* min_value, int64_t* max_value);
extern void carquet_dispatch_minmax_float(const float* values, int64_t count,
                                           float* min_value, float* max_value);
extern void carquet_dispatch_minmax_double(const double* values, int64_t count,
                                            double* min_value, double* max_value);

/* Copy+MinMax dispatch */
extern void carquet_dispatch_copy_minmax_i32(const int32_t* values, int64_t count,
                                              int32_t* output, int32_t* min_value, int32_t* max_value);
extern void carquet_dispatch_copy_minmax_i64(const int64_t* values, int64_t count,
                                              int64_t* output, int64_t* min_value, int64_t* max_value);
extern void carquet_dispatch_copy_minmax_float(const float* values, int64_t count,
                                                float* output, float* min_value, float* max_value);
extern void carquet_dispatch_copy_minmax_double(const double* values, int64_t count,
                                                 double* output, double* min_value, double* max_value);

/* Byte-stream-split dispatch */
extern void carquet_dispatch_byte_split_encode_float(const float* values, int64_t count, uint8_t* output);
extern void carquet_dispatch_byte_split_decode_float(const uint8_t* data, int64_t count, float* values);
extern void carquet_dispatch_byte_split_encode_double(const double* values, int64_t count, uint8_t* output);
extern void carquet_dispatch_byte_split_decode_double(const uint8_t* data, int64_t count, double* values);

/* Prefix sum dispatch */
extern void carquet_dispatch_prefix_sum_i32(int32_t* values, int64_t count, int32_t initial);
extern void carquet_dispatch_prefix_sum_i64(int64_t* values, int64_t count, int64_t initial);

/* LZ4 internal functions */
extern carquet_status_t carquet_lz4_compress(const uint8_t* src, size_t src_len,
                                              uint8_t* dst, size_t dst_capacity,
                                              size_t* dst_size);
extern carquet_status_t carquet_lz4_decompress(const uint8_t* src, size_t src_len,
                                                uint8_t* dst, size_t dst_capacity,
                                                size_t* dst_size);

/* Snappy internal functions */
extern size_t carquet_snappy_compress_bound(size_t source_length);
extern carquet_status_t carquet_snappy_compress(const uint8_t* input, size_t input_length,
                                                 uint8_t* output, size_t output_capacity,
                                                 size_t* output_length);
extern carquet_status_t carquet_snappy_decompress(const uint8_t* input, size_t input_length,
                                                    uint8_t* output, size_t output_capacity,
                                                    size_t* output_length);

/* ZSTD internal functions */
extern int carquet_zstd_compress(const uint8_t* src, size_t src_size,
                                  uint8_t* dst, size_t dst_capacity,
                                  size_t* dst_size, int level);
extern int carquet_zstd_decompress(const uint8_t* src, size_t src_size,
                                    uint8_t* dst, size_t dst_capacity,
                                    size_t* dst_size);
extern size_t carquet_zstd_compress_bound(size_t src_size);

/* GZIP internal functions */
extern int carquet_gzip_compress(const uint8_t* src, size_t src_size,
                                  uint8_t* dst, size_t dst_capacity,
                                  size_t* dst_size, int level);
extern int carquet_gzip_decompress(const uint8_t* src, size_t src_size,
                                    uint8_t* dst, size_t dst_capacity,
                                    size_t* dst_size);
extern size_t carquet_gzip_compress_bound(size_t src_size);

/* Delta encoding */
extern carquet_status_t carquet_delta_encode_int32(const int32_t* values, int32_t num_values,
                                                     uint8_t* data, size_t data_capacity,
                                                     size_t* bytes_written);
extern carquet_status_t carquet_delta_decode_int32(const uint8_t* data, size_t data_size,
                                                     int32_t* values, int32_t num_values,
                                                     size_t* bytes_consumed);

/* Dictionary encoding */
extern carquet_status_t carquet_dictionary_encode_int32(const int32_t* values, int64_t count,
                                                          carquet_buffer_t* dict_output,
                                                          carquet_buffer_t* indices_output);

/* Byte-stream-split encoding */
extern carquet_status_t carquet_byte_stream_split_encode_float(const float* values, int64_t count,
                                                                 uint8_t* output, size_t output_capacity,
                                                                 size_t* bytes_written);
extern carquet_status_t carquet_byte_stream_split_decode_float(const uint8_t* data, size_t data_size,
                                                                 float* values, int64_t count);

/* Bloom filter */
extern carquet_bloom_filter_t* carquet_bloom_filter_create_with_ndv(int64_t ndv, double fpp);
extern void carquet_bloom_filter_destroy(carquet_bloom_filter_t* filter);
extern void carquet_bloom_filter_insert_i32(carquet_bloom_filter_t* filter, int32_t value);
extern void carquet_bloom_filter_insert_i64(carquet_bloom_filter_t* filter, int64_t value);
extern bool carquet_bloom_filter_check_i32(const carquet_bloom_filter_t* filter, int32_t value);
extern bool carquet_bloom_filter_check_i64(const carquet_bloom_filter_t* filter, int64_t value);

/* CRC32 */
extern uint32_t carquet_crc32(const uint8_t* data, size_t length);

/* XXHash */
extern uint64_t carquet_xxhash64(const void* data, size_t length, uint64_t seed);

/* ============================================================================
 * Timing
 * ============================================================================ */

#ifdef _WIN32
static double get_time_ns(void) {
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart * 1e9 / (double)freq.QuadPart;
}
#else
static double get_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}
#endif

#define BENCH_START() double _bench_start = get_time_ns()
#define BENCH_END() (get_time_ns() - _bench_start)

/* Prevent compiler from optimizing away */
volatile int64_t g_sink;

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#define ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define NOINLINE
#define ALWAYS_INLINE inline
#endif

/* ============================================================================
 * RLE Micro-benchmark
 * ============================================================================ */

/*
 * Generates RLE-encoded data that mimics real Parquet definition levels.
 * Pattern: mostly non-null (1) with occasional null (0) values.
 */
static void generate_rle_test_data(uint8_t** out_data, size_t* out_size,
                                    int64_t num_values, double null_ratio) {
    /* Encode using the RLE encoder */
    carquet_buffer_t buf;
    carquet_buffer_init_capacity(&buf, 4096);

    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, &buf, 1); /* bit_width=1 for def levels */

    uint32_t seed = 42;
    for (int64_t i = 0; i < num_values; i++) {
        seed = seed * 1103515245 + 12345;
        double r = (double)(seed >> 16) / 32768.0;
        uint32_t val = (r >= null_ratio) ? 1 : 0;
        carquet_rle_encoder_put(&enc, val);
    }
    carquet_rle_encoder_flush(&enc);

    *out_data = malloc(buf.size);
    memcpy(*out_data, buf.data, buf.size);
    *out_size = buf.size;

    carquet_buffer_destroy(&buf);
}

NOINLINE static void bench_rle_decode_single(const uint8_t* data, size_t size,
                                              int bit_width, int64_t count,
                                              int64_t iterations) {
    int16_t* output = malloc(count * sizeof(int16_t));

    printf("  RLE single-value decode: ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_rle_decoder_t dec;
        carquet_rle_decoder_init(&dec, data, size, bit_width);

        for (int64_t i = 0; i < count; i++) {
            output[i] = (int16_t)carquet_rle_decoder_get(&dec);
        }
        g_sink = output[count / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    double values_per_sec = 1e9 / ns_per_value;

    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, values_per_sec / 1e6);

    free(output);
}

NOINLINE static void bench_rle_decode_batch(const uint8_t* data, size_t size,
                                             int bit_width, int64_t count,
                                             int64_t iterations) {
    uint32_t* output = malloc(count * sizeof(uint32_t));

    printf("  RLE batch decode:        ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_rle_decoder_t dec;
        carquet_rle_decoder_init(&dec, data, size, bit_width);

        int64_t decoded = carquet_rle_decoder_get_batch(&dec, output, count);
        g_sink = output[decoded / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    double values_per_sec = 1e9 / ns_per_value;

    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, values_per_sec / 1e6);

    free(output);
}

NOINLINE static void bench_rle_decode_levels(const uint8_t* data, size_t size,
                                              int bit_width, int64_t count,
                                              int64_t iterations) {
    int16_t* output = malloc(count * sizeof(int16_t));

    printf("  RLE decode_levels API:   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t decoded = carquet_rle_decode_levels(data, size, bit_width,
                                                     output, count);
        g_sink = output[decoded / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    double values_per_sec = 1e9 / ns_per_value;

    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, values_per_sec / 1e6);

    free(output);
}

static void run_rle_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== RLE Decoding Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    /* Generate test data with 10% nulls */
    uint8_t* rle_data;
    size_t rle_size;
    generate_rle_test_data(&rle_data, &rle_size, count, 0.10);

    printf("RLE data size: %zu bytes (%.2f bytes/value)\n\n",
           rle_size, (double)rle_size / count);

    bench_rle_decode_single(rle_data, rle_size, 1, count, iterations);
    bench_rle_decode_batch(rle_data, rle_size, 1, count, iterations);
    bench_rle_decode_levels(rle_data, rle_size, 1, count, iterations);

    free(rle_data);
}

/* ============================================================================
 * Dictionary Gather Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_gather_i32_scalar(const int32_t* dict, const uint32_t* indices,
                                              int64_t count, int32_t* output,
                                              int64_t iterations) {
    printf("  Gather i32 (scalar):     ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        for (int64_t i = 0; i < count; i++) {
            output[i] = dict[indices[i]];
        }
        g_sink = output[count / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_gather_i32_dispatch(const int32_t* dict, const uint32_t* indices,
                                                int64_t count, int32_t* output,
                                                int64_t iterations) {
    printf("  Gather i32 (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_gather_i32(dict, indices, count, output);
        g_sink = output[count / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_gather_i64_dispatch(const int64_t* dict, const uint32_t* indices,
                                                int64_t count, int64_t* output,
                                                int64_t iterations) {
    printf("  Gather i64 (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_gather_i64(dict, indices, count, output);
        g_sink = output[count / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

static void run_gather_benchmark_sized(int64_t count, int64_t iterations, int dict_size) {
    int32_t* dict_i32 = malloc((size_t)dict_size * sizeof(int32_t));
    int64_t* dict_i64 = malloc((size_t)dict_size * sizeof(int64_t));
    for (int i = 0; i < dict_size; i++) {
        dict_i32[i] = i * 7 + 13;
        dict_i64[i] = (int64_t)i * 7 + 13;
    }

    /* Generate random indices (xorshift for full 32-bit range) */
    uint32_t* indices = malloc((size_t)count * sizeof(uint32_t));
    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        indices[i] = seed % (uint32_t)dict_size;
    }

    int32_t* output_i32 = malloc((size_t)count * sizeof(int32_t));
    int64_t* output_i64 = malloc((size_t)count * sizeof(int64_t));

    double dict_mb = (double)dict_size * sizeof(int32_t) / (1024.0 * 1024.0);
    printf("\nDictionary size: %d entries (%.2f MB)\n", dict_size, dict_mb);

    bench_gather_i32_scalar(dict_i32, indices, count, output_i32, iterations);
    bench_gather_i32_dispatch(dict_i32, indices, count, output_i32, iterations);
    bench_gather_i64_dispatch(dict_i64, indices, count, output_i64, iterations);

    free(dict_i32);
    free(dict_i64);
    free(indices);
    free(output_i32);
    free(output_i64);
}

static void run_gather_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Dictionary Gather Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n", (long)count, (long)iterations);

    /* Test with different dictionary sizes to measure cache effects:
     * - 1K entries (4KB) - fits in L1 cache
     * - 100K entries (400KB) - fits in L2/L3 cache
     * - 10M entries (40MB) - exceeds L3 cache (tests memory bandwidth)
     */
    run_gather_benchmark_sized(count, iterations, 1000);      /* L1 cache */
    run_gather_benchmark_sized(count, iterations, 100000);    /* L2/L3 cache */
    run_gather_benchmark_sized(count, iterations / 10, 10000000); /* Main memory */
}

/* ============================================================================
 * Null Bitmap Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_count_nulls_scalar(const int16_t* def_levels, int64_t count,
                                               int64_t iterations) {
    printf("  Count non-nulls (scalar):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t non_null = 0;
        for (int64_t i = 0; i < count; i++) {
            if (def_levels[i] == 1) non_null++;
        }
        g_sink = non_null;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_count_nulls_dispatch(const int16_t* def_levels, int64_t count,
                                                  int64_t iterations) {
    printf("  Count non-nulls (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t non_null = carquet_dispatch_count_non_nulls(def_levels, count, 1);
        g_sink = non_null;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_build_bitmap_scalar(const int16_t* def_levels, int64_t count,
                                                uint8_t* bitmap, int64_t iterations) {
    printf("  Build bitmap (scalar):      ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t bytes = (count + 7) / 8;
        memset(bitmap, 0, bytes);
        for (int64_t i = 0; i < count; i++) {
            if (def_levels[i] < 1) {
                bitmap[i / 8] |= (1 << (i % 8));
            }
        }
        g_sink = bitmap[0];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_build_bitmap_dispatch(const int16_t* def_levels, int64_t count,
                                                   uint8_t* bitmap, int64_t iterations) {
    printf("  Build bitmap (dispatch):    ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_build_null_bitmap(def_levels, count, 1, bitmap);
        g_sink = bitmap[0];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

static void run_null_bitmap_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Null Bitmap Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    /* Generate def levels with 10% nulls */
    int16_t* def_levels = malloc(count * sizeof(int16_t));
    uint32_t seed = 42;
    int64_t null_count = 0;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        double r = (double)(seed >> 16) / 32768.0;
        def_levels[i] = (r >= 0.10) ? 1 : 0;
        if (def_levels[i] == 0) null_count++;
    }

    uint8_t* bitmap = malloc((count + 7) / 8);

    printf("Null ratio: %.1f%%\n\n", (double)null_count / count * 100);

    bench_count_nulls_scalar(def_levels, count, iterations);
    bench_count_nulls_dispatch(def_levels, count, iterations);
    bench_build_bitmap_scalar(def_levels, count, bitmap, iterations);
    bench_build_bitmap_dispatch(def_levels, count, bitmap, iterations);

    free(def_levels);
    free(bitmap);
}

/* ============================================================================
 * Compression Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_lz4_compress(const uint8_t* input, size_t input_size,
                                         uint8_t* output, size_t output_size,
                                         int64_t iterations) {
    printf("  LZ4 compress:    ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t compressed_size = 0;
        carquet_lz4_compress(input, input_size, output, output_size, &compressed_size);
        g_sink = (int64_t)compressed_size;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)input_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_lz4_decompress(const uint8_t* compressed, size_t comp_size,
                                           uint8_t* output, size_t output_size,
                                           int64_t iterations) {
    printf("  LZ4 decompress:  ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t decompressed_size = 0;
        carquet_lz4_decompress(compressed, comp_size, output, output_size, &decompressed_size);
        g_sink = (int64_t)decompressed_size;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)output_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_snappy_compress(const uint8_t* input, size_t input_size,
                                            uint8_t* output, size_t output_size,
                                            int64_t iterations) {
    printf("  Snappy compress: ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_snappy_compress(input, input_size, output, output_size, &out_len);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)input_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_snappy_decompress(const uint8_t* compressed, size_t comp_size,
                                              uint8_t* output, size_t output_size,
                                              int64_t iterations) {
    printf("  Snappy decompress: ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_snappy_decompress(compressed, comp_size, output, output_size, &out_len);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)output_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_zstd_compress(const uint8_t* input, size_t input_size,
                                          uint8_t* output, size_t output_size,
                                          int level, int64_t iterations) {
    printf("  ZSTD compress (L%d): ", level);
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_zstd_compress(input, input_size, output, output_size, &out_len, level);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)input_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_zstd_decompress(const uint8_t* compressed, size_t comp_size,
                                            uint8_t* output, size_t output_size,
                                            int64_t iterations) {
    printf("  ZSTD decompress:   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_zstd_decompress(compressed, comp_size, output, output_size, &out_len);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)output_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_gzip_compress(const uint8_t* input, size_t input_size,
                                          uint8_t* output, size_t output_size,
                                          int level, int64_t iterations) {
    printf("  GZIP compress (L%d): ", level);
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_gzip_compress(input, input_size, output, output_size, &out_len, level);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)input_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_gzip_decompress(const uint8_t* compressed, size_t comp_size,
                                            uint8_t* output, size_t output_size,
                                            int64_t iterations) {
    printf("  GZIP decompress:   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t out_len = 0;
        carquet_gzip_decompress(compressed, comp_size, output, output_size, &out_len);
        g_sink = (int64_t)out_len;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)output_size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

static void run_compression_benchmarks(size_t size, int64_t iterations) {
    printf("\n=== Compression Benchmarks ===\n");
    printf("Data size: %zu bytes, Iterations: %ld\n\n", size, (long)iterations);

    /* Generate compressible data (simulating Parquet column data) */
    uint8_t* input = malloc(size);
    uint32_t seed = 42;
    for (size_t i = 0; i < size; i++) {
        seed = seed * 1103515245 + 12345;
        /* Mix of patterns: some repetition for compressibility */
        input[i] = (uint8_t)((seed >> 16) % 64);
    }

    size_t max_compressed = size * 2;
    uint8_t* compressed = malloc(max_compressed);
    uint8_t* decompressed = malloc(size);

    /* LZ4 */
    size_t comp_size = 0;
    carquet_lz4_compress(input, size, compressed, max_compressed, &comp_size);
    printf("LZ4 ratio: %.2fx\n", (double)size / (double)comp_size);
    bench_lz4_compress(input, size, compressed, max_compressed, iterations);
    bench_lz4_decompress(compressed, comp_size, decompressed, size, iterations);

    printf("\n");

    /* Snappy */
    size_t snappy_size = 0;
    carquet_snappy_compress(input, size, compressed, max_compressed, &snappy_size);
    printf("Snappy ratio: %.2fx\n", (double)size / snappy_size);
    bench_snappy_compress(input, size, compressed, max_compressed, iterations);
    bench_snappy_decompress(compressed, snappy_size, decompressed, size, iterations);

    printf("\n");

    /* ZSTD */
    size_t zstd_size = 0;
    carquet_zstd_compress(input, size, compressed, max_compressed, &zstd_size, 1);
    printf("ZSTD L1 ratio: %.2fx\n", (double)size / zstd_size);
    bench_zstd_compress(input, size, compressed, max_compressed, 1, iterations);
    bench_zstd_decompress(compressed, zstd_size, decompressed, size, iterations);

    size_t zstd_size3 = 0;
    carquet_zstd_compress(input, size, compressed, max_compressed, &zstd_size3, 3);
    printf("  ZSTD L3 ratio: %.2fx\n", (double)size / zstd_size3);
    bench_zstd_compress(input, size, compressed, max_compressed, 3, iterations);

    printf("\n");

    /* GZIP */
    size_t gzip_size = 0;
    carquet_gzip_compress(input, size, compressed, max_compressed, &gzip_size, 6);
    printf("GZIP L6 ratio: %.2fx\n", (double)size / gzip_size);
    bench_gzip_compress(input, size, compressed, max_compressed, 6, iterations);
    bench_gzip_decompress(compressed, gzip_size, decompressed, size, iterations);

    free(input);
    free(compressed);
    free(decompressed);
}

/* ============================================================================
 * Dispatch Overhead Benchmark
 * ============================================================================ */

static void run_dispatch_overhead_benchmark(int64_t iterations) {
    printf("\n=== Dispatch Overhead Benchmark ===\n");
    printf("Iterations: %ld\n\n", (long)iterations);

    /* Minimal data to isolate function call overhead */
    int32_t dict[4] = {1, 2, 3, 4};
    uint32_t indices[4] = {0, 1, 2, 3};
    int32_t output[4];

    double elapsed_direct, elapsed_dispatch;

    printf("  Direct call (4 values): ");
    fflush(stdout);

    {
        BENCH_START();
        for (int64_t iter = 0; iter < iterations; iter++) {
            for (int i = 0; i < 4; i++) {
                output[i] = dict[indices[i]];
            }
            g_sink = output[0];
        }
        elapsed_direct = BENCH_END();
    }
    printf("%.2f ns/call\n", elapsed_direct / iterations);

    printf("  Dispatch call (4 values): ");
    fflush(stdout);

    {
        BENCH_START();
        for (int64_t iter = 0; iter < iterations; iter++) {
            carquet_dispatch_gather_i32(dict, indices, 4, output);
            g_sink = output[0];
        }
        elapsed_dispatch = BENCH_END();
    }
    printf("%.2f ns/call\n", elapsed_dispatch / iterations);

    printf("  Overhead: %.2f ns (%.1fx)\n",
           (elapsed_dispatch - elapsed_direct) / iterations,
           elapsed_dispatch / elapsed_direct);
}

/* ============================================================================
 * RLE Encoding Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_rle_encode(int64_t count, int64_t iterations) {
    printf("  RLE encode (10%% nulls):  ");
    fflush(stdout);

    /* Generate def level values */
    int16_t* levels = malloc(count * sizeof(int16_t));
    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        levels[i] = ((double)(seed >> 16) / 32768.0 >= 0.10) ? 1 : 0;
    }

    carquet_buffer_t buf;
    carquet_buffer_init_capacity(&buf, count);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        buf.size = 0;
        carquet_rle_encode_levels(levels, count, 1, &buf);
        g_sink = (int64_t)buf.size;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    carquet_buffer_destroy(&buf);
    free(levels);
}

NOINLINE static void bench_rle_encode_repeat(int64_t count, int64_t iterations) {
    printf("  RLE encode (all non-null): ");
    fflush(stdout);

    /* All 1s - exercises run-length path heavily */
    int16_t* levels = malloc(count * sizeof(int16_t));
    for (int64_t i = 0; i < count; i++) levels[i] = 1;

    carquet_buffer_t buf;
    carquet_buffer_init_capacity(&buf, count);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        buf.size = 0;
        carquet_rle_encode_levels(levels, count, 1, &buf);
        g_sink = (int64_t)buf.size;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    carquet_buffer_destroy(&buf);
    free(levels);
}

static void run_rle_encode_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== RLE Encoding Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    bench_rle_encode(count, iterations);
    bench_rle_encode_repeat(count, iterations);
}

/* ============================================================================
 * MinMax Micro-benchmark (statistics computation on write path)
 * ============================================================================ */

NOINLINE static void bench_minmax_i32_scalar(const int32_t* values, int64_t count,
                                              int64_t iterations) {
    printf("  MinMax i32 (scalar):     ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int32_t mn = values[0], mx = values[0];
        for (int64_t i = 1; i < count; i++) {
            if (values[i] < mn) mn = values[i];
            if (values[i] > mx) mx = values[i];
        }
        g_sink = mn + mx;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_minmax_i32_dispatch(const int32_t* values, int64_t count,
                                                int64_t iterations) {
    printf("  MinMax i32 (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int32_t mn, mx;
        carquet_dispatch_minmax_i32(values, count, &mn, &mx);
        g_sink = mn + mx;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_minmax_i64_dispatch(const int64_t* values, int64_t count,
                                                int64_t iterations) {
    printf("  MinMax i64 (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t mn, mx;
        carquet_dispatch_minmax_i64(values, count, &mn, &mx);
        g_sink = mn + mx;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_minmax_float_dispatch(const float* values, int64_t count,
                                                  int64_t iterations) {
    printf("  MinMax float (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        float mn, mx;
        carquet_dispatch_minmax_float(values, count, &mn, &mx);
        g_sink = (int64_t)(mn + mx);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_minmax_double_dispatch(const double* values, int64_t count,
                                                   int64_t iterations) {
    printf("  MinMax double (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        double mn, mx;
        carquet_dispatch_minmax_double(values, count, &mn, &mx);
        g_sink = (int64_t)(mn + mx);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_copy_minmax_i32(const int32_t* values, int64_t count,
                                            int32_t* output, int64_t iterations) {
    printf("  CopyMinMax i32 (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int32_t mn, mx;
        carquet_dispatch_copy_minmax_i32(values, count, output, &mn, &mx);
        g_sink = mn + mx;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_copy_minmax_i64(const int64_t* values, int64_t count,
                                            int64_t* output, int64_t iterations) {
    printf("  CopyMinMax i64 (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t mn, mx;
        carquet_dispatch_copy_minmax_i64(values, count, output, &mn, &mx);
        g_sink = mn + mx;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_copy_minmax_float(const float* values, int64_t count,
                                              float* output, int64_t iterations) {
    printf("  CopyMinMax float (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        float mn, mx;
        carquet_dispatch_copy_minmax_float(values, count, output, &mn, &mx);
        g_sink = (int64_t)(mn + mx);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_copy_minmax_double(const double* values, int64_t count,
                                               double* output, int64_t iterations) {
    printf("  CopyMinMax double (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        double mn, mx;
        carquet_dispatch_copy_minmax_double(values, count, output, &mn, &mx);
        g_sink = (int64_t)(mn + mx);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

static void run_minmax_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== MinMax Benchmarks (Statistics) ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    int32_t* i32_data = malloc(count * sizeof(int32_t));
    int64_t* i64_data = malloc(count * sizeof(int64_t));
    float* float_data = malloc(count * sizeof(float));
    double* double_data = malloc(count * sizeof(double));
    int32_t* i32_output = malloc(count * sizeof(int32_t));
    int64_t* i64_output = malloc(count * sizeof(int64_t));
    float* float_output = malloc(count * sizeof(float));
    double* double_output = malloc(count * sizeof(double));
    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        i32_data[i] = (int32_t)(seed % 1000000) - 500000;
        i64_data[i] = (int64_t)i32_data[i] * 1000;
        float_data[i] = (float)i32_data[i] * 0.1f;
        double_data[i] = (double)i32_data[i] * 0.01;
    }

    bench_minmax_i32_scalar(i32_data, count, iterations);
    bench_minmax_i32_dispatch(i32_data, count, iterations);
    bench_minmax_i64_dispatch(i64_data, count, iterations);
    bench_minmax_float_dispatch(float_data, count, iterations);
    bench_minmax_double_dispatch(double_data, count, iterations);
    bench_copy_minmax_i32(i32_data, count, i32_output, iterations);
    bench_copy_minmax_i64(i64_data, count, i64_output, iterations);
    bench_copy_minmax_float(float_data, count, float_output, iterations);
    bench_copy_minmax_double(double_data, count, double_output, iterations);

    free(i32_data);
    free(i64_data);
    free(float_data);
    free(double_data);
    free(i32_output);
    free(i64_output);
    free(float_output);
    free(double_output);
}

/* ============================================================================
 * Byte-Stream-Split Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_bss_encode_float(const float* values, int64_t count,
                                             uint8_t* output, int64_t iterations) {
    printf("  BSS encode float (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_byte_split_encode_float(values, count, output);
        g_sink = output[0];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_bss_decode_float(const uint8_t* encoded, int64_t count,
                                             float* output, int64_t iterations) {
    printf("  BSS decode float (dispatch): ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_byte_split_decode_float(encoded, count, output);
        g_sink = (int64_t)(output[0] * 100);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_bss_encode_double(const double* values, int64_t count,
                                              uint8_t* output, int64_t iterations) {
    printf("  BSS encode dbl (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_byte_split_encode_double(values, count, output);
        g_sink = output[0];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

NOINLINE static void bench_bss_decode_double(const uint8_t* encoded, int64_t count,
                                              double* output, int64_t iterations) {
    printf("  BSS decode dbl (dispatch):   ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_dispatch_byte_split_decode_double(encoded, count, output);
        g_sink = (int64_t)(output[0] * 100);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
}

static void run_bss_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Byte-Stream-Split Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    float* float_data = malloc(count * sizeof(float));
    double* double_data = malloc(count * sizeof(double));
    uint8_t* encoded_f = malloc(count * sizeof(float));
    uint8_t* encoded_d = malloc(count * sizeof(double));
    float* decoded_f = malloc(count * sizeof(float));
    double* decoded_d = malloc(count * sizeof(double));

    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        float_data[i] = (float)(seed % 10000) * 0.01f;
        double_data[i] = (double)(seed % 10000) * 0.001;
    }

    /* Pre-encode for decode benchmarks */
    carquet_dispatch_byte_split_encode_float(float_data, count, encoded_f);
    carquet_dispatch_byte_split_encode_double(double_data, count, encoded_d);

    bench_bss_encode_float(float_data, count, encoded_f, iterations);
    bench_bss_decode_float(encoded_f, count, decoded_f, iterations);
    bench_bss_encode_double(double_data, count, encoded_d, iterations);
    bench_bss_decode_double(encoded_d, count, decoded_d, iterations);

    free(float_data);
    free(double_data);
    free(encoded_f);
    free(encoded_d);
    free(decoded_f);
    free(decoded_d);
}

/* ============================================================================
 * Prefix Sum Micro-benchmark (used in delta decoding)
 * ============================================================================ */

NOINLINE static void bench_prefix_sum_i32_scalar(int32_t* values, int64_t count,
                                                   int64_t iterations) {
    printf("  PrefixSum i32 (scalar):  ");
    fflush(stdout);

    /* Save original for re-use */
    int32_t* original = malloc(count * sizeof(int32_t));
    memcpy(original, values, count * sizeof(int32_t));

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        memcpy(values, original, count * sizeof(int32_t));
        int32_t running = 0;
        for (int64_t i = 0; i < count; i++) {
            running += values[i];
            values[i] = running;
        }
        g_sink = values[count - 1];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    memcpy(values, original, count * sizeof(int32_t));
    free(original);
}

NOINLINE static void bench_prefix_sum_i32_dispatch(int32_t* values, int64_t count,
                                                     int64_t iterations) {
    printf("  PrefixSum i32 (dispatch):");
    fflush(stdout);

    int32_t* original = malloc(count * sizeof(int32_t));
    memcpy(original, values, count * sizeof(int32_t));

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        memcpy(values, original, count * sizeof(int32_t));
        carquet_dispatch_prefix_sum_i32(values, count, 0);
        g_sink = values[count - 1];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf(" %.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    memcpy(values, original, count * sizeof(int32_t));
    free(original);
}

static void run_prefix_sum_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Prefix Sum Benchmarks (Delta Decoding) ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    int32_t* data = malloc(count * sizeof(int32_t));
    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        data[i] = (int32_t)((seed >> 16) % 100) - 50; /* small deltas */
    }

    bench_prefix_sum_i32_scalar(data, count, iterations);
    bench_prefix_sum_i32_dispatch(data, count, iterations);

    free(data);
}

/* ============================================================================
 * CRC32 Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_crc32(const uint8_t* data, size_t size,
                                  int64_t iterations) {
    printf("  CRC32 (IEEE, zlib):      ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        uint32_t crc = carquet_crc32(data, size);
        g_sink = crc;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

NOINLINE static void bench_xxhash64(const uint8_t* data, size_t size,
                                     int64_t iterations) {
    printf("  XXHash64:                ");
    fflush(stdout);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        uint64_t hash = carquet_xxhash64(data, size, 0);
        g_sink = (int64_t)hash;
    }
    double elapsed = BENCH_END();

    double mb_per_sec = (double)size * iterations / elapsed * 1e3;
    printf("%.2f MB/sec\n", mb_per_sec);
}

static void run_hash_benchmarks(size_t size, int64_t iterations) {
    printf("\n=== Hash / CRC Benchmarks ===\n");
    printf("Data size: %zu bytes, Iterations: %ld\n\n", size, (long)iterations);

    uint8_t* data = malloc(size);
    uint32_t seed = 42;
    for (size_t i = 0; i < size; i++) {
        seed = seed * 1103515245 + 12345;
        data[i] = (uint8_t)(seed >> 16);
    }

    bench_crc32(data, size, iterations);
    bench_xxhash64(data, size, iterations);

    free(data);
}

/* ============================================================================
 * Bloom Filter Micro-benchmark
 * ============================================================================ */

NOINLINE static void bench_bloom_insert(int64_t count, int64_t iterations) {
    printf("  Bloom insert i32:        ");
    fflush(stdout);

    carquet_bloom_filter_t* filter = carquet_bloom_filter_create_with_ndv(count, 0.01);
    if (!filter) {
        printf("SKIP (create failed)\n");
        return;
    }

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        for (int64_t i = 0; i < count; i++) {
            carquet_bloom_filter_insert_i32(filter, (int32_t)i);
        }
        g_sink = count;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    carquet_bloom_filter_destroy(filter);
}

NOINLINE static void bench_bloom_check(int64_t count, int64_t iterations) {
    printf("  Bloom check i32:         ");
    fflush(stdout);

    carquet_bloom_filter_t* filter = carquet_bloom_filter_create_with_ndv(count, 0.01);
    if (!filter) {
        printf("SKIP (create failed)\n");
        return;
    }

    /* Insert values first */
    for (int64_t i = 0; i < count; i++) {
        carquet_bloom_filter_insert_i32(filter, (int32_t)i);
    }

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        int64_t found = 0;
        for (int64_t i = 0; i < count; i++) {
            if (carquet_bloom_filter_check_i32(filter, (int32_t)i))
                found++;
        }
        g_sink = found;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    carquet_bloom_filter_destroy(filter);
}

static void run_bloom_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Bloom Filter Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    bench_bloom_insert(count, iterations);
    bench_bloom_check(count, iterations);
}

/* ============================================================================
 * Encoding Micro-benchmarks (encode-side)
 * ============================================================================ */

NOINLINE static void bench_encode_plain_i32(const int32_t* values, int64_t count,
                                             int64_t iterations) {
    printf("  Plain encode i32:        ");
    fflush(stdout);

    carquet_buffer_t buf;
    carquet_buffer_init_capacity(&buf, count * sizeof(int32_t) + 64);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        buf.size = 0;
        carquet_encode_plain_int32(values, count, &buf);
        g_sink = (int64_t)buf.size;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
    carquet_buffer_destroy(&buf);
}

NOINLINE static void bench_encode_delta_i32(const int32_t* values, int64_t count,
                                             int64_t iterations) {
    printf("  Delta encode i32:        ");
    fflush(stdout);

    size_t cap = count * sizeof(int32_t) + 1024;
    uint8_t* output = malloc(cap);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t written = 0;
        carquet_delta_encode_int32(values, (int32_t)count, output, cap, &written);
        g_sink = (int64_t)written;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    /* Show compression ratio */
    size_t written = 0;
    carquet_delta_encode_int32(values, (int32_t)count, output, cap, &written);
    printf("    delta ratio: %.2fx (%ld -> %ld bytes)\n",
           (double)(count * sizeof(int32_t)) / written, (long)(count * sizeof(int32_t)), (long)written);

    free(output);
}

NOINLINE static void bench_decode_delta_i32(const uint8_t* encoded, size_t encoded_size,
                                             int64_t count, int64_t iterations) {
    printf("  Delta decode i32:        ");
    fflush(stdout);

    int32_t* output = malloc(count * sizeof(int32_t));

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t consumed = 0;
        carquet_delta_decode_int32(encoded, encoded_size, output, (int32_t)count, &consumed);
        g_sink = output[count / 2];
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
    free(output);
}

NOINLINE static void bench_encode_dict_i32(const int32_t* values, int64_t count,
                                            int64_t iterations) {
    printf("  Dict encode i32:         ");
    fflush(stdout);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init_capacity(&dict_buf, 64 * 1024);
    carquet_buffer_init_capacity(&indices_buf, count * 4 + 1024);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        dict_buf.size = 0;
        indices_buf.size = 0;
        carquet_dictionary_encode_int32(values, count, &dict_buf, &indices_buf);
        g_sink = (int64_t)dict_buf.size;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);

    /* Show dictionary stats */
    dict_buf.size = 0;
    indices_buf.size = 0;
    carquet_dictionary_encode_int32(values, count, &dict_buf, &indices_buf);
    int32_t dict_entries = (int32_t)(dict_buf.size / sizeof(int32_t));
    printf("    dict entries: %d, ratio: %.2fx\n", dict_entries,
           (double)(count * sizeof(int32_t)) / (dict_buf.size + indices_buf.size));

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
}

NOINLINE static void bench_encode_bss_float(const float* values, int64_t count,
                                             int64_t iterations) {
    printf("  BSS encode float:        ");
    fflush(stdout);

    size_t cap = count * sizeof(float) + 64;
    uint8_t* output = malloc(cap);

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        size_t written = 0;
        carquet_byte_stream_split_encode_float(values, count, output, cap, &written);
        g_sink = (int64_t)written;
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
    free(output);
}

NOINLINE static void bench_decode_bss_float(const uint8_t* encoded, size_t encoded_size,
                                             int64_t count, int64_t iterations) {
    printf("  BSS decode float:        ");
    fflush(stdout);

    float* output = malloc(count * sizeof(float));

    BENCH_START();
    for (int64_t iter = 0; iter < iterations; iter++) {
        carquet_byte_stream_split_decode_float(encoded, encoded_size, output, count);
        g_sink = (int64_t)(output[count / 2] * 100);
    }
    double elapsed = BENCH_END();

    double ns_per_value = elapsed / (iterations * count);
    printf("%.2f ns/value, %.2f M values/sec\n", ns_per_value, 1e9 / ns_per_value / 1e6);
    free(output);
}

static void run_encoding_benchmarks(int64_t count, int64_t iterations) {
    printf("\n=== Encoding Benchmarks ===\n");
    printf("Values: %ld, Iterations: %ld\n\n", (long)count, (long)iterations);

    /* Generate test data */
    int32_t* i32_data = malloc(count * sizeof(int32_t));
    int32_t* i32_lowcard = malloc(count * sizeof(int32_t)); /* low cardinality for dict */
    float* float_data = malloc(count * sizeof(float));
    uint32_t seed = 42;
    for (int64_t i = 0; i < count; i++) {
        seed = seed * 1103515245 + 12345;
        i32_data[i] = (int32_t)(seed % 1000000);
        i32_lowcard[i] = (int32_t)(seed % 1000);
        float_data[i] = (float)(seed % 10000) * 0.01f;
    }

    /* Plain encoding */
    bench_encode_plain_i32(i32_data, count, iterations);

    /* Delta encoding */
    bench_encode_delta_i32(i32_data, count, iterations);

    /* Delta decode (need to encode first) */
    size_t delta_cap = count * sizeof(int32_t) + 1024;
    uint8_t* delta_encoded = malloc(delta_cap);
    size_t delta_written = 0;
    carquet_delta_encode_int32(i32_data, (int32_t)count, delta_encoded, delta_cap, &delta_written);
    bench_decode_delta_i32(delta_encoded, delta_written, count, iterations);
    free(delta_encoded);

    printf("\n");

    /* Dictionary encoding (low cardinality) */
    printf("Low cardinality (1K unique values):\n");
    bench_encode_dict_i32(i32_lowcard, count, iterations);

    printf("High cardinality (1M unique values):\n");
    bench_encode_dict_i32(i32_data, count, iterations);

    printf("\n");

    /* BSS encode/decode */
    bench_encode_bss_float(float_data, count, iterations);

    size_t bss_cap = count * sizeof(float) + 64;
    uint8_t* bss_encoded = malloc(bss_cap);
    size_t bss_written = 0;
    carquet_byte_stream_split_encode_float(float_data, count, bss_encoded, bss_cap, &bss_written);
    bench_decode_bss_float(bss_encoded, bss_written, count, iterations);
    free(bss_encoded);

    free(i32_data);
    free(i32_lowcard);
    free(float_data);
}

/* ============================================================================
 * Main
 * ============================================================================ */

static void print_usage(const char* prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  --component NAME   Component to benchmark:\n");
    printf("                     rle, gather, null, compression, dispatch, minmax,\n");
    printf("                     bss, prefix_sum, hash, bloom, encoding, all\n");
    printf("  --count N          Number of values (default: 1000000)\n");
    printf("  --iterations N     Number of iterations (default: 100)\n");
    printf("  -h, --help         Show this help\n");
    printf("\nExamples:\n");
    printf("  %s --component encoding --iterations 50\n", prog);
    printf("  %s --component all\n", prog);
}

int main(int argc, char** argv) {
    const char* component = "all";
    int64_t count = 1000000;
    int64_t iterations = 100;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        else if (strcmp(argv[i], "--component") == 0 && i + 1 < argc) {
            component = argv[++i];
        }
        else if (strcmp(argv[i], "--count") == 0 && i + 1 < argc) {
            count = atol(argv[++i]);
        }
        else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            iterations = atol(argv[++i]);
        }
    }

    printf("=== Carquet Micro-benchmarks ===\n");
    printf("Component: %s\n", component);

    int all = strcmp(component, "all") == 0;

    if (all || strcmp(component, "rle") == 0) {
        run_rle_benchmarks(count, iterations);
        run_rle_encode_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "gather") == 0) {
        run_gather_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "null") == 0) {
        run_null_bitmap_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "compression") == 0) {
        run_compression_benchmarks(1024 * 1024, iterations / 10);  /* 1MB blocks */
    }
    if (all || strcmp(component, "dispatch") == 0) {
        run_dispatch_overhead_benchmark(iterations * 10000);
    }
    if (all || strcmp(component, "minmax") == 0) {
        run_minmax_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "bss") == 0) {
        run_bss_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "prefix_sum") == 0) {
        run_prefix_sum_benchmarks(count, iterations);
    }
    if (all || strcmp(component, "hash") == 0) {
        run_hash_benchmarks(1024 * 1024, iterations / 10);  /* 1MB blocks */
    }
    if (all || strcmp(component, "bloom") == 0) {
        run_bloom_benchmarks(count / 10, iterations);  /* smaller count for bloom */
    }
    if (all || strcmp(component, "encoding") == 0) {
        run_encoding_benchmarks(count, iterations);
    }

    printf("\nDone.\n");
    return 0;
}
