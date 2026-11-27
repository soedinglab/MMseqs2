/**
 * @file profile_write.c
 * @brief Comprehensive write path profiler
 *
 * Exercises all critical write paths:
 * - Plain, dictionary, delta, byte-stream-split encoding
 * - RLE level encoding
 * - Statistics computation (SIMD minmax)
 * - Compression codecs
 * - Bloom filter insertion
 * - CRC32 computation
 * - Page index building
 *
 * Build with debug symbols (-g) and optimization (-O2) for profiling.
 *
 * Usage:
 *   ./profile_write [options]
 *   sample ./profile_write -r 5000000          # macOS
 *   perf record -g ./profile_write -r 5000000  # Linux
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
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "carquet/carquet.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define DEFAULT_NUM_ROWS      10000000   /* 10M rows */
#define DEFAULT_ROW_GROUP     1000000    /* 1M rows per group */
#define DEFAULT_ITERATIONS    5
#define DEFAULT_WARMUP        1

typedef struct {
    int64_t num_rows;
    int64_t row_group_size;
    int iterations;
    int warmup;
    int use_dictionary;      /* 0=no, 1=yes */
    int use_nulls;           /* 0=no, 1=sparse (10%), 2=moderate (30%), 3=heavy (50%) */
    int compression;         /* 0=none, 1=snappy, 2=zstd, 3=lz4, 4=gzip */
    int write_statistics;    /* 0=no, 1=yes */
    int write_bloom;         /* 0=no, 1=yes */
    int write_page_index;    /* 0=no, 1=yes */
    int write_crc;           /* 0=no, 1=yes */
    int num_columns;         /* 5, 10, 20, 50 */
    int verbose;
    char output_file[512];
} profile_config_t;

/* ============================================================================
 * Timing
 * ============================================================================ */

#ifdef _WIN32
static double get_time_ms(void) {
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart * 1000.0 / (double)freq.QuadPart;
}
#else
static double get_time_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
}
#endif

static int64_t get_file_size(const char* path) {
#ifdef _WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    if (GetFileAttributesExA(path, GetFileExInfoStandard, &fad)) {
        return ((int64_t)fad.nFileSizeHigh << 32) | fad.nFileSizeLow;
    }
    return 0;
#else
    struct stat st;
    if (stat(path, &st) == 0) {
        return st.st_size;
    }
    return 0;
#endif
}

/* ============================================================================
 * Data Generation
 * ============================================================================ */

static uint64_t g_lcg_state = 42;

static void lcg_seed(uint64_t seed) { g_lcg_state = seed; }

static uint32_t lcg_rand(void) {
    g_lcg_state = g_lcg_state * 6364136223846793005ULL + 1;
    return (uint32_t)(g_lcg_state >> 32);
}

static double lcg_uniform(void) {
    return (double)lcg_rand() / (double)UINT32_MAX;
}

static double lcg_normal(double mean, double stddev) {
    double u1 = lcg_uniform();
    double u2 = lcg_uniform();
    if (u1 < 1e-10) u1 = 1e-10;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + z * stddev;
}

typedef struct {
    int64_t* int64_data;
    double* double_data;
    int32_t* int32_data;
    float* float_data;
    int16_t* def_levels;
    int64_t count;
} test_data_t;

static void generate_test_data(test_data_t* data, int64_t count, double null_ratio,
                                int dictionary_friendly) {
    lcg_seed(42);

    data->count = count;
    data->int64_data = malloc(count * sizeof(int64_t));
    data->double_data = malloc(count * sizeof(double));
    data->int32_data = malloc(count * sizeof(int32_t));
    data->float_data = malloc(count * sizeof(float));
    data->def_levels = null_ratio > 0 ? malloc(count * sizeof(int16_t)) : NULL;

    for (int64_t i = 0; i < count; i++) {
        if (dictionary_friendly) {
            /* Low cardinality data - exercises dictionary encoding path */
            data->int64_data[i] = (int64_t)(lcg_rand() % 1000);
            data->double_data[i] = (double)(lcg_rand() % 500) + 0.5;
            data->int32_data[i] = (int32_t)(lcg_rand() % 2000);
            data->float_data[i] = (float)(lcg_rand() % 250) + 0.25f;
        } else {
            /* High cardinality - exercises plain/delta encoding */
            data->int64_data[i] = (int64_t)(lcg_normal(50000.0, 25000.0));
            data->double_data[i] = lcg_normal(1000.0, 200.0);
            data->int32_data[i] = (int32_t)(lcg_rand() % 1000000);
            data->float_data[i] = (float)lcg_normal(100.0, 50.0);
        }

        if (data->def_levels) {
            data->def_levels[i] = (lcg_uniform() >= null_ratio) ? 1 : 0;
        }
    }
}

static void free_test_data(test_data_t* data) {
    free(data->int64_data);
    free(data->double_data);
    free(data->int32_data);
    free(data->float_data);
    free(data->def_levels);
    memset(data, 0, sizeof(*data));
}

/* ============================================================================
 * Schema Creation
 * ============================================================================ */

static carquet_schema_t* create_schema(int nullable, int num_columns, carquet_error_t* err) {
    carquet_schema_t* schema = carquet_schema_create(err);
    if (!schema) return NULL;

    carquet_field_repetition_t rep = nullable ? CARQUET_REPETITION_OPTIONAL : CARQUET_REPETITION_REQUIRED;

    /* Always add id column as required */
    if (carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    /* Add columns in a rotating pattern of types */
    char colname[64];
    for (int i = 1; i < num_columns; i++) {
        int type_idx = (i - 1) % 4;
        carquet_physical_type_t ptype;
        switch (type_idx) {
            case 0: ptype = CARQUET_PHYSICAL_INT64; break;
            case 1: ptype = CARQUET_PHYSICAL_DOUBLE; break;
            case 2: ptype = CARQUET_PHYSICAL_INT32; break;
            default: ptype = CARQUET_PHYSICAL_FLOAT; break;
        }
        snprintf(colname, sizeof(colname), "col_%d", i);
        if (carquet_schema_add_column(schema, colname, ptype, NULL, rep, 0, 0) != CARQUET_OK) {
            carquet_schema_free(schema);
            return NULL;
        }
    }

    return schema;
}

/* ============================================================================
 * Write Benchmark
 * ============================================================================ */

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

static carquet_compression_t get_compression(int mode) {
    switch (mode) {
        case 1: return CARQUET_COMPRESSION_SNAPPY;
        case 2: return CARQUET_COMPRESSION_ZSTD;
        case 3: return CARQUET_COMPRESSION_LZ4;
        case 4: return CARQUET_COMPRESSION_GZIP;
        default: return CARQUET_COMPRESSION_UNCOMPRESSED;
    }
}

NOINLINE static double write_test_file(const profile_config_t* config,
                                        const test_data_t* data,
                                        const char* filename,
                                        carquet_error_t* err) {
    carquet_schema_t* schema = create_schema(data->def_levels != NULL,
                                              config->num_columns, err);
    if (!schema) return -1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = get_compression(config->compression);
    opts.row_group_size = config->row_group_size;
    opts.write_statistics = config->write_statistics;
    opts.write_crc = config->write_crc;
    opts.write_page_index = config->write_page_index;
    opts.write_bloom_filters = config->write_bloom;

    if (!config->use_dictionary) {
        opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    }

    double start = get_time_ms();

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, err);
    if (!writer) {
        carquet_schema_free(schema);
        return -1;
    }

    int64_t chunk_size = config->row_group_size;
    int64_t remaining = data->count;
    int64_t offset = 0;

    while (remaining > 0) {
        int64_t to_write = remaining < chunk_size ? remaining : chunk_size;
        int16_t* def = data->def_levels ? data->def_levels + offset : NULL;

        /* Column 0: id (always required) */
        int64_t* ids = malloc(to_write * sizeof(int64_t));
        for (int64_t i = 0; i < to_write; i++) {
            ids[i] = offset + i;
        }
        if (carquet_writer_write_batch(writer, 0, ids, to_write, NULL, NULL) != CARQUET_OK) {
            free(ids);
            carquet_status_t close_st = carquet_writer_close(writer);
            (void)close_st;
            carquet_schema_free(schema);
            return -1;
        }
        free(ids);

        /* Write remaining columns in type rotation */
        for (int col = 1; col < config->num_columns; col++) {
            const void* col_data;
            int type_idx = (col - 1) % 4;
            switch (type_idx) {
                case 0: col_data = data->int64_data + offset; break;
                case 1: col_data = data->double_data + offset; break;
                case 2: col_data = data->int32_data + offset; break;
                default: col_data = data->float_data + offset; break;
            }

            if (carquet_writer_write_batch(writer, col, col_data, to_write, def, NULL) != CARQUET_OK) {
                carquet_status_t close_st = carquet_writer_close(writer);
                (void)close_st;
                carquet_schema_free(schema);
                return -1;
            }
        }

        offset += to_write;
        remaining -= to_write;
    }

    carquet_status_t status = carquet_writer_close(writer);
    double elapsed = get_time_ms() - start;

    carquet_schema_free(schema);
    return status == CARQUET_OK ? elapsed : -1;
}

/* ============================================================================
 * Feature Overhead Measurement
 * ============================================================================ */

typedef struct {
    double baseline_ms;       /* No optional features */
    double with_stats_ms;     /* +statistics */
    double with_crc_ms;       /* +CRC32 */
    double with_bloom_ms;     /* +bloom filters */
    double with_pageindex_ms; /* +page index */
    double with_all_ms;       /* All features */
} overhead_result_t;

NOINLINE static void measure_feature_overhead(const profile_config_t* base_config,
                                               const test_data_t* data,
                                               overhead_result_t* result) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    profile_config_t config = *base_config;
    char filename[512];

    /* Baseline: no optional features */
    config.write_statistics = 0;
    config.write_crc = 0;
    config.write_bloom = 0;
    config.write_page_index = 0;
    snprintf(filename, sizeof(filename), "%s_baseline", base_config->output_file);
    result->baseline_ms = write_test_file(&config, data, filename, &err);
    remove(filename);

    /* +statistics */
    config.write_statistics = 1;
    snprintf(filename, sizeof(filename), "%s_stats", base_config->output_file);
    result->with_stats_ms = write_test_file(&config, data, filename, &err);
    remove(filename);

    /* +CRC32 */
    config.write_statistics = 0;
    config.write_crc = 1;
    snprintf(filename, sizeof(filename), "%s_crc", base_config->output_file);
    result->with_crc_ms = write_test_file(&config, data, filename, &err);
    remove(filename);

    /* +bloom filters */
    config.write_crc = 0;
    config.write_bloom = 1;
    snprintf(filename, sizeof(filename), "%s_bloom", base_config->output_file);
    result->with_bloom_ms = write_test_file(&config, data, filename, &err);
    remove(filename);

    /* +page index */
    config.write_bloom = 0;
    config.write_page_index = 1;
    snprintf(filename, sizeof(filename), "%s_pidx", base_config->output_file);
    result->with_pageindex_ms = write_test_file(&config, data, filename, &err);
    remove(filename);

    /* All features */
    config.write_statistics = 1;
    config.write_crc = 1;
    config.write_bloom = 1;
    config.write_page_index = 1;
    snprintf(filename, sizeof(filename), "%s_all", base_config->output_file);
    result->with_all_ms = write_test_file(&config, data, filename, &err);
    remove(filename);
}

/* ============================================================================
 * Encoding Comparison
 * ============================================================================ */

typedef struct {
    double plain_ms;
    double dict_ms;
    int64_t plain_size;
    int64_t dict_size;
} encoding_result_t;

NOINLINE static void measure_encoding_comparison(const profile_config_t* base_config,
                                                   const test_data_t* data,
                                                   encoding_result_t* result) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    profile_config_t config = *base_config;
    char filename[512];

    /* Plain encoding */
    config.use_dictionary = 0;
    snprintf(filename, sizeof(filename), "%s_plain", base_config->output_file);
    result->plain_ms = write_test_file(&config, data, filename, &err);
    result->plain_size = get_file_size(filename);
    remove(filename);

    /* Dictionary encoding */
    config.use_dictionary = 1;
    snprintf(filename, sizeof(filename), "%s_dict", base_config->output_file);
    result->dict_ms = write_test_file(&config, data, filename, &err);
    result->dict_size = get_file_size(filename);
    remove(filename);
}

/* ============================================================================
 * Compression Comparison
 * ============================================================================ */

typedef struct {
    double time_ms;
    int64_t file_size;
} codec_result_t;

NOINLINE static void measure_compression_comparison(const profile_config_t* base_config,
                                                      const test_data_t* data) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    profile_config_t config = *base_config;
    char filename[512];

    const char* codec_names[] = {"none", "snappy", "zstd", "lz4", "gzip"};
    codec_result_t results[5];

    for (int c = 0; c <= 4; c++) {
        config.compression = c;
        snprintf(filename, sizeof(filename), "%s_codec%d", base_config->output_file, c);
        results[c].time_ms = write_test_file(&config, data, filename, &err);
        results[c].file_size = get_file_size(filename);
        remove(filename);
    }

    printf("\n=== Compression Codec Comparison ===\n");
    printf("  %-10s %10s %10s %10s %10s\n", "Codec", "Time (ms)", "Size (MB)", "Ratio", "Throughput");
    printf("  %-10s %10s %10s %10s %10s\n", "-----", "---------", "--------", "-----", "----------");

    int64_t uncompressed_size = results[0].file_size;
    for (int c = 0; c <= 4; c++) {
        if (results[c].time_ms < 0) continue;
        double size_mb = results[c].file_size / (1024.0 * 1024.0);
        double ratio = uncompressed_size > 0 ? (double)uncompressed_size / results[c].file_size : 1.0;
        double throughput = (base_config->num_rows / results[c].time_ms) * 1000.0 / 1e6;
        printf("  %-10s %10.2f %10.2f %9.2fx %8.2f M/s\n",
               codec_names[c], results[c].time_ms, size_mb, ratio, throughput);
    }
}

/* ============================================================================
 * Column Count Scaling
 * ============================================================================ */

NOINLINE static void measure_column_scaling(const profile_config_t* base_config,
                                              const test_data_t* data) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    profile_config_t config = *base_config;
    char filename[512];

    int col_counts[] = {5, 10, 20, 50};
    int n_configs = 4;

    printf("\n=== Column Count Scaling ===\n");
    printf("  %-10s %10s %10s %10s\n", "Columns", "Time (ms)", "Size (MB)", "ms/column");
    printf("  %-10s %10s %10s %10s\n", "-------", "---------", "--------", "---------");

    for (int i = 0; i < n_configs; i++) {
        config.num_columns = col_counts[i];
        snprintf(filename, sizeof(filename), "%s_cols%d", base_config->output_file, col_counts[i]);
        double time_ms = write_test_file(&config, data, filename, &err);
        int64_t fsize = get_file_size(filename);
        remove(filename);

        if (time_ms > 0) {
            printf("  %-10d %10.2f %10.2f %10.2f\n",
                   col_counts[i], time_ms, fsize / (1024.0 * 1024.0),
                   time_ms / col_counts[i]);
        }
    }
}

/* ============================================================================
 * Main
 * ============================================================================ */

static void print_usage(const char* prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  -r, --rows N          Number of rows (default: %d)\n", DEFAULT_NUM_ROWS);
    printf("  -g, --rowgroup N      Row group size (default: %d)\n", DEFAULT_ROW_GROUP);
    printf("  -i, --iterations N    Number of iterations (default: %d)\n", DEFAULT_ITERATIONS);
    printf("  -w, --warmup N        Warmup iterations (default: %d)\n", DEFAULT_WARMUP);
    printf("  -d, --dictionary      Enable dictionary encoding\n");
    printf("  -n, --nulls MODE      Null ratio: 0=none, 1=10%%, 2=30%%, 3=50%%\n");
    printf("  -c, --compression N   0=none, 1=snappy, 2=zstd, 3=lz4, 4=gzip\n");
    printf("  -C, --columns N       Number of columns (default: 5)\n");
    printf("  -S, --no-stats        Disable statistics\n");
    printf("  -B, --bloom           Enable bloom filters\n");
    printf("  -P, --page-index      Enable page index\n");
    printf("      --no-crc          Disable CRC32\n");
    printf("      --overhead        Measure feature overhead\n");
    printf("      --codecs          Compare compression codecs\n");
    printf("      --scaling         Measure column count scaling\n");
    printf("      --encoding        Compare plain vs dictionary encoding\n");
    printf("  -o, --output FILE     Output file path\n");
    printf("  -v, --verbose         Verbose output\n");
    printf("  -h, --help            Show this help\n");
}

static double null_ratio(int mode) {
    switch (mode) {
        case 1: return 0.10;
        case 2: return 0.30;
        case 3: return 0.50;
        default: return 0.0;
    }
}

static const char* compression_name(int mode) {
    switch (mode) {
        case 1: return "snappy";
        case 2: return "zstd";
        case 3: return "lz4";
        case 4: return "gzip";
        default: return "none";
    }
}

int main(int argc, char** argv) {
    profile_config_t config;
    config.num_rows = DEFAULT_NUM_ROWS;
    config.row_group_size = DEFAULT_ROW_GROUP;
    config.iterations = DEFAULT_ITERATIONS;
    config.warmup = DEFAULT_WARMUP;
    config.use_dictionary = 0;
    config.use_nulls = 0;
    config.compression = 0;
    config.write_statistics = 1;
    config.write_bloom = 0;
    config.write_page_index = 0;
    config.write_crc = 1;
    config.num_columns = 5;
    config.verbose = 0;
    snprintf(config.output_file, sizeof(config.output_file),
             "/tmp/carquet_profile_write_%d.parquet", getpid());

    int do_overhead = 0;
    int do_codecs = 0;
    int do_scaling = 0;
    int do_encoding = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        else if ((strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--rows") == 0) && i + 1 < argc)
            config.num_rows = atol(argv[++i]);
        else if ((strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--rowgroup") == 0) && i + 1 < argc)
            config.row_group_size = atol(argv[++i]);
        else if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--iterations") == 0) && i + 1 < argc)
            config.iterations = atoi(argv[++i]);
        else if ((strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--warmup") == 0) && i + 1 < argc)
            config.warmup = atoi(argv[++i]);
        else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--dictionary") == 0)
            config.use_dictionary = 1;
        else if ((strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nulls") == 0) && i + 1 < argc)
            config.use_nulls = atoi(argv[++i]);
        else if ((strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--compression") == 0) && i + 1 < argc)
            config.compression = atoi(argv[++i]);
        else if ((strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--columns") == 0) && i + 1 < argc)
            config.num_columns = atoi(argv[++i]);
        else if (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--no-stats") == 0)
            config.write_statistics = 0;
        else if (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--bloom") == 0)
            config.write_bloom = 1;
        else if (strcmp(argv[i], "-P") == 0 || strcmp(argv[i], "--page-index") == 0)
            config.write_page_index = 1;
        else if (strcmp(argv[i], "--no-crc") == 0)
            config.write_crc = 0;
        else if (strcmp(argv[i], "--overhead") == 0)
            do_overhead = 1;
        else if (strcmp(argv[i], "--codecs") == 0)
            do_codecs = 1;
        else if (strcmp(argv[i], "--scaling") == 0)
            do_scaling = 1;
        else if (strcmp(argv[i], "--encoding") == 0)
            do_encoding = 1;
        else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc)
            strncpy(config.output_file, argv[++i], sizeof(config.output_file) - 1);
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0)
            config.verbose = 1;
    }

    printf("=== Carquet Write Path Profiler ===\n\n");
    printf("Configuration:\n");
    printf("  Rows:        %ld\n", (long)config.num_rows);
    printf("  Row groups:  %ld\n", (long)config.row_group_size);
    printf("  Columns:     %d\n", config.num_columns);
    printf("  Iterations:  %d (+%d warmup)\n", config.iterations, config.warmup);
    printf("  Dictionary:  %s\n", config.use_dictionary ? "yes" : "no");
    printf("  Nulls:       %.0f%%\n", null_ratio(config.use_nulls) * 100);
    printf("  Compression: %s\n", compression_name(config.compression));
    printf("  Statistics:  %s\n", config.write_statistics ? "yes" : "no");
    printf("  CRC32:       %s\n", config.write_crc ? "yes" : "no");
    printf("  Bloom:       %s\n", config.write_bloom ? "yes" : "no");
    printf("  Page index:  %s\n", config.write_page_index ? "yes" : "no");
    printf("  Output file: %s\n", config.output_file);
    printf("\n");

    /* Generate test data */
    printf("Generating %ld rows of test data...\n", (long)config.num_rows);
    test_data_t data;
    generate_test_data(&data, config.num_rows, null_ratio(config.use_nulls),
                       config.use_dictionary);

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* === Standard write benchmark === */
    printf("\nWrite benchmark phase:\n");

    /* Warmup */
    for (int w = 0; w < config.warmup; w++) {
        double t = write_test_file(&config, &data, config.output_file, &err);
        if (t < 0) {
            fprintf(stderr, "Write failed: %s\n", err.message);
            free_test_data(&data);
            return 1;
        }
        remove(config.output_file);
        if (config.verbose) {
            printf("  Warmup %d: %.2f ms\n", w + 1, t);
        }
    }

    /* Measured iterations */
    double total_time = 0;
    double min_time = 1e9;
    double max_time = 0;
    int64_t last_file_size = 0;

    for (int i = 0; i < config.iterations; i++) {
        double t = write_test_file(&config, &data, config.output_file, &err);
        if (t < 0) {
            fprintf(stderr, "Write failed: %s\n", err.message);
            free_test_data(&data);
            return 1;
        }

        last_file_size = get_file_size(config.output_file);
        remove(config.output_file);

        total_time += t;
        if (t < min_time) min_time = t;
        if (t > max_time) max_time = t;

        if (config.verbose) {
            printf("  Iteration %d: %.2f ms (%.2f M rows/sec)\n",
                   i + 1, t, (config.num_rows / t) * 1000.0 / 1e6);
        }
    }

    double avg_time = total_time / config.iterations;
    double throughput = (config.num_rows / avg_time) * 1000.0 / 1e6;

    printf("\n=== Results ===\n");
    printf("  Avg write time: %.2f ms\n", avg_time);
    printf("  Min write time: %.2f ms\n", min_time);
    printf("  Max write time: %.2f ms\n", max_time);
    printf("  Throughput:     %.2f M rows/sec\n", throughput);
    printf("  File size:      %.2f MB\n", last_file_size / 1024.0 / 1024.0);
    printf("  Per-column:     %.2f ms/col\n", avg_time / config.num_columns);

    /* === Optional analysis modes === */

    if (do_overhead) {
        printf("\n=== Feature Overhead Analysis ===\n");
        overhead_result_t overhead;
        measure_feature_overhead(&config, &data, &overhead);

        printf("  %-20s %10s %10s\n", "Feature", "Time (ms)", "Overhead");
        printf("  %-20s %10s %10s\n", "-------", "---------", "--------");
        printf("  %-20s %10.2f %10s\n", "baseline", overhead.baseline_ms, "-");

        if (overhead.baseline_ms > 0) {
            double base = overhead.baseline_ms;
            if (overhead.with_stats_ms > 0)
                printf("  %-20s %10.2f %+9.1f%%\n", "+statistics", overhead.with_stats_ms,
                       (overhead.with_stats_ms - base) / base * 100);
            if (overhead.with_crc_ms > 0)
                printf("  %-20s %10.2f %+9.1f%%\n", "+CRC32", overhead.with_crc_ms,
                       (overhead.with_crc_ms - base) / base * 100);
            if (overhead.with_bloom_ms > 0)
                printf("  %-20s %10.2f %+9.1f%%\n", "+bloom filters", overhead.with_bloom_ms,
                       (overhead.with_bloom_ms - base) / base * 100);
            if (overhead.with_pageindex_ms > 0)
                printf("  %-20s %10.2f %+9.1f%%\n", "+page index", overhead.with_pageindex_ms,
                       (overhead.with_pageindex_ms - base) / base * 100);
            if (overhead.with_all_ms > 0)
                printf("  %-20s %10.2f %+9.1f%%\n", "all features", overhead.with_all_ms,
                       (overhead.with_all_ms - base) / base * 100);
        }
    }

    if (do_encoding) {
        /* Test with dictionary-friendly data */
        test_data_t dict_data;
        printf("\n=== Encoding Comparison (low cardinality data) ===\n");
        generate_test_data(&dict_data, config.num_rows, null_ratio(config.use_nulls), 1);
        encoding_result_t enc_result;
        measure_encoding_comparison(&config, &dict_data, &enc_result);
        printf("  Plain:      %.2f ms  (%.2f MB)\n", enc_result.plain_ms,
               enc_result.plain_size / (1024.0 * 1024.0));
        printf("  Dictionary: %.2f ms  (%.2f MB)\n", enc_result.dict_ms,
               enc_result.dict_size / (1024.0 * 1024.0));
        if (enc_result.plain_ms > 0 && enc_result.dict_ms > 0)
            printf("  Speedup:    %.2fx\n", enc_result.plain_ms / enc_result.dict_ms);
        if (enc_result.plain_size > 0 && enc_result.dict_size > 0)
            printf("  Size ratio: %.2fx\n", (double)enc_result.plain_size / enc_result.dict_size);
        free_test_data(&dict_data);

        /* Also test with high cardinality data */
        test_data_t hc_data;
        printf("\n=== Encoding Comparison (high cardinality data) ===\n");
        generate_test_data(&hc_data, config.num_rows, null_ratio(config.use_nulls), 0);
        measure_encoding_comparison(&config, &hc_data, &enc_result);
        printf("  Plain:      %.2f ms  (%.2f MB)\n", enc_result.plain_ms,
               enc_result.plain_size / (1024.0 * 1024.0));
        printf("  Dictionary: %.2f ms  (%.2f MB)\n", enc_result.dict_ms,
               enc_result.dict_size / (1024.0 * 1024.0));
        if (enc_result.plain_ms > 0 && enc_result.dict_ms > 0)
            printf("  Speedup:    %.2fx\n", enc_result.plain_ms / enc_result.dict_ms);
        free_test_data(&hc_data);
    }

    if (do_codecs) {
        measure_compression_comparison(&config, &data);
    }

    if (do_scaling) {
        measure_column_scaling(&config, &data);
    }

    /* Cleanup */
    free_test_data(&data);

    printf("\nDone.\n");
    return 0;
}
