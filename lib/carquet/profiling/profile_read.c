/**
 * @file profile_read.c
 * @brief Comprehensive read path profiler for perf analysis
 *
 * Designed to exercise all critical read paths:
 * - Dictionary encoding with gather operations
 * - RLE level decoding
 * - Null bitmap construction
 * - Various compression codecs
 * - SIMD dispatch paths
 *
 * Build with debug symbols (-g) and optimization (-O2) for perf profiling.
 *
 * Usage:
 *   perf record -g ./profile_read [options]
 *   perf report
 *   perf annotate carquet_rle_decoder_get
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
#define DEFAULT_BATCH_SIZE    262144     /* 256K - same as benchmark */
#define DEFAULT_ROW_GROUP     1000000    /* 1M rows per group */
#define DEFAULT_ITERATIONS    10
#define DEFAULT_WARMUP        2

typedef struct {
    int64_t num_rows;
    int64_t batch_size;
    int64_t row_group_size;
    int iterations;
    int warmup;
    int use_dictionary;      /* 0=no, 1=yes, 2=auto */
    int use_nulls;           /* 0=no, 1=sparse (10%), 2=moderate (30%), 3=heavy (50%) */
    int compression;         /* 0=none, 1=snappy, 2=zstd, 3=lz4 */
    int verbose;
    int profile_mode;        /* 0=full, 1=read_only, 2=write_only */
    char output_file[512];
} profile_config_t;

/* ============================================================================
 * Timing Infrastructure
 * ============================================================================ */

typedef struct {
    double write_time_ms;
    double read_time_ms;
    double dict_gather_time_ms;
    double rle_decode_time_ms;
    double decompress_time_ms;
    int64_t bytes_read;
    int64_t values_read;
    int64_t file_size;
} profile_result_t;

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
 * Data Generation (matches benchmark for reproducibility)
 * ============================================================================ */

static uint64_t g_lcg_state = 42;

static void lcg_seed(uint64_t seed) {
    g_lcg_state = seed;
}

static uint32_t lcg_rand(void) {
    g_lcg_state = g_lcg_state * 6364136223846793005ULL + 1;
    return (uint32_t)(g_lcg_state >> 32);
}

static double lcg_uniform(void) {
    return (double)lcg_rand() / (double)UINT32_MAX;
}

/* Box-Muller transform for normal distribution */
static double lcg_normal(double mean, double stddev) {
    double u1 = lcg_uniform();
    double u2 = lcg_uniform();
    if (u1 < 1e-10) u1 = 1e-10;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + z * stddev;
}

/* Generate data with configurable patterns */
typedef struct {
    int64_t* int64_data;
    double* double_data;
    int32_t* int32_data;
    float* float_data;
    int16_t* def_levels;
    int64_t count;
    double null_ratio;
} test_data_t;

static void generate_test_data(test_data_t* data, int64_t count, double null_ratio) {
    lcg_seed(42);

    data->count = count;
    data->null_ratio = null_ratio;

    data->int64_data = malloc(count * sizeof(int64_t));
    data->double_data = malloc(count * sizeof(double));
    data->int32_data = malloc(count * sizeof(int32_t));
    data->float_data = malloc(count * sizeof(float));
    data->def_levels = null_ratio > 0 ? malloc(count * sizeof(int16_t)) : NULL;

    for (int64_t i = 0; i < count; i++) {
        /* Non-sequential patterns to avoid trivial caching */
        data->int64_data[i] = (int64_t)(lcg_normal(50000.0, 25000.0));
        data->double_data[i] = lcg_normal(1000.0, 200.0);
        data->int32_data[i] = (int32_t)(lcg_rand() % 1000000);
        data->float_data[i] = (float)lcg_normal(100.0, 50.0);

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

static carquet_schema_t* create_test_schema(int nullable, carquet_error_t* err) {
    carquet_schema_t* schema = carquet_schema_create(err);
    if (!schema) return NULL;

    carquet_field_repetition_t rep = nullable ? CARQUET_REPETITION_OPTIONAL : CARQUET_REPETITION_REQUIRED;

    /* Comprehensive set of column types to exercise all paths */
    if (carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "int64_col", CARQUET_PHYSICAL_INT64, NULL, rep, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "double_col", CARQUET_PHYSICAL_DOUBLE, NULL, rep, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "int32_col", CARQUET_PHYSICAL_INT32, NULL, rep, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "float_col", CARQUET_PHYSICAL_FLOAT, NULL, rep, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    return schema;
}

/* ============================================================================
 * Write Phase
 * ============================================================================ */

static carquet_compression_t get_compression(int compression_mode) {
    switch (compression_mode) {
        case 1: return CARQUET_COMPRESSION_SNAPPY;
        case 2: return CARQUET_COMPRESSION_ZSTD;
        case 3: return CARQUET_COMPRESSION_LZ4;
        default: return CARQUET_COMPRESSION_UNCOMPRESSED;
    }
}

static double write_test_file(const profile_config_t* config, const test_data_t* data,
                              const char* filename, carquet_error_t* err) {
    carquet_schema_t* schema = create_test_schema(data->def_levels != NULL, err);
    if (!schema) return -1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = get_compression(config->compression);
    opts.row_group_size = config->row_group_size;
    /* Dictionary encoding: PLAIN to disable, PLAIN_DICTIONARY to enable */
    if (!config->use_dictionary) {
        opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    }

    double start = get_time_ms();

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, err);
    if (!writer) {
        carquet_schema_free(schema);
        return -1;
    }

    /* Write in chunks to simulate realistic workload */
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
            if (carquet_writer_close(writer) != CARQUET_OK) { /* best-effort cleanup */ }
            carquet_schema_free(schema);
            return -1;
        }
        free(ids);

        /* Columns 1-4 */
        if (carquet_writer_write_batch(writer, 1, data->int64_data + offset, to_write, def, NULL) != CARQUET_OK ||
            carquet_writer_write_batch(writer, 2, data->double_data + offset, to_write, def, NULL) != CARQUET_OK ||
            carquet_writer_write_batch(writer, 3, data->int32_data + offset, to_write, def, NULL) != CARQUET_OK ||
            carquet_writer_write_batch(writer, 4, data->float_data + offset, to_write, def, NULL) != CARQUET_OK) {
            if (carquet_writer_close(writer) != CARQUET_OK) { /* best-effort cleanup */ }
            carquet_schema_free(schema);
            return -1;
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
 * Read Phase - The critical path to profile
 * ============================================================================ */

/* Volatile to prevent optimization */
static volatile int64_t g_checksum;

/*
 * NOINLINE attribute ensures these functions appear in perf output
 * even when called frequently
 */
#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

/*
 * This is the main read loop - designed to be visible in perf annotate.
 * The inner loop structure matches what PyArrow does, making comparison fair.
 */
NOINLINE static int64_t read_batch_loop(carquet_batch_reader_t* batch_reader,
                                         const profile_config_t* config,
                                         int64_t* checksum) {
    (void)config;  /* May be used for future options */
    carquet_row_batch_t* batch = NULL;
    int64_t total_values = 0;
    int64_t local_sum = 0;

    while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        /* Process each column - this exercises the gather/decode paths */
        for (int col = 0; col < 5; col++) {
            const void* data;
            const uint8_t* null_bitmap;
            int64_t count;

            if (carquet_row_batch_column(batch, col, &data, &null_bitmap, &count) != CARQUET_OK) {
                continue;
            }

            total_values += count;

            /* Sample every 1000th value to create checksum without adding overhead */
            /* This matches the benchmark pattern for fair comparison */
            switch (col) {
                case 0: /* id - int64 */
                case 1: /* int64_col */
                    for (int64_t i = 0; i < count; i += 1000) {
                        local_sum += ((const int64_t*)data)[i];
                    }
                    break;
                case 2: /* double_col */
                    for (int64_t i = 0; i < count; i += 1000) {
                        local_sum += (int64_t)(((const double*)data)[i] * 100);
                    }
                    break;
                case 3: /* int32_col */
                    for (int64_t i = 0; i < count; i += 1000) {
                        local_sum += ((const int32_t*)data)[i];
                    }
                    break;
                case 4: /* float_col */
                    for (int64_t i = 0; i < count; i += 1000) {
                        local_sum += (int64_t)(((const float*)data)[i] * 100);
                    }
                    break;
            }
        }

        carquet_row_batch_free(batch);
        batch = NULL;
    }

    *checksum = local_sum;
    return total_values;
}

NOINLINE static double read_test_file(const profile_config_t* config,
                                       const char* filename,
                                       profile_result_t* result,
                                       carquet_error_t* err) {
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = 1;           /* Use mmap for better perf */
    opts.verify_checksums = 0;   /* Disable for pure read profiling */

    carquet_reader_t* reader = carquet_reader_open(filename, &opts, err);
    if (!reader) return -1;

    carquet_batch_reader_config_t batch_config;
    carquet_batch_reader_config_init(&batch_config);
    batch_config.batch_size = config->batch_size;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &batch_config, err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        return -1;
    }

    double start = get_time_ms();

    int64_t checksum;
    int64_t values = read_batch_loop(batch_reader, config, &checksum);

    double elapsed = get_time_ms() - start;

    /* Store volatile to prevent optimization */
    g_checksum = checksum;

    result->values_read = values;
    result->bytes_read = values * 8;  /* Approximate */

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);

    return elapsed;
}

/* ============================================================================
 * Main Profiling Loop
 * ============================================================================ */

static void print_usage(const char* prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  -r, --rows N          Number of rows (default: %d)\n", DEFAULT_NUM_ROWS);
    printf("  -b, --batch N         Batch size (default: %d)\n", DEFAULT_BATCH_SIZE);
    printf("  -g, --rowgroup N      Row group size (default: %d)\n", DEFAULT_ROW_GROUP);
    printf("  -i, --iterations N    Number of iterations (default: %d)\n", DEFAULT_ITERATIONS);
    printf("  -w, --warmup N        Warmup iterations (default: %d)\n", DEFAULT_WARMUP);
    printf("  -d, --dictionary      Enable dictionary encoding\n");
    printf("  -n, --nulls MODE      Null ratio: 0=none, 1=10%%, 2=30%%, 3=50%%\n");
    printf("  -c, --compression N   0=none, 1=snappy, 2=zstd, 3=lz4\n");
    printf("  -m, --mode MODE       0=full, 1=read-only, 2=write-only\n");
    printf("  -o, --output FILE     Output file path\n");
    printf("  -v, --verbose         Verbose output\n");
    printf("  -h, --help            Show this help\n");
    printf("\nExample:\n");
    printf("  perf record -g %s -r 1000000 -d -n 1 -c 2\n", prog);
    printf("  perf report --hierarchy\n");
}

static void parse_args(int argc, char** argv, profile_config_t* config) {
    config->num_rows = DEFAULT_NUM_ROWS;
    config->batch_size = DEFAULT_BATCH_SIZE;
    config->row_group_size = DEFAULT_ROW_GROUP;
    config->iterations = DEFAULT_ITERATIONS;
    config->warmup = DEFAULT_WARMUP;
    config->use_dictionary = 0;
    config->use_nulls = 0;
    config->compression = 0;
    config->verbose = 0;
    config->profile_mode = 0;
    snprintf(config->output_file, sizeof(config->output_file),
             "/tmp/carquet_profile_%d.parquet", getpid());

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        }
        else if ((strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--rows") == 0) && i + 1 < argc) {
            config->num_rows = atol(argv[++i]);
        }
        else if ((strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--batch") == 0) && i + 1 < argc) {
            config->batch_size = atol(argv[++i]);
        }
        else if ((strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--rowgroup") == 0) && i + 1 < argc) {
            config->row_group_size = atol(argv[++i]);
        }
        else if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--iterations") == 0) && i + 1 < argc) {
            config->iterations = atoi(argv[++i]);
        }
        else if ((strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--warmup") == 0) && i + 1 < argc) {
            config->warmup = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--dictionary") == 0) {
            config->use_dictionary = 1;
        }
        else if ((strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nulls") == 0) && i + 1 < argc) {
            config->use_nulls = atoi(argv[++i]);
        }
        else if ((strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--compression") == 0) && i + 1 < argc) {
            config->compression = atoi(argv[++i]);
        }
        else if ((strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mode") == 0) && i + 1 < argc) {
            config->profile_mode = atoi(argv[++i]);
        }
        else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            strncpy(config->output_file, argv[++i], sizeof(config->output_file) - 1);
        }
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            config->verbose = 1;
        }
    }
}

static const char* compression_name(int mode) {
    switch (mode) {
        case 1: return "snappy";
        case 2: return "zstd";
        case 3: return "lz4";
        default: return "none";
    }
}

static double null_ratio(int mode) {
    switch (mode) {
        case 1: return 0.10;
        case 2: return 0.30;
        case 3: return 0.50;
        default: return 0.0;
    }
}

int main(int argc, char** argv) {
    profile_config_t config;
    parse_args(argc, argv, &config);

    printf("=== Carquet Read Path Profiler ===\n\n");
    printf("Configuration:\n");
    printf("  Rows:        %ld\n", (long)config.num_rows);
    printf("  Batch size:  %ld\n", (long)config.batch_size);
    printf("  Row groups:  %ld\n", (long)config.row_group_size);
    printf("  Iterations:  %d (+%d warmup)\n", config.iterations, config.warmup);
    printf("  Dictionary:  %s\n", config.use_dictionary ? "yes" : "no");
    printf("  Nulls:       %.0f%%\n", null_ratio(config.use_nulls) * 100);
    printf("  Compression: %s\n", compression_name(config.compression));
    printf("  Output file: %s\n", config.output_file);
    printf("\n");

    /* Generate test data */
    printf("Generating %ld rows of test data...\n", (long)config.num_rows);
    test_data_t data;
    generate_test_data(&data, config.num_rows, null_ratio(config.use_nulls));

    carquet_error_t err = CARQUET_ERROR_INIT;
    profile_result_t result = {0};

    /* Write phase */
    if (config.profile_mode != 1) {
        printf("Writing test file...\n");
        double write_time = write_test_file(&config, &data, config.output_file, &err);
        if (write_time < 0) {
            fprintf(stderr, "Write failed: %s\n", err.message);
            free_test_data(&data);
            return 1;
        }
        result.write_time_ms = write_time;
        result.file_size = get_file_size(config.output_file);
        printf("  Write time: %.2f ms\n", write_time);
        printf("  File size:  %.2f MB\n", result.file_size / 1024.0 / 1024.0);
    }

    /* Read phase - the main profiling target */
    if (config.profile_mode != 2) {
        printf("\nRead profiling phase:\n");

        /* Warmup */
        for (int w = 0; w < config.warmup; w++) {
            double t = read_test_file(&config, config.output_file, &result, &err);
            if (t < 0) {
                fprintf(stderr, "Read failed: %s\n", err.message);
                free_test_data(&data);
                return 1;
            }
            if (config.verbose) {
                printf("  Warmup %d: %.2f ms\n", w + 1, t);
            }
        }

        /* Measured iterations */
        double total_time = 0;
        double min_time = 1e9;
        double max_time = 0;

        for (int i = 0; i < config.iterations; i++) {
            double t = read_test_file(&config, config.output_file, &result, &err);
            if (t < 0) {
                fprintf(stderr, "Read failed: %s\n", err.message);
                free_test_data(&data);
                return 1;
            }

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
        printf("  Avg read time: %.2f ms\n", avg_time);
        printf("  Min read time: %.2f ms\n", min_time);
        printf("  Max read time: %.2f ms\n", max_time);
        printf("  Throughput:    %.2f M rows/sec\n", throughput);
        printf("  Values read:   %ld per iteration\n", (long)result.values_read);
        printf("  Checksum:      %ld\n", (long)g_checksum);

        /* Output comparison reference */
        printf("\n=== PyArrow Comparison ===\n");
        printf("  PyArrow typically achieves ~50-100 M rows/sec for this workload.\n");
        printf("  Current: %.2f M rows/sec (%.1fx difference)\n",
               throughput, 75.0 / throughput);
    }

    /* Cleanup */
    free_test_data(&data);

    if (config.profile_mode != 1) {
        remove(config.output_file);
    }

    printf("\nDone. Use 'perf report' or 'perf annotate' to analyze results.\n");

    return 0;
}
