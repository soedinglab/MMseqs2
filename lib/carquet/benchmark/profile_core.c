/**
 * @file profile_core.c
 * @brief Profiling benchmark for Carquet core functions
 *
 * This benchmark is designed for CPU profiling (sample, Instruments, perf).
 * It can run in two modes:
 *
 * 1. READ-ONLY MODE (recommended for profiling):
 *    Uses pre-generated Parquet files from generate_test_files
 *    ./profile_core --read-only <path_to_profile_data>
 *
 * 2. FULL MODE (write + read):
 *    Generates data and writes/reads files
 *    ./profile_core [iterations] [rows_per_iter]
 *
 * Usage for profiling:
 *    ./generate_test_files benchmark/profile_data 10000000
 *    sample ./profile_core 30 -wait -f profile.txt &
 *    ./profile_core --read-only benchmark/profile_data --iterations 20
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <dirent.h>

#include <carquet/carquet.h>

/* Default parameters */
#define DEFAULT_ITERATIONS 5
#define DEFAULT_ROWS 1000000

/* ============================================================================
 * Read-Only Benchmark (for profiling reads)
 * ============================================================================ */

static int benchmark_read_file(const char* filepath, int iterations) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    printf("  Reading %s (%d iterations)...\n", filepath, iterations);

    for (int iter = 0; iter < iterations; iter++) {
        carquet_reader_options_t opts;
        carquet_reader_options_init(&opts);
        opts.verify_checksums = true;

        carquet_reader_t* reader = carquet_reader_open(filepath, &opts, &err);
        if (!reader) {
            fprintf(stderr, "    Failed to open: %s\n", err.message);
            return -1;
        }

        int32_t num_cols = carquet_reader_num_columns(reader);

        /* Use batch reader for realistic workload */
        carquet_batch_reader_config_t config;
        carquet_batch_reader_config_init(&config);
        config.batch_size = 65536;

        carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
        if (!batch_reader) {
            fprintf(stderr, "    Failed to create batch reader: %s\n", err.message);
            carquet_reader_close(reader);
            return -1;
        }

        int64_t rows_read = 0;
        carquet_row_batch_t* batch = NULL;

        while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
            int64_t batch_rows = carquet_row_batch_num_rows(batch);
            rows_read += batch_rows;

            /* Access each column to force full decoding */
            for (int32_t col = 0; col < carquet_row_batch_num_columns(batch); col++) {
                const void* data;
                const uint8_t* null_bitmap;
                int64_t num_values;
                if (carquet_row_batch_column(batch, col, &data, &null_bitmap, &num_values) != CARQUET_OK)
                    continue;

                /* Touch data to prevent optimization */
                if (data && num_values > 0) {
                    volatile uint8_t x = ((const uint8_t*)data)[0];
                    (void)x;
                }
            }

            carquet_row_batch_free(batch);
            batch = NULL;
        }

        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);

        if (iter == 0) {
            printf("    %lld rows, %d columns\n", (long long)rows_read, num_cols);
        }
    }

    return 0;
}

static int benchmark_column_read_file(const char* filepath, int iterations) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    printf("  Column-level read %s (%d iterations)...\n", filepath, iterations);

    for (int iter = 0; iter < iterations; iter++) {
        carquet_reader_t* reader = carquet_reader_open(filepath, NULL, &err);
        if (!reader) return -1;

        int32_t num_row_groups = carquet_reader_num_row_groups(reader);
        int32_t num_columns = carquet_reader_num_columns(reader);

        /* Allocate read buffers */
        void* buffer = malloc(65536 * sizeof(double));
        int16_t* def_levels = malloc(65536 * sizeof(int16_t));

        if (!buffer || !def_levels) {
            carquet_reader_close(reader);
            free(buffer);
            free(def_levels);
            return -1;
        }

        /* Read each column from each row group */
        for (int32_t rg = 0; rg < num_row_groups; rg++) {
            for (int32_t col = 0; col < num_columns; col++) {
                carquet_column_reader_t* col_reader =
                    carquet_reader_get_column(reader, rg, col, &err);

                if (!col_reader) continue;

                int64_t values_read;
                while ((values_read = carquet_column_read_batch(
                            col_reader, buffer, 65536, def_levels, NULL)) > 0) {
                    /* Touch data */
                    volatile uint8_t x = ((const uint8_t*)buffer)[0];
                    (void)x;
                }

                carquet_column_reader_free(col_reader);
            }
        }

        free(buffer);
        free(def_levels);
        carquet_reader_close(reader);
    }

    return 0;
}

static int run_read_only_benchmark(const char* data_dir, int iterations) {
    printf("\n=== Read-Only Profiling Mode ===\n");
    printf("Data directory: %s\n", data_dir);
    printf("Iterations per file: %d\n\n", iterations);

    const char* files[] = {
        "test_none.parquet",
        "test_snappy.parquet",
        "test_zstd.parquet",
    };
    int num_files = sizeof(files) / sizeof(files[0]);

    char filepath[1024];

    /* Batch reader benchmark */
    printf("--- Batch Reader Benchmark ---\n");
    for (int f = 0; f < num_files; f++) {
        snprintf(filepath, sizeof(filepath), "%s/%s", data_dir, files[f]);

        FILE* test = fopen(filepath, "rb");
        if (!test) {
            printf("  Skipping %s (not found)\n", files[f]);
            continue;
        }
        fclose(test);

        if (benchmark_read_file(filepath, iterations) != 0) {
            fprintf(stderr, "  FAILED: %s\n", files[f]);
        }
    }

    /* Column reader benchmark */
    printf("\n--- Column Reader Benchmark ---\n");
    for (int f = 0; f < num_files; f++) {
        snprintf(filepath, sizeof(filepath), "%s/%s", data_dir, files[f]);

        FILE* test = fopen(filepath, "rb");
        if (!test) continue;
        fclose(test);

        if (benchmark_column_read_file(filepath, iterations) != 0) {
            fprintf(stderr, "  FAILED: %s\n", files[f]);
        }
    }

    printf("\nRead-only benchmark complete.\n");
    return 0;
}

/* ============================================================================
 * Full Benchmark (Write + Read)
 * ============================================================================ */

static const char* TEST_FILE = "/tmp/carquet_profile.parquet";

static int64_t* g_int64_data = NULL;
static int32_t* g_int32_data = NULL;
static double* g_double_data = NULL;
static float* g_float_data = NULL;
static int32_t* g_low_cardinality = NULL;
static int16_t* g_def_levels = NULL;
static int64_t g_num_rows = 0;

static void generate_test_data(int64_t num_rows) {
    printf("Generating %lld rows of test data...\n", (long long)num_rows);

    g_num_rows = num_rows;

    g_int64_data = malloc(num_rows * sizeof(int64_t));
    g_int32_data = malloc(num_rows * sizeof(int32_t));
    g_double_data = malloc(num_rows * sizeof(double));
    g_float_data = malloc(num_rows * sizeof(float));
    g_low_cardinality = malloc(num_rows * sizeof(int32_t));
    g_def_levels = malloc(num_rows * sizeof(int16_t));

    if (!g_int64_data || !g_int32_data || !g_double_data ||
        !g_float_data || !g_low_cardinality || !g_def_levels) {
        fprintf(stderr, "Failed to allocate test data\n");
        exit(1);
    }

    unsigned int seed = 42;
    for (int64_t i = 0; i < num_rows; i++) {
        g_int64_data[i] = i * 1000 + (rand_r(&seed) % 100);
        g_int32_data[i] = (int32_t)(rand_r(&seed) % 1000000);
        g_double_data[i] = (double)i * 0.001 + (rand_r(&seed) % 1000) * 0.0001;
        g_float_data[i] = (float)(rand_r(&seed) % 10000) * 0.01f;
        g_low_cardinality[i] = (int32_t)(rand_r(&seed) % 100);
        g_def_levels[i] = (rand_r(&seed) % 10 == 0) ? 0 : 1;
    }

    printf("Test data generated.\n\n");
}

static void free_test_data(void) {
    free(g_int64_data);
    free(g_int32_data);
    free(g_double_data);
    free(g_float_data);
    free(g_low_cardinality);
    free(g_def_levels);
}

static carquet_schema_t* create_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return NULL;

    if (carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64,
                                  NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "value_i32", CARQUET_PHYSICAL_INT32,
                                  NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "value_f64", CARQUET_PHYSICAL_DOUBLE,
                                  NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "value_f32", CARQUET_PHYSICAL_FLOAT,
                                  NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "category", CARQUET_PHYSICAL_INT32,
                                  NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "nullable_val", CARQUET_PHYSICAL_DOUBLE,
                                  NULL, CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    return schema;
}

static int benchmark_write(carquet_compression_t compression) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = create_schema();
    if (!schema) return -1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = compression;
    opts.compression_level = 3;
    opts.row_group_size = 64 * 1024 * 1024;
    opts.write_statistics = true;

    carquet_writer_t* writer = carquet_writer_create(TEST_FILE, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        return -1;
    }

    if (carquet_writer_write_batch(writer, 0, g_int64_data, g_num_rows, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 1, g_int32_data, g_num_rows, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 2, g_double_data, g_num_rows, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 3, g_float_data, g_num_rows, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 4, g_low_cardinality, g_num_rows, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 5, g_double_data, g_num_rows, g_def_levels, NULL) != CARQUET_OK) {
        if (carquet_writer_close(writer) != CARQUET_OK) { /* best-effort cleanup */ }
        carquet_schema_free(schema);
        return -1;
    }

    carquet_status_t status = carquet_writer_close(writer);
    carquet_schema_free(schema);

    return (status == CARQUET_OK) ? 0 : -1;
}

static int benchmark_read(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.verify_checksums = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &opts, &err);
    if (!reader) return -1;

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 65536;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        return -1;
    }

    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        for (int32_t col = 0; col < carquet_row_batch_num_columns(batch); col++) {
            const void* data;
            const uint8_t* null_bitmap;
            int64_t num_values;
            if (carquet_row_batch_column(batch, col, &data, &null_bitmap, &num_values) != CARQUET_OK)
                continue;
            if (data && num_values > 0) {
                volatile uint8_t x = ((const uint8_t*)data)[0];
                (void)x;
            }
        }
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);

    return 0;
}

static int run_full_benchmark(int iterations, int64_t num_rows) {
    printf("\n=== Full Benchmark Mode (Write + Read) ===\n");
    printf("Iterations: %d\n", iterations);
    printf("Rows per iteration: %lld\n\n", (long long)num_rows);

    generate_test_data(num_rows);

    struct {
        carquet_compression_t codec;
        const char* name;
    } codecs[] = {
        { CARQUET_COMPRESSION_UNCOMPRESSED, "NONE" },
        { CARQUET_COMPRESSION_SNAPPY, "SNAPPY" },
        { CARQUET_COMPRESSION_ZSTD, "ZSTD" },
    };
    int num_codecs = sizeof(codecs) / sizeof(codecs[0]);

    for (int iter = 0; iter < iterations; iter++) {
        printf("--- Iteration %d/%d ---\n", iter + 1, iterations);

        for (int c = 0; c < num_codecs; c++) {
            printf("  [%s] Write... ", codecs[c].name);
            fflush(stdout);
            if (benchmark_write(codecs[c].codec) != 0) {
                printf("FAILED\n");
                continue;
            }
            printf("OK  ");

            printf("Read... ");
            fflush(stdout);
            if (benchmark_read() != 0) {
                printf("FAILED\n");
                continue;
            }
            printf("OK\n");
        }
        printf("\n");
    }

    free_test_data();
    remove(TEST_FILE);

    printf("Full benchmark complete.\n");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================ */

static void print_usage(const char* prog) {
    printf("Usage:\n");
    printf("  %s --read-only <data_dir> [--iterations N]\n", prog);
    printf("  %s [iterations] [rows_per_iter]\n", prog);
    printf("\nExamples:\n");
    printf("  ./generate_test_files benchmark/profile_data 10000000\n");
    printf("  %s --read-only benchmark/profile_data --iterations 20\n", prog);
    printf("  %s 5 2000000\n", prog);
}

int main(int argc, char* argv[]) {
    printf("=== Carquet Core Profiling Benchmark ===\n");

    /* Initialize library */
    if (carquet_init() != CARQUET_OK) {
        fprintf(stderr, "Failed to initialize carquet\n");
        return 1;
    }

    /* Print CPU info */
    const carquet_cpu_info_t* cpu = carquet_get_cpu_info();
    printf("CPU Features: ");
    if (cpu->has_neon) printf("NEON ");
    if (cpu->has_sse42) printf("SSE4.2 ");
    if (cpu->has_avx2) printf("AVX2 ");
    if (cpu->has_avx512f) printf("AVX-512 ");
    printf("\n");

    /* Parse arguments */
    int read_only = 0;
    const char* data_dir = NULL;
    int iterations = DEFAULT_ITERATIONS;
    int64_t num_rows = DEFAULT_ROWS;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--read-only") == 0 && i + 1 < argc) {
            read_only = 1;
            data_dir = argv[++i];
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (argv[i][0] != '-') {
            if (!read_only) {
                if (iterations == DEFAULT_ITERATIONS) {
                    iterations = atoi(argv[i]);
                } else {
                    num_rows = atoll(argv[i]);
                }
            }
        }
    }

    int rc;
    if (read_only) {
        if (!data_dir) {
            fprintf(stderr, "Error: --read-only requires a data directory\n");
            print_usage(argv[0]);
            return 1;
        }
        rc = run_read_only_benchmark(data_dir, iterations);
    } else {
        rc = run_full_benchmark(iterations, num_rows);
    }

    carquet_cleanup();
    return rc;
}
