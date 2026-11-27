/**
 * @file generate_test_files.c
 * @brief Generate test Parquet files for profiling (replaces generate_test_files.py)
 *
 * Usage:
 *   generate_test_files [output_dir] [num_rows]
 *
 * Default: ./profile_data/ with 10M rows
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>   /* _mkdir */
#endif

#include <carquet/carquet.h>

#define DEFAULT_NUM_ROWS 10000000
#define ROW_GROUP_SIZE   1000000

/* Simple LCG PRNG for reproducible results */
static uint32_t rng_state = 42;
static void rng_seed(uint32_t seed) { rng_state = seed; }
static uint32_t rng_next(void) {
    rng_state = rng_state * 1103515245 + 12345;
    return (rng_state >> 16) & 0x7FFFFFFF;
}
static double rng_double(void) {
    return (double)rng_next() / 2147483647.0;
}

static void ensure_dir(const char* path) {
    struct stat st;
    if (stat(path, &st) != 0) {
        #ifdef _WIN32
        _mkdir(path);
        #else
        mkdir(path, 0755);
        #endif
    }
}

static long get_file_size(const char* path) {
    struct stat st;
    if (stat(path, &st) == 0) return (long)st.st_size;
    return 0;
}

static int write_test_file(const char* filename, int num_rows,
                            carquet_compression_t compression,
                            const char* comp_name) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    printf("Writing %s...\n", filename);

    /* Create schema matching generate_test_files.py */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        fprintf(stderr, "  Failed to create schema: %s\n", err.message);
        return 1;
    }

    (void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "value_i32", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "value_f64", CARQUET_PHYSICAL_DOUBLE,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "value_f32", CARQUET_PHYSICAL_FLOAT,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "category", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = compression;
    opts.row_group_size = ROW_GROUP_SIZE;

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "  Failed to create writer: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }

    /* Allocate data buffers - write in row-group-sized chunks */
    int chunk = ROW_GROUP_SIZE;
    int64_t* ids = malloc(chunk * sizeof(int64_t));
    int32_t* i32_vals = malloc(chunk * sizeof(int32_t));
    double* f64_vals = malloc(chunk * sizeof(double));
    float* f32_vals = malloc(chunk * sizeof(float));
    int32_t* categories = malloc(chunk * sizeof(int32_t));

    if (!ids || !i32_vals || !f64_vals || !f32_vals || !categories) {
        fprintf(stderr, "  Out of memory\n");
        carquet_writer_close(writer);
        carquet_schema_free(schema);
        return 1;
    }

    rng_seed(42);

    int rows_written = 0;
    while (rows_written < num_rows) {
        int batch = num_rows - rows_written;
        if (batch > chunk) batch = chunk;

        for (int i = 0; i < batch; i++) {
            int64_t row = rows_written + i;
            ids[i] = row * 1000 + (int64_t)(rng_next() % 100);
            i32_vals[i] = (int32_t)(rng_next() % 1000000);
            f64_vals[i] = (double)row * 0.001 + rng_double() * 0.01;
            f32_vals[i] = (float)(rng_double() * 100.0);
            categories[i] = (int32_t)(rng_next() % 100);
        }

        (void)carquet_writer_write_batch(writer, 0, ids, batch, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 1, i32_vals, batch, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 2, f64_vals, batch, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 3, f32_vals, batch, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 4, categories, batch, NULL, NULL);

        rows_written += batch;
    }

    (void)carquet_writer_close(writer);

    free(ids);
    free(i32_vals);
    free(f64_vals);
    free(f32_vals);
    free(categories);
    carquet_schema_free(schema);

    long fsize = get_file_size(filename);
    printf("  Size: %.1f MB (%s compression)\n",
           fsize / (1024.0 * 1024.0), comp_name);

    /* Verify readability */
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    carquet_reader_t* reader = carquet_reader_open(filename, &ropts, &err);
    if (!reader) {
        fprintf(stderr, "  VERIFY FAILED: %s\n", err.message);
        return 1;
    }
    int64_t total = carquet_reader_num_rows(reader);
    int32_t ncols = carquet_reader_num_columns(reader);
    carquet_reader_close(reader);

    if (total != num_rows) {
        fprintf(stderr, "  VERIFY FAILED: row count %lld vs %d\n", (long long)total, num_rows);
        return 1;
    }
    printf("  Verified: %lld rows, %d columns\n", (long long)total, ncols);

    return 0;
}

int main(int argc, char* argv[]) {
    const char* output_dir = argc > 1 ? argv[1] : "benchmark/profile_data";
    int num_rows = argc > 2 ? atoi(argv[2]) : DEFAULT_NUM_ROWS;

    if (num_rows <= 0) {
        fprintf(stderr, "Invalid num_rows: %d\n", num_rows);
        return 1;
    }

    ensure_dir(output_dir);

    printf("Generating %d rows of test data...\n\n", num_rows);

    struct { carquet_compression_t codec; const char* name; const char* suffix; } comps[] = {
        { CARQUET_COMPRESSION_UNCOMPRESSED, "NONE",   "none"   },
        { CARQUET_COMPRESSION_SNAPPY,       "SNAPPY", "snappy" },
        { CARQUET_COMPRESSION_ZSTD,         "ZSTD",   "zstd"   },
    };

    for (int i = 0; i < 3; i++) {
        char filename[512];
        snprintf(filename, sizeof(filename), "%s/test_%s.parquet", output_dir, comps[i].suffix);

        if (write_test_file(filename, num_rows, comps[i].codec, comps[i].name) != 0) {
            fprintf(stderr, "Failed to write %s\n", filename);
            return 1;
        }
        printf("\n");
    }

    /* Write file list */
    char list_path[512];
    snprintf(list_path, sizeof(list_path), "%s/files.txt", output_dir);
    FILE* f = fopen(list_path, "w");
    if (f) {
        for (int i = 0; i < 3; i++) {
            fprintf(f, "test_%s.parquet\n", comps[i].suffix);
        }
        fclose(f);
    }

    printf("Done! Files written to %s/\n", output_dir);
    return 0;
}
