/**
 * @file append_rows.c
 * @brief Append row groups to an existing Parquet file
 *
 * Demonstrates carquet_writer_open_append(): write an initial file, then
 * reopen it and add more row groups without a read-then-rewrite. The schema
 * passed to open_append must match the file's existing leaf columns.
 */

#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(expr) do { if ((expr) != CARQUET_OK) { \
    fprintf(stderr, "FAIL at %s:%d\n", __FILE__, __LINE__); return 1; } } while(0)

static carquet_schema_t* make_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_column(s, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    return s;
}

/* Write rows [base, base+n) as a fresh file. */
static int write_initial(const char* path, int64_t base, int64_t n) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { fprintf(stderr, "create: %s\n", err.message); carquet_schema_free(schema); return 1; }

    int64_t* ids  = malloc((size_t)n * sizeof(int64_t));
    double*  vals = malloc((size_t)n * sizeof(double));
    for (int64_t i = 0; i < n; i++) { ids[i] = base + i; vals[i] = (double)(base + i) * 0.5; }

    CHECK(carquet_writer_write_batch(w, 0, ids,  n, NULL, NULL));
    CHECK(carquet_writer_write_batch(w, 1, vals, n, NULL, NULL));
    CHECK(carquet_writer_close(w));

    free(ids); free(vals);
    carquet_schema_free(schema);
    return 0;
}

/* Append rows [base, base+n) as a new row group in an existing file. */
static int append_rows(const char* path, int64_t base, int64_t n) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    carquet_writer_t* w = carquet_writer_open_append(path, schema, &opts, &err);
    if (!w) { fprintf(stderr, "open_append: %s\n", err.message); carquet_schema_free(schema); return 1; }

    int64_t* ids  = malloc((size_t)n * sizeof(int64_t));
    double*  vals = malloc((size_t)n * sizeof(double));
    for (int64_t i = 0; i < n; i++) { ids[i] = base + i; vals[i] = (double)(base + i) * 0.5; }

    CHECK(carquet_writer_write_batch(w, 0, ids,  n, NULL, NULL));
    CHECK(carquet_writer_write_batch(w, 1, vals, n, NULL, NULL));
    CHECK(carquet_writer_close(w));

    free(ids); free(vals);
    carquet_schema_free(schema);
    return 0;
}

/* Read the whole file back and verify the contiguous id sequence. */
static int verify(const char* path, int64_t expected_rows) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { fprintf(stderr, "open: %s\n", err.message); return 1; }

    int64_t nrows = carquet_reader_num_rows(r);
    printf("  file now has %lld rows, %d row group(s)\n",
           (long long)nrows, carquet_reader_num_row_groups(r));

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);

    int64_t seen = 0;
    int64_t expect_id = 0;
    int rc = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nulls; int64_t n;
        carquet_row_batch_column(batch, 0, &data, &nulls, &n);
        const int64_t* ids = (const int64_t*)data;
        for (int64_t i = 0; i < n; i++) {
            if (ids[i] != expect_id) {
                fprintf(stderr, "  MISMATCH at row %lld: got %lld want %lld\n",
                        (long long)(seen + i), (long long)ids[i], (long long)expect_id);
                rc = 1;
            }
            expect_id++;
        }
        seen += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    if (seen != expected_rows) {
        fprintf(stderr, "  ROW COUNT MISMATCH: got %lld want %lld\n",
                (long long)seen, (long long)expected_rows);
        rc = 1;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    if (rc == 0) printf("  verified ids 0..%lld contiguous\n", (long long)(expected_rows - 1));
    return rc;
}

int main(void) {
    carquet_init();
    const char* path = "/tmp/carquet_append_example.parquet";

    printf("Writing initial file (1000 rows)...\n");
    if (write_initial(path, 0, 1000)) return 1;
    if (verify(path, 1000)) return 1;

    printf("\nAppending 500 rows...\n");
    if (append_rows(path, 1000, 500)) return 1;
    if (verify(path, 1500)) return 1;

    printf("\nAppending another 250 rows...\n");
    if (append_rows(path, 1500, 250)) return 1;
    if (verify(path, 1750)) return 1;

    remove(path);
    printf("\nDone.\n");
    return 0;
}
