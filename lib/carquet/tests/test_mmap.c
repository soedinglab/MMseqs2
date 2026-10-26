/**
 * @file test_mmap.c
 * @brief Tests for memory-mapped I/O and zero-copy reading
 */

#include <carquet/carquet.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* Batch column buffers may be zero-copy views into the mmap'd file and are
 * not guaranteed to satisfy the alignment of the element type. Read elements
 * via memcpy to avoid unaligned-load undefined behavior. */
static inline int64_t rd_i64(const void* base, int64_t i) {
    int64_t v; memcpy(&v, (const uint8_t*)base + i * (int64_t)sizeof v, sizeof v); return v;
}
static inline int32_t rd_i32(const void* base, int64_t i) {
    int32_t v; memcpy(&v, (const uint8_t*)base + i * (int64_t)sizeof v, sizeof v); return v;
}
static inline double rd_f64(const void* base, int64_t i) {
    double v; memcpy(&v, (const uint8_t*)base + i * (int64_t)sizeof v, sizeof v); return v;
}

/* ============================================================================
 * Helper: Create a test file with uncompressed data
 * ============================================================================
 */

static const char* TEST_FILE = "test_mmap_data.parquet";

static int create_test_file_with_page_size(int64_t num_rows, int64_t page_size) {
    carquet_error_t error = CARQUET_ERROR_INIT;

    /* Create schema with REQUIRED INT64 column (eligible for zero-copy) */
    carquet_schema_t* schema = carquet_schema_create(&error);
    if (!schema) return 1;

    /* Add required int64 column */
    carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64,
                               NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Add required double column */
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE,
                               NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Writer options - NO COMPRESSION for zero-copy */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.row_group_size = num_rows * 32;  /* Keep test data in a single row group */
    if (page_size > 0) {
        opts.page_size = page_size;
    }

    carquet_writer_t* writer = carquet_writer_create(TEST_FILE, schema, &opts, &error);
    if (!writer) {
        carquet_schema_free(schema);
        return 1;
    }

    /* Write data */
    int64_t* ids = malloc(sizeof(int64_t) * num_rows);
    double* values = malloc(sizeof(double) * num_rows);

    for (int64_t i = 0; i < num_rows; i++) {
        ids[i] = i * 100;
        values[i] = (double)i * 3.14159;
    }

    carquet_writer_write_batch(writer, 0, ids, num_rows, NULL, NULL);
    carquet_writer_write_batch(writer, 1, values, num_rows, NULL, NULL);

    carquet_writer_close(writer);
    carquet_schema_free(schema);
    free(ids);
    free(values);

    return 0;
}

static int create_test_file(int64_t num_rows) {
    return create_test_file_with_page_size(num_rows, 0);
}

static int create_compressed_optional_file(void) {
    carquet_error_t error = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&error);
    if (!schema) return 1;

    carquet_schema_add_column(schema, "maybe_id", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.row_group_size = 10;

    carquet_writer_t* writer = carquet_writer_create(TEST_FILE, schema, &opts, &error);
    if (!writer) {
        carquet_schema_free(schema);
        return 1;
    }

    for (int rg = 0; rg < 2; rg++) {
        int32_t values[5];
        int16_t def[10];
        int nn = 0;

        for (int i = 0; i < 10; i++) {
            int global = rg * 10 + i;
            if ((global % 2) == 0) {
                def[i] = 1;
                values[nn++] = global * 10;
            } else {
                def[i] = 0;
            }
        }

        if (carquet_writer_write_batch(writer, 0, values, 10, def, NULL) != CARQUET_OK) {
            carquet_writer_abort(writer);
            carquet_schema_free(schema);
            return 1;
        }

        if (rg == 0 && carquet_writer_new_row_group(writer) != CARQUET_OK) {
            carquet_writer_abort(writer);
            carquet_schema_free(schema);
            return 1;
        }
    }

    carquet_status_t status = carquet_writer_close(writer);
    carquet_schema_free(schema);
    return status == CARQUET_OK ? 0 : 1;
}

/* ============================================================================
 * Test: mmap reader opens correctly
 * ============================================================================
 */

static int test_mmap_open(void) {
    const char* name = "mmap_open";
    carquet_error_t error = CARQUET_ERROR_INIT;

    /* Create test file */
    if (create_test_file(1000) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Open with mmap */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader with mmap");
    }

    /* Verify mmap is active */
    if (!carquet_reader_is_mmap(reader)) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "mmap should be active");
    }

    /* Verify metadata */
    if (carquet_reader_num_rows(reader) != 1000) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Wrong row count");
    }

    if (carquet_reader_num_columns(reader) != 2) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Wrong column count");
    }

    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: zero-copy eligibility check
 * ============================================================================
 */

static int test_zero_copy_eligibility(void) {
    const char* name = "zero_copy_eligibility";
    carquet_error_t error = CARQUET_ERROR_INIT;

    /* Create test file */
    if (create_test_file(1000) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Open with mmap */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    /* Column 0 (INT64, REQUIRED, uncompressed) should be zero-copy eligible */
    if (!carquet_reader_can_zero_copy(reader, 0, 0)) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "INT64 column should be zero-copy eligible");
    }

    /* Column 1 (DOUBLE, REQUIRED, uncompressed) should also be eligible */
    if (!carquet_reader_can_zero_copy(reader, 0, 1)) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "DOUBLE column should be zero-copy eligible");
    }

    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: Read data via mmap
 * ============================================================================
 */

static int test_mmap_read_data(void) {
    const char* name = "mmap_read_data";
    carquet_error_t error = CARQUET_ERROR_INIT;
    int64_t num_rows = 1000;

    /* Create test file */
    if (create_test_file(num_rows) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Open with mmap */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    /* Read column 0 (INT64) */
    carquet_column_reader_t* col_reader = carquet_reader_get_column(reader, 0, 0, &error);
    if (!col_reader) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to get column reader");
    }

    int64_t* data = malloc(sizeof(int64_t) * num_rows);
    int64_t values_read = carquet_column_read_batch(col_reader, data, num_rows, NULL, NULL);

    if (values_read != num_rows) {
        free(data);
        carquet_column_reader_free(col_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Wrong number of values read");
    }

    /* Verify data */
    for (int64_t i = 0; i < num_rows; i++) {
        if (data[i] != i * 100) {
            free(data);
            carquet_column_reader_free(col_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Data mismatch");
        }
    }

    free(data);
    carquet_column_reader_free(col_reader);
    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: Batch reader with mmap
 * ============================================================================
 */

static int test_mmap_batch_reader(void) {
    const char* name = "mmap_batch_reader";
    carquet_error_t error = CARQUET_ERROR_INIT;
    int64_t num_rows = 1000;

    /* Create test file */
    if (create_test_file(num_rows) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Open with mmap */
    carquet_reader_options_t reader_opts;
    carquet_reader_options_init(&reader_opts);
    reader_opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &reader_opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    /* Create batch reader */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = num_rows;  /* Read all in one batch */

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &error);
    if (!batch_reader) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to create batch reader");
    }

    /* Read batch */
    carquet_row_batch_t* batch = NULL;
    carquet_status_t status = carquet_batch_reader_next(batch_reader, &batch);
    if (status != CARQUET_OK || !batch) {
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to read batch");
    }

    /* Verify row count */
    if (carquet_row_batch_num_rows(batch) != num_rows) {
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Wrong batch row count");
    }

    /* Get column data */
    const void* data;
    const uint8_t* null_bitmap;
    int64_t col_num_values;

    status = carquet_row_batch_column(batch, 0, &data, &null_bitmap, &col_num_values);
    if (status != CARQUET_OK) {
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to get column data");
    }

    /* Verify data (batch buffer may be an unaligned mmap view). */
    for (int64_t i = 0; i < num_rows; i++) {
        if (rd_i64(data, i) != i * 100) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Data mismatch in batch");
        }
    }

    carquet_row_batch_free(batch);
    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: Oversized mmap batch request stays page-aligned for zero-copy
 * ============================================================================
 */

static int test_mmap_batch_reader_page_aligned(void) {
    const char* name = "mmap_batch_reader_page_aligned";
    carquet_error_t error = CARQUET_ERROR_INIT;
    const int64_t page_size = 4096;
    const int64_t num_rows = 4096;

    if (create_test_file_with_page_size(num_rows, page_size) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    carquet_reader_options_t reader_opts;
    carquet_reader_options_init(&reader_opts);
    reader_opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &reader_opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = (int32_t)num_rows;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &error);
    if (!batch_reader) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to create batch reader");
    }

    int64_t total_rows = 0;
    int64_t expected_row = 0;
    int batch_count = 0;
    bool saw_split_batch = false;
    carquet_row_batch_t* batch = NULL;

    while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        const void* data = NULL;
        const uint8_t* null_bitmap = NULL;
        int64_t col_num_values = 0;
        const int64_t* ids;
        const double* values;
        int64_t batch_rows = carquet_row_batch_num_rows(batch);

        if (batch_count == 0) {
            if (batch_rows >= num_rows) {
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                TEST_FAIL(name, "Oversized batch request should be split across page boundaries");
            }
            saw_split_batch = true;
        }

        if (batch_rows <= 0) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Batch should contain rows");
        }

        if (carquet_row_batch_column(batch, 0, &data, &null_bitmap, &col_num_values) != CARQUET_OK ||
            !data || null_bitmap != NULL || col_num_values != batch_rows) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Failed to read INT64 batch column");
        }
        ids = (const int64_t*)data;  /* compared via rd_i64 (may be unaligned) */

        if (carquet_row_batch_column(batch, 1, &data, &null_bitmap, &col_num_values) != CARQUET_OK ||
            !data || null_bitmap != NULL || col_num_values != batch_rows) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Failed to read DOUBLE batch column");
        }
        values = (const double*)data;

        for (int64_t i = 0; i < batch_rows; i++) {
            if (rd_i64(ids, i) != expected_row * 100 ||
                fabs(rd_f64(values, i) - ((double)expected_row * 3.14159)) > 1e-9) {
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                TEST_FAIL(name, "Batch data mismatch");
            }
            expected_row++;
        }

        total_rows += batch_rows;
        batch_count++;
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    if (total_rows != num_rows || batch_count < 2 || expected_row != num_rows || !saw_split_batch) {
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Unexpected batch segmentation");
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: compressed mmap batch reader keeps optional nulls intact
 * ============================================================================
 */

static int test_mmap_compressed_optional_batch(void) {
    const char* name = "mmap_compressed_optional_batch";
    carquet_error_t error = CARQUET_ERROR_INIT;

    if (create_compressed_optional_file() != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    carquet_reader_options_t reader_opts;
    carquet_reader_options_init(&reader_opts);
    reader_opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &reader_opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 20;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &error);
    if (!batch_reader) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "Failed to create batch reader");
    }

    int64_t total = 0;
    carquet_status_t status;
    carquet_row_batch_t* batch = NULL;
    while ((status = carquet_batch_reader_next(batch_reader, &batch)) == CARQUET_OK && batch) {
        const void* data = NULL;
        const uint8_t* null_bitmap = NULL;
        int64_t n = 0;

        if (carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n) != CARQUET_OK ||
            !data || !null_bitmap || n <= 0) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            TEST_FAIL(name, "Failed to read optional batch column");
        }

        const int32_t* values = (const int32_t*)data;  /* via rd_i32 (may be unaligned) */
        for (int64_t i = 0; i < n; i++) {
            int64_t row = total + i;
            int is_present = (null_bitmap[i / 8] >> (i % 8)) & 1;
            if ((row % 2) == 0) {
                if (!is_present || rd_i32(values, i) != row * 10) {
                    carquet_row_batch_free(batch);
                    carquet_batch_reader_free(batch_reader);
                    carquet_reader_close(reader);
                    TEST_FAIL(name, "present optional value mismatch");
                }
            } else if (is_present || rd_i32(values, i) != 0) {
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                TEST_FAIL(name, "null optional value mismatch");
            }
        }

        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    if (status != CARQUET_ERROR_END_OF_DATA || total != 20) {
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        TEST_FAIL(name, "Unexpected optional batch result");
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: Compare mmap vs fread results
 * ============================================================================
 */

static int test_mmap_vs_fread(void) {
    const char* name = "mmap_vs_fread";
    carquet_error_t error = CARQUET_ERROR_INIT;
    int64_t num_rows = 5000;

    /* Create test file */
    if (create_test_file(num_rows) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Read with mmap */
    carquet_reader_options_t mmap_opts;
    carquet_reader_options_init(&mmap_opts);
    mmap_opts.use_mmap = true;

    carquet_reader_t* mmap_reader = carquet_reader_open(TEST_FILE, &mmap_opts, &error);
    if (!mmap_reader) {
        TEST_FAIL(name, "Failed to open mmap reader");
    }

    carquet_column_reader_t* mmap_col = carquet_reader_get_column(mmap_reader, 0, 0, &error);
    int64_t* mmap_data = malloc(sizeof(int64_t) * num_rows);
    carquet_column_read_batch(mmap_col, mmap_data, num_rows, NULL, NULL);

    carquet_column_reader_free(mmap_col);
    carquet_reader_close(mmap_reader);

    /* Read with fread */
    carquet_reader_options_t fread_opts;
    carquet_reader_options_init(&fread_opts);
    fread_opts.use_mmap = false;

    carquet_reader_t* fread_reader = carquet_reader_open(TEST_FILE, &fread_opts, &error);
    if (!fread_reader) {
        free(mmap_data);
        TEST_FAIL(name, "Failed to open fread reader");
    }

    carquet_column_reader_t* fread_col = carquet_reader_get_column(fread_reader, 0, 0, &error);
    int64_t* fread_data = malloc(sizeof(int64_t) * num_rows);
    carquet_column_read_batch(fread_col, fread_data, num_rows, NULL, NULL);

    carquet_column_reader_free(fread_col);
    carquet_reader_close(fread_reader);

    /* Compare results */
    int mismatch = 0;
    for (int64_t i = 0; i < num_rows; i++) {
        if (mmap_data[i] != fread_data[i]) {
            mismatch = 1;
            break;
        }
    }

    free(mmap_data);
    free(fread_data);

    if (mismatch) {
        TEST_FAIL(name, "mmap and fread results differ");
    }

    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Test: Fallback to fread when mmap fails
 * ============================================================================
 */

static int test_fread_fallback(void) {
    const char* name = "fread_fallback";
    carquet_error_t error = CARQUET_ERROR_INIT;

    /* Create test file */
    if (create_test_file(100) != 0) {
        TEST_FAIL(name, "Failed to create test file");
    }

    /* Open without mmap */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = false;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &opts, &error);
    if (!reader) {
        TEST_FAIL(name, "Failed to open reader");
    }

    /* Verify mmap is NOT active */
    if (carquet_reader_is_mmap(reader)) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "mmap should NOT be active");
    }

    /* Zero-copy should not be possible without mmap */
    if (carquet_reader_can_zero_copy(reader, 0, 0)) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "zero-copy should not be possible without mmap");
    }

    carquet_reader_close(reader);
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

/* Write `num_rgs` row groups of `rows_per_rg` REQUIRED INT64 values under ZSTD
 * (so the multi-row-group mmap pipeline engages). constant != 0 writes that
 * value everywhere, otherwise a pseudo-random sequence. */
static int write_pipeline_file(const char* path, int num_rgs, int64_t rows_per_rg,
                               int64_t constant) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;
    if (carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) return 1;
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { carquet_schema_free(schema); return 1; }
    if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN) != CARQUET_OK) return 1;
    int64_t* v = malloc((size_t)rows_per_rg * sizeof(int64_t));
    if (!v) return 1;
    uint64_t x = 0x243F6A8885A308D3ULL;
    int rc = 0;
    for (int g = 0; g < num_rgs && rc == 0; g++) {
        for (int64_t i = 0; i < rows_per_rg; i++) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            v[i] = constant ? constant : (int64_t)x;
        }
        if (carquet_writer_write_batch(w, 0, v, rows_per_rg, NULL, NULL) != CARQUET_OK) rc = 1;
        if (rc == 0 && g + 1 < num_rgs &&
            carquet_writer_new_row_group(w) != CARQUET_OK) rc = 1;
    }
    free(v);
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    carquet_schema_free(schema);
    return rc;
}

/* Read every batch through an mmap batch reader. Returns the last status and
 * the rows delivered; *sum accumulates the values. */
static carquet_status_t drain_mmap_batches(const char* path, int64_t* rows,
                                           uint64_t* sum, bool* had_error_detail) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = true;
    *rows = 0; *sum = 0; *had_error_detail = false;
    carquet_reader_t* r = carquet_reader_open(path, &ropts, &err);
    if (!r) return CARQUET_ERROR_FILE_OPEN;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); return err.code; }
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st;
    while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nulls; int64_t n;
        if (carquet_row_batch_column(batch, 0, &data, &nulls, &n) != CARQUET_OK) break;
        for (int64_t i = 0; i < n; i++) *sum += (uint64_t)rd_i64(data, i);
        *rows += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    *had_error_detail = carquet_batch_reader_last_error(br) != NULL;
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return st;
}

/* Constant data compresses to a few hundred bytes per row group, so the file
 * holds far more than 8 rows per byte. That is legitimate, and every row must
 * still come back: the pre-read pipeline's allocation guard used to treat it
 * as malformed and end the read with zero rows and END_OF_DATA. */
static int test_mmap_highly_compressible(void) {
    const char* path = "test_mmap_compressible.parquet";
    enum { RGS = 4, ROWS = 200000 };
    if (write_pipeline_file(path, RGS, ROWS, 7) != 0)
        TEST_FAIL("mmap_highly_compressible", "write failed");
    int64_t rows; uint64_t sum; bool detail;
    carquet_status_t st = drain_mmap_batches(path, &rows, &sum, &detail);
    remove(path);
    if (st != CARQUET_ERROR_END_OF_DATA && st != CARQUET_OK)
        TEST_FAIL("mmap_highly_compressible", "read failed");
    if (rows != (int64_t)RGS * ROWS || sum != (uint64_t)RGS * ROWS * 7)
        TEST_FAIL("mmap_highly_compressible", "rows were dropped");
    TEST_PASS("mmap_highly_compressible");
    return 0;
}

/* Damage a page payload in the middle of a multi-row-group file. The pipeline
 * decodes row groups ahead on worker threads; a column that failed there must
 * surface as an error from next(), with detail, and never as a normal end of
 * data or as a full-size batch over a buffer the failed task left unwritten. */
static int test_mmap_pipeline_corruption(void) {
    const char* path = "test_mmap_corrupt.parquet";
    enum { RGS = 4, ROWS = 50000 };
    if (write_pipeline_file(path, RGS, ROWS, 0) != 0)
        TEST_FAIL("mmap_pipeline_corruption", "write failed");

    int64_t rows; uint64_t good_sum; bool detail;
    carquet_status_t st = drain_mmap_batches(path, &rows, &good_sum, &detail);
    if ((st != CARQUET_ERROR_END_OF_DATA && st != CARQUET_OK) || rows != (int64_t)RGS * ROWS) {
        remove(path);
        TEST_FAIL("mmap_pipeline_corruption", "clean read failed");
    }

    FILE* f = fopen(path, "r+b");
    if (!f) { remove(path); TEST_FAIL("mmap_pipeline_corruption", "reopen failed"); }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    /* Random INT64 barely compresses, so the file is ~RGS * ROWS * 8 bytes of
     * page payload; 5/8 of the way in lands inside the third row group. */
    long at = size / 8 * 5;
    unsigned char junk[64];
    fseek(f, at, SEEK_SET);
    if (fread(junk, 1, sizeof(junk), f) != sizeof(junk)) { fclose(f); remove(path);
        TEST_FAIL("mmap_pipeline_corruption", "read for corruption failed"); }
    for (size_t i = 0; i < sizeof(junk); i++) junk[i] = (unsigned char)~junk[i];
    fseek(f, at, SEEK_SET);
    fwrite(junk, 1, sizeof(junk), f);
    fclose(f);

    uint64_t sum;
    st = drain_mmap_batches(path, &rows, &sum, &detail);
    remove(path);
    if (st == CARQUET_OK || st == CARQUET_ERROR_END_OF_DATA)
        TEST_FAIL("mmap_pipeline_corruption", "corruption was not reported");
    if (!detail)
        TEST_FAIL("mmap_pipeline_corruption", "no error detail");
    if (rows >= (int64_t)RGS * ROWS)
        TEST_FAIL("mmap_pipeline_corruption", "corrupt row group was delivered");
    TEST_PASS("mmap_pipeline_corruption");
    return 0;
}

/* ============================================================================
 * Page-range segment reads (intra-chunk parallel decode)
 * ============================================================================
 * A file with one large row group gives the pipeline only one task per
 * column, so it cuts each coalescable chunk into page-aligned segments that
 * decode concurrently. These tests pin down that path: every value must come
 * back exactly (including the BYTE_STREAM_SPLIT double column, whose segment
 * tasks stage through the shared scratch pool), and damage inside one
 * segment must surface as an error with detail, never as silently wrong or
 * missing rows.
 */

/* One row group of `rows` REQUIRED (INT64 id = 3*i, DOUBLE value = 0.5*i)
 * under Snappy. */
static int write_segment_file(const char* path, int64_t rows) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;
    if (carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) return 1;
    if (carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) return 1;
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;
    opts.row_group_size = rows * 64;   /* everything in one row group */
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { carquet_schema_free(schema); return 1; }
    int64_t* ids = malloc((size_t)rows * sizeof(int64_t));
    double* vals = malloc((size_t)rows * sizeof(double));
    if (!ids || !vals) { free(ids); free(vals); return 1; }
    for (int64_t i = 0; i < rows; i++) { ids[i] = 3 * i; vals[i] = 0.5 * (double)i; }
    int rc = 0;
    if (carquet_writer_write_batch(w, 0, ids, rows, NULL, NULL) != CARQUET_OK) rc = 1;
    if (carquet_writer_write_batch(w, 1, vals, rows, NULL, NULL) != CARQUET_OK) rc = 1;
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    carquet_schema_free(schema);
    free(ids); free(vals);
    return rc;
}

/* Read the segment file back through the mmap batch reader, checking every
 * value against its formula. *rows counts rows delivered; *bad counts
 * mismatches. */
static carquet_status_t verify_segment_file(const char* path, int64_t* rows,
                                            int64_t* bad, bool* detail) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = true;
    *rows = 0; *bad = 0; *detail = false;
    carquet_reader_t* r = carquet_reader_open(path, &ropts, &err);
    if (!r) return CARQUET_ERROR_FILE_OPEN;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 70000;   /* not a divisor of the page size: batches straddle segments */
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); return err.code; }
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st;
    while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
        const void* ids; const void* vals; const uint8_t* nulls; int64_t n, n2;
        if (carquet_row_batch_column(batch, 0, &ids, &nulls, &n) != CARQUET_OK ||
            carquet_row_batch_column(batch, 1, &vals, &nulls, &n2) != CARQUET_OK ||
            n != n2) { *bad += 1; break; }
        for (int64_t i = 0; i < n; i++) {
            int64_t row = *rows + i;
            if (rd_i64(ids, i) != 3 * row) (*bad)++;
            if (rd_f64(vals, i) != 0.5 * (double)row) (*bad)++;
        }
        *rows += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    *detail = carquet_batch_reader_last_error(br) != NULL;
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return st;
}

static int test_mmap_segment_read_exact(void) {
    const char* path = "test_mmap_segments.parquet";
    const int64_t ROWS = 600000;   /* ~5 default-size pages per column */
    if (write_segment_file(path, ROWS) != 0)
        TEST_FAIL("mmap_segment_read_exact", "write failed");
    int64_t rows, bad; bool detail;
    carquet_status_t st = verify_segment_file(path, &rows, &bad, &detail);
    remove(path);
    if (st != CARQUET_ERROR_END_OF_DATA && st != CARQUET_OK)
        TEST_FAIL("mmap_segment_read_exact", "read failed");
    if (rows != ROWS)
        TEST_FAIL("mmap_segment_read_exact", "row count mismatch");
    if (bad != 0)
        TEST_FAIL("mmap_segment_read_exact", "value mismatch");
    TEST_PASS("mmap_segment_read_exact");
    return 0;
}

/* Same as test_mmap_pipeline_corruption but with a single large row group,
 * so the damaged page sits inside one of several concurrently decoded
 * segments of the same chunk. The segment cannot recover on its own (its
 * siblings share the column reader); the consumer must detect the short
 * column, re-read it through the checked page reader, and report. */
static int test_mmap_segment_corruption(void) {
    const char* path = "test_mmap_segment_corrupt.parquet";
    const int64_t ROWS = 600000;
    if (write_pipeline_file(path, 1, ROWS, 0) != 0)
        TEST_FAIL("mmap_segment_corruption", "write failed");

    int64_t rows; uint64_t good_sum; bool detail;
    carquet_status_t st = drain_mmap_batches(path, &rows, &good_sum, &detail);
    if ((st != CARQUET_ERROR_END_OF_DATA && st != CARQUET_OK) || rows != ROWS) {
        remove(path);
        TEST_FAIL("mmap_segment_corruption", "clean read failed");
    }

    FILE* f = fopen(path, "r+b");
    if (!f) { remove(path); TEST_FAIL("mmap_segment_corruption", "reopen failed"); }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    long at = size / 2;   /* inside the chunk: the third of ~five pages */
    unsigned char junk[64];
    fseek(f, at, SEEK_SET);
    if (fread(junk, 1, sizeof(junk), f) != sizeof(junk)) { fclose(f); remove(path);
        TEST_FAIL("mmap_segment_corruption", "read for corruption failed"); }
    for (size_t i = 0; i < sizeof(junk); i++) junk[i] = (unsigned char)~junk[i];
    fseek(f, at, SEEK_SET);
    fwrite(junk, 1, sizeof(junk), f);
    fclose(f);

    uint64_t sum;
    st = drain_mmap_batches(path, &rows, &sum, &detail);
    remove(path);
    if (st == CARQUET_OK || st == CARQUET_ERROR_END_OF_DATA)
        TEST_FAIL("mmap_segment_corruption", "corruption was not reported");
    if (!detail)
        TEST_FAIL("mmap_segment_corruption", "no error detail");
    if (rows != 0)
        TEST_FAIL("mmap_segment_corruption", "rows of the corrupt row group were delivered");
    TEST_PASS("mmap_segment_corruption");
    return 0;
}

int main(void) {
    printf("=== Memory-Mapped I/O Tests ===\n\n");

    int failures = 0;

    failures += test_mmap_open();
    failures += test_zero_copy_eligibility();
    failures += test_mmap_read_data();
    failures += test_mmap_batch_reader();
    failures += test_mmap_batch_reader_page_aligned();
    failures += test_mmap_compressed_optional_batch();
    failures += test_mmap_vs_fread();
    failures += test_fread_fallback();
    failures += test_mmap_highly_compressible();
    failures += test_mmap_pipeline_corruption();
    failures += test_mmap_segment_read_exact();
    failures += test_mmap_segment_corruption();

    /* Cleanup */
    remove(TEST_FILE);

    printf("\n=== Results: %d failures ===\n", failures);
    return failures > 0 ? 1 : 0;
}
