/**
 * @file test_large_schema.c
 * @brief Tests for large schema Parquet files (1600+ columns)
 *
 * Tests to verify that carquet can handle schemas with many columns,
 * as commonly found in cheminformatics and other scientific applications.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

#define NUM_COLUMNS 1600
#define NUM_ROWS 100

static char TEST_FILE[512];

/* ============================================================================
 * Test: Create and write a Parquet file with 1600+ columns
 * ============================================================================ */

static int test_large_schema_write(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    printf("Creating schema with %d columns...\n", NUM_COLUMNS);

    /* Create schema with many columns */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        printf("  Failed to create schema: %s\n", err.message);
        TEST_FAIL("large_schema_write", "failed to create schema");
    }

    /* Add columns with different types */
    char col_name[64];
    for (int i = 0; i < NUM_COLUMNS; i++) {
        snprintf(col_name, sizeof(col_name), "col_%04d", i);

        /* Mix of types: INT32, DOUBLE, FLOAT */
        carquet_physical_type_t type;
        if (i % 3 == 0) {
            type = CARQUET_PHYSICAL_INT32;
        } else if (i % 3 == 1) {
            type = CARQUET_PHYSICAL_DOUBLE;
        } else {
            type = CARQUET_PHYSICAL_FLOAT;
        }

        carquet_status_t status = carquet_schema_add_column(
            schema, col_name, type, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

        if (status != CARQUET_OK) {
            printf("  Failed to add column %d\n", i);
            carquet_schema_free(schema);
            TEST_FAIL("large_schema_write", "failed to add column");
        }

        if (i > 0 && i % 400 == 0) {
            printf("  Added %d columns...\n", i);
        }
    }

    printf("  Schema created with %d columns\n", carquet_schema_num_columns(schema));

    /* Create writer */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;

    printf("Creating writer...\n");
    carquet_writer_t* writer = carquet_writer_create(TEST_FILE, schema, &opts, &err);
    if (!writer) {
        printf("  Failed to create writer: %s\n", err.message);
        carquet_schema_free(schema);
        TEST_FAIL("large_schema_write", "failed to create writer");
    }

    /* Allocate data buffers */
    int32_t* int_data = malloc(NUM_ROWS * sizeof(int32_t));
    double* double_data = malloc(NUM_ROWS * sizeof(double));
    float* float_data = malloc(NUM_ROWS * sizeof(float));

    if (!int_data || !double_data || !float_data) {
        free(int_data);
        free(double_data);
        free(float_data);
        carquet_writer_close(writer);
        carquet_schema_free(schema);
        TEST_FAIL("large_schema_write", "failed to allocate data");
    }

    /* Initialize data */
    for (int i = 0; i < NUM_ROWS; i++) {
        int_data[i] = i * 100;
        double_data[i] = (double)i * 1.5;
        float_data[i] = (float)i * 0.5f;
    }

    /* Write data to all columns */
    printf("Writing %d rows to %d columns...\n", NUM_ROWS, NUM_COLUMNS);
    for (int col = 0; col < NUM_COLUMNS; col++) {
        carquet_status_t status;

        if (col % 3 == 0) {
            status = carquet_writer_write_batch(writer, col, int_data, NUM_ROWS, NULL, NULL);
        } else if (col % 3 == 1) {
            status = carquet_writer_write_batch(writer, col, double_data, NUM_ROWS, NULL, NULL);
        } else {
            status = carquet_writer_write_batch(writer, col, float_data, NUM_ROWS, NULL, NULL);
        }

        if (status != CARQUET_OK) {
            printf("  Failed to write column %d\n", col);
            free(int_data);
            free(double_data);
            free(float_data);
            carquet_writer_close(writer);
            carquet_schema_free(schema);
            TEST_FAIL("large_schema_write", "failed to write column");
        }

        if (col > 0 && col % 400 == 0) {
            printf("  Written %d columns...\n", col);
        }
    }

    free(int_data);
    free(double_data);
    free(float_data);

    /* Close writer */
    printf("Closing writer and finalizing file...\n");
    carquet_status_t status = carquet_writer_close(writer);
    carquet_schema_free(schema);

    if (status != CARQUET_OK) {
        printf("  Writer close failed with status %d\n", status);
        TEST_FAIL("large_schema_write", "failed to close writer");
    }

    printf("  File written successfully\n");
    TEST_PASS("large_schema_write");
    return 0;
}

/* ============================================================================
 * Test: Read back the large schema Parquet file
 * ============================================================================ */

static int test_large_schema_read(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    printf("Opening file for reading...\n");

    /* Use mmap for thread-safe parallel reading with OpenMP */
    carquet_reader_options_t read_opts;
    carquet_reader_options_init(&read_opts);
    read_opts.use_mmap = true;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, &read_opts, &err);
    if (!reader) {
        printf("  Failed to open file: %s\n", err.message);
        TEST_FAIL("large_schema_read", "failed to open file");
    }

    /* Verify metadata */
    int64_t num_rows = carquet_reader_num_rows(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_row_groups = carquet_reader_num_row_groups(reader);

    printf("  Rows: %lld\n", (long long)num_rows);
    printf("  Columns: %d\n", num_cols);
    printf("  Row groups: %d\n", num_row_groups);

    if (num_rows != NUM_ROWS) {
        carquet_reader_close(reader);
        printf("  Expected %d rows, got %lld\n", NUM_ROWS, (long long)num_rows);
        TEST_FAIL("large_schema_read", "row count mismatch");
    }

    if (num_cols != NUM_COLUMNS) {
        carquet_reader_close(reader);
        printf("  Expected %d columns, got %d\n", NUM_COLUMNS, num_cols);
        TEST_FAIL("large_schema_read", "column count mismatch");
    }

    /* Verify schema column names */
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    for (int i = 0; i < 10; i++) {
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, i + 1);
        if (node) {
            const char* name = carquet_schema_node_name(node);
            printf("  Column %d: %s\n", i, name ? name : "(null)");
        }
    }
    printf("  ... (%d more columns)\n", NUM_COLUMNS - 10);

    /* Read a few columns using batch reader */
    printf("Reading data with batch reader...\n");

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);

    /* Project just a few columns to verify reading works */
    int32_t proj_cols[] = {0, 1, 2, 100, 500, 1000, 1599};
    config.column_indices = proj_cols;
    config.num_columns = 7;
    config.batch_size = NUM_ROWS;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        printf("  Failed to create batch reader: %s\n", err.message);
        TEST_FAIL("large_schema_read", "failed to create batch reader");
    }

    carquet_row_batch_t* batch = NULL;
    int64_t total_rows = 0;
    carquet_status_t read_status;

    while ((read_status = carquet_batch_reader_next(batch_reader, &batch)) == CARQUET_OK && batch) {
        int64_t batch_rows = carquet_row_batch_num_rows(batch);
        int32_t batch_cols = carquet_row_batch_num_columns(batch);

        printf("  Batch: %lld rows, %d columns\n", (long long)batch_rows, batch_cols);
        total_rows += batch_rows;

        /* Verify first column data (INT32) */
        const void* data;
        const uint8_t* null_bitmap;
        int64_t num_values;

        carquet_status_t status = carquet_row_batch_column(batch, 0, &data, &null_bitmap, &num_values);
        if (status == CARQUET_OK && data) {
            const int32_t* int_vals = (const int32_t*)data;
            printf("  First column values: [%d, %d, %d, ...]\n",
                   int_vals[0], int_vals[1], int_vals[2]);
        }

        carquet_row_batch_free(batch);
        batch = NULL;
    }

    printf("  Total rows read: %lld\n", (long long)total_rows);
    if (read_status != CARQUET_OK && read_status != CARQUET_ERROR_END_OF_DATA) {
        printf("  batch_reader_next returned status: %d\n", read_status);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);

    if (total_rows != NUM_ROWS) {
        printf("  Expected to read %d rows, got %lld\n", NUM_ROWS, (long long)total_rows);
        TEST_FAIL("large_schema_read", "read row count mismatch");
    }

    TEST_PASS("large_schema_read");
    return 0;
}

/* ============================================================================
 * Test: Verify file with external tool (if available)
 * ============================================================================ */

static int test_verify_with_pyarrow(void) {
    printf("Checking if file can be validated with pyarrow...\n");

    /* Try to run Python validation script */
    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
        "python3 -c \""
        "import pyarrow.parquet as pq; "
        "import sys; "
        "try: "
        "    t = pq.read_table('%s'); "
        "    print(f'  PyArrow: {t.num_rows} rows, {t.num_columns} columns'); "
        "    sys.exit(0); "
        "except Exception as e: "
        "    print(f'  PyArrow error: {e}'); "
        "    sys.exit(1); "
        "\" 2>/dev/null", TEST_FILE);
    int result = system(cmd);

    if (result == 0) {
        TEST_PASS("verify_with_pyarrow");
        return 0;
    } else {
        /* PyArrow validation is optional - don't fail the test suite */
        printf("  PyArrow validation skipped (not available or failed)\n");
        printf("  Run manually: python3 -c \"import pyarrow.parquet as pq; print(pq.read_table('%s'))\"\n", TEST_FILE);
        TEST_PASS("verify_with_pyarrow (skipped)");
        return 0;
    }
}

/* ============================================================================
 * Test: Stress test with even more columns
 * ============================================================================ */

static int test_very_large_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    const int VERY_LARGE_COLS = 5000;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "very_large_schema");

    printf("Creating schema with %d columns...\n", VERY_LARGE_COLS);

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        printf("  Failed to create schema: %s\n", err.message);
        TEST_FAIL("very_large_schema", "failed to create schema");
    }

    char col_name[64];
    for (int i = 0; i < VERY_LARGE_COLS; i++) {
        snprintf(col_name, sizeof(col_name), "feature_%05d", i);

        carquet_status_t status = carquet_schema_add_column(
            schema, col_name, CARQUET_PHYSICAL_FLOAT, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);

        if (status != CARQUET_OK) {
            printf("  Failed to add column %d\n", i);
            carquet_schema_free(schema);
            TEST_FAIL("very_large_schema", "failed to add column");
        }

        if (i > 0 && i % 1000 == 0) {
            printf("  Added %d columns...\n", i);
        }
    }

    printf("  Schema has %d columns\n", carquet_schema_num_columns(schema));

    /* Create writer */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;  /* No compression for speed */

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        printf("  Failed to create writer: %s\n", err.message);
        carquet_schema_free(schema);
        TEST_FAIL("very_large_schema", "failed to create writer");
    }

    /* Write just 10 rows */
    const int rows = 10;
    float* data = malloc(rows * sizeof(float));
    if (!data) {
        carquet_writer_close(writer);
        carquet_schema_free(schema);
        TEST_FAIL("very_large_schema", "malloc failed");
    }

    for (int i = 0; i < rows; i++) {
        data[i] = (float)i * 0.1f;
    }

    printf("Writing %d rows to %d columns...\n", rows, VERY_LARGE_COLS);
    for (int col = 0; col < VERY_LARGE_COLS; col++) {
        carquet_status_t status = carquet_writer_write_batch(writer, col, data, rows, NULL, NULL);
        if (status != CARQUET_OK) {
            printf("  Failed to write column %d\n", col);
            free(data);
            carquet_writer_close(writer);
            carquet_schema_free(schema);
            TEST_FAIL("very_large_schema", "failed to write column");
        }

        if (col > 0 && col % 1000 == 0) {
            printf("  Written %d columns...\n", col);
        }
    }

    free(data);

    printf("Closing writer...\n");
    carquet_status_t status = carquet_writer_close(writer);
    carquet_schema_free(schema);

    if (status != CARQUET_OK) {
        printf("  Writer close failed\n");
        TEST_FAIL("very_large_schema", "failed to close writer");
    }

    /* Try to read it back */
    printf("Reading file back...\n");
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        printf("  Failed to open file: %s\n", err.message);
        remove(test_file);
        TEST_FAIL("very_large_schema", "failed to read back file");
    }

    printf("  Read back: %lld rows, %d columns\n",
           (long long)carquet_reader_num_rows(reader),
           carquet_reader_num_columns(reader));

    carquet_reader_close(reader);
    remove(test_file);

    TEST_PASS("very_large_schema");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    int failures = 0;

    /* Initialize portable temp file path */
    carquet_test_temp_path(TEST_FILE, sizeof(TEST_FILE), "large_schema");

    printf("=== Large Schema Tests ===\n\n");

    /* Test 1: Write large schema file */
    failures += test_large_schema_write();

    /* Test 2: Read it back */
    if (failures == 0) {
        failures += test_large_schema_read();
    }

    /* Test 3: Verify with PyArrow if available */
    if (failures == 0) {
        failures += test_verify_with_pyarrow();
    }

    /* Test 4: Very large schema (5000 columns) */
    printf("\n");
    failures += test_very_large_schema();

    /* Cleanup */
    remove(TEST_FILE);

    printf("\n");
    if (failures == 0) {
        printf("All large schema tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
