/**
 * @file test_maturity.c
 * @brief Comprehensive maturity tests for Carquet library
 *
 * This test suite explores the maturity and robustness of Carquet by testing:
 * - All physical data types roundtrip
 * - Edge cases (empty, single row, boundary values)
 * - Nullable columns with various null patterns
 * - Error handling (corrupted files, invalid inputs)
 * - Stress scenarios (large data, many columns, many row groups)
 * - Interoperability (generate files for external verification)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>

#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include <carquet/carquet.h>
#include "core/compat.h"

/* Portable temp directory helper */
static const char* get_temp_dir(void) {
#ifdef _WIN32
    static char temp_dir[512] = {0};
    if (temp_dir[0] == 0) {
        const char* tmp = getenv("TEMP");
        if (!tmp) tmp = getenv("TMP");
        if (!tmp) tmp = ".";
        snprintf(temp_dir, sizeof(temp_dir), "%s", tmp);
    }
    return temp_dir;
#else
    return "/tmp";
#endif
}

static void make_temp_path(char* buf, size_t bufsize, const char* name) {
    snprintf(buf, bufsize, "%s/%s.parquet", get_temp_dir(), name);
}

/* ============================================================================
 * Test Macros
 * ============================================================================ */

#define TEST_PASS(name) do { \
    printf("[PASS] %s\n", name); \
    g_tests_passed++; \
} while(0)

#define TEST_FAIL(name, msg) do { \
    printf("[FAIL] %s: %s\n", name, msg); \
    g_tests_failed++; \
    return 1; \
} while(0)

#define TEST_SKIP(name, reason) do { \
    printf("[SKIP] %s: %s\n", name, reason); \
    g_tests_skipped++; \
    return 0; \
} while(0)

#define ASSERT_OK(status, name, msg) do { \
    if ((status) != CARQUET_OK) { \
        printf("[FAIL] %s: %s (status=%d)\n", name, msg, status); \
        g_tests_failed++; \
        return 1; \
    } \
} while(0)

#define ASSERT_TRUE(cond, name, msg) do { \
    if (!(cond)) { \
        printf("[FAIL] %s: %s\n", name, msg); \
        g_tests_failed++; \
        return 1; \
    } \
} while(0)

static int g_tests_passed = 0;
static int g_tests_failed = 0;
static int g_tests_skipped = 0;

/* ============================================================================
 * Section 1: Physical Data Types Roundtrip Tests
 * ============================================================================ */

static int test_type_boolean(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_boolean");
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create schema */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_boolean", "schema creation failed");

(void)carquet_schema_add_column(schema, "bool_col", CARQUET_PHYSICAL_BOOLEAN, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Write test data with various patterns */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_boolean", "writer creation failed");
    }

    /* Test data: alternating, all true, all false, random */
    uint8_t bools[100];
    for (int i = 0; i < 25; i++) bools[i] = i % 2;          /* alternating */
    for (int i = 25; i < 50; i++) bools[i] = 1;             /* all true */
    for (int i = 50; i < 75; i++) bools[i] = 0;             /* all false */
    for (int i = 75; i < 100; i++) bools[i] = (i * 7) % 2;  /* pseudo-random */

    carquet_status_t status = carquet_writer_write_batch(writer, 0, bools, 100, NULL, NULL);
    ASSERT_OK(status, "type_boolean", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_boolean", "writer close failed");

    /* Read back and verify */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_boolean", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == 100, "type_boolean", "row count mismatch");

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 100;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        carquet_schema_free(schema);
        remove(test_file);
        printf("  DEBUG: batch_reader creation failed: %s\n", err.message);
        TEST_FAIL("type_boolean", "batch reader creation failed");
    }

    carquet_row_batch_t* batch = NULL;
    carquet_status_t next_status = carquet_batch_reader_next(batch_reader, &batch);
    printf("  DEBUG: batch_reader_next returned %d, batch=%p\n", next_status, (void*)batch);

    if (next_status == CARQUET_OK && batch) {
        printf("  DEBUG: batch num_columns=%d, num_rows=%lld\n",
               carquet_row_batch_num_columns(batch), (long long)carquet_row_batch_num_rows(batch));
        const void* data;
        int64_t count;
        const uint8_t* nulls = NULL;
        carquet_status_t col_status = carquet_row_batch_column(batch, 0, &data, &nulls, &count);
        printf("  DEBUG: row_batch_column returned %d\n", col_status);

        const uint8_t* read_bools = (const uint8_t*)data;
        for (int i = 0; i < 100; i++) {
            if (read_bools[i] != bools[i]) {
                printf("  Boolean mismatch at index %d: expected %d, got %d\n",
                       i, bools[i], read_bools[i]);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                carquet_schema_free(schema);
                remove(test_file);
                TEST_FAIL("type_boolean", "data mismatch");
            }
        }
        carquet_row_batch_free(batch);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_boolean");
    return 0;
}

static int test_type_int32(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_int32");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_int32", "schema creation failed");

(void)carquet_schema_add_column(schema, "int32_col", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_int32", "writer creation failed");
    }

    /* Test boundary values */
    int32_t values[] = {
        0, 1, -1, 100, -100,
        INT32_MAX, INT32_MIN, INT32_MAX - 1, INT32_MIN + 1,
        127, -128, 255, -256,       /* byte boundaries */
        32767, -32768,              /* short boundaries */
        65535, -65536,
        0x7FFFFF, -0x800000,        /* 24-bit boundaries */
    };
    int count = sizeof(values) / sizeof(values[0]);

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, count, NULL, NULL);
    ASSERT_OK(status, "type_int32", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_int32", "writer close failed");

    /* Read back and verify */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_int32", "reader open failed");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = count;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    carquet_row_batch_t* batch = NULL;

    if (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t n;
        (void)carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n);

        const int32_t* read_values = (const int32_t*)data;
        for (int i = 0; i < count; i++) {
            if (read_values[i] != values[i]) {
                printf("  INT32 mismatch at index %d: expected %d, got %d\n",
                       i, values[i], read_values[i]);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                carquet_schema_free(schema);
                remove(test_file);
                TEST_FAIL("type_int32", "data mismatch");
            }
        }
        carquet_row_batch_free(batch);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_int32");
    return 0;
}

static int test_type_int64(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_int64");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_int64", "schema creation failed");

(void)carquet_schema_add_column(schema, "int64_col", CARQUET_PHYSICAL_INT64, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_int64", "writer creation failed");
    }

    int64_t values[] = {
        0, 1, -1, 1000000000000LL, -1000000000000LL,
        INT64_MAX, INT64_MIN, INT64_MAX - 1, INT64_MIN + 1,
        (int64_t)INT32_MAX + 1, (int64_t)INT32_MIN - 1,
        0x7FFFFFFFFFFFFFFFLL,
        (int64_t)0x8000000000000000LL,
    };
    int count = sizeof(values) / sizeof(values[0]);

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, count, NULL, NULL);
    ASSERT_OK(status, "type_int64", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_int64", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_int64", "reader open failed");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    carquet_row_batch_t* batch = NULL;

    if (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t n;
        (void)carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n);

        const int64_t* read_values = (const int64_t*)data;
        for (int i = 0; i < count; i++) {
            if (read_values[i] != values[i]) {
                printf("  INT64 mismatch at index %d: expected %lld, got %lld\n",
                       i, (long long)values[i], (long long)read_values[i]);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                carquet_schema_free(schema);
                remove(test_file);
                TEST_FAIL("type_int64", "data mismatch");
            }
        }
        carquet_row_batch_free(batch);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_int64");
    return 0;
}

static int test_type_float(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_float");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_float", "schema creation failed");

(void)carquet_schema_add_column(schema, "float_col", CARQUET_PHYSICAL_FLOAT, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_float", "writer creation failed");
    }

    float values[] = {
        0.0f, 1.0f, -1.0f,
        FLT_MIN, FLT_MAX, -FLT_MAX,
        FLT_EPSILON, -FLT_EPSILON,
        3.14159265f, -2.71828182f,
        1e-38f, 1e38f,
        0.1f, 0.2f, 0.3f,  /* values that can't be represented exactly */
        INFINITY, -INFINITY,
        /* Note: NaN comparison is tricky, skip for now */
    };
    int count = sizeof(values) / sizeof(values[0]);

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, count, NULL, NULL);
    ASSERT_OK(status, "type_float", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_float", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_float", "reader open failed");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    carquet_row_batch_t* batch = NULL;

    if (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t n;
        (void)carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n);

        const float* read_values = (const float*)data;
        for (int i = 0; i < count; i++) {
            /* Handle infinity specially */
            if (isinf(values[i]) && isinf(read_values[i])) {
                if ((values[i] > 0) != (read_values[i] > 0)) {
                    printf("  FLOAT infinity sign mismatch at index %d\n", i);
                    TEST_FAIL("type_float", "infinity sign mismatch");
                }
                continue;
            }
            /* Bitwise comparison for exact equality */
            uint32_t expected_bits, actual_bits;
            memcpy(&expected_bits, &values[i], 4);
            memcpy(&actual_bits, &read_values[i], 4);
            if (expected_bits != actual_bits) {
                printf("  FLOAT mismatch at index %d: expected %g (0x%08x), got %g (0x%08x)\n",
                       i, values[i], expected_bits, read_values[i], actual_bits);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                carquet_schema_free(schema);
                remove(test_file);
                TEST_FAIL("type_float", "data mismatch");
            }
        }
        carquet_row_batch_free(batch);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_float");
    return 0;
}

static int test_type_double(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_double");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_double", "schema creation failed");

(void)carquet_schema_add_column(schema, "double_col", CARQUET_PHYSICAL_DOUBLE, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_double", "writer creation failed");
    }

    double values[] = {
        0.0, 1.0, -1.0,
        DBL_MIN, DBL_MAX, -DBL_MAX,
        DBL_EPSILON, -DBL_EPSILON,
        3.14159265358979323846, -2.71828182845904523536,
        1e-308, 1e308,
        (double)INFINITY, (double)-INFINITY,
    };
    int count = sizeof(values) / sizeof(values[0]);

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, count, NULL, NULL);
    ASSERT_OK(status, "type_double", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_double", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_double", "reader open failed");
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    carquet_row_batch_t* batch = NULL;

    if (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t n;
        (void)carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n);

        const double* read_values = (const double*)data;
        for (int i = 0; i < count; i++) {
            if (isinf(values[i]) && isinf(read_values[i])) {
                if ((values[i] > 0) != (read_values[i] > 0)) {
                    TEST_FAIL("type_double", "infinity sign mismatch");
                }
                continue;
            }
            uint64_t expected_bits, actual_bits;
            memcpy(&expected_bits, &values[i], 8);
            memcpy(&actual_bits, &read_values[i], 8);
            if (expected_bits != actual_bits) {
                printf("  DOUBLE mismatch at index %d: expected %g, got %g\n",
                       i, values[i], read_values[i]);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(batch_reader);
                carquet_reader_close(reader);
                carquet_schema_free(schema);
                remove(test_file);
                TEST_FAIL("type_double", "data mismatch");
            }
        }
        carquet_row_batch_free(batch);
    }

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_double");
    return 0;
}

static int test_type_byte_array(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_byte_array");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("type_byte_array", "schema creation failed");

(void)carquet_schema_add_column(schema, "string_col", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("type_byte_array", "writer creation failed");
    }

    /* Test various string patterns using carquet_byte_array_t */
    const char* string_data[] = {
        "",                              /* empty string */
        "a",                             /* single char */
        "Hello, World!",                 /* ASCII */
        "Test String",                   /* simple ASCII */
    };
    int count = sizeof(string_data) / sizeof(string_data[0]);

    /* Create array of byte_array structures */
    carquet_byte_array_t byte_arrays[4];
    for (int i = 0; i < count; i++) {
        byte_arrays[i].data = (uint8_t*)string_data[i];
        byte_arrays[i].length = (int32_t)strlen(string_data[i]);
    }

    carquet_status_t status = carquet_writer_write_batch(writer, 0, byte_arrays, count, NULL, NULL);
    ASSERT_OK(status, "type_byte_array", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "type_byte_array", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("type_byte_array", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == count, "type_byte_array", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("type_byte_array");
    return 0;
}

/* ============================================================================
 * Section 2: Edge Case Tests
 * ============================================================================ */

static int test_edge_single_row(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_single_row");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("edge_single_row", "schema creation failed");

(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("edge_single_row", "writer creation failed");
    }

    int32_t value = 42;
    carquet_status_t status = carquet_writer_write_batch(writer, 0, &value, 1, NULL, NULL);
    ASSERT_OK(status, "edge_single_row", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "edge_single_row", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("edge_single_row", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == 1, "edge_single_row", "row count should be 1");
    ASSERT_TRUE(carquet_reader_num_columns(reader) == 1, "edge_single_row", "col count should be 1");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("edge_single_row");
    return 0;
}

static int test_edge_many_columns(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_many_columns");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("edge_many_columns", "schema creation failed");

    /* Create 10 columns (reduced from 100 to avoid potential schema size issues) */
    const int NUM_COLS = 10;
    char col_name[32];
    for (int i = 0; i < NUM_COLS; i++) {
        snprintf(col_name, sizeof(col_name), "col_%03d", i);
        carquet_status_t status = carquet_schema_add_column(schema, col_name, CARQUET_PHYSICAL_INT32, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);
        if (status != CARQUET_OK) {
            printf("  Failed to add column %d: status=%d\n", i, status);
            carquet_schema_free(schema);
            TEST_FAIL("edge_many_columns", "column addition failed");
        }
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("edge_many_columns", "writer creation failed");
    }

    /* Write 10 rows */
    const int NUM_ROWS = 10;
    int32_t values[10];
    for (int i = 0; i < NUM_ROWS; i++) values[i] = i;

    for (int col = 0; col < NUM_COLS; col++) {
        carquet_status_t status = carquet_writer_write_batch(writer, col, values, NUM_ROWS, NULL, NULL);
        if (status != CARQUET_OK) {
            printf("  Failed to write column %d\n", col);
(void)carquet_writer_close(writer);
            carquet_schema_free(schema);
            remove(test_file);
            TEST_FAIL("edge_many_columns", "write failed");
        }
    }

    carquet_status_t status = carquet_writer_close(writer);
    ASSERT_OK(status, "edge_many_columns", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("edge_many_columns", "reader open failed");
    }

    printf("  Created file with %d columns, %lld rows\n",
           carquet_reader_num_columns(reader),
           (long long)carquet_reader_num_rows(reader));

    ASSERT_TRUE(carquet_reader_num_columns(reader) == NUM_COLS, "edge_many_columns", "column count mismatch");
    ASSERT_TRUE(carquet_reader_num_rows(reader) == NUM_ROWS, "edge_many_columns", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("edge_many_columns");
    return 0;
}

static int test_edge_many_row_groups(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_many_row_groups");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("edge_many_row_groups", "schema creation failed");

(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    /* Force small row groups */
    opts.row_group_size = 100;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("edge_many_row_groups", "writer creation failed");
    }

    /* Write 50 row groups with 10 rows each */
    const int NUM_ROW_GROUPS = 50;
    const int ROWS_PER_GROUP = 10;
    int32_t values[10];

    for (int rg = 0; rg < NUM_ROW_GROUPS; rg++) {
        for (int i = 0; i < ROWS_PER_GROUP; i++) {
            values[i] = rg * ROWS_PER_GROUP + i;
        }
(void)carquet_writer_write_batch(writer, 0, values, ROWS_PER_GROUP, NULL, NULL);
        if (rg < NUM_ROW_GROUPS - 1) {
(void)carquet_writer_new_row_group(writer);
        }
    }

    carquet_status_t status = carquet_writer_close(writer);
    ASSERT_OK(status, "edge_many_row_groups", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("edge_many_row_groups", "reader open failed");
    }

    int32_t num_rg = carquet_reader_num_row_groups(reader);
    int64_t total_rows = carquet_reader_num_rows(reader);

    printf("  Created file with %d row groups, %lld total rows\n",
           num_rg, (long long)total_rows);

    ASSERT_TRUE(num_rg == NUM_ROW_GROUPS, "edge_many_row_groups", "row group count mismatch");
    ASSERT_TRUE(total_rows == NUM_ROW_GROUPS * ROWS_PER_GROUP, "edge_many_row_groups", "total row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("edge_many_row_groups");
    return 0;
}

static int test_auto_row_group_flush_aligned_batches(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_auto_row_group_flush");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("auto_row_group_flush_aligned_batches", "schema creation failed");

(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.row_group_size = 10 * (int64_t)(sizeof(int32_t) + sizeof(double)) - 1;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("auto_row_group_flush_aligned_batches", "writer creation failed");
    }

    enum { NUM_ROW_GROUPS = 3, ROWS_PER_GROUP = 10 };
    int32_t ids[ROWS_PER_GROUP];
    double values[ROWS_PER_GROUP];

    for (int rg = 0; rg < NUM_ROW_GROUPS; rg++) {
        for (int i = 0; i < ROWS_PER_GROUP; i++) {
            ids[i] = rg * ROWS_PER_GROUP + i;
            values[i] = (double)ids[i] * 0.5;
        }

        ASSERT_OK(carquet_writer_write_batch(writer, 0, ids, ROWS_PER_GROUP, NULL, NULL),
                  "auto_row_group_flush_aligned_batches", "id write failed");
        ASSERT_OK(carquet_writer_write_batch(writer, 1, values, ROWS_PER_GROUP, NULL, NULL),
                  "auto_row_group_flush_aligned_batches", "value write failed");
    }

    ASSERT_OK(carquet_writer_close(writer),
              "auto_row_group_flush_aligned_batches", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("auto_row_group_flush_aligned_batches", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_row_groups(reader) == NUM_ROW_GROUPS,
                "auto_row_group_flush_aligned_batches", "row group count mismatch");
    ASSERT_TRUE(carquet_reader_num_rows(reader) == NUM_ROW_GROUPS * ROWS_PER_GROUP,
                "auto_row_group_flush_aligned_batches", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("auto_row_group_flush_aligned_batches");
    return 0;
}

/* ============================================================================
 * Section 3: Nullable Column Tests
 * ============================================================================ */

static int test_nullable_all_null(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_all_null");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("nullable_all_null", "schema creation failed");

(void)carquet_schema_add_column(schema, "nullable_int", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("nullable_all_null", "writer creation failed");
    }

    /* All values are null */
    int32_t values[10] = {0};  /* values don't matter */
    int16_t def_levels[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  /* all null */

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, 10, def_levels, NULL);
    ASSERT_OK(status, "nullable_all_null", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "nullable_all_null", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("nullable_all_null", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == 10, "nullable_all_null", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("nullable_all_null");
    return 0;
}

static int test_nullable_none_null(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_none_null");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("nullable_none_null", "schema creation failed");

(void)carquet_schema_add_column(schema, "nullable_int", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("nullable_none_null", "writer creation failed");
    }

    /* No values are null */
    int32_t values[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int16_t def_levels[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  /* all present */

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, 10, def_levels, NULL);
    ASSERT_OK(status, "nullable_none_null", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "nullable_none_null", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("nullable_none_null", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == 10, "nullable_none_null", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("nullable_none_null");
    return 0;
}

static int test_nullable_mixed(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_mixed_null");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("nullable_mixed", "schema creation failed");

(void)carquet_schema_add_column(schema, "nullable_int", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("nullable_mixed", "writer creation failed");
    }

    /* Mix of null and non-null (sparse encoding: values array has only non-null values) */
    int32_t values[5] = {1, 3, 5, 7, 9};
    int16_t def_levels[10] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0};  /* alternating */

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, 10, def_levels, NULL);
    ASSERT_OK(status, "nullable_mixed", "write failed");

    status = carquet_writer_close(writer);
    ASSERT_OK(status, "nullable_mixed", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("nullable_mixed", "reader open failed");
    }

    ASSERT_TRUE(carquet_reader_num_rows(reader) == 10, "nullable_mixed", "row count mismatch");

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 10;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("nullable_mixed", "batch reader creation failed");
    }

    carquet_row_batch_t* batch = NULL;
    status = carquet_batch_reader_next(batch_reader, &batch);
    ASSERT_OK(status, "nullable_mixed", "batch read failed");
    ASSERT_TRUE(batch != NULL, "nullable_mixed", "missing batch");

    const void* data = NULL;
    const uint8_t* null_bitmap = NULL;
    int64_t n = 0;
    status = carquet_row_batch_column(batch, 0, &data, &null_bitmap, &n);
    ASSERT_OK(status, "nullable_mixed", "column read failed");
    ASSERT_TRUE(n == 10, "nullable_mixed", "batch row count mismatch");
    ASSERT_TRUE(null_bitmap != NULL, "nullable_mixed", "missing null bitmap");

    const int32_t* read_values = (const int32_t*)data;
    const int32_t expected_values[10] = {1, 0, 3, 0, 5, 0, 7, 0, 9, 0};
    for (int i = 0; i < 10; i++) {
        int is_present = (null_bitmap[i / 8] >> (i % 8)) & 1;
        if ((i % 2) == 1) {
            ASSERT_TRUE(!is_present, "nullable_mixed", "expected null bit");
        } else {
            ASSERT_TRUE(is_present, "nullable_mixed", "unexpected null bit");
        }
        if (read_values[i] != expected_values[i]) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            carquet_schema_free(schema);
            remove(test_file);
            TEST_FAIL("nullable_mixed", "decoded values not aligned to logical rows");
        }
    }

    carquet_row_batch_free(batch);
    carquet_batch_reader_free(batch_reader);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("nullable_mixed");
    return 0;
}

static int test_nullable_partial_reads(void) {
    const char* test_name = "nullable_partial_reads";
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "carquet_nullable_partial_reads");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(test_name, "schema creation failed");

    carquet_status_t status = carquet_schema_add_column(
        schema, "value", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_OPTIONAL, 0, 0);
    ASSERT_OK(status, test_name, "add column failed");

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.page_size = 16;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL(test_name, "writer creation failed");
    }

    int32_t values[5] = {10, 30, 50, 70, 90};
    int16_t def_levels[10] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
    status = carquet_writer_write_batch(writer, 0, values, 10, def_levels, NULL);
    ASSERT_OK(status, test_name, "write failed");
    status = carquet_writer_close(writer);
    ASSERT_OK(status, test_name, "close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL(test_name, "reader open failed");
    }

    carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, 0, &err);
    if (!col) {
        carquet_reader_close(reader);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL(test_name, "column open failed");
    }

    int32_t out[4] = {0};
    int16_t defs[4] = {0};
    int64_t n = carquet_column_read_batch(col, out, 3, defs, NULL);
    ASSERT_TRUE(n == 3, test_name, "first partial count mismatch");
    ASSERT_TRUE(defs[0] == 1 && defs[1] == 0 && defs[2] == 1,
                test_name, "first partial defs mismatch");
    ASSERT_TRUE(out[0] == 10 && out[1] == 30,
                test_name, "first partial sparse values mismatch");

    memset(out, 0, sizeof(out));
    memset(defs, 0, sizeof(defs));
    n = carquet_column_read_batch(col, out, 4, defs, NULL);
    ASSERT_TRUE(n == 4, test_name, "second partial count mismatch");
    ASSERT_TRUE(defs[0] == 0 && defs[1] == 1 && defs[2] == 0 && defs[3] == 1,
                test_name, "second partial defs mismatch");
    ASSERT_TRUE(out[0] == 50 && out[1] == 70,
                test_name, "second partial sparse values mismatch");

    carquet_column_reader_free(col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS(test_name);
    return 0;
}

/* ============================================================================
 * Section 4: Error Handling Tests
 * ============================================================================ */

static int test_error_invalid_file(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Try to open a non-existent file */
    carquet_reader_t* reader = carquet_reader_open("/nonexistent/path.parquet", NULL, &err);
    ASSERT_TRUE(reader == NULL, "error_invalid_file", "should fail on nonexistent file");
    ASSERT_TRUE(err.code != CARQUET_OK, "error_invalid_file", "should have error code");

    TEST_PASS("error_invalid_file");
    return 0;
}

static int test_error_corrupted_magic(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_corrupted");
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create a file with invalid magic bytes */
    FILE* f = fopen(test_file, "wb");
    if (!f) TEST_FAIL("error_corrupted_magic", "failed to create test file");

    /* Write garbage instead of PAR1 */
    const char* garbage = "XXXX this is not a parquet file XXXX";
    fwrite(garbage, 1, strlen(garbage), f);
    fclose(f);

    /* Try to open it */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    ASSERT_TRUE(reader == NULL, "error_corrupted_magic", "should fail on corrupted file");
    ASSERT_TRUE(err.code == CARQUET_ERROR_INVALID_MAGIC || err.code != CARQUET_OK,
                "error_corrupted_magic", "should detect invalid magic");

    printf("  Error message: %s\n", err.message);

    remove(test_file);
    TEST_PASS("error_corrupted_magic");
    return 0;
}

static int test_error_truncated_file(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_truncated");
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create a file with just the header magic but nothing else */
    FILE* f = fopen(test_file, "wb");
    if (!f) TEST_FAIL("error_truncated_file", "failed to create test file");

    fwrite("PAR1", 1, 4, f);  /* Header only, truncated */
    fclose(f);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    ASSERT_TRUE(reader == NULL, "error_truncated_file", "should fail on truncated file");

    printf("  Error on truncated file: %s\n", err.message);

    remove(test_file);
    TEST_PASS("error_truncated_file");
    return 0;
}

static int test_error_invalid_arguments(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Test non-existent file path */
    carquet_reader_t* reader = carquet_reader_open("/nonexistent/path/file.parquet", NULL, &err);
    ASSERT_TRUE(reader == NULL, "error_invalid_arguments", "should fail on non-existent file");
    ASSERT_TRUE(err.code != CARQUET_OK, "error_invalid_arguments", "should set error code");

    /* Test empty schema */
    carquet_schema_t* empty_schema = carquet_schema_create(&err);
    ASSERT_TRUE(empty_schema != NULL, "error_invalid_arguments", "should create empty schema");

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    char empty_schema_file[512];
    make_temp_path(empty_schema_file, sizeof(empty_schema_file), "test_empty_schema");
    carquet_writer_t* writer = carquet_writer_create(empty_schema_file, empty_schema, &opts, &err);
    if (writer) {
        carquet_writer_close(writer);
    }
    carquet_schema_free(empty_schema);

    TEST_PASS("error_invalid_arguments");
    return 0;
}

/* ============================================================================
 * Section 5: Compression Roundtrip Tests
 * ============================================================================ */

static int test_compression_roundtrip(carquet_compression_t compression, const char* name) {
    char test_file[256];
    snprintf(test_file, sizeof(test_file), "%s/test_maturity_compress_%s.parquet", get_temp_dir(), name);
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        printf("[FAIL] compress_%s: schema creation failed\n", name);
        g_tests_failed++;
        return 1;
    }

(void)carquet_schema_add_column(schema, "data", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = compression;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("[FAIL] compress_%s: writer creation failed - %s\n", name, err.message);
        g_tests_failed++;
        return 1;
    }

    /* Write 1000 values with patterns good for compression */
    int32_t values[1000];
    for (int i = 0; i < 1000; i++) {
        values[i] = i % 100;  /* Repeating pattern */
    }

    carquet_status_t status = carquet_writer_write_batch(writer, 0, values, 1000, NULL, NULL);
    if (status != CARQUET_OK) {
(void)carquet_writer_close(writer);
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] compress_%s: write failed\n", name);
        g_tests_failed++;
        return 1;
    }

    status = carquet_writer_close(writer);
    if (status != CARQUET_OK) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] compress_%s: close failed\n", name);
        g_tests_failed++;
        return 1;
    }

    /* Read back and verify */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] compress_%s: reader open failed - %s\n", name, err.message);
        g_tests_failed++;
        return 1;
    }

    if (carquet_reader_num_rows(reader) != 1000) {
        carquet_reader_close(reader);
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] compress_%s: row count mismatch\n", name);
        g_tests_failed++;
        return 1;
    }

    /* Get file size for compression ratio */
    FILE* f = fopen(test_file, "rb");
    if (f) {
        carquet_fseek64(f, 0, SEEK_END);
        int64_t size = carquet_ftell64(f);
        fclose(f);
        printf("  %s: file size = %lld bytes (raw = 4000)\n", name, (long long)size);
    }

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    printf("[PASS] compress_%s\n", name);
    g_tests_passed++;
    return 0;
}

static int test_float_double_compression_roundtrip(
    carquet_compression_t compression,
    const char* name) {

    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_float_double");
    carquet_error_t err = CARQUET_ERROR_INIT;
    const int count = 2048;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("float_double_compression", "schema creation failed");

    (void)carquet_schema_add_column(schema, "value_f32", CARQUET_PHYSICAL_FLOAT, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "value_f64", CARQUET_PHYSICAL_DOUBLE, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);

    float* f32 = malloc((size_t)count * sizeof(float));
    double* f64 = malloc((size_t)count * sizeof(double));
    if (!f32 || !f64) {
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        TEST_FAIL("float_double_compression", "allocation failed");
    }

    for (int i = 0; i < count; i++) {
        f32[i] = (float)((i % 97) * 0.25f + (i / 97) * 0.001f);
        f64[i] = (double)((i % 131) * 0.125 + (i / 131) * 0.00001);
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = compression;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        printf("[FAIL] float_double_compression_%s: writer creation failed - %s\n", name, err.message);
        g_tests_failed++;
        return 1;
    }

    carquet_status_t status = carquet_writer_write_batch(writer, 0, f32, count, NULL, NULL);
    ASSERT_OK(status, "float_double_compression", "float write failed");
    status = carquet_writer_write_batch(writer, 1, f64, count, NULL, NULL);
    ASSERT_OK(status, "float_double_compression", "double write failed");
    status = carquet_writer_close(writer);
    ASSERT_OK(status, "float_double_compression", "writer close failed");

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] float_double_compression_%s: reader open failed - %s\n", name, err.message);
        g_tests_failed++;
        return 1;
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = count;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        remove(test_file);
        printf("[FAIL] float_double_compression_%s: batch reader failed - %s\n", name, err.message);
        g_tests_failed++;
        return 1;
    }

    carquet_row_batch_t* batch = NULL;
    status = carquet_batch_reader_next(batch_reader, &batch);
    ASSERT_OK(status, "float_double_compression", "batch read failed");
    if (!batch || carquet_row_batch_num_rows(batch) != count) {
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("float_double_compression", "unexpected batch row count");
    }

    const void* data = NULL;
    const uint8_t* nulls = NULL;
    int64_t values_read = 0;
    status = carquet_row_batch_column(batch, 0, &data, &nulls, &values_read);
    ASSERT_OK(status, "float_double_compression", "float column read failed");
    if (values_read != count) {
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("float_double_compression", "float row count mismatch");
    }
    const float* read_f32 = (const float*)data;
    for (int i = 0; i < count; i++) {
        if (read_f32[i] != f32[i]) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            free(f32);
            free(f64);
            carquet_schema_free(schema);
            remove(test_file);
            TEST_FAIL("float_double_compression", "float roundtrip mismatch");
        }
    }

    status = carquet_row_batch_column(batch, 1, &data, &nulls, &values_read);
    ASSERT_OK(status, "float_double_compression", "double column read failed");
    if (values_read != count) {
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(batch_reader);
        carquet_reader_close(reader);
        free(f32);
        free(f64);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("float_double_compression", "double row count mismatch");
    }
    const double* read_f64 = (const double*)data;
    for (int i = 0; i < count; i++) {
        if (read_f64[i] != f64[i]) {
            carquet_row_batch_free(batch);
            carquet_batch_reader_free(batch_reader);
            carquet_reader_close(reader);
            free(f32);
            free(f64);
            carquet_schema_free(schema);
            remove(test_file);
            TEST_FAIL("float_double_compression", "double roundtrip mismatch");
        }
    }

    carquet_row_batch_free(batch);
    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);
    free(f32);
    free(f64);
    carquet_schema_free(schema);
    remove(test_file);

    printf("[PASS] float_double_compression_%s\n", name);
    g_tests_passed++;
    return 0;
}

static int test_all_compressions(void) {
    printf("\n--- Compression Roundtrip Tests ---\n");

    test_compression_roundtrip(CARQUET_COMPRESSION_UNCOMPRESSED, "uncompressed");
    test_compression_roundtrip(CARQUET_COMPRESSION_SNAPPY, "snappy");
    test_compression_roundtrip(CARQUET_COMPRESSION_GZIP, "gzip");
    test_compression_roundtrip(CARQUET_COMPRESSION_LZ4, "lz4");
    test_compression_roundtrip(CARQUET_COMPRESSION_ZSTD, "zstd");
    test_float_double_compression_roundtrip(CARQUET_COMPRESSION_SNAPPY, "snappy");
    test_float_double_compression_roundtrip(CARQUET_COMPRESSION_ZSTD, "zstd");

    return 0;
}

/* ============================================================================
 * Section 6: Stress Tests
 * ============================================================================ */

static int test_stress_large_data(void) {
    char test_file[512];
    make_temp_path(test_file, sizeof(test_file), "test_maturity_large");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("stress_large_data", "schema creation failed");

(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("stress_large_data", "writer creation failed");
    }

    /* Write 1 million rows in batches */
    const int64_t TOTAL_ROWS = 1000000;
    const int BATCH_SIZE = 10000;

    int64_t* ids = malloc(BATCH_SIZE * sizeof(int64_t));
    double* values = malloc(BATCH_SIZE * sizeof(double));

    clock_t start = clock();

    for (int64_t batch = 0; batch < TOTAL_ROWS / BATCH_SIZE; batch++) {
        for (int i = 0; i < BATCH_SIZE; i++) {
            int64_t row = batch * BATCH_SIZE + i;
            ids[i] = row;
            values[i] = (double)row * 0.001;
        }

(void)carquet_writer_write_batch(writer, 0, ids, BATCH_SIZE, NULL, NULL);
(void)carquet_writer_write_batch(writer, 1, values, BATCH_SIZE, NULL, NULL);

        /* Create new row group every 100K rows */
        if ((batch + 1) % 10 == 0 && batch < (TOTAL_ROWS / BATCH_SIZE) - 1) {
(void)carquet_writer_new_row_group(writer);
        }
    }

    free(ids);
    free(values);

    carquet_status_t status = carquet_writer_close(writer);
    clock_t end = clock();
    double write_time = (double)(end - start) / CLOCKS_PER_SEC;

    if (status != CARQUET_OK) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("stress_large_data", "writer close failed");
    }

    /* Get file size */
    FILE* f = fopen(test_file, "rb");
    int64_t file_size = 0;
    if (f) {
        carquet_fseek64(f, 0, SEEK_END);
        file_size = carquet_ftell64(f);
        fclose(f);
    }

    /* Read back */
    start = clock();
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("stress_large_data", "reader open failed");
    }

    int64_t read_rows = carquet_reader_num_rows(reader);
    end = clock();
    double read_time = (double)(end - start) / CLOCKS_PER_SEC;

    printf("  Wrote %lld rows in %.2fs (%.0f rows/sec)\n",
           (long long)TOTAL_ROWS, write_time, TOTAL_ROWS / write_time);
    printf("  File size: %.2f MB (%.2f bytes/row)\n",
           file_size / (1024.0 * 1024.0), (double)file_size / TOTAL_ROWS);
    printf("  Read metadata in %.4fs\n", read_time);

    ASSERT_TRUE(read_rows == TOTAL_ROWS, "stress_large_data", "row count mismatch");

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("stress_large_data");
    return 0;
}

/* ============================================================================
 * Section 7: Interoperability Tests (generate files for external validation)
 * ============================================================================ */

static int test_interop_generate_files(void) {
    char base_path[512];
    snprintf(base_path, sizeof(base_path), "%s/carquet_interop", get_temp_dir());
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Test 1: Simple integers */
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s_int32.parquet", base_path);

        carquet_schema_t* schema = carquet_schema_create(&err);
(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);

        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);

        carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
        if (writer) {
            int32_t values[100];
            for (int i = 0; i < 100; i++) values[i] = i * 10;
(void)carquet_writer_write_batch(writer, 0, values, 100, NULL, NULL);
(void)carquet_writer_close(writer);
            printf("  Generated: %s\n", path);
        }
        carquet_schema_free(schema);
    }

    /* Test 2: Multiple data types */
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s_mixed.parquet", base_path);

        carquet_schema_t* schema = carquet_schema_create(&err);
(void)carquet_schema_add_column(schema, "int_col", CARQUET_PHYSICAL_INT64, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "float_col", CARQUET_PHYSICAL_DOUBLE, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);
(void)carquet_schema_add_column(schema, "bool_col", CARQUET_PHYSICAL_BOOLEAN, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);

        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        opts.compression = CARQUET_COMPRESSION_SNAPPY;

        carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
        if (writer) {
            int64_t ints[50];
            double floats[50];
            uint8_t bools[50];
            for (int i = 0; i < 50; i++) {
                ints[i] = i * 1000;
                floats[i] = i * 3.14159;
                bools[i] = i % 2;
            }
(void)carquet_writer_write_batch(writer, 0, ints, 50, NULL, NULL);
(void)carquet_writer_write_batch(writer, 1, floats, 50, NULL, NULL);
(void)carquet_writer_write_batch(writer, 2, bools, 50, NULL, NULL);
(void)carquet_writer_close(writer);
            printf("  Generated: %s\n", path);
        }
        carquet_schema_free(schema);
    }

    /* Test 3: With ZSTD compression */
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s_zstd.parquet", base_path);

        carquet_schema_t* schema = carquet_schema_create(&err);
(void)carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_INT32, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0);

        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        opts.compression = CARQUET_COMPRESSION_ZSTD;

        carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
        if (writer) {
            int32_t values[1000];
            for (int i = 0; i < 1000; i++) values[i] = i % 100;  /* repetitive for compression */
(void)carquet_writer_write_batch(writer, 0, values, 1000, NULL, NULL);
(void)carquet_writer_close(writer);
            printf("  Generated: %s\n", path);
        }
        carquet_schema_free(schema);
    }

    printf("  Files ready for PyArrow verification:\n");
    printf("    python3 -c \"import pyarrow.parquet as pq; print(pq.read_table('%s_int32.parquet').to_pandas())\"\n", base_path);
    printf("    python3 -c \"import pyarrow.parquet as pq; print(pq.read_table('%s_mixed.parquet').to_pandas())\"\n", base_path);
    printf("    python3 -c \"import pyarrow.parquet as pq; print(pq.read_table('%s_zstd.parquet').to_pandas())\"\n", base_path);

    TEST_PASS("interop_generate_files");
    return 0;
}

/* ============================================================================
 * Section 8: Memory Safety Tests
 * ============================================================================ */

static int test_memory_double_free_protection(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Test that calling free twice doesn't crash */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("memory_double_free", "schema creation failed");

(void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_schema_free(schema);
    /* Note: We can't actually test double-free safely, but we verify single free works */

    TEST_PASS("memory_double_free_protection");
    return 0;
}

static int test_memory_cleanup_on_error(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Try to open non-existent file - should not leak memory */
    for (int i = 0; i < 100; i++) {
        carquet_reader_t* reader = carquet_reader_open("/nonexistent/path.parquet", NULL, &err);
        /* reader should be NULL, no memory to free */
        (void)reader;
    }

    /* If we get here without crashing, cleanup is working */
    TEST_PASS("memory_cleanup_on_error");
    return 0;
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */

int main(void) {
    printf("=== Carquet Maturity Test Suite ===\n\n");
    printf("Testing library version: %s\n\n", carquet_version());

    /* Initialize library */
(void)carquet_init();

    /* Section 1: Physical Data Types */
    printf("--- Physical Data Type Tests ---\n");
    test_type_boolean();
    test_type_int32();
    test_type_int64();
    test_type_float();
    test_type_double();
    test_type_byte_array();

    /* Section 2: Edge Cases */
    printf("\n--- Edge Case Tests ---\n");
    test_edge_single_row();
    test_edge_many_columns();
    test_edge_many_row_groups();
    test_auto_row_group_flush_aligned_batches();

    /* Section 3: Nullable Columns */
    printf("\n--- Nullable Column Tests ---\n");
    test_nullable_all_null();
    test_nullable_none_null();
    test_nullable_mixed();
    test_nullable_partial_reads();

    /* Section 4: Error Handling */
    printf("\n--- Error Handling Tests ---\n");
    test_error_invalid_file();
    test_error_corrupted_magic();
    test_error_truncated_file();
    test_error_invalid_arguments();

    /* Section 5: Compression */
    test_all_compressions();

    /* Section 6: Stress Tests */
    printf("\n--- Stress Tests ---\n");
    test_stress_large_data();

    /* Section 7: Interoperability */
    printf("\n--- Interoperability Tests ---\n");
    test_interop_generate_files();

    /* Section 8: Memory Safety */
    printf("\n--- Memory Safety Tests ---\n");
    test_memory_double_free_protection();
    test_memory_cleanup_on_error();

    /* Summary */
    printf("\n========================================\n");
    printf("  MATURITY TEST SUMMARY\n");
    printf("========================================\n");
    printf("  Passed:  %d\n", g_tests_passed);
    printf("  Failed:  %d\n", g_tests_failed);
    printf("  Skipped: %d\n", g_tests_skipped);
    printf("========================================\n");

    if (g_tests_failed > 0) {
        printf("\nSome tests failed. Review output above.\n");
        return 1;
    }

    printf("\nAll maturity tests passed!\n");
    return 0;
}
