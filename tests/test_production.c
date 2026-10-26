/**
 * @file test_production.c
 * @brief Tests for production-ready features
 *
 * Tests for:
 * - Column projection (batch reader)
 * - Row group statistics
 * - Predicate pushdown / row group filtering
 * - Memory-mapped I/O
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

#include <carquet/carquet.h>
#include "core/compat.h"
#include "test_helpers.h"

#define NUM_ROWS 10000
#define NUM_ROW_GROUPS 10

static char TEST_FILE[512];

/* ============================================================================
 * Helper: Create test file with multiple row groups
 * ============================================================================
 */

static int create_test_file(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create schema with multiple columns */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return -1;

    (void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "category", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "score", CARQUET_PHYSICAL_FLOAT, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Writer options - small row groups for testing */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;
    opts.row_group_size = (NUM_ROWS / NUM_ROW_GROUPS) * 32;  /* Force multiple row groups */

    carquet_writer_t* writer = carquet_writer_create(TEST_FILE, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        return -1;
    }

    /* Generate test data */
    int32_t* ids = malloc(NUM_ROWS * sizeof(int32_t));
    double* values = malloc(NUM_ROWS * sizeof(double));
    int32_t* categories = malloc(NUM_ROWS * sizeof(int32_t));
    float* scores = malloc(NUM_ROWS * sizeof(float));

    for (int i = 0; i < NUM_ROWS; i++) {
        ids[i] = i;
        values[i] = (double)i * 1.5;
        categories[i] = i % 10;  /* 0-9 categories */
        scores[i] = (float)(i % 100) / 10.0f;  /* 0.0 - 9.9 */
    }

    /* Write data in chunks to create multiple row groups */
    int rows_per_group = NUM_ROWS / NUM_ROW_GROUPS;
    for (int g = 0; g < NUM_ROW_GROUPS; g++) {
        int offset = g * rows_per_group;

        (void)carquet_writer_write_batch(writer, 0, ids + offset, rows_per_group, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 1, values + offset, rows_per_group, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 2, categories + offset, rows_per_group, NULL, NULL);
        (void)carquet_writer_write_batch(writer, 3, scores + offset, rows_per_group, NULL, NULL);

        if (g < NUM_ROW_GROUPS - 1) {
            (void)carquet_writer_new_row_group(writer);
        }
    }

    free(ids);
    free(values);
    free(categories);
    free(scores);

    carquet_status_t status = carquet_writer_close(writer);
    carquet_schema_free(schema);

    return (status == CARQUET_OK) ? 0 : -1;
}

/* ============================================================================
 * Test: Column Projection
 * ============================================================================
 */

static int test_column_projection(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open file */
    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, NULL, &err);
    if (!reader) {
        printf("  Failed to open file: %s\n", err.message);
        TEST_FAIL("column_projection", "failed to open file");
    }

    /* Verify file has expected columns */
    assert(carquet_reader_num_columns(reader) == 4);

    /* Test 1: Project only 2 columns by index */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);

    int32_t proj_cols[] = {0, 2};  /* id and category only */
    config.column_indices = proj_cols;
    config.num_columns = 2;
    config.batch_size = 1000;

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        printf("  Failed to create batch reader: %s\n", err.message);
        TEST_FAIL("column_projection", "failed to create batch reader");
    }

    /* Read all batches and verify */
    int64_t total_rows = 0;
    carquet_row_batch_t* batch = NULL;

    while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        /* Verify batch has only 2 columns */
        assert(carquet_row_batch_num_columns(batch) == 2);

        int64_t batch_rows = carquet_row_batch_num_rows(batch);
        total_rows += batch_rows;

        /* Verify we can access the projected columns */
        const void* data;
        const uint8_t* null_bitmap;
        int64_t num_values;

        assert(carquet_row_batch_column(batch, 0, &data, &null_bitmap, &num_values) == CARQUET_OK);
        assert(data != NULL);
        assert(num_values == batch_rows);

        carquet_status_t status = carquet_row_batch_column(batch, 1, &data, &null_bitmap, &num_values);
        assert(status == CARQUET_OK);
        (void)status;
        assert(data != NULL);

        carquet_row_batch_free(batch);
        batch = NULL;
    }

    printf("  Read %lld rows with 2-column projection\n", (long long)total_rows);
    assert(total_rows == NUM_ROWS);

    carquet_batch_reader_free(batch_reader);

    /* Test 2: Project by column names */
    const char* col_names[] = {"value", "score"};
    carquet_batch_reader_config_init(&config);
    config.column_names = col_names;
    config.num_column_names = 2;
    config.batch_size = 2000;

    batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (!batch_reader) {
        carquet_reader_close(reader);
        printf("  Failed to create batch reader by name: %s\n", err.message);
        TEST_FAIL("column_projection", "failed to create batch reader by name");
    }

    total_rows = 0;
    while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
        assert(carquet_row_batch_num_columns(batch) == 2);
        total_rows += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    printf("  Read %lld rows with name-based projection\n", (long long)total_rows);
    assert(total_rows == NUM_ROWS);

    carquet_batch_reader_free(batch_reader);
    carquet_reader_close(reader);

    TEST_PASS("column_projection");
    return 0;
}

/* ============================================================================
 * Test: Row Group Statistics
 * ============================================================================
 */

static int test_row_group_statistics(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, NULL, &err);
    if (!reader) {
        TEST_FAIL("row_group_statistics", "failed to open file");
    }

    int32_t num_row_groups = carquet_reader_num_row_groups(reader);
    printf("  File has %d row groups\n", num_row_groups);

    if (num_row_groups != NUM_ROW_GROUPS) {
        TEST_FAIL("row_group_statistics", "unexpected row group count");
    }

    int rows_per_group = NUM_ROWS / NUM_ROW_GROUPS;

    /* Statistics must be present for every row group across every column. */
    for (int32_t rg = 0; rg < num_row_groups; rg++) {
        carquet_column_statistics_t stats;
        carquet_status_t status = carquet_reader_column_statistics(reader, rg, 0, &stats);
        if (status != CARQUET_OK) {
            TEST_FAIL("row_group_statistics", "id stats not readable");
        }
        if (!stats.has_min_max) {
            TEST_FAIL("row_group_statistics", "id min/max not written");
        }
        if (!stats.has_null_count || stats.null_count != 0) {
            TEST_FAIL("row_group_statistics", "id null_count not zero/missing");
        }
        if (stats.num_values != rows_per_group) {
            TEST_FAIL("row_group_statistics", "id num_values mismatch");
        }

        int32_t min_val = *(const int32_t*)stats.min_value;
        int32_t max_val = *(const int32_t*)stats.max_value;
        int32_t expected_min = rg * rows_per_group;
        int32_t expected_max = expected_min + rows_per_group - 1;
        if (min_val != expected_min || max_val != expected_max) {
            printf("  Row group %d: got [%d, %d] expected [%d, %d]\n",
                   rg, min_val, max_val, expected_min, expected_max);
            TEST_FAIL("row_group_statistics", "id min/max range incorrect");
        }
        printf("  Row group %d: id range [%d, %d], %lld values\n",
               rg, min_val, max_val, (long long)stats.num_values);

        /* Float / double columns must also carry stats. */
        carquet_column_statistics_t fstats;
        if (carquet_reader_column_statistics(reader, rg, 1, &fstats) != CARQUET_OK ||
            !fstats.has_min_max) {
            TEST_FAIL("row_group_statistics", "value (DOUBLE) stats missing");
        }
        if (carquet_reader_column_statistics(reader, rg, 3, &fstats) != CARQUET_OK ||
            !fstats.has_min_max) {
            TEST_FAIL("row_group_statistics", "score (FLOAT) stats missing");
        }
    }

    carquet_reader_close(reader);
    TEST_PASS("row_group_statistics");
    return 0;
}

/* ============================================================================
 * Test: write_statistics=false suppresses stats
 * ============================================================================
 */

static int test_statistics_toggle(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "nostats");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) {
        TEST_FAIL("statistics_toggle", "schema create failed");
    }
    (void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = false;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("statistics_toggle", "writer create failed");
    }

    int32_t ids[100];
    for (int i = 0; i < 100; i++) ids[i] = i;
    (void)carquet_writer_write_batch(writer, 0, ids, 100, NULL, NULL);
    (void)carquet_writer_close(writer);
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) {
        remove(path);
        TEST_FAIL("statistics_toggle", "reader open failed");
    }

    carquet_column_statistics_t stats;
    carquet_status_t status = carquet_reader_column_statistics(reader, 0, 0, &stats);
    carquet_reader_close(reader);
    remove(path);

    if (status == CARQUET_OK && stats.has_min_max) {
        TEST_FAIL("statistics_toggle", "stats present with write_statistics=false");
    }

    TEST_PASS("statistics_toggle");
    return 0;
}

/* ============================================================================
 * Statistics functional test suite
 *
 * These tests exercise the per-column-physical-type stat code paths end-to-end
 * (writer → file → reader). They are designed to fail loudly when:
 *   - stats aren't emitted at all (the original bug)
 *   - min/max aggregation across pages / row-groups is wrong
 *   - byte-array lex ordering or truncation rules are violated
 *   - per-column overrides are not honored
 *   - null counts don't aggregate
 * ============================================================================
 */

/* Compare a stat's binary payload to an expected primitive value. */
#define ASSERT_PRIMITIVE_EQ(name, stats, field, type, expected) do { \
    if ((stats).field##_size != (int32_t)sizeof(type)) { \
        TEST_FAIL(name, "stat " #field "_size != sizeof(" #type ")"); \
    } \
    type _v; \
    memcpy(&_v, (stats).field, sizeof(_v)); \
    if (_v != (expected)) { \
        TEST_FAIL(name, "stat " #field " mismatch"); \
    } \
} while (0)

#define ASSERT_BYTES_EQ(name, ptr, size, expected_ptr, expected_size, msg) do { \
    if ((int32_t)(size) != (int32_t)(expected_size) || \
        memcmp((ptr), (expected_ptr), (size_t)(expected_size)) != 0) { \
        TEST_FAIL(name, msg); \
    } \
} while (0)

/* Numeric stats: positive AND negative ranges (signed compare must hold). */
static int test_stats_numeric_signed(void) {
    const char* tname = "stats_numeric_signed";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_numsigned");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "i32", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "i64", CARQUET_PHYSICAL_INT64, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "f32", CARQUET_PHYSICAL_FLOAT, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "f64", CARQUET_PHYSICAL_DOUBLE, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    int32_t i32[] = {  10, -7, 100, -2147483648, 0, 2147483647 };
    int64_t i64[] = { 1000, -1, 0, INT64_MIN, INT64_MAX, 42 };
    float   f32[] = { 1.5f, -3.25f, 0.0f, -0.0f, 100.0f, -100.5f };
    double  f64[] = { 1e10, -1e10, 0.0, 0.5, -0.5, 7.0 };
    const int N = 6;

    if (carquet_writer_write_batch(writer, 0, i32, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 1, i64, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 2, f32, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 3, f64, N, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write_batch");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "INT32 stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, int32_t, INT32_MIN);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, int32_t, INT32_MAX);

    if (carquet_reader_column_statistics(reader, 0, 1, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "INT64 stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, int64_t, INT64_MIN);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, int64_t, INT64_MAX);

    if (carquet_reader_column_statistics(reader, 0, 2, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "FLOAT stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, float, -100.5f);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, float, 100.0f);

    if (carquet_reader_column_statistics(reader, 0, 3, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "DOUBLE stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, double, -1e10);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, double, 1e10);

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

static int test_stats_unsigned_logical_order(void) {
    const char* tname = "stats_unsigned_logical_order";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_unsigned");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");

    carquet_logical_type_t u32_type = {
        .id = CARQUET_LOGICAL_INTEGER,
        .params.integer = { .bit_width = 32, .is_signed = false }
    };
    carquet_logical_type_t u64_type = {
        .id = CARQUET_LOGICAL_INTEGER,
        .params.integer = { .bit_width = 64, .is_signed = false }
    };

    (void)carquet_schema_add_column(schema, "u32", CARQUET_PHYSICAL_INT32, &u32_type,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "u64", CARQUET_PHYSICAL_INT64, &u64_type,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    uint32_t u32[] = { UINT32_MAX, 1u, 42u };
    uint64_t u64[] = { UINT64_MAX, 1ull, 42ull };
    if (carquet_writer_write_batch(writer, 0, u32, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 1, u64, 3, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write_batch");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "UINT32 stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, uint32_t, 1u);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, uint32_t, UINT32_MAX);

    if (carquet_reader_column_statistics(reader, 0, 1, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "UINT64 stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, uint64_t, 1ull);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, uint64_t, UINT64_MAX);

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* BOOLEAN stats: all-false / all-true / mixed across three row groups. */
static int test_stats_boolean(void) {
    const char* tname = "stats_boolean";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_bool");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BOOLEAN, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    uint8_t all_false[8] = {0,0,0,0,0,0,0,0};
    uint8_t all_true[8]  = {1,1,1,1,1,1,1,1};
    uint8_t mixed[8]     = {0,1,0,1,1,0,1,0};

    if (carquet_writer_write_batch(writer, 0, all_false, 8, NULL, NULL) != CARQUET_OK ||
        carquet_writer_new_row_group(writer) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 0, all_true, 8, NULL, NULL) != CARQUET_OK ||
        carquet_writer_new_row_group(writer) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 0, mixed, 8, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    /* all-false → min=0, max=0 */
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max || s.min_value_size != 1 || s.max_value_size != 1 ||
        ((const uint8_t*)s.min_value)[0] != 0 ||
        ((const uint8_t*)s.max_value)[0] != 0) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "all-false RG min/max wrong");
    }
    /* all-true → min=1, max=1 */
    if (carquet_reader_column_statistics(reader, 1, 0, &s) != CARQUET_OK ||
        !s.has_min_max ||
        ((const uint8_t*)s.min_value)[0] != 1 ||
        ((const uint8_t*)s.max_value)[0] != 1) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "all-true RG min/max wrong");
    }
    /* mixed → min=0, max=1 */
    if (carquet_reader_column_statistics(reader, 2, 0, &s) != CARQUET_OK ||
        !s.has_min_max ||
        ((const uint8_t*)s.min_value)[0] != 0 ||
        ((const uint8_t*)s.max_value)[0] != 1) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "mixed RG min/max wrong");
    }

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* BYTE_ARRAY: varying-length strings, lex min/max, null counts. */
static int test_stats_byte_array(void) {
    const char* tname = "stats_byte_array";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_ba");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "name", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    /* "a" is the lex min (shortest), "zz" is the lex max. Multiple lengths. */
    const char* strs[] = {"delta", "a", "echo", "zz", "bravo", "ch", "foxtrot"};
    const int N = (int)(sizeof(strs)/sizeof(strs[0]));
    carquet_byte_array_t arr[7];
    for (int i = 0; i < N; i++) {
        arr[i].data = (uint8_t*)(uintptr_t)strs[i];
        arr[i].length = (int32_t)strlen(strs[i]);
    }

    if (carquet_writer_write_batch(writer, 0, arr, N, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "stats missing");
    }
    ASSERT_BYTES_EQ(tname, s.min_value, s.min_value_size, "a", 1, "min != \"a\"");
    ASSERT_BYTES_EQ(tname, s.max_value, s.max_value_size, "zz", 2, "max != \"zz\"");

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* FIXED_LEN_BYTE_ARRAY with type_length=16 (e.g., UUID-ish). */
static int test_stats_flba(void) {
    const char* tname = "stats_flba";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_flba");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "uuid",
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 16, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    const int N = 5;
    uint8_t buf[5 * 16];
    memset(buf, 0, sizeof(buf));
    /* Hand-crafted to ensure a non-trivial lex ordering. */
    memcpy(buf + 0*16,  "AAAA000000000000", 16);  /* should be MIN */
    memcpy(buf + 1*16,  "BBBB000000000000", 16);
    memcpy(buf + 2*16,  "BAAA000000000000", 16);
    memcpy(buf + 3*16,  "ZZZZ999999999999", 16);  /* should be MAX */
    memcpy(buf + 4*16,  "MMMM000000000000", 16);

    if (carquet_writer_write_batch(writer, 0, buf, N, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "stats missing");
    }
    ASSERT_BYTES_EQ(tname, s.min_value, s.min_value_size,
        "AAAA000000000000", 16, "FLBA min wrong");
    ASSERT_BYTES_EQ(tname, s.max_value, s.max_value_size,
        "ZZZZ999999999999", 16, "FLBA max wrong");

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* BYTE_ARRAY > 32-byte truncation: must increment last byte to keep the
 * stored max >= the real max, and min must be exact when short enough. */
static int test_stats_byte_array_truncation(void) {
    const char* tname = "stats_byte_array_truncation";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_trunc");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "s", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    /* The lex max is 64 'a's ending in 'b' — truncation must increment the
     * 32nd byte (still 'a') so the stored bound is lex >= the real one. */
    char actual_max[64];
    memset(actual_max, 'a', 64);
    actual_max[63] = 'b';

    carquet_byte_array_t arr[2];
    arr[0].data = (uint8_t*)(uintptr_t)"aaa";
    arr[0].length = 3;
    arr[1].data = (uint8_t*)actual_max;
    arr[1].length = 64;

    if (carquet_writer_write_batch(writer, 0, arr, 2, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "stats missing");
    }
    /* min is short → must be exact. */
    ASSERT_BYTES_EQ(tname, s.min_value, s.min_value_size, "aaa", 3, "min not exact");

    /* max truncated to <=32 bytes, lex >= actual_max prefix. */
    if (s.max_value_size > 32 || s.max_value_size == 0) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "max length not truncated to 32");
    }
    /* The stored max must lex-compare >= the first max_value_size bytes of
     * the actual max. Compare against the original full max. */
    size_t cmp_n = (size_t)s.max_value_size < 64 ? (size_t)s.max_value_size : 64;
    int c = memcmp(s.max_value, actual_max, cmp_n);
    if (c < 0) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "truncated max is lex < actual max (invalid)");
    }
    /* The 32nd byte should have been incremented to 'b'. */
    if (s.max_value_size != 32 ||
        ((const uint8_t*)s.max_value)[31] != 'b') {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "expected last byte incremented from 'a' to 'b'");
    }

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* BYTE_ARRAY truncation with all-0xFF prefix: increment overflows so max
 * must be DROPPED to remain a valid upper bound. */
static int test_stats_byte_array_truncation_overflow(void) {
    const char* tname = "stats_byte_array_truncation_overflow";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_truncFF");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "s", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    /* 64 bytes all 0xFF: truncation increment cannot represent a valid
     * upper bound, so the stored max must be absent. */
    uint8_t max_FF[64];
    memset(max_FF, 0xFF, 64);

    carquet_byte_array_t arr[2];
    arr[0].data = (uint8_t*)(uintptr_t)"abc"; arr[0].length = 3;
    arr[1].data = max_FF; arr[1].length = 64;

    if (carquet_writer_write_batch(writer, 0, arr, 2, NULL, NULL) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "write");
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    /* The reader reports has_min_max only when BOTH bounds are present.
     * With max dropped due to overflow, we expect has_min_max=false. */
    carquet_status_t st = carquet_reader_column_statistics(reader, 0, 0, &s);
    if (st == CARQUET_OK && s.has_min_max && s.max_value_size > 0) {
        /* If a max IS stored, it must be lex >= 0xFF...0xFF, which means it
         * must be at least as long and not less anywhere. Verify it doesn't
         * lie. */
        size_t cmp_n = (size_t)s.max_value_size < 64 ? (size_t)s.max_value_size : 64;
        if (memcmp(s.max_value, max_FF, cmp_n) < 0) {
            carquet_reader_close(reader); carquet_test_cleanup(path);
            TEST_FAIL(tname, "max < 0xFF... (overflow truncation produced "
                             "invalid upper bound)");
        }
    }
    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* Stats must aggregate across multiple pages within a single column chunk. */
static int test_stats_multi_page_aggregation(void) {
    const char* tname = "stats_multi_page_aggregation";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_multipage");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "n", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    /* Tiny page size to force many pages within one column chunk. */
    opts.page_size = 1024;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    /* 100000 ints from -50000..49999 written in 10 batches. The smallest
     * value (-49997) appears in the middle batch; the largest (49999) at
     * the end. Without proper aggregation only the last page's stats
     * would survive. */
    const int total = 100000;
    int32_t* data = malloc((size_t)total * sizeof(int32_t));
    if (!data) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "alloc");
    }
    for (int i = 0; i < total; i++) data[i] = i - 50000;
    /* Force the extreme min into a middle batch (already there by linear
     * fill); add a sentinel guaranteed-low value to a middle position. */
    data[total / 2] = -123456;

    int written = 0;
    int batch_sz = total / 10;
    while (written < total) {
        int chunk = (total - written) < batch_sz ? (total - written) : batch_sz;
        if (carquet_writer_write_batch(writer, 0, data + written, chunk,
                                       NULL, NULL) != CARQUET_OK) {
            free(data); (void)carquet_writer_close(writer);
            carquet_schema_free(schema); carquet_test_cleanup(path);
            TEST_FAIL(tname, "write batch");
        }
        written += chunk;
    }
    free(data);

    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "writer close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    int32_t num_rgs = carquet_reader_num_row_groups(reader);
    /* For each row group find the row-group-local min/max from the row
     * group's index range, then verify reader returns matching stats. */
    int32_t global_min = INT32_MAX, global_max = INT32_MIN;
    for (int32_t rg = 0; rg < num_rgs; rg++) {
        carquet_column_statistics_t s;
        if (carquet_reader_column_statistics(reader, rg, 0, &s) != CARQUET_OK ||
            !s.has_min_max || s.min_value_size != 4 || s.max_value_size != 4) {
            carquet_reader_close(reader); carquet_test_cleanup(path);
            TEST_FAIL(tname, "row-group stats missing across pages");
        }
        int32_t mn, mx;
        memcpy(&mn, s.min_value, 4);
        memcpy(&mx, s.max_value, 4);
        if (mn > mx) {
            carquet_reader_close(reader); carquet_test_cleanup(path);
            TEST_FAIL(tname, "row-group min > max");
        }
        if (mn < global_min) global_min = mn;
        if (mx > global_max) global_max = mx;
    }
    /* Sentinel must show up in the global min; the highest value must be
     * 49999 (the largest data[i] before the sentinel-injection). */
    if (global_min != -123456) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "aggregated min did not include sentinel");
    }
    if (global_max != 49999) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "aggregated max wrong");
    }

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* Null counts must be aggregated across pages and reflected in stats. */
static int test_stats_null_counts(void) {
    const char* tname = "stats_null_counts";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_nulls");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "opt", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    /* 100 logical rows: every 3rd value is NULL. Non-null values are
     * sparse-packed at the front per the writer contract. */
    const int N = 100;
    int16_t def_levels[100];
    int32_t values[100];
    int n_non_null = 0;
    int expected_nulls = 0;
    for (int i = 0; i < N; i++) {
        if (i % 3 == 0) {
            def_levels[i] = 0;
            expected_nulls++;
        } else {
            def_levels[i] = 1;
            values[n_non_null++] = i;
        }
    }

    if (carquet_writer_write_batch(writer, 0, values, N, def_levels, NULL) != CARQUET_OK ||
        carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "write/close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "stats fetch failed");
    }
    if (!s.has_null_count || s.null_count != expected_nulls) {
        printf("    null_count=%lld expected=%d\n",
               (long long)s.null_count, expected_nulls);
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "null_count wrong");
    }
    /* min should be the smallest non-null value = 1 (i % 3 != 0 starts at 1). */
    if (!s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "min/max missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, int32_t, 1);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, int32_t, 98);

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* Per-column stats override: global ON, but column 1 disabled. */
static int test_stats_per_column_override(void) {
    const char* tname = "stats_per_column_override";
    carquet_error_t err = CARQUET_ERROR_INIT;
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "stats_override");

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL(tname, "schema create");
    (void)carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) { carquet_schema_free(schema); TEST_FAIL(tname, "writer create"); }

    if (carquet_writer_set_column_statistics(writer, 1, false) != CARQUET_OK) {
        (void)carquet_writer_close(writer); carquet_schema_free(schema);
        carquet_test_cleanup(path); TEST_FAIL(tname, "set_column_statistics failed");
    }

    int32_t a[] = {3, 1, 4, 1, 5, 9, 2, 6};
    int32_t b[] = {7, 0, 8, 1, 8, 2, 4, 5};
    if (carquet_writer_write_batch(writer, 0, a, 8, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(writer, 1, b, 8, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(writer) != CARQUET_OK) {
        carquet_schema_free(schema); carquet_test_cleanup(path);
        TEST_FAIL(tname, "write/close");
    }
    carquet_schema_free(schema);

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { carquet_test_cleanup(path); TEST_FAIL(tname, "reader open"); }

    /* Column 0: stats present. */
    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, 0, 0, &s) != CARQUET_OK ||
        !s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "col 0 stats missing");
    }
    ASSERT_PRIMITIVE_EQ(tname, s, min_value, int32_t, 1);
    ASSERT_PRIMITIVE_EQ(tname, s, max_value, int32_t, 9);

    /* Column 1: override disabled, stats must be absent. */
    carquet_status_t st = carquet_reader_column_statistics(reader, 0, 1, &s);
    if (st == CARQUET_OK && s.has_min_max) {
        carquet_reader_close(reader); carquet_test_cleanup(path);
        TEST_FAIL(tname, "col 1 stats present despite override=false");
    }

    carquet_reader_close(reader);
    carquet_test_cleanup(path);
    TEST_PASS(tname);
    return 0;
}

/* ============================================================================
 * Test: Predicate Pushdown (Row Group Filtering)
 * ============================================================================
 */

static int test_predicate_pushdown(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, NULL, &err);
    if (!reader) {
        TEST_FAIL("predicate_pushdown", "failed to open file");
    }

    int32_t num_row_groups = carquet_reader_num_row_groups(reader);

    /* Test: Find row groups where id > 5000 */
    int32_t search_value = 5000;
    int32_t matching[100];

    int32_t num_matching = carquet_reader_filter_row_groups(
        reader,
        0,  /* column 0 = id */
        CARQUET_COMPARE_GT,
        &search_value,
        sizeof(int32_t),
        matching,
        100);

    printf("  Row groups with id > %d: %d (of %d total)\n",
           search_value, num_matching, num_row_groups);

    /* Should filter out roughly half the row groups */
    assert(num_matching > 0);
    assert(num_matching <= num_row_groups);

    /* Test: Find row groups where id == 100 (should match only 1 or few) */
    search_value = 100;
    num_matching = carquet_reader_filter_row_groups(
        reader,
        0,
        CARQUET_COMPARE_EQ,
        &search_value,
        sizeof(int32_t),
        matching,
        100);

    printf("  Row groups that might contain id == %d: %d\n", search_value, num_matching);
    assert(num_matching >= 1);  /* At least one should match */

    /* Test: Find row groups where id < 0 (should match none) */
    search_value = 0;
    num_matching = carquet_reader_filter_row_groups(
        reader,
        0,
        CARQUET_COMPARE_LT,
        &search_value,
        sizeof(int32_t),
        matching,
        100);

    printf("  Row groups with id < 0: %d (should be 0)\n", num_matching);
    /* All IDs start from 0, so no row group should have id < 0 */

    carquet_reader_close(reader);
    TEST_PASS("predicate_pushdown");
    return 0;
}

/* ============================================================================
 * Test: Buffer-based Reading (simulates mmap)
 * ============================================================================
 */

static int test_buffer_reading(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Read entire file into memory */
    FILE* f = fopen(TEST_FILE, "rb");
    if (!f) {
        TEST_FAIL("buffer_reading", "failed to open file");
    }

    carquet_fseek64(f, 0, SEEK_END);
    size_t size = (size_t)carquet_ftell64(f);
    carquet_fseek64(f, 0, SEEK_SET);

    uint8_t* buffer = malloc(size);
    if (!buffer) {
        fclose(f);
        TEST_FAIL("buffer_reading", "failed to allocate buffer");
    }

    size_t read = fread(buffer, 1, size, f);
    fclose(f);

    if (read != size) {
        free(buffer);
        TEST_FAIL("buffer_reading", "failed to read file");
    }

    printf("  File size: %zu bytes\n", size);

    /* Open from buffer */
    carquet_reader_t* reader = carquet_reader_open_buffer(buffer, size, NULL, &err);
    if (!reader) {
        free(buffer);
        printf("  Failed to open from buffer: %s\n", err.message);
        TEST_FAIL("buffer_reading", "failed to open from buffer");
    }

    /* Verify we can read metadata */
    int64_t num_rows = carquet_reader_num_rows(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);

    printf("  Buffer read: %lld rows, %d columns\n", (long long)num_rows, num_cols);
    assert(num_rows == NUM_ROWS);
    assert(num_cols == 4);

    carquet_reader_close(reader);
    free(buffer);

    TEST_PASS("buffer_reading");
    return 0;
}

/* ============================================================================
 * Test: Full Pipeline (projection + filtering)
 * ============================================================================
 */

static int test_full_pipeline(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open(TEST_FILE, NULL, &err);
    if (!reader) {
        TEST_FAIL("full_pipeline", "failed to open file");
    }

    /* Step 1: Filter row groups where category column might have value 5 */
    int32_t category_value = 5;
    int32_t matching_rgs[100];
    int32_t num_matching = carquet_reader_filter_row_groups(
        reader,
        2,  /* category column */
        CARQUET_COMPARE_EQ,
        &category_value,
        sizeof(int32_t),
        matching_rgs,
        100);

    printf("  Row groups that might contain category=5: %d\n", num_matching);

    /* Step 2: Read only id and category columns from matching row groups */
    int32_t proj_cols[] = {0, 2};
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.column_indices = proj_cols;
    config.num_columns = 2;
    config.batch_size = 1000;

    int64_t total_matching_rows = 0;

    /* In a real implementation, we'd only read the matching row groups */
    /* For now, we read all and count matches */
    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (batch_reader) {
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
            const void* cat_data;
            const uint8_t* null_bitmap;
            int64_t num_values;
            (void)carquet_row_batch_column(batch, 1, &cat_data, &null_bitmap, &num_values);
            (void)null_bitmap;

            const int32_t* categories = (const int32_t*)cat_data;
            for (int64_t i = 0; i < num_values; i++) {
                if (categories[i] == category_value) {
                    total_matching_rows++;
                }
            }

            carquet_row_batch_free(batch);
            batch = NULL;
        }
        carquet_batch_reader_free(batch_reader);
    }

    printf("  Rows with category=5: %lld (expected ~%d)\n",
           (long long)total_matching_rows, NUM_ROWS / 10);

    /* Should be approximately 10% of rows (category 0-9) */
    assert(total_matching_rows > 0);
    assert(total_matching_rows <= NUM_ROWS);

    carquet_reader_close(reader);
    TEST_PASS("full_pipeline");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    int failures = 0;

    /* Initialize portable temp file path */
    carquet_test_temp_path(TEST_FILE, sizeof(TEST_FILE), "production");

    printf("=== Production Feature Tests ===\n\n");

    /* Create test file */
    printf("Creating test file with %d rows in ~%d row groups...\n", NUM_ROWS, NUM_ROW_GROUPS);
    if (create_test_file() != 0) {
        printf("FATAL: Failed to create test file\n");
        return 1;
    }
    printf("Test file created: %s\n\n", TEST_FILE);

    /* Run tests */
    failures += test_column_projection();
    failures += test_row_group_statistics();
    failures += test_statistics_toggle();
    failures += test_stats_numeric_signed();
    failures += test_stats_unsigned_logical_order();
    failures += test_stats_boolean();
    failures += test_stats_byte_array();
    failures += test_stats_flba();
    failures += test_stats_byte_array_truncation();
    failures += test_stats_byte_array_truncation_overflow();
    failures += test_stats_multi_page_aggregation();
    failures += test_stats_null_counts();
    failures += test_stats_per_column_override();
    failures += test_predicate_pushdown();
    failures += test_buffer_reading();
    failures += test_full_pipeline();

    /* Cleanup */
    remove(TEST_FILE);

    printf("\n");
    if (failures == 0) {
        printf("All production tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
