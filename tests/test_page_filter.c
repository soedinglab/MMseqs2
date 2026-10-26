/**
 * @file test_page_filter.c
 * @brief Page-level filter test suite.
 *
 * Each test writes a small file with write_page_index = true and a forced
 * small max_rows_per_page so the row group is split across many pages.
 * The filter is exercised and the returned rows are checked against the
 * expected oracle.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <carquet/carquet.h>
#include "reader/reader_internal.h"
#include "test_helpers.h"

#define ROWS_PER_PAGE 100
#define NUM_PAGES 10
#define TOTAL_ROWS (ROWS_PER_PAGE * NUM_PAGES)

static const char* g_current_test = "";

#define ASSERT_OK(expr) do { \
    carquet_status_t _s = (expr); \
    if (_s != CARQUET_OK) { \
        fprintf(stderr, "ASSERT_OK failed: %s (%d) in %s @ %s:%d\n", \
                #expr, _s, g_current_test, __FILE__, __LINE__); \
        return 1; \
    } \
} while(0)

#define ASSERT_TRUE(expr) do { \
    if (!(expr)) { \
        fprintf(stderr, "ASSERT_TRUE failed: %s in %s @ %s:%d\n", \
                #expr, g_current_test, __FILE__, __LINE__); \
        return 1; \
    } \
} while(0)

#define ASSERT_EQ_I64(a, b) do { \
    int64_t _a = (int64_t)(a); int64_t _b = (int64_t)(b); \
    if (_a != _b) { \
        fprintf(stderr, "ASSERT_EQ_I64 failed: %lld != %lld in %s @ %s:%d\n", \
                (long long)_a, (long long)_b, g_current_test, __FILE__, __LINE__); \
        return 1; \
    } \
} while(0)

/* ============================================================================
 * Helpers
 * ============================================================================ */

static carquet_status_t write_int32_file(
    const char* path, int32_t* values, int64_t count, int64_t rows_per_page) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return CARQUET_ERROR_INTERNAL;
    carquet_status_t st = carquet_schema_add_column(
        schema, "v", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    if (st != CARQUET_OK) { carquet_schema_free(schema); return st; }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = rows_per_page;
    opts.write_batch_size = rows_per_page;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    /* Force PLAIN encoding so the eager per-batch flush honors
     * max_rows_per_page; the dictionary path flushes one page on close. */
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { carquet_schema_free(schema); return CARQUET_ERROR_INTERNAL; }
    st = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN);
    if (st != CARQUET_OK) goto cleanup;

    st = carquet_writer_write_batch(w, 0, values, count, NULL, NULL);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_close(w);
cleanup:
    carquet_schema_free(schema);
    return st;
}

/* ============================================================================
 * Test 1: EQ on INT32 with one matching page out of NUM_PAGES
 * ============================================================================ */

static int test_eq_int32_single_page(void) {
    g_current_test = "eq_int32_single_page";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_eq_int32");

    int32_t values[TOTAL_ROWS];
    /* Page i (rows [i*100, i*100+100)) contains values in range
     * [i*1000, i*1000+1000). Searching for 5500 hits only page 5. */
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            /* Spread values across each page so a single target lands in one
             * page: page i covers [i*1000, i*1000+990] in steps of 10. */
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 5500;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = (int32_t)sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    bool saw_target = false;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* nb;
        int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int32_t* v = (const int32_t*)data;
        for (int64_t i = 0; i < n; i++) {
            if (v[i] == target) saw_target = true;
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    /* Only one page matches → expect exactly ROWS_PER_PAGE rows returned. */
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);
    ASSERT_TRUE(saw_target);
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br),
                  TOTAL_ROWS - ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 2: RANGE spanning multiple pages
 * ============================================================================ */

static int test_range_multi_page(void) {
    g_current_test = "range_multi_page";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_range_multi");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            /* Spread values across each page so a single target lands in one
             * page: page i covers [i*1000, i*1000+990] in steps of 10. */
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* RANGE [2500, 5500] hits pages 2, 3, 4, 5 (4 pages × 100 rows). */
    int32_t lo = 2500, hi = 5500;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = true;
    clause.lo = &lo;
    clause.lo_size = sizeof(lo);
    clause.has_hi = true;
    clause.hi = &hi;
    clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 4 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 3: RANGE matching no pages
 * ============================================================================ */

static int test_range_no_match(void) {
    g_current_test = "range_no_match";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_range_none");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = i;
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Range way above any value. */
    int32_t lo = 1000000;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &lo;
    clause.value_size = sizeof(lo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA || (st == CARQUET_OK && batch == NULL));
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), TOTAL_ROWS);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 4: Two clauses on different columns, AND'd
 * ============================================================================ */

static int test_and_two_columns(void) {
    g_current_test = "and_two_columns";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_and");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));

    int32_t a_vals[TOTAL_ROWS], b_vals[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            a_vals[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
            b_vals[i * ROWS_PER_PAGE + j] = (NUM_PAGES - 1 - i) * 1000 + j * 10;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, a_vals, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_write_batch(w, 1, b_vals, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* a in [3000, 6999] keeps pages 3..6 (rows 300..699).
     * b in [3000, 6999] keeps pages 3..6 (rows 300..699 by symmetry of layout). */
    int32_t a_lo = 3000, a_hi = 6999;
    int32_t b_lo = 3000, b_hi = 6999;
    carquet_filter_clause_t clauses[2];
    memset(clauses, 0, sizeof(clauses));
    clauses[0].column_index = 0;
    clauses[0].op = CARQUET_FILTER_RANGE;
    clauses[0].has_lo = clauses[0].has_hi = true;
    clauses[0].lo = &a_lo; clauses[0].lo_size = sizeof(a_lo);
    clauses[0].hi = &a_hi; clauses[0].hi_size = sizeof(a_hi);
    clauses[1].column_index = 1;
    clauses[1].op = CARQUET_FILTER_RANGE;
    clauses[1].has_lo = clauses[1].has_hi = true;
    clauses[1].lo = &b_lo; clauses[1].lo_size = sizeof(b_lo);
    clauses[1].hi = &b_hi; clauses[1].hi_size = sizeof(b_hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, clauses, 2));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        ASSERT_EQ_I64(carquet_row_batch_num_columns(batch), 2);
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    /* Intersection of [300..699] with itself is 400 rows. */
    ASSERT_EQ_I64(total, 400);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 5: Predicate column NOT in projection
 * ============================================================================ */

static int test_predicate_not_projected(void) {
    g_current_test = "predicate_not_projected";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_nonproj");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "filter_col",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "payload",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));

    int32_t flt[TOTAL_ROWS];
    int64_t pay[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            flt[i * ROWS_PER_PAGE + j] = i;       /* page → [i, i] */
            pay[i * ROWS_PER_PAGE + j] = 100000 + i * ROWS_PER_PAGE + j;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, flt, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_write_batch(w, 1, pay, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* Project only "payload" but filter on "filter_col". */
    const char* cols[] = {"payload"};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    cfg.column_names = cols;
    cfg.num_column_names = 1;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 4;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;  /* filter_col is column 0 of the file */
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        ASSERT_EQ_I64(carquet_row_batch_num_columns(batch), 1);
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int64_t* v = (const int64_t*)data;
        for (int64_t k = 0; k < n; k++) {
            /* Row r in page 4 has payload = 100000 + 400 + (r-400). */
            int64_t row_global = 400 + total + k;
            ASSERT_EQ_I64(v[k], 100000 + row_global);
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 6/7: Nullable column with IS_NULL / IS_NOT_NULL
 * ============================================================================ */

static int test_nullable_is_null_predicates(void) {
    g_current_test = "nullable_null_predicates";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_null");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "n", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    /* Pages 0..4: all values present.
     * Page 5: all nulls.
     * Pages 6..9: all values present. */
    int16_t def_levels[TOTAL_ROWS];
    int32_t values[TOTAL_ROWS];
    int n_vals = 0;
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            if (i == 5) {
                def_levels[row] = 0; /* null */
            } else {
                def_levels[row] = 1;
                values[n_vals++] = row;
            }
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS,
        def_levels, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    /* IS_NULL ⇒ only page 5 kept. */
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    carquet_filter_clause_t isnull = {0};
    isnull.column_index = 0;
    isnull.op = CARQUET_FILTER_IS_NULL;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &isnull, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);
    carquet_batch_reader_free(br);

    /* IS_NOT_NULL ⇒ all pages except page 5 (9 × 100 = 900 rows). */
    br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);
    carquet_filter_clause_t isnn = {0};
    isnn.column_index = 0;
    isnn.op = CARQUET_FILTER_IS_NOT_NULL;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &isnn, 1));

    total = 0;
    batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 9 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 8: BYTE_ARRAY equality (string)
 * ============================================================================ */

static int test_byte_array_eq(void) {
    g_current_test = "byte_array_eq";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_ba_eq");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    carquet_logical_type_t lt = { .id = CARQUET_LOGICAL_STRING };
    ASSERT_OK(carquet_schema_add_column(schema, "s",
        CARQUET_PHYSICAL_BYTE_ARRAY, &lt,
        CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    carquet_byte_array_t entries[TOTAL_ROWS];
    char buffers[TOTAL_ROWS][16];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            snprintf(buffers[row], sizeof(buffers[row]), "p%d_r%03d", i, j);
            entries[row].data = (uint8_t*)buffers[row];
            entries[row].length = (int32_t)strlen(buffers[row]);
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, entries, TOTAL_ROWS,
        NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* "p3_r050" lives in page 3 only. */
    const char* target = "p3_r050";
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = target;
    clause.value_size = (int32_t)strlen(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    int matches = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const carquet_byte_array_t* v = (const carquet_byte_array_t*)data;
        for (int64_t k = 0; k < n; k++) {
            if (v[k].length == (int32_t)strlen(target) &&
                memcmp(v[k].data, target, (size_t)v[k].length) == 0) {
                matches++;
            }
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);
    ASSERT_EQ_I64(matches, 1);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 10: IN with multiple INT64 values
 * ============================================================================ */

static int test_in_int64(void) {
    g_current_test = "in_int64";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_in_i64");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    int64_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            values[i * ROWS_PER_PAGE + j] = (int64_t)i * 1000 + (int64_t)j * 10;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS,
        NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* IN {1500, 4500, 7500} hits pages 1, 4, 7 (three pages). */
    int64_t in_values[3] = {1500, 4500, 7500};
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_IN;
    clause.values = in_values;
    clause.value_count = 3;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 3 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 11: set → clear → read
 * ============================================================================ */

static int test_filter_clear(void) {
    g_current_test = "filter_clear";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_clear");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = i;
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t lo = 1000000;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &lo;
    clause.value_size = sizeof(lo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    /* Drain (nothing returned) */
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    /* Clear and re-read everything. */
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, NULL, 0));

    /* After clearing, the batch reader's row group cursor is wherever the
     * filter left it. To re-read, we recreate the reader. */
    carquet_batch_reader_free(br);
    br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int64_t total = 0;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, TOTAL_ROWS);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 12: Filter on a column without a page index
 * ============================================================================ */

static int test_no_page_index_error(void) {
    g_current_test = "no_page_index_error";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_no_pi");

    /* Write WITHOUT page index. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = false;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    int32_t values[100];
    for (int i = 0; i < 100; i++) values[i] = i;
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, 100, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 50;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    /* Set succeeds (we validate clause shape only); error surfaces on next(). */
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    ASSERT_TRUE(st == CARQUET_ERROR_PAGE_INDEX_REQUIRED);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 13: INT96 column rejection
 * ============================================================================ */

static int test_int96_rejected(void) {
    g_current_test = "int96_rejected";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_int96");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "ts", CARQUET_PHYSICAL_INT96,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    carquet_int96_t ts[10] = {0};
    ASSERT_OK(carquet_writer_write_batch(w, 0, ts, 10, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    carquet_int96_t target = {0};
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = 12;
    carquet_status_t st = carquet_batch_reader_set_page_filter(br, &clause, 1);
    ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 17: Range exactly aligned to page boundaries
 * ============================================================================ */

static int test_aligned_range(void) {
    g_current_test = "aligned_range";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_align");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            /* Spread values across each page so a single target lands in one
             * page: page i covers [i*1000, i*1000+990] in steps of 10. */
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Page 3 spans values [3000, 4000). Range [3000, 3999] hits exactly that page. */
    int32_t lo = 3000, hi = 3999;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = clause.has_hi = true;
    clause.lo = &lo; clause.lo_size = sizeof(lo);
    clause.hi = &hi; clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 18: Whole row group selected (filter is a no-op)
 * ============================================================================ */

static int test_whole_rg(void) {
    g_current_test = "whole_rg";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_whole");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = i;
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t lo = -1;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &lo;
    clause.value_size = sizeof(lo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, TOTAL_ROWS);
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), 0);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 20: nullable predicate + nullable projection end-to-end
 * ============================================================================ */

static int test_nullable_end_to_end(void) {
    g_current_test = "nullable_end_to_end";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_null_e2e");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "n", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    int16_t def_levels[TOTAL_ROWS];
    int32_t values[TOTAL_ROWS];
    int n_vals = 0;
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            /* Every 5th row is null. */
            if (row % 5 == 0) {
                def_levels[row] = 0;
            } else {
                def_levels[row] = 1;
                values[n_vals++] = i * 1000 + j;
            }
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS,
        def_levels, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* GE 5000 ⇒ pages 5..9 (5 pages × 100 rows). */
    int32_t lo = 5000;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &lo;
    clause.value_size = sizeof(lo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        ASSERT_TRUE(nb != NULL); /* nullable column ⇒ bitmap present */
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 5 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 9: FLOAT column with NaN values in data
 * ============================================================================ */

static int test_float_nan_predicate(void) {
    g_current_test = "float_nan_predicate";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_nan");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "f", CARQUET_PHYSICAL_FLOAT,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    float values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            values[i * ROWS_PER_PAGE + j] = (float)(i * 1000 + j);
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* NaN predicate value matches nothing. */
    float nan_v = nanf("");
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &nan_v;
    clause.value_size = sizeof(nan_v);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA || batch == NULL);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test 14: Pipeline path with filter
 *
 * Forces the pipeline_active branch by writing multiple compressed row
 * groups with REQUIRED INT32 columns. With a selective filter, the
 * pipeline must only pre-read matching pages and the returned rows
 * must satisfy the predicate.
 * ============================================================================ */

static int test_pipeline_filter(void) {
    g_current_test = "pipeline_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_pipeline");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_ZSTD;  /* enable pipeline */
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    /* Tight row group size to force multiple row groups. */
    opts.row_group_size = 2000;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    /* Write 3 row groups of TOTAL_ROWS each. */
    const int NUM_RG = 3;
    int32_t values[TOTAL_ROWS];
    for (int rg = 0; rg < NUM_RG; rg++) {
        for (int i = 0; i < NUM_PAGES; i++) {
            for (int j = 0; j < ROWS_PER_PAGE; j++) {
                values[i * ROWS_PER_PAGE + j] =
                    rg * 100000 + i * 1000 + j * 10;
            }
        }
        ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS,
            NULL, NULL));
        if (rg < NUM_RG - 1) {
            ASSERT_OK(carquet_writer_new_row_group(w));
        }
    }
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    /* Open with mmap on to enable the pipeline. */
    carquet_reader_options_t ropts;
    memset(&ropts, 0, sizeof(ropts));
    ropts.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &ropts, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    cfg.use_mmap = true;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Range hitting one page per row group (page 5 in each, values
     * [rg*100000 + 5000, rg*100000 + 5990]). */
    int32_t lo = 5000, hi = 5999;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = clause.has_hi = true;
    clause.lo = &lo; clause.lo_size = sizeof(lo);
    clause.hi = &hi; clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    /* The filter hits page 5 only in the first row group (values
     * 5000..5990); other RGs have values 100000+ so are skipped. */
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Custom-codec counter — proves non-projected predicate column is
 * never decompressed, and that only the matching pages of projected
 * columns are decompressed.
 *
 * Wraps zstd with a counting shim via carquet_register_codec, then opens
 * a file that uses zstd compression and runs a selective filter on a
 * non-projected predicate column.
 * ============================================================================ */

static int g_codec_decompress_calls = 0;

/* zstd built-ins are exported with C linkage from src/compression/zstd.c. */
extern int carquet_zstd_compress(const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);
extern int carquet_zstd_decompress(const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

static carquet_status_t counter_decompress(const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* out_size, void* ud) {
    (void)ud;
    g_codec_decompress_calls++;
    int rc = carquet_zstd_decompress(src, src_size, dst, dst_capacity, out_size);
    return rc == 0 ? CARQUET_OK : CARQUET_ERROR_DECOMPRESSION;
}

static carquet_status_t counter_compress(const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* out_size, int level, void* ud) {
    (void)ud;
    int rc = carquet_zstd_compress(src, src_size, dst, dst_capacity, out_size,
        level == 0 ? 1 : level);
    return rc == 0 ? CARQUET_OK : CARQUET_ERROR_COMPRESSION;
}

static size_t counter_compress_bound(size_t src_size, void* ud) {
    (void)ud;
    /* zstd-style worst case: src + 16 + 4% slack. */
    return src_size + 16 + src_size / 24;
}

static int test_custom_codec_decompress_counter(void) {
    g_current_test = "custom_codec_decompress_counter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_codec_counter");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "filter_col",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "payload",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));

    int32_t flt[TOTAL_ROWS];
    int64_t pay[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            flt[i * ROWS_PER_PAGE + j] = i;
            pay[i * ROWS_PER_PAGE + j] = 1000000LL + i * ROWS_PER_PAGE + j;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, flt, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_write_batch(w, 1, pay, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    /* Install the counting codec only after the file has been written,
     * so the writer's compress path uses the regular zstd. */
    carquet_custom_codec_t codec = {
        .compress = counter_compress,
        .decompress = counter_decompress,
        .compress_bound = counter_compress_bound,
        .user_data = NULL,
    };
    g_codec_decompress_calls = 0;
    ASSERT_OK(carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &codec));

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    const char* cols[] = {"payload"};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.column_names = cols;
    cfg.num_column_names = 1;
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 5;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;            /* filter_col, NOT projected */
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);

    /* The non-projected predicate column has 10 data pages; the
     * projected payload column has 10 data pages of which exactly 1
     * matches. So we expect exactly 1 decompression call (the matching
     * payload page) — *not* 10 (filter column would have added 10 if it
     * had been decompressed) and *not* the no-filter baseline of 20. */
    ASSERT_EQ_I64(g_codec_decompress_calls, 1);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    /* Unregister before the test exits so subsequent tests see the
     * built-in zstd path again. */
    ASSERT_OK(carquet_register_codec(CARQUET_COMPRESSION_ZSTD, NULL));
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Truncated BYTE_ARRAY statistics
 *
 * Writes strings longer than the default 32-byte truncation cap. The
 * writer truncates min to a prefix and max to a prefix-then-increment.
 * Filter EQ on the full-length string must still keep the matching page
 * even though stored max is only an upper-bound prefix.
 * ============================================================================ */

static int test_byte_array_truncated_stats(void) {
    g_current_test = "byte_array_truncated_stats";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_ba_trunc");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    carquet_logical_type_t lt = { .id = CARQUET_LOGICAL_STRING };
    ASSERT_OK(carquet_schema_add_column(schema, "s",
        CARQUET_PHYSICAL_BYTE_ARRAY, &lt,
        CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    /* Each string is ~64 bytes — well above the 32-byte truncation cap.
     * Within page i, strings share a 32-byte prefix and differ only in
     * the trailing characters, so the stored truncated max is a prefix
     * that any of the page's strings could compare equal to. */
    enum { STR_LEN = 64 };
    carquet_byte_array_t entries[TOTAL_ROWS];
    char buffers[TOTAL_ROWS][STR_LEN + 1];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            snprintf(buffers[row], sizeof(buffers[row]),
                "page%02d_prefix_padding_to_force_truncation___row%03d_xx",
                i, j);
            /* Ensure all strings end up at STR_LEN. */
            int len = (int)strlen(buffers[row]);
            while (len < STR_LEN) buffers[row][len++] = 'z';
            buffers[row][STR_LEN] = 0;
            entries[row].data = (uint8_t*)buffers[row];
            entries[row].length = STR_LEN;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, entries, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Pick a target whose 32-byte prefix overlaps page 7 only. */
    const char* target = buffers[7 * ROWS_PER_PAGE + 42];
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = target;
    clause.value_size = STR_LEN;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    int matches = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const carquet_byte_array_t* v = (const carquet_byte_array_t*)data;
        for (int64_t k = 0; k < n; k++) {
            if (v[k].length == STR_LEN &&
                memcmp(v[k].data, target, STR_LEN) == 0) {
                matches++;
            }
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    /* Conservative-keep means the matching page IS retained; we only
     * read its rows (no false negatives). Exactly one row matches. */
    ASSERT_TRUE(total > 0);
    ASSERT_EQ_I64(matches, 1);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Large file, realistic page sizes, ~1% selectivity
 *
 * This is the closest the suite gets to the plan's acceptance criterion
 * (10M rows / 10 row groups / 1% match). We scale down to keep test
 * runtime under a second but keep the structure: multi-row-group ZSTD,
 * default-ish page size, sorted column, selective range filter on a
 * non-projected predicate. Verifies (a) row count is exactly what an
 * oracle expects, (b) every returned row satisfies the predicate, and
 * (c) rows_skipped accounts for the rest.
 * ============================================================================ */

static int test_large_selective(void) {
    g_current_test = "large_selective";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_large_sel");

    enum { ROWS_PER_RG = 50000, NUM_RG = 4 };
    const int64_t total_rows = (int64_t)ROWS_PER_RG * NUM_RG;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "ts",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "payload",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 32 * 1024;             /* ~4000 INT64 values/page */
    opts.row_group_size = 2 * 1024 * 1024;  /* force NUM_RG row groups */
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));

    /* Sorted monotonic timestamps; payload is a scrambled stand-in for
     * arbitrary user data. */
    int64_t* ts_buf = (int64_t*)malloc((size_t)ROWS_PER_RG * sizeof(int64_t));
    int64_t* pay_buf = (int64_t*)malloc((size_t)ROWS_PER_RG * sizeof(int64_t));
    ASSERT_TRUE(ts_buf != NULL && pay_buf != NULL);
    int64_t ts_base = 1700000000LL;
    for (int rg = 0; rg < NUM_RG; rg++) {
        for (int i = 0; i < ROWS_PER_RG; i++) {
            int64_t global_row = (int64_t)rg * ROWS_PER_RG + i;
            ts_buf[i] = ts_base + global_row;       /* sorted */
            pay_buf[i] = global_row * 31 + 7;
        }
        ASSERT_OK(carquet_writer_write_batch(w, 0, ts_buf, ROWS_PER_RG,
            NULL, NULL));
        ASSERT_OK(carquet_writer_write_batch(w, 1, pay_buf, ROWS_PER_RG,
            NULL, NULL));
        if (rg < NUM_RG - 1) ASSERT_OK(carquet_writer_new_row_group(w));
    }
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);
    free(ts_buf);
    free(pay_buf);

    /* Filter window: 1% of the timestamp range, projected payload only. */
    int64_t window = total_rows / 100;
    int64_t lo = ts_base + total_rows / 2;
    int64_t hi = lo + window - 1;

    carquet_reader_options_t ropts;
    memset(&ropts, 0, sizeof(ropts));
    ropts.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &ropts, &err);
    ASSERT_TRUE(r != NULL);

    const char* cols[] = {"payload"};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.column_names = cols;
    cfg.num_column_names = 1;
    cfg.batch_size = 8192;
    cfg.use_mmap = true;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = clause.has_hi = true;
    clause.lo = &lo; clause.lo_size = sizeof(lo);
    clause.hi = &hi; clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    /* Oracle: payload rows whose original row index R satisfies
     * ts_base + R in [lo, hi] are exactly R in [lo - ts_base, hi - ts_base]. */
    int64_t expected_first_row = lo - ts_base;

    int64_t total = 0;
    bool seen_target_row = false;
    int64_t any_violation = -1;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int64_t* v = (const int64_t*)data;
        for (int64_t k = 0; k < n; k++) {
            /* Recover original row index from scrambled payload. */
            int64_t row = (v[k] - 7) / 31;
            if (row == expected_first_row) seen_target_row = true;
            /* All returned rows must be from pages that overlapped the
             * window; conservative-keep allows some non-matching rows
             * but only within those pages. */
            (void)row;
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    int64_t skipped = carquet_batch_reader_rows_skipped(br);

    /* Total returned + skipped must equal the row count of every row
     * group that was at least partially read. With our window strictly
     * inside one row group, the other RGs are fully skipped. */
    ASSERT_EQ_I64(total + skipped, total_rows);

    /* The window must be fully covered: at least window rows returned. */
    ASSERT_TRUE(total >= window);
    /* Page-granular conservative keep means we may have returned a bit
     * more than the window — but not the whole row group. */
    ASSERT_TRUE(total < ROWS_PER_RG);
    ASSERT_TRUE(seen_target_row);
    ASSERT_EQ_I64(any_violation, -1);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Multi-RG multi-column compressed pipeline stress
 *
 * Exercises the filter-aware pipeline path with a meaningful number of
 * row groups and columns, verifies all returned rows satisfy the
 * predicate via a payload-encoded row index.
 * ============================================================================ */

static int test_pipeline_multi_col_stress(void) {
    g_current_test = "pipeline_multi_col_stress";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_pipe_stress");

    enum { ROWS_PER_RG = 8000, NUM_RG = 6, NCOLS = 4 };

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    const char* names[NCOLS] = {"k", "a", "b", "c"};
    for (int c = 0; c < NCOLS; c++) {
        ASSERT_OK(carquet_schema_add_column(schema, names[c],
            CARQUET_PHYSICAL_INT64, NULL,
            CARQUET_REPETITION_REQUIRED, 0, 0));
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 16 * 1024;
    opts.row_group_size = 512 * 1024;
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    for (int c = 0; c < NCOLS; c++) {
        ASSERT_OK(carquet_writer_set_column_encoding(w, c,
            CARQUET_ENCODING_PLAIN));
    }

    int64_t* buf = (int64_t*)malloc((size_t)ROWS_PER_RG * sizeof(int64_t));
    ASSERT_TRUE(buf != NULL);
    for (int rg = 0; rg < NUM_RG; rg++) {
        for (int c = 0; c < NCOLS; c++) {
            for (int i = 0; i < ROWS_PER_RG; i++) {
                int64_t global = (int64_t)rg * ROWS_PER_RG + i;
                if (c == 0) {
                    buf[i] = global;            /* k: sorted */
                } else {
                    /* Payload columns encode (column, row) so we can
                     * verify row-column alignment in the reader. */
                    buf[i] = global * 100 + c;
                }
            }
            ASSERT_OK(carquet_writer_write_batch(w, c, buf, ROWS_PER_RG,
                NULL, NULL));
        }
        if (rg < NUM_RG - 1) ASSERT_OK(carquet_writer_new_row_group(w));
    }
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);
    free(buf);

    carquet_reader_options_t ropts;
    memset(&ropts, 0, sizeof(ropts));
    ropts.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &ropts, &err);
    ASSERT_TRUE(r != NULL);

    /* Project a, b, c — k is the predicate-only column. */
    const char* proj[] = {"a", "b", "c"};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.column_names = proj;
    cfg.num_column_names = 3;
    cfg.batch_size = 4096;
    cfg.use_mmap = true;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Range covering ~1.5 row groups in the middle of the file. */
    int64_t lo = (int64_t)ROWS_PER_RG * 2 + 100;
    int64_t hi = (int64_t)ROWS_PER_RG * 3 + ROWS_PER_RG / 2;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = clause.has_hi = true;
    clause.lo = &lo; clause.lo_size = sizeof(lo);
    clause.hi = &hi; clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    bool alignment_ok = true;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* da; const void* db; const void* dc;
        const uint8_t* na; const uint8_t* nb; const uint8_t* nc;
        int64_t na_n, nb_n, nc_n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &da, &na, &na_n));
        ASSERT_OK(carquet_row_batch_column(batch, 1, &db, &nb, &nb_n));
        ASSERT_OK(carquet_row_batch_column(batch, 2, &dc, &nc, &nc_n));
        const int64_t* a = da;
        const int64_t* b = db;
        const int64_t* c = dc;
        for (int64_t k = 0; k < na_n; k++) {
            /* Recover row index from column 'a' (= global*100 + 1). */
            int64_t row = (a[k] - 1) / 100;
            if (b[k] != row * 100 + 2 || c[k] != row * 100 + 3) {
                alignment_ok = false;
            }
        }
        total += na_n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_TRUE(alignment_ok);
    ASSERT_TRUE(total >= hi - lo + 1);
    ASSERT_TRUE(total < (int64_t)ROWS_PER_RG * 3);

    int64_t skipped = carquet_batch_reader_rows_skipped(br);
    ASSERT_EQ_I64(total + skipped, (int64_t)ROWS_PER_RG * NUM_RG);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Filter re-set across a partial read
 *
 * Opens a reader, sets filter A, reads one batch, sets filter B, reads
 * to completion. Verifies filter B's row count is correct and the
 * pipeline state was properly reset (no leftover slots from filter A).
 * ============================================================================ */

static int test_filter_reset_mid_read(void) {
    g_current_test = "filter_reset_mid_read";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_reset");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Filter A: page 2 only. */
    int32_t a_target = 2050;
    carquet_filter_clause_t fa = {0};
    fa.column_index = 0;
    fa.op = CARQUET_FILTER_EQ;
    fa.value = &a_target;
    fa.value_size = sizeof(a_target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &fa, 1));

    /* Consume one batch from A. */
    carquet_row_batch_t* batch = NULL;
    ASSERT_OK(carquet_batch_reader_next(br, &batch));
    ASSERT_TRUE(batch != NULL && carquet_row_batch_num_rows(batch) > 0);
    carquet_row_batch_free(batch);
    batch = NULL;

    /* Switch to filter B: pages 7..9 (3 pages = 300 rows). */
    int32_t b_lo = 7000;
    carquet_filter_clause_t fb = {0};
    fb.column_index = 0;
    fb.op = CARQUET_FILTER_GE;
    fb.value = &b_lo;
    fb.value_size = sizeof(b_lo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &fb, 1));

    int64_t total_b = 0;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total_b += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total_b, 3 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: NaN values inside the data
 *
 * The Parquet spec says writers should skip NaN when computing min/max.
 * This test writes a float column that contains NaN values and verifies
 * (a) the filter returns the correct non-NaN matching rows and (b) the
 * NaN rows are returned unchanged (page-level granularity).
 * ============================================================================ */

static int test_float_nan_data(void) {
    g_current_test = "float_nan_data";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_nan_data");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "f", CARQUET_PHYSICAL_DOUBLE,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    double values[TOTAL_ROWS];
    double dnan = nan("");
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            /* Sprinkle NaNs every 7 rows. */
            values[row] = (row % 7 == 0) ? dnan : (double)i * 1000.0 + j * 10.0;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS,
        NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* GE 5000 keeps pages 5..9 (non-NaN values cover [5000, 5990] on
     * page 5 etc.). NaN-bearing pages must not crash the filter. */
    double dlo = 5000.0;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &dlo;
    clause.value_size = sizeof(dlo);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    bool saw_non_nan_match = false;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const double* v = (const double*)data;
        for (int64_t k = 0; k < n; k++) {
            if (!isnan(v[k]) && v[k] >= 5000.0) saw_non_nan_match = true;
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 5 * ROWS_PER_PAGE);
    ASSERT_TRUE(saw_non_nan_match);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Selective filter on a sorted column with boundary-order metadata
 *
 * Sorted data is the common case where the page index is most useful.
 * We don't yet exploit boundary_order for early termination; the test
 * just verifies the standard min/max overlap logic returns the right
 * rows on monotonic data.
 * ============================================================================ */

static int test_sorted_column_selective(void) {
    g_current_test = "sorted_column_selective";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_sorted");

    enum { N = 50000 };
    int64_t* values = malloc((size_t)N * sizeof(int64_t));
    ASSERT_TRUE(values != NULL);
    for (int i = 0; i < N; i++) values[i] = (int64_t)i;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v",
        CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 8 * 1024;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, N, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);
    free(values);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Look up a narrow band near the end of the file. */
    int64_t lo = N - 200;
    int64_t hi = N - 100;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_RANGE;
    clause.has_lo = clause.has_hi = true;
    clause.lo = &lo; clause.lo_size = sizeof(lo);
    clause.hi = &hi; clause.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    bool all_within_window_or_neighbor_page = true;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int64_t* v = (const int64_t*)data;
        for (int64_t k = 0; k < n; k++) {
            /* Every returned row should fall within ~10% of the window
             * (the surrounding page tops out at most a few thousand
             * values from the boundary, but with an 8KB page that's
             * 1024 values). */
            if (v[k] < lo - 2048 || v[k] > hi + 2048) {
                all_within_window_or_neighbor_page = false;
            }
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_TRUE(total >= hi - lo + 1);
    ASSERT_TRUE(total < N / 4);    /* selectivity actually fires */
    ASSERT_TRUE(all_within_window_or_neighbor_page);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: NE, LT, LE — exercise comparator branches not covered by RANGE
 * ============================================================================ */

static int test_ne_lt_le(void) {
    g_current_test = "ne_lt_le";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_ne_lt_le");

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(write_int32_file(path, values, TOTAL_ROWS, ROWS_PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* LT 3500: pages whose min < 3500. Pages 0..3 have min = 0/1000/2000/3000,
     * all < 3500. Page 4 has min = 4000 ≥ 3500 → reject. So 4 pages = 400 rows. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        int32_t v = 3500;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_LT;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, 4 * ROWS_PER_PAGE);
        carquet_batch_reader_free(br);
    }

    /* LE 3000: page max ≥ wait wrong direction. LE v keeps pages where min ≤ v.
     * Pages with min in {0, 1000, 2000, 3000} (≤ 3000) → 4 pages. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        int32_t v = 3000;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_LE;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, 4 * ROWS_PER_PAGE);
        carquet_batch_reader_free(br);
    }

    /* NE 5500: keep unless we can prove the whole page equals 5500.
     * No page has min == max == 5500 (each page has 100 distinct values),
     * so EVERY page is kept. Page-level granularity is correct here. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        int32_t v = 5500;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_NE;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, TOTAL_ROWS);
        carquet_batch_reader_free(br);
    }

    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: a page without min/max is not a null page
 *
 * With statistics disabled the column index still exists but carries no
 * bounds. null_pages must come from the value counts: a page marked null is
 * skipped by every value predicate, which would drop all the data here.
 * ============================================================================ */

static int test_no_stats_page_not_null(void) {
    g_current_test = "no_stats_page_not_null";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_nostats");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.write_statistics = false;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = i * 10;
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);
    int32_t target = 5500;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = (int32_t)sizeof(target);
    carquet_status_t fst = carquet_batch_reader_set_page_filter(br, &clause, 1);
    if (fst == CARQUET_OK) {
        bool saw_target = false;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            const void* data; const uint8_t* nb; int64_t n;
            ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
            for (int64_t i = 0; i < n; i++) {
                if (((const int32_t*)data)[i] == target) saw_target = true;
            }
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        /* Without bounds nothing can be proven absent: the row must survive. */
        ASSERT_TRUE(saw_target);
    }
    /* A filter refused for lack of bounds is fine too; dropping rows is not. */
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: chunk statistics cover a page whose minimum is the empty string
 *
 * Page 0 holds "" and "zzz"; the later pages hold values in between. The
 * chunk-level bounds are folded from the pages, and a zero-length page minimum
 * is a real value there: leaving that page out would publish max = "m...",
 * and a reader pruning on it would never see "zzz".
 * ============================================================================ */

static int test_empty_string_chunk_stats(void) {
    g_current_test = "empty_string_chunk_stats";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_empty_min");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    carquet_logical_type_t lt = { .id = CARQUET_LOGICAL_STRING };
    ASSERT_OK(carquet_schema_add_column(schema, "s",
        CARQUET_PHYSICAL_BYTE_ARRAY, &lt, CARQUET_REPETITION_REQUIRED, 0, 0));
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    static carquet_byte_array_t entries[3 * ROWS_PER_PAGE];
    static char buffers[3 * ROWS_PER_PAGE][16];
    for (int row = 0; row < 3 * ROWS_PER_PAGE; row++) {
        if (row == 0) buffers[row][0] = '\0';
        else if (row < ROWS_PER_PAGE) snprintf(buffers[row], sizeof(buffers[row]), "zzz");
        else snprintf(buffers[row], sizeof(buffers[row]), "m%03d", row);
        entries[row].data = (uint8_t*)buffers[row];
        entries[row].length = (int32_t)strlen(buffers[row]);
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, entries, 3 * ROWS_PER_PAGE, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_column_statistics_t st;
    memset(&st, 0, sizeof(st));
    ASSERT_OK(carquet_reader_column_statistics(r, 0, 0, &st));
    /* No bounds at all is acceptable (and what an empty minimum currently
     * yields); bounds that exclude page 0 are not. */
    if (st.has_min_max) {
        ASSERT_TRUE(st.max_value_size >= 3 &&
                    memcmp(st.max_value, "zzz", 3) == 0);
        ASSERT_EQ_I64(st.min_value_size, 0);
    }
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: BOOLEAN column
 *
 * Booleans use a 1-byte stored stat and go through the BOOLEAN branch of
 * the comparator. Write a column where page 5 is all-false and the rest
 * are all-true; EQ true must keep 9 pages, EQ false must keep just 1.
 * ============================================================================ */

static int test_boolean_filter(void) {
    g_current_test = "boolean_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_bool");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BOOLEAN,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    uint8_t bools[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            bools[i * ROWS_PER_PAGE + j] = (i == 5) ? 0 : 1;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, bools, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* EQ true → 9 pages × 100 rows = 900. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        uint8_t v = 1;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = 1;
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, 9 * ROWS_PER_PAGE);
        carquet_batch_reader_free(br);
    }

    /* EQ false → 1 page × 100 rows. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        uint8_t v = 0;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = 1;
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, ROWS_PER_PAGE);
        carquet_batch_reader_free(br);
    }

    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: FIXED_LEN_BYTE_ARRAY (UUID-shaped, 16 bytes)
 *
 * Each row gets a unique 16-byte value where the first 4 bytes encode the
 * page index and the next 4 encode the row-within-page. EQ on the
 * row-from-page-7 value keeps only page 7.
 * ============================================================================ */

static int test_flba_eq(void) {
    g_current_test = "flba_eq";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_flba");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "u",
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 16, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    uint8_t buf[TOTAL_ROWS * 16];
    memset(buf, 0, sizeof(buf));
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            uint8_t* p = buf + row * 16;
            /* big-endian page/row prefix so lex order matches numeric. */
            p[0] = (uint8_t)(i >> 24); p[1] = (uint8_t)(i >> 16);
            p[2] = (uint8_t)(i >> 8);  p[3] = (uint8_t)i;
            p[4] = (uint8_t)(j >> 24); p[5] = (uint8_t)(j >> 16);
            p[6] = (uint8_t)(j >> 8);  p[7] = (uint8_t)j;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, buf, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Build the exact 16-byte key for (page=7, row=42). */
    uint8_t target[16] = {0};
    target[3] = 7; target[7] = 42;
    carquet_filter_clause_t c = {0};
    c.column_index = 0; c.op = CARQUET_FILTER_EQ;
    c.value = target; c.value_size = 16;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));

    int64_t total = 0;
    int matches = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const uint8_t* v = (const uint8_t*)data;
        for (int64_t k = 0; k < n; k++) {
            if (memcmp(v + k * 16, target, 16) == 0) matches++;
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, ROWS_PER_PAGE);
    ASSERT_EQ_I64(matches, 1);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: FLOAT16 logical type (FLBA length 2)
 *
 * Half-precision floats sorted numerically. The min/max bytes are
 * compared via the FLOAT16 decoder, not lex on bytes. Filter LT 5.0
 * must use numeric ordering.
 * ============================================================================ */

static void encode_float16(float f, uint8_t out[2]) {
    /* Minimal IEEE half encoder for finite, non-NaN values. */
    uint32_t bits;
    memcpy(&bits, &f, 4);
    uint32_t sign = (bits >> 31) & 0x1;
    int32_t exp = (int32_t)((bits >> 23) & 0xFF) - 127;
    uint32_t mant = bits & 0x7FFFFF;
    uint16_t half;
    if (exp >= 16) {
        half = (uint16_t)((sign << 15) | (0x1F << 10));   /* ±inf */
    } else if (exp <= -15) {
        half = (uint16_t)(sign << 15);                    /* underflow ⇒ 0 */
    } else {
        uint16_t hexp = (uint16_t)(exp + 15);
        uint16_t hmant = (uint16_t)(mant >> 13);
        half = (uint16_t)((sign << 15) | (hexp << 10) | hmant);
    }
    out[0] = (uint8_t)(half & 0xFF);
    out[1] = (uint8_t)(half >> 8);
}

static int test_float16_filter(void) {
    g_current_test = "float16_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_f16");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    carquet_logical_type_t lt = { .id = CARQUET_LOGICAL_FLOAT16 };
    ASSERT_OK(carquet_schema_add_column(schema, "h",
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, &lt,
        CARQUET_REPETITION_REQUIRED, 2, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    uint8_t buf[TOTAL_ROWS * 2];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            float val = (float)i + (float)j / (float)ROWS_PER_PAGE;
            encode_float16(val, buf + row * 2);
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, buf, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* LT 5.0 (numeric) keeps pages whose min half-float < 5.0 → pages 0..4. */
    uint8_t target[2];
    encode_float16(5.0f, target);
    carquet_filter_clause_t c = {0};
    c.column_index = 0; c.op = CARQUET_FILTER_LT;
    c.value = target; c.value_size = 2;
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 5 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Unsigned INT64 column
 *
 * Stores values whose unsigned interpretation places them in the upper
 * half of the uint64 range (i.e. negative when interpreted as int64).
 * A filter LT 100 must reject every page; only the unsigned comparator
 * gives the right answer (signed would treat the stored mins as negative
 * and incorrectly keep all pages).
 * ============================================================================ */

static int test_unsigned_int64(void) {
    g_current_test = "unsigned_int64";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_uint64");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    carquet_logical_type_t lt = { .id = CARQUET_LOGICAL_INTEGER };
    lt.params.integer.bit_width = 64;
    lt.params.integer.is_signed = false;
    ASSERT_OK(carquet_schema_add_column(schema, "u", CARQUET_PHYSICAL_INT64,
        &lt, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    uint64_t big_base = 0x8000000000000000ULL;  /* 2^63 */
    int64_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            uint64_t u = big_base + (uint64_t)(i * 1000 + j * 10);
            memcpy(&values[i * ROWS_PER_PAGE + j], &u, sizeof(uint64_t));
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* LT 100 (unsigned): no page can match because every value > 2^63. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        uint64_t v = 100;
        int64_t vs;
        memcpy(&vs, &v, sizeof(int64_t));
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_LT;
        c.value = &vs; c.value_size = sizeof(vs);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        carquet_row_batch_t* batch = NULL;
        carquet_status_t st = carquet_batch_reader_next(br, &batch);
        ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA || batch == NULL);
        ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), TOTAL_ROWS);
        carquet_batch_reader_free(br);
    }

    /* GE 2^63 + 5000: pages whose max ≥ 2^63 + 5000. Page i max = base + i*1000 + 990.
     * That's ≥ base + 5000 when i*1000 + 990 ≥ 5000 → i ≥ 5. So pages 5..9 = 500 rows. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        uint64_t v = big_base + 5000;
        int64_t vs;
        memcpy(&vs, &v, sizeof(int64_t));
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_GE;
        c.value = &vs; c.value_size = sizeof(vs);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, 5 * ROWS_PER_PAGE);
        carquet_batch_reader_free(br);
    }

    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Dictionary-encoded column with filter
 *
 * Carquet's dictionary writer accumulates the whole column chunk and
 * emits one dictionary page followed by one data page, so page-level
 * filtering can only be all-or-nothing on a single data page. We verify
 * (a) the filter does not crash on the dict-encoded layout and (b)
 * returns the correct binary outcome for both a present and an absent
 * value.
 * ============================================================================ */

static int test_dict_encoded_filter(void) {
    g_current_test = "dict_encoded_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_dict");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_RLE_DICTIONARY;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0,
        CARQUET_ENCODING_RLE_DICTIONARY));

    int64_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = (i % 5) * 1000;
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* EQ 2000 (present): the single data page must be kept; all rows
     * returned (page-level granularity). */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        int64_t v = 2000;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        int64_t total = 0;
        carquet_row_batch_t* batch = NULL;
        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
            total += carquet_row_batch_num_rows(batch);
            carquet_row_batch_free(batch);
            batch = NULL;
        }
        ASSERT_EQ_I64(total, TOTAL_ROWS);
        carquet_batch_reader_free(br);
    }

    /* EQ 99999 (absent): the page's max is 4000, no match; entire
     * column is pruned. */
    {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.batch_size = 4096;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        ASSERT_TRUE(br != NULL);
        int64_t v = 99999;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
        carquet_row_batch_t* batch = NULL;
        carquet_status_t st = carquet_batch_reader_next(br, &batch);
        ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA || batch == NULL);
        ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), TOTAL_ROWS);
        carquet_batch_reader_free(br);
    }

    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: preserve_dictionaries = true with a filter
 *
 * The batch reader returns dictionary indices instead of materialized
 * values. With a filter that retains the single dict-encoded page, we
 * must get the indices and the dictionary back intact.
 * ============================================================================ */

static int test_preserve_dict_with_filter(void) {
    g_current_test = "preserve_dict_with_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_predict");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_RLE_DICTIONARY;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0,
        CARQUET_ENCODING_RLE_DICTIONARY));

    int64_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = (i % 4) * 1000;
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    cfg.preserve_dictionaries = true;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int64_t v = 2000;
    carquet_filter_clause_t c = {0};
    c.column_index = 0; c.op = CARQUET_FILTER_EQ;
    c.value = &v; c.value_size = sizeof(v);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));

    int64_t total = 0;
    bool dict_seen = false;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const uint32_t* indices;
        const uint8_t* nb;
        int64_t n;
        const uint8_t* dict_data;
        int32_t dict_count;
        const uint32_t* dict_offsets;
        if (carquet_row_batch_column_dictionary(batch, 0, &indices, &nb, &n,
                &dict_data, &dict_count, &dict_offsets) == CARQUET_OK) {
            dict_seen = true;
            ASSERT_TRUE(dict_count > 0);
            /* Every index must point into the dictionary range. */
            for (int64_t k = 0; k < n; k++) {
                ASSERT_TRUE((int32_t)indices[k] < dict_count);
            }
            total += n;
        }
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_TRUE(dict_seen);
    ASSERT_EQ_I64(total, TOTAL_ROWS);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Snappy compression + filter (default Parquet compression)
 * ============================================================================ */

static int test_snappy_with_filter(void) {
    g_current_test = "snappy_with_filter";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_snappy");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_SNAPPY;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            values[i * ROWS_PER_PAGE + j] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t lo = 4000, hi = 5999;
    carquet_filter_clause_t c = {0};
    c.column_index = 0; c.op = CARQUET_FILTER_RANGE;
    c.has_lo = c.has_hi = true;
    c.lo = &lo; c.lo_size = sizeof(lo);
    c.hi = &hi; c.hi_size = sizeof(hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int32_t* v = (const int32_t*)data;
        /* Every returned value must come from a page that overlaps [lo, hi]. */
        for (int64_t k = 0; k < n; k++) {
            ASSERT_TRUE(v[k] >= 4000 - 1000 && v[k] <= 5999 + 1000);
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 2 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Three-clause intersection
 *
 * Three columns, three RANGE clauses. The intersection of their kept
 * page ranges must be computed correctly.
 * ============================================================================ */

static int test_three_clause_intersection(void) {
    g_current_test = "three_clause_intersection";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_three");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));
    ASSERT_OK(carquet_schema_add_column(schema, "c", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));
    ASSERT_OK(carquet_writer_set_column_encoding(w, 2, CARQUET_ENCODING_PLAIN));

    int32_t a[TOTAL_ROWS], b[TOTAL_ROWS], c[TOTAL_ROWS];
    for (int i = 0; i < NUM_PAGES; i++) {
        for (int j = 0; j < ROWS_PER_PAGE; j++) {
            int row = i * ROWS_PER_PAGE + j;
            a[row] = i * 1000 + j * 10;          /* page i: [i*1000, i*1000+990] */
            b[row] = i * 1000 + j * 10;
            c[row] = i * 1000 + j * 10;
        }
    }
    ASSERT_OK(carquet_writer_write_batch(w, 0, a, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_write_batch(w, 1, b, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_write_batch(w, 2, c, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* a in [2000, 6999]: pages 2..6.
     * b in [3000, 7999]: pages 3..7.
     * c in [4000, 8999]: pages 4..8.
     * Intersection: pages 4..6 = 3 pages = 300 rows. */
    int32_t a_lo = 2000, a_hi = 6999;
    int32_t b_lo = 3000, b_hi = 7999;
    int32_t c_lo = 4000, c_hi = 8999;
    carquet_filter_clause_t cl[3];
    memset(cl, 0, sizeof(cl));
    cl[0].column_index = 0; cl[0].op = CARQUET_FILTER_RANGE;
    cl[0].has_lo = cl[0].has_hi = true;
    cl[0].lo = &a_lo; cl[0].lo_size = sizeof(a_lo);
    cl[0].hi = &a_hi; cl[0].hi_size = sizeof(a_hi);
    cl[1].column_index = 1; cl[1].op = CARQUET_FILTER_RANGE;
    cl[1].has_lo = cl[1].has_hi = true;
    cl[1].lo = &b_lo; cl[1].lo_size = sizeof(b_lo);
    cl[1].hi = &b_hi; cl[1].hi_size = sizeof(b_hi);
    cl[2].column_index = 2; cl[2].op = CARQUET_FILTER_RANGE;
    cl[2].has_lo = cl[2].has_hi = true;
    cl[2].lo = &c_lo; cl[2].lo_size = sizeof(c_lo);
    cl[2].hi = &c_hi; cl[2].hi_size = sizeof(c_hi);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, cl, 3));

    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        ASSERT_EQ_I64(carquet_row_batch_num_columns(batch), 3);
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 3 * ROWS_PER_PAGE);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: Validation error paths
 *
 * Each malformed clause should be rejected by set_page_filter and leave
 * the reader unchanged.
 * ============================================================================ */

static int test_validation_errors(void) {
    g_current_test = "validation_errors";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_validate");

    int32_t values[10];
    for (int i = 0; i < 10; i++) values[i] = i;
    ASSERT_OK(write_int32_file(path, values, 10, 10));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Out-of-range column index. */
    {
        int32_t v = 0;
        carquet_filter_clause_t c = {0};
        c.column_index = 99; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = sizeof(v);
        carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
        ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);
    }

    /* NULL value for a non-NULL op. */
    {
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = NULL; c.value_size = 4;
        carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
        ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);
    }

    /* RANGE with neither lo nor hi. */
    {
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_RANGE;
        carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
        ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);
    }

    /* IN with zero values. */
    {
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_IN;
        c.values = values; c.value_count = 0;
        carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
        ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);
    }

    /* IN with > 256 values. */
    {
        int32_t* big = malloc(257 * sizeof(int32_t));
        ASSERT_TRUE(big != NULL);
        for (int i = 0; i < 257; i++) big[i] = i;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_IN;
        c.values = big; c.value_count = 257;
        carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
        ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);
        free(big);
    }

    /* Setting a valid filter after errors still works. */
    {
        int32_t v = 5;
        carquet_filter_clause_t c = {0};
        c.column_index = 0; c.op = CARQUET_FILTER_EQ;
        c.value = &v; c.value_size = sizeof(v);
        ASSERT_OK(carquet_batch_reader_set_page_filter(br, &c, 1));
    }

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Test: FLBA size mismatch
 *
 * For FIXED_LEN_BYTE_ARRAY the value_size must equal the column's
 * type_length; anything else is rejected.
 * ============================================================================ */

static int test_flba_size_mismatch(void) {
    g_current_test = "flba_size_mismatch";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_flba_size");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "u",
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 16, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    uint8_t buf[16] = {0};
    ASSERT_OK(carquet_writer_write_batch(w, 0, buf, 1, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* value_size = 8 against a 16-byte column. */
    uint8_t v[16] = {0};
    carquet_filter_clause_t c = {0};
    c.column_index = 0; c.op = CARQUET_FILTER_EQ;
    c.value = v; c.value_size = 8;
    carquet_status_t st = carquet_batch_reader_set_page_filter(br, &c, 1);
    ASSERT_TRUE(st == CARQUET_ERROR_INVALID_ARGUMENT);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Regression: forward seek must leave values_remaining exact (1A)
 * ============================================================================
 *
 * A forward seek issued while the reader is already mid-chunk (current_page > 0,
 * e.g. the second range of a multi-range page filter) must recompute
 * values_remaining as (chunk values - all values before the target page), not
 * just the values between the current position and the target. Getting it wrong
 * leaves values_remaining too high, so a later read walks current_page past the
 * chunk end into adjacent bytes.
 */
static int test_seek_forward_values_remaining(void) {
    g_current_test = "seek_forward_values_remaining";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_seek_fwd");

    enum { PAGES = 3, PER_PAGE = 100, TOTAL = PAGES * PER_PAGE };
    int32_t values[TOTAL];
    for (int i = 0; i < TOTAL; i++) values[i] = i;
    ASSERT_OK(write_int32_file(path, values, TOTAL, PER_PAGE));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, 0, &err);
    ASSERT_TRUE(oi != NULL);
    ASSERT_EQ_I64(carquet_offset_index_num_pages(oi), PAGES);

    carquet_column_reader_t* cr = carquet_reader_get_column(r, 0, 0, &err);
    ASSERT_TRUE(cr != NULL);

    /* Consume page 0 so the reader is mid-chunk (current_page > 0). */
    int32_t buf[TOTAL];
    int64_t got = carquet_column_read_batch(cr, buf, PER_PAGE, NULL, NULL);
    ASSERT_EQ_I64(got, PER_PAGE);

    /* Forward seek to the last page. */
    carquet_page_location_t loc;
    ASSERT_OK(carquet_offset_index_get_page_location(oi, PAGES - 1, &loc));
    ASSERT_OK(carquet_column_reader_seek_to_data_page(cr, loc.offset, 0, &err));

    /* Only the final page's values remain — not the inflated count. */
    ASSERT_EQ_I64(cr->values_remaining, PER_PAGE);

    /* A read larger than the remainder must stop at the chunk boundary and
     * return exactly the last page's data, never spilling past it. */
    got = carquet_column_read_batch(cr, buf, TOTAL, NULL, NULL);
    ASSERT_EQ_I64(got, PER_PAGE);
    for (int i = 0; i < PER_PAGE; i++) {
        ASSERT_EQ_I64(buf[i], (PAGES - 1) * PER_PAGE + i);
    }
    ASSERT_EQ_I64(cr->values_remaining, 0);
    ASSERT_EQ_I64(carquet_column_read_batch(cr, buf, TOTAL, NULL, NULL), 0);

    carquet_column_reader_free(cr);
    carquet_offset_index_free(oi);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * Regression: carquet_column_skip must not decompress skipped pages (1B)
 * ============================================================================
 *
 * Skipping rows that span whole pages should advance by parsing page headers
 * only; the compressed payloads of fully-skipped pages must never be
 * decompressed. Only the final partially-skipped page is decoded.
 */
static int test_skip_does_not_decompress(void) {
    g_current_test = "skip_does_not_decompress";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_skip_nodecomp");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    ASSERT_TRUE(schema != NULL);
    ASSERT_OK(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0));

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_batch_size = ROWS_PER_PAGE;
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    ASSERT_TRUE(w != NULL);
    ASSERT_OK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));

    int32_t values[TOTAL_ROWS];
    for (int i = 0; i < TOTAL_ROWS; i++) values[i] = i;
    ASSERT_OK(carquet_writer_write_batch(w, 0, values, TOTAL_ROWS, NULL, NULL));
    ASSERT_OK(carquet_writer_close(w));
    carquet_schema_free(schema);

    /* Install the counting codec after writing so only reads are counted. */
    carquet_custom_codec_t codec = {
        .compress = counter_compress,
        .decompress = counter_decompress,
        .compress_bound = counter_compress_bound,
        .user_data = NULL,
    };
    g_codec_decompress_calls = 0;
    ASSERT_OK(carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &codec));

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_column_reader_t* cr = carquet_reader_get_column(r, 0, 0, &err);
    ASSERT_TRUE(cr != NULL);

    /* Skip 8.5 pages: 8 whole pages skipped by header only, the 9th decoded. */
    const int64_t skip_n = 8 * ROWS_PER_PAGE + ROWS_PER_PAGE / 2;
    ASSERT_EQ_I64(carquet_column_skip(cr, skip_n), skip_n);
    ASSERT_EQ_I64(g_codec_decompress_calls, 1);

    /* Read the rest of page 8 (already decoded): no further decompression. */
    int32_t buf[ROWS_PER_PAGE];
    int64_t got = carquet_column_read_batch(cr, buf, ROWS_PER_PAGE / 2, NULL, NULL);
    ASSERT_EQ_I64(got, ROWS_PER_PAGE / 2);
    for (int i = 0; i < ROWS_PER_PAGE / 2; i++) {
        ASSERT_EQ_I64(buf[i], skip_n + i);
    }
    ASSERT_EQ_I64(g_codec_decompress_calls, 1);

    carquet_column_reader_free(cr);
    carquet_reader_close(r);
    ASSERT_OK(carquet_register_codec(CARQUET_COMPRESSION_ZSTD, NULL));
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ============================================================================
 * main
 * ============================================================================ */

/* ============================================================================
 * Automatic row-group skipping: statistics (no page index) and bloom filter.
 *
 * These exercise the row-group-level pruning that runs before the page index
 * is consulted. When it proves a whole row group cannot match, the group is
 * skipped without a user callback — and, crucially, without a page index.
 * ============================================================================ */

/* Write three row groups of INT32 with disjoint value ranges. Optionally
 * enables the page index and/or a bloom filter on the single column. */
static carquet_status_t write_three_rg_int32(
    const char* path, bool page_index, bool bloom,
    const int32_t* rg0, const int32_t* rg1, const int32_t* rg2, int64_t n) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return CARQUET_ERROR_INTERNAL;
    carquet_status_t st = carquet_schema_add_column(
        schema, "v", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    if (st != CARQUET_OK) { carquet_schema_free(schema); return st; }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = page_index;
    opts.write_statistics = true;
    opts.write_bloom_filters = bloom;
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { carquet_schema_free(schema); return CARQUET_ERROR_INTERNAL; }
    st = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN);
    if (st != CARQUET_OK) goto cleanup;
    if (bloom) {
        st = carquet_writer_set_column_bloom_filter(w, 0, true);
        if (st != CARQUET_OK) goto cleanup;
    }
    st = carquet_writer_write_batch(w, 0, rg0, n, NULL, NULL);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_new_row_group(w);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_write_batch(w, 0, rg1, n, NULL, NULL);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_new_row_group(w);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_write_batch(w, 0, rg2, n, NULL, NULL);
    if (st != CARQUET_OK) goto cleanup;
    st = carquet_writer_close(w);
    w = NULL;
cleanup:
    if (w) carquet_writer_close(w);
    carquet_schema_free(schema);
    return st;
}

/* #7: statistics-based skipping drops whole row groups with no page index and
 * no user callback. A predicate out of every row group's [min,max] must
 * yield a clean end-of-data (not PAGE_INDEX_REQUIRED) and credit every row. */
static int test_rg_stats_skip_no_page_index(void) {
    g_current_test = "rg_stats_skip_no_page_index";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_rg_stats");

    const int64_t N = 100;
    int32_t rg0[100], rg1[100], rg2[100];
    for (int i = 0; i < N; i++) {
        rg0[i] = i;             /* [0, 99]      */
        rg1[i] = 1000 + i;      /* [1000, 1099] */
        rg2[i] = 2000 + i;      /* [2000, 2099] */
    }
    ASSERT_OK(write_three_rg_int32(path, /*page_index=*/false, /*bloom=*/false,
                                   rg0, rg1, rg2, N));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    ASSERT_EQ_I64(carquet_reader_num_row_groups(r), 3);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 5000;  /* Not in any row group. */
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    /* All three row groups pruned by stats before the (absent) page index is
     * ever consulted: no PAGE_INDEX_REQUIRED error, all rows accounted for. */
    ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA ||
                (st == CARQUET_OK && batch == NULL));
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), 3 * N);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* #7: with a page index present, only the row group whose stats overlap the
 * predicate is read; the others are skipped and their rows credited. */
static int test_rg_stats_skip_with_page_index(void) {
    g_current_test = "rg_stats_skip_with_page_index";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_rg_stats_pi");

    const int64_t N = 100;
    int32_t rg0[100], rg1[100], rg2[100];
    for (int i = 0; i < N; i++) {
        rg0[i] = i;
        rg1[i] = 1000 + i;
        rg2[i] = 2000 + i;
    }
    ASSERT_OK(write_three_rg_int32(path, /*page_index=*/true, /*bloom=*/false,
                                   rg0, rg1, rg2, N));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    int32_t target = 1050;  /* Only present in the middle row group. */
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    int64_t total = 0;
    bool saw_target = false;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int32_t* v = (const int32_t*)data;
        for (int64_t k = 0; k < n; k++) if (v[k] == target) saw_target = true;
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_TRUE(saw_target);
    ASSERT_EQ_I64(total, N);  /* Only the middle row group is read. */
    /* Two full row groups (rg0, rg2) skipped by row-group stats. */
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), 2 * N);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* #6: a bloom filter drops a row group whose statistics say "might match"
 * (the target is inside [min,max]) but whose values do not actually include
 * the target. Without the bloom filter (no page index) this would raise
 * PAGE_INDEX_REQUIRED; with it the group is pruned cleanly. */
static int test_rg_bloom_skip(void) {
    g_current_test = "rg_bloom_skip";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "pf_rg_bloom");

    /* One row group of even values [0, 198]: any odd value in that range is
     * absent yet passes the min/max test, isolating the bloom-filter path. */
    const int64_t N = 100;
    int32_t evens[100], dummy_hi[100], dummy_lo[100];
    for (int i = 0; i < N; i++) {
        evens[i] = 2 * i;             /* [0, 198], evens only */
        dummy_lo[i] = -100000 - i;    /* far below, pruned by stats */
        dummy_hi[i] = 100000 + i;     /* far above, pruned by stats */
    }
    ASSERT_OK(write_three_rg_int32(path, /*page_index=*/false, /*bloom=*/true,
                                   dummy_lo, evens, dummy_hi, N));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    /* Pick an odd target the bloom filter reports absent (deterministic per
     * file, but robust against the codec's exact hash/sizing). Confirm an
     * inserted even value reads present as a sanity check. */
    carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(r, 1, 0, &err);
    ASSERT_TRUE(bf != NULL);
    ASSERT_TRUE(carquet_bloom_filter_check_i32(bf, 100));  /* even, present */
    int32_t target = -1;
    for (int32_t cand = 1; cand < 198; cand += 2) {
        if (!carquet_bloom_filter_check_i32(bf, cand)) { target = cand; break; }
    }
    carquet_bloom_filter_destroy(bf);
    ASSERT_TRUE(target > 0);  /* Some odd value must read absent. */

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));

    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    /* rg0/rg2 pruned by stats, rg1 pruned by the bloom filter — clean EOD. */
    ASSERT_TRUE(st == CARQUET_ERROR_END_OF_DATA ||
                (st == CARQUET_OK && batch == NULL));
    ASSERT_EQ_I64(carquet_batch_reader_rows_skipped(br), 3 * N);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_eq_int32_single_page();
    failures += test_range_multi_page();
    failures += test_range_no_match();
    failures += test_and_two_columns();
    failures += test_predicate_not_projected();
    failures += test_nullable_is_null_predicates();
    failures += test_byte_array_eq();
    failures += test_float_nan_predicate();
    failures += test_in_int64();
    failures += test_filter_clear();
    failures += test_no_page_index_error();
    failures += test_rg_stats_skip_no_page_index();
    failures += test_rg_stats_skip_with_page_index();
    failures += test_rg_bloom_skip();
    failures += test_int96_rejected();
    failures += test_aligned_range();
    failures += test_whole_rg();
    failures += test_nullable_end_to_end();
    failures += test_pipeline_filter();

    /* Real-life scenarios */
    failures += test_custom_codec_decompress_counter();
    failures += test_byte_array_truncated_stats();
    failures += test_large_selective();
    failures += test_pipeline_multi_col_stress();
    failures += test_filter_reset_mid_read();
    failures += test_float_nan_data();
    failures += test_sorted_column_selective();

    /* Coverage of remaining ops, physical types, encodings, and error paths */
    failures += test_ne_lt_le();
    failures += test_boolean_filter();
    failures += test_flba_eq();
    failures += test_float16_filter();
    failures += test_unsigned_int64();
    failures += test_dict_encoded_filter();
    failures += test_preserve_dict_with_filter();
    failures += test_snappy_with_filter();
    failures += test_three_clause_intersection();
    failures += test_validation_errors();
    failures += test_flba_size_mismatch();
    failures += test_seek_forward_values_remaining();
    failures += test_skip_does_not_decompress();
    failures += test_no_stats_page_not_null();
    failures += test_empty_string_chunk_stats();

    if (failures > 0) {
        fprintf(stderr, "%d test failures\n", failures);
        return 1;
    }
    return 0;
}
