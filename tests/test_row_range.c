/**
 * @file test_row_range.c
 * @brief Tests for the public row-range [offset, limit) read API
 *        (carquet_batch_reader_set_row_range).
 *
 * Files are written so that the value in row i equals i (its global row
 * index), which makes every window trivially verifiable. Multi-row-group
 * files are produced with carquet_writer_new_row_group() so windows that
 * straddle row-group boundaries are exercised.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

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

/* Write `rg_count` row groups of `rows_per_rg` INT32 rows, value == global row
 * index. rows_per_page controls page splitting (0 = writer default). */
static carquet_status_t write_indexed_file(
    const char* path, int rg_count, int rows_per_rg,
    int64_t rows_per_page, carquet_compression_t codec,
    carquet_encoding_t enc) {

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return CARQUET_ERROR_INTERNAL;
    carquet_status_t st = carquet_schema_add_column(
        schema, "v", CARQUET_PHYSICAL_INT32, NULL,
        CARQUET_REPETITION_REQUIRED, 0, 0);
    if (st != CARQUET_OK) { carquet_schema_free(schema); return st; }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codec;
    if (rows_per_page > 0) {
        opts.write_page_index = true;
        opts.max_rows_per_page = rows_per_page;
        opts.write_batch_size = rows_per_page;
    }
    opts.dictionary_encoding = enc;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { carquet_schema_free(schema); return CARQUET_ERROR_INTERNAL; }
    st = carquet_writer_set_column_encoding(w, 0, enc);
    if (st != CARQUET_OK) goto cleanup;

    int32_t* vals = malloc(sizeof(int32_t) * (size_t)rows_per_rg);
    if (!vals) { st = CARQUET_ERROR_OUT_OF_MEMORY; goto cleanup; }
    for (int g = 0; g < rg_count; g++) {
        for (int j = 0; j < rows_per_rg; j++) {
            vals[j] = g * rows_per_rg + j;
        }
        st = carquet_writer_write_batch(w, 0, vals, rows_per_rg, NULL, NULL);
        if (st != CARQUET_OK) { free(vals); goto cleanup; }
        if (g + 1 < rg_count) {
            st = carquet_writer_new_row_group(w);
            if (st != CARQUET_OK) { free(vals); goto cleanup; }
        }
    }
    free(vals);
    st = carquet_writer_close(w);
cleanup:
    carquet_schema_free(schema);
    return st;
}

/* Read [offset, limit) with the given batch size and assert the returned
 * values are exactly the contiguous run [offset, offset+expect_count). */
static int check_window(const char* path, int64_t offset, int64_t limit,
                        int32_t batch_size, int64_t expect_count) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = batch_size;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    ASSERT_OK(carquet_batch_reader_set_row_range(br, offset, limit));

    int64_t total = 0;
    int64_t expected_val = offset;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data; const uint8_t* nb; int64_t n;
        ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
        const int32_t* v = (const int32_t*)data;
        for (int64_t i = 0; i < n; i++) {
            if (v[i] != (int32_t)expected_val) {
                fprintf(stderr, "value mismatch at row %lld: got %d want %lld\n",
                        (long long)(total + i), v[i], (long long)expected_val);
                carquet_row_batch_free(batch);
                carquet_batch_reader_free(br);
                carquet_reader_close(r);
                return 1;
            }
            expected_val++;
        }
        total += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, expect_count);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return 0;
}

/* ---- Test 1: single row group, several windows -------------------------- */
static int test_single_rg_windows(void) {
    g_current_test = "single_rg_windows";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_single");
    ASSERT_OK(write_indexed_file(path, 1, 1000, 0,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    /* interior window, batch smaller than window (multi-batch) */
    if (check_window(path, 100, 250, 64, 250)) return 1;
    /* window from row 0 */
    if (check_window(path, 0, 10, 4096, 10)) return 1;
    /* window to the last row */
    if (check_window(path, 990, 10, 4096, 10)) return 1;
    /* limit past end clamps to available rows */
    if (check_window(path, 900, 500, 4096, 100)) return 1;
    /* limit < 0 means "to end" */
    if (check_window(path, 700, -1, 4096, 300)) return 1;

    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 2: cross-row-group windows ------------------------------------ */
static int test_multi_rg_windows(void) {
    g_current_test = "multi_rg_windows";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_multi");
    /* 4 row groups × 250 rows = 1000 rows total. */
    ASSERT_OK(write_indexed_file(path, 4, 250, 0,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    /* window straddling RG0/RG1 boundary (row 250) */
    if (check_window(path, 200, 100, 4096, 100)) return 1;
    /* window spanning three row groups */
    if (check_window(path, 240, 520, 128, 520)) return 1;
    /* whole file via to-end */
    if (check_window(path, 0, -1, 300, 1000)) return 1;
    /* window entirely inside RG3 */
    if (check_window(path, 800, 50, 4096, 50)) return 1;

    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 3: with page index (offset-index seek path) ------------------- */
static int test_page_index_seek(void) {
    g_current_test = "page_index_seek";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_pidx");
    /* one RG of 1000 rows split into 100-row pages, page index written. */
    ASSERT_OK(write_indexed_file(path, 1, 1000, 100,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    if (check_window(path, 350, 200, 64, 200)) return 1;   /* spans pages 3-5 */
    if (check_window(path, 0, 100, 4096, 100)) return 1;
    if (check_window(path, 950, 100, 4096, 50)) return 1;   /* clamp */

    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 4: compressed file (dictionary column) ------------------------ */
static int test_compressed_dict(void) {
    g_current_test = "compressed_dict";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_zstd");
    ASSERT_OK(write_indexed_file(path, 3, 400, 0,
                                 CARQUET_COMPRESSION_ZSTD,
                                 CARQUET_ENCODING_RLE_DICTIONARY));

    if (check_window(path, 350, 200, 128, 200)) return 1;   /* crosses RG0/RG1 */
    if (check_window(path, 0, -1, 256, 1200)) return 1;
    if (check_window(path, 1100, 200, 4096, 100)) return 1; /* clamp in last RG */

    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 5: boundary cases -------------------------------------------- */
static int test_boundaries(void) {
    g_current_test = "boundaries";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_bound");
    ASSERT_OK(write_indexed_file(path, 1, 100, 0,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    /* limit == 0 → zero rows */
    if (check_window(path, 10, 0, 4096, 0)) return 1;
    /* offset at total → zero rows */
    if (check_window(path, 100, 50, 4096, 0)) return 1;
    /* offset beyond total → zero rows */
    if (check_window(path, 500, 50, 4096, 0)) return 1;

    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 6: clearing + re-positioning a reader ------------------------- */
static int test_clear_and_reposition(void) {
    g_current_test = "clear_and_reposition";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_clear");
    ASSERT_OK(write_indexed_file(path, 1, 500, 0,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* First a window [10,20). */
    ASSERT_OK(carquet_batch_reader_set_row_range(br, 10, 20));
    carquet_row_batch_t* batch = NULL;
    ASSERT_OK(carquet_batch_reader_next(br, &batch));
    ASSERT_TRUE(batch != NULL);
    ASSERT_EQ_I64(carquet_row_batch_num_rows(batch), 20);
    carquet_row_batch_free(batch);

    /* Re-position to a different window on the same reader. */
    ASSERT_OK(carquet_batch_reader_set_row_range(br, 400, 25));
    batch = NULL;
    ASSERT_OK(carquet_batch_reader_next(br, &batch));
    ASSERT_TRUE(batch != NULL);
    const void* data; const uint8_t* nb; int64_t n;
    ASSERT_OK(carquet_row_batch_column(batch, 0, &data, &nb, &n));
    ASSERT_EQ_I64(n, 25);
    ASSERT_EQ_I64(((const int32_t*)data)[0], 400);
    carquet_row_batch_free(batch);

    /* Clear → full sequential read of all 500 rows. */
    ASSERT_OK(carquet_batch_reader_set_row_range(br, 0, -1));
    int64_t total = 0;
    batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    ASSERT_EQ_I64(total, 500);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

/* ---- Test 7: rejection / error paths ------------------------------------ */
static int test_rejections(void) {
    g_current_test = "rejections";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "rr_reject");
    ASSERT_OK(write_indexed_file(path, 1, 100, 100,
                                 CARQUET_COMPRESSION_UNCOMPRESSED,
                                 CARQUET_ENCODING_PLAIN));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    /* Negative offset is rejected. */
    ASSERT_TRUE(carquet_batch_reader_set_row_range(br, -1, 10) ==
                CARQUET_ERROR_INVALID_ARGUMENT);

    /* Mutual exclusion: a page filter blocks a row range. */
    int32_t target = 50;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &target;
    clause.value_size = (int32_t)sizeof(target);
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, &clause, 1));
    ASSERT_TRUE(carquet_batch_reader_set_row_range(br, 0, 10) ==
                CARQUET_ERROR_INVALID_ARGUMENT);
    /* Clear the filter, then a row range is accepted. */
    ASSERT_OK(carquet_batch_reader_set_page_filter(br, NULL, 0));
    ASSERT_OK(carquet_batch_reader_set_row_range(br, 0, 10));
    /* And now installing a filter while a range is active is rejected. */
    ASSERT_TRUE(carquet_batch_reader_set_page_filter(br, &clause, 1) ==
                CARQUET_ERROR_INVALID_ARGUMENT);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(g_current_test);
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_single_rg_windows();
    failures += test_multi_rg_windows();
    failures += test_page_index_seek();
    failures += test_compressed_dict();
    failures += test_boundaries();
    failures += test_clear_and_reposition();
    failures += test_rejections();

    if (failures) {
        fprintf(stderr, "%d row-range test(s) failed\n", failures);
        return 1;
    }
    printf("All row-range tests passed\n");
    return 0;
}
