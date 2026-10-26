/**
 * @file test_fixes_0_7_2.c
 * @brief Regression tests for the v0.7.2 bug fixes and optimisations.
 *
 * One test per fix, each failing on the v0.7.1 code:
 *  - required_list_empty_not_null   batch reader reported an empty list in a
 *                                   REQUIRED LIST as a null list
 *  - dict_nan_first_stats           dictionary FLOAT/DOUBLE chunk whose first
 *                                   value is NaN recorded NaN as min/max
 *  - empty_string_chunk_stats       BYTE_ARRAY chunk whose minimum is "" lost
 *                                   its footer statistics (and, on the
 *                                   dictionary path, its page bounds)
 *  - dict_fallback_not_sticky       a dictionary fallback in one row group made
 *                                   every later row group PLAIN
 *  - list_offset_index_rows         OffsetIndex.first_row_index counted levels,
 *                                   and pages were cut inside a list
 *  - column_index_withheld          a ColumnIndex with empty bounds was written
 *                                   for pages without statistics
 *  - byte_array_page_size           page_size was ignored for BYTE_ARRAY
 *                                   columns written in one call
 *  - rle_bitpacked_grouping         RLE emitted one header per 8 bit-packed
 *                                   values; now up to 63 groups per header
 *  - batch_size_zero_config         a zero-initialised config spun forever
 *  - dispatch_match_copy            NEON/SSE match_copy wrong for offsets 3,5,6,7
 *  - bss_large_tiled_decode         BYTE_STREAM_SPLIT tile moved off the stack
 */

#include <carquet/carquet.h>
#include "core/buffer.h"
#include "encoding/rle.h"
#include "test_helpers.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void carquet_dispatch_match_copy(uint8_t* dst, const uint8_t* src,
                                        size_t len, size_t offset);

/* ------------------------------------------------------------------------ */
static carquet_writer_t* open_writer(const char* path, carquet_schema_t* s,
                                     carquet_writer_options_t* opts,
                                     carquet_error_t* err) {
    carquet_writer_t* w = carquet_writer_create(path, s, opts, err);
    return w;
}

/* ------------------------------------------------------------------------ */
static int test_required_list_empty_not_null(void) {
    const char* name = "required_list_empty_not_null";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_reqlist");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    assert(carquet_schema_add_list(s, "vals", CARQUET_PHYSICAL_INT32, NULL,
                                   CARQUET_REPETITION_REQUIRED, 0, 0) >= 0);
    carquet_writer_t* w = open_writer(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);

    /* rows: [1], [], [2]  -> REQUIRED list of REQUIRED ints: max_def = 1 */
    int32_t values[2] = {1, 2};
    int16_t def[3] = {1, 0, 1};
    int16_t rep[3] = {0, 0, 0};
    if (carquet_writer_write_batch(w, 0, values, 3, def, rep) != CARQUET_OK)
        TEST_FAIL(name, "write_batch");
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) TEST_FAIL(name, err.message);
    carquet_row_batch_t* batch = NULL;
    if (carquet_batch_reader_next(br, &batch) != CARQUET_OK || !batch)
        TEST_FAIL(name, "next");

    const int32_t* offsets = NULL; int64_t num_lists = 0;
    const void* vals = NULL; const uint8_t* vvalid = NULL; int64_t nvals = 0;
    const uint8_t* list_validity = NULL;
    if (carquet_row_batch_column_list(batch, 0, &offsets, &num_lists, &vals,
                                      &vvalid, &nvals, &list_validity) != CARQUET_OK)
        TEST_FAIL(name, "column_list");
    assert(num_lists == 3 && nvals == 2);
    assert(offsets[0] == 0 && offsets[1] == 1 && offsets[2] == 1 && offsets[3] == 2);
    /* The container is REQUIRED: no list can be null, so the empty list in
     * row 1 must not be reported as one. */
    if (list_validity != NULL && !(list_validity[0] & 2))
        TEST_FAIL(name, "empty list in a REQUIRED LIST reported as null");

    carquet_row_batch_free(batch);
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_dict_nan_first_stats(void) {
    const char* name = "dict_nan_first_stats";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_dictnan");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "f", CARQUET_PHYSICAL_FLOAT, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "z", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.dictionary_encoding = CARQUET_ENCODING_RLE_DICTIONARY;
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    for (int c = 0; c < 3; c++) {
        assert(carquet_writer_set_column_encoding(w, c, CARQUET_ENCODING_RLE_DICTIONARY) == CARQUET_OK);
    }
    /* NaN first, with repeats so the dictionary is kept (not all-unique). */
    double d[6] = {NAN, 5.0, 1.0, 9.0, 5.0, 1.0};
    float f[6] = {NAN, 5.0f, 1.0f, 9.0f, 5.0f, 1.0f};
    double z[6] = {0.0, -0.0, 0.0, -0.0, 0.0, 0.0};
    assert(carquet_writer_write_batch(w, 0, d, 6, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 1, f, 6, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 2, z, 6, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_column_chunk_metadata_t cm;
    assert(carquet_reader_column_chunk_metadata(r, 0, 0, &cm) == CARQUET_OK);
    if (!cm.has_dictionary_page) TEST_FAIL(name, "column was not dictionary-encoded");

    carquet_column_statistics_t st;
    assert(carquet_reader_column_statistics(r, 0, 0, &st) == CARQUET_OK);
    if (!st.has_min_max) TEST_FAIL(name, "double: no min/max");
    double mn, mx;
    memcpy(&mn, st.min_value, 8); memcpy(&mx, st.max_value, 8);
    if (!(mn == 1.0 && mx == 9.0)) TEST_FAIL(name, "double: NaN-first dictionary chunk has wrong min/max");

    assert(carquet_reader_column_statistics(r, 0, 1, &st) == CARQUET_OK);
    if (!st.has_min_max) TEST_FAIL(name, "float: no min/max");
    float fmn, fmx;
    memcpy(&fmn, st.min_value, 4); memcpy(&fmx, st.max_value, 4);
    if (!(fmn == 1.0f && fmx == 9.0f)) TEST_FAIL(name, "float: NaN-first dictionary chunk has wrong min/max");

    /* Signed zero normalisation: min is -0.0, max is +0.0. */
    assert(carquet_reader_column_statistics(r, 0, 2, &st) == CARQUET_OK);
    if (!st.has_min_max) TEST_FAIL(name, "zero: no min/max");
    memcpy(&mn, st.min_value, 8); memcpy(&mx, st.max_value, 8);
    if (!(mn == 0.0 && signbit(mn))) TEST_FAIL(name, "zero: min is not -0.0");
    if (!(mx == 0.0 && !signbit(mx))) TEST_FAIL(name, "zero: max is not +0.0");

    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_empty_string_chunk_stats(void) {
    const char* name = "empty_string_chunk_stats";
    for (int dict = 0; dict < 2; dict++) {
        char path[512];
        carquet_test_temp_path(path, sizeof(path), dict ? "fix072_emptydict" : "fix072_emptyplain");
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* s = carquet_schema_create(&err);
        carquet_logical_type_t str = {0};
        str.id = CARQUET_LOGICAL_STRING;
        carquet_schema_add_column(s, "s", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                                  CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        opts.write_page_index = true;
        opts.dictionary_encoding = dict ? CARQUET_ENCODING_RLE_DICTIONARY : CARQUET_ENCODING_PLAIN;
        carquet_writer_t* w = open_writer(path, s, &opts, &err);
        if (!w) TEST_FAIL(name, err.message);
        assert(carquet_writer_set_column_encoding(w, 0,
            dict ? CARQUET_ENCODING_RLE_DICTIONARY : CARQUET_ENCODING_PLAIN) == CARQUET_OK);
        carquet_byte_array_t v[6] = {
            {(uint8_t*)"", 0}, {(uint8_t*)"zzz", 3}, {(uint8_t*)"m", 1},
            {NULL, 0}, {(uint8_t*)"zzz", 3}, {(uint8_t*)"m", 1} };
        assert(carquet_writer_write_batch(w, 0, v, 6, NULL, NULL) == CARQUET_OK);
        if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
        carquet_schema_free(s);

        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL(name, err.message);
        carquet_column_statistics_t st;
        assert(carquet_reader_column_statistics(r, 0, 0, &st) == CARQUET_OK);
        if (!st.has_min_max) TEST_FAIL(name, "chunk with \"\" minimum lost its statistics");
        if (st.min_value_size != 0) TEST_FAIL(name, "min is not the empty string");
        if (st.max_value_size != 3 || memcmp(st.max_value, "zzz", 3) != 0)
            TEST_FAIL(name, "max is not \"zzz\"");

        /* Page bounds must say the same. */
        carquet_column_index_t* ci = carquet_reader_get_column_index(r, 0, 0, &err);
        if (!ci) TEST_FAIL(name, "no column index");
        assert(carquet_column_index_num_pages(ci) >= 1);
        carquet_page_stats_t ps;
        assert(carquet_column_index_get_page_stats(ci, 0, &ps) == CARQUET_OK);
        if (ps.is_null_page) TEST_FAIL(name, "page flagged null");
        if (ps.max_value_size != 3 || memcmp(ps.max_value, "zzz", 3) != 0)
            TEST_FAIL(name, "page max is not \"zzz\"");
        carquet_column_index_free(ci);

        /* A predicate above "m" must still reach the page holding "zzz". */
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        carquet_filter_clause_t clause = {0};
        clause.column_index = 0;
        clause.op = CARQUET_FILTER_GT;
        clause.value = "y";
        clause.value_size = 1;
        assert(carquet_batch_reader_set_page_filter(br, &clause, 1) == CARQUET_OK);
        int64_t total = 0;
        carquet_row_batch_t* b = NULL;
        while (carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
            total += carquet_row_batch_num_rows(b);
            carquet_row_batch_free(b);
            b = NULL;
        }
        if (total != 6) TEST_FAIL(name, "page with \"zzz\" was pruned");
        carquet_batch_reader_free(br);
        carquet_reader_close(r);
        remove(path);
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_dict_fallback_not_sticky(void) {
    const char* name = "dict_fallback_not_sticky";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_sticky");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.dictionary_encoding = CARQUET_ENCODING_RLE_DICTIONARY;
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    assert(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY) == CARQUET_OK);

    /* Row group 0: all unique -> dictionary gives nothing, falls back to PLAIN. */
    int64_t uniq[64];
    for (int i = 0; i < 64; i++) uniq[i] = i;
    assert(carquet_writer_write_batch(w, 0, uniq, 64, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_new_row_group(w) == CARQUET_OK);
    /* Row group 1: 2 distinct values -> must get a dictionary page again. */
    int64_t rep[64];
    for (int i = 0; i < 64; i++) rep[i] = i & 1;
    assert(carquet_writer_write_batch(w, 0, rep, 64, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    assert(carquet_reader_num_row_groups(r) == 2);
    carquet_column_chunk_metadata_t cm0, cm1;
    assert(carquet_reader_column_chunk_metadata(r, 0, 0, &cm0) == CARQUET_OK);
    assert(carquet_reader_column_chunk_metadata(r, 1, 0, &cm1) == CARQUET_OK);
    if (cm0.has_dictionary_page) TEST_FAIL(name, "all-unique row group kept a dictionary");
    if (!cm1.has_dictionary_page)
        TEST_FAIL(name, "dictionary fallback in row group 0 stuck to row group 1");

    /* Values round-trip in both groups. */
    for (int rg = 0; rg < 2; rg++) {
        carquet_column_reader_t* c = carquet_reader_get_column(r, rg, 0, &err);
        int64_t out[64];
        assert(carquet_column_read_batch(c, out, 64, NULL, NULL) == 64);
        for (int i = 0; i < 64; i++) assert(out[i] == (rg == 0 ? uniq[i] : rep[i]));
        carquet_column_reader_free(c);
    }
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_list_offset_index_rows(void) {
    const char* name = "list_offset_index_rows";
    for (int v2 = 0; v2 < 2; v2++) {
        char path[512];
        carquet_test_temp_path(path, sizeof(path), v2 ? "fix072_listoi2" : "fix072_listoi1");
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* s = carquet_schema_create(&err);
        assert(carquet_schema_add_list(s, "l", CARQUET_PHYSICAL_INT32, NULL,
                                       CARQUET_REPETITION_OPTIONAL, 0, 0) >= 0);
        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        opts.write_page_index = true;
        opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
        opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
        /* page_size / 4 = 1000 values, clamped to the 1024 chunk minimum:
         * the first size check happens after 1024 levels, which is not a
         * multiple of 3. The page must extend to the next record start, so
         * page 0 holds 1026 levels = 342 rows. */
        opts.page_size = 4000;
        opts.data_page_version = v2 ? 2 : 1;
        carquet_writer_t* w = open_writer(path, s, &opts, &err);
        if (!w) TEST_FAIL(name, err.message);
        assert(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN) == CARQUET_OK);

        const int rows = 2000;
        int32_t* vals = malloc(sizeof(int32_t) * rows * 3);
        int16_t* def = malloc(sizeof(int16_t) * rows * 3);
        int16_t* rep = malloc(sizeof(int16_t) * rows * 3);
        for (int i = 0; i < rows * 3; i++) {
            vals[i] = i; def[i] = 3; rep[i] = (i % 3) ? 1 : 0;  /* OPTIONAL list of REQUIRED ints: max_def 2? see below */
        }
        /* max_def for optional list<required int> is 2 (list present, element
         * present); clamp the levels to the schema's maximum. */
        int16_t max_def = carquet_schema_max_def_level(s, 0);
        for (int i = 0; i < rows * 3; i++) def[i] = max_def;
        if (carquet_writer_write_batch(w, 0, vals, rows * 3, def, rep) != CARQUET_OK)
            TEST_FAIL(name, "write_batch");
        if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
        carquet_schema_free(s);
        free(vals); free(def); free(rep);

        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL(name, err.message);
        carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, 0, &err);
        if (!oi) TEST_FAIL(name, "no offset index");
        int32_t n = carquet_offset_index_num_pages(oi);
        if (n < 2) TEST_FAIL(name, "expected several pages");
        carquet_page_location_t loc;
        assert(carquet_offset_index_get_page_location(oi, 0, &loc) == CARQUET_OK);
        assert(loc.first_row_index == 0);
        assert(carquet_offset_index_get_page_location(oi, 1, &loc) == CARQUET_OK);
        if (loc.first_row_index != 342) {
            printf("  page 1 first_row_index = %lld (expected 342)\n", (long long)loc.first_row_index);
            TEST_FAIL(name, "first_row_index counts levels or page cut mid-record");
        }
        int64_t prev = -1;
        for (int32_t p = 0; p < n; p++) {
            assert(carquet_offset_index_get_page_location(oi, p, &loc) == CARQUET_OK);
            assert(loc.first_row_index > prev && loc.first_row_index < rows);
            prev = loc.first_row_index;
        }
        carquet_offset_index_free(oi);

        /* Data still round-trips through the nested batch reader. */
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        int64_t total_lists = 0, total_vals = 0;
        carquet_row_batch_t* b = NULL;
        while (carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
            const int32_t* offsets; int64_t nl; const void* vv; const uint8_t* vva;
            int64_t nv; const uint8_t* lv;
            assert(carquet_row_batch_column_list(b, 0, &offsets, &nl, &vv, &vva, &nv, &lv) == CARQUET_OK);
            total_lists += nl; total_vals += nv;
            carquet_row_batch_free(b); b = NULL;
        }
        assert(total_lists == rows && total_vals == rows * 3);
        carquet_batch_reader_free(br);
        carquet_reader_close(r);
        remove(path);
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_column_index_withheld(void) {
    const char* name = "column_index_withheld";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_nostats_ci");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "w", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.write_statistics = false;   /* no page bounds for any column */
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    int32_t v[16];
    for (int i = 0; i < 16; i++) v[i] = i;
    assert(carquet_writer_write_batch(w, 0, v, 16, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 1, v, 16, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_column_chunk_metadata_t cm;
    assert(carquet_reader_column_chunk_metadata(r, 0, 0, &cm) == CARQUET_OK);
    if (cm.has_column_index) TEST_FAIL(name, "ColumnIndex written for a column without statistics");
    if (!cm.has_offset_index) TEST_FAIL(name, "OffsetIndex missing");
    assert(carquet_reader_column_chunk_metadata(r, 0, 1, &cm) == CARQUET_OK);
    if (cm.has_column_index) TEST_FAIL(name, "ColumnIndex written for a column without statistics");
    if (!cm.has_offset_index) TEST_FAIL(name, "OffsetIndex missing");

    /* A page filter on the stats-less column keeps every row (nothing can
     * be proven absent) rather than failing or pruning. */
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    int32_t target = 7;
    carquet_filter_clause_t clause = {0};
    clause.column_index = 0; clause.op = CARQUET_FILTER_EQ;
    clause.value = &target; clause.value_size = 4;
    assert(carquet_batch_reader_set_page_filter(br, &clause, 1) == CARQUET_OK);
    int64_t total = 0;
    carquet_row_batch_t* b = NULL;
    carquet_status_t st;
    while ((st = carquet_batch_reader_next(br, &b)) == CARQUET_OK && b) {
        total += carquet_row_batch_num_rows(b);
        carquet_row_batch_free(b); b = NULL;
    }
    if (st != CARQUET_ERROR_END_OF_DATA && st != CARQUET_OK) TEST_FAIL(name, "page filter failed");
    if (total != 16) TEST_FAIL(name, "rows dropped without page bounds");
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_byte_array_page_size(void) {
    const char* name = "byte_array_page_size";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_bapage");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "s", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 64 * 1024;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    assert(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN) == CARQUET_OK);

    const int n = 20000;
    char* blob = malloc((size_t)n * 50);
    carquet_byte_array_t* v = malloc(sizeof(carquet_byte_array_t) * n);
    for (int i = 0; i < n; i++) {
        memset(blob + (size_t)i * 50, 'a' + (i % 26), 50);
        v[i].data = (uint8_t*)blob + (size_t)i * 50;
        v[i].length = 50;
    }
    /* One call for ~1 MB of strings: previously one page. */
    assert(carquet_writer_write_batch(w, 0, v, n, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);
    free(blob); free(v);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, 0, &err);
    if (!oi) TEST_FAIL(name, "no offset index");
    int32_t pages = carquet_offset_index_num_pages(oi);
    carquet_offset_index_free(oi);
    if (pages < 8) {
        printf("  %d page(s) for 1 MB of strings at page_size 64 KB\n", pages);
        TEST_FAIL(name, "page_size ignored for BYTE_ARRAY");
    }
    assert(carquet_reader_num_rows(r) == n);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_rle_bitpacked_grouping(void) {
    const char* name = "rle_bitpacked_grouping";
    /* 4096 alternating bits: no run of 8, so everything is bit-packed. One
     * header per 8 values costs 1024 bytes; grouping 63 × 8 values per
     * header costs 512 + 9 = 521. */
    uint32_t in[4096];
    for (int i = 0; i < 4096; i++) in[i] = (uint32_t)(i & 1);
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    assert(carquet_rle_encode_all(in, 4096, 1, &buf) == CARQUET_OK);
    if (buf.size > 530) {
        printf("  encoded size %zu\n", buf.size);
        TEST_FAIL(name, "bit-packed groups are not sharing headers");
    }
    uint32_t out[4096];
    int64_t got = carquet_rle_decode_all(buf.data, buf.size, 1, out, 4096);
    assert(got == 4096 && memcmp(in, out, sizeof(in)) == 0);
    carquet_buffer_destroy(&buf);

    /* Mixed streams: runs interleaved with bit-packed stretches, every
     * length and width, including stretches that end in a partial group. */
    uint32_t seed = 7;
    for (int bw = 1; bw <= 32; bw += (bw < 8 ? 1 : 7)) {
        for (int len = 1; len <= 1200; len += 37) {
            uint32_t* a = malloc(sizeof(uint32_t) * len);
            uint32_t* b = malloc(sizeof(uint32_t) * len);
            uint32_t mask = bw == 32 ? 0xFFFFFFFFu : ((1u << bw) - 1);
            int i = 0;
            while (i < len) {
                seed = seed * 1103515245u + 12345u;
                int run = (int)((seed >> 16) % 20) + 1;
                seed = seed * 1103515245u + 12345u;
                uint32_t v = (seed >> 8) & mask;
                for (int k = 0; k < run && i < len; k++, i++) {
                    if (run >= 10) a[i] = v;
                    else { seed = seed * 1103515245u + 12345u; a[i] = (seed >> 8) & mask; }
                }
            }
            carquet_buffer_init(&buf);
            assert(carquet_rle_encode_all(a, len, bw, &buf) == CARQUET_OK);
            got = carquet_rle_decode_all(buf.data, buf.size, bw, b, len);
            if (got != len || memcmp(a, b, sizeof(uint32_t) * len) != 0) {
                printf("  bw=%d len=%d got=%lld\n", bw, len, (long long)got);
                TEST_FAIL(name, "mixed RLE/bit-packed stream does not round-trip");
            }
            carquet_buffer_destroy(&buf);
            free(a); free(b);
        }
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_batch_size_zero_config(void) {
    const char* name = "batch_size_zero_config";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_bs0");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_t* w = open_writer(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);
    int32_t v[1000];
    for (int i = 0; i < 1000; i++) v[i] = i;
    assert(carquet_writer_write_batch(w, 0, v, 1000, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_batch_reader_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));   /* batch_size == 0 */
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) TEST_FAIL(name, err.message);
    int64_t total = 0; int iters = 0;
    carquet_row_batch_t* b = NULL;
    while (carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
        total += carquet_row_batch_num_rows(b);
        carquet_row_batch_free(b); b = NULL;
        if (++iters > 10000) TEST_FAIL(name, "batch_size 0 loops forever on empty batches");
    }
    if (total != 1000) TEST_FAIL(name, "row count mismatch");
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_dispatch_match_copy(void) {
    const char* name = "dispatch_match_copy";
    uint8_t buf[256], ref[256];
    for (size_t offset = 1; offset <= 20; offset++) {
        for (size_t len = 0; len <= 100; len++) {
            for (size_t i = 0; i < sizeof(buf); i++) buf[i] = ref[i] = (uint8_t)(i * 37 + 11);
            /* dst starts at `offset`; src is the preceding `offset` bytes */
            carquet_dispatch_match_copy(buf + offset, buf, len, offset);
            for (size_t i = 0; i < len; i++) ref[offset + i] = ref[i];
            if (memcmp(buf, ref, sizeof(buf)) != 0) {
                printf("  offset=%zu len=%zu\n", offset, len);
                TEST_FAIL(name, "match_copy differs from the overlapping byte copy");
            }
        }
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_bss_large_tiled_decode(void) {
    const char* name = "bss_large_tiled_decode";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_bss");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "f", CARQUET_PHYSICAL_FLOAT, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.page_size = 8 * 1024 * 1024;  /* one page above both tile sizes */
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    assert(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_BYTE_STREAM_SPLIT) == CARQUET_OK);
    assert(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_BYTE_STREAM_SPLIT) == CARQUET_OK);
    const int n = 150000;  /* > 65536 floats and > 32768 doubles: tiled path */
    float* f = malloc(sizeof(float) * n);
    double* d = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) { f[i] = (float)i * 0.5f; d[i] = (double)i * 0.25; }
    assert(carquet_writer_write_batch(w, 0, f, n, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 1, d, n, NULL, NULL) == CARQUET_OK);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    float* fo = malloc(sizeof(float) * n);
    double* dout = malloc(sizeof(double) * n);
    carquet_column_reader_t* c0 = carquet_reader_get_column(r, 0, 0, &err);
    carquet_column_reader_t* c1 = carquet_reader_get_column(r, 0, 1, &err);
    assert(carquet_column_read_batch(c0, fo, n, NULL, NULL) == n);
    assert(carquet_column_read_batch(c1, dout, n, NULL, NULL) == n);
    assert(memcmp(f, fo, sizeof(float) * n) == 0);
    assert(memcmp(d, dout, sizeof(double) * n) == 0);
    carquet_column_reader_free(c0);
    carquet_column_reader_free(c1);
    carquet_reader_close(r);
    free(f); free(d); free(fo); free(dout);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* BROTLI (codec 4) round trip through the built-in codec; skipped when the
 * library was built without libbrotli. */
static int test_brotli_roundtrip(void) {
    const char* name = "brotli_roundtrip";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix072_brotli");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t str = {0};
    str.id = CARQUET_LOGICAL_STRING;
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "t", CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_BROTLI;
    opts.compression_level = 5;
    carquet_writer_t* w = open_writer(path, s, &opts, &err);
    if (!w) TEST_FAIL(name, err.message);
    const int n = 20000;
    int64_t* v = malloc(sizeof(int64_t) * n);
    carquet_byte_array_t* t = malloc(sizeof(carquet_byte_array_t) * n);
    int16_t* def = malloc(sizeof(int16_t) * n);
    static const char* words[] = {"alpha", "beta", "gamma", "delta", ""};
    int present = 0;
    for (int i = 0; i < n; i++) {
        v[i] = (int64_t)i * 31 % 977;
        def[i] = (i % 7) != 3;
        if (def[i]) { const char* wd = words[i % 5]; t[present].data = (uint8_t*)wd; t[present].length = (int32_t)strlen(wd); present++; }
    }
    assert(carquet_writer_write_batch(w, 0, v, n, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 1, t, n, def, NULL) == CARQUET_OK);
    carquet_status_t st = carquet_writer_close(w);
    if (st == CARQUET_ERROR_UNSUPPORTED_CODEC) {
        free(v); free(t); free(def); carquet_schema_free(s); remove(path);
        printf("[SKIP] %s (built without libbrotli)\n", name);
        return 0;
    }
    if (st != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_column_chunk_metadata_t cm;
    assert(carquet_reader_column_chunk_metadata(r, 0, 0, &cm) == CARQUET_OK);
    if (cm.codec != CARQUET_COMPRESSION_BROTLI) TEST_FAIL(name, "codec not recorded as BROTLI");
    if (cm.total_compressed_size >= cm.total_uncompressed_size) TEST_FAIL(name, "no compression happened");
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    carquet_row_batch_t* b = NULL;
    int64_t row = 0; int seen = 0;
    while (carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
        const void* d; const uint8_t* nb; int64_t cnt;
        assert(carquet_row_batch_column(b, 0, &d, &nb, &cnt) == CARQUET_OK);
        for (int64_t i = 0; i < cnt; i++) assert(((const int64_t*)d)[i] == v[row + i]);
        const void* td; const uint8_t* tv; int64_t tcnt;
        assert(carquet_row_batch_column(b, 1, &td, &tv, &tcnt) == CARQUET_OK);
        /* The batch reader is row-aligned for nullable columns: one slot per
         * logical row, null slots zero-filled, the bitmap indexed by row. */
        const carquet_byte_array_t* ba = (const carquet_byte_array_t*)td;
        for (int64_t i = 0; i < cnt; i++) {
            bool pres = !tv || (tv[i >> 3] & (1u << (i & 7)));
            assert(pres == (def[row + i] != 0));
            if (pres) {
                const char* wd = words[(row + i) % 5];
                if (ba[i].length != (int32_t)strlen(wd)) {
                    printf("  row %lld: got length %d, expected %zu\n",
                           (long long)(row + i), ba[i].length, strlen(wd));
                    TEST_FAIL(name, "string length mismatch after BROTLI round trip");
                }
                assert(ba[i].length == 0 || memcmp(ba[i].data, wd, (size_t)ba[i].length) == 0);
                seen++;
            }
        }
        row += cnt;
        carquet_row_batch_free(b); b = NULL;
    }
    assert(row == n && seen == present);
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    free(v); free(t); free(def);
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* An encrypted-footer file ("PARE" trailing magic) is reported as
 * NOT_IMPLEMENTED with an explicit message, not as a corrupt file. */
static int test_encrypted_file_detected(void) {
    const char* name = "encrypted_file_detected";
    uint8_t buf[64];
    memset(buf, 0, sizeof(buf));
    memcpy(buf, "PAR1", 4);
    memcpy(buf + sizeof(buf) - 4, "PARE", 4);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sizeof(buf), NULL, &err);
    if (r) { carquet_reader_close(r); TEST_FAIL(name, "opened an encrypted file"); }
    if (err.code != CARQUET_ERROR_NOT_IMPLEMENTED) TEST_FAIL(name, err.message);
    if (!strstr(err.message, "encryption")) TEST_FAIL(name, "message does not mention encryption");
    TEST_PASS(name);
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_required_list_empty_not_null();
    rc |= test_dict_nan_first_stats();
    rc |= test_empty_string_chunk_stats();
    rc |= test_dict_fallback_not_sticky();
    rc |= test_list_offset_index_rows();
    rc |= test_column_index_withheld();
    rc |= test_byte_array_page_size();
    rc |= test_rle_bitpacked_grouping();
    rc |= test_batch_size_zero_config();
    rc |= test_dispatch_match_copy();
    rc |= test_bss_large_tiled_decode();
    rc |= test_brotli_roundtrip();
    rc |= test_encrypted_file_detected();
    if (rc == 0) printf("\nAll 0.7.2 regression tests passed.\n");
    return rc;
}
