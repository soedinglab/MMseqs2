/**
 * @file test_fixes_0_7_4.c
 * @brief Regression tests for the bugs fixed after v0.7.3.
 *
 * One test per fix, each failing on the v0.7.3 code:
 *  - dense_cursor_after_skip      a skip into the next page reused the dense
 *                                 cursor of the previous page (OPTIONAL
 *                                 column returned shifted values)
 *  - list_filter_null_count       page filter treated a LIST leaf as all-null
 *                                 when null_count (level slots) >= rows
 *  - column_index_list_lengths    ColumnIndex num_pages was clamped to the
 *                                 optional null_counts list, dropping pages
 *  - bloom_signed_zero            bloom pruning missed -0.0 for EQ +0.0
 *  - deprecated_stats_unsigned    deprecated signed min/max trusted for UINT
 *  - append_shorter_footer        append left the old footer tail in place
 *  - schema_preorder              adding to a closed group corrupted the schema
 *  - geo_stats_per_row_group      GeospatialStatistics accumulated across
 *                                 row groups
 *  - rle_encoder_oom              RLE encoder returned OK on a failed append
 *  - max_rows_per_page_one_call   one write_batch call ignored the row cap
 *  - arrow_dict_null_non_nullable null dictionary slot on a non-nullable field
 *                                 desynchronised levels and values (heap read)
 *  - arrow_schema_wide_struct     ARROW:schema of a struct with > 64 fields
 *                                 overflowed a stack array
 *  - preserved_dict_index_range   preserve_dictionaries handed out-of-range
 *                                 indices to the caller (V1 and V2 pages)
 *  - v2_empty_compressed_section  GZIP V2 all-null page with a 0-byte data
 *                                 section could not be read
 */

#include <carquet/carquet.h>
#include "core/buffer.h"
#include "encoding/rle.h"
#include "reader/reader_internal.h"
#include "thrift/thrift_encode.h"
#include "test_helpers.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern carquet_column_index_t* carquet_column_index_parse(const uint8_t* data, size_t size);

/* Count the rows a batch reader returns with one page-filter clause. */
static int64_t filtered_rows(const char* path, const carquet_filter_clause_t* clause,
                             carquet_status_t* status_out) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { *status_out = CARQUET_ERROR_FILE_OPEN; return -1; }
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); *status_out = CARQUET_ERROR_INVALID_STATE; return -1; }
    carquet_status_t st = carquet_batch_reader_set_page_filter(br, clause, 1);
    int64_t rows = 0;
    carquet_row_batch_t* b = NULL;
    while (st == CARQUET_OK && (st = carquet_batch_reader_next(br, &b)) == CARQUET_OK && b) {
        rows += carquet_row_batch_num_rows(b);
        carquet_row_batch_free(b);
        b = NULL;
    }
    if (st == CARQUET_ERROR_END_OF_DATA) st = CARQUET_OK;
    *status_out = st;
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return rows;
}

/* ------------------------------------------------------------------------ */
static int test_dense_cursor_after_skip(void) {
    const char* name = "dense_cursor_after_skip";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_cursor");
    carquet_error_t err = CARQUET_ERROR_INIT;

    enum { N = 3000 };
    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.max_rows_per_page = 1000;
    o.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    carquet_writer_t* w = carquet_writer_create(path, s, &o, &err);
    if (!w) TEST_FAIL(name, err.message);

    /* Rows 0..299 null, row i >= 300 holds i. Values are dense. */
    static int32_t vals[N];
    static int16_t def[N];
    int nv = 0;
    for (int i = 0; i < N; i++) {
        def[i] = i < 300 ? 0 : 1;
        if (def[i]) vals[nv++] = i;
    }
    /* One call per page, so the layout does not depend on the writer's
     * chunking (see max_rows_per_page_one_call). */
    int off = 0;
    for (int p = 0; p < N / 1000; p++) {
        int present = 0;
        for (int i = p * 1000; i < (p + 1) * 1000; i++) present += def[i];
        if (carquet_writer_write_batch(w, 0, vals + off, 1000, def + p * 1000, NULL)
                != CARQUET_OK)
            TEST_FAIL(name, "write");
        off += present;
    }
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    if (!c) TEST_FAIL(name, err.message);

    int32_t out[300];
    int16_t dl[300];
    /* 300 null slots of page 0: leaves the dense cursor at slot 300. */
    int64_t got = carquet_column_read_batch(c, out, 300, dl, NULL);
    int64_t skipped = carquet_column_skip(c, 1000);   /* -> row 1300, page 1 */
    int64_t got2 = carquet_column_read_batch(c, out, 5, dl, NULL);
    int ok = got == 300 && skipped == 1000 && got2 == 5;
    for (int i = 0; ok && i < 5; i++) ok = dl[i] == 1 && out[i] == 1300 + i;
    if (!ok) {
        printf("  got=%lld skipped=%lld got2=%lld first=%d (want 1300)\n",
               (long long)got, (long long)skipped, (long long)got2, out[0]);
    }
    carquet_column_reader_free(c);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "values after skip shifted");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_list_filter_null_count(void) {
    const char* name = "list_filter_null_count";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_listfilter");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_list(s, "l", CARQUET_PHYSICAL_INT32, NULL,
                                CARQUET_REPETITION_OPTIONAL, 0, 0) < 0)
        TEST_FAIL(name, "schema");
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.write_page_index = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &o, &err);
    if (!w) TEST_FAIL(name, err.message);

    /* 10 rows [r, null, null, null]: 30 null slots for 10 rows. */
    int32_t vals[10];
    int16_t def[40], rep[40];
    for (int r = 0; r < 10; r++) {
        vals[r] = r;
        for (int k = 0; k < 4; k++) {
            rep[r * 4 + k] = k == 0 ? 0 : 1;
            def[r * 4 + k] = k == 0 ? 3 : 2;
        }
    }
    if (carquet_writer_write_batch(w, 0, vals, 40, def, rep) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(s);

    int32_t five = 5;
    carquet_filter_clause_t eq;
    memset(&eq, 0, sizeof(eq));
    eq.column_index = 0;
    eq.op = CARQUET_FILTER_EQ;
    eq.value = &five;
    eq.value_size = 4;
    carquet_filter_clause_t nn;
    memset(&nn, 0, sizeof(nn));
    nn.column_index = 0;
    nn.op = CARQUET_FILTER_IS_NOT_NULL;

    carquet_status_t st1, st2;
    int64_t rows_eq = filtered_rows(path, &eq, &st1);
    int64_t rows_nn = filtered_rows(path, &nn, &st2);
    carquet_test_cleanup(path);
    /* The filter works at page granularity: the single page must survive. */
    if (st1 != CARQUET_OK || st2 != CARQUET_OK || rows_eq != 10 || rows_nn != 10) {
        printf("  EQ: st=%d rows=%lld  IS_NOT_NULL: st=%d rows=%lld (want 10)\n",
               st1, (long long)rows_eq, st2, (long long)rows_nn);
        TEST_FAIL(name, "row group pruned on a slot null count");
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static carquet_column_index_t* parse_ci(int with_null_pages, int n_null_counts) {
    carquet_buffer_t b;
    carquet_buffer_init(&b);
    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &b);
    int32_t mn[3] = {1, 10, 20}, mx[3] = {5, 15, 25};
    thrift_write_struct_begin(&enc);
    if (with_null_pages) {
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 1);
        thrift_write_list_begin(&enc, THRIFT_TYPE_TRUE, 3);
        for (int i = 0; i < 3; i++) thrift_write_byte(&enc, 2); /* false */
    }
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 2);
    thrift_write_list_begin(&enc, THRIFT_TYPE_BINARY, 3);
    for (int i = 0; i < 3; i++) thrift_write_binary(&enc, (const uint8_t*)&mn[i], 4);
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 3);
    thrift_write_list_begin(&enc, THRIFT_TYPE_BINARY, 3);
    for (int i = 0; i < 3; i++) thrift_write_binary(&enc, (const uint8_t*)&mx[i], 4);
    thrift_write_field_header(&enc, THRIFT_TYPE_I32, 4);
    thrift_write_i32(&enc, 1);
    if (n_null_counts > 0) {
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 5);
        thrift_write_list_begin(&enc, THRIFT_TYPE_I64, n_null_counts);
        for (int i = 0; i < n_null_counts; i++) thrift_write_i64(&enc, 0);
    }
    thrift_write_struct_end(&enc);
    carquet_column_index_t* ci =
        carquet_column_index_parse(carquet_buffer_data(&b), carquet_buffer_size(&b));
    carquet_buffer_destroy(&b);
    return ci;
}

static int test_column_index_list_lengths(void) {
    const char* name = "column_index_list_lengths";
    int ok = 1;

    /* 3 pages, 1 null_count: the optional list is dropped, 3 pages remain. */
    carquet_column_index_t* ci = parse_ci(1, 1);
    carquet_page_stats_t ps;
    if (!ci || carquet_column_index_num_pages(ci) != 3) {
        printf("  short null_counts: num_pages=%d (want 3)\n",
               ci ? carquet_column_index_num_pages(ci) : -1);
        ok = 0;
    } else if (carquet_column_index_get_page_stats(ci, 2, &ps) != CARQUET_OK ||
               ps.has_null_count || !ps.min_value) {
        printf("  page 2 stats wrong\n");
        ok = 0;
    }
    if (ci) carquet_column_index_free(ci);

    /* Consistent lists still parse. */
    ci = parse_ci(1, 3);
    if (!ci || carquet_column_index_num_pages(ci) != 3) ok = 0;
    if (ci) carquet_column_index_free(ci);

    /* No null_pages (required): the index cannot be mapped, reject it. It
     * used to parse with 0 pages, which dropped the whole row group. */
    ci = parse_ci(0, 3);
    if (ci) {
        printf("  ColumnIndex without null_pages accepted (%d pages)\n",
               carquet_column_index_num_pages(ci));
        carquet_column_index_free(ci);
        ok = 0;
    }
    if (!ok) TEST_FAIL(name, "ColumnIndex page count");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_bloom_signed_zero(void) {
    const char* name = "bloom_signed_zero";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_bloomzero");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_DOUBLE, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(s, "f", CARQUET_PHYSICAL_FLOAT, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.write_page_index = true;
    o.write_bloom_filters = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &o, &err);
    if (!w) TEST_FAIL(name, err.message);
    double d[3] = {-0.0, 1.0, 2.0};
    float f[3] = {-0.0f, 1.0f, 2.0f};
    if (carquet_writer_write_batch(w, 0, d, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, f, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(s);

    double dz = 0.0;
    float fz = 0.0f;
    carquet_filter_clause_t c;
    memset(&c, 0, sizeof(c));
    c.op = CARQUET_FILTER_EQ;
    c.column_index = 0;
    c.value = &dz;
    c.value_size = 8;
    carquet_status_t st1, st2;
    int64_t rows_d = filtered_rows(path, &c, &st1);
    c.column_index = 1;
    c.value = &fz;
    c.value_size = 4;
    int64_t rows_f = filtered_rows(path, &c, &st2);
    carquet_test_cleanup(path);
    if (st1 != CARQUET_OK || st2 != CARQUET_OK || rows_d != 3 || rows_f != 3) {
        printf("  double: st=%d rows=%lld  float: st=%d rows=%lld (want 3)\n",
               st1, (long long)rows_d, st2, (long long)rows_f);
        TEST_FAIL(name, "EQ +0.0 pruned a row group holding -0.0");
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
/* Turn the footer min/max of chunk (0, col) into the deprecated fields, the
 * way pre-1.10 parquet-mr wrote them (signed order). */
static void make_deprecated(carquet_reader_t* r, int col, const void* min, const void* max) {
    parquet_statistics_t* ps = &r->metadata.row_groups[0].columns[col].metadata.statistics;
    ps->has_min_value = false;
    ps->has_max_value = false;
    memcpy(ps->min_value, min, 4);
    memcpy(ps->max_value, max, 4);
    ps->min_deprecated = ps->min_value;
    ps->min_deprecated_len = 4;
    ps->max_deprecated = ps->max_value;
    ps->max_deprecated_len = 4;
}

static int test_deprecated_stats_unsigned(void) {
    const char* name = "deprecated_stats_unsigned";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_depstats");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t u32;
    memset(&u32, 0, sizeof(u32));
    u32.id = CARQUET_LOGICAL_INTEGER;
    u32.params.integer.bit_width = 32;
    u32.params.integer.is_signed = false;
    if (carquet_schema_add_column(s, "u", CARQUET_PHYSICAL_INT32, &u32,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(s, "i", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);
    uint32_t u[2] = {1u, 0x80000000u};
    int32_t iv[2] = {-5, 7};
    if (carquet_writer_write_batch(w, 0, u, 2, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, iv, 2, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    /* Signed order: min = 0x80000000 (negative), max = 1. */
    uint32_t smin = 0x80000000u, smax = 1u;
    make_deprecated(r, 0, &smin, &smax);
    int32_t imin = -5, imax = 7;
    make_deprecated(r, 1, &imin, &imax);

    carquet_column_statistics_t cs;
    int ok = carquet_reader_column_statistics(r, 0, 0, &cs) == CARQUET_OK && !cs.has_min_max;
    if (!ok) printf("  UINT_32 deprecated min/max trusted\n");

    /* Pruning must keep the group: 1 is in the column. */
    bool might = false;
    uint32_t one = 1u;
    if (carquet_reader_row_group_matches(r, 0, 0, CARQUET_COMPARE_EQ, &one, 4, &might)
            != CARQUET_OK || !might) {
        printf("  EQ 1 pruned the row group\n");
        ok = 0;
    }

    /* Signed INT32 keeps using the deprecated pair. */
    if (carquet_reader_column_statistics(r, 0, 1, &cs) != CARQUET_OK || !cs.has_min_max) {
        printf("  signed INT32 deprecated min/max dropped\n");
        ok = 0;
    }
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "deprecated statistics");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_append_shorter_footer(void) {
    const char* name = "append_shorter_footer";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_append");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");

    /* Footer parts that append does not restore byte for byte (empty-string
     * minimum, encoding stats) make the rewritten footer shorter. */
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);
    carquet_byte_array_t v[2] = {{(uint8_t*)"", 0}, {(uint8_t*)"zz", 2}};
    for (int rg = 0; rg < 8; rg++) {
        if (carquet_writer_write_batch(w, 0, v, 2, NULL, NULL) != CARQUET_OK ||
            carquet_writer_new_row_group(w) != CARQUET_OK)
            TEST_FAIL(name, "write");
    }
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");

    w = carquet_writer_open_append(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);
    carquet_byte_array_t one = {(uint8_t*)"a", 1};
    if (carquet_writer_write_batch(w, 0, &one, 1, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "append");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    int ok = r && carquet_reader_num_rows(r) == 17;
    if (!ok) printf("  reopen: %s\n", r ? "wrong row count" : err.message);
    if (r) carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "stale footer tail after append");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_schema_preorder(void) {
    const char* name = "schema_preorder";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    int32_t a = carquet_schema_add_group(s, "a", CARQUET_REPETITION_REQUIRED, 0);
    int32_t b = carquet_schema_add_group(s, "b", CARQUET_REPETITION_REQUIRED, 0);
    int ok = a > 0 && b > 0;
    /* a is closed once b follows it: x would serialise under b. */
    ok = ok && carquet_schema_add_column(s, "x", CARQUET_PHYSICAL_INT32, NULL,
                                         CARQUET_REPETITION_REQUIRED, 0, a)
               == CARQUET_ERROR_INVALID_ARGUMENT;
    ok = ok && carquet_schema_add_group(s, "g", CARQUET_REPETITION_REQUIRED, a) == -1;
    /* Open groups (the last element and its ancestors, root) still accept. */
    ok = ok && carquet_schema_add_column(s, "y", CARQUET_PHYSICAL_INT32, NULL,
                                         CARQUET_REPETITION_REQUIRED, 0, b) == CARQUET_OK;
    ok = ok && carquet_schema_add_column(s, "z", CARQUET_PHYSICAL_INT32, NULL,
                                         CARQUET_REPETITION_REQUIRED, 0, b) == CARQUET_OK;
    ok = ok && carquet_schema_add_list(s, "l", CARQUET_PHYSICAL_INT32, NULL,
                                       CARQUET_REPETITION_OPTIONAL, 0, 0) >= 0;
    ok = ok && carquet_schema_add_column(s, "t", CARQUET_PHYSICAL_INT32, NULL,
                                         CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
    carquet_schema_free(s);
    if (!ok) TEST_FAIL(name, "child of a closed group accepted");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static void wkb_point(uint8_t out[21], double x, double y) {
    out[0] = 1; /* little endian */
    out[1] = 1; out[2] = 0; out[3] = 0; out[4] = 0; /* Point */
    uint8_t t[8];
    memcpy(t, &x, 8);
    for (int i = 0; i < 8; i++) out[5 + i] = t[i];
    memcpy(t, &y, 8);
    for (int i = 0; i < 8; i++) out[13 + i] = t[i];
}

static int test_geo_stats_per_row_group(void) {
    const char* name = "geo_stats_per_row_group";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_geo");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t geo;
    memset(&geo, 0, sizeof(geo));
    geo.id = CARQUET_LOGICAL_GEOMETRY;
    if (carquet_schema_add_column(s, "g", CARQUET_PHYSICAL_BYTE_ARRAY, &geo,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);
    uint8_t p0[21], p1[21];
    wkb_point(p0, 1.0, 1.0);
    wkb_point(p1, 10.0, 20.0);
    carquet_byte_array_t a = {p0, 21}, b = {p1, 21};
    if (carquet_writer_write_batch(w, 0, &a, 1, NULL, NULL) != CARQUET_OK ||
        carquet_writer_new_row_group(w) != CARQUET_OK ||
        carquet_writer_write_batch(w, 0, &b, 1, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_geospatial_statistics_t g0, g1;
    int ok = carquet_reader_num_row_groups(r) == 2 &&
             carquet_reader_geospatial_statistics(r, 0, 0, &g0) == CARQUET_OK &&
             carquet_reader_geospatial_statistics(r, 1, 0, &g1) == CARQUET_OK &&
             g0.has_bbox && g0.xmin == 1.0 && g0.xmax == 1.0 &&
             g1.has_bbox && g1.xmin == 10.0 && g1.xmax == 10.0 &&
             g1.ymin == 20.0 && g1.ymax == 20.0;
    if (!ok) printf("  row group 1 bbox x=[%g, %g] (want [10, 10])\n", g1.xmin, g1.xmax);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "geospatial statistics leaked across row groups");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int g_fail_big_realloc = 0;
static void* t_malloc(size_t n, void* ctx) { (void)ctx; return malloc(n); }
static void* t_realloc(void* p, size_t n, void* ctx) {
    (void)ctx;
    if (g_fail_big_realloc && n >= 8192) return NULL;
    return realloc(p, n);
}
static void t_free(void* p, void* ctx) { (void)ctx; free(p); }

static int test_rle_encoder_oom(void) {
    const char* name = "rle_encoder_oom";
    enum { N = 20000 };
    static int16_t levels[N];
    for (int i = 0; i < N; i++) levels[i] = (int16_t)((i * 7) % 200); /* bit-packed */

    carquet_allocator_t a = {t_malloc, t_realloc, t_free, NULL};
    carquet_set_allocator(&a);
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    g_fail_big_realloc = 1;
    carquet_status_t st = carquet_rle_encode_levels(levels, N, 8, &buf);
    g_fail_big_realloc = 0;
    carquet_buffer_destroy(&buf);
    carquet_set_allocator(NULL);

    if (st == CARQUET_OK) TEST_FAIL(name, "truncated stream reported as CARQUET_OK");

    /* Same input without the failure still round-trips. */
    carquet_buffer_init(&buf);
    st = carquet_rle_encode_levels(levels, N, 8, &buf);
    static int16_t back[N];
    int64_t n = st == CARQUET_OK
        ? carquet_rle_decode_levels(carquet_buffer_data(&buf), carquet_buffer_size(&buf),
                                    8, back, N)
        : -1;
    int ok = n == N && memcmp(back, levels, sizeof(levels)) == 0;
    carquet_buffer_destroy(&buf);
    if (!ok) TEST_FAIL(name, "round trip");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int page_starts(const char* path, int col, int64_t* starts, int cap) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) return -1;
    carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, col, &err);
    int n = -1;
    if (oi) {
        n = carquet_offset_index_num_pages(oi);
        for (int i = 0; i < n && i < cap; i++) {
            carquet_page_location_t loc;
            if (carquet_offset_index_get_page_location(oi, i, &loc) != CARQUET_OK) {
                n = -1;
                break;
            }
            starts[i] = loc.first_row_index;
        }
        carquet_offset_index_free(oi);
    }
    carquet_reader_close(r);
    return n;
}

static int test_max_rows_per_page_one_call(void) {
    const char* name = "max_rows_per_page_one_call";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_maxrows");
    carquet_error_t err = CARQUET_ERROR_INIT;

    enum { N = 3000 };
    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_column(s, "a", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_DOUBLE, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.max_rows_per_page = 1000;
    o.write_page_index = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &o, &err);
    if (!w) TEST_FAIL(name, err.message);
    static int32_t a[N];
    static double b[N];
    static int16_t def[N];
    int nb = 0;
    for (int i = 0; i < N; i++) {
        a[i] = i * 7919;
        def[i] = (int16_t)(i % 5 != 0);
        if (def[i]) b[nb++] = i * 0.5;
    }
    if (carquet_writer_write_batch(w, 0, a, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, b, N, def, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(s);

    int ok = 1;
    for (int col = 0; col < 2; col++) {
        int64_t st[8] = {0};
        int n = page_starts(path, col, st, 8);
        if (n != 3 || st[0] != 0 || st[1] != 1000 || st[2] != 2000) {
            printf("  column %d: %d pages (want 3 of 1000 rows)\n", col, n);
            ok = 0;
        }
    }
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "row cap ignored within one write_batch");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static void noop_schema_release(struct ArrowSchema* s) { s->release = NULL; }
static void noop_array_release(struct ArrowArray* a) { a->release = NULL; }

static int test_arrow_dict_null_non_nullable(void) {
    const char* name = "arrow_dict_null_non_nullable";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_dictnull");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* sc = carquet_schema_create(&err);
    if (carquet_schema_add_column(sc, "s", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_t* w = carquet_writer_create(path, sc, NULL, &err);
    if (!w) TEST_FAIL(name, err.message);

    /* Dictionary utf8 ["a", null], indices [1, 1, 1, 1], field non-nullable. */
    static const uint8_t dvalid = 0x1;
    static const int32_t doff[3] = {0, 1, 1};
    static const char dbytes[] = "a";
    const void* dbufs[3] = {&dvalid, doff, dbytes};
    struct ArrowArray dict;
    memset(&dict, 0, sizeof(dict));
    dict.length = 2; dict.null_count = 1; dict.n_buffers = 3; dict.buffers = dbufs;
    dict.release = noop_array_release;
    static const int32_t idx[4] = {1, 1, 1, 1};
    const void* cbufs[2] = {NULL, idx};
    struct ArrowArray child;
    memset(&child, 0, sizeof(child));
    child.length = 4; child.n_buffers = 2; child.buffers = cbufs;
    child.dictionary = &dict; child.release = noop_array_release;
    struct ArrowSchema dsch, csch, top;
    memset(&dsch, 0, sizeof(dsch));
    memset(&csch, 0, sizeof(csch));
    memset(&top, 0, sizeof(top));
    dsch.format = "u"; dsch.name = ""; dsch.flags = ARROW_FLAG_NULLABLE;
    dsch.release = noop_schema_release;
    csch.format = "i"; csch.name = "s"; csch.flags = 0; csch.dictionary = &dsch;
    csch.release = noop_schema_release;
    struct ArrowSchema* kids[1] = {&csch};
    struct ArrowArray* akids[1] = {&child};
    top.format = "+s"; top.name = ""; top.n_children = 1; top.children = kids;
    top.release = noop_schema_release;
    const void* tb[1] = {NULL};
    struct ArrowArray topa;
    memset(&topa, 0, sizeof(topa));
    topa.length = 4; topa.n_buffers = 1; topa.buffers = tb; topa.n_children = 1;
    topa.children = akids; topa.release = noop_array_release;

    carquet_status_t st = carquet_writer_write_arrow(w, &topa, &top, &err);
    carquet_writer_abort(w);
    carquet_schema_free(sc);
    carquet_test_cleanup(path);
    if (st != CARQUET_ERROR_INVALID_ARGUMENT) {
        printf("  write_arrow returned %d (want INVALID_ARGUMENT)\n", st);
        TEST_FAIL(name, "null in a non-nullable dictionary field accepted");
    }
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_arrow_schema_wide_struct(void) {
    const char* name = "arrow_schema_wide_struct";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_widestruct");
    carquet_error_t err = CARQUET_ERROR_INIT;

    enum { K = 80 };
    carquet_schema_t* sc = carquet_schema_create(&err);
    int32_t g = carquet_schema_add_group(sc, "g", CARQUET_REPETITION_REQUIRED, 0);
    if (g < 0) TEST_FAIL(name, "group");
    char nm[16];
    for (int i = 0; i < K; i++) {
        snprintf(nm, sizeof(nm), "c%d", i);
        if (carquet_schema_add_column(sc, nm, CARQUET_PHYSICAL_INT32, NULL,
                                      CARQUET_REPETITION_REQUIRED, 0, g) != CARQUET_OK)
            TEST_FAIL(name, "schema");
    }
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.write_arrow_schema = true;
    carquet_writer_t* w = carquet_writer_create(path, sc, &o, &err);
    if (!w) TEST_FAIL(name, err.message);
    int32_t v = 1;
    for (int i = 0; i < K; i++) {
        if (carquet_writer_write_batch(w, i, &v, 1, NULL, NULL) != CARQUET_OK)
            TEST_FAIL(name, "write");
    }
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(sc);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    int ok = r && carquet_reader_find_metadata(r, "ARROW:schema") != NULL &&
             carquet_reader_num_columns(r) == K;
    if (r) carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "ARROW:schema of an 80-field struct");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int patch_file(const char* path, const uint8_t* find, size_t flen,
                      size_t at, uint8_t value, long from, long span) {
    FILE* f = fopen(path, "r+b");
    if (!f) return -1;
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    uint8_t* buf = malloc((size_t)size);
    rewind(f);
    if (!buf || fread(buf, 1, (size_t)size, f) != (size_t)size) {
        free(buf);
        fclose(f);
        return -1;
    }
    long end = span > 0 && from + span < size ? from + span : size;
    long hit = -1;
    int hits = 0;
    for (long i = from; i + (long)flen <= end; i++) {
        if (memcmp(buf + i, find, flen) == 0) { hit = i; hits++; }
    }
    int rc = -1;
    if (hits == 1) {
        fseek(f, hit + (long)at, SEEK_SET);
        rc = fputc(value, f) == value ? 0 : -1;
    }
    free(buf);
    fclose(f);
    return rc;
}

static int64_t data_page_offset(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) return -1;
    int64_t off = r->metadata.row_groups[0].columns[0].metadata.data_page_offset;
    carquet_reader_close(r);
    return off;
}

static int test_preserved_dict_index_range(void) {
    const char* name = "preserved_dict_index_range";
    int ok = 1;
    for (int version = 1; version <= 2; version++) {
        char path[512];
        carquet_test_temp_path(path, sizeof(path), version == 1 ? "fix074_dictv1"
                                                                : "fix074_dictv2");
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* sc = carquet_schema_create(&err);
        if (carquet_schema_add_column(sc, "v", CARQUET_PHYSICAL_INT32, NULL,
                                      CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
            TEST_FAIL(name, "schema");
        carquet_writer_options_t o;
        carquet_writer_options_init(&o);
        o.write_crc = false;
        o.write_statistics = false;
        o.data_page_version = version;
        carquet_writer_t* w = carquet_writer_create(path, sc, &o, &err);
        if (!w) TEST_FAIL(name, err.message);
        /* 3 dictionary entries; indices (100 x 0)(100 x 1)(100 x 2) are three
         * RLE runs: bit width 2, then C8 01 <index> per run. */
        static int32_t vals[300];
        for (int i = 0; i < 300; i++) vals[i] = (i / 100 + 1) * 10;
        if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY)
                != CARQUET_OK ||
            carquet_writer_write_batch(w, 0, vals, 300, NULL, NULL) != CARQUET_OK ||
            carquet_writer_close(w) != CARQUET_OK)
            TEST_FAIL(name, "write");
        carquet_schema_free(sc);

        /* Last run's index 2 -> 3: past the 3-entry dictionary. */
        static const uint8_t run[4] = {0xC8, 0x01, 0x01, 0xC8};
        int64_t dpo = data_page_offset(path);
        if (dpo < 0 || patch_file(path, run, 4, 5, 0x03, (long)dpo, 0) != 0) {
            carquet_test_cleanup(path);
            TEST_FAIL(name, "could not locate the index runs");
        }

        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL(name, err.message);
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.preserve_dictionaries = true;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        carquet_status_t st = CARQUET_ERROR_INVALID_STATE;
        uint32_t max_idx = 0;
        int32_t dict_count = 0;
        if (br) {
            carquet_row_batch_t* b = NULL;
            while ((st = carquet_batch_reader_next(br, &b)) == CARQUET_OK && b) {
                const uint32_t* ix;
                const uint8_t* nb;
                int64_t n;
                const uint8_t* dd;
                const uint32_t* doffs;
                if (carquet_row_batch_column_dictionary(b, 0, &ix, &nb, &n, &dd,
                                                        &dict_count, &doffs) == CARQUET_OK) {
                    for (int64_t i = 0; i < n; i++) if (ix[i] > max_idx) max_idx = ix[i];
                }
                carquet_row_batch_free(b);
                b = NULL;
            }
            carquet_batch_reader_free(br);
        }
        carquet_reader_close(r);
        carquet_test_cleanup(path);
        if (st == CARQUET_OK || st == CARQUET_ERROR_END_OF_DATA) {
            printf("  V%d: read succeeded, max index %u for %d entries\n",
                   version, max_idx, dict_count);
            ok = 0;
        }
    }
    if (!ok) TEST_FAIL(name, "out-of-range dictionary index handed to the caller");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------ */
static int test_v2_empty_compressed_section(void) {
    const char* name = "v2_empty_compressed_section";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fix074_v2empty");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* sc = carquet_schema_create(&err);
    if (carquet_schema_add_column(sc, "v", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "schema");
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.compression = CARQUET_COMPRESSION_GZIP;
    o.data_page_version = 2;
    o.write_crc = false;
    o.write_statistics = false;
    carquet_writer_t* w = carquet_writer_create(path, sc, &o, &err);
    if (!w) TEST_FAIL(name, err.message);
    int16_t def[5] = {0, 0, 0, 0, 0};
    int32_t none = 0;
    if (carquet_writer_write_batch(w, 0, &none, 5, def, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL(name, "write");
    carquet_schema_free(sc);

    /* carquet writes is_compressed = false for an empty data section; flip
     * it to true like parquet-java (field 7 header 0x12 -> 0x11). */
    static const uint8_t f7_false[1] = {0x12};
    int64_t dpo = data_page_offset(path);
    if (dpo < 0 || patch_file(path, f7_false, 1, 0, 0x11, (long)dpo, 24) != 0) {
        carquet_test_cleanup(path);
        TEST_FAIL(name, "could not locate is_compressed");
    }

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, err.message);
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int32_t out[5];
    int16_t dl[5] = {9, 9, 9, 9, 9};
    int64_t got = c ? carquet_column_read_batch(c, out, 5, dl, NULL) : -1;
    if (c) carquet_column_reader_free(c);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (got != 5 || dl[0] != 0 || dl[4] != 0) {
        printf("  read %lld values\n", (long long)got);
        TEST_FAIL(name, "all-null V2 page with an empty GZIP section unreadable");
    }
    TEST_PASS(name);
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_dense_cursor_after_skip();
    rc |= test_list_filter_null_count();
    rc |= test_column_index_list_lengths();
    rc |= test_bloom_signed_zero();
    rc |= test_deprecated_stats_unsigned();
    rc |= test_append_shorter_footer();
    rc |= test_schema_preorder();
    rc |= test_geo_stats_per_row_group();
    rc |= test_rle_encoder_oom();
    rc |= test_max_rows_per_page_one_call();
    rc |= test_arrow_dict_null_non_nullable();
    rc |= test_arrow_schema_wide_struct();
    rc |= test_preserved_dict_index_range();
    rc |= test_v2_empty_compressed_section();
    if (rc == 0) printf("\nAll post-0.7.3 regression tests passed.\n");
    return rc;
}
