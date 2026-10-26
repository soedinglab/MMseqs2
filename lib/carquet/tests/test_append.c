/**
 * @file test_append.c
 * @brief Tests for carquet_writer_open_append() that complement the basic
 *        multi-append roundtrip already covered in test_writer_extensions.c.
 *
 * Focus here is on the v0.6.0 append guarantees the basic test does not assert:
 *   - key-value metadata is carried over and accumulates across appends
 *   - existing bloom filters and page indexes survive an append (they sit
 *     between the row-group data and the footer the writer overwrites)
 *   - OPTIONAL (nullable) columns keep their null pattern across row groups
 *   - multi-column files append correctly for every leaf
 *   - structural-mismatch appends are rejected and leave the file intact
 *     (wrong column count, wrong column name, and a missing target file)
 *
 * Everything is read back through carquet's own reader and asserted exact.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"
#include <stdio.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Test 1: metadata carryover + accumulation across appends            */
/* ------------------------------------------------------------------ */
static int test_metadata_carryover(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "append_meta");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);

    int64_t a[100], b[100];
    for (int i = 0; i < 100; i++) { a[i] = i; b[i] = 100 + i; }

    /* Initial file with two metadata keys. */
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("meta_carryover", "create"); }
    carquet_writer_add_metadata(w, "origin", "initial");
    carquet_writer_add_metadata(w, "schema_id", "42");
    carquet_writer_write_batch(w, 0, a, 100, NULL, NULL);
    carquet_writer_close(w);

    /* Append a row group and add one more metadata key. */
    w = carquet_writer_open_append(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); carquet_test_cleanup(path);
              TEST_FAIL("meta_carryover", "open_append"); }
    carquet_writer_add_metadata(w, "appended", "true");
    carquet_writer_write_batch(w, 0, b, 100, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    /* All three keys must be present after the append. */
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("meta_carryover", "open"); }
    const char* origin    = carquet_reader_find_metadata(r, "origin");
    const char* schema_id = carquet_reader_find_metadata(r, "schema_id");
    const char* appended  = carquet_reader_find_metadata(r, "appended");
    int ok = origin && strcmp(origin, "initial") == 0
          && schema_id && strcmp(schema_id, "42") == 0
          && appended && strcmp(appended, "true") == 0
          && carquet_reader_num_rows(r) == 200
          && carquet_reader_num_row_groups(r) == 2;
    carquet_reader_close(r);
    carquet_test_cleanup(path);

    if (!ok) TEST_FAIL("meta_carryover", "metadata not carried over/accumulated");
    TEST_PASS("meta_carryover");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 2: bloom filters + page indexes on row group 0 survive append  */
/* ------------------------------------------------------------------ */
static int test_bloom_pageindex_preserved(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "append_bloom");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "k", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);

    int64_t a[500], b[500];
    for (int i = 0; i < 500; i++) { a[i] = i; b[i] = 1000 + i; }

    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.write_statistics    = true;
    wo.write_bloom_filters = true;
    wo.write_page_index    = true;

    /* Initial file: row group 0 carries bloom filter + page index. */
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("append_bloom", "create"); }
    carquet_writer_write_batch(w, 0, a, 500, NULL, NULL);
    carquet_writer_close(w);

    /* Append row group 1 (also with bloom + page index). */
    w = carquet_writer_open_append(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); carquet_test_cleanup(path);
              TEST_FAIL("append_bloom", "open_append"); }
    carquet_writer_write_batch(w, 0, b, 500, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("append_bloom", "open"); }

    int ok = carquet_reader_num_row_groups(r) == 2;

    /* Row group 0's bloom filter must still resolve and answer correctly:
     * a value it contains -> maybe-present; a value far outside -> absent. */
    if (ok) {
        carquet_bloom_filter_t* bf0 = carquet_reader_get_bloom_filter(r, 0, 0, &err);
        ok = bf0 != NULL;
        if (bf0) {
            ok = ok && carquet_bloom_filter_check_i64(bf0, 250)    /* in rg0 */
                    && !carquet_bloom_filter_check_i64(bf0, 999999);
            carquet_bloom_filter_destroy(bf0);
        }
    }
    /* The appended row group 1 must also carry a usable bloom filter. */
    if (ok) {
        carquet_bloom_filter_t* bf1 = carquet_reader_get_bloom_filter(r, 1, 0, &err);
        ok = bf1 != NULL;
        if (bf1) {
            ok = ok && carquet_bloom_filter_check_i64(bf1, 1250)   /* in rg1 */
                    && !carquet_bloom_filter_check_i64(bf1, 250);  /* only in rg0 */
            carquet_bloom_filter_destroy(bf1);
        }
    }
    /* Page index on row group 0 must still be readable. */
    if (ok) {
        carquet_column_index_t* ci = carquet_reader_get_column_index(r, 0, 0, &err);
        ok = ci != NULL && carquet_column_index_num_pages(ci) >= 1;
        if (ci) carquet_column_index_free(ci);
    }
    /* Chunk metadata flags should report both features on both row groups. */
    if (ok) {
        for (int32_t rg = 0; rg < 2 && ok; rg++) {
            carquet_column_chunk_metadata_t m;
            ok = carquet_reader_column_chunk_metadata(r, rg, 0, &m) == CARQUET_OK
                 && m.has_bloom_filter && m.has_column_index;
        }
    }
    carquet_reader_close(r);
    carquet_test_cleanup(path);

    if (!ok) TEST_FAIL("append_bloom", "bloom/page-index not preserved across append");
    TEST_PASS("append_bloom");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 3: nullable (OPTIONAL) column keeps null pattern across append */
/* ------------------------------------------------------------------ */
static int test_nullable_append(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "append_null");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "n", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);

    /* Row group 0: 6 rows, every other one null. Dense values for def==1. */
    int32_t v0[] = {10, 20, 30};                 /* present values only */
    int16_t d0[] = {1, 0, 1, 0, 1, 0};           /* 6 rows, 3 present    */
    /* Row group 1: 4 rows, first two null. */
    int32_t v1[] = {40, 50};
    int16_t d1[] = {0, 0, 1, 1};                 /* 4 rows, 2 present     */

    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("append_null", "create"); }
    carquet_writer_write_batch(w, 0, v0, 6, d0, NULL);
    carquet_writer_close(w);

    w = carquet_writer_open_append(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); carquet_test_cleanup(path);
              TEST_FAIL("append_null", "open_append"); }
    carquet_writer_write_batch(w, 0, v1, 4, d1, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("append_null", "open"); }
    int ok = carquet_reader_num_rows(r) == 10 && carquet_reader_num_row_groups(r) == 2;

    /* Read each row group's def levels back and confirm the null pattern. */
    if (ok) {
        int16_t exp_def[2][6] = { {1,0,1,0,1,0}, {0,0,1,1,0,0} };
        int     exp_n[2]      = { 6, 4 };
        int32_t exp_val[2][3] = { {10,20,30}, {40,50,0} };
        for (int32_t rg = 0; rg < 2 && ok; rg++) {
            carquet_column_reader_t* c = carquet_reader_get_column(r, rg, 0, &err);
            if (!c) { ok = 0; break; }
            int32_t vals[16]; int16_t def[16];
            int64_t n = carquet_column_read_batch(c, vals, 16, def, NULL);
            carquet_column_reader_free(c);
            ok = (n == exp_n[rg]);
            int vi = 0;
            for (int i = 0; i < exp_n[rg] && ok; i++) {
                ok = (def[i] == exp_def[rg][i]);
                if (ok && def[i] == 1) { ok = (vals[vi] == exp_val[rg][vi]); vi++; }
            }
        }
    }
    carquet_reader_close(r);
    carquet_test_cleanup(path);

    if (!ok) TEST_FAIL("append_null", "null pattern not preserved across append");
    TEST_PASS("append_null");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 4: multi-column file appends every leaf correctly              */
/* ------------------------------------------------------------------ */
static int test_multicolumn_append(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "append_multicol");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "id",   CARQUET_PHYSICAL_INT64,      NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "val",  CARQUET_PHYSICAL_DOUBLE,     NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "name", CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    enum { N = 50 };
    int64_t id0[N], id1[N];
    double  val0[N], val1[N];
    carquet_byte_array_t nm0[N], nm1[N];
    static char buf0[N][16], buf1[N][16];
    for (int i = 0; i < N; i++) {
        id0[i] = i;       val0[i] = i * 1.5;
        id1[i] = N + i;   val1[i] = (N + i) * 1.5;
        snprintf(buf0[i], sizeof(buf0[i]), "a%d", i);
        snprintf(buf1[i], sizeof(buf1[i]), "b%d", i);
        nm0[i].data = (uint8_t*)buf0[i]; nm0[i].length = (int64_t)strlen(buf0[i]);
        nm1[i].data = (uint8_t*)buf1[i]; nm1[i].length = (int64_t)strlen(buf1[i]);
    }

    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("append_multicol", "create"); }
    carquet_writer_write_batch(w, 0, id0,  N, NULL, NULL);
    carquet_writer_write_batch(w, 1, val0, N, NULL, NULL);
    carquet_writer_write_batch(w, 2, nm0,  N, NULL, NULL);
    carquet_writer_close(w);

    w = carquet_writer_open_append(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); carquet_test_cleanup(path);
              TEST_FAIL("append_multicol", "open_append"); }
    carquet_writer_write_batch(w, 0, id1,  N, NULL, NULL);
    carquet_writer_write_batch(w, 1, val1, N, NULL, NULL);
    carquet_writer_write_batch(w, 2, nm1,  N, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("append_multicol", "open"); }
    int ok = carquet_reader_num_rows(r) == 2 * N
          && carquet_reader_num_columns(r) == 3
          && carquet_reader_num_row_groups(r) == 2;

    /* Spot-check the three columns across both row groups. */
    for (int32_t rg = 0; rg < 2 && ok; rg++) {
        carquet_column_reader_t* ci = carquet_reader_get_column(r, rg, 0, &err);
        carquet_column_reader_t* cv = carquet_reader_get_column(r, rg, 1, &err);
        carquet_column_reader_t* cn = carquet_reader_get_column(r, rg, 2, &err);
        if (!ci || !cv || !cn) { ok = 0; }
        if (ok) {
            int64_t ids[N]; double vals[N]; carquet_byte_array_t names[N];
            int64_t ni = carquet_column_read_batch(ci, ids,   N, NULL, NULL);
            int64_t nv = carquet_column_read_batch(cv, vals,  N, NULL, NULL);
            int64_t nn = carquet_column_read_batch(cn, names, N, NULL, NULL);
            ok = (ni == N && nv == N && nn == N);
            for (int i = 0; i < N && ok; i++) {
                int64_t base = (int64_t)rg * N;
                char expname[16];
                snprintf(expname, sizeof(expname), "%c%d", rg == 0 ? 'a' : 'b', i);
                ok = ids[i] == base + i
                  && vals[i] == (base + i) * 1.5
                  && names[i].length == (int64_t)strlen(expname)
                  && memcmp(names[i].data, expname, (size_t)names[i].length) == 0;
            }
        }
        carquet_column_reader_free(ci);
        carquet_column_reader_free(cv);
        carquet_column_reader_free(cn);
    }
    carquet_reader_close(r);
    carquet_test_cleanup(path);

    if (!ok) TEST_FAIL("append_multicol", "multi-column append mismatch");
    TEST_PASS("append_multicol");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 5: structural-mismatch appends are rejected, file left intact  */
/* ------------------------------------------------------------------ */
static int test_mismatch_rejected(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "append_mismatch");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "a", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    int32_t a[10], b[10];
    for (int i = 0; i < 10; i++) { a[i] = i; b[i] = i * 2; }

    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("append_mismatch", "create"); }
    carquet_writer_write_batch(w, 0, a, 10, NULL, NULL);
    carquet_writer_write_batch(w, 1, b, 10, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    int ok = 1;

    /* (a) Wrong column count: only one column. Must be rejected. */
    {
        carquet_schema_t* bad = carquet_schema_create(&err);
        carquet_schema_add_column(bad, "a", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_writer_t* bw = carquet_writer_open_append(path, bad, NULL, &err);
        ok = ok && (bw == NULL);
        if (bw) carquet_writer_close(bw);
        carquet_schema_free(bad);
    }

    /* (b) Wrong column name: right count/types, "a" renamed to "x". */
    {
        carquet_schema_t* bad = carquet_schema_create(&err);
        carquet_schema_add_column(bad, "x", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_schema_add_column(bad, "b", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_writer_t* bw = carquet_writer_open_append(path, bad, NULL, &err);
        ok = ok && (bw == NULL);
        if (bw) carquet_writer_close(bw);
        carquet_schema_free(bad);
    }

    /* (c) Append to a non-existent file must fail, not create one. */
    {
        char missing[512]; carquet_test_temp_path(missing, sizeof(missing), "append_does_not_exist");
        carquet_test_cleanup(missing); /* ensure absent */
        carquet_schema_t* s2 = carquet_schema_create(&err);
        carquet_schema_add_column(s2, "a", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_schema_add_column(s2, "b", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_writer_t* bw = carquet_writer_open_append(missing, s2, NULL, &err);
        ok = ok && (bw == NULL);
        if (bw) carquet_writer_close(bw);
        carquet_schema_free(s2);
        carquet_test_cleanup(missing);
    }

    /* After all rejected attempts the original file must be intact: 1 row
     * group, 2 columns, 10 rows, values unchanged. */
    if (ok) {
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        ok = r != NULL
          && carquet_reader_num_row_groups(r) == 1
          && carquet_reader_num_columns(r) == 2
          && carquet_reader_num_rows(r) == 10;
        if (ok) {
            carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
            int32_t got[10];
            int64_t n = c ? carquet_column_read_batch(c, got, 10, NULL, NULL) : -1;
            ok = (n == 10);
            for (int i = 0; i < 10 && ok; i++) ok = (got[i] == i);
            if (c) carquet_column_reader_free(c);
        }
        if (r) carquet_reader_close(r);
    }
    carquet_test_cleanup(path);

    if (!ok) TEST_FAIL("append_mismatch", "mismatch not rejected or file corrupted");
    TEST_PASS("append_mismatch");
    return 0;
}

int main(void) {
    carquet_init();
    int failures = 0;
    failures += test_metadata_carryover();
    failures += test_bloom_pageindex_preserved();
    failures += test_nullable_append();
    failures += test_multicolumn_append();
    failures += test_mismatch_rejected();

    if (failures == 0) {
        printf("\nAll append tests passed.\n");
        return 0;
    }
    printf("\n%d append test(s) failed.\n", failures);
    return 1;
}
