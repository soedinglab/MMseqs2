/**
 * @file test_arrow_c_data.c
 * @brief Tests for the Arrow C Data Interface bridge (export + import).
 *
 * Covers: schema/array export with exact buffer-byte assertions (validity
 * bitmaps, bit-packed booleans, int32 offsets), full export -> import
 * round-trip through a second writer, schema import, and error / boundary
 * paths (NULL args, nested rejection, unsupported format, column mismatch,
 * sliced arrays, all-null columns).
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Dataset (5 rows):
 *   id    INT64 REQUIRED  = 10,20,30,40,50
 *   score DOUBLE OPTIONAL = 1.5, NULL, 3.5, NULL, 5.5
 *   name  STRING OPTIONAL = "aa","b", NULL, "dddd", ""
 *   flag  BOOL   OPTIONAL = true,false, NULL, true, false
 */

static carquet_schema_t* make_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_logical_type_t str = {0};
    str.id = CARQUET_LOGICAL_STRING;
    carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "score", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_schema_add_column(s, "name", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_schema_add_column(s, "flag", CARQUET_PHYSICAL_BOOLEAN, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    return s;
}

/* Write the dataset to an in-memory Parquet buffer. Caller frees *buf. */
static int write_dataset(void** buf, size_t* sz) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) { carquet_schema_free(schema); return -1; }

    int64_t id[5] = {10, 20, 30, 40, 50};
    if (carquet_writer_write_batch(w, 0, id, 5, NULL, NULL) != CARQUET_OK) return -1;

    int16_t score_def[5] = {1, 0, 1, 0, 1};
    double score[3] = {1.5, 3.5, 5.5};
    if (carquet_writer_write_batch(w, 1, score, 5, score_def, NULL) != CARQUET_OK) return -1;

    int16_t name_def[5] = {1, 1, 0, 1, 1};
    carquet_byte_array_t name[4] = {
        {(uint8_t*)"aa", 2}, {(uint8_t*)"b", 1}, {(uint8_t*)"dddd", 4}, {(uint8_t*)"", 0}
    };
    if (carquet_writer_write_batch(w, 2, name, 5, name_def, NULL) != CARQUET_OK) return -1;

    int16_t flag_def[5] = {1, 1, 0, 1, 1};
    uint8_t flag[4] = {1, 0, 1, 0};
    if (carquet_writer_write_batch(w, 3, flag, 5, flag_def, NULL) != CARQUET_OK) return -1;

    if (carquet_writer_close(w) != CARQUET_OK) return -1;
    if (carquet_writer_get_buffer(w, buf, sz) != CARQUET_OK) return -1;
    carquet_schema_free(schema);
    return 0;
}

/* Read the single batch of the dataset. The batch data stays valid until the
 * batch reader is freed, so we return it too; the Arrow export copies out. */
static int read_one_batch(void* buf, size_t sz,
                          carquet_reader_t** out_reader,
                          carquet_batch_reader_t** out_br,
                          carquet_row_batch_t** out_batch) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) return -1;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 100;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); return -1; }
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st = carquet_batch_reader_next(br, &batch);
    if (st != CARQUET_OK || !batch) { carquet_batch_reader_free(br); carquet_reader_close(r); return -1; }
    *out_reader = r;
    *out_br = br;
    *out_batch = batch;
    return 0;
}

static void free_one_batch(carquet_reader_t* r, carquet_batch_reader_t* br) {
    if (br) carquet_batch_reader_free(br);
    if (r) carquet_reader_close(r);
}

/* ------------------------------------------------------------------ */

static int test_export_bytes(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz) != 0) TEST_FAIL("export_bytes", "write failed");

    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("export_bytes", "read failed"); }

    const carquet_schema_t* schema = carquet_reader_schema(r);
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_status_t st = carquet_arrow_export_batch(batch, schema, &aschema, &aarray, &err);
    if (st != CARQUET_OK) { free_one_batch(r, br); free(buf); TEST_FAIL("export_bytes", err.message); }

    /* Schema: struct with 4 children, expected formats + nullability. */
    assert(strcmp(aschema.format, "+s") == 0);
    assert(aschema.n_children == 4);
    assert(strcmp(aschema.children[0]->format, "l") == 0);
    assert(strcmp(aschema.children[0]->name, "id") == 0);
    assert((aschema.children[0]->flags & ARROW_FLAG_NULLABLE) == 0);
    assert(strcmp(aschema.children[1]->format, "g") == 0);
    assert((aschema.children[1]->flags & ARROW_FLAG_NULLABLE) != 0);
    assert(strcmp(aschema.children[2]->format, "u") == 0);
    assert(strcmp(aschema.children[3]->format, "b") == 0);

    /* Array: struct length 5, 4 children. */
    assert(aarray.length == 5);
    assert(aarray.n_children == 4);

    /* id: no nulls, direct int64 data. */
    const struct ArrowArray* c_id = aarray.children[0];
    assert(c_id->null_count == 0);
    assert(c_id->buffers[0] == NULL);
    const int64_t* idv = (const int64_t*)c_id->buffers[1];
    for (int i = 0; i < 5; i++) assert(idv[i] == (int64_t)(10 * (i + 1)));

    /* score: nulls at 1,3. validity byte = 0b00010101 = 0x15. */
    const struct ArrowArray* c_score = aarray.children[1];
    assert(c_score->null_count == 2);
    assert(((const uint8_t*)c_score->buffers[0])[0] == 0x15);
    const double* sv = (const double*)c_score->buffers[1];
    assert(sv[0] == 1.5 && sv[2] == 3.5 && sv[4] == 5.5);

    /* name: null at 2. validity = 0b00011011 = 0x1B. offsets/data exact. */
    const struct ArrowArray* c_name = aarray.children[2];
    assert(c_name->null_count == 1);
    assert(((const uint8_t*)c_name->buffers[0])[0] == 0x1B);
    const int32_t* off = (const int32_t*)c_name->buffers[1];
    int32_t exp_off[6] = {0, 2, 3, 3, 7, 7};
    for (int i = 0; i < 6; i++) assert(off[i] == exp_off[i]);
    assert(memcmp(c_name->buffers[2], "aabdddd", 7) == 0);

    /* flag: bit-packed bool. validity 0x1B. data bits: row0=1,row3=1 -> 0x09. */
    const struct ArrowArray* c_flag = aarray.children[3];
    assert(c_flag->null_count == 1);
    assert(((const uint8_t*)c_flag->buffers[0])[0] == 0x1B);
    assert(((const uint8_t*)c_flag->buffers[1])[0] == 0x09);

    aschema.release(&aschema);
    aarray.release(&aarray);
    free_one_batch(r, br);
    free(buf);
    TEST_PASS("export_bytes");
    return 0;
}

/* Verify the exported struct survives freeing the source batch (owned copies). */
static int test_export_owns_memory(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz) != 0) TEST_FAIL("export_owns", "write failed");
    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("export_owns", "read failed"); }

    const carquet_schema_t* schema = carquet_reader_schema(r);
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_arrow_export_batch(batch, schema, &aschema, &aarray, &err) != CARQUET_OK)
        { free_one_batch(r, br); free(buf); TEST_FAIL("export_owns", err.message); }

    /* Free source data first, then read the copies. */
    free_one_batch(r, br);
    free(buf);

    const int64_t* idv = (const int64_t*)aarray.children[0]->buffers[1];
    assert(idv[0] == 10 && idv[4] == 50);
    assert(memcmp(aarray.children[2]->buffers[2], "aabdddd", 7) == 0);

    aschema.release(&aschema);
    aarray.release(&aarray);
    TEST_PASS("export_owns");
    return 0;
}

static int test_roundtrip(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz) != 0) TEST_FAIL("roundtrip", "write failed");
    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("roundtrip", "read failed"); }

    const carquet_schema_t* schema = carquet_reader_schema(r);
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_arrow_export_batch(batch, schema, &aschema, &aarray, &err) != CARQUET_OK)
        { free_one_batch(r, br); free(buf); TEST_FAIL("roundtrip", err.message); }
    free_one_batch(r, br);
    free(buf);

    /* Write the Arrow array back through a fresh writer (consumes aschema/aarray). */
    carquet_schema_t* schema2 = make_schema();
    carquet_writer_t* w = carquet_writer_create_buffer(schema2, NULL, &err);
    if (!w) { carquet_schema_free(schema2); aschema.release(&aschema); aarray.release(&aarray);
              TEST_FAIL("roundtrip", "writer create failed"); }
    if (carquet_writer_write_arrow(w, &aarray, &aschema, &err) != CARQUET_OK) {
        carquet_schema_free(schema2);
        TEST_FAIL("roundtrip", err.message);
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(schema2); TEST_FAIL("roundtrip", "close2"); }
    void* buf2 = NULL; size_t sz2 = 0;
    if (carquet_writer_get_buffer(w, &buf2, &sz2) != CARQUET_OK) { carquet_schema_free(schema2); TEST_FAIL("roundtrip", "getbuf2"); }
    carquet_schema_free(schema2);

    /* Read back and compare against the original dataset. */
    carquet_reader_t* r2 = carquet_reader_open_buffer(buf2, sz2, NULL, &err);
    if (!r2) { free(buf2); TEST_FAIL("roundtrip", "open2"); }
    assert(carquet_reader_num_rows(r2) == 5);

    carquet_column_reader_t* c0 = carquet_reader_get_column(r2, 0, 0, NULL);
    int64_t idv[5]; assert(carquet_column_read_batch(c0, idv, 5, NULL, NULL) == 5);
    for (int i = 0; i < 5; i++) assert(idv[i] == (int64_t)(10 * (i + 1)));
    carquet_column_reader_free(c0);

    carquet_column_reader_t* c1 = carquet_reader_get_column(r2, 0, 1, NULL);
    double sv[5]; int16_t sd[5];
    int64_t sn = carquet_column_read_batch(c1, sv, 5, sd, NULL);
    assert(sn == 5);
    assert(sd[0] == 1 && sd[1] == 0 && sd[2] == 1 && sd[3] == 0 && sd[4] == 1);
    assert(sv[0] == 1.5 && sv[1] == 3.5 && sv[2] == 5.5);  /* dense non-null */
    carquet_column_reader_free(c1);

    carquet_column_reader_t* c2 = carquet_reader_get_column(r2, 0, 2, NULL);
    carquet_byte_array_t nv[5]; int16_t nd[5];
    int64_t nn = carquet_column_read_batch(c2, nv, 5, nd, NULL);
    assert(nn == 5);
    assert(nd[0] == 1 && nd[1] == 1 && nd[2] == 0 && nd[3] == 1 && nd[4] == 1);
    assert(nv[0].length == 2 && memcmp(nv[0].data, "aa", 2) == 0);
    assert(nv[1].length == 1 && nv[1].data[0] == 'b');
    assert(nv[2].length == 4 && memcmp(nv[2].data, "dddd", 4) == 0);
    assert(nv[3].length == 0);
    carquet_column_reader_free(c2);

    carquet_column_reader_t* c3 = carquet_reader_get_column(r2, 0, 3, NULL);
    uint8_t fv[5]; int16_t fd[5];
    int64_t fn = carquet_column_read_batch(c3, fv, 5, fd, NULL);
    assert(fn == 5);
    assert(fd[0] == 1 && fd[1] == 1 && fd[2] == 0 && fd[3] == 1 && fd[4] == 1);
    assert(fv[0] == 1 && fv[1] == 0 && fv[2] == 1 && fv[3] == 0);  /* dense */
    carquet_column_reader_free(c3);

    carquet_reader_close(r2);
    free(buf2);
    TEST_PASS("roundtrip");
    return 0;
}

static int test_import_schema(void) {
    carquet_schema_t* s = make_schema();
    struct ArrowSchema aschema;
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_arrow_export_schema(s, &aschema, &err) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("import_schema", err.message); }
    carquet_schema_free(s);

    carquet_schema_t* imported = NULL;
    if (carquet_arrow_import_schema(&aschema, &imported, &err) != CARQUET_OK)
        TEST_FAIL("import_schema", err.message);  /* aschema consumed */

    assert(carquet_schema_num_columns(imported) == 4);
    assert(strcmp(carquet_schema_column_name(imported, 0), "id") == 0);
    assert(carquet_schema_column_type(imported, 0) == CARQUET_PHYSICAL_INT64);
    assert(carquet_schema_max_def_level(imported, 0) == 0);   /* required */
    assert(carquet_schema_column_type(imported, 1) == CARQUET_PHYSICAL_DOUBLE);
    assert(carquet_schema_max_def_level(imported, 1) == 1);   /* optional */
    assert(carquet_schema_column_type(imported, 2) == CARQUET_PHYSICAL_BYTE_ARRAY);
    assert(carquet_schema_column_type(imported, 3) == CARQUET_PHYSICAL_BOOLEAN);
    carquet_schema_free(imported);
    TEST_PASS("import_schema");
    return 0;
}

/* ---- Error and boundary paths ---- */

static void noop_schema_release(struct ArrowSchema* s) { s->release = NULL; }

static int test_errors(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    struct ArrowSchema aschema;

    /* NULL args. */
    assert(carquet_arrow_export_schema(NULL, &aschema, &err) == CARQUET_ERROR_INVALID_ARGUMENT);

    /* Nested/repeated column export is now supported: a LIST<int32> yields a
     * top-level struct child with Arrow format "+l" whose element is "i". */
    carquet_schema_t* ls = carquet_schema_create(NULL);
    carquet_schema_add_list(ls, "vals", CARQUET_PHYSICAL_INT32, NULL,
                            CARQUET_REPETITION_OPTIONAL, 0, 0);
    assert(carquet_arrow_export_schema(ls, &aschema, &err) == CARQUET_OK);
    assert(aschema.n_children == 1);
    assert(strcmp(aschema.children[0]->format, "+l") == 0);
    assert(aschema.children[0]->n_children == 1);
    assert(strcmp(aschema.children[0]->children[0]->format, "i") == 0);
    aschema.release(&aschema);
    carquet_schema_free(ls);

    /* Import: non-struct top-level rejected. */
    struct ArrowSchema flat = {0};
    flat.format = "i";
    flat.release = noop_schema_release;
    carquet_schema_t* out = NULL;
    assert(carquet_arrow_import_schema(&flat, &out, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    assert(out == NULL);

    /* Import: unsupported child format rejected. */
    struct ArrowSchema bad_child = {0};
    bad_child.format = "ZZZ";
    bad_child.name = "x";
    bad_child.release = noop_schema_release;
    struct ArrowSchema* kids[1] = {&bad_child};
    struct ArrowSchema root = {0};
    root.format = "+s";
    root.n_children = 1;
    root.children = kids;
    root.release = noop_schema_release;
    assert(carquet_arrow_import_schema(&root, &out, &err) == CARQUET_ERROR_NOT_IMPLEMENTED);

    /* Truncated timestamp formats must be rejected without reading past the
     * NUL terminator (regression: "ts" over-read fmt[3], found by fuzzing). */
    const char* truncated[] = { "ts", "tt", "td", "t", "w" };
    for (size_t i = 0; i < sizeof(truncated) / sizeof(truncated[0]); i++) {
        struct ArrowSchema tchild = {0};
        tchild.format = truncated[i];
        tchild.name = "x";
        tchild.release = noop_schema_release;
        struct ArrowSchema* tkids[1] = {&tchild};
        struct ArrowSchema troot = {0};
        troot.format = "+s";
        troot.n_children = 1;
        troot.children = tkids;
        troot.release = noop_schema_release;
        carquet_status_t rc = carquet_arrow_import_schema(&troot, &out, &err);
        assert(rc != CARQUET_OK);  /* rejected, no crash / over-read */
        assert(out == NULL);
    }

    TEST_PASS("errors");
    return 0;
}

static int test_write_arrow_errors(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz) != 0) TEST_FAIL("write_errors", "write failed");
    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("write_errors", "read failed"); }
    const carquet_schema_t* schema = carquet_reader_schema(r);

    /* Column-count mismatch: 4-col array into a 2-col writer. */
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_arrow_export_batch(batch, schema, &aschema, &aarray, &err) != CARQUET_OK)
        { free_one_batch(r, br); free(buf); TEST_FAIL("write_errors", err.message); }

    carquet_schema_t* small = carquet_schema_create(NULL);
    carquet_schema_add_column(small, "id", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(small, "score", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_writer_t* w = carquet_writer_create_buffer(small, NULL, &err);
    carquet_status_t st = carquet_writer_write_arrow(w, &aarray, &aschema, &err);
    assert(st == CARQUET_ERROR_INVALID_ARGUMENT);  /* consumes aarray/aschema */
    /* Abandon the writer (never closed): free via close path. */
    (void)carquet_writer_close(w);
    void* tmp = NULL; size_t tsz = 0;
    (void)carquet_writer_get_buffer(w, &tmp, &tsz);
    free(tmp);
    carquet_schema_free(small);

    /* A child shorter than the struct's offset + length is rejected (it would
     * be read past its end). */
    struct ArrowSchema aschema2; struct ArrowArray aarray2;
    if (carquet_arrow_export_batch(batch, schema, &aschema2, &aarray2, &err) != CARQUET_OK)
        { free_one_batch(r, br); free(buf); TEST_FAIL("write_errors", "export2"); }
    aarray2.offset = 1;  /* slice [1, 6) of a 5-row struct: child 5 < 1 + 5 */
    carquet_schema_t* full = make_schema();
    carquet_writer_t* w2 = carquet_writer_create_buffer(full, NULL, &err);
    st = carquet_writer_write_arrow(w2, &aarray2, &aschema2, &err);
    assert(st == CARQUET_ERROR_INVALID_ARGUMENT);
    (void)carquet_writer_close(w2);
    tmp = NULL; tsz = 0;
    (void)carquet_writer_get_buffer(w2, &tmp, &tsz);
    free(tmp);
    carquet_schema_free(full);

    free_one_batch(r, br);
    free(buf);
    TEST_PASS("write_errors");
    return 0;
}

/* Sliced arrays (offset != 0): a struct slice (offset on the top-level array,
 * applied to every child) and a child slice (offset on one leaf) must write
 * exactly the rows of the slice. */
static int test_sliced_arrays(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz) != 0) TEST_FAIL("sliced", "write failed");
    carquet_error_t err = CARQUET_ERROR_INIT;

    for (int mode = 0; mode < 2; mode++) {
        carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
        if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("sliced", "read failed"); }
        struct ArrowSchema as; struct ArrowArray aa;
        if (carquet_arrow_export_batch(batch, carquet_reader_schema(r), &as, &aa, &err) != CARQUET_OK)
            { free_one_batch(r, br); free(buf); TEST_FAIL("sliced", err.message); }
        free_one_batch(r, br);

        if (mode == 0) {
            /* struct slice: rows 1..4 of every column */
            aa.offset = 1; aa.length = 4;
        } else {
            /* child slice: only `id` is shifted by one row; the struct keeps
             * rows 0..3 of the other columns */
            aa.length = 4;
            aa.children[0]->offset = 1; aa.children[0]->length = 4;
        }

        carquet_schema_t* schema2 = make_schema();
        carquet_writer_t* w = carquet_writer_create_buffer(schema2, NULL, &err);
        if (carquet_writer_write_arrow(w, &aa, &as, &err) != CARQUET_OK)
            { carquet_schema_free(schema2); free(buf); TEST_FAIL("sliced", err.message); }
        if (carquet_writer_close(w) != CARQUET_OK) { free(buf); TEST_FAIL("sliced", "close"); }
        void* buf2 = NULL; size_t sz2 = 0;
        if (carquet_writer_get_buffer(w, &buf2, &sz2) != CARQUET_OK) { free(buf); TEST_FAIL("sliced", "getbuf"); }
        carquet_schema_free(schema2);

        carquet_reader_t* r2 = carquet_reader_open_buffer(buf2, sz2, NULL, &err);
        if (!r2) { free(buf2); free(buf); TEST_FAIL("sliced", "open2"); }
        assert(carquet_reader_num_rows(r2) == 4);

        carquet_column_reader_t* c0 = carquet_reader_get_column(r2, 0, 0, NULL);
        int64_t idv[4]; assert(carquet_column_read_batch(c0, idv, 4, NULL, NULL) == 4);
        for (int i = 0; i < 4; i++) assert(idv[i] == (int64_t)(10 * (i + 2)));  /* 20..50 both modes */
        carquet_column_reader_free(c0);

        carquet_column_reader_t* c1 = carquet_reader_get_column(r2, 0, 1, NULL);
        double sv[4]; int16_t sd[4];
        assert(carquet_column_read_batch(c1, sv, 4, sd, NULL) == 4);
        if (mode == 0) {
            assert(sd[0] == 0 && sd[1] == 1 && sd[2] == 0 && sd[3] == 1);
            assert(sv[0] == 3.5 && sv[1] == 5.5);
        } else {
            assert(sd[0] == 1 && sd[1] == 0 && sd[2] == 1 && sd[3] == 0);
            assert(sv[0] == 1.5 && sv[1] == 3.5);
        }
        carquet_column_reader_free(c1);

        carquet_column_reader_t* c2 = carquet_reader_get_column(r2, 0, 2, NULL);
        carquet_byte_array_t nv[4]; int16_t nd[4];
        assert(carquet_column_read_batch(c2, nv, 4, nd, NULL) == 4);
        if (mode == 0) {
            assert(nd[0] == 1 && nd[1] == 0 && nd[2] == 1 && nd[3] == 1);
            assert(nv[0].length == 1 && nv[0].data[0] == 'b');
            assert(nv[1].length == 4 && memcmp(nv[1].data, "dddd", 4) == 0);
            assert(nv[2].length == 0);
        } else {
            assert(nd[0] == 1 && nd[1] == 1 && nd[2] == 0 && nd[3] == 1);
            assert(nv[0].length == 2 && memcmp(nv[0].data, "aa", 2) == 0);
            assert(nv[1].length == 1 && nv[1].data[0] == 'b');
            assert(nv[2].length == 4 && memcmp(nv[2].data, "dddd", 4) == 0);
        }
        carquet_column_reader_free(c2);

        carquet_column_reader_t* c3 = carquet_reader_get_column(r2, 0, 3, NULL);
        uint8_t fv[4]; int16_t fd[4];
        assert(carquet_column_read_batch(c3, fv, 4, fd, NULL) == 4);
        if (mode == 0) {
            assert(fd[0] == 1 && fd[1] == 0 && fd[2] == 1 && fd[3] == 1);
            assert(fv[0] == 0 && fv[1] == 1 && fv[2] == 0);
        } else {
            assert(fd[0] == 1 && fd[1] == 1 && fd[2] == 0 && fd[3] == 1);
            assert(fv[0] == 1 && fv[1] == 0 && fv[2] == 1);
        }
        carquet_column_reader_free(c3);
        carquet_reader_close(r2);
        free(buf2);
    }
    free(buf);
    TEST_PASS("sliced_arrays");
    return 0;
}

/* ---- minimal Arrow C struct builders for the conversion tests ---- */
static void rel_schema_b(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format); free((void*)s->name);
    for (int64_t i = 0; i < s->n_children; i++)
        if (s->children[i]) { if (s->children[i]->release) s->children[i]->release(s->children[i]); free(s->children[i]); }
    if (s->dictionary) { if (s->dictionary->release) s->dictionary->release(s->dictionary); free(s->dictionary); }
    free(s->children); s->release = NULL;
}
static void rel_array_b(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) { for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]); free(a->buffers); }
    if (a->dictionary) { if (a->dictionary->release) a->dictionary->release(a->dictionary); free(a->dictionary); }
    for (int64_t i = 0; i < a->n_children; i++)
        if (a->children[i]) { if (a->children[i]->release) a->children[i]->release(a->children[i]); free(a->children[i]); }
    free(a->children); a->release = NULL;
}
static char* dup_b(const char* s) { char* o = malloc(strlen(s) + 1); strcpy(o, s); return o; }
static struct ArrowSchema* SB(const char* fmt, const char* name, int nullable, int nch) {
    struct ArrowSchema* s = calloc(1, sizeof(*s));
    s->format = dup_b(fmt); s->name = dup_b(name); s->flags = nullable ? ARROW_FLAG_NULLABLE : 0;
    s->n_children = nch; s->children = nch ? calloc(nch, sizeof(void*)) : NULL; s->release = rel_schema_b;
    return s;
}
static struct ArrowArray* AB(int64_t len, int nbuf, int nch) {
    struct ArrowArray* a = calloc(1, sizeof(*a));
    a->length = len; a->null_count = -1; a->n_buffers = nbuf; a->n_children = nch;
    a->buffers = nbuf ? calloc(nbuf, sizeof(void*)) : NULL;
    a->children = nch ? calloc(nch, sizeof(void*)) : NULL;
    a->release = rel_array_b;
    return a;
}
static void* dupmem(const void* src, size_t n) { void* o = malloc(n ? n : 1); memcpy(o, src, n); return o; }

/* Write a one-field struct (schema `sc` / array `top`, both heap nodes
 * released by the call) to a buffer, then read it back with read_arrow. */
static int write_read_arrow(struct ArrowSchema* sc, struct ArrowArray* top,
                            struct ArrowSchema* sc2, void** buf, size_t* sz,
                            carquet_error_t* err) {
    carquet_schema_t* cs = NULL;
    if (carquet_arrow_import_schema(sc, &cs, err) != CARQUET_OK) return -1;
    free(sc);
    carquet_writer_t* w = carquet_writer_create_buffer(cs, NULL, err);
    if (!w) return -1;
    carquet_status_t rc = carquet_writer_write_arrow(w, top, sc2, err);
    free(top); free(sc2);
    if (rc != CARQUET_OK) return -1;
    if (carquet_writer_close(w) != CARQUET_OK) return -1;
    if (carquet_writer_get_buffer(w, buf, sz) != CARQUET_OK) return -1;
    carquet_schema_free(cs);
    return 0;
}

/* Arrow decimal128 "d:p,s" -> FIXED_LEN_BYTE_ARRAY DECIMAL(p, s) with the
 * parquet-cpp byte length, big-endian two's complement. */
static int test_decimal_import(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    /* DECIMAL(10,2): 5 bytes. Values 123456, -1, 0, -123456 as 16-byte LE. */
    uint8_t le[4][16] = {{0}};
    int64_t v[4] = {123456, -1, 0, -123456};
    for (int i = 0; i < 4; i++) {
        uint64_t lo = (uint64_t)v[i], hi = v[i] < 0 ? UINT64_MAX : 0;
        memcpy(le[i], &lo, 8); memcpy(le[i] + 8, &hi, 8);
    }
    struct ArrowSchema* sc = SB("+s", "schema", 0, 1);
    sc->children[0] = SB("d:10,2", "amount", 0, 0);
    struct ArrowSchema* sc2 = SB("+s", "schema", 0, 1);
    sc2->children[0] = SB("d:10,2", "amount", 0, 0);
    struct ArrowArray* leaf = AB(4, 2, 0); leaf->buffers[1] = dupmem(le, sizeof(le));
    struct ArrowArray* top = AB(4, 1, 1); top->children[0] = leaf;

    void* buf = NULL; size_t sz = 0;
    if (write_read_arrow(sc, top, sc2, &buf, &sz, &err) != 0) TEST_FAIL("decimal_import", err.message);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("decimal_import", "open");
    const carquet_schema_t* rs = carquet_reader_schema(r);
    const carquet_schema_node_t* col = carquet_schema_get_element(rs, 1);  /* 0 is the root */
    assert(col && carquet_schema_node_physical_type(col) == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY);
    assert(carquet_schema_node_type_length(col) == 5);
    const carquet_logical_type_t* clt = carquet_schema_node_logical_type(col);
    assert(clt && clt->id == CARQUET_LOGICAL_DECIMAL);
    assert(clt->params.decimal.precision == 10 && clt->params.decimal.scale == 2);

    struct ArrowSchema as; struct ArrowArray aa;
    if (carquet_reader_read_arrow(r, 0, &as, &aa, &err) != CARQUET_OK) TEST_FAIL("decimal_import", err.message);
    assert(aa.length == 4 && strcmp(as.children[0]->format, "w:5") == 0);
    const uint8_t* be = (const uint8_t*)aa.children[0]->buffers[1];
    static const uint8_t exp[4][5] = {
        {0x00, 0x00, 0x01, 0xE2, 0x40}, {0xFF, 0xFF, 0xFF, 0xFF, 0xFF},
        {0x00, 0x00, 0x00, 0x00, 0x00}, {0xFF, 0xFF, 0xFE, 0x1D, 0xC0} };
    assert(memcmp(be, exp, sizeof(exp)) == 0);
    as.release(&as); aa.release(&aa);
    carquet_reader_close(r);
    free(buf);

    /* decimal256 and decimal64 forms; precision beyond the bit width rejected. */
    carquet_schema_t* cs = NULL;
    struct ArrowSchema* s3 = SB("+s", "schema", 0, 1);
    s3->children[0] = SB("d:50,3,256", "big", 1, 0);
    assert(carquet_arrow_import_schema(s3, &cs, &err) == CARQUET_OK);
    free(s3);
    col = carquet_schema_get_element(cs, 1);
    assert(carquet_schema_node_type_length(col) == 21);
    assert(carquet_schema_node_logical_type(col)->params.decimal.precision == 50);
    carquet_schema_free(cs); cs = NULL;
    struct ArrowSchema* s4 = SB("+s", "schema", 0, 1);
    s4->children[0] = SB("d:18,0,64", "mid", 1, 0);
    assert(carquet_arrow_import_schema(s4, &cs, &err) == CARQUET_OK);
    free(s4);
    assert(carquet_schema_node_type_length(carquet_schema_get_element(cs, 1)) == 8);
    carquet_schema_free(cs); cs = NULL;
    struct ArrowSchema* s5 = SB("+s", "schema", 0, 1);
    s5->children[0] = SB("d:19,0,64", "toobig", 1, 0);
    assert(carquet_arrow_import_schema(s5, &cs, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    free(s5);
    assert(cs == NULL);

    TEST_PASS("decimal_import");
    return 0;
}

/* date64 / time32[s] / timestamp[s] are coerced to DATE / TIME(MILLIS) /
 * TIMESTAMP(MILLIS), as pyarrow does. */
static int test_temporal_coercions(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    int64_t d64[3] = { 3 * 86400000LL + 5, -1, 0 };          /* days 3, -1, 0 */
    int32_t t32[3] = { 3661, 0, 86399 };                      /* ms 3661000, 0, 86399000 */
    int64_t ts[3]  = { 1700000000LL, -1, 0 };                 /* ms *1000 */

    struct ArrowSchema* sc = SB("+s", "schema", 0, 3);
    sc->children[0] = SB("tdm", "d", 0, 0);
    sc->children[1] = SB("tts", "t", 0, 0);
    sc->children[2] = SB("tss:UTC", "ts", 0, 0);
    struct ArrowSchema* sc2 = SB("+s", "schema", 0, 3);
    sc2->children[0] = SB("tdm", "d", 0, 0);
    sc2->children[1] = SB("tts", "t", 0, 0);
    sc2->children[2] = SB("tss:UTC", "ts", 0, 0);
    struct ArrowArray* top = AB(3, 1, 3);
    top->children[0] = AB(3, 2, 0); top->children[0]->buffers[1] = dupmem(d64, sizeof(d64));
    top->children[1] = AB(3, 2, 0); top->children[1]->buffers[1] = dupmem(t32, sizeof(t32));
    top->children[2] = AB(3, 2, 0); top->children[2]->buffers[1] = dupmem(ts, sizeof(ts));

    void* buf = NULL; size_t sz = 0;
    if (write_read_arrow(sc, top, sc2, &buf, &sz, &err) != 0) TEST_FAIL("temporal_coercions", err.message);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("temporal_coercions", "open");
    struct ArrowSchema as; struct ArrowArray aa;
    if (carquet_reader_read_arrow(r, 0, &as, &aa, &err) != CARQUET_OK) TEST_FAIL("temporal_coercions", err.message);
    assert(strcmp(as.children[0]->format, "tdD") == 0);
    assert(strcmp(as.children[1]->format, "ttm") == 0);
    assert(strncmp(as.children[2]->format, "tsm:", 4) == 0);
    const int32_t* dv = (const int32_t*)aa.children[0]->buffers[1];
    assert(dv[0] == 3 && dv[1] == -1 && dv[2] == 0);
    const int32_t* tv = (const int32_t*)aa.children[1]->buffers[1];
    assert(tv[0] == 3661000 && tv[1] == 0 && tv[2] == 86399000);
    const int64_t* sv = (const int64_t*)aa.children[2]->buffers[1];
    assert(sv[0] == 1700000000000LL && sv[1] == -1000 && sv[2] == 0);
    as.release(&as); aa.release(&aa);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS("temporal_coercions");
    return 0;
}

/* fixed_size_list "+w:n" is written as a LIST; a null slot still occupies n
 * child values. Also exercises a sliced list array. */
static int test_fixed_size_list(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    /* [[1,2],[3,4],NULL,[7,8]] sliced to [1, 4): [[3,4],NULL,[7,8]] */
    int32_t vals[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    uint8_t valid = 0x0B;  /* slots 0,1,3 present; slot 2 null */
    struct ArrowSchema* sc = SB("+s", "schema", 0, 1);
    sc->children[0] = SB("+w:2", "pair", 1, 1);
    sc->children[0]->children[0] = SB("i", "item", 0, 0);
    struct ArrowSchema* sc2 = SB("+s", "schema", 0, 1);
    sc2->children[0] = SB("+w:2", "pair", 1, 1);
    sc2->children[0]->children[0] = SB("i", "item", 0, 0);
    struct ArrowArray* leaf = AB(8, 2, 0); leaf->buffers[1] = dupmem(vals, sizeof(vals));
    struct ArrowArray* lst = AB(3, 1, 1); lst->buffers[0] = dupmem(&valid, 1);
    lst->offset = 1; lst->children[0] = leaf;
    struct ArrowArray* top = AB(3, 1, 1); top->children[0] = lst;

    void* buf = NULL; size_t sz = 0;
    if (write_read_arrow(sc, top, sc2, &buf, &sz, &err) != 0) TEST_FAIL("fixed_size_list", err.message);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("fixed_size_list", "open");
    assert(carquet_reader_num_rows(r) == 3);
    struct ArrowSchema as; struct ArrowArray aa;
    if (carquet_reader_read_arrow(r, 0, &as, &aa, &err) != CARQUET_OK) TEST_FAIL("fixed_size_list", err.message);
    assert(strcmp(as.children[0]->format, "+l") == 0);
    struct ArrowArray* l = aa.children[0];
    assert(l->length == 3);
    const int32_t* off = (const int32_t*)l->buffers[1];
    assert(off[0] == 0 && off[1] == 2 && off[2] == 2 && off[3] == 4);
    const uint8_t* lv = (const uint8_t*)l->buffers[0];
    assert(lv && (lv[0] & 7) == 5);  /* present, null, present */
    const int32_t* iv = (const int32_t*)l->children[0]->buffers[1];
    assert(iv[0] == 3 && iv[1] == 4 && iv[2] == 7 && iv[3] == 8);
    as.release(&as); aa.release(&aa);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS("fixed_size_list");
    return 0;
}

/* A dictionary whose value slot is itself null yields a null row. */
static int test_dictionary_null_slot(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    /* dictionary values: [10, NULL, 30]; indices: [2, 1, 0, 1] -> 30, NULL, 10, NULL */
    int32_t dvals[3] = {10, 0, 30};
    uint8_t dvalid = 0x05;
    int8_t idx[4] = {2, 1, 0, 1};
    struct ArrowSchema* sc = SB("+s", "schema", 0, 1);
    sc->children[0] = SB("c", "cat", 1, 0);
    sc->children[0]->dictionary = SB("i", "", 1, 0);
    struct ArrowSchema* sc2 = SB("+s", "schema", 0, 1);
    sc2->children[0] = SB("c", "cat", 1, 0);
    sc2->children[0]->dictionary = SB("i", "", 1, 0);
    struct ArrowArray* dict = AB(3, 2, 0);
    dict->buffers[0] = dupmem(&dvalid, 1); dict->buffers[1] = dupmem(dvals, sizeof(dvals));
    struct ArrowArray* ind = AB(4, 2, 0); ind->buffers[1] = dupmem(idx, sizeof(idx));
    ind->dictionary = dict;
    struct ArrowArray* top = AB(4, 1, 1); top->children[0] = ind;

    void* buf = NULL; size_t sz = 0;
    if (write_read_arrow(sc, top, sc2, &buf, &sz, &err) != 0) TEST_FAIL("dictionary_null_slot", err.message);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("dictionary_null_slot", "open");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, NULL);
    int32_t v[4]; int16_t d[4];
    assert(carquet_column_read_batch(c, v, 4, d, NULL) == 4);
    assert(d[0] == 1 && d[1] == 0 && d[2] == 1 && d[3] == 0);
    assert(v[0] == 30 && v[1] == 10);
    carquet_column_reader_free(c);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS("dictionary_null_slot");
    return 0;
}

/* Boundary: a single all-null optional column. */
static int test_all_null(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_column(s, "x", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_writer_t* w = carquet_writer_create_buffer(s, NULL, &err);
    int16_t def[4] = {0, 0, 0, 0};
    int32_t no_values = 0;  /* 0 present values; pointer is unused but non-NULL */
    if (carquet_writer_write_batch(w, 0, &no_values, 4, def, NULL) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("all_null", "write failed"); }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); TEST_FAIL("all_null", "close"); }
    void* buf = NULL; size_t sz = 0;
    carquet_writer_get_buffer(w, &buf, &sz);
    carquet_schema_free(s);

    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("all_null", "read failed"); }
    const carquet_schema_t* rs = carquet_reader_schema(r);
    struct ArrowArray aarray;
    if (carquet_arrow_export_batch(batch, rs, NULL, &aarray, &err) != CARQUET_OK)
        { free_one_batch(r, br); free(buf); TEST_FAIL("all_null", err.message); }
    assert(aarray.length == 4);
    assert(aarray.children[0]->null_count == 4);
    /* All validity bits clear. */
    assert((((const uint8_t*)aarray.children[0]->buffers[0])[0] & 0x0F) == 0);
    aarray.release(&aarray);
    free_one_batch(r, br);
    free(buf);
    TEST_PASS("all_null");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Nested batch export: single-level LIST<int32> becomes Arrow "+l".    */
static int test_export_list(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_list(s, "tags", CARQUET_PHYSICAL_INT32, NULL,
                            CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_writer_t* w = carquet_writer_create_buffer(s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("export_list", "writer"); }
    /* rows: [10,20], [], [30]  ->  offsets {0,2,2,3}, values {10,20,30} */
    int32_t offsets[4] = {0, 2, 2, 3};
    int32_t values[3] = {10, 20, 30};
    if (carquet_writer_write_list_column(w, 0, 3, offsets, NULL, values, NULL, &err) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("export_list", err.message); }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); TEST_FAIL("export_list", "close"); }
    void* buf = NULL; size_t sz = 0;
    carquet_writer_get_buffer(w, &buf, &sz);
    carquet_schema_free(s);

    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("export_list", "read failed"); }
    const carquet_schema_t* rs = carquet_reader_schema(r);
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_status_t st = carquet_arrow_export_batch(batch, rs, &aschema, &aarray, &err);
    if (st != CARQUET_OK) { free_one_batch(r, br); free(buf); TEST_FAIL("export_list", err.message); }

    /* Schema: struct { tags: list<element:int32> }. */
    assert(strcmp(aschema.format, "+s") == 0);
    assert(aschema.n_children == 1);
    assert(strcmp(aschema.children[0]->format, "+l") == 0);
    assert(strcmp(aschema.children[0]->name, "tags") == 0);
    assert((aschema.children[0]->flags & ARROW_FLAG_NULLABLE) != 0);
    assert(aschema.children[0]->n_children == 1);
    assert(strcmp(aschema.children[0]->children[0]->format, "i") == 0);

    /* Array: 3 lists, offsets {0,2,2,3}, flattened element values {10,20,30}. */
    assert(aarray.length == 3);
    assert(aarray.n_children == 1);
    const struct ArrowArray* list = aarray.children[0];
    assert(list->length == 3);
    assert(list->n_buffers == 2);
    const int32_t* off = (const int32_t*)list->buffers[1];
    int32_t exp_off[4] = {0, 2, 2, 3};
    for (int i = 0; i < 4; i++) assert(off[i] == exp_off[i]);
    assert(list->n_children == 1);
    const struct ArrowArray* elem = list->children[0];
    assert(elem->length == 3);
    const int32_t* ev = (const int32_t*)elem->buffers[1];
    assert(ev[0] == 10 && ev[1] == 20 && ev[2] == 30);

    aarray.release(&aarray);
    aschema.release(&aschema);
    free_one_batch(r, br);
    free(buf);
    TEST_PASS("export_list");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Nested batch export: MAP<int32,int64> becomes Arrow "+m".            */
static int test_export_map(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_map(s, "m", CARQUET_PHYSICAL_INT32, NULL, 0,
                           CARQUET_PHYSICAL_INT64, NULL, 0,
                           CARQUET_REPETITION_OPTIONAL, 0);
    carquet_writer_t* w = carquet_writer_create_buffer(s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("export_map", "writer"); }
    /* rows: {1:100, 2:200}, {3:300}  ->  offsets {0,2,3} */
    int32_t offsets[3] = {0, 2, 3};
    int32_t keys[3] = {1, 2, 3};
    int64_t vals[3] = {100, 200, 300};
    /* leaf 0 = key (REQUIRED), leaf 1 = value; both share the entry offsets. */
    if (carquet_writer_write_list_column(w, 0, 2, offsets, NULL, keys, NULL, &err) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("export_map", err.message); }
    if (carquet_writer_write_list_column(w, 1, 2, offsets, NULL, vals, NULL, &err) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("export_map", err.message); }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); TEST_FAIL("export_map", "close"); }
    void* buf = NULL; size_t sz = 0;
    carquet_writer_get_buffer(w, &buf, &sz);
    carquet_schema_free(s);

    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("export_map", "read failed"); }
    const carquet_schema_t* rs = carquet_reader_schema(r);
    struct ArrowSchema aschema; struct ArrowArray aarray;
    carquet_status_t st = carquet_arrow_export_batch(batch, rs, &aschema, &aarray, &err);
    if (st != CARQUET_OK) { free_one_batch(r, br); free(buf); TEST_FAIL("export_map", err.message); }

    /* Schema: struct { m: map<entries: struct<key:int32, value:int64>> }. */
    assert(aschema.n_children == 1);
    assert(strcmp(aschema.children[0]->format, "+m") == 0);
    const struct ArrowSchema* entries = aschema.children[0]->children[0];
    assert(strcmp(entries->format, "+s") == 0);
    assert(entries->n_children == 2);
    assert(strcmp(entries->children[0]->format, "i") == 0);   /* key int32 */
    assert(strcmp(entries->children[1]->format, "l") == 0);   /* value int64 */

    /* Array: 2 maps, entry offsets {0,2,3}; entries struct holds 3 kv pairs. */
    assert(aarray.length == 2);
    const struct ArrowArray* map = aarray.children[0];
    assert(map->length == 2);
    assert(map->n_buffers == 2);
    const int32_t* off = (const int32_t*)map->buffers[1];
    assert(off[0] == 0 && off[1] == 2 && off[2] == 3);
    const struct ArrowArray* ent = map->children[0];
    assert(ent->length == 3);
    assert(ent->n_children == 2);
    const int32_t* ek = (const int32_t*)ent->children[0]->buffers[1];
    const int64_t* evv = (const int64_t*)ent->children[1]->buffers[1];
    assert(ek[0] == 1 && ek[1] == 2 && ek[2] == 3);
    assert(evv[0] == 100 && evv[1] == 200 && evv[2] == 300);

    aarray.release(&aarray);
    aschema.release(&aschema);
    free_one_batch(r, br);
    free(buf);
    TEST_PASS("export_map");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Rejection path: a STRUCT top-level field is not exportable from a    */
/* row batch (use carquet_reader_read_arrow) — must return NOT_IMPLEMENTED. */
static int test_export_struct_rejected(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(NULL);
    int32_t g = carquet_schema_add_group(s, "point", CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(s, "x", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, g);
    carquet_schema_add_column(s, "y", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, g);
    carquet_writer_t* w = carquet_writer_create_buffer(s, NULL, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("export_struct_rejected", "writer"); }
    int32_t xs[3] = {1, 2, 3}, ys[3] = {4, 5, 6};
    if (carquet_writer_write_batch(w, 0, xs, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, ys, 3, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        { carquet_schema_free(s); TEST_FAIL("export_struct_rejected", "write"); }
    void* buf = NULL; size_t sz = 0;
    carquet_writer_get_buffer(w, &buf, &sz);
    carquet_schema_free(s);

    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) { free(buf); TEST_FAIL("export_struct_rejected", "read"); }
    const carquet_schema_t* rs = carquet_reader_schema(r);
    struct ArrowArray aarray;
    carquet_status_t st = carquet_arrow_export_batch(batch, rs, NULL, &aarray, &err);
    assert(st == CARQUET_ERROR_NOT_IMPLEMENTED);
    free_one_batch(r, br);
    free(buf);
    TEST_PASS("export_struct_rejected");
    return 0;
}

/* ------------------------------------------------------------------ */
/* INTEGER(8|16) columns are INT32 on the Parquet side but int8/uint8/int16/
 * uint16 on the Arrow side: export must narrow the values buffer to the width
 * its format string declares, and import must widen from that width rather than
 * read the Arrow buffer with the INT32 stride. */
static carquet_schema_t* make_narrow_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    static const struct { const char* name; int bw; bool sgn; int rep; } cols[4] = {
        {"i8", 8, true, CARQUET_REPETITION_REQUIRED},
        {"u8", 8, false, CARQUET_REPETITION_REQUIRED},
        {"i16", 16, true, CARQUET_REPETITION_REQUIRED},
        {"u16", 16, false, CARQUET_REPETITION_OPTIONAL},
    };
    for (int i = 0; i < 4; i++) {
        carquet_logical_type_t lt = {0};
        lt.id = CARQUET_LOGICAL_INTEGER;
        lt.params.integer.bit_width = (int8_t)cols[i].bw;
        lt.params.integer.is_signed = cols[i].sgn;
        carquet_schema_add_column(s, cols[i].name, CARQUET_PHYSICAL_INT32, &lt,
                                  (carquet_field_repetition_t)cols[i].rep, 0, 0);
    }
    return s;
}

static const int32_t NARROW_I8[5]  = {-128, -1, 0, 1, 127};
static const int32_t NARROW_U8[5]  = {0, 1, 128, 200, 255};
static const int32_t NARROW_I16[5] = {-32768, -2, 0, 3, 32767};
static const int32_t NARROW_U16[4] = {0, 40000, 65535, 7};      /* dense */
static const int16_t NARROW_U16_DEF[5] = {1, 0, 1, 1, 1};

static void check_narrow_array(const struct ArrowSchema* sc, const struct ArrowArray* ar) {
    assert(sc->n_children == 4 && ar->n_children == 4);
    assert(strcmp(sc->children[0]->format, "c") == 0);
    assert(strcmp(sc->children[1]->format, "C") == 0);
    assert(strcmp(sc->children[2]->format, "s") == 0);
    assert(strcmp(sc->children[3]->format, "S") == 0);
    const int8_t*   a = (const int8_t*)ar->children[0]->buffers[1];
    const uint8_t*  b = (const uint8_t*)ar->children[1]->buffers[1];
    const int16_t*  c = (const int16_t*)ar->children[2]->buffers[1];
    const uint16_t* d = (const uint16_t*)ar->children[3]->buffers[1];
    const uint8_t*  dv = (const uint8_t*)ar->children[3]->buffers[0];
    for (int i = 0; i < 5; i++) {
        assert(a[i] == NARROW_I8[i]);
        assert(b[i] == NARROW_U8[i]);
        assert(c[i] == NARROW_I16[i]);
    }
    assert(dv && dv[0] == 0x1D);                    /* row 1 null */
    assert(d[0] == 0 && d[2] == 40000 && d[3] == 65535 && d[4] == 7);
}

static int test_narrow_ints(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_narrow_schema();
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) TEST_FAIL("narrow_ints", "writer create failed");
    if (carquet_writer_write_batch(w, 0, NARROW_I8, 5, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, NARROW_U8, 5, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 2, NARROW_I16, 5, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 3, NARROW_U16, 5, NARROW_U16_DEF, NULL) != CARQUET_OK)
        TEST_FAIL("narrow_ints", "write failed");
    void* buf = NULL; size_t sz = 0;
    if (carquet_writer_close(w) != CARQUET_OK ||
        carquet_writer_get_buffer(w, &buf, &sz) != CARQUET_OK)
        TEST_FAIL("narrow_ints", "close failed");
    carquet_schema_free(schema);

    /* Whole-row-group read: buffers at the declared Arrow width. */
    carquet_reader_t* r0 = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r0) TEST_FAIL("narrow_ints", "open failed");
    struct ArrowSchema rs; struct ArrowArray ra;
    if (carquet_reader_read_arrow(r0, 0, &rs, &ra, &err) != CARQUET_OK)
        TEST_FAIL("narrow_ints", err.message);
    check_narrow_array(&rs, &ra);
    rs.release(&rs); ra.release(&ra);
    carquet_reader_close(r0);

    /* Batch export: same layout. */
    carquet_reader_t* r = NULL; carquet_batch_reader_t* br = NULL; carquet_row_batch_t* batch = NULL;
    if (read_one_batch(buf, sz, &r, &br, &batch) != 0) TEST_FAIL("narrow_ints", "read failed");
    struct ArrowSchema aschema; struct ArrowArray aarray;
    if (carquet_arrow_export_batch(batch, carquet_reader_schema(r), &aschema, &aarray, &err) != CARQUET_OK)
        TEST_FAIL("narrow_ints", err.message);
    check_narrow_array(&aschema, &aarray);
    free_one_batch(r, br);
    free(buf);

    /* Import the narrow Arrow arrays (consumes them) and read INT32 back. */
    carquet_schema_t* schema2 = make_narrow_schema();
    carquet_writer_t* w2 = carquet_writer_create_buffer(schema2, NULL, &err);
    if (!w2) TEST_FAIL("narrow_ints", "writer2 create failed");
    if (carquet_writer_write_arrow(w2, &aarray, &aschema, &err) != CARQUET_OK)
        TEST_FAIL("narrow_ints", err.message);
    void* buf2 = NULL; size_t sz2 = 0;
    if (carquet_writer_close(w2) != CARQUET_OK ||
        carquet_writer_get_buffer(w2, &buf2, &sz2) != CARQUET_OK)
        TEST_FAIL("narrow_ints", "close2 failed");
    carquet_schema_free(schema2);

    carquet_reader_t* r2 = carquet_reader_open_buffer(buf2, sz2, NULL, &err);
    if (!r2) TEST_FAIL("narrow_ints", "open2 failed");
    const int32_t* expect[3] = {NARROW_I8, NARROW_U8, NARROW_I16};
    for (int col = 0; col < 3; col++) {
        carquet_column_reader_t* c = carquet_reader_get_column(r2, 0, col, NULL);
        int32_t v[5];
        assert(carquet_column_read_batch(c, v, 5, NULL, NULL) == 5);
        assert(memcmp(v, expect[col], sizeof(v)) == 0);
        carquet_column_reader_free(c);
    }
    carquet_column_reader_t* c3 = carquet_reader_get_column(r2, 0, 3, NULL);
    int32_t v3[5]; int16_t d3[5];
    assert(carquet_column_read_batch(c3, v3, 5, d3, NULL) == 5);
    assert(memcmp(d3, NARROW_U16_DEF, sizeof(d3)) == 0);
    assert(memcmp(v3, NARROW_U16, sizeof(NARROW_U16)) == 0);
    carquet_column_reader_free(c3);
    carquet_reader_close(r2);
    free(buf2);
    TEST_PASS("narrow_ints");
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_export_bytes();
    rc |= test_export_owns_memory();
    rc |= test_roundtrip();
    rc |= test_narrow_ints();
    rc |= test_import_schema();
    rc |= test_errors();
    rc |= test_write_arrow_errors();
    rc |= test_sliced_arrays();
    rc |= test_decimal_import();
    rc |= test_temporal_coercions();
    rc |= test_fixed_size_list();
    rc |= test_dictionary_null_slot();
    rc |= test_all_null();
    rc |= test_export_list();
    rc |= test_export_map();
    rc |= test_export_struct_rejected();
    if (rc == 0) printf("\nAll Arrow C Data Interface tests passed.\n");
    return rc;
}
