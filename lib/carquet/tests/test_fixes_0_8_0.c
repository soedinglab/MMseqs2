/**
 * @file test_fixes_0_8_0.c
 * @brief One regression test per v0.8.0 fix in the Arrow stream code
 *        (arrow_restore.c, arrow_c_stream.c, arrow_filter.c, arrow_ref.c).
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

static char g_dir[512];

static const char* tmp_path(const char* name) {
    static char buf[4][700];
    static int k;
    char* p = buf[k++ & 3];
    snprintf(p, 700, "%s/carquet_fix080_%s.parquet", g_dir, name);
    return p;
}

static carquet_logical_type_t string_type(void) {
    carquet_logical_type_t t;
    memset(&t, 0, sizeof(t));
    t.id = CARQUET_LOGICAL_STRING;
    return t;
}

/* A leaf ArrowSchema; `dict_format` makes it a dictionary of that value type. */
static struct ArrowSchema leaf_schema(const char* format, const char* name, const char* dict_format) {
    struct ArrowSchema out, dict;
    memset(&dict, 0, sizeof(dict));
    if (dict_format) {
        assert(carquet_arrow_schema_build(dict_format, "", ARROW_FLAG_NULLABLE, NULL, 0, NULL, NULL, &dict) == CARQUET_OK);
    }
    assert(carquet_arrow_schema_build(format, name, ARROW_FLAG_NULLABLE, NULL, 0, NULL,
                                      dict_format ? &dict : NULL, &out) == CARQUET_OK);
    return out;
}

static struct ArrowSchema struct_schema(struct ArrowSchema* kids, int32_t n) {
    struct ArrowSchema out;
    assert(carquet_arrow_schema_build("+s", "", 0, NULL, n, kids, NULL, &out) == CARQUET_OK);
    return out;
}

static bool valid_at(const struct ArrowArray* a, int64_t i) {
    const uint8_t* v = a->n_buffers > 0 ? (const uint8_t*)a->buffers[0] : NULL;
    return !v || ((v[(a->offset + i) >> 3] >> ((a->offset + i) & 7)) & 1);
}

/* Structural agreement of an array with its schema: what a consumer relies on
 * before touching a buffer. */
static void check_matches(const struct ArrowSchema* s, const struct ArrowArray* a) {
    assert((s->dictionary != NULL) == (a->dictionary != NULL));
    const char* f = s->format;
    int64_t want_buffers = -1;
    if (strcmp(f, "n") == 0) want_buffers = 0;
    else if (strcmp(f, "+s") == 0 || strncmp(f, "+w:", 3) == 0) want_buffers = 1;
    else if (strcmp(f, "u") == 0 || strcmp(f, "z") == 0 || strcmp(f, "U") == 0 || strcmp(f, "Z") == 0) want_buffers = 3;
    else if (f[0] == 'v') want_buffers = -1;
    else want_buffers = 2;
    if (want_buffers >= 0) assert(a->n_buffers == want_buffers);
    if (f[0] == '+') {
        assert(a->n_children == s->n_children);
        for (int64_t i = 0; i < s->n_children; i++) check_matches(s->children[i], a->children[i]);
    }
    if (s->dictionary) check_matches(s->dictionary, a->dictionary);
}

/* ARROW:schema says dictionary<int8, string>, the column holds 300 distinct
 * values: the batch used to come back as plain strings under a schema that
 * announces int8 indices and a dictionary. */
static int test_dictionary_index_overflow(void) {
    const char* path = tmp_path("dict_overflow");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t str = string_type();
    assert(s && carquet_schema_add_column(s, "s", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                                          CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    carquet_schema_free(s);
    assert(w);
    struct ArrowSchema kid = leaf_schema("c", "s", "u");
    struct ArrowSchema top = struct_schema(&kid, 1);
    assert(carquet_writer_set_arrow_schema(w, &top) == CARQUET_OK);
    top.release(&top);

    enum { N = 300 };
    static char text[N][8];
    carquet_byte_array_t ba[N];
    int16_t def[N];
    for (int i = 0; i < N; i++) {
        snprintf(text[i], sizeof(text[i]), "v%03d", i);
        ba[i].data = (uint8_t*)text[i];
        ba[i].length = 4;
        def[i] = 1;
    }
    assert(carquet_writer_write_batch(w, 0, ba, N, def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(w) == CARQUET_OK);

    /* Small batches fit int8 codes: dictionary arrays that agree with the schema. */
    for (int pass = 0; pass < 2; pass++) {
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL("dictionary_index_overflow", err.message);
        carquet_arrow_stream_options_t o;
        carquet_arrow_stream_options_init(&o);
        o.restore_arrow_schema = true;
        o.batch_size = pass == 0 ? 100 : 1000;
        struct ArrowArrayStream st;
        assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
        struct ArrowSchema sc;
        assert(st.get_schema(&st, &sc) == 0);
        assert(strcmp(sc.children[0]->format, "c") == 0 && sc.children[0]->dictionary);
        int64_t rows = 0;
        int failed = 0;
        for (;;) {
            struct ArrowArray b;
            int rc = st.get_next(&st, &b);
            if (rc != 0) {
                /* 300 distinct values do not fit int8: an error, not a batch
                 * of another type. */
                assert(st.get_last_error(&st) != NULL && b.release == NULL);
                failed = 1;
                break;
            }
            if (!b.release) break;
            check_matches(&sc, &b);
            const struct ArrowArray* c = b.children[0];
            for (int64_t i = 0; i < b.length; i++) {
                int8_t code = ((const int8_t*)c->buffers[1])[c->offset + i];
                const int32_t* off = (const int32_t*)c->dictionary->buffers[1];
                const char* data = (const char*)c->dictionary->buffers[2];
                assert(code >= 0 && code < c->dictionary->length);
                assert(off[code + 1] - off[code] == 4 &&
                       memcmp(data + off[code], text[rows + i], 4) == 0);
            }
            rows += b.length;
            b.release(&b);
        }
        assert(pass == 0 ? (!failed && rows == N) : failed);
        sc.release(&sc);
        st.release(&st);
        carquet_reader_close(r);
    }
    remove(path);
    TEST_PASS("dictionary_index_overflow");
    return 0;
}

/* ARROW:schema gives a decimal a bit width its precision cannot have
 * (decimal32 with 38 digits): the values used to be cut to 4 bytes. */
static int test_decimal_width_from_stored_schema(void) {
    const char* path = tmp_path("decimal_width");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t dec;
    memset(&dec, 0, sizeof(dec));
    dec.id = CARQUET_LOGICAL_DECIMAL;
    dec.params.decimal.precision = 38;
    dec.params.decimal.scale = 2;
    assert(s && carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, &dec,
                                          CARQUET_REPETITION_REQUIRED, 16, 0) == CARQUET_OK);
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    carquet_schema_free(s);
    assert(w);
    struct ArrowSchema kid = leaf_schema("d:38,2,32", "d", NULL);
    struct ArrowSchema top = struct_schema(&kid, 1);
    assert(carquet_writer_set_arrow_schema(w, &top) == CARQUET_OK);
    top.release(&top);
    uint8_t v[2][16];
    memset(v, 0, sizeof(v));
    v[0][15] = 7;                 /* 7 */
    v[1][3] = 1;                  /* 2^96: does not fit 4 bytes */
    assert(carquet_writer_write_batch(w, 0, v, 2, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_close(w) == CARQUET_OK);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("decimal_width_from_stored_schema", err.message);
    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    o.restore_arrow_schema = true;
    struct ArrowArrayStream st;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
    struct ArrowSchema sc;
    struct ArrowArray b;
    assert(st.get_schema(&st, &sc) == 0);
    assert(strcmp(sc.children[0]->format, "d:38,2") == 0);
    assert(st.get_next(&st, &b) == 0 && b.release && b.length == 2);
    const uint8_t* d = (const uint8_t*)b.children[0]->buffers[1];
    assert(d[0] == 7 && d[16 + 12] == 1);
    b.release(&b);
    sc.release(&sc);
    st.release(&st);
    carquet_reader_close(r);
    remove(path);
    TEST_PASS("decimal_width_from_stored_schema");
    return 0;
}

/* id INT32 | name STRING OPTIONAL | tag STRING OPTIONAL | pt STRUCT{x INT32}
 * (the struct sends the stream down the row-group path), 3 row groups of 6
 * rows; stored Arrow types: name large_string, tag dictionary<int16, string>.
 * A row window, a predicate and the restore together. */
static int test_window_predicate_restore(void) {
    const char* path = tmp_path("window_restore");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t str = string_type();
    int32_t pt;
    assert(s);
    assert(carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);
    assert(carquet_schema_add_column(s, "name", CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);
    assert(carquet_schema_add_column(s, "tag", CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);
    assert((pt = carquet_schema_add_group(s, "pt", CARQUET_REPETITION_REQUIRED, 0)) >= 0);
    assert(carquet_schema_add_column(s, "x", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, pt) == CARQUET_OK);
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    carquet_schema_free(s);
    assert(w);
    {
        struct ArrowSchema kids[4], x;
        kids[0] = leaf_schema("i", "id", NULL);
        kids[1] = leaf_schema("U", "name", NULL);
        kids[2] = leaf_schema("s", "tag", "u");
        x = leaf_schema("i", "x", NULL);
        assert(carquet_arrow_schema_build("+s", "pt", 0, NULL, 1, &x, NULL, &kids[3]) == CARQUET_OK);
        struct ArrowSchema top = struct_schema(kids, 4);
        assert(carquet_writer_set_arrow_schema(w, &top) == CARQUET_OK);
        top.release(&top);
    }
    enum { RG = 3, PER = 6, TOTAL = RG * PER };
    static char names[TOTAL][16];
    static const char* tags[3] = {"red", "green", "blue"};
    for (int g = 0; g < RG; g++) {
        int32_t id[PER], x[PER];
        carquet_byte_array_t nm[PER], tg[PER];
        int16_t dn[PER], dt[PER];
        int nn = 0, nt = 0;
        for (int i = 0; i < PER; i++) {
            int r = g * PER + i;
            id[i] = r;
            x[i] = r * 10;
            dn[i] = (r % 4 == 2) ? 0 : 1;      /* name null at 2, 6, 10, 14 */
            dt[i] = (r % 5 == 0) ? 0 : 1;      /* tag null at 0, 5, 10, 15 */
            if (dn[i]) {
                snprintf(names[r], sizeof(names[r]), "name-%d", r);
                nm[nn].data = (uint8_t*)names[r];
                nm[nn++].length = (int32_t)strlen(names[r]);
            }
            if (dt[i]) {
                tg[nt].data = (uint8_t*)tags[r % 3];
                tg[nt++].length = (int32_t)strlen(tags[r % 3]);
            }
        }
        if (g > 0) assert(carquet_writer_new_row_group(w) == CARQUET_OK);
        assert(carquet_writer_write_batch(w, 0, id, PER, NULL, NULL) == CARQUET_OK);
        assert(carquet_writer_write_batch(w, 1, nm, PER, dn, NULL) == CARQUET_OK);
        assert(carquet_writer_write_batch(w, 2, tg, PER, dt, NULL) == CARQUET_OK);
        assert(carquet_writer_write_batch(w, 3, x, PER, NULL, NULL) == CARQUET_OK);
    }
    assert(carquet_writer_close(w) == CARQUET_OK);

    /* with_pred 0: window only; 1: window + id >= 9 (on a field not returned
     * when `project`). */
    for (int variant = 0; variant < 4; variant++) {
        bool with_pred = variant & 1, project = variant & 2;
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL("window_predicate_restore", err.message);
        carquet_arrow_stream_options_t o;
        carquet_arrow_stream_options_init(&o);
        o.restore_arrow_schema = true;
        o.row_offset = 4;
        o.row_limit = 10;                       /* rows 4..13 */
        int32_t all[4] = {0, 1, 2, 3}, some[3] = {2, 1, 3};
        o.field_indices = project ? some : all;
        o.num_fields = project ? 3 : 4;
        int32_t nine = 9, one = 1;
        carquet_filter_clause_t cl;
        memset(&cl, 0, sizeof(cl));
        cl.column_index = 0;
        cl.op = CARQUET_FILTER_GE;
        cl.value = &nine;
        cl.value_size = 4;
        carquet_predicate_t pred = { &cl, &one, 1 };
        if (with_pred) o.predicate = &pred;
        struct ArrowArrayStream st;
        if (carquet_reader_export_arrow_stream(r, &o, &st, &err) != CARQUET_OK)
            TEST_FAIL("window_predicate_restore", err.message);
        struct ArrowSchema sc;
        assert(st.get_schema(&st, &sc) == 0);
        int ci_name = 1, ci_tag = project ? 0 : 2, ci_pt = project ? 2 : 3;
        assert(sc.n_children == (project ? 3 : 4));
        assert(strcmp(sc.children[ci_name]->format, "U") == 0);
        assert(strcmp(sc.children[ci_tag]->format, "s") == 0 && sc.children[ci_tag]->dictionary);
        int expect = with_pred ? 9 : 4;
        for (;;) {
            struct ArrowArray b;
            if (st.get_next(&st, &b) != 0) TEST_FAIL("window_predicate_restore", st.get_last_error(&st));
            if (!b.release) break;
            check_matches(&sc, &b);
            /* The window is in the fields: a batch with an offset of its own
             * cannot be imported as a record batch by Arrow C++ / PyArrow. */
            assert(b.offset == 0);
            for (int64_t i = 0; i < b.length; i++, expect++) {
                int64_t row = b.offset + i;       /* a struct's offset applies to its children */
                const struct ArrowArray* nm = b.children[ci_name];
                const struct ArrowArray* tg = b.children[ci_tag];
                const struct ArrowArray* x = b.children[ci_pt]->children[0];
                assert(row < nm->length && row < tg->length && row < b.children[ci_pt]->length);
                if (!project) {
                    const struct ArrowArray* id = b.children[0];
                    assert(((const int32_t*)id->buffers[1])[id->offset + row] == expect);
                }
                int64_t xrow = b.children[ci_pt]->offset + row;
                assert(((const int32_t*)x->buffers[1])[x->offset + xrow] == expect * 10);
                if (expect % 4 == 2) {
                    assert(!valid_at(nm, row));
                } else {
                    const int64_t* off = (const int64_t*)nm->buffers[1];
                    int64_t a = off[nm->offset + row], e = off[nm->offset + row + 1];
                    assert(valid_at(nm, row) && e - a == (int64_t)strlen(names[expect]) &&
                           memcmp((const char*)nm->buffers[2] + a, names[expect], (size_t)(e - a)) == 0);
                }
                if (expect % 5 == 0) {
                    assert(!valid_at(tg, row));
                } else {
                    int16_t code;
                    memcpy(&code, (const uint8_t*)tg->buffers[1] + 2 * (tg->offset + row), 2);
                    const struct ArrowArray* d = tg->dictionary;
                    const int32_t* off = (const int32_t*)d->buffers[1];
                    assert(valid_at(tg, row) && code >= 0 && code < d->length);
                    const char* want = tags[expect % 3];
                    assert(off[code + 1] - off[code] == (int32_t)strlen(want) &&
                           memcmp((const char*)d->buffers[2] + off[code], want, strlen(want)) == 0);
                }
            }
            b.release(&b);
        }
        assert(expect == 14);
        sc.release(&sc);
        st.release(&st);
        carquet_reader_close(r);
    }
    remove(path);
    TEST_PASS("window_predicate_restore");
    return 0;
}

/* Minimised fuzz_arrow_stream input: a footer with no column and a row group
 * of -21 rows. With no field to read there was no decoded column to check the
 * count against, and the stream returned a batch of length -21. */
static int test_negative_row_count(void) {
    static const uint8_t file[] = {
        0x50, 0x41, 0x52, 0x31, 0x15, 0x06, 0x15, 0x18, 0x2c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00,
        0x16, 0x04, 0x19, 0x1c, 0x19, 0x2c, 0x06, 0x10, 0x19, 0x18, 0x01, 0x64, 0x15, 0x00, 0x16, 0x04,
        0x16, 0x8a, 0x01, 0x8a, 0x26, 0x01, 0x2a, 0x08, 0x16, 0x26, 0x1c, 0x18, 0x03, 0x01, 0x3f, 0x00,
        0x00, 0x02, 0x00, 0x00, 0x0f, 0x16, 0x02, 0x28, 0x03, 0xff, 0x04, 0x0f, 0x18, 0x03, 0x00, 0x00,
        0x0f, 0x11, 0x11, 0x00, 0x19, 0x2c, 0x15, 0x04, 0x15, 0x00, 0x02, 0x10, 0x15, 0x00, 0x15, 0x15,
        0x00, 0x15, 0x02, 0x00, 0x3c, 0x29, 0x06, 0x19, 0x26, 0x02, 0x02, 0x00, 0x00, 0x00, 0x37, 0x02,
        0x00, 0x00, 0x50, 0x3d, 0x3d, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
        0x00, 0x50, 0x41, 0x52, 0x31, 0x65, 0x74, 0x2d, 0x63, 0x70, 0x70, 0x2d, 0x61, 0x72, 0x72, 0x6f,
        0x77, 0x20, 0x76, 0x65, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x50,
        0x41, 0x52, 0x31, 0x50, 0x41, 0x52, 0x31, 0x41, 0x41, 0x42, 0xce, 0x87, 0x00, 0x00, 0x40, 0x61,
        0x00, 0x00, 0x67, 0x41, 0x41, 0x41, 0x41, 0x55, 0x41, 0x41, 0x41, 0x41, 0x42, 0x41, 0x41, 0x41,
        0x41, 0x49, 0x41, 0x41, 0x41, 0x41, 0x41, 0x3d, 0x3d, 0x00, 0x18, 0x20, 0x70, 0x61, 0x72, 0x71,
        0x75, 0x65, 0x74, 0x2d, 0x63, 0x70, 0x70, 0x2d, 0x61, 0x72, 0x72, 0x6f, 0x77, 0x20, 0x76, 0x65,
        0x72, 0x73, 0x69, 0x6f, 0x6e, 0x20, 0x32, 0x36, 0x2e, 0x30, 0x2e, 0x3d, 0x00, 0x02, 0x00, 0x00,
        0x50, 0x41, 0x52, 0x31, 0x2e, 0x30, 0x2e, 0x31, 0x19, 0x2c, 0x1c, 0x00, 0x00, 0x1c, 0x00, 0x00,
        0x00, 0x37, 0x02, 0x00, 0x00, 0x50, 0x3d, 0x3d, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x01, 0x00, 0x00, 0x50, 0x41, 0x52, 0x31,
    };
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(file, sizeof(file), NULL, &err);
    if (r) {
        struct ArrowArrayStream st;
        if (carquet_reader_export_arrow_stream(r, NULL, &st, &err) == CARQUET_OK) {
            for (;;) {
                struct ArrowArray b;
                if (st.get_next(&st, &b) != 0 || !b.release) break;
                assert(b.length >= 0);
                b.release(&b);
            }
            st.release(&st);
        }
        for (int32_t g = 0; g < carquet_reader_num_row_groups(r); g++) {
            struct ArrowArray a;
            if (carquet_reader_read_arrow(r, g, NULL, &a, &err) == CARQUET_OK) {
                assert(a.length >= 0);
                a.release(&a);
            }
        }
        carquet_reader_close(r);
    }
    TEST_PASS("negative_row_count");
    return 0;
}

/* f LIST<INT32> OPTIONAL (optional elements), stored Arrow type
 * fixed_size_list<int32>[size]. `spec` gives one character per row: 'n' a null
 * list, 'f' a full one (elements 10*row + k, the odd rows' last one null). */
static const char* write_fsl_file(const char* name, const char* spec, const char* fsl_format, int size) {
    const char* path = tmp_path(name);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    assert(s);
    int32_t inner = carquet_schema_add_list_group(s, "f", CARQUET_REPETITION_OPTIONAL, 0);
    assert(inner >= 0);
    assert(carquet_schema_add_column(s, "element", CARQUET_PHYSICAL_INT32, NULL,
                                     CARQUET_REPETITION_OPTIONAL, 0, inner) == CARQUET_OK);
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    carquet_schema_free(s);
    assert(w);
    struct ArrowSchema item = leaf_schema("i", "item", NULL), f, top;
    assert(carquet_arrow_schema_build(fsl_format, "f", ARROW_FLAG_NULLABLE, NULL, 1, &item, NULL, &f) == CARQUET_OK);
    top = struct_schema(&f, 1);
    assert(carquet_writer_set_arrow_schema(w, &top) == CARQUET_OK);
    top.release(&top);

    int32_t values[64];
    int16_t def[64], rep[64];
    int nl = 0, nv = 0;
    for (int r = 0; spec[r]; r++) {
        if (spec[r] == 'n') { def[nl] = 0; rep[nl++] = 0; continue; }
        for (int k = 0; k < size; k++) {
            bool null_elem = (r % 2 == 1) && k == size - 1;
            rep[nl] = k ? 1 : 0;
            def[nl++] = null_elem ? 2 : 3;
            if (!null_elem) values[nv++] = 10 * r + k;
        }
    }
    assert(carquet_writer_write_batch(w, 0, values, nl, def, rep) == CARQUET_OK);
    assert(carquet_writer_close(w) == CARQUET_OK);
    return path;
}

/* Restoring a fixed-size list: full lists are served from the list's own
 * child (no copy), null lists stored empty get made-up child slots, and a
 * stored size that would need gigabytes of them is refused instead of
 * allocated (fuzz_arrow_stream: 3 null lists of size 2^30 asked for 25 GB). */
static int test_fixed_size_list_restore(void) {
    static const struct { const char* spec; const char* format; int size; bool fails; } cases[] = {
        {"ffff", "+w:2", 2, false},
        {"fnff", "+w:2", 2, false},
        {"nfnnf", "+w:3", 3, false},
        {"nnn", "+w:2", 2, false},
        {"nnn", "+w:1073741824", 2, true},
    };
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
        const char* path = write_fsl_file("fsl", cases[c].spec, cases[c].format, cases[c].size);
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL("fixed_size_list_restore", err.message);
        carquet_arrow_stream_options_t o;
        carquet_arrow_stream_options_init(&o);
        o.restore_arrow_schema = true;
        struct ArrowArrayStream st;
        if (carquet_reader_export_arrow_stream(r, &o, &st, &err) != CARQUET_OK)
            TEST_FAIL("fixed_size_list_restore", err.message);
        struct ArrowSchema sc;
        struct ArrowArray b;
        assert(st.get_schema(&st, &sc) == 0);
        assert(strcmp(sc.children[0]->format, cases[c].format) == 0);
        int rc = st.get_next(&st, &b);
        if (cases[c].fails) {
            assert(rc != 0 && b.release == NULL && st.get_last_error(&st) != NULL);
        } else {
            int size = cases[c].size;
            assert(rc == 0 && b.release && b.length == (int64_t)strlen(cases[c].spec));
            check_matches(&sc, &b);
            const struct ArrowArray* f = b.children[0];
            const struct ArrowArray* item = f->children[0];
            assert(item->length >= (f->offset + f->length) * size);
            for (int64_t row = 0; row < b.length; row++) {
                if (cases[c].spec[row] == 'n') { assert(!valid_at(f, row)); continue; }
                assert(valid_at(f, row));
                for (int k = 0; k < size; k++) {
                    int64_t slot = (f->offset + row) * size + k;
                    if ((row % 2 == 1) && k == size - 1) {
                        assert(!valid_at(item, slot));
                    } else {
                        assert(valid_at(item, slot));
                        assert(((const int32_t*)item->buffers[1])[item->offset + slot] == 10 * (int32_t)row + k);
                    }
                }
            }
            b.release(&b);
        }
        sc.release(&sc);
        st.release(&st);
        carquet_reader_close(r);
        remove(path);
    }
    TEST_PASS("fixed_size_list_restore");
    return 0;
}

/* ARROW:schema says the null type for an INT32 column (how files from before
 * the Null annotation carry it). All-null: restored. With a value in it: the
 * value used to be dropped without a word. */
static int test_null_type_over_values(void) {
    for (int with_value = 0; with_value < 2; with_value++) {
        const char* path = tmp_path("null_type");
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* s = carquet_schema_create(&err);
        assert(s && carquet_schema_add_column(s, "n", CARQUET_PHYSICAL_INT32, NULL,
                                              CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);
        carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
        carquet_schema_free(s);
        assert(w);
        struct ArrowSchema kid = leaf_schema("n", "n", NULL);
        struct ArrowSchema top = struct_schema(&kid, 1);
        assert(carquet_writer_set_arrow_schema(w, &top) == CARQUET_OK);
        top.release(&top);
        int32_t v[1] = {42};
        int16_t def[3] = {0, (int16_t)with_value, 0};
        assert(carquet_writer_write_batch(w, 0, v, 3, def, NULL) == CARQUET_OK);
        assert(carquet_writer_close(w) == CARQUET_OK);

        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) TEST_FAIL("null_type_over_values", err.message);
        carquet_arrow_stream_options_t o;
        carquet_arrow_stream_options_init(&o);
        o.restore_arrow_schema = true;
        struct ArrowArrayStream st;
        assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
        struct ArrowSchema sc;
        struct ArrowArray b;
        assert(st.get_schema(&st, &sc) == 0 && strcmp(sc.children[0]->format, "n") == 0);
        int rc = st.get_next(&st, &b);
        if (with_value) {
            assert(rc != 0 && b.release == NULL);
        } else {
            assert(rc == 0 && b.release && b.length == 3);
            check_matches(&sc, &b);
            assert(b.children[0]->null_count == 3);
            b.release(&b);
        }
        sc.release(&sc);
        st.release(&st);
        carquet_reader_close(r);
        remove(path);
    }
    TEST_PASS("null_type_over_values");
    return 0;
}

int main(void) {
    assert(carquet_init() == CARQUET_OK);
    const char* t = getenv("TMPDIR");
#ifdef _WIN32
    if (!t) t = getenv("TEMP");
    if (!t) t = ".";
#else
    if (!t) t = "/tmp";
#endif
    snprintf(g_dir, sizeof(g_dir), "%s", t);
    size_t n = strlen(g_dir);
    if (n > 1 && (g_dir[n - 1] == '/' || g_dir[n - 1] == '\\')) g_dir[n - 1] = '\0';

    int failures = 0;
    failures += test_dictionary_index_overflow();
    failures += test_decimal_width_from_stored_schema();
    failures += test_window_predicate_restore();
    failures += test_negative_row_count();
    failures += test_fixed_size_list_restore();
    failures += test_null_type_over_values();
    if (failures) {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
    printf("All v0.8.0 fix tests passed\n");
    return 0;
}
