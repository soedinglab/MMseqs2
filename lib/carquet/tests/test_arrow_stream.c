/**
 * @file test_arrow_stream.c
 * @brief carquet_reader_export_arrow_stream() and the pyarrow-compatible
 *        Arrow types of DECIMAL / INT96 (arrow_leaf.c).
 *
 *  - stream_flat            batch path: batch sizes, row groups not spanned,
 *                           schema == batch types, values, arrays outliving
 *                           the stream and the reader
 *  - stream_decimal_int96   DECIMAL over INT32 / INT64 / FLBA / BYTE_ARRAY as
 *                           decimal128 / decimal256, INT96 as timestamp[ns],
 *                           null slots, on both the batch and nested paths
 *  - stream_projection      selected fields reordered, nested (STRUCT) field
 *                           -> one batch per row group
 *  - stream_errors          bad / duplicate field index, BYTE_ARRAY decimal
 *                           too wide for its precision (get_next error +
 *                           get_last_error), close_reader
 */
#include <carquet/carquet.h>

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_helpers.h"

#define ROWS_PER_RG 5
#define NUM_RG 2

static carquet_logical_type_t decimal(int32_t p, int32_t s) {
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_DECIMAL;
    lt.params.decimal.precision = p;
    lt.params.decimal.scale = s;
    return lt;
}

/* Little-endian 128/256-bit two's complement value of a decimal slot, as int64
 * (the tests only use values that fit), checking the sign extension. */
static int64_t dec_at(const struct ArrowArray* a, size_t width, int64_t i) {
    const uint8_t* p = (const uint8_t*)a->buffers[1] + (size_t)(a->offset + i) * width;
    int64_t v;
    memcpy(&v, p, 8);
    for (size_t k = 8; k < width; k++) assert(p[k] == (v < 0 ? 0xFF : 0x00));
    return v;
}

static bool is_valid(const struct ArrowArray* a, int64_t i) {
    const uint8_t* v = (const uint8_t*)a->buffers[0];
    int64_t j = a->offset + i;
    return !v || (v[j >> 3] >> (j & 7)) & 1;
}

/* id INT32 | d32 DECIMAL(9,2) INT32 OPTIONAL | d64 DECIMAL(18,3) INT64 |
 * dfl DECIMAL(10,2) FLBA(5) | dba DECIMAL(40,0) BYTE_ARRAY | t96 INT96 |
 * name STRING OPTIONAL | pt STRUCT{x INT32} */
static int write_file(const char* path, bool too_wide) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t d32 = decimal(9, 2), d64 = decimal(18, 3), dfl = decimal(10, 2);
    carquet_logical_type_t dba = decimal(40, 0), str;
    memset(&str, 0, sizeof(str));
    str.id = CARQUET_LOGICAL_STRING;
    int32_t pt;
    if (!s ||
        carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) ||
        carquet_schema_add_column(s, "d32", CARQUET_PHYSICAL_INT32, &d32, CARQUET_REPETITION_OPTIONAL, 0, 0) ||
        carquet_schema_add_column(s, "d64", CARQUET_PHYSICAL_INT64, &d64, CARQUET_REPETITION_REQUIRED, 0, 0) ||
        carquet_schema_add_column(s, "dfl", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, &dfl,
                                  CARQUET_REPETITION_REQUIRED, 5, 0) ||
        carquet_schema_add_column(s, "dba", CARQUET_PHYSICAL_BYTE_ARRAY, &dba, CARQUET_REPETITION_REQUIRED, 0, 0) ||
        carquet_schema_add_column(s, "t96", CARQUET_PHYSICAL_INT96, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) ||
        carquet_schema_add_column(s, "name", CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_OPTIONAL, 0, 0) ||
        (pt = carquet_schema_add_group(s, "pt", CARQUET_REPETITION_REQUIRED, 0)) < 0 ||
        carquet_schema_add_column(s, "x", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, pt)) {
        if (s) carquet_schema_free(s);
        return 1;
    }
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    carquet_schema_free(s);
    if (!w) return 1;

    int rc = 0;
    for (int g = 0; g < NUM_RG && rc == 0; g++) {
        int32_t id[ROWS_PER_RG], x[ROWS_PER_RG], d32v[ROWS_PER_RG];
        int64_t d64v[ROWS_PER_RG];
        uint8_t flba[ROWS_PER_RG][5];
        uint8_t bytes[ROWS_PER_RG][8];
        carquet_byte_array_t ba[ROWS_PER_RG], names[ROWS_PER_RG];
        carquet_int96_t t96[ROWS_PER_RG];
        int16_t def[ROWS_PER_RG];
        int nd = 0, nn = 0;
        static const char* words[ROWS_PER_RG] = {"a", "bb", "ccc", "dddd", "eeeee"};
        for (int i = 0; i < ROWS_PER_RG; i++) {
            int r = g * ROWS_PER_RG + i;
            int64_t sv = (r % 2) ? -(int64_t)(r * 1000 + 7) : (int64_t)(r * 1000 + 7);
            id[i] = r;
            x[i] = r * 10;
            def[i] = (r % 3 == 1) ? 0 : 1;          /* rows 1, 4, 7 null */
            if (def[i]) {
                d32v[nd++] = (int32_t)sv;
                names[nn].data = (uint8_t*)words[i];
                names[nn++].length = (int32_t)strlen(words[i]);
            }
            d64v[i] = sv * 3;
            uint64_t u = (uint64_t)sv;
            for (int k = 0; k < 5; k++) flba[i][k] = (uint8_t)(u >> (8 * (4 - k)));
            /* dba: the low 1..8 bytes of sv, big-endian (see dba_expected) */
            int len = 1 + (r % 8);
            for (int k = 0; k < len; k++) bytes[i][k] = (uint8_t)(u >> (8 * (len - 1 - k)));
            ba[i].data = bytes[i];
            ba[i].length = len;
            t96[i].value[0] = (uint32_t)(r * 1000);          /* nanos of day, low */
            t96[i].value[1] = 0;
            t96[i].value[2] = (uint32_t)(2440588 + r - 3);   /* Julian day */
        }
        static uint8_t wide[33];
        if (too_wide && g == 1) {
            /* 33 significant bytes: more than decimal256 holds */
            memset(wide, 0x11, sizeof(wide));
            ba[2].data = wide;
            ba[2].length = (int32_t)sizeof(wide);
        }
        if (g > 0 && carquet_writer_new_row_group(w) != CARQUET_OK) rc = 1;
        if (rc == 0 &&
            (carquet_writer_write_batch(w, 0, id, ROWS_PER_RG, NULL, NULL) ||
             carquet_writer_write_batch(w, 1, d32v, ROWS_PER_RG, def, NULL) ||
             carquet_writer_write_batch(w, 2, d64v, ROWS_PER_RG, NULL, NULL) ||
             carquet_writer_write_batch(w, 3, flba, ROWS_PER_RG, NULL, NULL) ||
             carquet_writer_write_batch(w, 4, ba, ROWS_PER_RG, NULL, NULL) ||
             carquet_writer_write_batch(w, 5, t96, ROWS_PER_RG, NULL, NULL) ||
             carquet_writer_write_batch(w, 6, names, ROWS_PER_RG, def, NULL) ||
             carquet_writer_write_batch(w, 7, x, ROWS_PER_RG, NULL, NULL)))
            rc = 1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    return rc;
}

/* The value `dba` holds for row r: its low (1 + r % 8) bytes, sign-extended. */
static int64_t dba_expected(int r) {
    int64_t sv = (r % 2) ? -(int64_t)(r * 1000 + 7) : (int64_t)(r * 1000 + 7);
    int len = 1 + (r % 8);
    if (len >= 8) return sv;
    uint64_t u = (uint64_t)sv & ((1ULL << (8 * len)) - 1);
    if (u >> (8 * len - 1)) u |= ~((1ULL << (8 * len)) - 1);   /* sign-extend */
    return (int64_t)u;
}

/* Check the converted columns of `b` (fields id, d32, d64, dfl, dba, t96 at
 * children 0..5), rows starting at file row `first`. */
static void check_converted(const struct ArrowArray* b, int first) {
    for (int64_t i = 0; i < b->length; i++) {
        int r = first + (int)i;
        int64_t sv = (r % 2) ? -(int64_t)(r * 1000 + 7) : (int64_t)(r * 1000 + 7);
        assert(((const int32_t*)b->children[0]->buffers[1])[b->children[0]->offset + i] == r);
        if (r % 3 == 1) {
            assert(!is_valid(b->children[1], i));
        } else {
            assert(is_valid(b->children[1], i) && dec_at(b->children[1], 16, i) == sv);
        }
        assert(dec_at(b->children[2], 16, i) == sv * 3);
        /* 5-byte FLBA keeps the low 40 bits, sign-extended */
        int64_t f = (int64_t)((uint64_t)sv << 24) >> 24;
        assert(dec_at(b->children[3], 16, i) == f);
        assert(dec_at(b->children[4], 32, i) == dba_expected(r));
        int64_t ns;
        memcpy(&ns, (const uint8_t*)b->children[5]->buffers[1] + (size_t)(b->children[5]->offset + i) * 8, 8);
        assert(ns == (int64_t)(r - 3) * 86400000000000LL + r * 1000);
    }
}

static int test_stream_flat_and_types(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("stream_flat", err.message);

    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    int32_t fields[7] = {0, 1, 2, 3, 4, 5, 6};   /* everything but the struct: batch path */
    o.field_indices = fields;
    o.num_fields = 7;
    o.batch_size = 3;
    struct ArrowArrayStream st;
    if (carquet_reader_export_arrow_stream(r, &o, &st, &err) != CARQUET_OK) TEST_FAIL("stream_flat", err.message);

    struct ArrowSchema sc;
    assert(st.get_schema(&st, &sc) == 0);
    static const char* fmt[7] = {"i", "d:9,2", "d:18,3", "d:10,2", "d:40,0,256", "tsn:", "u"};
    assert(strcmp(sc.format, "+s") == 0 && sc.n_children == 7);
    for (int i = 0; i < 7; i++) assert(strcmp(sc.children[i]->format, fmt[i]) == 0);
    assert(sc.children[1]->flags & ARROW_FLAG_NULLABLE);
    sc.release(&sc);

    /* 5 rows per row group, batch_size 3: 3, 2, 3, 2 (never across groups). */
    struct ArrowArray kept[4];
    int nb = 0, row = 0;
    static const int64_t expect_len[4] = {3, 2, 3, 2};
    for (;;) {
        struct ArrowArray b;
        assert(st.get_next(&st, &b) == 0);
        if (!b.release) break;
        assert(nb < 4 && b.length == expect_len[nb] && b.n_children == 7);
        check_converted(&b, row);
        row += (int)b.length;
        kept[nb++] = b;
    }
    assert(nb == 4 && row == NUM_RG * ROWS_PER_RG);
    st.release(&st);
    assert(st.release == NULL);
    carquet_reader_close(r);

    /* Arrays stay valid after the stream and the reader are gone. */
    row = 0;
    for (int i = 0; i < nb; i++) {
        check_converted(&kept[i], row);
        const int32_t* off = (const int32_t*)kept[i].children[6]->buffers[1];
        assert(off[kept[i].length] >= 0);
        row += (int)kept[i].length;
        kept[i].release(&kept[i]);
    }

    /* The nested path (read_arrow) types and converts the same way. */
    r = carquet_reader_open(path, NULL, &err);
    struct ArrowSchema ns; struct ArrowArray na;
    assert(r && carquet_reader_read_arrow(r, 1, &ns, &na, &err) == CARQUET_OK);
    for (int i = 0; i < 7; i++) assert(strcmp(ns.children[i]->format, fmt[i]) == 0);
    check_converted(&na, ROWS_PER_RG);
    ns.release(&ns); na.release(&na);
    carquet_reader_close(r);

    TEST_PASS("stream_flat_and_types");
    return 0;
}

static int test_stream_projection(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("stream_projection", err.message);

    /* Batch path, reordered: name, id. */
    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    int32_t f1[2] = {6, 0};
    o.field_indices = f1; o.num_fields = 2;
    struct ArrowArrayStream st;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
    struct ArrowSchema sc;
    assert(st.get_schema(&st, &sc) == 0 && sc.n_children == 2);
    assert(strcmp(sc.children[0]->name, "name") == 0 && strcmp(sc.children[1]->name, "id") == 0);
    sc.release(&sc);
    struct ArrowArray b;
    assert(st.get_next(&st, &b) == 0 && b.release && b.length == ROWS_PER_RG && b.n_children == 2);
    assert(((const int32_t*)b.children[1]->buffers[1])[2] == 2);
    assert(!is_valid(b.children[0], 1) && is_valid(b.children[0], 2));
    b.release(&b);
    st.release(&st);

    /* A STRUCT field: one batch per row group, whatever batch_size. */
    int32_t f2[2] = {7, 0};
    o.field_indices = f2; o.num_fields = 2; o.batch_size = 2;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
    assert(st.get_schema(&st, &sc) == 0 && sc.n_children == 2);
    assert(strcmp(sc.children[0]->format, "+s") == 0 && strcmp(sc.children[1]->name, "id") == 0);
    sc.release(&sc);
    int nb = 0;
    for (;;) {
        assert(st.get_next(&st, &b) == 0);
        if (!b.release) break;
        assert(b.length == ROWS_PER_RG && b.n_children == 2);
        const struct ArrowArray* x = b.children[0]->children[0];
        assert(((const int32_t*)x->buffers[1])[1] == (nb * ROWS_PER_RG + 1) * 10);
        b.release(&b);
        nb++;
    }
    assert(nb == NUM_RG);
    st.release(&st);

    /* Every field (options NULL): the struct sends it down the nested path. */
    assert(carquet_reader_export_arrow_stream(r, NULL, &st, &err) == CARQUET_OK);
    assert(st.get_schema(&st, &sc) == 0 && sc.n_children == 8);
    sc.release(&sc);
    assert(st.get_next(&st, &b) == 0 && b.release && b.n_children == 8);
    check_converted(&b, 0);
    b.release(&b);
    st.release(&st);

    carquet_reader_close(r);
    TEST_PASS("stream_projection");
    return 0;
}

static int test_stream_errors(const char* path, const char* wide_path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("stream_errors", err.message);
    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    struct ArrowArrayStream st;

    int32_t bad[1] = {8};
    o.field_indices = bad; o.num_fields = 1;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    assert(st.release == NULL);
    int32_t dup[2] = {0, 0};
    o.field_indices = dup; o.num_fields = 2;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    assert(st.release == NULL);
    carquet_reader_close(r);   /* still the caller's after a failed export */

    /* 33 significant bytes do not fit decimal256: an error, not a wrong value. */
    for (int nested = 0; nested < 2; nested++) {
        r = carquet_reader_open(wide_path, NULL, &err);
        if (!r) TEST_FAIL("stream_errors", err.message);
        carquet_arrow_stream_options_init(&o);
        int32_t f[2] = {4, 7};
        o.field_indices = f; o.num_fields = nested ? 2 : 1;
        o.close_reader = true;     /* the stream closes it */
        assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_OK);
        struct ArrowArray b;
        assert(st.get_next(&st, &b) == 0 && b.release);   /* row group 0 is fine */
        b.release(&b);
        int got = 0;
        for (;;) {
            got = st.get_next(&st, &b);
            if (got != 0 || !b.release) break;
            b.release(&b);
        }
        assert(got == EIO);
        const char* msg = st.get_last_error(&st);
        assert(msg && strlen(msg) > 0);
        st.release(&st);
    }

    TEST_PASS("stream_errors");
    return 0;
}

/* ---- row selection: row groups, window, predicate ------------------------ */

/* Ids (field "id" at child `idcol`) of every row the stream returns. */
static int collect_ids(carquet_reader_t* r, const carquet_arrow_stream_options_t* o, int idcol,
                       int32_t* ids, int* n, int32_t* nchildren) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    struct ArrowArrayStream st;
    *n = 0;
    if (carquet_reader_export_arrow_stream(r, o, &st, &err) != CARQUET_OK) return 1;
    for (;;) {
        struct ArrowArray b;
        if (st.get_next(&st, &b) != 0) { st.release(&st); return 1; }
        if (!b.release) break;
        assert(b.length > 0);
        if (nchildren) *nchildren = (int32_t)b.n_children;
        const struct ArrowArray* c = b.children[idcol];
        for (int64_t i = 0; i < b.length; i++)
            ids[(*n)++] = ((const int32_t*)c->buffers[1])[b.offset + c->offset + i];
        b.release(&b);
    }
    st.release(&st);
    return 0;
}

static bool same(const int32_t* got, int n, const int32_t* want, int m) {
    if (n != m) return false;
    for (int i = 0; i < n; i++) if (got[i] != want[i]) return false;
    return true;
}

static carquet_filter_clause_t clause(int32_t col, carquet_filter_op_t op, const void* v, int32_t size) {
    carquet_filter_clause_t c;
    memset(&c, 0, sizeof(c));
    c.column_index = col;
    c.op = op;
    c.value = v;
    c.value_size = size;
    return c;
}

static int test_stream_selection(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("stream_selection", err.message);
    int32_t ids[64]; int n = 0; int32_t nch = 0;
    carquet_arrow_stream_options_t o;

    /* Row groups in the order asked, both paths. */
    for (int nested = 0; nested < 2; nested++) {
        carquet_arrow_stream_options_init(&o);
        int32_t rgs[2] = {1, 0};
        int32_t f[2] = {0, 7};
        o.row_groups = rgs; o.num_row_groups = 2;
        o.field_indices = f; o.num_fields = nested ? 2 : 1;
        assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
        const int32_t want[10] = {5, 6, 7, 8, 9, 0, 1, 2, 3, 4};
        assert(same(ids, n, want, 10));
    }
    /* A window, both paths, across the row-group boundary. */
    for (int nested = 0; nested < 2; nested++) {
        carquet_arrow_stream_options_init(&o);
        int32_t f[2] = {0, 7};
        o.field_indices = f; o.num_fields = nested ? 2 : 1;
        o.row_offset = 3; o.row_limit = 4;
        assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
        const int32_t want[4] = {3, 4, 5, 6};
        assert(same(ids, n, want, 4));
    }
    int32_t bad_rgs[2] = {0, 0};
    carquet_arrow_stream_options_init(&o);
    o.row_groups = bad_rgs; o.num_row_groups = 2;
    struct ArrowArrayStream st;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    o.num_row_groups = 1; o.row_limit = 3;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);

    /* Predicates: output field id only, the predicate's fields read and dropped. */
    int32_t only_id[1] = {0};
    int32_t v7 = 7;
    carquet_filter_clause_t c1 = clause(0, CARQUET_FILTER_GE, &v7, 4);
    int32_t one = 1;
    carquet_predicate_t p = { &c1, &one, 1 };
    carquet_arrow_stream_options_init(&o);
    o.field_indices = only_id; o.num_fields = 1; o.predicate = &p;
    assert(collect_ids(r, &o, 0, ids, &n, &nch) == 0);
    { const int32_t w[3] = {7, 8, 9}; assert(same(ids, n, w, 3)); }

    /* IS NULL on a string; a negative decimal (nulls never match). */
    carquet_filter_clause_t cn = clause(6, CARQUET_FILTER_IS_NULL, NULL, 0);
    p.clauses = &cn;
    assert(collect_ids(r, &o, 0, ids, &n, &nch) == 0 && nch == 1);
    { const int32_t w[3] = {1, 4, 7}; assert(same(ids, n, w, 3)); }
    int32_t zero = 0;
    carquet_filter_clause_t cd = clause(1, CARQUET_FILTER_LT, &zero, 4);
    p.clauses = &cd;
    assert(collect_ids(r, &o, 0, ids, &n, &nch) == 0 && nch == 1);
    { const int32_t w[3] = {3, 5, 9}; assert(same(ids, n, w, 3)); }

    /* OR: id == 0, or name == "ccc" (row 2; row 7's name is null). */
    carquet_filter_clause_t dnf[2] = { clause(0, CARQUET_FILTER_EQ, &zero, 4),
                                       clause(6, CARQUET_FILTER_EQ, "ccc", 3) };
    int32_t sizes[2] = {1, 1};
    carquet_predicate_t p2 = { dnf, sizes, 2 };
    o.predicate = &p2;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[2] = {0, 2}; assert(same(ids, n, w, 2)); }

    /* IN / NOT IN. */
    carquet_byte_array_t names[2] = { {(uint8_t*)"a", 1}, {(uint8_t*)"eeeee", 5} };
    carquet_filter_clause_t cin = clause(6, CARQUET_FILTER_IN, NULL, 0);
    cin.values = names; cin.value_count = 2;
    p.clauses = &cin; o.predicate = &p;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[3] = {0, 5, 9}; assert(same(ids, n, w, 3)); }
    int32_t small[3] = {0, 1, 2};
    carquet_filter_clause_t cnin = clause(0, CARQUET_FILTER_NOT_IN, NULL, 0);
    cnin.values = small; cnin.value_count = 3;
    p.clauses = &cnin;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[7] = {3, 4, 5, 6, 7, 8, 9}; assert(same(ids, n, w, 7)); }

    /* INT96 compared by instant: after row 5's. */
    carquet_int96_t t5;
    t5.value[0] = 5000; t5.value[1] = 0; t5.value[2] = (uint32_t)(2440588 + 5 - 3);
    carquet_filter_clause_t ct = clause(5, CARQUET_FILTER_GT, &t5, 12);
    p.clauses = &ct;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[4] = {6, 7, 8, 9}; assert(same(ids, n, w, 4)); }

    /* BYTE_ARRAY decimal: row 6 holds 6007 in 7 bytes; give it in 2. */
    uint8_t be6007[2] = {0x17, 0x77};
    carquet_filter_clause_t cb = clause(4, CARQUET_FILTER_EQ, be6007, 2);
    p.clauses = &cb;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[1] = {6}; assert(same(ids, n, w, 1)); }

    /* A struct leaf (pt.x): the row-group path. */
    int32_t fifty = 50;
    carquet_filter_clause_t cx = clause(7, CARQUET_FILTER_GT, &fifty, 4);
    p.clauses = &cx;
    assert(collect_ids(r, &o, 0, ids, &n, &nch) == 0 && nch == 1);
    { const int32_t w[4] = {6, 7, 8, 9}; assert(same(ids, n, w, 4)); }

    /* Statistics prune the first row group entirely; a window still applies. */
    o.row_offset = 0; o.row_limit = 8;
    p.clauses = &c1;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t w[1] = {7}; assert(same(ids, n, w, 1)); }

    /* Malformed predicates. */
    carquet_arrow_stream_options_init(&o);
    carquet_filter_clause_t cbad = clause(99, CARQUET_FILTER_EQ, &zero, 4);
    p.clauses = &cbad; o.predicate = &p;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);
    int32_t empty = 0;
    carquet_predicate_t p0 = { &c1, &empty, 1 };
    o.predicate = &p0;
    assert(carquet_reader_export_arrow_stream(r, &o, &st, &err) == CARQUET_ERROR_INVALID_ARGUMENT);

    carquet_reader_close(r);
    TEST_PASS("stream_selection");
    return 0;
}

/* Floating point semantics and the page filter, on a DOUBLE column. */
static int test_stream_float_predicates(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "arrow_stream_f");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    assert(s && carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);
    assert(carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.write_page_index = true;
    wo.max_rows_per_page = 2;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    carquet_schema_free(s);
    assert(w);
    /* rows: 1.0, NaN, -0.0, 3.0, null */
    int32_t id[5] = {0, 1, 2, 3, 4};
    double v[4] = {1.0, NAN, -0.0, 3.0};
    int16_t def[5] = {1, 1, 1, 1, 0};
    assert(carquet_writer_write_batch(w, 0, id, 5, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(w, 1, v, 5, def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(w) == CARQUET_OK);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    assert(r);
    int32_t ids[8]; int n;
    int32_t one = 1, f0[1] = {0};
    double d1 = 1.0, d3 = 3.0, dz = 0.0, dnan = NAN;
    carquet_filter_clause_t c;
    carquet_predicate_t p = { &c, &one, 1 };
    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    o.field_indices = f0; o.num_fields = 1; o.predicate = &p;

    c = clause(1, CARQUET_FILTER_EQ, &d1, 8);
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[1] = {0}; assert(same(ids, n, want, 1)); }
    c = clause(1, CARQUET_FILTER_NE, &d1, 8);           /* NaN != 1 holds, null does not */
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[3] = {1, 2, 3}; assert(same(ids, n, want, 3)); }
    c = clause(1, CARQUET_FILTER_LT, &d3, 8);
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[2] = {0, 2}; assert(same(ids, n, want, 2)); }
    c = clause(1, CARQUET_FILTER_EQ, &dz, 8);           /* -0.0 == 0.0 */
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[1] = {2}; assert(same(ids, n, want, 1)); }
    c = clause(1, CARQUET_FILTER_GT, &dnan, 8);
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0 && n == 0);
    c = clause(1, CARQUET_FILTER_IN, NULL, 0);
    c.values = &dnan; c.value_count = 1;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[1] = {1}; assert(same(ids, n, want, 1)); }
    c = clause(1, CARQUET_FILTER_NOT_IN, NULL, 0);
    c.values = &d1; c.value_count = 1;
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[3] = {1, 2, 3}; assert(same(ids, n, want, 3)); }
    /* Page index: id >= 3 prunes the first page, the exact filter does the rest. */
    int32_t three = 3;
    c = clause(0, CARQUET_FILTER_GE, &three, 4);
    assert(collect_ids(r, &o, 0, ids, &n, NULL) == 0);
    { const int32_t want[2] = {3, 4}; assert(same(ids, n, want, 2)); }

    carquet_reader_close(r);
    carquet_test_cleanup(path);
    TEST_PASS("stream_float_predicates");
    return 0;
}

/* carquet_arrow_ref_*: several independent exports, slices, any release order. */
static int test_arrow_ref(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    struct ArrowArray a;
    assert(r && carquet_reader_read_arrow(r, 1, NULL, &a, &err) == CARQUET_OK);
    carquet_reader_close(r);
    carquet_arrow_ref_t* ref = carquet_arrow_ref_create(&a);
    assert(ref && a.release == NULL);
    const struct ArrowArray* held = carquet_arrow_ref_array(ref);
    assert(held->length == ROWS_PER_RG);
    struct ArrowArray e1, e2;
    assert(carquet_arrow_ref_export(ref, 0, -1, &e1) == CARQUET_OK);
    assert(carquet_arrow_ref_export(ref, 2, 2, &e2) == CARQUET_OK);
    struct ArrowArray e3;
    assert(carquet_arrow_ref_export(ref, 6, 1, &e3) == CARQUET_ERROR_INVALID_ARGUMENT);
    carquet_arrow_ref_release(ref);            /* exports keep the data */
    /* A record batch takes the slice in its children (offset 0 itself). */
    assert(e2.length == 2 && e2.offset == 0 && e2.children[0]->offset == 2);
    const struct ArrowArray* id = e2.children[0];
    assert(((const int32_t*)id->buffers[1])[e2.offset + id->offset] == 7);
    /* Move a child out, release the parent, then the child. */
    struct ArrowArray moved = *e1.children[6];
    e1.children[6]->release = NULL;
    e1.release(&e1);
    const int32_t* off = (const int32_t*)moved.buffers[1];
    assert(off[moved.length] >= off[0]);
    moved.release(&moved);
    e2.release(&e2);

    /* Schema copy. */
    r = carquet_reader_open(path, NULL, &err);
    struct ArrowSchema sch, cp;
    assert(r && carquet_arrow_export_schema(carquet_reader_schema(r), &sch, &err) == CARQUET_OK);
    assert(carquet_arrow_schema_copy(&sch, &cp) == CARQUET_OK);
    sch.release(&sch);
    assert(cp.n_children == 8 && strcmp(cp.children[1]->format, "d:9,2") == 0);
    assert(strcmp(cp.children[7]->name, "pt") == 0);
    cp.release(&cp);
    carquet_reader_close(r);
    TEST_PASS("arrow_ref");
    return 0;
}

int main(void) {
    char path[512], wide[512];
    carquet_test_temp_path(path, sizeof(path), "arrow_stream");
    carquet_test_temp_path(wide, sizeof(wide), "arrow_stream_wide");
    if (write_file(path, false) || write_file(wide, true)) {
        printf("[FAIL] arrow_stream: write\n");
        return 1;
    }
    int failures = 0;
    failures += test_stream_flat_and_types(path);
    failures += test_stream_projection(path);
    failures += test_stream_errors(path, wide);
    failures += test_stream_selection(path);
    failures += test_stream_float_predicates();
    failures += test_arrow_ref(path);
    carquet_test_cleanup(path);
    carquet_test_cleanup(wide);
    if (failures) return 1;
    printf("\nAll Arrow stream tests passed.\n");
    return 0;
}
