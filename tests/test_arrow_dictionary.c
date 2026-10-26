/**
 * @file test_arrow_dictionary.c
 * @brief Regression test for writing dictionary-encoded Arrow arrays.
 *
 * A dictionary-encoded Arrow array (e.g. a pandas categorical / pyarrow
 * DictionaryArray) carries integer indices in its primary buffers and the
 * actual values in `array->dictionary` (typed by `schema->dictionary`). The
 * write bridge previously ignored the dictionary and wrote the raw index buffer
 * verbatim as an integer column — wrong data, silently. It now resolves each
 * index to its dictionary value and writes the value column. This test writes a
 * categorical string array and asserts the materialized string values on read.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- minimal heap-owned Arrow C builders (freed via release) ------------- */
static void rel_schema(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format); free((void*)s->name);
    for (int64_t i = 0; i < s->n_children; i++)
        if (s->children[i]) { if (s->children[i]->release) s->children[i]->release(s->children[i]); free(s->children[i]); }
    if (s->dictionary) { if (s->dictionary->release) s->dictionary->release(s->dictionary); free(s->dictionary); }
    free(s->children); s->release = NULL;
}
static void rel_array(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) { for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]); free(a->buffers); }
    for (int64_t i = 0; i < a->n_children; i++)
        if (a->children[i]) { if (a->children[i]->release) a->children[i]->release(a->children[i]); free(a->children[i]); }
    if (a->dictionary) { if (a->dictionary->release) a->dictionary->release(a->dictionary); free(a->dictionary); }
    free(a->children); a->release = NULL;
}
static char* dups(const char* s) { char* o = malloc(strlen(s) + 1); strcpy(o, s); return o; }
static struct ArrowSchema* S(const char* fmt, const char* name, int nullable, int nch) {
    struct ArrowSchema* s = calloc(1, sizeof(*s));
    s->format = dups(fmt); s->name = dups(name); s->flags = nullable ? ARROW_FLAG_NULLABLE : 0;
    s->n_children = nch; s->children = nch ? calloc(nch, sizeof(void*)) : NULL; s->release = rel_schema;
    return s;
}
static struct ArrowArray* A(int64_t len, int nbuf, int nch) {
    struct ArrowArray* a = calloc(1, sizeof(*a));
    a->length = len; a->null_count = 0; a->n_buffers = nbuf; a->n_children = nch;
    a->buffers = nbuf ? calloc(nbuf, sizeof(void*)) : NULL;
    a->children = nch ? calloc(nch, sizeof(void*)) : NULL;
    a->release = rel_array;
    return a;
}
static void* dupb(const void* p, size_t n) { void* o = malloc(n ? n : 1); if (n) memcpy(o, p, n); return o; }

/* Build the dictionary VALUES array: utf8 ["red","green","blue"]. */
static struct ArrowArray* make_dict_values(void) {
    static const char data[] = "redgreenblue";
    static const int32_t off[4] = {0, 3, 8, 12};
    struct ArrowArray* v = A(3, 3, 0);
    v->buffers[0] = NULL;                      /* validity: all valid */
    v->buffers[1] = dupb(off, sizeof(off));    /* int32 offsets (n+1) */
    v->buffers[2] = dupb(data, 12);            /* char data */
    return v;
}

static int test_categorical_string(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Schema: struct { color: dictionary<int8, utf8> (non-null) }. */
    struct ArrowSchema* sc = S("+s", "schema", 0, 1);
    sc->children[0] = S("c", "color", 0, 0);              /* index type int8 */
    sc->children[0]->dictionary = S("u", "item", 0, 0);   /* value type utf8 */

    carquet_schema_t* cs = NULL;
    if (carquet_arrow_import_schema(sc, &cs, &err) != CARQUET_OK)
        TEST_FAIL("categorical_string", err.message);
    free(sc);  /* import released internals; free the top node */

    /* Data: indices [0,2,1,0,2] -> ["red","blue","green","red","blue"]. */
    static const int8_t idx[5] = {0, 2, 1, 0, 2};
    struct ArrowSchema* sc2 = S("+s", "schema", 0, 1);
    sc2->children[0] = S("c", "color", 0, 0);
    sc2->children[0]->dictionary = S("u", "item", 0, 0);

    struct ArrowArray* icol = A(5, 2, 0);
    icol->buffers[0] = NULL;                          /* validity: all valid */
    icol->buffers[1] = dupb(idx, sizeof(idx));        /* int8 indices */
    icol->dictionary = make_dict_values();
    struct ArrowArray* top = A(5, 1, 1);
    top->buffers[0] = NULL;
    top->children[0] = icol;

    char path[512]; carquet_test_temp_path(path, sizeof(path), "arrow_dict");
    carquet_writer_t* w = carquet_writer_create(path, cs, NULL, &err);
    if (!w) { carquet_schema_free(cs); TEST_FAIL("categorical_string", "create"); }
    if (carquet_writer_write_arrow(w, top, sc2, &err) != CARQUET_OK) {
        carquet_schema_free(cs); TEST_FAIL("categorical_string", err.message);
    }
    free(top); free(sc2);  /* write_arrow released internals */
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(cs); TEST_FAIL("categorical_string", "close"); }
    carquet_schema_free(cs);

    /* Read the column back and verify the resolved string values. */
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("categorical_string", "open"); }
    carquet_batch_reader_config_t cfg; carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    carquet_row_batch_t* batch = NULL;
    if (!br || carquet_batch_reader_next(br, &batch) != CARQUET_OK || !batch) {
        if (br) carquet_batch_reader_free(br); carquet_reader_close(r); carquet_test_cleanup(path);
        TEST_FAIL("categorical_string", "read");
    }
    const void* data = NULL; const uint8_t* nulls = NULL; int64_t n = 0;
    if (carquet_row_batch_column(batch, 0, &data, &nulls, &n) != CARQUET_OK || n != 5) {
        carquet_batch_reader_free(br); carquet_reader_close(r); carquet_test_cleanup(path);
        TEST_FAIL("categorical_string", "column");
    }
    const carquet_byte_array_t* ba = (const carquet_byte_array_t*)data;
    const char* expect[5] = {"red", "blue", "green", "red", "blue"};
    int ok = 1;
    for (int i = 0; i < 5; i++) {
        size_t elen = strlen(expect[i]);
        if ((size_t)ba[i].length != elen || memcmp(ba[i].data, expect[i], elen) != 0) { ok = 0; break; }
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("categorical_string", "dictionary values not materialized");
    TEST_PASS("categorical_string");
    return 0;
}

/* Nested (non-primitive) dictionary values are explicitly rejected, not
 * silently corrupted. */
static int test_nested_dict_rejected(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    struct ArrowSchema* sc = S("+s", "schema", 0, 1);
    sc->children[0] = S("c", "x", 0, 0);
    sc->children[0]->dictionary = S("+s", "item", 0, 0);  /* struct value: unsupported */
    carquet_schema_t* cs = NULL;
    carquet_status_t st = carquet_arrow_import_schema(sc, &cs, &err);
    free(sc);
    if (cs) carquet_schema_free(cs);
    if (st != CARQUET_ERROR_NOT_IMPLEMENTED) TEST_FAIL("nested_dict_rejected", "expected NOT_IMPLEMENTED");
    TEST_PASS("nested_dict_rejected");
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_categorical_string();
    rc |= test_nested_dict_rejected();
    if (rc == 0) printf("\nAll arrow-dictionary tests passed.\n");
    return rc;
}
