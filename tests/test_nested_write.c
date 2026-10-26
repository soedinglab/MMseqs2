/**
 * @file test_nested_write.c
 * @brief Tests for the nested write helper (carquet_writer_write_list_column).
 *
 * Auto-shredding takes Arrow-style offsets + validity bitmaps and produces the
 * definition/repetition levels the low-level writer needs. These tests:
 *   - round-trip a LIST<int32> through the batch reader,
 *   - assert the shredded output is BYTE-IDENTICAL to a hand-computed
 *     carquet_writer_write_batch() with the same levels,
 *   - round-trip a MAP<int32,int64> (two leaves sharing the map offsets),
 *   - exercise the error/rejection and boundary paths.
 *
 * Dataset (5 lists), mirroring test_batch_nested.c:
 *   row0: [10, 20]      row1: []       row2: NULL
 *   row3: [NULL]        row4: [30, 40]
 * Offsets = {0,2,2,2,3,5}; child array (5) = {10,20,<null>,30,40};
 * list validity = 0b11011 (row2 null); child validity = 0b11011 (child2 null).
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

/* Shared dataset. */
static const int32_t DS_OFFSETS[6] = {0, 2, 2, 2, 3, 5};
static const int32_t DS_CHILD[5]   = {10, 20, 0 /*null*/, 30, 40};
static const uint8_t DS_LIST_VALID  = 0x1B; /* bits 0,1,3,4 set; bit2 clear */
static const uint8_t DS_CHILD_VALID = 0x1B; /* child 2 is null              */

/* ------------------------------------------------------------------------ */

static carquet_schema_t* make_list_schema(void) {
    carquet_schema_t* schema = carquet_schema_create(NULL);
    if (!schema) return NULL;
    carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_INT32, NULL,
                            CARQUET_REPETITION_OPTIONAL, 0, 0);
    return schema;
}

/* Write the dataset via the auto-shredding helper. */
static int write_via_helper(void** buf, size_t* sz) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_list_schema();
    if (!schema) return -1;
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) { carquet_schema_free(schema); return -1; }
    carquet_status_t st = carquet_writer_write_list_column(
        w, 0, 5, DS_OFFSETS, &DS_LIST_VALID, DS_CHILD, &DS_CHILD_VALID, &err);
    if (st != CARQUET_OK) { carquet_writer_abort(w); carquet_schema_free(schema); return -1; }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(schema); return -1; }
    st = carquet_writer_get_buffer(w, buf, sz);
    carquet_schema_free(schema);
    return st == CARQUET_OK ? 0 : -1;
}

/* Write the dataset via hand-computed levels (the oracle). */
static int write_via_manual(void** buf, size_t* sz) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_list_schema();
    if (!schema) return -1;
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) { carquet_schema_free(schema); return -1; }
    int32_t vals[4] = {10, 20, 30, 40};
    int16_t def[7] = {3, 3, 1, 0, 2, 3, 3};
    int16_t rep[7] = {0, 1, 0, 0, 0, 0, 1};
    if (carquet_writer_write_batch(w, 0, vals, 7, def, rep) != CARQUET_OK) {
        carquet_writer_abort(w); carquet_schema_free(schema); return -1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(schema); return -1; }
    carquet_status_t st = carquet_writer_get_buffer(w, buf, sz);
    carquet_schema_free(schema);
    return st == CARQUET_OK ? 0 : -1;
}

/* Round-trip the list and validate against the dataset. */
static int test_list_roundtrip(void) {
    g_current_test = "list_roundtrip";
    void* buf = NULL; size_t sz = 0;
    ASSERT_TRUE(write_via_helper(&buf, &sz) == 0);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    ASSERT_TRUE(r != NULL);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);

    carquet_row_batch_t* batch = NULL;
    ASSERT_OK(carquet_batch_reader_next(br, &batch));
    ASSERT_TRUE(batch != NULL);

    const int32_t* offsets = NULL;
    int64_t num_lists = 0, num_values = 0;
    const void* values = NULL;
    const uint8_t* vvalid = NULL;
    const uint8_t* lvalid = NULL;
    ASSERT_OK(carquet_row_batch_column_list(batch, 0, &offsets, &num_lists,
                                            &values, &vvalid, &num_values, &lvalid));
    ASSERT_TRUE(num_lists == 5);
    ASSERT_TRUE(num_values == 5);
    for (int i = 0; i < 6; i++) ASSERT_TRUE(offsets[i] == DS_OFFSETS[i]);
    const int32_t* v = (const int32_t*)values;
    ASSERT_TRUE(v[0] == 10 && v[1] == 20 && v[3] == 30 && v[4] == 40);
    ASSERT_TRUE(vvalid != NULL && lvalid != NULL);
    ASSERT_TRUE(!(vvalid[0] & (1 << 2)));   /* child 2 null */
    ASSERT_TRUE(!(lvalid[0] & (1 << 2)));   /* list 2 null  */
    ASSERT_TRUE((lvalid[0] & (1 << 0)) && (lvalid[0] & (1 << 4)));

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS(g_current_test);
    return 0;
}

/* The shredded bytes must exactly equal the hand-computed-levels write. */
static int test_byte_identical_to_manual(void) {
    g_current_test = "byte_identical_to_manual";
    void* a = NULL; size_t asz = 0;
    void* b = NULL; size_t bsz = 0;
    ASSERT_TRUE(write_via_helper(&a, &asz) == 0);
    ASSERT_TRUE(write_via_manual(&b, &bsz) == 0);
    ASSERT_TRUE(asz == bsz);
    ASSERT_TRUE(memcmp(a, b, asz) == 0);
    free(a);
    free(b);
    TEST_PASS(g_current_test);
    return 0;
}

/* MAP<int32,int64>: key (REQUIRED) and value (OPTIONAL) share the entry
 * offsets and the map-level validity; only value carries element validity. */
static int test_map_roundtrip(void) {
    g_current_test = "map_roundtrip";
    /* 3 maps:
     *   row0: {1:100, 2:NULL}   row1: {}   row2: {3:300}
     * entry offsets = {0,2,2,3}; keys = {1,2,3}; values = {100,<null>,300};
     * value validity = 0b101 (entry1 null); map validity = all present. */
    const int32_t offsets[4] = {0, 2, 2, 3};
    const int32_t keys[3]    = {1, 2, 3};
    const int64_t vals[3]    = {100, 0 /*null*/, 300};
    const uint8_t vvalid     = 0x05;  /* entries 0,2 present; entry1 null */

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    ASSERT_TRUE(schema != NULL);
    ASSERT_TRUE(carquet_schema_add_map(schema, "m",
        CARQUET_PHYSICAL_INT32, NULL, 0,
        CARQUET_PHYSICAL_INT64, NULL, 0,
        CARQUET_REPETITION_OPTIONAL, 0) >= 0);

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    ASSERT_TRUE(w != NULL);
    /* key leaf = column 0, value leaf = column 1. */
    ASSERT_OK(carquet_writer_write_list_column(w, 0, 3, offsets, NULL,
                                               keys, NULL, &err));
    ASSERT_OK(carquet_writer_write_list_column(w, 1, 3, offsets, NULL,
                                               vals, &vvalid, &err));
    ASSERT_OK(carquet_writer_close(w));
    void* buf = NULL; size_t sz = 0;
    ASSERT_OK(carquet_writer_get_buffer(w, &buf, &sz));
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    ASSERT_TRUE(r != NULL);
    ASSERT_TRUE(carquet_reader_num_rows(r) == 3);
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    ASSERT_TRUE(br != NULL);
    carquet_row_batch_t* batch = NULL;
    ASSERT_OK(carquet_batch_reader_next(br, &batch));
    ASSERT_TRUE(batch != NULL);

    /* Key column. */
    const int32_t* koff = NULL; int64_t knl = 0, knv = 0;
    const void* kv = NULL; const uint8_t* kvv = NULL; const uint8_t* klv = NULL;
    ASSERT_OK(carquet_row_batch_column_list(batch, 0, &koff, &knl, &kv, &kvv, &knv, &klv));
    ASSERT_TRUE(knl == 3);
    ASSERT_TRUE(koff[0] == 0 && koff[1] == 2 && koff[2] == 2 && koff[3] == 3);
    const int32_t* kd = (const int32_t*)kv;
    ASSERT_TRUE(kd[0] == 1 && kd[1] == 2 && kd[2] == 3);

    /* Value column. */
    const int32_t* voff = NULL; int64_t vnl = 0, vnv = 0;
    const void* vv = NULL; const uint8_t* vvv = NULL; const uint8_t* vlv = NULL;
    ASSERT_OK(carquet_row_batch_column_list(batch, 1, &voff, &vnl, &vv, &vvv, &vnv, &vlv));
    ASSERT_TRUE(vnl == 3);
    ASSERT_TRUE(voff[0] == 0 && voff[1] == 2 && voff[2] == 2 && voff[3] == 3);
    const int64_t* vd = (const int64_t*)vv;
    ASSERT_TRUE(vd[0] == 100 && vd[2] == 300);
    ASSERT_TRUE(vvv != NULL && !(vvv[0] & (1 << 1)));  /* entry 1 value null */

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS(g_current_test);
    return 0;
}

/* Error and boundary paths. */
static int test_error_paths(void) {
    g_current_test = "error_paths";
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Flat (non-repeated) column ⇒ NOT_IMPLEMENTED. */
    {
        carquet_schema_t* schema = carquet_schema_create(NULL);
        ASSERT_TRUE(schema != NULL);
        ASSERT_OK(carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
            NULL, CARQUET_REPETITION_OPTIONAL, 0, 0));
        carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
        ASSERT_TRUE(w != NULL);
        int32_t off[2] = {0, 1};
        int32_t val[1] = {7};
        ASSERT_TRUE(carquet_writer_write_list_column(w, 0, 1, off, NULL, val, NULL, &err)
                    == CARQUET_ERROR_NOT_IMPLEMENTED);
        carquet_writer_abort(w);
        carquet_schema_free(schema);
    }

    /* Column index out of range. */
    {
        carquet_schema_t* schema = make_list_schema();
        ASSERT_TRUE(schema != NULL);
        carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
        ASSERT_TRUE(w != NULL);
        int32_t off[2] = {0, 1};
        int32_t val[1] = {7};
        ASSERT_TRUE(carquet_writer_write_list_column(w, 9, 1, off, NULL, val, NULL, &err)
                    == CARQUET_ERROR_INVALID_ARGUMENT);
        ASSERT_TRUE(carquet_writer_write_list_column(w, -1, 1, off, NULL, val, NULL, &err)
                    == CARQUET_ERROR_INVALID_ARGUMENT);
        carquet_writer_abort(w);
        carquet_schema_free(schema);
    }

    /* offsets[0] != 0 and non-monotonic offsets are rejected. */
    {
        carquet_schema_t* schema = make_list_schema();
        ASSERT_TRUE(schema != NULL);
        carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
        ASSERT_TRUE(w != NULL);
        int32_t bad0[3]  = {1, 2, 3};      /* offsets[0] != 0 */
        int32_t val[3]   = {1, 2, 3};
        ASSERT_TRUE(carquet_writer_write_list_column(w, 0, 2, bad0, NULL, val, NULL, &err)
                    == CARQUET_ERROR_INVALID_ARGUMENT);
        int32_t nonmono[3] = {0, 2, 1};    /* decreasing */
        ASSERT_TRUE(carquet_writer_write_list_column(w, 0, 2, nonmono, NULL, val, NULL, &err)
                    == CARQUET_ERROR_INVALID_ARGUMENT);
        carquet_writer_abort(w);
        carquet_schema_free(schema);
    }

    /* Null element written into a REQUIRED leaf is rejected. add_map keys are
     * REQUIRED: clearing a value-validity bit for the key leaf must fail. */
    {
        carquet_schema_t* schema = carquet_schema_create(NULL);
        ASSERT_TRUE(schema != NULL);
        ASSERT_TRUE(carquet_schema_add_map(schema, "m",
            CARQUET_PHYSICAL_INT32, NULL, 0,
            CARQUET_PHYSICAL_INT64, NULL, 0,
            CARQUET_REPETITION_OPTIONAL, 0) >= 0);
        carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
        ASSERT_TRUE(w != NULL);
        int32_t off[2] = {0, 1};
        int32_t keys[1] = {5};
        uint8_t none = 0x00;  /* key element 0 marked null */
        ASSERT_TRUE(carquet_writer_write_list_column(w, 0, 1, off, NULL, keys, &none, &err)
                    == CARQUET_ERROR_INVALID_ARGUMENT);
        carquet_writer_abort(w);
        carquet_schema_free(schema);
    }

    /* Boundary: num_lists == 0 is a no-op that still produces a valid file. */
    {
        carquet_schema_t* schema = make_list_schema();
        ASSERT_TRUE(schema != NULL);
        carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
        ASSERT_TRUE(w != NULL);
        ASSERT_OK(carquet_writer_write_list_column(w, 0, 0, NULL, NULL, NULL, NULL, &err));
        ASSERT_OK(carquet_writer_close(w));
        void* buf = NULL; size_t sz = 0;
        ASSERT_OK(carquet_writer_get_buffer(w, &buf, &sz));
        carquet_schema_free(schema);
        carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
        ASSERT_TRUE(r != NULL);
        ASSERT_TRUE(carquet_reader_num_rows(r) == 0);
        carquet_reader_close(r);
        free(buf);
    }

    TEST_PASS(g_current_test);
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_list_roundtrip();
    failures += test_byte_identical_to_manual();
    failures += test_map_roundtrip();
    failures += test_error_paths();
    if (failures > 0) {
        fprintf(stderr, "%d test failures\n", failures);
        return 1;
    }
    return 0;
}
