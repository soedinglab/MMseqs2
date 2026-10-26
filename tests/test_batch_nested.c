/**
 * @file test_batch_nested.c
 * @brief Tests for LIST reconstruction in the batch reader.
 *
 * Verifies carquet_row_batch_column_list(): Arrow List<T> offsets, flattened
 * child values, child validity (null elements), and list-level validity (null
 * lists). Covers a flat + list mix, a list-only projection (nested column
 * driving row-group advance), and the rejection of a list column through the
 * flat accessor.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Dataset (5 rows):
 *   id   INT32 REQUIRED         = 1,2,3,4,5
 *   tags LIST<int32 optional>:
 *     row0: [10, 20]
 *     row1: []            (empty list)
 *     row2: NULL          (null list)
 *     row3: [NULL]        (one null element)
 *     row4: [30, 40]
 *
 *   Leaf slots (def/rep), 7 entries; packed non-null values = {10,20,30,40}:
 *     def = {3,3, 1, 0, 2, 3,3}
 *     rep = {0,1, 0, 0, 0, 0,1}
 */
static int write_dataset(void** buf, size_t* sz, bool with_id) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    int32_t elem_col;
    if (with_id) {
        carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_INT32, NULL,
                                CARQUET_REPETITION_OPTIONAL, 0, 0);
        elem_col = 1;
    } else {
        carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_INT32, NULL,
                                CARQUET_REPETITION_OPTIONAL, 0, 0);
        elem_col = 0;
    }

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) { carquet_schema_free(schema); return -1; }

    if (with_id) {
        int32_t ids[5] = {1, 2, 3, 4, 5};
        if (carquet_writer_write_batch(w, 0, ids, 5, NULL, NULL) != CARQUET_OK) return -1;
    }
    int32_t tag_vals[4] = {10, 20, 30, 40};
    int16_t tag_def[7] = {3, 3, 1, 0, 2, 3, 3};
    int16_t tag_rep[7] = {0, 1, 0, 0, 0, 0, 1};
    if (carquet_writer_write_batch(w, elem_col, tag_vals, 7, tag_def, tag_rep) != CARQUET_OK)
        return -1;

    if (carquet_writer_close(w) != CARQUET_OK) return -1;
    if (carquet_writer_get_buffer(w, buf, sz) != CARQUET_OK) return -1;
    carquet_schema_free(schema);
    return 0;
}

/* Assert the reconstructed "tags" list column matches the dataset. */
static void check_tags(carquet_row_batch_t* batch, int32_t col) {
    const int32_t* offsets = NULL;
    int64_t num_lists = 0, num_values = 0;
    const void* values = NULL;
    const uint8_t* vvalid = NULL;
    const uint8_t* lvalid = NULL;
    carquet_status_t st = carquet_row_batch_column_list(
        batch, col, &offsets, &num_lists, &values, &vvalid, &num_values, &lvalid);
    assert(st == CARQUET_OK);

    assert(num_lists == 5);
    assert(num_values == 5);
    int32_t exp_off[6] = {0, 2, 2, 2, 3, 5};
    for (int i = 0; i < 6; i++) assert(offsets[i] == exp_off[i]);

    const int32_t* v = (const int32_t*)values;
    assert(v[0] == 10 && v[1] == 20 && v[3] == 30 && v[4] == 40);
    /* child[2] is a null element (value zeroed). */

    /* Child validity: elements 0,1,3,4 present; element 2 null. */
    assert(vvalid != NULL);
    assert(vvalid[0] & (1 << 0));
    assert(vvalid[0] & (1 << 1));
    assert(!(vvalid[0] & (1 << 2)));
    assert(vvalid[0] & (1 << 3));
    assert(vvalid[0] & (1 << 4));

    /* List validity: rows 0,1,3,4 present; row 2 is a null list. */
    assert(lvalid != NULL);
    assert(lvalid[0] & (1 << 0));
    assert(lvalid[0] & (1 << 1));
    assert(!(lvalid[0] & (1 << 2)));
    assert(lvalid[0] & (1 << 3));
    assert(lvalid[0] & (1 << 4));

    /* Row semantics: row0 [10,20]; row1 empty; row2 null; row3 [null]; row4 [30,40]. */
    assert(offsets[1] - offsets[0] == 2);   /* row0: 2 elements */
    assert(offsets[2] - offsets[1] == 0);   /* row1: empty */
    assert(offsets[3] - offsets[2] == 0);   /* row2: null list, 0 elements */
    assert(offsets[4] - offsets[3] == 1);   /* row3: 1 (null) element */
    assert(offsets[5] - offsets[4] == 2);   /* row4: 2 elements */

    /* The flat accessor must reject a list column. */
    const void* d; const uint8_t* nb; int64_t nv;
    assert(carquet_row_batch_column(batch, col, &d, &nb, &nv) == CARQUET_ERROR_INVALID_ARGUMENT);
}

static int test_flat_plus_list(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz, true) != 0) TEST_FAIL("flat_plus_list", "write failed");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) { free(buf); TEST_FAIL("flat_plus_list", "open failed"); }
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 100;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); free(buf); TEST_FAIL("flat_plus_list", "br create"); }

    carquet_row_batch_t* batch = NULL;
    if (carquet_batch_reader_next(br, &batch) != CARQUET_OK || !batch)
        { carquet_batch_reader_free(br); carquet_reader_close(r); free(buf); TEST_FAIL("flat_plus_list", "next"); }

    assert(carquet_row_batch_num_rows(batch) == 5);
    assert(carquet_row_batch_num_columns(batch) == 2);

    /* Flat id column. */
    const void* id_data; const uint8_t* id_nb; int64_t id_nv;
    assert(carquet_row_batch_column(batch, 0, &id_data, &id_nb, &id_nv) == CARQUET_OK);
    assert(id_nv == 5);
    const int32_t* ids = (const int32_t*)id_data;
    for (int i = 0; i < 5; i++) assert(ids[i] == i + 1);

    /* Reconstructed tags list column. */
    check_tags(batch, 1);

    /* One row group -> next returns END_OF_DATA. */
    carquet_row_batch_t* b2 = NULL;
    assert(carquet_batch_reader_next(br, &b2) == CARQUET_ERROR_END_OF_DATA);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS("flat_plus_list");
    return 0;
}

static int test_list_only(void) {
    void* buf = NULL; size_t sz = 0;
    if (write_dataset(&buf, &sz, false) != 0) TEST_FAIL("list_only", "write failed");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) { free(buf); TEST_FAIL("list_only", "open failed"); }
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 100;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); free(buf); TEST_FAIL("list_only", "br create"); }

    carquet_row_batch_t* batch = NULL;
    if (carquet_batch_reader_next(br, &batch) != CARQUET_OK || !batch)
        { carquet_batch_reader_free(br); carquet_reader_close(r); free(buf); TEST_FAIL("list_only", "next"); }

    assert(carquet_row_batch_num_rows(batch) == 5);
    assert(carquet_row_batch_num_columns(batch) == 1);
    check_tags(batch, 0);

    carquet_row_batch_t* b2 = NULL;
    assert(carquet_batch_reader_next(br, &b2) == CARQUET_ERROR_END_OF_DATA);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    free(buf);
    TEST_PASS("list_only");
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_flat_plus_list();
    rc |= test_list_only();
    if (rc == 0) printf("\nAll batch-reader nested tests passed.\n");
    return rc;
}
