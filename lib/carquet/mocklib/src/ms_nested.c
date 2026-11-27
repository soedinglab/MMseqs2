/**
 * @file ms_nested.c
 * @brief Nested-data self-check: a LIST<STRING> column round-tripped and
 *        verified, plus the repetition-level helper functions.
 *
 * Exercises carquet_schema_add_list, carquet_writer_write_list_column,
 * carquet_row_batch_column_list, carquet_list_offsets and carquet_count_rows.
 */
#include "ms_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ms_nested_selfcheck(const char* path, char* errbuf, size_t errlen) {
    if (!path) { ms_set_err(errbuf, errlen, "path is NULL"); return -1; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_logical_type_t str = { .id = CARQUET_LOGICAL_STRING };

    /* Schema: id INT64 (leaf 0) + tags LIST<STRING> (element leaf 1). */
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) { ms_set_err(errbuf, errlen, err.message); return -1; }
    carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                            CARQUET_REPETITION_OPTIONAL, 0, 0);
    const int32_t TAGS_COL = 1;

    /* Four rows. tags: ["a","bb"], [], ["ccc"], ["d","e","f"]. */
    int64_t ids[4] = { 10, 20, 30, 40 };
    const char* strs[6] = { "a", "bb", "ccc", "d", "e", "f" };
    carquet_byte_array_t vals[6];
    for (int i = 0; i < 6; i++) {
        vals[i].data = (uint8_t*)strs[i];
        vals[i].length = (int32_t)strlen(strs[i]);
    }
    int32_t offsets[5] = { 0, 2, 2, 3, 6 };

    carquet_writer_t* w = carquet_writer_create(path, schema, NULL, &err);
    if (!w) { ms_set_err(errbuf, errlen, err.message); carquet_schema_free(schema); return -1; }

    int rc = -1;
    if (carquet_writer_write_batch(w, 0, ids, 4, NULL, NULL) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, "write id failed"); goto wfail;
    }
    if (carquet_writer_write_list_column(w, TAGS_COL, 4, offsets,
                                         /*list_validity=*/NULL, vals,
                                         /*value_validity=*/NULL, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message); goto wfail;
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, "close failed"); carquet_schema_free(schema); return -1;
    }
    carquet_schema_free(schema);
    goto readback;

wfail:
    carquet_writer_abort(w);
    carquet_schema_free(schema);
    return -1;

readback:;
    /* Read the list column back and verify boundaries + values. */
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { ms_set_err(errbuf, errlen, err.message); return -1; }

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 16;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { ms_set_err(errbuf, errlen, err.message); carquet_reader_close(r); return -1; }

    carquet_row_batch_t* batch = NULL;
    if (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const int32_t* off = NULL;
        int64_t num_lists = 0, num_values = 0;
        const void* values = NULL;
        const uint8_t* vvalid = NULL;
        const uint8_t* lvalid = NULL;
        if (carquet_row_batch_column_list(batch, TAGS_COL, &off, &num_lists,
                                          &values, &vvalid, &num_values,
                                          &lvalid) == CARQUET_OK) {
            const carquet_byte_array_t* v = (const carquet_byte_array_t*)values;
            bool ok = (num_lists == 4 && num_values == 6 &&
                       off[0] == 0 && off[1] == 2 && off[2] == 2 &&
                       off[3] == 3 && off[4] == 6);
            for (int i = 0; ok && i < 6; i++)
                ok = ((int)v[i].length == (int)strlen(strs[i]) &&
                      memcmp(v[i].data, strs[i], v[i].length) == 0);
            rc = ok ? 0 : -1;
            if (!ok) ms_set_err(errbuf, errlen, "list reconstruction mismatch");
        } else {
            ms_set_err(errbuf, errlen, "row_batch_column_list failed");
        }
        carquet_row_batch_free(batch);
    } else {
        ms_set_err(errbuf, errlen, "no batch returned");
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    if (rc != 0) return rc;

    /* Exercise the repetition-level helpers on a documented synthetic array:
       rep==0 marks a new logical row. 5 zeros => 5 rows; offsets bound lists. */
    int16_t rep[7]  = { 0, 1, 0, 0, 0, 0, 1 };
    int64_t reconstructed[8];
    int64_t rows = carquet_count_rows(rep, 7);
    int64_t nl = carquet_list_offsets(rep, 7, /*list_rep_level=*/1, reconstructed, 8);
    if (rows != 5 || nl != 5) {
        ms_set_err(errbuf, errlen, "rep-level helpers gave unexpected counts");
        return -1;
    }
    return 0;
}
