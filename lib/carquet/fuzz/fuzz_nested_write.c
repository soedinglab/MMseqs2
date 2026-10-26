/**
 * @file fuzz_nested_write.c
 * @brief Fuzz target for the nested write helper (auto-shredding).
 *
 * carquet_writer_write_list_column() shreds caller-supplied Arrow offsets and
 * validity bitmaps into definition/repetition levels and allocates internal
 * buffers sized from offsets[num_lists]. This target drives it with arbitrary
 * offsets, validity, and child values over both LIST and MAP schemas, then
 * reads the result back so the reconstruction path is exercised too. The point
 * is that no input — malformed offsets, oversized child counts, mismatched
 * validity — may crash, leak, or read out of bounds; it may only return an
 * error status.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <carquet/carquet.h>

#define MAX_LISTS 256
#define MAX_CHILDREN 1024

typedef struct { const uint8_t* data; size_t size; size_t pos; } fuzz_input_t;

static uint8_t take8(fuzz_input_t* in) {
    return in->pos < in->size ? in->data[in->pos++] : 0;
}
static uint32_t take32(fuzz_input_t* in) {
    uint32_t a = take8(in), b = take8(in), c = take8(in), d = take8(in);
    return a | (b << 8) | (c << 16) | (d << 24);
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 4 || size > (1u << 20)) return 0;
    fuzz_input_t in = {data, size, 0};

    uint8_t sel = take8(&in);
    bool is_map = (sel & 1) != 0;
    bool list_optional = (sel & 2) != 0;
    bool provide_list_valid = (sel & 4) != 0;
    bool provide_value_valid = (sel & 8) != 0;

    int32_t num_lists = (int32_t)(take8(&in) | ((uint32_t)take8(&in) << 8));
    num_lists %= (MAX_LISTS + 1);

    /* Build monotonically *non-decreasing-ish* offsets, but allow the fuzzer
     * to inject non-monotonic / out-of-range steps that the API must reject. */
    int32_t* offsets = (int32_t*)malloc(sizeof(int32_t) * (size_t)(num_lists + 1));
    if (!offsets) return 0;
    offsets[0] = (int32_t)(take8(&in) & 0x3);  /* usually 0, sometimes not */
    for (int32_t i = 0; i < num_lists; i++) {
        int32_t step = (int32_t)(take8(&in) & 0x7);
        /* Occasionally emit a negative step to test rejection. */
        if ((take8(&in) & 0x1f) == 0) step = -step;
        int64_t next = (int64_t)offsets[i] + step;
        if (next > MAX_CHILDREN) next = MAX_CHILDREN;
        offsets[i + 1] = (int32_t)next;
    }
    int32_t total_children = offsets[num_lists] > 0 ? offsets[num_lists] : 0;
    if (total_children > MAX_CHILDREN) total_children = MAX_CHILDREN;

    /* Child values (int32 for LIST/map-key, int64 for map-value). */
    size_t vslot = is_map ? sizeof(int64_t) : sizeof(int32_t);
    void* values = calloc((size_t)total_children + 1, vslot);
    if (!values) { free(offsets); return 0; }
    for (int32_t i = 0; i < total_children; i++) {
        if (is_map) ((int64_t*)values)[i] = (int64_t)take32(&in);
        else        ((int32_t*)values)[i] = (int32_t)take32(&in);
    }

    size_t lv_bytes = (size_t)(num_lists / 8 + 1);
    size_t cv_bytes = (size_t)(total_children / 8 + 1);
    uint8_t* list_valid = (uint8_t*)malloc(lv_bytes);
    uint8_t* value_valid = (uint8_t*)malloc(cv_bytes);
    if (!list_valid || !value_valid) { free(offsets); free(values); free(list_valid); free(value_valid); return 0; }
    for (size_t i = 0; i < lv_bytes; i++) list_valid[i] = take8(&in);
    for (size_t i = 0; i < cv_bytes; i++) value_valid[i] = take8(&in);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    if (!schema) { free(offsets); free(values); free(list_valid); free(value_valid); return 0; }
    carquet_field_repetition_t rep =
        list_optional ? CARQUET_REPETITION_OPTIONAL : CARQUET_REPETITION_REQUIRED;
    int32_t leaf_col = 0;
    if (is_map) {
        carquet_schema_add_map(schema, "m",
            CARQUET_PHYSICAL_INT32, NULL, 0,
            CARQUET_PHYSICAL_INT64, NULL, 0, rep, 0);
        leaf_col = 1;  /* value leaf (OPTIONAL) */
    } else {
        carquet_schema_add_list(schema, "l", CARQUET_PHYSICAL_INT32, NULL, rep, 0, 0);
        leaf_col = 0;
    }

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (w) {
        carquet_status_t st = carquet_writer_write_list_column(
            w, leaf_col, num_lists, offsets,
            provide_list_valid ? list_valid : NULL,
            values,
            provide_value_valid ? value_valid : NULL, &err);
        if (st == CARQUET_OK && carquet_writer_close(w) == CARQUET_OK) {
            void* buf = NULL; size_t bsz = 0;
            /* get_buffer frees the writer on success; abort otherwise. */
            if (carquet_writer_get_buffer(w, &buf, &bsz) == CARQUET_OK && buf) {
                carquet_error_t rerr = CARQUET_ERROR_INIT;
                carquet_reader_t* r = carquet_reader_open_buffer(buf, bsz, NULL, &rerr);
                if (r) {
                    carquet_batch_reader_config_t cfg;
                    carquet_batch_reader_config_init(&cfg);
                    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &rerr);
                    if (br) {
                        carquet_row_batch_t* batch = NULL;
                        while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
                            const int32_t* off = NULL; int64_t nl = 0, nv = 0;
                            const void* vv = NULL; const uint8_t* cvld = NULL; const uint8_t* lvld = NULL;
                            (void)carquet_row_batch_column_list(batch, leaf_col, &off, &nl,
                                                                &vv, &cvld, &nv, &lvld);
                            carquet_row_batch_free(batch);
                            batch = NULL;
                        }
                        carquet_batch_reader_free(br);
                    }
                    carquet_reader_close(r);
                }
                free(buf);
            } else {
                carquet_writer_abort(w);
            }
        } else {
            carquet_writer_abort(w);
        }
    }

    carquet_schema_free(schema);
    free(offsets);
    free(values);
    free(list_valid);
    free(value_valid);
    return 0;
}
