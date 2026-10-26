/**
 * @file fuzz_page_filter.c
 * @brief Fuzz target for the v0.6.0 page-level filter evaluator.
 *
 * The page filter (src/reader/page_filter.c) consumes attacker-controlled
 * column-index / offset-index statistics and compares them against caller
 * predicate values. Unlike fuzz_reader, this target installs a page filter
 * so that carquet_page_filter_eval_row_group() — and in particular the
 * fixed-width typed comparator over file-supplied page min/max — runs on
 * malformed input.
 *
 * The fuzzer derives well-typed clauses from the leaf schema so the
 * predicate side always validates; any out-of-bounds access therefore comes
 * from the file-supplied stats, i.e. a genuine library bug rather than a
 * malformed predicate from the harness.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <carquet/carquet.h>

#define MAX_LEAVES 64

/* Build leaf-index -> (physical type, type_length) in column order by
 * walking the schema elements and picking out leaves. This order matches
 * the leaf ordering the filter uses internally. */
static int32_t collect_leaf_types(const carquet_schema_t* schema,
                                  carquet_physical_type_t* ptypes,
                                  int32_t* type_lengths) {
    if (!schema) return 0;
    int32_t n_elem = carquet_schema_num_elements(schema);
    int32_t nleaf = 0;
    for (int32_t i = 0; i < n_elem && nleaf < MAX_LEAVES; i++) {
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, i);
        if (!node) continue;
        if (!carquet_schema_node_is_leaf(node)) continue;
        ptypes[nleaf] = carquet_schema_node_physical_type(node);
        type_lengths[nleaf] = carquet_schema_node_type_length(node);
        nleaf++;
    }
    return nleaf;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 16) return 0;
    (void)carquet_init();

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open_buffer(data, size, NULL, &err);
    if (!reader) return 0;

    int32_t num_cols = carquet_reader_num_columns(reader);
    if (num_cols <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    carquet_physical_type_t ptypes[MAX_LEAVES];
    int32_t type_lengths[MAX_LEAVES];
    int32_t nleaf = collect_leaf_types(schema, ptypes, type_lengths);
    if (nleaf <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    /* Use the first bytes of input to pick a column + op; the remaining
     * bytes provide a scratch predicate value buffer. */
    uint8_t sel_col = data[0];
    uint8_t sel_op = data[1];
    /* A 16-byte predicate value covers every numeric width and supplies
     * raw bytes for FLBA/BYTE_ARRAY. Comes from the fuzz input. */
    const uint8_t* val = &data[size >= 32 ? size - 16 : 0];

    int32_t col = (int32_t)(sel_col % (uint8_t)nleaf);
    carquet_physical_type_t pt = ptypes[col];
    int32_t tlen = type_lengths[col];

    /* INT96 is rejected by validate; skip to avoid a guaranteed error
     * (keeps the interesting paths reachable). */
    if (pt == CARQUET_PHYSICAL_INT96) {
        carquet_reader_close(reader);
        return 0;
    }

    carquet_filter_clause_t clause;
    memset(&clause, 0, sizeof(clause));
    clause.column_index = col;

    carquet_filter_op_t ops[] = {
        CARQUET_FILTER_EQ, CARQUET_FILTER_NE, CARQUET_FILTER_LT,
        CARQUET_FILTER_LE, CARQUET_FILTER_GT, CARQUET_FILTER_GE,
        CARQUET_FILTER_RANGE, CARQUET_FILTER_IN,
        CARQUET_FILTER_IS_NULL, CARQUET_FILTER_IS_NOT_NULL,
    };
    carquet_filter_op_t op = ops[sel_op % (sizeof(ops) / sizeof(ops[0]))];
    clause.op = op;

    /* value_size: for FLBA must equal type_length; for BYTE_ARRAY a byte
     * length; for numeric it is ignored. */
    int32_t vsize;
    switch (pt) {
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: vsize = tlen; break;
        case CARQUET_PHYSICAL_BYTE_ARRAY:           vsize = 4;    break;
        default:                                    vsize = 8;    break;
    }
    /* Guard: an FLBA with an implausibly large type_length would make our
     * 16-byte scratch value too small. Skip those. */
    if (pt == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && (tlen <= 0 || tlen > 16)) {
        carquet_reader_close(reader);
        return 0;
    }

    switch (op) {
        case CARQUET_FILTER_RANGE:
            clause.has_lo = true;
            clause.lo = val;
            clause.lo_size = vsize;
            clause.has_hi = true;
            clause.hi = val;
            clause.hi_size = vsize;
            break;
        case CARQUET_FILTER_IN:
            /* IN on BYTE_ARRAY expects carquet_byte_array_t entries; only
             * exercise IN for fixed-width numeric/FLBA where the packed
             * scalar layout matches what the evaluator reads. */
            if (pt == CARQUET_PHYSICAL_BYTE_ARRAY) {
                clause.op = CARQUET_FILTER_EQ;
                clause.value = val;
                clause.value_size = vsize;
            } else {
                clause.values = val;
                clause.value_count = 1;
            }
            break;
        case CARQUET_FILTER_IS_NULL:
        case CARQUET_FILTER_IS_NOT_NULL:
            /* No value needed. */
            break;
        default:
            clause.value = val;
            clause.value_size = vsize;
            break;
    }

    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = (int32_t)((data[2] % 64) + 1) * 32;
    config.num_threads = data[3] % 4;
    config.use_mmap = false;  /* buffer-backed reader */

    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &config, &err);
    if (br) {
        if (carquet_batch_reader_set_page_filter(br, &clause, 1) == CARQUET_OK) {
            carquet_row_batch_t* batch = NULL;
            int batch_count = 0;
            while (batch_count < 16 &&
                   carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
                for (int32_t c = 0; c < num_cols && c < 64; c++) {
                    const void* dp;
                    const uint8_t* nulls;
                    int64_t cnt;
                    (void)carquet_row_batch_column(batch, c, &dp, &nulls, &cnt);
                }
                carquet_row_batch_free(batch);
                batch = NULL;
                batch_count++;
            }
            (void)carquet_batch_reader_rows_skipped(br);
        }
        carquet_batch_reader_free(br);
    }

    carquet_reader_close(reader);
    return 0;
}

#ifdef AFL_MAIN
#include <stdio.h>
#include <sys/stat.h>
int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <input_file>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 1; }
    struct stat st; fstat(fileno(f), &st);
    uint8_t* d = malloc((size_t)st.st_size);
    if (!d) { fclose(f); return 1; }
    fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, (size_t)st.st_size);
    free(d); return r;
}
#endif
