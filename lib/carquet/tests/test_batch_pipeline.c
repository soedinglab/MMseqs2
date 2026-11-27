/**
 * @file test_batch_pipeline.c
 * @brief Batch-reader pipeline coverage beyond REQUIRED + mmap: flat OPTIONAL
 *        columns (row-aligned values + validity bitmap), dictionary chunks
 *        split into segments decoded through their own column readers, and
 *        files read without mmap (positional I/O).
 *
 * Every value is a function of its global row index, so each batch is checked
 * row by row against the generator, nulls included, for several batch sizes
 * (multiples of 8 serve bitmap views, others shifted copies), with and
 * without mmap, and once with a page filter installed after creation (which
 * moves a pipeline with OPTIONAL columns to the sequential filtered path),
 * and with preserve_dictionaries (dictionary indices must reach the caller).
 * A second file covers BYTE_ARRAY columns (OPTIONAL dictionary and PLAIN,
 * REQUIRED), whose values point into pages owned by the pipeline's column
 * readers: batches are checked before the next call, exported Arrow arrays
 * after the readers are gone, and a page filter that seeks between ranges.
 * Last, OPTIONAL columns of every value width through null shapes that
 * stress the definition-level bitmap decode (see "Definition levels decoded
 * as a validity bitmap").
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include <carquet/carquet.h>
#include "test_helpers.h"
#include "reader/reader_internal.h"
#include "core/allocator.h"

#define NUM_RG      3
#define RG_ROWS     200000   /* >= 2 x 65536: chunks are split into segments */
#define TOTAL_ROWS  ((int64_t)NUM_RG * RG_ROWS)
#define NUM_COLS    5
#define PAGE_ROWS   10007

/* Column generators. a: INT64 OPTIONAL PLAIN; b: DOUBLE OPTIONAL dictionary;
 * c: INT32 REQUIRED dictionary; d: INT32 REQUIRED PLAIN; e: BOOLEAN OPTIONAL. */
static bool a_null(int64_t r) { return r % 7 == 3 || (r / 1000) % 11 == 5; }
static int64_t a_val(int64_t r) { return r * 3 - 7; }
static bool b_null(int64_t r) { return r % 5 == 0; }
static double b_val(int64_t r) { return (double)(r % 37) * 0.5; }
static int32_t c_val(int64_t r) { return (int32_t)(r % 113); }
static int32_t d_val(int64_t r) { return (int32_t)(r ^ 0x55); }
static bool e_null(int64_t r) { return r % 13 == 1; }
static uint8_t e_val(int64_t r) { return r % 3 == 0; }

static int write_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;
    if (carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_DOUBLE, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "c", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "d", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "e", CARQUET_PHYSICAL_BOOLEAN, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return 1;
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    /* Many pages per chunk (dictionary ones included), holding row counts
     * that are not multiples of 8, so segments share bitmap bytes with their
     * neighbours. The row cap applies between write calls, so each row group
     * is written in PAGE_ROWS chunks. */
    opts.max_rows_per_page = PAGE_ROWS;
    opts.write_page_index = true;        /* for the page-filter case */
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) {
        carquet_schema_free(schema);
        return 1;
    }
    int rc = 0;
    if (carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK ||
        carquet_writer_set_column_encoding(w, 2, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK) {
        rc = 1;
    }

    int64_t* a = malloc(RG_ROWS * sizeof(int64_t));
    double* b = malloc(RG_ROWS * sizeof(double));
    int32_t* c = malloc(RG_ROWS * sizeof(int32_t));
    int32_t* d = malloc(RG_ROWS * sizeof(int32_t));
    uint8_t* e = malloc(RG_ROWS);
    int16_t* da = malloc(RG_ROWS * sizeof(int16_t));
    int16_t* db = malloc(RG_ROWS * sizeof(int16_t));
    int16_t* de = malloc(RG_ROWS * sizeof(int16_t));
    if (!a || !b || !c || !d || !e || !da || !db || !de) rc = 1;

    for (int rg = 0; rg < NUM_RG && rc == 0; rg++) {
        int64_t na = 0, nb = 0, ne = 0;
        for (int64_t k = 0; k < RG_ROWS; k++) {
            int64_t r = (int64_t)rg * RG_ROWS + k;
            da[k] = a_null(r) ? 0 : 1;
            if (da[k]) a[na++] = a_val(r);
            db[k] = b_null(r) ? 0 : 1;
            if (db[k]) b[nb++] = b_val(r);
            c[k] = c_val(r);
            d[k] = d_val(r);
            de[k] = e_null(r) ? 0 : 1;
            if (de[k]) e[ne++] = e_val(r);
        }
        /* Dense value cursors per OPTIONAL column. */
        int64_t pa = 0, pb = 0, pe = 0;
        for (int64_t k = 0; k < RG_ROWS && rc == 0; k += PAGE_ROWS) {
            int64_t m = RG_ROWS - k < PAGE_ROWS ? RG_ROWS - k : PAGE_ROWS;
            int64_t ca = 0, cb = 0, ce = 0;
            for (int64_t q = k; q < k + m; q++) {
                ca += da[q];
                cb += db[q];
                ce += de[q];
            }
            if (carquet_writer_write_batch(w, 0, a + pa, m, da + k, NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 1, b + pb, m, db + k, NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 2, c + k, m, NULL, NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 3, d + k, m, NULL, NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 4, e + pe, m, de + k, NULL) != CARQUET_OK) {
                rc = 1;
            }
            pa += ca;
            pb += cb;
            pe += ce;
        }
        (void)na; (void)nb; (void)ne;
        if (rc == 0 && rg + 1 < NUM_RG && carquet_writer_new_row_group(w) != CARQUET_OK) rc = 1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    free(a); free(b); free(c); free(d); free(e); free(da); free(db); free(de);
    carquet_schema_free(schema);
    return rc;
}

static bool bit(const uint8_t* bm, int64_t j) {
    return (bm[j >> 3] >> (j & 7)) & 1;
}

/* Check one batch whose first row is global row `base`. */
static int check_batch(carquet_row_batch_t* batch, int64_t base, const char** why) {
    int64_t n = carquet_row_batch_num_rows(batch);
    const void* data[NUM_COLS];
    const uint8_t* nulls[NUM_COLS];
    for (int c = 0; c < NUM_COLS; c++) {
        int64_t count = 0;
        if (carquet_row_batch_column(batch, c, &data[c], &nulls[c], &count) != CARQUET_OK ||
            count != n) {
            *why = "column count mismatch";
            return 1;
        }
    }
    if (nulls[2] || nulls[3]) { *why = "REQUIRED column has a bitmap"; return 1; }
    if (!nulls[0] || !nulls[1] || !nulls[4]) { *why = "OPTIONAL column lacks a bitmap"; return 1; }

    const int64_t* a = data[0];
    const double* b = data[1];
    const int32_t* c = data[2];
    const int32_t* d = data[3];
    const uint8_t* e = data[4];
    for (int64_t j = 0; j < n; j++) {
        int64_t r = base + j;
        if (bit(nulls[0], j) != !a_null(r)) { *why = "a: validity"; return 1; }
        if (a[j] != (a_null(r) ? 0 : a_val(r))) { *why = "a: value"; return 1; }
        if (bit(nulls[1], j) != !b_null(r)) { *why = "b: validity"; return 1; }
        if (b[j] != (b_null(r) ? 0.0 : b_val(r))) { *why = "b: value"; return 1; }
        if (c[j] != c_val(r)) { *why = "c: value"; return 1; }
        if (d[j] != d_val(r)) { *why = "d: value"; return 1; }
        if (bit(nulls[4], j) != !e_null(r)) { *why = "e: validity"; return 1; }
        if (e[j] != (e_null(r) ? 0 : e_val(r))) { *why = "e: value"; return 1; }
    }
    return 0;
}

static int read_and_check(const char* path, bool use_mmap, int32_t batch_size,
                          bool with_filter, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = batch_size;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }

    /* c >= 0 holds for every row: the filter keeps everything. */
    int32_t zero = 0;
    carquet_filter_clause_t clause;
    memset(&clause, 0, sizeof(clause));
    clause.column_index = 2;
    clause.op = CARQUET_FILTER_GE;
    clause.value = &zero;
    clause.value_size = (int32_t)sizeof(zero);
    if (with_filter &&
        carquet_batch_reader_set_page_filter(br, &clause, 1) != CARQUET_OK) {
        carquet_batch_reader_free(br);
        carquet_reader_close(reader);
        TEST_FAIL(name, "set page filter");
    }

    int64_t seen = 0;
    const char* why = NULL;
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st;
    while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
        int bad = check_batch(batch, seen, &why);
        seen += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
        if (bad) break;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);

    if (why) {
        char msg[160];
        snprintf(msg, sizeof(msg), "%s (rows seen %lld)", why, (long long)seen);
        TEST_FAIL(name, msg);
    }
    if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) TEST_FAIL(name, "next failed");
    if (seen != TOTAL_ROWS) TEST_FAIL(name, "row count");
    TEST_PASS(name);
    return 0;
}

/* Export every batch to Arrow and check all arrays only once the batch
 * reader and the file reader are gone: fixed-width columns the pipeline
 * served are exported as views holding a reference on its row-group
 * buffers, which the reader hands over instead of reusing. */
static int read_arrow_kept(const char* path, bool use_mmap, int32_t batch_size,
                           bool expect_views, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = batch_size;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }
    struct ArrowArray* arrays = NULL;
    int64_t* bases = NULL;
    int32_t n_arr = 0, cap = 0, views = 0;
    const char* why = NULL;
    int64_t seen = 0;
    carquet_row_batch_t* batch = NULL;
    while (!why && carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        if (n_arr == cap) {
            int32_t ncap = cap ? cap * 2 : 16;
            struct ArrowArray* ga = realloc(arrays, (size_t)ncap * sizeof(*ga));
            if (ga) arrays = ga;
            int64_t* gb = realloc(bases, (size_t)ncap * sizeof(*gb));
            if (gb) bases = gb;
            if (!ga || !gb) { why = "alloc"; break; }
            cap = ncap;
        }
        if (carquet_arrow_export_batch(batch, carquet_reader_schema(reader), NULL,
                                       &arrays[n_arr], &err) != CARQUET_OK) {
            why = "arrow export";
        } else {
            /* A view (not a copy) of column c: the zero-copy export. */
            const void* cdata = NULL;
            const uint8_t* cnulls = NULL;
            int64_t cn = 0;
            if (carquet_row_batch_column(batch, 2, &cdata, &cnulls, &cn) == CARQUET_OK &&
                arrays[n_arr].children[2]->buffers[1] == cdata) {
                views++;
            }
            bases[n_arr++] = seen;
        }
        seen += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);

    int64_t checked = 0;
    for (int32_t k = 0; !why && k < n_arr; k++) {
        const struct ArrowArray* top = &arrays[k];
        if (top->n_children != NUM_COLS) { why = "child count"; break; }
        const struct ArrowArray* ca = top->children[0];
        const struct ArrowArray* cb = top->children[1];
        const struct ArrowArray* cc = top->children[2];
        const struct ArrowArray* cd = top->children[3];
        const struct ArrowArray* ce = top->children[4];
        const int64_t* a = ca->buffers[1];
        const double* b = cb->buffers[1];
        const int32_t* c = cc->buffers[1];
        const int32_t* d = cd->buffers[1];
        const uint8_t* e = ce->buffers[1];
        int64_t na = 0, nb = 0, ne = 0;
        for (int64_t j = 0; j < top->length && !why; j++) {
            int64_t r = bases[k] + j;
            if (bit(ca->buffers[0], j) != !a_null(r)) why = "a: validity";
            else if (!a_null(r) && a[j] != a_val(r)) why = "a: value";
            else if (bit(cb->buffers[0], j) != !b_null(r)) why = "b: validity";
            else if (!b_null(r) && b[j] != b_val(r)) why = "b: value";
            else if (c[j] != c_val(r)) why = "c: value";
            else if (d[j] != d_val(r)) why = "d: value";
            else if (bit(ce->buffers[0], j) != !e_null(r)) why = "e: validity";
            else if (!e_null(r) && bit(e, j) != (e_val(r) != 0)) why = "e: value";
            na += a_null(r);
            nb += b_null(r);
            ne += e_null(r);
        }
        if (!why && (ca->null_count != na || cb->null_count != nb || ce->null_count != ne ||
                     cc->null_count != 0 || cd->null_count != 0)) {
            why = "null_count";
        }
        checked += top->length;
    }
    for (int32_t k = 0; k < n_arr; k++) arrays[k].release(&arrays[k]);
    free(arrays);
    free(bases);
    if (why) TEST_FAIL(name, why);
    if (seen != TOTAL_ROWS || checked != TOTAL_ROWS) TEST_FAIL(name, "row count");
    /* Batches of at least a quarter row group are exported as views; smaller
     * ones are copied so they do not pin a whole row group's buffer. */
    if (expect_views ? views == 0 : views != 0) TEST_FAIL(name, "zero-copy export use");
    TEST_PASS(name);
    return 0;
}

/* preserve_dictionaries must expose indices for a dictionary column, with
 * or without mmap (the pipeline only holds physical values, so such a
 * reader takes the sequential path). Column c is REQUIRED dictionary. */
static int check_preserved_dictionary(const char* path, bool use_mmap, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");

    int32_t cols[1] = {2};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.column_indices = cols;
    cfg.num_columns = 1;
    cfg.preserve_dictionaries = true;
    cfg.batch_size = 65536;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }

    const char* why = NULL;
    int64_t seen = 0;
    carquet_row_batch_t* batch = NULL;
    while (!why && carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const uint32_t* idx = NULL;
        const uint8_t* nulls = NULL;
        int64_t n = 0;
        const uint8_t* dict = NULL;
        int32_t dict_count = 0;
        const uint32_t* dict_offsets = NULL;
        if (carquet_row_batch_column_dictionary(batch, 0, &idx, &nulls, &n, &dict,
                                                &dict_count, &dict_offsets) != CARQUET_OK) {
            why = "dictionary not preserved";
        } else {
            const int32_t* values = (const int32_t*)(const void*)dict;
            for (int64_t j = 0; j < n && !why; j++) {
                if (idx[j] >= (uint32_t)dict_count) why = "index out of range";
                else if (values[idx[j]] != c_val(seen + j)) why = "wrong value";
            }
        }
        seen += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);
    if (why) TEST_FAIL(name, why);
    if (seen != TOTAL_ROWS) TEST_FAIL(name, "row count");
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * BYTE_ARRAY columns in the pipeline
 * ============================================================================
 * s0: STRING OPTIONAL dictionary; s1: STRING OPTIONAL PLAIN, lengths 0..299
 * (empty strings included) and null for the whole second row group; s2:
 * STRING REQUIRED PLAIN; k: INT32 REQUIRED, constant per page (0, 1, 2, ...)
 * so an EQ filter on it selects every third page: several row ranges per
 * row group.
 */
#define S_COLS 4

static bool s0_null(int64_t r) { return r % 11 == 2; }
static bool s1_null(int64_t r) { return r % 9 == 4 || (r / RG_ROWS) == 1; }
static int32_t k_val(int64_t r) { return (int32_t)(((r % RG_ROWS) / PAGE_ROWS) % 3); }

/* Value bytes of row r of column `col` (0..2) into buf; returns the length. */
static int32_t s_val(int col, int64_t r, char* buf) {
    if (col == 0) return (int32_t)snprintf(buf, 32, "k%d", (int)(r % 97));
    if (col == 2) return (int32_t)snprintf(buf, 32, "req-%lld", (long long)r);
    /* Up to 299 bytes: page header statistics then outgrow the first window
     * a header is read with without mmap (re-read larger). */
    int32_t len = (int32_t)((r * 7) % 300);
    for (int32_t i = 0; i < len; i++) buf[i] = (char)('a' + (r + i) % 26);
    return len;
}

static bool s_null(int col, int64_t r) {
    return col == 0 ? s0_null(r) : col == 1 ? s1_null(r) : false;
}

static int write_string_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;
    carquet_logical_type_t str;
    memset(&str, 0, sizeof(str));
    str.id = CARQUET_LOGICAL_STRING;
    if (carquet_schema_add_column(schema, "s0", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "s1", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "s2", CARQUET_PHYSICAL_BYTE_ARRAY, &str,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "k", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return 1;
    }
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.max_rows_per_page = PAGE_ROWS;
    opts.write_page_index = true;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) {
        carquet_schema_free(schema);
        return 1;
    }
    int rc = 0;
    if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK ||
        carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN) != CARQUET_OK ||
        carquet_writer_set_column_encoding(w, 2, CARQUET_ENCODING_PLAIN) != CARQUET_OK) {
        rc = 1;
    }
    char* text = malloc((size_t)PAGE_ROWS * 3 * 320);
    carquet_byte_array_t* v[3];
    int16_t* def[3];
    for (int c = 0; c < 3; c++) {
        v[c] = malloc(PAGE_ROWS * sizeof(carquet_byte_array_t));
        def[c] = malloc(PAGE_ROWS * sizeof(int16_t));
        if (!v[c] || !def[c]) rc = 1;
    }
    int32_t* k = malloc(PAGE_ROWS * sizeof(int32_t));
    if (!text || !k) rc = 1;

    for (int rg = 0; rg < NUM_RG && rc == 0; rg++) {
        for (int64_t s = 0; s < RG_ROWS && rc == 0; s += PAGE_ROWS) {
            int64_t m = RG_ROWS - s < PAGE_ROWS ? RG_ROWS - s : PAGE_ROWS;
            char* tp = text;
            int64_t nv[3] = {0, 0, 0};
            for (int64_t q = 0; q < m; q++) {
                int64_t r = (int64_t)rg * RG_ROWS + s + q;
                k[q] = k_val(r);
                for (int c = 0; c < 3; c++) {
                    def[c][q] = s_null(c, r) ? 0 : 1;
                    if (!def[c][q]) continue;
                    int32_t len = s_val(c, r, tp);
                    v[c][nv[c]].data = (uint8_t*)tp;
                    v[c][nv[c]].length = len;
                    nv[c]++;
                    tp += len;
                }
            }
            if (carquet_writer_write_batch(w, 0, v[0], m, def[0], NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 1, v[1], m, def[1], NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 2, v[2], m, NULL, NULL) != CARQUET_OK ||
                carquet_writer_write_batch(w, 3, k, m, NULL, NULL) != CARQUET_OK) {
                rc = 1;
            }
        }
        if (rc == 0 && rg + 1 < NUM_RG && carquet_writer_new_row_group(w) != CARQUET_OK) rc = 1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    for (int c = 0; c < 3; c++) { free(v[c]); free(def[c]); }
    free(text);
    free(k);
    carquet_schema_free(schema);
    return rc;
}

/* Batch column `bc` holds string column `col`; check `n` rows whose global
 * row indices are rows[0..n) (or base + j when rows is NULL). */
static int check_string_column(const carquet_row_batch_t* batch, int bc, int col,
                               int64_t base, const int64_t* rows, int64_t n,
                               const char** why) {
    const void* data = NULL;
    const uint8_t* nulls = NULL;
    int64_t count = 0;
    if (carquet_row_batch_column(batch, bc, &data, &nulls, &count) != CARQUET_OK || count != n) {
        *why = "string column count";
        return 1;
    }
    if ((col == 2) != (nulls == NULL)) { *why = "string bitmap presence"; return 1; }
    const carquet_byte_array_t* ba = data;
    char buf[320];
    for (int64_t j = 0; j < n; j++) {
        int64_t r = rows ? rows[j] : base + j;
        bool present = nulls ? bit(nulls, j) : true;
        if (present == s_null(col, r)) { *why = "string validity"; return 1; }
        if (!present) {
            if (ba[j].length != 0) { *why = "null string slot not empty"; return 1; }
            continue;
        }
        int32_t len = s_val(col, r, buf);
        if (ba[j].length != len || (len > 0 && memcmp(ba[j].data, buf, (size_t)len) != 0)) {
            *why = "string value";
            return 1;
        }
    }
    return 0;
}

/* Arrow arrays kept past the batch reader's lifetime. */
typedef struct {
    struct ArrowArray arr;
    int64_t base;
} kept_array_t;

static int check_kept_arrays(const kept_array_t* kept, int32_t nkept, const char** why) {
    char buf[320];
    for (int32_t a = 0; a < nkept; a++) {
        const struct ArrowArray* top = &kept[a].arr;
        if (top->n_children != S_COLS) { *why = "arrow: child count"; return 1; }
        for (int c = 0; c < 3; c++) {
            const struct ArrowArray* ch = top->children[c];
            const uint8_t* valid = ch->buffers[0];
            const int32_t* off = ch->buffers[1];
            const uint8_t* bytes = ch->buffers[2];
            int64_t nulls = 0;
            for (int64_t j = 0; j < ch->length; j++) {
                int64_t r = kept[a].base + j;
                bool present = valid ? bit(valid, j) : true;
                if (present == s_null(c, r)) { *why = "arrow: validity"; return 1; }
                if (!present) {
                    nulls++;
                    if (off[j + 1] != off[j]) { *why = "arrow: null slot has bytes"; return 1; }
                    continue;
                }
                int32_t len = s_val(c, r, buf);
                if (off[j + 1] - off[j] != len ||
                    (len > 0 && memcmp(bytes + off[j], buf, (size_t)len) != 0)) {
                    *why = "arrow: string value";
                    return 1;
                }
            }
            if (ch->null_count != nulls) { *why = "arrow: null_count"; return 1; }
        }
        const struct ArrowArray* kc = top->children[3];
        const int32_t* kv = kc->buffers[1];
        for (int64_t j = 0; j < kc->length; j++) {
            if (kv[j] != k_val(kept[a].base + j)) { *why = "arrow: k value"; return 1; }
        }
    }
    return 0;
}

/* Read every row of the string file and check each batch before the next
 * call; with `arrow`, also export each batch and check all exported arrays
 * only after the batch reader and the file reader are gone. */
static int read_strings(const char* path, bool use_mmap, int32_t batch_size, bool arrow,
                        const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = batch_size;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }
    /* The writer may cut row groups smaller than asked, so batches can be
     * smaller than batch_size: grow the kept-array list as needed. */
    int32_t cap = 0;
    kept_array_t* kept = NULL;
    int32_t nkept = 0;

    const char* why = NULL;
    int64_t seen = 0;
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st;
    while (!why && (st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
        int64_t n = carquet_row_batch_num_rows(batch);
        for (int c = 0; c < 3 && !why; c++) {
            check_string_column(batch, c, c, seen, NULL, n, &why);
        }
        if (!why && arrow && nkept == cap) {
            int32_t ncap = cap ? cap * 2 : 16;
            kept_array_t* grown = realloc(kept, (size_t)ncap * sizeof(*kept));
            if (grown) {
                kept = grown;
                cap = ncap;
            } else {
                why = "alloc";
            }
        }
        if (!why && arrow) {
            if (carquet_arrow_export_batch(batch, carquet_reader_schema(reader), NULL,
                                           &kept[nkept].arr, &err) != CARQUET_OK) {
                why = "arrow export";
            } else {
                kept[nkept++].base = seen;
            }
        }
        seen += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);
    if (!why && arrow) check_kept_arrays(kept, nkept, &why);
    for (int32_t a = 0; a < nkept; a++) kept[a].arr.release(&kept[a].arr);
    free(kept);
    if (why) {
        char msg[160];
        snprintf(msg, sizeof(msg), "%s (rows seen %lld)", why, (long long)seen);
        TEST_FAIL(name, msg);
    }
    if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) TEST_FAIL(name, "next failed");
    if (seen != TOTAL_ROWS) TEST_FAIL(name, "row count");
    TEST_PASS(name);
    return 0;
}

/* REQUIRED-only projection (s2, k) under a page filter that keeps every
 * third page: the pipeline's filtered pre-read seeks between ranges, which
 * must not free the pages the strings of earlier ranges point into. */
static int read_strings_filtered(const char* path, bool use_mmap, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");
    int32_t cols[2] = {2, 3};
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.column_indices = cols;
    cfg.num_columns = 2;
    cfg.batch_size = 1 << 20;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }
    int32_t zero = 0;
    carquet_filter_clause_t clause;
    memset(&clause, 0, sizeof(clause));
    clause.column_index = 3;
    clause.op = CARQUET_FILTER_EQ;
    clause.value = &zero;
    clause.value_size = (int32_t)sizeof(zero);
    if (carquet_batch_reader_set_page_filter(br, &clause, 1) != CARQUET_OK) {
        carquet_batch_reader_free(br);
        carquet_reader_close(reader);
        TEST_FAIL(name, "set page filter");
    }
    int64_t* rows = malloc((size_t)TOTAL_ROWS * sizeof(int64_t));
    int64_t nrows = 0;
    for (int64_t r = 0; rows && r < TOTAL_ROWS; r++) {
        if (k_val(r) == 0) rows[nrows++] = r;
    }
    const char* why = rows ? NULL : "alloc";
    int64_t seen = 0;
    carquet_row_batch_t* batch = NULL;
    while (!why && carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        int64_t n = carquet_row_batch_num_rows(batch);
        if (seen + n > nrows) {
            why = "too many rows";
        } else {
            check_string_column(batch, 0, 2, 0, rows + seen, n, &why);
        }
        seen += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);
    free(rows);
    if (why) TEST_FAIL(name, why);
    if (seen != nrows) TEST_FAIL(name, "row count");
    TEST_PASS(name);
    return 0;
}

/* ============================================================================
 * Definition levels decoded as a validity bitmap
 * ============================================================================
 * Flat OPTIONAL columns (max_def == 1) decode their 1-bit level stream
 * straight into a page bitmap. First the bitmap helpers against naive
 * references over arbitrary (unaligned) offsets, in buffers of the exact
 * size they cover; then a file whose OPTIONAL columns of every value width
 * cycle through null shapes (all null, none, sparse, dense, alternating,
 * runs) read through the batch reader and through the column reader API
 * (int16 levels expanded from the bitmap, skips, odd read sizes), with V1
 * and V2 pages, compressed and not.
 */

static uint64_t ls_mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

static int test_bitmap_helpers(void) {
    const char* name = "bitmap_helpers_unaligned";
    uint64_t seed = 1;
    for (int iter = 0; iter < 20000; iter++) {
        int64_t src_at = (int64_t)(ls_mix(seed++) % 80);
        int64_t dst_at = (int64_t)(ls_mix(seed++) % 80);
        int64_t n = (int64_t)(ls_mix(seed++) % (iter % 4 == 0 ? 600 : 140));
        size_t src_bytes = (size_t)((src_at + n + 7) / 8);
        size_t dst_bytes = (size_t)((dst_at + n + 7) / 8);
        if (src_bytes == 0) src_bytes = 1;
        if (dst_bytes == 0) dst_bytes = 1;
        uint8_t* src = malloc(src_bytes);        /* exact: ASan flags over-reads */
        uint8_t* dst = malloc(dst_bytes);
        int16_t* lv = malloc(((size_t)n + 1) * sizeof(int16_t));
        if (!src || !dst || !lv) { free(src); free(dst); free(lv); TEST_FAIL(name, "alloc"); }
        int density = (int)(ls_mix(seed++) % 4);
        for (size_t b = 0; b < src_bytes; b++) {
            uint64_t r = ls_mix(seed++);
            src[b] = density == 0 ? 0 : density == 1 ? 0xFF : (uint8_t)r;
        }
        /* dst: random bits below dst_at, zero from dst_at on. */
        for (size_t b = 0; b < dst_bytes; b++) dst[b] = (uint8_t)ls_mix(seed++);
        for (int64_t i = dst_at; i < (int64_t)dst_bytes * 8; i++) dst[i >> 3] &= (uint8_t)~(1u << (i & 7));
        uint8_t* orig = malloc(dst_bytes);
        if (!orig) { free(src); free(dst); free(lv); TEST_FAIL(name, "alloc"); }
        memcpy(orig, dst, dst_bytes);
        int64_t ref_count = 0;
        for (int64_t i = 0; i < n; i++) ref_count += bit(src, src_at + i);

        if (carquet_bitmap_count_range(src, src_at, n) != ref_count) {
            free(src); free(dst); free(lv); TEST_FAIL(name, "count_range");
        }
        carquet_bitmap_copy_range(dst, dst_at, src, src_at, n);
        for (int64_t i = 0; i < (int64_t)dst_bytes * 8; i++) {
            bool want = i < dst_at ? bit(orig, i)
                      : i < dst_at + n ? bit(src, src_at + i - dst_at) : false;
            if (bit(dst, i) != want) {
                free(src); free(dst); free(lv); free(orig);
                TEST_FAIL(name, "copy_range");
            }
        }
        free(orig);
        carquet_bitmap_to_levels(src, src_at, n, lv);
        for (int64_t i = 0; i < n; i++) {
            if (lv[i] != (int16_t)bit(src, src_at + i)) {
                free(src); free(dst); free(lv); TEST_FAIL(name, "to_levels");
            }
        }
        for (int64_t i = 0; n > 0 && i < n; i += 1 + (int64_t)(ls_mix(seed++) % 64)) {
            int k = (int)(n - i < 64 ? n - i : 64);
            k = 1 + (int)(ls_mix(seed++) % (uint64_t)k);
            uint64_t w = carquet_bitmap_load_bits(src, src_at + i, k);
            for (int j = 0; j < 64; j++) {
                bool want = j < k && bit(src, src_at + i + j);
                if ((bool)((w >> j) & 1) != want) {
                    free(src); free(dst); free(lv); TEST_FAIL(name, "load_bits");
                }
            }
        }
        free(src); free(dst); free(lv);
    }
    TEST_PASS(name);
    return 0;
}

#define LS_RGS 2
#define LS_RG_ROWS 150000          /* >= 2 x 65536: pipeline segments */
#define LS_ROWS ((int64_t)LS_RGS * LS_RG_ROWS)
#define LS_PAGE_ROWS 3001
#define LS_COLS 7

/* Columns: 0 INT64 PLAIN, 1 INT32 dictionary, 2 DOUBLE PLAIN, 3 FLBA(16)
 * PLAIN, 4 FLBA(5) PLAIN, 5 BOOLEAN, 6 BYTE_ARRAY PLAIN; all OPTIONAL. */
static const size_t ls_width[LS_COLS] = {8, 4, 8, 16, 5, 1, sizeof(carquet_byte_array_t)};

/* Null shape of row r in column c: 4999-row segments cycling (shifted per
 * column) through all null, no null, sparse (7%), dense (50%), alternating
 * and runs of 37. */
static bool ls_null(int c, int64_t r) {
    uint64_t h = ls_mix((uint64_t)r * 8 + (uint64_t)c);
    switch ((r / 4999 + c) % 6) {
    case 0: return true;
    case 1: return false;
    case 2: return h % 100 < 7;
    case 3: return h & 1;
    case 4: return r & 1;
    default: return (r / 37) & 1;
    }
}

/* Value bytes of row r of fixed-width column c (c != 6). */
static void ls_val(int c, int64_t r, uint8_t* out) {
    switch (c) {
    case 0: { int64_t v = r * 3 + 11; memcpy(out, &v, 8); break; }
    case 1: { int32_t v = (int32_t)(r % 1000) - 500; memcpy(out, &v, 4); break; }
    case 2: { double v = (double)r * 0.25; memcpy(out, &v, 8); break; }
    case 3: case 4:
        for (size_t i = 0; i < ls_width[c]; i++) out[i] = (uint8_t)(ls_mix((uint64_t)r + i) >> 7);
        break;
    case 5: out[0] = r % 3 == 0; break;
    default: break;
    }
}

static int32_t ls_str(int64_t r, char* buf) {
    int32_t len = (int32_t)(r % 23);
    for (int32_t i = 0; i < len; i++) buf[i] = (char)('A' + (r + i) % 26);
    return len;
}

static int ls_write(const char* path, int page_version, carquet_compression_t codec) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;
    const carquet_physical_type_t types[LS_COLS] = {
        CARQUET_PHYSICAL_INT64, CARQUET_PHYSICAL_INT32, CARQUET_PHYSICAL_DOUBLE,
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
        CARQUET_PHYSICAL_BOOLEAN, CARQUET_PHYSICAL_BYTE_ARRAY};
    const char* names[LS_COLS] = {"i64", "i32", "dbl", "f16", "f5", "bool", "str"};
    for (int c = 0; c < LS_COLS; c++) {
        int32_t tl = types[c] == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY ? (int32_t)ls_width[c] : 0;
        if (carquet_schema_add_column(schema, names[c], types[c], NULL,
                                      CARQUET_REPETITION_OPTIONAL, tl, 0) != CARQUET_OK) {
            carquet_schema_free(schema);
            return 1;
        }
    }
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codec;
    opts.data_page_version = page_version;
    opts.max_rows_per_page = LS_PAGE_ROWS;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) {
        carquet_schema_free(schema);
        return 1;
    }
    int rc = 0;
    for (int c = 0; c < LS_COLS; c++) {
        carquet_encoding_t e = c == 1 ? CARQUET_ENCODING_RLE_DICTIONARY : CARQUET_ENCODING_PLAIN;
        if (carquet_writer_set_column_encoding(w, c, e) != CARQUET_OK) rc = 1;
    }
    uint8_t* vals = malloc((size_t)LS_PAGE_ROWS * 16);
    carquet_byte_array_t* ba = malloc((size_t)LS_PAGE_ROWS * sizeof(*ba));
    char* text = malloc((size_t)LS_PAGE_ROWS * 32);
    int16_t* def = malloc((size_t)LS_PAGE_ROWS * sizeof(int16_t));
    if (!vals || !ba || !text || !def) rc = 1;
    for (int rg = 0; rg < LS_RGS && rc == 0; rg++) {
        for (int64_t s = 0; s < LS_RG_ROWS && rc == 0; s += LS_PAGE_ROWS) {
            int64_t m = LS_RG_ROWS - s < LS_PAGE_ROWS ? LS_RG_ROWS - s : LS_PAGE_ROWS;
            for (int c = 0; c < LS_COLS && rc == 0; c++) {
                int64_t nd = 0;
                size_t tpos = 0;
                for (int64_t k = 0; k < m; k++) {
                    int64_t r = (int64_t)rg * LS_RG_ROWS + s + k;
                    def[k] = ls_null(c, r) ? 0 : 1;
                    if (!def[k]) continue;
                    if (c == 6) {
                        int32_t len = ls_str(r, text + tpos);
                        ba[nd].data = (uint8_t*)text + tpos;
                        ba[nd].length = len;
                        tpos += (size_t)len;
                    } else {
                        ls_val(c, r, vals + (size_t)nd * ls_width[c]);
                    }
                    nd++;
                }
                const void* src = c == 6 ? (const void*)ba : (const void*)vals;
                if (carquet_writer_write_batch(w, c, src, m, def, NULL) != CARQUET_OK) rc = 1;
            }
        }
        if (rc == 0 && rg + 1 < LS_RGS && carquet_writer_new_row_group(w) != CARQUET_OK) rc = 1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    free(vals); free(ba); free(text); free(def);
    carquet_schema_free(schema);
    return rc;
}

/* Row r of column c is `cell` (width ls_width[c]); `present` its validity. */
static bool ls_cell_ok(int c, int64_t r, const uint8_t* cell, bool present) {
    if (present == ls_null(c, r)) return false;
    uint8_t want[32];
    memset(want, 0, sizeof(want));
    if (!present) {
        /* Null rows read as zero bytes in batches. */
        return memcmp(cell, want, ls_width[c]) == 0;
    }
    if (c == 6) {
        const carquet_byte_array_t* v = (const carquet_byte_array_t*)(const void*)cell;
        char buf[32];
        int32_t len = ls_str(r, buf);
        return v->length == len && (len == 0 || memcmp(v->data, buf, (size_t)len) == 0);
    }
    ls_val(c, r, want);
    return memcmp(cell, want, ls_width[c]) == 0;
}

static int ls_read_batches(const char* path, bool use_mmap, int32_t batch_size, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = batch_size;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        carquet_reader_close(reader);
        TEST_FAIL(name, "batch reader create");
    }
    const char* why = NULL;
    int64_t seen = 0;
    carquet_row_batch_t* batch = NULL;
    while (!why && carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        int64_t n = carquet_row_batch_num_rows(batch);
        for (int c = 0; c < LS_COLS && !why; c++) {
            const void* data = NULL;
            const uint8_t* nulls = NULL;
            int64_t cnt = 0;
            if (carquet_row_batch_column(batch, c, &data, &nulls, &cnt) != CARQUET_OK ||
                cnt != n || !nulls) {
                why = "column";
                break;
            }
            for (int64_t j = 0; j < n; j++) {
                if (!ls_cell_ok(c, seen + j, (const uint8_t*)data + (size_t)j * ls_width[c],
                                bit(nulls, j))) {
                    why = "cell";
                    printf("  col %d row %lld\n", c, (long long)(seen + j));
                    break;
                }
            }
        }
        seen += n;
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(reader);
    if (why) TEST_FAIL(name, why);
    if (seen != LS_ROWS) TEST_FAIL(name, "row count");
    TEST_PASS(name);
    return 0;
}

/* carquet_column_read_batch with int16 levels (expanded from the page
 * bitmaps) and dense values, in odd read sizes interleaved with skips that
 * land mid-page. */
static int ls_read_columns(const char* path, bool use_mmap, const char* name) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = use_mmap;
    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) TEST_FAIL(name, "open");
    enum { MAXR = 5000 };
    uint8_t* vals = malloc((size_t)MAXR * 16);
    int16_t* def = malloc((size_t)MAXR * sizeof(int16_t));
    const char* why = vals && def ? NULL : "alloc";
    uint64_t seed = 99;
    for (int rg = 0; rg < LS_RGS && !why; rg++) {
        for (int c = 0; c < LS_COLS && !why; c++) {
            carquet_column_reader_t* cr = carquet_reader_get_column(reader, rg, c, &err);
            if (!cr) { why = "get column"; break; }
            int64_t r = 0;
            while (r < LS_RG_ROWS && !why) {
                uint64_t op = ls_mix(seed++);
                int64_t k = 1 + (int64_t)((op >> 8) % (op & 1 ? 37 : MAXR));
                if (k > LS_RG_ROWS - r) k = LS_RG_ROWS - r;
                if (op % 5 == 0) {
                    if (carquet_column_skip(cr, k) != k) why = "skip";
                    r += k;
                    continue;
                }
                int64_t got = carquet_column_read_batch(cr, vals, k, def, NULL);
                if (got != k) { why = "read count"; break; }
                int64_t d = 0;
                for (int64_t j = 0; j < k; j++) {
                    int64_t row = (int64_t)rg * LS_RG_ROWS + r + j;
                    if (def[j] != (ls_null(c, row) ? 0 : 1)) { why = "level"; break; }
                    if (!def[j]) continue;
                    if (!ls_cell_ok(c, row, vals + (size_t)d * ls_width[c], true)) {
                        why = "value";
                        break;
                    }
                    d++;
                }
                r += k;
            }
            carquet_column_reader_free(cr);
        }
    }
    free(vals);
    free(def);
    carquet_reader_close(reader);
    if (why) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* Page decode, bitmap levels vs int16 levels: the same V1 (RLE and legacy
 * MSB-first BIT_PACKED level sections) and V2 pages decoded once into int16
 * levels (def_levels passed) and once into the reader's def bitmap, for
 * random null shapes and sizes, intact and with a truncated level section. */
extern carquet_status_t carquet_read_data_page_v1(
    carquet_column_reader_t* reader, const uint8_t* page_data, size_t page_size,
    const parquet_data_page_header_t* header, void* values, int64_t max_values,
    int16_t* def_levels, int16_t* rep_levels, int64_t* values_read, carquet_error_t* error);
extern carquet_status_t carquet_read_data_page_v2(
    carquet_column_reader_t* reader, const uint8_t* page_data, size_t page_size,
    const parquet_data_page_header_v2_t* header, void* values, int64_t max_values,
    int16_t* def_levels, int16_t* rep_levels, int64_t* values_read, carquet_error_t* error);

/* Hybrid RLE/bit-packed stream of 1-bit `lv[0, n)`, runs and groups mixed at
 * random. Returns its length. */
static size_t pg_hybrid(const int16_t* lv, int32_t n, uint8_t* out, uint64_t* seed) {
    size_t pos = 0;
    int32_t i = 0;
    while (i < n) {
        uint64_t r = ls_mix((*seed)++);
        if (r % 3 == 0) {
            int32_t g = 1 + (int32_t)((r >> 8) % 4);
            out[pos++] = (uint8_t)((g << 1) | 1);
            for (int32_t q = 0; q < g * 8; q += 8) {
                uint8_t b = 0;
                for (int j = 0; j < 8; j++) {
                    if (i + q + j < n && lv[i + q + j]) b |= (uint8_t)(1u << j);
                }
                out[pos++] = b;
            }
            i += g * 8;
        } else {
            int32_t len = 1;
            while (i + len < n && lv[i + len] == lv[i] && len < 1000) len++;
            uint32_t h = (uint32_t)len << 1;
            while (h >= 0x80) { out[pos++] = (uint8_t)(h | 0x80); h >>= 7; }
            out[pos++] = (uint8_t)h;
            out[pos++] = (uint8_t)lv[i];
            i += len;
        }
    }
    return pos;
}

static int test_page_decode_bits_vs_levels(const char* path) {
    const char* name = "page_decode_bits_vs_levels";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) TEST_FAIL(name, "open");
    carquet_column_reader_t* cr = carquet_reader_get_column(reader, 0, 1, &err);  /* INT32 OPTIONAL */
    if (!cr || !carquet_column_def_as_bits(cr)) {
        if (cr) carquet_column_reader_free(cr);
        carquet_reader_close(reader);
        TEST_FAIL(name, "column reader");
    }
    enum { MAXN = 3000 };
    uint8_t* page = malloc(16 + MAXN * 4 + MAXN * 4);
    int16_t* lv = malloc(MAXN * sizeof(int16_t));
    int16_t* ref_lv = malloc(MAXN * sizeof(int16_t));
    int32_t* ref_vals = malloc(MAXN * sizeof(int32_t));
    int32_t* vals = malloc(MAXN * sizeof(int32_t));
    cr->decoded_def_bits = carquet_mem_malloc(MAXN / 8 + 16);
    cr->def_bits_capacity = cr->decoded_def_bits ? MAXN / 8 + 16 : 0;
    const char* why = (page && lv && ref_lv && ref_vals && vals && cr->decoded_def_bits) ? NULL : "alloc";
    uint64_t seed = 4242;
    for (int iter = 0; iter < 3000 && !why; iter++) {
        int32_t n = 1 + (int32_t)(ls_mix(seed++) % MAXN);
        int shape = (int)(ls_mix(seed++) % 4);
        int32_t nn = 0;
        for (int32_t i = 0; i < n; i++) {
            uint64_t h = ls_mix(seed++);
            lv[i] = shape == 0 ? (h % 100 >= 7) : shape == 1 ? (int16_t)(h & 1)
                  : shape == 2 ? (int16_t)((i / 50) & 1) : (int16_t)(h % 100 < 3);
            nn += lv[i];
        }
        int kind = iter % 3;   /* 0: V1 RLE, 1: V1 BIT_PACKED, 2: V2 */
        size_t pos = 0;
        size_t lvl_bytes;
        if (kind == 1) {
            lvl_bytes = ((size_t)n + 7) / 8;
            memset(page, 0, lvl_bytes);
            for (int32_t i = 0; i < n; i++) {
                if (lv[i]) page[i >> 3] |= (uint8_t)(0x80u >> (i & 7));
            }
            pos = lvl_bytes;
        } else if (kind == 0) {
            lvl_bytes = pg_hybrid(lv, n, page + 4, &seed);
            page[0] = (uint8_t)lvl_bytes; page[1] = (uint8_t)(lvl_bytes >> 8);
            page[2] = 0; page[3] = 0;
            pos = 4 + lvl_bytes;
        } else {
            lvl_bytes = pg_hybrid(lv, n, page, &seed);
            pos = lvl_bytes;
        }
        for (int32_t i = 0; i < nn; i++) {
            int32_t v = (int32_t)(ls_mix((uint64_t)i + seed) & 0x7FFFFFFF);
            memcpy(page + pos + (size_t)i * 4, &v, 4);
        }
        size_t page_size = pos + (size_t)nn * 4;
        /* Second pass: the level section cut short (a decode error). */
        for (int cut = 0; cut < 2 && !why; cut++) {
            parquet_data_page_header_t h1;
            parquet_data_page_header_v2_t h2;
            memset(&h1, 0, sizeof(h1));
            memset(&h2, 0, sizeof(h2));
            h1.num_values = n;
            h1.encoding = CARQUET_ENCODING_PLAIN;
            h1.definition_level_encoding = kind == 1 ? CARQUET_ENCODING_BIT_PACKED : CARQUET_ENCODING_RLE;
            h1.repetition_level_encoding = CARQUET_ENCODING_RLE;
            h2.num_values = n;
            h2.num_rows = n;
            h2.num_nulls = n - nn;
            h2.encoding = CARQUET_ENCODING_PLAIN;
            h2.definition_levels_byte_length = (int32_t)lvl_bytes;
            size_t psize = page_size;
            if (cut) {
                if (kind == 0) {
                    page[0] = (uint8_t)(lvl_bytes / 2); page[1] = (uint8_t)((lvl_bytes / 2) >> 8);
                } else if (kind == 1) {
                    psize = lvl_bytes - 1;   /* fewer bytes than ceil(n/8) */
                } else {
                    h2.definition_levels_byte_length = (int32_t)(lvl_bytes / 2);
                }
            }
            int64_t got_ref = -1, got = -1;
            carquet_error_t e1 = CARQUET_ERROR_INIT, e2 = CARQUET_ERROR_INIT;
            carquet_status_t s1, s2;
            if (kind == 2) {
                s1 = carquet_read_data_page_v2(cr, page, psize, &h2, ref_vals, n, ref_lv, NULL, &got_ref, &e1);
                s2 = carquet_read_data_page_v2(cr, page, psize, &h2, vals, n, NULL, NULL, &got, &e2);
            } else {
                s1 = carquet_read_data_page_v1(cr, page, psize, &h1, ref_vals, n, ref_lv, NULL, &got_ref, &e1);
                s2 = carquet_read_data_page_v1(cr, page, psize, &h1, vals, n, NULL, NULL, &got, &e2);
            }
            if (cut && kind == 0) {
                page[0] = (uint8_t)lvl_bytes; page[1] = (uint8_t)(lvl_bytes >> 8);
            }
            if ((s1 == CARQUET_OK) != (s2 == CARQUET_OK)) { why = "status differs"; break; }
            if (s1 != CARQUET_OK) {
                if (!cut) why = "intact page failed";
                continue;
            }
            if (got != got_ref || got != n) { why = "count"; break; }
            int32_t dense = 0;
            for (int32_t i = 0; i < n; i++) {
                if (ref_lv[i] != lv[i]) { why = "reference levels"; break; }
                if (bit(cr->decoded_def_bits, i) != (ref_lv[i] == 1)) { why = "bit vs level"; break; }
                dense += ref_lv[i];
            }
            if (!why && memcmp(vals, ref_vals, (size_t)dense * 4) != 0) why = "values";
        }
    }
    free(page); free(lv); free(ref_lv); free(ref_vals); free(vals);
    carquet_column_reader_free(cr);
    carquet_reader_close(reader);
    if (why) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

static int run_level_shape_tests(void) {
    int failures = 0;
    struct { int version; carquet_compression_t codec; const char* tag; } files[] = {
        {1, CARQUET_COMPRESSION_UNCOMPRESSED, "v1_none"},
        {1, CARQUET_COMPRESSION_ZSTD, "v1_zstd"},
        {2, CARQUET_COMPRESSION_UNCOMPRESSED, "v2_none"},
        {2, CARQUET_COMPRESSION_SNAPPY, "v2_snappy"},
    };
    for (size_t f = 0; f < sizeof(files) / sizeof(files[0]); f++) {
        char path[512];
        char name[96];
        carquet_test_temp_path(path, sizeof(path), "batch_pipeline_levels");
        if (ls_write(path, files[f].version, files[f].codec) != 0) {
            printf("[FAIL] write level-shape file %s\n", files[f].tag);
            remove(path);
            failures++;
            continue;
        }
        const int32_t sizes[] = {65536, 999, 64, 7};
        for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
            for (int mm = 1; mm >= 0; mm--) {
                if (sizes[s] == 7 && mm == 0) continue;   /* slow, one mode is enough */
                snprintf(name, sizeof(name), "levels_%s_%s_batch_%d", files[f].tag,
                         mm ? "mmap" : "pread", (int)sizes[s]);
                failures += ls_read_batches(path, mm != 0, sizes[s], name);
            }
        }
        if (f == 0) failures += test_page_decode_bits_vs_levels(path);
        snprintf(name, sizeof(name), "levels_%s_column_reader_mmap", files[f].tag);
        failures += ls_read_columns(path, true, name);
        snprintf(name, sizeof(name), "levels_%s_column_reader_pread", files[f].tag);
        failures += ls_read_columns(path, false, name);
        remove(path);
    }
    return failures;
}

int main(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "batch_pipeline");
    if (write_file(path) != 0) {
        printf("[FAIL] write test file\n");
        remove(path);
        return 1;
    }

    int failures = 0;
    failures += read_and_check(path, true, 65536, false, "mmap_batch_65536");
    failures += read_and_check(path, true, 999, false, "mmap_batch_999_shifted_bitmaps");
    failures += read_and_check(path, false, 65536, false, "pread_batch_65536");
    failures += read_and_check(path, false, 1001, false, "pread_batch_1001");
    failures += read_and_check(path, true, 4096, true, "mmap_optional_with_page_filter");
    failures += check_preserved_dictionary(path, true, "mmap_preserve_dictionaries");
    failures += check_preserved_dictionary(path, false, "pread_preserve_dictionaries");
    failures += read_arrow_kept(path, true, 1 << 20, true, "mmap_arrow_kept_row_group_batches");
    failures += read_arrow_kept(path, true, 65536, true, "mmap_arrow_kept_batch_65536");
    failures += read_arrow_kept(path, true, 999, false, "mmap_arrow_kept_batch_999");
    failures += read_arrow_kept(path, false, 100000, true, "pread_arrow_kept_batch_100000");
    remove(path);

    carquet_test_temp_path(path, sizeof(path), "batch_pipeline_strings");
    if (write_string_file(path) != 0) {
        printf("[FAIL] write string test file\n");
        remove(path);
        return 1;
    }
    failures += read_strings(path, true, 65536, false, "strings_mmap_batch_65536");
    failures += read_strings(path, true, 999, false, "strings_mmap_batch_999");
    failures += read_strings(path, false, 65536, false, "strings_pread_batch_65536");
    failures += read_strings(path, false, 4093, false, "strings_pread_batch_4093");
    failures += read_strings(path, true, 1 << 20, true, "strings_mmap_arrow_kept_past_reader");
    failures += read_strings(path, false, 50000, true, "strings_pread_arrow_kept_past_reader");
    failures += read_strings_filtered(path, true, "strings_required_filtered_mmap");
    failures += read_strings_filtered(path, false, "strings_required_filtered_pread");

    remove(path);

    failures += test_bitmap_helpers();
    failures += run_level_shape_tests();
    return failures ? 1 : 0;
}
