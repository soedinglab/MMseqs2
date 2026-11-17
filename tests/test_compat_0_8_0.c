/**
 * @file test_compat_0_8_0.c
 * @brief Regression tests for the cross-library compatibility fixes of v0.8.0
 *        (found by reading carquet's files with pyarrow / DuckDB / polars /
 *        fastparquet and theirs, plus apache/parquet-testing, with carquet).
 *
 *  - decimal_scale_zero        DECIMAL(p, 0) wrote `precision` without `scale`
 *                              in the SchemaElement; DuckDB, polars and
 *                              fastparquet refuse such a file
 *  - dict_page_no_encodings    a chunk with an empty `encodings` list whose
 *                              dictionary page leads the data pages (no
 *                              dictionary_page_offset) could not be read
 *  - required_leaf_null_struct a REQUIRED leaf under a null OPTIONAL struct
 *                              was exported with nulls in a non-nullable field
 *  - legacy_lists              repeated fields outside a LIST group, LIST
 *                              groups whose repeated node is the element
 *                              (repeated primitive, "array" group) and a MAP
 *                              without values failed to assemble
 *  - ref_export_batch_offset   a whole-array view of a batch that carries an
 *                              offset kept the offset on the top-level struct
 *  - arrow_extension_logical   arrow.uuid / arrow.json Arrow fields were
 *                              written without the UUID / JSON logical type
 *  - utf8_truncated_stats      truncated STRING min / max were cut inside a
 *                              character (pyarrow cannot decode them)
 */

#include <carquet/carquet.h>
#include "reader/reader_internal.h"
#include "test_helpers.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const uint8_t* find_bytes(const uint8_t* hay, size_t n, const uint8_t* needle, size_t m) {
    if (m == 0 || n < m) return NULL;
    for (size_t i = 0; i + m <= n; i++) {
        if (hay[i] == needle[0] && memcmp(hay + i, needle, m) == 0) return hay + i;
    }
    return NULL;
}

static bool bit_set(const void* bitmap, int64_t i) {
    return (((const uint8_t*)bitmap)[i >> 3] >> (i & 7)) & 1;
}

/* ------------------------------------------------------------------------- */

static int test_decimal_scale_zero(void) {
    const char* name = "decimal_scale_zero";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_DECIMAL;
    lt.params.decimal.precision = 5;
    lt.params.decimal.scale = 0;
    if (carquet_schema_add_column(schema, "d", CARQUET_PHYSICAL_INT32, &lt,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "add_column");

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    if (!w) TEST_FAIL(name, "create");
    int32_t v[3] = {1, 2, 3};
    if (carquet_writer_write_batch(w, 0, v, 3, NULL, NULL) != CARQUET_OK) TEST_FAIL(name, "write");
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    void* buf = NULL; size_t size = 0;
    if (carquet_writer_get_buffer(w, &buf, &size) != CARQUET_OK) TEST_FAIL(name, "get_buffer");

    /* SchemaElement, Thrift compact: field 6 converted_type = DECIMAL (5),
     * field 7 scale = 0, field 8 precision = 5. Without the fix the scale
     * field is absent: 25 0A 25 0A. */
    static const uint8_t with_scale[] = {0x25, 0x0A, 0x15, 0x00, 0x15, 0x0A};
    static const uint8_t without_scale[] = {0x25, 0x0A, 0x25, 0x0A};
    bool ok = find_bytes((const uint8_t*)buf, size, with_scale, sizeof(with_scale)) != NULL;
    bool bad = find_bytes((const uint8_t*)buf, size, without_scale, sizeof(without_scale)) != NULL;

    carquet_reader_t* r = carquet_reader_open_buffer(buf, size, NULL, &err);
    bool reads = false;
    if (r) {
        const carquet_schema_t* rs = carquet_reader_schema(r);
        const parquet_schema_element_t* el = &rs->elements[1];
        reads = el->has_logical_type && el->logical_type.id == CARQUET_LOGICAL_DECIMAL &&
                el->logical_type.params.decimal.precision == 5 &&
                el->logical_type.params.decimal.scale == 0;
        carquet_reader_close(r);
    }
    carquet_free(buf);
    carquet_schema_free(schema);
    if (!ok || bad) TEST_FAIL(name, "scale 0 is not written next to precision");
    if (!reads) TEST_FAIL(name, "DECIMAL(5, 0) does not read back");
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

static int test_dict_page_no_encodings(void) {
    const char* name = "dict_page_no_encodings";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "carquet_compat_dict_noenc.parquet");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (carquet_schema_add_column(schema, "s", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "add_column");
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) TEST_FAIL(name, "create");
    if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK)
        TEST_FAIL(name, "set encoding");
    enum { N = 600 };
    static const char* words[3] = {"alpha", "beta", "gamma"};
    carquet_byte_array_t* vals = (carquet_byte_array_t*)malloc(N * sizeof(*vals));
    for (int i = 0; i < N; i++) {
        vals[i].data = (uint8_t*)words[i % 3];
        vals[i].length = (int32_t)strlen(words[i % 3]);
    }
    if (carquet_writer_write_batch(w, 0, vals, N, NULL, NULL) != CARQUET_OK) TEST_FAIL(name, "write");
    free(vals);
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    carquet_schema_free(schema);

    for (int use_mmap = 0; use_mmap < 2; use_mmap++) {
        carquet_reader_options_t ro;
        carquet_reader_options_init(&ro);
        ro.use_mmap = use_mmap != 0;
        carquet_reader_t* r = carquet_reader_open(path, &ro, &err);
        if (!r) TEST_FAIL(name, "open");
        /* The layout of parquet-testing's nation.dict-malformed.parquet (old
         * parquet-mr): no dictionary_page_offset, data_page_offset on the
         * dictionary page, and no encodings listed. */
        parquet_column_metadata_t* m = &r->metadata.row_groups[0].columns[0].metadata;
        if (!m->has_dictionary_page_offset) {
            carquet_reader_close(r);
            TEST_FAIL(name, "expected a dictionary-encoded chunk");
        }
        m->data_page_offset = m->dictionary_page_offset;
        m->has_dictionary_page_offset = false;
        m->dictionary_page_offset = 0;
        m->num_encodings = 0;

        carquet_column_reader_t* cr = carquet_reader_get_column(r, 0, 0, &err);
        if (!cr) { carquet_reader_close(r); TEST_FAIL(name, "get_column"); }
        carquet_byte_array_t out[N];
        int64_t got = 0;
        while (got < N) {
            int64_t n = carquet_column_read_batch(cr, out + got, N - got, NULL, NULL);
            if (n <= 0) break;
            for (int64_t i = got; i < got + n; i++) {
                const char* want = words[i % 3];
                if (out[i].length != (int32_t)strlen(want) ||
                    memcmp(out[i].data, want, (size_t)out[i].length) != 0) {
                    carquet_column_reader_free(cr);
                    carquet_reader_close(r);
                    TEST_FAIL(name, "wrong value");
                }
            }
            got += n;
        }
        carquet_column_reader_free(cr);
        carquet_reader_close(r);
        if (got != N) TEST_FAIL(name, "chunk with an empty encodings list could not be read");
    }
    remove(path);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

/* Read row group 0 of an in-memory file as a nested Arrow array. */
static carquet_status_t read_arrow_buffer(const void* buf, size_t size, struct ArrowArray* arr,
                                          struct ArrowSchema* sch, carquet_reader_t** reader_out) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(buf, size, NULL, &err);
    if (!r) return CARQUET_ERROR_FILE_OPEN;
    carquet_status_t st = carquet_reader_read_arrow(r, 0, sch, arr, &err);
    *reader_out = r;
    return st;
}

static int test_required_leaf_null_struct(void) {
    const char* name = "required_leaf_null_struct";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    /* optional group s { required int32 a; required binary b; } */
    int32_t g = carquet_schema_add_group(schema, "s", CARQUET_REPETITION_OPTIONAL, 0);
    if (g < 0) TEST_FAIL(name, "add_group");
    if (carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, g) != CARQUET_OK ||
        carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, g) != CARQUET_OK)
        TEST_FAIL(name, "add_column");
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    if (!w) TEST_FAIL(name, "create");
    /* rows: {7,"x"}, null, {9,"yz"} */
    int16_t def[3] = {1, 0, 1};
    int32_t a[2] = {7, 9};
    carquet_byte_array_t b[2] = {{(uint8_t*)"x", 1}, {(uint8_t*)"yz", 2}};
    if (carquet_writer_write_batch(w, 0, a, 3, def, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, b, 3, def, NULL) != CARQUET_OK)
        TEST_FAIL(name, "write");
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    void* buf = NULL; size_t size = 0;
    if (carquet_writer_get_buffer(w, &buf, &size) != CARQUET_OK) TEST_FAIL(name, "get_buffer");

    struct ArrowArray arr; struct ArrowSchema sch; carquet_reader_t* r = NULL;
    if (read_arrow_buffer(buf, size, &arr, &sch, &r) != CARQUET_OK) TEST_FAIL(name, "read_arrow");
    const struct ArrowArray* s = arr.children[0];
    const struct ArrowSchema* ss = sch.children[0];
    int failed = 0;
    const char* why = "";
    if (!(ss->flags & ARROW_FLAG_NULLABLE) || s->null_count != 1 || bit_set(s->buffers[0], 1)) {
        failed = 1; why = "struct validity";
    }
    for (int c = 0; c < 2 && !failed; c++) {
        if (ss->children[c]->flags & ARROW_FLAG_NULLABLE) { failed = 1; why = "child declared nullable"; }
        else if (s->children[c]->null_count != 0 || s->children[c]->buffers[0] != NULL) {
            failed = 1; why = "non-nullable child carries nulls";
        } else if (s->children[c]->length != 3) { failed = 1; why = "child length"; }
    }
    if (!failed) {
        const int32_t* av = (const int32_t*)s->children[0]->buffers[1];
        const int32_t* bo = (const int32_t*)s->children[1]->buffers[1];
        if (av[0] != 7 || av[2] != 9 || bo[0] != 0 || bo[1] != 1 || bo[2] != 1 || bo[3] != 3) {
            failed = 1; why = "values";
        }
    }
    arr.release(&arr); sch.release(&sch);
    carquet_reader_close(r);
    carquet_free(buf);
    carquet_schema_free(schema);
    if (failed) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

static void annotate(carquet_schema_t* s, int32_t e, carquet_logical_type_id_t id) {
    memset(&s->elements[e].logical_type, 0, sizeof(s->elements[e].logical_type));
    s->elements[e].logical_type.id = id;
    s->elements[e].has_logical_type = true;
}

static int test_legacy_lists(void) {
    const char* name = "legacy_lists";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    /*
     * repeated int32 bare;                               (leaf 0)
     * optional group ph { repeated group phone {         (unannotated repeated group)
     *     required int64 number; } }                     (leaf 1)
     * required group ll (LIST) { repeated group array (LIST) {
     *     repeated int32 array; } }                      (leaf 2, old_list_structure)
     * optional group prim (LIST) { repeated int32 item; }   (leaf 3)
     * required group set (MAP) { repeated group key_value {
     *     required int32 key; } }                        (leaf 4, map without values)
     */
    int ok = carquet_schema_add_column(schema, "bare", CARQUET_PHYSICAL_INT32, NULL,
                                       CARQUET_REPETITION_REPEATED, 0, 0) == CARQUET_OK;
    int32_t ph = carquet_schema_add_group(schema, "ph", CARQUET_REPETITION_OPTIONAL, 0);
    int32_t phone = carquet_schema_add_group(schema, "phone", CARQUET_REPETITION_REPEATED, ph);
    ok = ok && ph > 0 && phone > 0 &&
         carquet_schema_add_column(schema, "number", CARQUET_PHYSICAL_INT64, NULL,
                                   CARQUET_REPETITION_REQUIRED, 0, phone) == CARQUET_OK;
    int32_t ll = carquet_schema_add_group(schema, "ll", CARQUET_REPETITION_REQUIRED, 0);
    int32_t ll_arr = carquet_schema_add_group(schema, "array", CARQUET_REPETITION_REPEATED, ll);
    ok = ok && ll > 0 && ll_arr > 0 &&
         carquet_schema_add_column(schema, "array", CARQUET_PHYSICAL_INT32, NULL,
                                   CARQUET_REPETITION_REPEATED, 0, ll_arr) == CARQUET_OK;
    int32_t prim = carquet_schema_add_group(schema, "prim", CARQUET_REPETITION_OPTIONAL, 0);
    ok = ok && prim > 0 &&
         carquet_schema_add_column(schema, "item", CARQUET_PHYSICAL_INT32, NULL,
                                   CARQUET_REPETITION_REPEATED, 0, prim) == CARQUET_OK;
    int32_t set = carquet_schema_add_group(schema, "set", CARQUET_REPETITION_REQUIRED, 0);
    int32_t kv = carquet_schema_add_group(schema, "key_value", CARQUET_REPETITION_REPEATED, set);
    ok = ok && set > 0 && kv > 0 &&
         carquet_schema_add_column(schema, "key", CARQUET_PHYSICAL_INT32, NULL,
                                   CARQUET_REPETITION_REQUIRED, 0, kv) == CARQUET_OK;
    if (!ok) TEST_FAIL(name, "schema");
    annotate(schema, ll, CARQUET_LOGICAL_LIST);
    annotate(schema, ll_arr, CARQUET_LOGICAL_LIST);
    annotate(schema, prim, CARQUET_LOGICAL_LIST);
    annotate(schema, set, CARQUET_LOGICAL_MAP);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    if (!w) TEST_FAIL(name, "create");

    /* 3 rows.
     * bare: [1, 2], [], [3] */
    { int32_t v[] = {1, 2, 3}; int16_t d[] = {1, 1, 0, 1}; int16_t r[] = {0, 1, 0, 0};
      if (carquet_writer_write_batch(w, 0, v, 4, d, r) != CARQUET_OK) TEST_FAIL(name, "write bare"); }
    /* ph: null, {phone: []}, {phone: [{10}, {20}]}      (max_def 2) */
    { int64_t v[] = {10, 20}; int16_t d[] = {0, 1, 2, 2}; int16_t r[] = {0, 0, 0, 1};
      if (carquet_writer_write_batch(w, 1, v, 4, d, r) != CARQUET_OK) TEST_FAIL(name, "write ph"); }
    /* ll: [[1, 2], [3]], [], [[], [4]]                  (max_def 2, max_rep 2) */
    { int32_t v[] = {1, 2, 3, 4}; int16_t d[] = {2, 2, 2, 0, 1, 2}; int16_t r[] = {0, 2, 1, 0, 0, 1};
      if (carquet_writer_write_batch(w, 2, v, 6, d, r) != CARQUET_OK) TEST_FAIL(name, "write ll"); }
    /* prim: [5, 6], null, []                            (max_def 2) */
    { int32_t v[] = {5, 6}; int16_t d[] = {2, 2, 0, 1}; int16_t r[] = {0, 1, 0, 0};
      if (carquet_writer_write_batch(w, 3, v, 4, d, r) != CARQUET_OK) TEST_FAIL(name, "write prim"); }
    /* set: {7, 8}, {}, {9}                              (max_def 1) */
    { int32_t v[] = {7, 8, 9}; int16_t d[] = {1, 1, 0, 1}; int16_t r[] = {0, 1, 0, 0};
      if (carquet_writer_write_batch(w, 4, v, 4, d, r) != CARQUET_OK) TEST_FAIL(name, "write set"); }
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    void* buf = NULL; size_t size = 0;
    if (carquet_writer_get_buffer(w, &buf, &size) != CARQUET_OK) TEST_FAIL(name, "get_buffer");

    struct ArrowArray arr; struct ArrowSchema sch; carquet_reader_t* r = NULL;
    carquet_status_t st = read_arrow_buffer(buf, size, &arr, &sch, &r);
    if (st != CARQUET_OK) {
        if (r) carquet_reader_close(r);
        carquet_free(buf);
        carquet_schema_free(schema);
        TEST_FAIL(name, "legacy list shapes could not be assembled");
    }
    int failed = 0;
    const char* why = "";
#define CHECK(cond, msg) do { if (!failed && !(cond)) { failed = 1; why = msg; } } while (0)
    CHECK(arr.length == 3 && arr.n_children == 5, "batch shape");
    if (!failed) {
        /* bare: list<bare: int32 not null> not null */
        const struct ArrowSchema* fs = sch.children[0];
        const struct ArrowArray* fa = arr.children[0];
        CHECK(strcmp(fs->format, "+l") == 0 && !(fs->flags & ARROW_FLAG_NULLABLE), "bare: type");
        CHECK(strcmp(fs->children[0]->format, "i") == 0 &&
              !(fs->children[0]->flags & ARROW_FLAG_NULLABLE) &&
              strcmp(fs->children[0]->name, "bare") == 0, "bare: element");
        const int32_t* o = (const int32_t*)fa->buffers[1];
        const int32_t* v = (const int32_t*)fa->children[0]->buffers[1];
        CHECK(fa->null_count == 0 && o[0] == 0 && o[1] == 2 && o[2] == 2 && o[3] == 3, "bare: offsets");
        CHECK(fa->children[0]->length == 3 && v[0] == 1 && v[1] == 2 && v[2] == 3, "bare: values");
    }
    if (!failed) {
        /* ph: struct<phone: list<phone: struct<number: int64 not null> not null> not null> */
        const struct ArrowSchema* fs = sch.children[1];
        const struct ArrowArray* fa = arr.children[1];
        CHECK(strcmp(fs->format, "+s") == 0 && fa->null_count == 1 && !bit_set(fa->buffers[0], 0), "ph: struct");
        const struct ArrowSchema* ls = fs->children[0];
        const struct ArrowArray* la = fa->children[0];
        CHECK(strcmp(ls->format, "+l") == 0 && !(ls->flags & ARROW_FLAG_NULLABLE), "ph: list type");
        CHECK(strcmp(ls->children[0]->format, "+s") == 0 &&
              !(ls->children[0]->flags & ARROW_FLAG_NULLABLE), "ph: element type");
        const int32_t* o = (const int32_t*)la->buffers[1];
        CHECK(la->length == 3 && o[0] == 0 && o[1] == 0 && o[2] == 0 && o[3] == 2, "ph: offsets");
        const struct ArrowArray* num = la->children[0]->children[0];
        const int64_t* v = (const int64_t*)num->buffers[1];
        CHECK(num->length == 2 && num->null_count == 0 && v[0] == 10 && v[1] == 20, "ph: values");
    }
    if (!failed) {
        /* ll: list<array: list<array: int32 not null> not null> not null */
        const struct ArrowSchema* fs = sch.children[2];
        const struct ArrowArray* fa = arr.children[2];
        CHECK(strcmp(fs->format, "+l") == 0 && strcmp(fs->children[0]->format, "+l") == 0 &&
              strcmp(fs->children[0]->children[0]->format, "i") == 0, "ll: type");
        CHECK(!(fs->flags & ARROW_FLAG_NULLABLE) && !(fs->children[0]->flags & ARROW_FLAG_NULLABLE) &&
              !(fs->children[0]->children[0]->flags & ARROW_FLAG_NULLABLE), "ll: nullability");
        const int32_t* o = (const int32_t*)fa->buffers[1];
        const struct ArrowArray* in = fa->children[0];
        const int32_t* io = (const int32_t*)in->buffers[1];
        const int32_t* v = (const int32_t*)in->children[0]->buffers[1];
        CHECK(o[0] == 0 && o[1] == 2 && o[2] == 2 && o[3] == 4, "ll: outer offsets");
        CHECK(in->length == 4 && io[0] == 0 && io[1] == 2 && io[2] == 3 && io[3] == 3 && io[4] == 4,
              "ll: inner offsets");
        CHECK(in->children[0]->length == 4 && v[0] == 1 && v[1] == 2 && v[2] == 3 && v[3] == 4, "ll: values");
    }
    if (!failed) {
        /* prim: list<item: int32 not null> (nullable) */
        const struct ArrowSchema* fs = sch.children[3];
        const struct ArrowArray* fa = arr.children[3];
        CHECK(strcmp(fs->format, "+l") == 0 && (fs->flags & ARROW_FLAG_NULLABLE) &&
              !(fs->children[0]->flags & ARROW_FLAG_NULLABLE), "prim: type");
        const int32_t* o = (const int32_t*)fa->buffers[1];
        const int32_t* v = (const int32_t*)fa->children[0]->buffers[1];
        CHECK(fa->null_count == 1 && !bit_set(fa->buffers[0], 1) && bit_set(fa->buffers[0], 2), "prim: validity");
        CHECK(o[0] == 0 && o[1] == 2 && o[2] == 2 && o[3] == 2 && v[0] == 5 && v[1] == 6, "prim: values");
    }
    if (!failed) {
        /* set: list<key: int32 not null> not null */
        const struct ArrowSchema* fs = sch.children[4];
        const struct ArrowArray* fa = arr.children[4];
        CHECK(strcmp(fs->format, "+l") == 0 && strcmp(fs->children[0]->format, "i") == 0 &&
              strcmp(fs->children[0]->name, "key") == 0, "set: type");
        const int32_t* o = (const int32_t*)fa->buffers[1];
        const int32_t* v = (const int32_t*)fa->children[0]->buffers[1];
        CHECK(o[0] == 0 && o[1] == 2 && o[2] == 2 && o[3] == 3 && v[0] == 7 && v[1] == 8 && v[2] == 9,
              "set: values");
    }
#undef CHECK
    arr.release(&arr); sch.release(&sch);

    /* The same file through the stream (row-group path). */
    if (!failed) {
        struct ArrowArrayStream stream;
        carquet_arrow_stream_options_t so;
        carquet_arrow_stream_options_init(&so);
        if (carquet_reader_export_arrow_stream(r, &so, &stream, &err) != CARQUET_OK) {
            failed = 1; why = "stream export";
        } else {
            struct ArrowArray b;
            memset(&b, 0, sizeof(b));
            if (stream.get_next(&stream, &b) != 0 || !b.release || b.length != 3) {
                failed = 1; why = "stream batch";
            }
            if (b.release) b.release(&b);
            stream.release(&stream);
        }
    }
    carquet_reader_close(r);
    carquet_free(buf);
    carquet_schema_free(schema);
    if (failed) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

static void noop_release_array(struct ArrowArray* a) { a->release = NULL; }

static int test_ref_export_batch_offset(void) {
    const char* name = "ref_export_batch_offset";
    /* struct<x: int32> batch of 4 rows carrying offset 2 (rows 2..5 of x). */
    static int32_t xv[6] = {0, 1, 2, 3, 4, 5};
    static const void* xbufs[2] = {NULL, xv};
    static const void* sbufs[1] = {NULL};
    struct ArrowArray* x = (struct ArrowArray*)calloc(1, sizeof(*x));
    struct ArrowArray** kids = (struct ArrowArray**)calloc(1, sizeof(*kids));
    x->length = 6; x->n_buffers = 2; x->buffers = xbufs; x->release = noop_release_array;
    kids[0] = x;
    struct ArrowArray top;
    memset(&top, 0, sizeof(top));
    top.length = 4; top.offset = 2; top.n_buffers = 1; top.buffers = sbufs;
    top.n_children = 1; top.children = kids; top.release = noop_release_array;

    carquet_arrow_ref_t* ref = carquet_arrow_ref_create(&top);
    if (!ref) { free(x); free(kids); TEST_FAIL(name, "ref_create"); }
    struct ArrowArray view;
    memset(&view, 0, sizeof(view));
    int failed = 0;
    const char* why = "";
    if (carquet_arrow_ref_export(ref, 0, 4, &view) != CARQUET_OK) {
        failed = 1; why = "export";
    } else {
        /* PyArrow refuses a record batch whose struct has an offset: the view
         * must carry it in the children. */
        if (view.offset != 0 || view.length != 4) { failed = 1; why = "whole view keeps the batch offset"; }
        else if (view.children[0]->offset != 2 || view.children[0]->length != 4) {
            failed = 1; why = "child slice";
        } else {
            const int32_t* v = (const int32_t*)view.children[0]->buffers[1] + view.children[0]->offset;
            if (v[0] != 2 || v[3] != 5) { failed = 1; why = "values"; }
        }
        view.release(&view);
    }
    carquet_arrow_ref_release(ref);
    free(x); free(kids);
    if (failed) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

/* Arrow field metadata blob with one key/value pair. */
static size_t metadata_blob(char* out, const char* key, const char* value) {
    int32_t n = 1, kl = (int32_t)strlen(key), vl = (int32_t)strlen(value);
    char* p = out;
    memcpy(p, &n, 4); p += 4;
    memcpy(p, &kl, 4); p += 4; memcpy(p, key, (size_t)kl); p += kl;
    memcpy(p, &vl, 4); p += 4; memcpy(p, value, (size_t)vl); p += vl;
    return (size_t)(p - out);
}

static void noop_release_schema(struct ArrowSchema* s) { s->release = NULL; }

static int test_arrow_extension_logical(void) {
    const char* name = "arrow_extension_logical";
    char uuid_md[128], json_md[128], other_md[128];
    metadata_blob(uuid_md, "ARROW:extension:name", "arrow.uuid");
    metadata_blob(json_md, "ARROW:extension:name", "arrow.json");
    metadata_blob(other_md, "ARROW:extension:name", "my.ext");

    struct ArrowSchema f[4];
    struct ArrowSchema* kids[4];
    memset(f, 0, sizeof(f));
    const char* names[4] = {"u", "j", "plain16", "other"};
    const char* fmts[4] = {"w:16", "u", "w:16", "w:16"};
    const char* mds[4] = {uuid_md, json_md, NULL, other_md};
    for (int i = 0; i < 4; i++) {
        f[i].format = fmts[i]; f[i].name = names[i]; f[i].metadata = mds[i];
        f[i].flags = ARROW_FLAG_NULLABLE; f[i].release = noop_release_schema;
        kids[i] = &f[i];
    }
    struct ArrowSchema top;
    memset(&top, 0, sizeof(top));
    top.format = "+s"; top.name = ""; top.n_children = 4; top.children = kids;
    top.release = noop_release_schema;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* cs = NULL;
    if (carquet_arrow_import_schema(&top, &cs, &err) != CARQUET_OK || !cs)
        TEST_FAIL(name, "import_schema");
    const parquet_schema_element_t* e = cs->elements;
    int failed = 0;
    const char* why = "";
    if (!(e[1].has_logical_type && e[1].logical_type.id == CARQUET_LOGICAL_UUID &&
          e[1].type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && e[1].type_length == 16)) {
        failed = 1; why = "arrow.uuid is not written as UUID";
    } else if (!(e[2].has_logical_type && e[2].logical_type.id == CARQUET_LOGICAL_JSON)) {
        failed = 1; why = "arrow.json is not written as JSON";
    } else if (e[3].has_logical_type || e[4].has_logical_type) {
        failed = 1; why = "plain fixed_size_binary(16) / unknown extension got a logical type";
    }
    carquet_schema_free(cs);
    if (failed) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

static bool utf8_valid(const uint8_t* s, size_t n) {
    size_t i = 0;
    while (i < n) {
        uint8_t b = s[i];
        size_t len = b < 0x80 ? 1 : (b & 0xE0) == 0xC0 ? 2 : (b & 0xF0) == 0xE0 ? 3 :
                     (b & 0xF8) == 0xF0 ? 4 : 0;
        if (len == 0 || i + len > n) return false;
        for (size_t k = 1; k < len; k++) if ((s[i + k] & 0xC0) != 0x80) return false;
        i += len;
    }
    return true;
}

static int bytes_cmp(const uint8_t* a, size_t na, const uint8_t* b, size_t nb) {
    size_t n = na < nb ? na : nb;
    int c = n ? memcmp(a, b, n) : 0;
    if (c) return c;
    return na < nb ? -1 : (na > nb ? 1 : 0);
}

/* One STRING column (and the same bytes in a plain BYTE_ARRAY column) holding
 * `lo` and `hi`; checks the chunk statistics written for them. */
static int check_string_stats(const char* name, const char* lo, const char* hi,
                              bool expect_max) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_STRING;
    if (carquet_schema_add_column(schema, "s", CARQUET_PHYSICAL_BYTE_ARRAY, &lt,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL(name, "add_column");
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    if (!w) TEST_FAIL(name, "create");
    carquet_byte_array_t v[2] = {{(uint8_t*)lo, (int32_t)strlen(lo)}, {(uint8_t*)hi, (int32_t)strlen(hi)}};
    if (carquet_writer_write_batch(w, 0, v, 2, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, v, 2, NULL, NULL) != CARQUET_OK)
        TEST_FAIL(name, "write");
    if (carquet_writer_close(w) != CARQUET_OK) TEST_FAIL(name, "close");
    void* buf = NULL; size_t size = 0;
    if (carquet_writer_get_buffer(w, &buf, &size) != CARQUET_OK) TEST_FAIL(name, "get_buffer");
    carquet_reader_t* r = carquet_reader_open_buffer(buf, size, NULL, &err);
    if (!r) TEST_FAIL(name, "open");

    int failed = 0;
    const char* why = "";
    const parquet_statistics_t* st = &r->metadata.row_groups[0].columns[0].metadata.statistics;
    const parquet_statistics_t* sb = &r->metadata.row_groups[0].columns[1].metadata.statistics;
    size_t nlo = strlen(lo), nhi = strlen(hi);
    if (!st->min_value || st->min_value_len > 32 ||
        !utf8_valid(st->min_value, (size_t)st->min_value_len)) {
        failed = 1; why = "STRING min is not valid UTF-8";
    } else if (bytes_cmp(st->min_value, (size_t)st->min_value_len, (const uint8_t*)lo, nlo) > 0) {
        failed = 1; why = "STRING min is not a lower bound";
    } else if (expect_max != (st->max_value != NULL)) {
        failed = 1; why = expect_max ? "STRING max missing" : "STRING max should be dropped";
    } else if (st->max_value) {
        if (st->max_value_len > 32 || !utf8_valid(st->max_value, (size_t)st->max_value_len)) {
            failed = 1; why = "STRING max is not valid UTF-8";
        } else if (bytes_cmp(st->max_value, (size_t)st->max_value_len, (const uint8_t*)hi, nhi) < 0) {
            failed = 1; why = "STRING max is not an upper bound";
        }
    }
    /* The untyped column keeps the byte-wise truncation. */
    if (!failed && sb->min_value && nlo > 32 &&
        (sb->min_value_len != 32 || memcmp(sb->min_value, lo, 32) != 0)) {
        failed = 1; why = "BYTE_ARRAY min truncation changed";
    }
    if (!failed && sb->max_value &&
        bytes_cmp(sb->max_value, (size_t)sb->max_value_len, (const uint8_t*)hi, nhi) < 0) {
        failed = 1; why = "BYTE_ARRAY max is not an upper bound";
    }
    carquet_reader_close(r);
    carquet_free(buf);
    carquet_schema_free(schema);
    if (failed) TEST_FAIL(name, why);
    return 0;
}

static int test_utf8_truncated_stats(void) {
    const char* name = "utf8_truncated_stats";
    char two[81], three[61], four[81], top[81], ascii_then[64];
    /* U+00FF (C3 BF) x 40: byte 32 falls on a boundary, the byte increment
     * BF -> C0 does not give UTF-8. */
    for (int i = 0; i < 40; i++) memcpy(two + 2 * i, "\xC3\xBF", 2);
    two[80] = 0;
    /* U+FFFF (EF BF BF) x 20: byte 32 is inside a character, and the last
     * 3-byte code point has no successor of its length. */
    for (int i = 0; i < 20; i++) memcpy(three + 3 * i, "\xEF\xBF\xBF", 3);
    three[60] = 0;
    /* U+1F600 (F0 9F 98 80) x 20 */
    for (int i = 0; i < 20; i++) memcpy(four + 4 * i, "\xF0\x9F\x98\x80", 4);
    four[80] = 0;
    /* U+10FFFF (F4 8F BF BF) x 20: no successor at all, max is dropped */
    for (int i = 0; i < 20; i++) memcpy(top + 4 * i, "\xF4\x8F\xBF\xBF", 4);
    top[80] = 0;
    /* 31 ASCII bytes then 2-byte characters: the cut lands inside the first one */
    memset(ascii_then, 'a', 31);
    for (int i = 0; i < 16; i++) memcpy(ascii_then + 31 + 2 * i, "\xC3\xA9", 2);
    ascii_then[63] = 0;

    if (check_string_stats(name, "a", two, true)) return 1;
    if (check_string_stats(name, three, three, true)) return 1;
    if (check_string_stats(name, four, four, true)) return 1;
    if (check_string_stats(name, ascii_then, ascii_then, true)) return 1;
    if (check_string_stats(name, "a", top, false)) return 1;
    TEST_PASS(name);
    return 0;
}

/* ------------------------------------------------------------------------- */

int main(void) {
    int failures = 0;
    printf("=== v0.8.0 compatibility regression tests ===\n\n");
    failures += test_decimal_scale_zero();
    failures += test_dict_page_no_encodings();
    failures += test_required_leaf_null_struct();
    failures += test_legacy_lists();
    failures += test_ref_export_batch_offset();
    failures += test_arrow_extension_logical();
    failures += test_utf8_truncated_stats();
    printf("\n");
    if (failures == 0) {
        printf("All v0.8.0 compatibility tests passed!\n");
        return 0;
    }
    printf("%d test(s) failed\n", failures);
    return 1;
}
