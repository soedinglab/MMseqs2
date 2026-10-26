/**
 * @file roundtrip_writer.c
 * @brief Comprehensive roundtrip writer for the interop suite.
 *
 * Writes carquet Parquet files that other libraries (PyArrow, DuckDB) read
 * back, and emits a JSON manifest of expected values on stdout for
 * run_interop.py to verify against.
 *
 * Coverage:
 *   - Basic physical types (BOOLEAN/INT32/INT64/FLOAT/DOUBLE/BYTE_ARRAY) with
 *     nullability, across all 5 supported compression codecs.
 *   - The Arrow-broadening encoding paths (opt-in via set_column_encoding):
 *     RLE_DICTIONARY, DELTA_BINARY_PACKED, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT.
 *   - Data Page V2 framing and "ARROW:schema" footer metadata.
 *   - Logical types: STRING, DATE, TIMESTAMP(MICROS, UTC), DECIMAL(9,2),
 *     UUID (FIXED_LEN_BYTE_ARRAY[16]).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <carquet/carquet.h>

#define NUM_ROWS 5000

/* DATE / TIMESTAMP anchors at 2021-01-01T00:00:00Z. */
#define LOGICAL_DATE_BASE 18628                  /* days since 1970-01-01 */
#define LOGICAL_TS_BASE   1609459200000000LL     /* micros since epoch    */

static const char* SAMPLE_STRINGS[10] = {
    "hello", "world", "carquet", "parquet", "test",
    "alpha", "beta", "gamma", "delta", "epsilon"
};

/* ── Basic file: 7 columns, configurable encoding / page version / arrow ───── */

typedef enum {
    ENC_DEFAULT,   /* automatic (PLAIN / auto-BSS)        */
    ENC_DICT,      /* RLE_DICTIONARY on every value column */
    ENC_DELTA,     /* DELTA_BINARY_PACKED + DELTA_BYTE_ARRAY */
    ENC_BSS,       /* BYTE_STREAM_SPLIT where supported    */
} enc_mode_t;

typedef struct {
    carquet_compression_t codec;
    enc_mode_t enc;
    int data_page_version;   /* 1 or 2 */
    bool write_arrow_schema;
    bool write_page_index;   /* emit column/offset index (v0.6.0 offset fix)   */
    int64_t page_size;       /* 0 = default; small value forces multi-page     */
} basic_cfg_t;

/* Column indices for the basic schema. */
enum {
    COL_BOOL = 0, COL_INT32, COL_INT64, COL_FLOAT, COL_DOUBLE,
    COL_STRING, COL_NULLABLE_INT
};

static void generate_basic_data(
    uint8_t* bools, int32_t* int32s, int64_t* int64s,
    float* floats, double* doubles,
    carquet_byte_array_t* strings, int16_t* string_def_levels,
    int32_t* nullable_ints, int16_t* nullable_def_levels,
    int n
) {
    int string_value_count = 0;
    int nullable_int_count = 0;

    for (int i = 0; i < n; i++) {
        bools[i] = (i % 2 == 0) ? 1 : 0;
        int32s[i] = i * 10 - 5000;
        int64s[i] = (int64_t)i * 1000000LL - 2500000000LL;
        floats[i] = (float)i * 0.5f - 1250.0f;
        doubles[i] = (double)i * 0.125 - 312.5;

        if (i % 7 == 0) {
            string_def_levels[i] = 0;          /* every 7th null */
        } else {
            string_def_levels[i] = 1;
            const char* s = SAMPLE_STRINGS[i % 10];
            strings[string_value_count].length = (int32_t)strlen(s);
            strings[string_value_count].data = (uint8_t*)s;
            string_value_count++;
        }

        if (i % 5 == 0) {
            nullable_def_levels[i] = 0;        /* every 5th null */
        } else {
            nullable_def_levels[i] = 1;
            nullable_ints[nullable_int_count] = i * 100;
            nullable_int_count++;
        }
    }
}

static int apply_encoding(carquet_writer_t* w, enc_mode_t enc) {
    if (enc == ENC_DEFAULT) return 0;

    if (enc == ENC_DICT) {
        const int cols[] = {COL_INT32, COL_INT64, COL_FLOAT, COL_DOUBLE,
                            COL_STRING, COL_NULLABLE_INT};
        for (size_t i = 0; i < sizeof(cols) / sizeof(cols[0]); i++)
            if (carquet_writer_set_column_encoding(w, cols[i],
                    CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK)
                return 1;
    } else if (enc == ENC_DELTA) {
        const int dbp[] = {COL_INT32, COL_INT64, COL_NULLABLE_INT};
        for (size_t i = 0; i < sizeof(dbp) / sizeof(dbp[0]); i++)
            if (carquet_writer_set_column_encoding(w, dbp[i],
                    CARQUET_ENCODING_DELTA_BINARY_PACKED) != CARQUET_OK)
                return 1;
        if (carquet_writer_set_column_encoding(w, COL_STRING,
                CARQUET_ENCODING_DELTA_BYTE_ARRAY) != CARQUET_OK)
            return 1;
    } else if (enc == ENC_BSS) {
        /* BYTE_STREAM_SPLIT for INT32/INT64 is a newer Parquet addition that
         * DuckDB 1.4.3 cannot read; keep this variant portable (FLOAT/DOUBLE). */
        const int bss[] = {COL_FLOAT, COL_DOUBLE};
        for (size_t i = 0; i < sizeof(bss) / sizeof(bss[0]); i++)
            if (carquet_writer_set_column_encoding(w, bss[i],
                    CARQUET_ENCODING_BYTE_STREAM_SPLIT) != CARQUET_OK)
                return 1;
    }
    return 0;
}

/* Write one full copy of the basic dataset (NUM_ROWS rows) to a writer. */
static carquet_status_t populate_basic(carquet_writer_t* writer) {
    uint8_t* bools = malloc(NUM_ROWS * sizeof(uint8_t));
    int32_t* int32s = malloc(NUM_ROWS * sizeof(int32_t));
    int64_t* int64s = malloc(NUM_ROWS * sizeof(int64_t));
    float* floats = malloc(NUM_ROWS * sizeof(float));
    double* doubles = malloc(NUM_ROWS * sizeof(double));
    carquet_byte_array_t* strings = malloc(NUM_ROWS * sizeof(carquet_byte_array_t));
    int16_t* string_def_levels = malloc(NUM_ROWS * sizeof(int16_t));
    int32_t* nullable_ints = malloc(NUM_ROWS * sizeof(int32_t));
    int16_t* nullable_def_levels = malloc(NUM_ROWS * sizeof(int16_t));

    generate_basic_data(bools, int32s, int64s, floats, doubles,
                        strings, string_def_levels,
                        nullable_ints, nullable_def_levels, NUM_ROWS);

    carquet_status_t st = CARQUET_OK;
    #define WB(col, vals, defs) do { \
        carquet_status_t _s = carquet_writer_write_batch(writer, (col), (vals), \
            NUM_ROWS, (defs), NULL); \
        if (_s != CARQUET_OK) st = _s; \
    } while (0)
    WB(COL_BOOL, bools, NULL);
    WB(COL_INT32, int32s, NULL);
    WB(COL_INT64, int64s, NULL);
    WB(COL_FLOAT, floats, NULL);
    WB(COL_DOUBLE, doubles, NULL);
    WB(COL_STRING, strings, string_def_levels);
    WB(COL_NULLABLE_INT, nullable_ints, nullable_def_levels);
    #undef WB

    free(bools);
    free(int32s);
    free(int64s);
    free(floats);
    free(doubles);
    free(strings);
    free(string_def_levels);
    free(nullable_ints);
    free(nullable_def_levels);
    return st;
}

/* Build the shared 7-column basic schema (also used by the append path). */
static carquet_schema_t* make_basic_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return NULL;

    (void)carquet_schema_add_column(schema, "bool_col", CARQUET_PHYSICAL_BOOLEAN,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "int32_col", CARQUET_PHYSICAL_INT32,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "int64_col", CARQUET_PHYSICAL_INT64,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "float_col", CARQUET_PHYSICAL_FLOAT,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "double_col", CARQUET_PHYSICAL_DOUBLE,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "string_col", CARQUET_PHYSICAL_BYTE_ARRAY,
                                    NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    (void)carquet_schema_add_column(schema, "nullable_int", CARQUET_PHYSICAL_INT32,
                                    NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    return schema;
}

static int write_basic_file(const char* path, basic_cfg_t cfg) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = make_basic_schema();
    if (!schema) return 1;

    /* When emitting ARROW:schema, attach per-field custom_metadata (variable
     * labels) so pyarrow reads them back as field metadata. Elements 1..2 are
     * bool_col / int32_col (root is element 0). */
    if (cfg.write_arrow_schema) {
        (void)carquet_schema_set_field_metadata(schema, 1, "Label",
                                                "Boolean Flag Variable");
        (void)carquet_schema_set_field_metadata(schema, 2, "Label",
                                                "Signed 32-bit Measure");
        (void)carquet_schema_set_field_metadata(schema, 2, "Units", "count");
    }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = cfg.codec;
    opts.row_group_size = 2000 * 100;  /* force multiple row groups */
    opts.data_page_version = cfg.data_page_version;
    opts.write_arrow_schema = cfg.write_arrow_schema;
    opts.write_page_index = cfg.write_page_index;
    if (cfg.page_size > 0) opts.page_size = cfg.page_size;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "Failed to create writer: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }

    if (apply_encoding(writer, cfg.enc) != 0) {
        fprintf(stderr, "Failed to set column encoding for %s\n", path);
        carquet_writer_close(writer);
        carquet_schema_free(schema);
        return 1;
    }

    carquet_status_t status = populate_basic(writer);
    if (status == CARQUET_OK) status = carquet_writer_close(writer);
    else (void)carquet_writer_close(writer);

    carquet_schema_free(schema);

    return (status == CARQUET_OK) ? 0 : 1;
}

/* ── Append file: create, then open_append a second copy of the dataset ─────── */

static int write_append_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Phase 1: create the file with the first copy. */
    carquet_schema_t* schema = make_basic_schema();
    if (!schema) return 1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.row_group_size = 2000 * 100;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "append: create failed: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }
    carquet_status_t status = populate_basic(writer);
    if (status == CARQUET_OK) status = carquet_writer_close(writer);
    else (void)carquet_writer_close(writer);
    if (status != CARQUET_OK) {
        fprintf(stderr, "append: first write failed\n");
        carquet_schema_free(schema);
        return 1;
    }

    /* Phase 2: reopen for append and add a second copy as new row group(s). */
    carquet_writer_t* aw = carquet_writer_open_append(path, schema, &opts, &err);
    if (!aw) {
        fprintf(stderr, "append: open_append failed: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }
    status = populate_basic(aw);
    if (status == CARQUET_OK) status = carquet_writer_close(aw);
    else (void)carquet_writer_close(aw);

    carquet_schema_free(schema);
    return (status == CARQUET_OK) ? 0 : 1;
}

/* ── Nested file: LIST + STRUCT (locks in carquet→Arrow nested write parity) ──
 *
 * Schema:
 *   id    (required INT32)
 *   tags  (optional LIST<INT32>)
 *   info  (required group / struct) { name (required STRING), age (required INT32) }
 *
 * Expected logical content (4 rows), verified by run_interop.py:
 *   id   = [1, 2, 3, 4]
 *   tags = [[100, 200], null, [300], [400, 500, 600]]
 *   info = [{a,10}, {b,20}, {c,30}, {d,40}]
 */
static int write_nested_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;

    carquet_logical_type_t str_lt = { .id = CARQUET_LOGICAL_STRING };

    /* leaf 0: id */
    (void)carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    /* leaf 1: tags list element */
    (void)carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_INT32,
                                  NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    /* struct group with two required leaves (leaf 2: info.name, leaf 3: info.age) */
    int32_t info = carquet_schema_add_group(schema, "info",
                                            CARQUET_REPETITION_REQUIRED, 0);
    (void)carquet_schema_add_column(schema, "name", CARQUET_PHYSICAL_BYTE_ARRAY,
                                    &str_lt, CARQUET_REPETITION_REQUIRED, 0, info);
    (void)carquet_schema_add_column(schema, "age", CARQUET_PHYSICAL_INT32,
                                    NULL, CARQUET_REPETITION_REQUIRED, 0, info);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    /* Emit the nested Arrow schema (LIST<int32> + STRUCT) into the footer so a
     * reader recovers the exact Arrow type tree, not just the Parquet schema. */
    opts.write_arrow_schema = true;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) {
        fprintf(stderr, "nested: create failed: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }

    carquet_status_t st = CARQUET_OK;

    int32_t ids[] = {1, 2, 3, 4};
    if (carquet_writer_write_batch(w, 0, ids, 4, NULL, NULL) != CARQUET_OK) st = 1;

    /* tags: [[100,200], null, [300], [400,500,600]] */
    int32_t tag_vals[] = {100, 200, 300, 400, 500, 600};
    int16_t tag_def[]  = {3, 3, 0, 3, 3, 3, 3};
    int16_t tag_rep[]  = {0, 1, 0, 0, 0, 1, 1};
    if (carquet_writer_write_batch(w, 1, tag_vals, 7, tag_def, tag_rep) != CARQUET_OK) st = 1;

    const char* names[] = {"a", "b", "c", "d"};
    carquet_byte_array_t name_vals[4];
    for (int i = 0; i < 4; i++) {
        name_vals[i].data = (uint8_t*)names[i];
        name_vals[i].length = 1;
    }
    if (carquet_writer_write_batch(w, 2, name_vals, 4, NULL, NULL) != CARQUET_OK) st = 1;

    int32_t ages[] = {10, 20, 30, 40};
    if (carquet_writer_write_batch(w, 3, ages, 4, NULL, NULL) != CARQUET_OK) st = 1;

    if (st == CARQUET_OK) st = carquet_writer_close(w);
    else (void)carquet_writer_close(w);

    carquet_schema_free(schema);
    return (st == CARQUET_OK) ? 0 : 1;
}

static void print_basic_columns(void) {
    printf("      \"columns\": {\n");
    printf("        \"bool_col\": { \"first\": [true, false, true, false, true], \"type\": \"bool\" },\n");
    printf("        \"int32_col\": { \"first\": [-5000, -4990, -4980, -4970, -4960], \"type\": \"int32\" },\n");
    printf("        \"int64_col\": { \"first\": [-2500000000, -2499000000, -2498000000, -2497000000, -2496000000], \"type\": \"int64\" },\n");
    printf("        \"float_col\": { \"first\": [-1250.0, -1249.5, -1249.0, -1248.5, -1248.0], \"type\": \"float\" },\n");
    printf("        \"double_col\": { \"first\": [-312.5, -312.375, -312.25, -312.125, -312.0], \"type\": \"double\" },\n");
    printf("        \"string_col\": { \"first\": [null, \"world\", \"carquet\", \"parquet\", \"test\"], \"null_pattern\": \"every_7th\", \"type\": \"string\" },\n");
    printf("        \"nullable_int\": { \"first\": [null, 100, 200, 300, 400], \"null_pattern\": \"every_5th\", \"type\": \"int32\" }\n");
    printf("      }\n");
}

/* ── Logical-types file ────────────────────────────────────────────────────── */

enum {
    LCOL_STRING = 0, LCOL_DATE, LCOL_TS, LCOL_DECIMAL, LCOL_UUID
};

static int write_logical_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return 1;

    carquet_logical_type_t str_lt = { .id = CARQUET_LOGICAL_STRING };
    carquet_logical_type_t date_lt = { .id = CARQUET_LOGICAL_DATE };
    carquet_logical_type_t ts_lt = { .id = CARQUET_LOGICAL_TIMESTAMP };
    ts_lt.params.timestamp.unit = CARQUET_TIME_UNIT_MICROS;
    ts_lt.params.timestamp.is_adjusted_to_utc = true;
    carquet_logical_type_t dec_lt = { .id = CARQUET_LOGICAL_DECIMAL };
    dec_lt.params.decimal.precision = 9;
    dec_lt.params.decimal.scale = 2;
    carquet_logical_type_t uuid_lt = { .id = CARQUET_LOGICAL_UUID };

    (void)carquet_schema_add_column(schema, "str_col", CARQUET_PHYSICAL_BYTE_ARRAY,
                                    &str_lt, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "date_col", CARQUET_PHYSICAL_INT32,
                                    &date_lt, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "ts_col", CARQUET_PHYSICAL_INT64,
                                    &ts_lt, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "dec_col", CARQUET_PHYSICAL_INT32,
                                    &dec_lt, CARQUET_REPETITION_REQUIRED, 0, 0);
    (void)carquet_schema_add_column(schema, "uuid_col", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
                                    &uuid_lt, CARQUET_REPETITION_REQUIRED, 16, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.row_group_size = 2000 * 100;

    carquet_writer_t* writer = carquet_writer_create(path, schema, &opts, &err);
    if (!writer) {
        fprintf(stderr, "Failed to create logical writer: %s\n", err.message);
        carquet_schema_free(schema);
        return 1;
    }

    carquet_byte_array_t* strings = malloc(NUM_ROWS * sizeof(carquet_byte_array_t));
    int32_t* dates = malloc(NUM_ROWS * sizeof(int32_t));
    int64_t* tss = malloc(NUM_ROWS * sizeof(int64_t));
    int32_t* decs = malloc(NUM_ROWS * sizeof(int32_t));
    uint8_t* uuids = malloc((size_t)NUM_ROWS * 16);

    for (int i = 0; i < NUM_ROWS; i++) {
        const char* s = SAMPLE_STRINGS[i % 10];
        strings[i].data = (uint8_t*)s;
        strings[i].length = (int32_t)strlen(s);
        dates[i] = LOGICAL_DATE_BASE + i;
        tss[i] = LOGICAL_TS_BASE + (int64_t)i * 1000000LL;
        decs[i] = i % 100000;   /* unscaled; value = decs[i] / 100 */
        for (int j = 0; j < 16; j++)
            uuids[(size_t)i * 16 + j] = (uint8_t)((i + j) & 0xFF);
    }

    (void)carquet_writer_write_batch(writer, LCOL_STRING, strings, NUM_ROWS, NULL, NULL);
    (void)carquet_writer_write_batch(writer, LCOL_DATE, dates, NUM_ROWS, NULL, NULL);
    (void)carquet_writer_write_batch(writer, LCOL_TS, tss, NUM_ROWS, NULL, NULL);
    (void)carquet_writer_write_batch(writer, LCOL_DECIMAL, decs, NUM_ROWS, NULL, NULL);
    (void)carquet_writer_write_batch(writer, LCOL_UUID, uuids, NUM_ROWS, NULL, NULL);

    carquet_status_t status = carquet_writer_close(writer);

    free(strings);
    free(dates);
    free(tss);
    free(decs);
    free(uuids);
    carquet_schema_free(schema);

    return (status == CARQUET_OK) ? 0 : 1;
}

static void print_logical_columns(void) {
    printf("      \"columns\": {\n");
    printf("        \"str_col\": { \"first\": [\"hello\", \"world\", \"carquet\", \"parquet\", \"test\"], \"type\": \"string\" },\n");

    printf("        \"date_col\": { \"first\": [");
    for (int i = 0; i < 5; i++)
        printf("%s%d", i ? ", " : "", LOGICAL_DATE_BASE + i);
    printf("], \"type\": \"date\" },\n");

    printf("        \"ts_col\": { \"first\": [");
    for (int i = 0; i < 5; i++)
        printf("%s%lld", i ? ", " : "", (long long)(LOGICAL_TS_BASE + (int64_t)i * 1000000LL));
    printf("], \"type\": \"timestamp_us\" },\n");

    printf("        \"dec_col\": { \"first\": [");
    for (int i = 0; i < 5; i++)
        printf("%s\"%d.%02d\"", i ? ", " : "", i / 100, i % 100);
    printf("], \"type\": \"decimal\" },\n");

    printf("        \"uuid_col\": { \"first\": [");
    for (int i = 0; i < 5; i++) {
        printf("%s\"", i ? ", " : "");
        for (int j = 0; j < 16; j++)
            printf("%02x", (unsigned)((i + j) & 0xFF));
        printf("\"");
    }
    printf("], \"type\": \"uuid\" }\n");
    printf("      }\n");
}

/* ── Manifest ──────────────────────────────────────────────────────────────── */

/* Build a filesystem-safe name from a display tag (e.g. "zstd+dict"). */
static void safe_filename(char* out, size_t out_sz, const char* dir, const char* tag) {
    char clean[64];
    size_t k = 0;
    for (size_t i = 0; tag[i] && k < sizeof(clean) - 1; i++)
        clean[k++] = (tag[i] == '+' || tag[i] == ' ') ? '_' : tag[i];
    clean[k] = '\0';
    snprintf(out, out_sz, "%s/carquet_%s.parquet", dir, clean);
}

/* ---- Deep-nested file via the generic Arrow shredder (write_arrow) ----
 * matrix: list<list<int32>> = [[[1,2],[3]], [], [[4]]]  (verifies arbitrary
 * depth on the wire; read back by PyArrow / DuckDB in run_interop.py). */
static void nn_rel_schema(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format); free((void*)s->name);
    for (int64_t i = 0; i < s->n_children; i++)
        if (s->children[i]) { if (s->children[i]->release) s->children[i]->release(s->children[i]); free(s->children[i]); }
    free(s->children); s->release = NULL;
}
static void nn_rel_array(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) { for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]); free(a->buffers); }
    for (int64_t i = 0; i < a->n_children; i++)
        if (a->children[i]) { if (a->children[i]->release) a->children[i]->release(a->children[i]); free(a->children[i]); }
    free(a->children); a->release = NULL;
}
static char* nn_dup(const char* s){ char* o=(char*)malloc(strlen(s)+1); strcpy(o,s); return o; }
static struct ArrowSchema* nn_S(const char* fmt, const char* name, int nullable, int nch){
    struct ArrowSchema* s=(struct ArrowSchema*)calloc(1,sizeof(*s));
    s->format=nn_dup(fmt); s->name=nn_dup(name); s->flags=nullable?ARROW_FLAG_NULLABLE:0;
    s->n_children=nch; s->children=nch?(struct ArrowSchema**)calloc(nch,sizeof(void*)):NULL;
    s->release=nn_rel_schema; return s;
}
static struct ArrowArray* nn_A(int64_t len, int nbuf, int nch){
    struct ArrowArray* a=(struct ArrowArray*)calloc(1,sizeof(*a));
    a->length=len; a->null_count=-1; a->n_buffers=nbuf; a->n_children=nch;
    a->buffers=nbuf?(const void**)calloc(nbuf,sizeof(void*)):NULL;
    a->children=nch?(struct ArrowArray**)calloc(nch,sizeof(void*)):NULL;
    a->release=nn_rel_array; return a;
}
static int32_t* nn_i32(const int32_t* x, int64_t n){ int32_t* o=(int32_t*)malloc((n?n:1)*4); memcpy(o,x,n*4); return o; }

static struct ArrowSchema* nn_build_schema(void) {
    struct ArrowSchema* t=nn_S("+s","schema",0,1);
    t->children[0]=nn_S("+l","matrix",1,1);
    t->children[0]->children[0]=nn_S("+l","element",1,1);
    t->children[0]->children[0]->children[0]=nn_S("i","element",0,0);
    return t;
}

static int write_deep_nested_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* cs = NULL;
    if (carquet_arrow_import_schema(nn_build_schema(), &cs, &err) != CARQUET_OK) return 1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.write_arrow_schema = true;
    carquet_writer_t* w = carquet_writer_create(path, cs, &opts, &err);
    if (!w) { carquet_schema_free(cs); return 1; }

    int32_t oo[4]={0,2,2,3}, io[4]={0,2,3,4}, iv[4]={1,2,3,4};
    struct ArrowArray* ii=nn_A(4,2,0); ii->buffers[1]=nn_i32(iv,4);
    struct ArrowArray* il=nn_A(3,2,1); il->buffers[1]=nn_i32(io,4); il->children[0]=ii;
    struct ArrowArray* ol=nn_A(3,2,1); ol->buffers[1]=nn_i32(oo,4); ol->children[0]=il;
    struct ArrowArray* top=nn_A(3,1,1); top->children[0]=ol;

    int rc = carquet_writer_write_arrow(w, top, nn_build_schema(), &err) == CARQUET_OK ? 0 : 1;
    if (carquet_writer_close(w) != CARQUET_OK) rc = 1;
    carquet_schema_free(cs);
    return rc;
}

int main(int argc, char** argv) {
    const char* output_dir = (argc > 1) ? argv[1] : "/tmp";

    if (carquet_init() != CARQUET_OK) {
        fprintf(stderr, "Failed to init carquet\n");
        return 1;
    }

    struct {
        const char* tag;
        basic_cfg_t cfg;
    } basics[] = {
        {"uncompressed", {CARQUET_COMPRESSION_UNCOMPRESSED, ENC_DEFAULT, 1, false, false, 0}},
        {"snappy",       {CARQUET_COMPRESSION_SNAPPY,       ENC_DEFAULT, 1, false, false, 0}},
        {"gzip",         {CARQUET_COMPRESSION_GZIP,         ENC_DEFAULT, 1, false, false, 0}},
        {"lz4_raw",      {CARQUET_COMPRESSION_LZ4_RAW,      ENC_DEFAULT, 1, false, false, 0}},
        /* v0.7.0: codec 5 (deprecated Hadoop-framed LZ4), distinct from LZ4_RAW. */
        {"lz4",          {CARQUET_COMPRESSION_LZ4,          ENC_DEFAULT, 1, false, false, 0}},
        {"zstd",         {CARQUET_COMPRESSION_ZSTD,         ENC_DEFAULT, 1, false, false, 0}},
        /* v0.7.2: built-in BROTLI (codec 4). */
        {"brotli",       {CARQUET_COMPRESSION_BROTLI,       ENC_DEFAULT, 1, false, false, 0}},
        {"zstd+dict",    {CARQUET_COMPRESSION_ZSTD,         ENC_DICT,    1, false, false, 0}},
        {"zstd+delta",   {CARQUET_COMPRESSION_ZSTD,         ENC_DELTA,   1, false, false, 0}},
        {"zstd+bss",     {CARQUET_COMPRESSION_ZSTD,         ENC_BSS,     1, false, false, 0}},
        {"zstd+pagev2",  {CARQUET_COMPRESSION_ZSTD,         ENC_DEFAULT, 2, false, false, 0}},
        {"zstd+arrow",   {CARQUET_COMPRESSION_ZSTD,         ENC_DEFAULT, 1, true,  false, 0}},
        /* v0.6.0: page index (exercises the absolute-offset fix). */
        {"zstd+pageidx", {CARQUET_COMPRESSION_ZSTD,         ENC_DEFAULT, 1, false, true, 0}},
        /* v0.6.0: tiny pages force OPTIONAL columns across multiple page chunks,
         * exercising the nullable batched-write value mis-indexing fix. */
        {"zstd+multipage", {CARQUET_COMPRESSION_ZSTD,       ENC_DEFAULT, 1, false, false, 4096}},
    };
    const int num_basics = (int)(sizeof(basics) / sizeof(basics[0]));

    printf("{\n");
    printf("  \"num_rows\": %d,\n", NUM_ROWS);
    printf("  \"files\": [\n");

    int first = 1;

    for (int i = 0; i < num_basics; i++) {
        char path[512];
        safe_filename(path, sizeof(path), output_dir, basics[i].tag);

        if (write_basic_file(path, basics[i].cfg) != 0) {
            fprintf(stderr, "Failed to write %s\n", path);
            continue;
        }

        if (!first) printf(",\n");
        first = 0;

        printf("    {\n");
        printf("      \"path\": \"%s\",\n", path);
        printf("      \"compression\": \"%s\",\n", basics[i].tag);
        if (basics[i].cfg.write_arrow_schema) {
            /* Expected Arrow per-field custom_metadata (see write_basic_file). */
            printf("      \"field_metadata\": {\n");
            printf("        \"bool_col\": {\"Label\": \"Boolean Flag Variable\"},\n");
            printf("        \"int32_col\": {\"Label\": \"Signed 32-bit Measure\", \"Units\": \"count\"}\n");
            printf("      },\n");
        }
        print_basic_columns();
        printf("    }");
    }

    /* Logical-types file. */
    {
        char path[512];
        safe_filename(path, sizeof(path), output_dir, "logical");
        if (write_logical_file(path) == 0) {
            if (!first) printf(",\n");
            first = 0;
            printf("    {\n");
            printf("      \"path\": \"%s\",\n", path);
            printf("      \"compression\": \"%s\",\n", "logical");
            print_logical_columns();
            printf("    }");
        } else {
            fprintf(stderr, "Failed to write %s\n", path);
        }
    }

    /* Append file: two copies of the dataset (v0.6.0 open_append). Carries its
     * own num_rows / verification since the totals differ from the basics. */
    {
        char path[512];
        safe_filename(path, sizeof(path), output_dir, "append");
        if (write_append_file(path) == 0) {
            const long long basic_sum =
                (long long)((NUM_ROWS - 1) * NUM_ROWS / 2 * 10 - 5000LL * NUM_ROWS);
            if (!first) printf(",\n");
            first = 0;
            printf("    {\n");
            printf("      \"path\": \"%s\",\n", path);
            printf("      \"compression\": \"%s\",\n", "append");
            printf("      \"num_rows\": %d,\n", NUM_ROWS * 2);
            printf("      \"verification\": {\n");
            printf("        \"null_count_string_col\": %d,\n", 2 * ((NUM_ROWS + 6) / 7));
            printf("        \"null_count_nullable_int\": %d,\n", 2 * ((NUM_ROWS + 4) / 5));
            printf("        \"int32_sum\": %lld,\n", 2 * basic_sum);
            printf("        \"last_int32\": %d\n", (NUM_ROWS - 1) * 10 - 5000);
            printf("      },\n");
            print_basic_columns();
            printf("    }");
        } else {
            fprintf(stderr, "Failed to write %s\n", path);
        }
    }

    /* Nested file: LIST + STRUCT. Verified specially by run_interop.py against
     * a fixed expected structure (the column model is nested, not flat). */
    {
        char path[512];
        safe_filename(path, sizeof(path), output_dir, "nested");
        if (write_nested_file(path) == 0) {
            if (!first) printf(",\n");
            first = 0;
            printf("    {\n");
            printf("      \"path\": \"%s\",\n", path);
            printf("      \"compression\": \"%s\",\n", "nested");
            printf("      \"nested\": true,\n");
            printf("      \"num_rows\": 4,\n");
            printf("      \"columns\": {}\n");
            printf("    }");
        } else {
            fprintf(stderr, "Failed to write %s\n", path);
        }
    }

    /* Deep-nested file (list<list<int32>>) written via the generic shredder. */
    {
        char path[512];
        safe_filename(path, sizeof(path), output_dir, "nested_deep");
        if (write_deep_nested_file(path) == 0) {
            if (!first) printf(",\n");
            first = 0;
            printf("    {\n");
            printf("      \"path\": \"%s\",\n", path);
            printf("      \"compression\": \"%s\",\n", "nested_deep");
            printf("      \"nested\": true,\n");
            printf("      \"num_rows\": 3,\n");
            printf("      \"columns\": {}\n");
            printf("    }");
        } else {
            fprintf(stderr, "Failed to write %s\n", path);
        }
    }

    printf("\n  ],\n");
    printf("  \"verification\": {\n");
    printf("    \"row_counts\": %d,\n", NUM_ROWS);
    printf("    \"null_count_string_col\": %d,\n", (NUM_ROWS + 6) / 7);
    printf("    \"null_count_nullable_int\": %d,\n", (NUM_ROWS + 4) / 5);
    printf("    \"bool_true_count\": %d,\n", (NUM_ROWS + 1) / 2);
    printf("    \"int32_sum\": %lld,\n", (long long)((NUM_ROWS - 1) * NUM_ROWS / 2 * 10 - 5000LL * NUM_ROWS));
    printf("    \"last_int32\": %d\n", (NUM_ROWS - 1) * 10 - 5000);
    printf("  }\n");
    printf("}\n");

    return 0;
}
