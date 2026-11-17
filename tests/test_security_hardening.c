/**
 * @file test_security_hardening.c
 * @brief Regression tests for memory-safety fixes on untrusted-input and
 *        error paths (cross-platform + sanitizer audit, Oct 2026).
 *
 * One test per fix:
 *  - geospatial_types list repeated in one ColumnMetaData appended past
 *    types[64] (and over num_types itself)
 *  - schema validation at open: num_children not covering the element list,
 *    FLBA type_length <= 0, FLOAT16 on a non-FLBA(2) leaf
 *  - DELTA_BYTE_ARRAY prefix lengths not bounded by the previous value
 *    (attacker-chosen reconstruction size) and DELTA_LENGTH_BYTE_ARRAY sum
 *  - batch-reader coalesced split planner trusting page num_values (heap
 *    write past the slot buffer on mmap + multi-threaded reads)
 *  - pipelined writer: a background write error left current_row_group
 *    aliased by the pipeline slot (double free / use-after-free on close)
 *  - Arrow import: struct child-count mismatch, leaf/writer column type
 *    mismatch, NULL data buffer with non-empty values
 */

#include <carquet/carquet.h>
#include "thrift/thrift_encode.h"
#include "thrift/thrift_decode.h"
#include "thrift/parquet_types.h"
#include "core/buffer.h"
#include "core/arena.h"
#include "test_helpers.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if !defined(_WIN32)
#include <signal.h>
#include <sys/resource.h>
#endif

extern carquet_status_t carquet_delta_encode_int32(const int32_t* values, int32_t count,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);
extern carquet_status_t carquet_delta_strings_decode_lengths(const uint8_t* data,
    size_t data_size, int32_t num_values, int32_t** prefix_out, int32_t** suffix_out,
    size_t* header_bytes, size_t* total_size);
extern carquet_status_t carquet_delta_length_decode(const uint8_t* data, size_t data_size,
    carquet_byte_array_t* values, int32_t num_values, size_t* bytes_consumed);
extern void carquet_mem_free(void* ptr);

/* ------------------------------------------------------------------------- */
/* Thrift footer crafting helpers                                             */
/* ------------------------------------------------------------------------- */

static void w_i32(thrift_encoder_t* e, int16_t id, int32_t v) {
    thrift_write_field_header(e, THRIFT_TYPE_I32, id);
    thrift_write_i32(e, v);
}
static void w_i64(thrift_encoder_t* e, int16_t id, int64_t v) {
    thrift_write_field_header(e, THRIFT_TYPE_I64, id);
    thrift_write_i64(e, v);
}
static void w_str(thrift_encoder_t* e, int16_t id, const char* s) {
    thrift_write_field_header(e, THRIFT_TYPE_BINARY, id);
    thrift_write_string(e, s);
}

typedef struct {
    const char* name;
    int32_t type;          /* -1: group (no type field) */
    int32_t type_length;   /* 0: absent */
    int32_t num_children;  /* -1: absent */
    int float16;           /* emit LogicalType FLOAT16 */
} elem_spec_t;

/* Note: thrift_write_struct_end() emits the STOP byte itself. */
static void w_schema_element(thrift_encoder_t* e, const elem_spec_t* s) {
    thrift_write_struct_begin(e);
    if (s->type >= 0) w_i32(e, 1, s->type);
    if (s->type_length) w_i32(e, 2, s->type_length);
    if (s->type >= 0) w_i32(e, 3, 0);  /* REQUIRED */
    w_str(e, 4, s->name);
    if (s->num_children >= 0) w_i32(e, 5, s->num_children);
    if (s->float16) {
        thrift_write_field_header(e, THRIFT_TYPE_STRUCT, 10);   /* logicalType */
        thrift_write_struct_begin(e);
        thrift_write_field_header(e, THRIFT_TYPE_STRUCT, 15);   /* FLOAT16 */
        thrift_write_struct_begin(e);
        thrift_write_struct_end(e);
        thrift_write_struct_end(e);
    }
    thrift_write_struct_end(e);
}

/* FileMetaData with the given schema and no row groups. */
static void w_file_metadata(thrift_encoder_t* e, const elem_spec_t* elems, int n) {
    thrift_write_struct_begin(e);
    w_i32(e, 1, 1);
    thrift_write_field_header(e, THRIFT_TYPE_LIST, 2);
    thrift_write_list_begin(e, THRIFT_TYPE_STRUCT, n);
    for (int i = 0; i < n; i++) w_schema_element(e, &elems[i]);
    w_i64(e, 3, 0);
    thrift_write_field_header(e, THRIFT_TYPE_LIST, 4);
    thrift_write_list_begin(e, THRIFT_TYPE_STRUCT, 0);
    thrift_write_struct_end(e);
}

/* Wrap a serialized footer into a minimal Parquet file image. */
static uint8_t* make_file(const carquet_buffer_t* footer, size_t* out_size) {
    size_t fl = carquet_buffer_size(footer);
    size_t total = 4 + fl + 4 + 4;
    uint8_t* f = (uint8_t*)malloc(total);
    memcpy(f, "PAR1", 4);
    memcpy(f + 4, carquet_buffer_data_const(footer), fl);
    uint32_t len = (uint32_t)fl;
    f[4 + fl + 0] = (uint8_t)(len);
    f[4 + fl + 1] = (uint8_t)(len >> 8);
    f[4 + fl + 2] = (uint8_t)(len >> 16);
    f[4 + fl + 3] = (uint8_t)(len >> 24);
    memcpy(f + 4 + fl + 4, "PAR1", 4);
    *out_size = total;
    return f;
}

/* Returns 1 if a file with this schema opens, 0 if it is rejected. */
static int schema_opens(const elem_spec_t* elems, int n) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);
    w_file_metadata(&enc, elems, n);
    size_t size = 0;
    uint8_t* file = make_file(&buf, &size);
    carquet_buffer_destroy(&buf);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(file, size, NULL, &err);
    int ok = r != NULL;
    if (r) carquet_reader_close(r);
    free(file);
    return ok;
}

/* ------------------------------------------------------------------------- */

static int test_geo_types_repeated_field(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    thrift_encoder_t e;
    thrift_encoder_init(&e, &buf);

    elem_spec_t schema[2] = {
        { "schema", -1, 0, 1, 0 },
        { "g", CARQUET_PHYSICAL_BYTE_ARRAY, 0, -1, 0 },
    };

    thrift_write_struct_begin(&e);
    w_i32(&e, 1, 1);
    thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
    thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 2);
    for (int i = 0; i < 2; i++) w_schema_element(&e, &schema[i]);
    w_i64(&e, 3, 0);
    thrift_write_field_header(&e, THRIFT_TYPE_LIST, 4);      /* row_groups */
    thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 1);
      thrift_write_struct_begin(&e);                        /* RowGroup */
      thrift_write_field_header(&e, THRIFT_TYPE_LIST, 1);    /* columns */
      thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 1);
        thrift_write_struct_begin(&e);                      /* ColumnChunk */
        w_i64(&e, 2, 4);
        thrift_write_field_header(&e, THRIFT_TYPE_STRUCT, 3); /* ColumnMetaData */
        thrift_write_struct_begin(&e);
          w_i32(&e, 1, CARQUET_PHYSICAL_BYTE_ARRAY);
          thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
          thrift_write_list_begin(&e, THRIFT_TYPE_I32, 1);
          thrift_write_i32(&e, 0);
          thrift_write_field_header(&e, THRIFT_TYPE_LIST, 3);
          thrift_write_list_begin(&e, THRIFT_TYPE_BINARY, 1);
          thrift_write_string(&e, "g");
          w_i32(&e, 4, 0);
          w_i64(&e, 5, 0);
          w_i64(&e, 6, 0);
          w_i64(&e, 7, 0);
          w_i64(&e, 9, 4);
          thrift_write_field_header(&e, THRIFT_TYPE_STRUCT, 17); /* geo stats */
          thrift_write_struct_begin(&e);
            /* geospatial_types twice: 64 codes, then 2 more. */
            thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
            thrift_write_list_begin(&e, THRIFT_TYPE_I32, CARQUET_GEO_MAX_TYPES);
            for (int i = 0; i < CARQUET_GEO_MAX_TYPES; i++) thrift_write_i32(&e, 1000 + i);
            /* Same field id again (delta 0 -> long-form header). */
            thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
            thrift_write_list_begin(&e, THRIFT_TYPE_I32, 2);
            thrift_write_i32(&e, 5);
            thrift_write_i32(&e, 7);
          thrift_write_struct_end(&e);
        thrift_write_struct_end(&e);
        thrift_write_struct_end(&e);
      w_i64(&e, 2, 0);
      w_i64(&e, 3, 0);
      thrift_write_struct_end(&e);
    thrift_write_struct_end(&e);

    if (thrift_encoder_has_error(&e)) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("geo_types_repeated_field", "encoder error");
    }

    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("geo_types_repeated_field", "arena init");
    }
    parquet_file_metadata_t md;
    memset(&md, 0, sizeof(md));
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_status_t st = parquet_parse_file_metadata(
        carquet_buffer_data(&buf), carquet_buffer_size(&buf), &arena, &md, &err);
    carquet_buffer_destroy(&buf);

    int rc = 0;
    if (st != CARQUET_OK || md.num_row_groups != 1 ||
        md.row_groups[0].num_columns != 1) {
        printf("[FAIL] geo_types_repeated_field: parse status %d\n", (int)st);
        rc = 1;
    } else {
        const parquet_geospatial_statistics_t* g =
            &md.row_groups[0].columns[0].metadata.geospatial_statistics;
        /* Before the fix: types[64] (== num_types) was overwritten with 5, so
         * num_types ended at 7 and types[6] was clobbered with 7. */
        if (g->num_types != CARQUET_GEO_MAX_TYPES || g->types[6] != 1006) {
            printf("[FAIL] geo_types_repeated_field: num_types=%d types[6]=%d\n",
                   (int)g->num_types, (int)g->types[6]);
            rc = 1;
        }
    }
    carquet_arena_destroy(&arena);
    if (!rc) TEST_PASS("geo_types_repeated_field");
    return rc;
}

static int test_schema_validation(void) {
    /* Baseline: a well-formed one-column schema opens. */
    elem_spec_t good[2] = {
        { "schema", -1, 0, 1, 0 },
        { "a", CARQUET_PHYSICAL_INT64, 0, -1, 0 },
    };
    if (!schema_opens(good, 2)) TEST_FAIL("schema_validation", "valid schema rejected");

    elem_spec_t good16[2] = {
        { "schema", -1, 0, 1, 0 },
        { "h", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 2, -1, 1 },
    };
    if (!schema_opens(good16, 2)) TEST_FAIL("schema_validation", "FLOAT16 FLBA(2) rejected");

    /* Root claims 1 child but 3 leaves follow: leaves b, c were never
     * reached and aliased the root element (levels 0, no type). */
    elem_spec_t short_root[4] = {
        { "schema", -1, 0, 1, 0 },
        { "a", CARQUET_PHYSICAL_INT64, 0, -1, 0 },
        { "b", CARQUET_PHYSICAL_INT64, 0, -1, 0 },
        { "c", CARQUET_PHYSICAL_INT64, 0, -1, 0 },
    };
    if (schema_opens(short_root, 4))
        TEST_FAIL("schema_validation", "inconsistent num_children accepted");

    /* FIXED_LEN_BYTE_ARRAY with type_length 0. */
    elem_spec_t flba0[2] = {
        { "schema", -1, 0, 1, 0 },
        { "f", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 0, -1, 0 },
    };
    if (schema_opens(flba0, 2)) TEST_FAIL("schema_validation", "FLBA(0) accepted");

    /* FLOAT16 on FLBA(1): 2-byte decoders would over-read 1-byte values/stats. */
    elem_spec_t f16_1[2] = {
        { "schema", -1, 0, 1, 0 },
        { "h", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 1, -1, 1 },
    };
    if (schema_opens(f16_1, 2)) TEST_FAIL("schema_validation", "FLOAT16 on FLBA(1) accepted");

    TEST_PASS("schema_validation");
    return 0;
}

static size_t build_dba(uint8_t* buf, size_t cap, const int32_t* pfx,
                        const int32_t* sfx, int n, size_t payload) {
    size_t w1 = 0, w2 = 0;
    if (carquet_delta_encode_int32(pfx, n, buf, cap, &w1) != CARQUET_OK) return 0;
    if (carquet_delta_encode_int32(sfx, n, buf + w1, cap - w1, &w2) != CARQUET_OK) return 0;
    if (w1 + w2 + payload > cap) return 0;
    memset(buf + w1 + w2, 'A', payload);
    return w1 + w2 + payload;
}

static int test_delta_byte_array_lengths(void) {
    static uint8_t buf[1 << 14];
    int32_t *p = NULL, *s = NULL;
    size_t hdr = 0, total = 0;

    /* Every prefix = INT32_MAX with no previous value: used to be accepted
     * with a 16 GiB reconstruction size from a 24-byte page. */
    int32_t pA[8], sA[8];
    for (int i = 0; i < 8; i++) { pA[i] = INT32_MAX; sA[i] = 0; }
    size_t n = build_dba(buf, sizeof(buf), pA, sA, 8, 0);
    if (!n) TEST_FAIL("delta_byte_array_lengths", "encode A");
    carquet_status_t st = carquet_delta_strings_decode_lengths(buf, n, 8, &p, &s, &hdr, &total);
    carquet_mem_free(p); carquet_mem_free(s); p = s = NULL;
    if (st == CARQUET_OK) TEST_FAIL("delta_byte_array_lengths", "prefix > previous length accepted");

    /* Total of 0xFFFFFFF8: wrapped the retain-page allocation on 32-bit. */
    int32_t pB[3] = {0, INT32_MAX, 0x7FFFFC11}, sB[3] = {1000, 0, 0};
    n = build_dba(buf, sizeof(buf), pB, sB, 3, 1000);
    if (!n) TEST_FAIL("delta_byte_array_lengths", "encode B");
    st = carquet_delta_strings_decode_lengths(buf, n, 3, &p, &s, &hdr, &total);
    carquet_mem_free(p); carquet_mem_free(s); p = s = NULL;
    if (st == CARQUET_OK) TEST_FAIL("delta_byte_array_lengths", "wrapping total accepted");

    /* A valid stream still decodes: "abc", "abd" (prefix 2), "x". */
    int32_t pC[3] = {0, 2, 0}, sC[3] = {3, 1, 1};
    n = build_dba(buf, sizeof(buf), pC, sC, 3, 5);
    if (!n) TEST_FAIL("delta_byte_array_lengths", "encode C");
    st = carquet_delta_strings_decode_lengths(buf, n, 3, &p, &s, &hdr, &total);
    carquet_mem_free(p); carquet_mem_free(s);
    if (st != CARQUET_OK || total != 3 + 3 + 1)
        TEST_FAIL("delta_byte_array_lengths", "valid stream rejected");

    /* DELTA_LENGTH_BYTE_ARRAY lengths summing to 2^32 (0 on 32-bit). */
    int32_t lD[3] = {INT32_MAX, INT32_MAX, 2};
    size_t w = 0;
    if (carquet_delta_encode_int32(lD, 3, buf, sizeof(buf), &w) != CARQUET_OK)
        TEST_FAIL("delta_byte_array_lengths", "encode D");
    carquet_byte_array_t out[3];
    size_t consumed = 0;
    if (carquet_delta_length_decode(buf, w, out, 3, &consumed) == CARQUET_OK)
        TEST_FAIL("delta_byte_array_lengths", "DELTA_LENGTH oversized lengths accepted");

    TEST_PASS("delta_byte_array_lengths");
    return 0;
}

/* ------------------------------------------------------------------------- */

static int test_coalesced_split_inflated_num_values(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "sec_split");
    const int64_t nrows = 500000;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema || carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64, NULL,
                                             CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL("coalesced_split", "schema");
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_ZSTD;
    wo.write_crc = false;
    wo.row_group_size = (int64_t)1 << 40;
    wo.max_rows_per_page = 131072;
    carquet_writer_t* w = carquet_writer_create(path, schema, &wo, &err);
    if (!w) { carquet_schema_free(schema); TEST_FAIL("coalesced_split", "writer"); }
    if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN) != CARQUET_OK) {
        carquet_writer_abort(w); carquet_schema_free(schema);
        TEST_FAIL("coalesced_split", "encoding");
    }
    int64_t* vals = (int64_t*)malloc((size_t)nrows * sizeof(int64_t));
    for (int64_t i = 0; i < nrows; i++) vals[i] = i;
    carquet_status_t st = carquet_writer_write_batch(w, 0, vals, nrows, NULL, NULL);
    st = st == CARQUET_OK ? carquet_writer_close(w) : (carquet_writer_abort(w), st);
    free(vals);
    carquet_schema_free(schema);
    if (st != CARQUET_OK) TEST_FAIL("coalesced_split", "write");

    /* Load the file, locate the first data page header's num_values varint
     * (DataPageHeader field 1, compact header byte 0x15) and inflate it to the
     * largest value of the same varint width, so offsets do not move. */
    FILE* f = fopen(path, "rb");
    if (!f) TEST_FAIL("coalesced_split", "reopen");
    fseek(f, 0, SEEK_END);
    long fsz = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t* file = (uint8_t*)malloc((size_t)fsz);
    size_t got = fread(file, 1, (size_t)fsz, f);
    fclose(f);
    if (got != (size_t)fsz) { free(file); TEST_FAIL("coalesced_split", "read"); }

    long page_off = 4;  /* first data page follows the magic */
    parquet_page_header_t hdr;
    size_t hsz = 0;
    if (parquet_parse_page_header(file + page_off, (size_t)(fsz - page_off), &hdr, &hsz,
                                  NULL) != CARQUET_OK || hdr.type != CARQUET_PAGE_DATA) {
        free(file);
        TEST_FAIL("coalesced_split", "first page header");
    }
    uint32_t zz = (uint32_t)hdr.data_page_header.num_values << 1;
    uint8_t want[5];
    int wl = 0;
    for (uint32_t t = zz;;) { uint8_t b = t & 0x7F; t >>= 7; if (t) b |= 0x80; want[wl++] = b; if (!t) break; }
    long found = -1;
    for (long q = page_off; q + 1 + wl <= page_off + (long)hsz; q++) {
        if (file[q] == 0x15 && memcmp(file + q + 1, want, (size_t)wl) == 0) { found = q + 1; break; }
    }
    if (found < 0 || wl >= 5) { free(file); TEST_FAIL("coalesced_split", "num_values varint"); }
    uint32_t zz_new = (((1u << (7 * wl)) - 1u) >> 1) << 1;
    for (int k = 0; k < wl; k++) {
        file[found + k] = (uint8_t)((zz_new >> (7 * k)) & 0x7F) | (k + 1 < wl ? 0x80 : 0);
    }
    f = fopen(path, "wb");
    if (!f) { free(file); TEST_FAIL("coalesced_split", "rewrite"); }
    fwrite(file, 1, (size_t)fsz, f);
    fclose(f);
    free(file);

    /* mmap + 4 threads activates the pipelined, split-column read. */
    carquet_reader_options_t ro;
    carquet_reader_options_init(&ro);
    ro.use_mmap = true;
    carquet_reader_t* r = carquet_reader_open(path, &ro, &err);
    if (!r) { remove(path); TEST_FAIL("coalesced_split", "open"); }
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 65536;
    cfg.num_threads = 4;
    cfg.use_mmap = true;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); remove(path); TEST_FAIL("coalesced_split", "batch reader"); }
    carquet_row_batch_t* batch = NULL;
    int64_t rows = 0;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        rows += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(path);
    /* The corrupt page must not yield more rows than the row group has; the
     * real assertion is "no out-of-bounds write" (ASan / crash). */
    if (rows > nrows) TEST_FAIL("coalesced_split", "more rows than written");
    TEST_PASS("coalesced_split_inflated_num_values");
    return 0;
}

/* ------------------------------------------------------------------------- */

static int test_pipeline_write_error(void) {
#if defined(_WIN32)
    TEST_PASS("pipeline_write_error (skipped: needs RLIMIT_FSIZE)");
    return 0;
#else
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "sec_pipe");

    /* Cap the file size so a background row-group write fails mid-stream. */
    struct rlimit old_rl;
    if (getrlimit(RLIMIT_FSIZE, &old_rl) != 0) {
        TEST_PASS("pipeline_write_error (skipped: getrlimit)");
        return 0;
    }
    void (*old_sig)(int) = signal(SIGXFSZ, SIG_IGN);
    struct rlimit rl = old_rl;
    rl.rlim_cur = 256 * 1024;
    if (setrlimit(RLIMIT_FSIZE, &rl) != 0) {
        signal(SIGXFSZ, old_sig);
        TEST_PASS("pipeline_write_error (skipped: setrlimit)");
        return 0;
    }

    int rc = 0;
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s ||
        carquet_schema_add_column(s, "a", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        rc = 1;
    }
    carquet_writer_t* w = NULL;
    if (!rc) {
        carquet_writer_options_t o;
        carquet_writer_options_init(&o);
        o.compression = CARQUET_COMPRESSION_SNAPPY;
        o.dictionary_encoding = CARQUET_ENCODING_PLAIN;
        o.async_io = true;
        w = carquet_writer_create(path, s, &o, &err);
        if (!w) rc = 1;
    }
    carquet_status_t failed = CARQUET_OK;
    if (!rc) {
        const int64_t n = 200000;
        int64_t* v = (int64_t*)malloc((size_t)n * sizeof(int64_t));
        uint64_t x = 88172645463325252ull;
        for (int64_t i = 0; i < n; i++) {  /* incompressible */
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            v[i] = (int64_t)x;
        }
        for (int rg = 0; rg < 6 && failed == CARQUET_OK; rg++) {
            carquet_status_t st = carquet_writer_write_batch(w, 0, v, n, NULL, NULL);
            if (st == CARQUET_OK) st = carquet_writer_write_batch(w, 1, v, n, NULL, NULL);
            if (st == CARQUET_OK) st = carquet_writer_new_row_group(w);
            failed = st;
        }
        /* Keep writing after the failure: must keep failing, not append to a
         * row group the encode tasks still read. */
        if (failed != CARQUET_OK &&
            carquet_writer_write_batch(w, 0, v, n, NULL, NULL) == CARQUET_OK) {
            printf("[FAIL] pipeline_write_error: write after error succeeded\n");
            rc = 1;
        }
        free(v);
        /* Close after the error used to free the in-flight row group twice. */
        carquet_status_t cst = carquet_writer_close(w);
        if (failed == CARQUET_OK && cst == CARQUET_OK) {
            printf("[FAIL] pipeline_write_error: size limit never hit\n");
            rc = 1;
        } else if (failed != CARQUET_OK && cst == CARQUET_OK) {
            printf("[FAIL] pipeline_write_error: close hid the write error\n");
            rc = 1;
        }
    }
    if (s) carquet_schema_free(s);
    setrlimit(RLIMIT_FSIZE, &old_rl);
    signal(SIGXFSZ, old_sig);
    remove(path);
    if (!rc) TEST_PASS("pipeline_write_error");
    return rc;
#endif
}

/* ------------------------------------------------------------------------- */
/* Arrow import validation (stack-built Arrow structs, no-op release)         */
/* ------------------------------------------------------------------------- */

static void noop_release_schema(struct ArrowSchema* s) { s->release = NULL; }
static void noop_release_array(struct ArrowArray* a) { a->release = NULL; }

static void init_schema(struct ArrowSchema* s, const char* fmt, const char* name,
                        int64_t n_children, struct ArrowSchema** children) {
    memset(s, 0, sizeof(*s));
    s->format = fmt;
    s->name = name;
    s->n_children = n_children;
    s->children = children;
    s->release = noop_release_schema;
}

static void init_array(struct ArrowArray* a, int64_t length, int64_t n_buffers,
                       const void** buffers, int64_t n_children,
                       struct ArrowArray** children) {
    memset(a, 0, sizeof(*a));
    a->length = length;
    a->n_buffers = n_buffers;
    a->buffers = buffers;
    a->n_children = n_children;
    a->children = children;
    a->release = noop_release_array;
}

static carquet_writer_t* make_buffer_writer(carquet_physical_type_t pt,
                                            carquet_field_repetition_t rep,
                                            carquet_schema_t** out_schema) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) return NULL;
    if (carquet_schema_add_column(s, "x", pt, NULL, rep, 0, 0) != CARQUET_OK) {
        carquet_schema_free(s);
        return NULL;
    }
    carquet_writer_t* w = carquet_writer_create_buffer(s, NULL, &err);
    if (!w) { carquet_schema_free(s); return NULL; }
    *out_schema = s;
    return w;
}

static int test_arrow_import_validation(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    int32_t ints[4] = {1, 2, 3, 4};

    /* (1) a leaf into a column of another physical type: write_batch used to
     *     read 8 bytes per value from a 4-byte-per-value buffer. An int32 leaf
     *     into an INT64 column is now widened (0.8.0); a double leaf into an
     *     INT64 column must still be refused. */
    {
        carquet_schema_t* cs = NULL;
        carquet_writer_t* w = make_buffer_writer(CARQUET_PHYSICAL_INT64,
                                                 CARQUET_REPETITION_REQUIRED, &cs);
        if (!w) TEST_FAIL("arrow_import_validation", "writer (1)");
        struct ArrowSchema leaf_s, top_s;
        struct ArrowSchema* kids_s[1] = { &leaf_s };
        init_schema(&leaf_s, "g", "x", 0, NULL);
        init_schema(&top_s, "+s", "", 1, kids_s);
        const void* leaf_bufs[2] = { NULL, ints };
        struct ArrowArray leaf_a, top_a;
        struct ArrowArray* kids_a[1] = { &leaf_a };
        init_array(&leaf_a, 4, 2, leaf_bufs, 0, NULL);
        const void* top_bufs[1] = { NULL };
        init_array(&top_a, 4, 1, top_bufs, 1, kids_a);
        carquet_status_t st = carquet_writer_write_arrow(w, &top_a, &top_s, &err);
        carquet_writer_abort(w);
        carquet_schema_free(cs);
        if (st == CARQUET_OK) TEST_FAIL("arrow_import_validation", "double leaf into INT64 accepted");
    }

    /* (2) nested struct array with fewer children than its schema: the
     *     shredder indexed a->children[1] past the array. */
    {
        carquet_error_t e2 = CARQUET_ERROR_INIT;
        carquet_schema_t* cs = carquet_schema_create(&e2);
        int32_t g = cs ? carquet_schema_add_group(cs, "st", CARQUET_REPETITION_REQUIRED, 0) : -1;
        if (g < 0 ||
            carquet_schema_add_column(cs, "p", CARQUET_PHYSICAL_INT32, NULL,
                                      CARQUET_REPETITION_REQUIRED, 0, g) != CARQUET_OK ||
            carquet_schema_add_column(cs, "q", CARQUET_PHYSICAL_INT32, NULL,
                                      CARQUET_REPETITION_REQUIRED, 0, g) != CARQUET_OK) {
            if (cs) carquet_schema_free(cs);
            TEST_FAIL("arrow_import_validation", "nested schema");
        }
        carquet_writer_t* w = carquet_writer_create_buffer(cs, NULL, &e2);
        if (!w) { carquet_schema_free(cs); TEST_FAIL("arrow_import_validation", "writer (2)"); }

        struct ArrowSchema p_s, q_s, st_s, top_s;
        struct ArrowSchema* st_kids_s[2] = { &p_s, &q_s };
        struct ArrowSchema* top_kids_s[1] = { &st_s };
        init_schema(&p_s, "i", "p", 0, NULL);
        init_schema(&q_s, "i", "q", 0, NULL);
        init_schema(&st_s, "+s", "st", 2, st_kids_s);
        init_schema(&top_s, "+s", "", 1, top_kids_s);

        const void* p_bufs[2] = { NULL, ints };
        struct ArrowArray p_a, st_a, top_a;
        init_array(&p_a, 4, 2, p_bufs, 0, NULL);
        /* children array holds only one entry while the schema says 2 */
        struct ArrowArray** st_kids_a = (struct ArrowArray**)malloc(sizeof(*st_kids_a));
        st_kids_a[0] = &p_a;
        const void* nobuf[1] = { NULL };
        init_array(&st_a, 4, 1, nobuf, 1, st_kids_a);
        struct ArrowArray* top_kids_a[1] = { &st_a };
        init_array(&top_a, 4, 1, nobuf, 1, top_kids_a);
        carquet_status_t st = carquet_writer_write_arrow(w, &top_a, &top_s, &e2);
        free(st_kids_a);
        carquet_writer_abort(w);
        carquet_schema_free(cs);
        if (st == CARQUET_OK) TEST_FAIL("arrow_import_validation", "struct child-count mismatch accepted");
    }

    /* (3) utf8 value with a non-empty length but no data buffer. */
    {
        carquet_schema_t* cs = NULL;
        carquet_writer_t* w = make_buffer_writer(CARQUET_PHYSICAL_BYTE_ARRAY,
                                                 CARQUET_REPETITION_REQUIRED, &cs);
        if (!w) TEST_FAIL("arrow_import_validation", "writer (3)");
        struct ArrowSchema leaf_s, top_s;
        struct ArrowSchema* kids_s[1] = { &leaf_s };
        init_schema(&leaf_s, "z", "x", 0, NULL);
        init_schema(&top_s, "+s", "", 1, kids_s);
        int32_t offs[3] = {0, 5, 9};
        const void* leaf_bufs[3] = { NULL, offs, NULL };
        struct ArrowArray leaf_a, top_a;
        struct ArrowArray* kids_a[1] = { &leaf_a };
        init_array(&leaf_a, 2, 3, leaf_bufs, 0, NULL);
        const void* top_bufs[1] = { NULL };
        init_array(&top_a, 2, 1, top_bufs, 1, kids_a);
        carquet_status_t st = carquet_writer_write_arrow(w, &top_a, &top_s, &err);
        carquet_writer_abort(w);
        carquet_schema_free(cs);
        if (st == CARQUET_OK) TEST_FAIL("arrow_import_validation", "NULL data buffer accepted");
    }

    /* (4) the matching case still works. */
    {
        carquet_schema_t* cs = NULL;
        carquet_writer_t* w = make_buffer_writer(CARQUET_PHYSICAL_INT32,
                                                 CARQUET_REPETITION_REQUIRED, &cs);
        if (!w) TEST_FAIL("arrow_import_validation", "writer (4)");
        struct ArrowSchema leaf_s, top_s;
        struct ArrowSchema* kids_s[1] = { &leaf_s };
        init_schema(&leaf_s, "i", "x", 0, NULL);
        init_schema(&top_s, "+s", "", 1, kids_s);
        const void* leaf_bufs[2] = { NULL, ints };
        struct ArrowArray leaf_a, top_a;
        struct ArrowArray* kids_a[1] = { &leaf_a };
        init_array(&leaf_a, 4, 2, leaf_bufs, 0, NULL);
        const void* top_bufs[1] = { NULL };
        init_array(&top_a, 4, 1, top_bufs, 1, kids_a);
        carquet_status_t st = carquet_writer_write_arrow(w, &top_a, &top_s, &err);
        carquet_writer_abort(w);
        carquet_schema_free(cs);
        if (st != CARQUET_OK) {
            printf("[FAIL] arrow_import_validation: valid import failed: %s\n", err.message);
            return 1;
        }
    }

    TEST_PASS("arrow_import_validation");
    return 0;
}

int main(void) {
    printf("=== Security hardening regression tests ===\n");
    int rc = 0;
    rc |= test_geo_types_repeated_field();
    rc |= test_schema_validation();
    rc |= test_delta_byte_array_lengths();
    rc |= test_coalesced_split_inflated_num_values();
    rc |= test_pipeline_write_error();
    rc |= test_arrow_import_validation();
    if (rc) {
        printf("=== FAILURES ===\n");
        return 1;
    }
    printf("=== All security hardening tests passed ===\n");
    return 0;
}
