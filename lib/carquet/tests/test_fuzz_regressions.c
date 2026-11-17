/**
 * @file test_fuzz_regressions.c
 * @brief File-level regressions for inputs found by the fuzz targets and the
 *        interop suite (interop/run_interop.py full-value corpus oracle).
 *
 * Fuzz cases embed a minimised fuzzer input and check the library's behaviour
 * on it (terminates, rejects cleanly, never spins). Interop cases rebuild the
 * legacy-writer quirk found in apache/parquet-testing from a carquet-written
 * file by rewriting its footer, so the test needs no external data.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <carquet/carquet.h>
#include "test_helpers.h"
#include "core/arena.h"
#include "core/buffer.h"
#include "thrift/parquet_types.h"
#include "thrift/thrift_encode.h"

/* Codec entry points (src/compression/, src/reader/page_reader.c). */
int carquet_gzip_compress(const uint8_t* src, size_t src_size, uint8_t* dst,
                          size_t dst_capacity, size_t* dst_size, int level);
int carquet_gzip_decompress(const uint8_t* src, size_t src_size, uint8_t* dst,
                            size_t dst_capacity, size_t* dst_size);
carquet_status_t carquet_lz4_compress(const uint8_t* src, size_t src_size,
                                      uint8_t* dst, size_t dst_capacity,
                                      size_t* dst_size);
size_t carquet_lz4_compress_bound(size_t src_size);
carquet_status_t carquet_decompress_page(carquet_compression_t codec,
                                         const uint8_t* compressed,
                                         size_t compressed_size,
                                         uint8_t* decompressed,
                                         size_t decompressed_capacity,
                                         size_t* decompressed_size);

/* ----------------------------------------------------------------------------
 * fuzz_reader_file timeout (Oct 2026): a LIST/MAP file whose only row group
 * declares num_rows = 0 while its column chunks declare num_values = 2. The
 * nested batch path returned an empty batch for the empty row group without
 * consuming the column readers, then saw the first column still "has_next"
 * and re-emitted the same empty batch forever — carquet_validate_file() and
 * any `while (next() == OK)` loop hung. Minimised by libFuzzer, then the
 * footer re-framed behind a "PAR1" header (573 bytes).
 * -------------------------------------------------------------------------- */
static const uint8_t k_nested_empty_rg_spin[] = {
    0x50, 0x41, 0x52, 0x31, 0x15, 0x04, 0x19, 0x5c, 0x35, 0x00, 0x18, 0x06, 0x73, 0x63, 0x68, 0x65,
    0x6d, 0x61, 0x15, 0x02, 0x00, 0x35, 0x02, 0x18, 0x03, 0x61, 0x72, 0x72, 0x15, 0x02, 0x15, 0x02,
    0x4c, 0x2c, 0x00, 0x00, 0x00, 0x35, 0x04, 0x18, 0x09, 0x6b, 0x65, 0x79, 0x5f, 0x76, 0x61, 0x6c,
    0x75, 0x65, 0x15, 0x04, 0x00, 0x15, 0x0c, 0x25, 0x00, 0x18, 0x03, 0x6b, 0x65, 0x79, 0x25, 0x00,
    0x4c, 0x1c, 0x00, 0x00, 0x00, 0x15, 0x02, 0x25, 0x02, 0x18, 0x05, 0x76, 0x61, 0x6c, 0x75, 0x65,
    0x00, 0x16, 0x04, 0x19, 0x1c, 0x19, 0x2c, 0x26, 0xac, 0x34, 0x1c, 0x15, 0x0c, 0x19, 0x45, 0x10,
    0x00, 0x06, 0x00, 0x19, 0x38, 0x03, 0x61, 0x72, 0x72, 0x09, 0x6b, 0x65, 0x79, 0x5f, 0x76, 0x61,
    0x6c, 0x75, 0x65, 0x03, 0x6b, 0x65, 0x79, 0x15, 0x08, 0x16, 0x04, 0x16, 0xca, 0x81, 0x80, 0x80,
    0x10, 0x16, 0xa4, 0x34, 0x26, 0xe4, 0x19, 0x26, 0x08, 0x1c, 0x36, 0x00, 0x00, 0x19, 0x3c, 0x15,
    0x04, 0x15, 0x00, 0x15, 0x02, 0x00, 0x15, 0x00, 0x15, 0x00, 0x15, 0x02, 0x00, 0x15, 0x00, 0x15,
    0x10, 0x15, 0x02, 0x00, 0x00, 0x00, 0x26, 0xec, 0x36, 0x1c, 0x15, 0x02, 0x19, 0x35, 0x10, 0x00,
    0x06, 0x19, 0x38, 0x03, 0x61, 0x72, 0x72, 0x09, 0x6b, 0x65, 0x79, 0x5f, 0x76, 0x61, 0x6c, 0x75,
    0x65, 0x05, 0x76, 0x61, 0x6c, 0x75, 0x65, 0x15, 0x08, 0x16, 0x04, 0x16, 0x9c, 0x01, 0x16, 0xa2,
    0x01, 0x26, 0xf6, 0x35, 0x26, 0xca, 0x35, 0x1c, 0x18, 0x04, 0x01, 0x00, 0x00, 0x00, 0x18, 0x04,
    0x01, 0x00, 0x00, 0x00, 0x16, 0x00, 0xe1, 0x04, 0x01, 0x00, 0x00, 0x00, 0x18, 0x04, 0x01, 0x00,
    0x00, 0x00, 0x00, 0x19, 0x2c, 0x15, 0x04, 0x15, 0x00, 0x15, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x0c, 0x42, 0x68, 0x63, 0x6e, 0x49, 0x41, 0x56, 0x50, 0x2f, 0x2f, 0x2f, 0x34,
    0x6a, 0x2f, 0x2f, 0x2f, 0x38, 0x41, 0x41, 0x41, 0x41, 0x4e, 0x47, 0x41, 0x41, 0x41, 0x41, 0x43,
    0x41, 0x41, 0x41, 0x41, 0x41, 0x45, 0x41, 0x41, 0x41, 0x41, 0x41, 0x67, 0x41, 0x41, 0x41, 0x48,
    0x41, 0x41, 0x41, 0x41, 0x41, 0x6b, 0x41, 0x41, 0x41, 0x41, 0x42, 0x77, 0x41, 0x41, 0x41, 0x47,
    0x56, 0x75, 0x41, 0x41, 0x41, 0x3d, 0x3d, 0x00, 0x18, 0x20, 0x70, 0x61, 0x72, 0x71, 0x75, 0x65,
    0x74, 0x2d, 0x63, 0x70, 0x70, 0x2d, 0x61, 0x72, 0x72, 0x6f, 0x77, 0x20, 0x76, 0x65, 0x72, 0x73,
    0x69, 0x6f, 0x6e, 0x20, 0x31, 0x31, 0x2e, 0x30, 0x2e, 0x30, 0x19, 0x2c, 0x1c, 0x00, 0x00, 0x1c,
    0x00, 0x00, 0x00, 0xc8, 0x02, 0x00, 0x00, 0x50, 0x41, 0x52, 0x31, 0x41, 0x41, 0x41, 0x41, 0x41,
    0x41, 0x41, 0x41, 0x41, 0x41, 0x46, 0x41, 0x41, 0x41, 0x41, 0x64, 0x6d, 0x46, 0x73, 0x64, 0x57,
    0x55, 0x41, 0x41, 0x41, 0x41, 0x49, 0x41, 0x41, 0x77, 0x41, 0x43, 0x41, 0x41, 0x48, 0x41, 0x41,
    0x67, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x42, 0x49, 0x41, 0x41, 0x41, 0x41, 0x42,
    0x41, 0x41, 0x46, 0x41, 0x41, 0x49, 0x41, 0x41, 0x41, 0x41, 0x42, 0x77, 0x41, 0x4d, 0x41, 0x41,
    0x41, 0x41, 0x45, 0x41, 0x41, 0x51, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x42, 0x52,
    0x41, 0x41, 0x41, 0x41, 0x41, 0x59, 0x41, 0x41, 0x41, 0x41, 0x42, 0x41, 0x41, 0x41, 0x41, 0x41,
    0x41, 0x41, 0x41, 0x41, 0x41, 0x44, 0x41, 0x41, 0x41, 0x41, 0x61, 0x32, 0x56, 0x35, 0x41, 0x41,
    0x41, 0x42, 0x51, 0x41, 0x41, 0x45, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41,
    0xff, 0xff, 0xff, 0x13, 0x20, 0x70, 0x61, 0x72, 0x71, 0x75, 0x65, 0x74, 0x2d, 0x63, 0x70, 0x70,
    0x2d, 0x61, 0x72, 0x72, 0x6f, 0x77, 0x20, 0x76, 0x65, 0x72, 0x73, 0x69, 0x6f, 0x6e, 0x20, 0x31,
    0x31, 0x2e, 0x30, 0x2e, 0x30, 0x19, 0x2c, 0x1c, 0x00, 0x00, 0x1c, 0x00, 0x00, 0x00, 0xc8, 0x02,
    0x00, 0x00, 0x50, 0x41, 0x52, 0x31, 0x02, 0x00, 0x00, 0x50, 0x41, 0x52, 0x31,
};

static int drain_batches(carquet_reader_t* reader, int* iterations) {
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) return 0;  /* rejecting the file up front is fine too */
    int n = 0;
    for (; n < 1000; n++) {
        carquet_row_batch_t* batch = NULL;
        carquet_status_t st = carquet_batch_reader_next(br, &batch);
        if (st != CARQUET_OK) break;
        carquet_row_batch_free(batch);
    }
    carquet_batch_reader_free(br);
    *iterations = n;
    return n < 1000 ? 0 : 1;
}

static int test_nested_empty_row_group_terminates(void) {
    const char* name = "nested_empty_row_group_terminates";
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open_buffer(
        k_nested_empty_rg_spin, sizeof(k_nested_empty_rg_spin), NULL, &err);
    if (!reader) TEST_FAIL(name, "repro no longer opens (update the test input)");
    int iters = 0;
    int spun = drain_batches(reader, &iters);
    carquet_reader_close(reader);
    if (spun) TEST_FAIL(name, "batch reader did not terminate after 1000 batches");
    if (iters != 0) TEST_FAIL(name, "empty row group should yield no batch");

    /* Same input through the path-based validator (it drives the batch
     * reader internally and used to hang). */
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fuzz_regress_nested_empty_rg");
    FILE* f = fopen(path, "wb");
    if (!f) TEST_FAIL(name, "cannot create temp file");
    size_t w = fwrite(k_nested_empty_rg_spin, 1, sizeof(k_nested_empty_rg_spin), f);
    fclose(f);
    if (w != sizeof(k_nested_empty_rg_spin)) { remove(path); TEST_FAIL(name, "write failed"); }
    err = (carquet_error_t)CARQUET_ERROR_INIT;
    (void)carquet_validate_file(path, &err);  /* any verdict; it must return */
    remove(path);

    TEST_PASS(name);
    return 0;
}

/* ----------------------------------------------------------------------------
 * Found alongside the case above: the minimised input had no leading "PAR1"
 * (its footer started at offset 0) yet carquet_validate_file() accepted it,
 * because the fread open path only reads the file tail. The validator now
 * checks the header magic like the buffer / mmap open paths do.
 * -------------------------------------------------------------------------- */
static int test_validate_rejects_missing_header_magic(void) {
    const char* name = "validate_rejects_missing_header_magic";
    uint8_t buf[sizeof(k_nested_empty_rg_spin)];
    memcpy(buf, k_nested_empty_rg_spin, sizeof(buf));
    memcpy(buf, "XXXX", 4);  /* footer and trailer untouched */

    char path[512];
    carquet_test_temp_path(path, sizeof(path), "fuzz_regress_no_header_magic");
    FILE* f = fopen(path, "wb");
    if (!f) TEST_FAIL(name, "cannot create temp file");
    size_t w = fwrite(buf, 1, sizeof(buf), f);
    fclose(f);
    if (w != sizeof(buf)) { remove(path); TEST_FAIL(name, "write failed"); }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_status_t st = carquet_validate_file(path, &err);
    remove(path);
    if (st != CARQUET_ERROR_INVALID_MAGIC) TEST_FAIL(name, "expected CARQUET_ERROR_INVALID_MAGIC");
    TEST_PASS(name);
    return 0;
}

/* ----------------------------------------------------------------------------
 * Footer rewriting helpers (interop cases)
 * -------------------------------------------------------------------------- */

#define DICT_ROWS 1000

/* One required INT64 column, RLE_DICTIONARY, values i % 7, one row group. */
static int write_dict_file(uint8_t** out, size_t* out_size) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) return 1;
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    carquet_writer_t* w = carquet_writer_create_buffer(s, &o, &err);
    if (!w) { carquet_schema_free(s); return 1; }
    int64_t v[DICT_ROWS];
    for (int i = 0; i < DICT_ROWS; i++) v[i] = i % 7;
    void* buf = NULL;
    size_t sz = 0;
    int rc = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK ||
             carquet_writer_write_batch(w, 0, v, DICT_ROWS, NULL, NULL) != CARQUET_OK ||
             carquet_writer_close(w) != CARQUET_OK ||
             carquet_writer_get_buffer(w, &buf, &sz) != CARQUET_OK;
    carquet_schema_free(s);
    if (rc) { free(buf); return 1; }
    *out = buf;
    *out_size = sz;
    return 0;
}

typedef void (*footer_mutator_t)(parquet_file_metadata_t* md);

/* Parse the footer of `in`, apply `mut`, re-serialize it behind the same
 * column-chunk bytes. Returns a malloc'd file or NULL. */
static uint8_t* rewrite_footer(const uint8_t* in, size_t in_size,
                               footer_mutator_t mut, size_t* out_size) {
    if (in_size < 12) return NULL;
    uint32_t flen = (uint32_t)in[in_size - 8] | ((uint32_t)in[in_size - 7] << 8) |
                    ((uint32_t)in[in_size - 6] << 16) | ((uint32_t)in[in_size - 5] << 24);
    if ((size_t)flen + 12 > in_size) return NULL;
    size_t fstart = in_size - 8 - flen;

    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) return NULL;
    parquet_file_metadata_t md;
    carquet_error_t err = CARQUET_ERROR_INIT;
    uint8_t* out = NULL;
    if (parquet_parse_file_metadata(in + fstart, flen, &arena, &md, &err) == CARQUET_OK) {
        mut(&md);
        carquet_buffer_t b;
        carquet_buffer_init(&b);
        if (parquet_write_file_metadata(&md, &b, &err) == CARQUET_OK) {
            size_t nf = carquet_buffer_size(&b);
            *out_size = fstart + nf + 8;
            out = malloc(*out_size);
            if (out) {
                memcpy(out, in, fstart);
                memcpy(out + fstart, carquet_buffer_data(&b), nf);
                for (int k = 0; k < 4; k++) out[fstart + nf + k] = (uint8_t)(nf >> (8 * k));
                memcpy(out + fstart + nf + 4, "PAR1", 4);
            }
        }
        carquet_buffer_destroy(&b);
        parquet_file_metadata_free(&md);
    }
    carquet_arena_destroy(&arena);
    return out;
}

/* Legacy parquet-mr: no dictionary_page_offset; data_page_offset points at
 * the dictionary page that leads the chunk. */
static void drop_dictionary_offset(parquet_file_metadata_t* md) {
    parquet_column_metadata_t* m = &md->row_groups[0].columns[0].metadata;
    m->data_page_offset = m->dictionary_page_offset;
    m->has_dictionary_page_offset = false;
    m->dictionary_page_offset = 0;
}

/* Dremio / parquet-mr fork: dictionary_page_offset written as 0. */
static void zero_dictionary_offset(parquet_file_metadata_t* md) {
    parquet_column_metadata_t* m = &md->row_groups[0].columns[0].metadata;
    m->data_page_offset = m->dictionary_page_offset;
    m->dictionary_page_offset = 0;
}

/* Decode the whole column through the batch reader and the low-level column
 * reader (with a skip); every value must be i % 7. */
static int check_dict_values(carquet_reader_t* r, const char** why) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { *why = "batch reader create"; return 1; }
    int64_t row = 0;
    carquet_row_batch_t* b = NULL;
    while (carquet_batch_reader_next(br, &b) == CARQUET_OK && b) {
        const void* d; const uint8_t* nulls; int64_t n;
        if (carquet_row_batch_column(b, 0, &d, &nulls, &n) != CARQUET_OK) {
            carquet_row_batch_free(b); carquet_batch_reader_free(br);
            *why = "batch column"; return 1;
        }
        for (int64_t i = 0; i < n; i++, row++) {
            if (((const int64_t*)d)[i] != row % 7) {
                carquet_row_batch_free(b); carquet_batch_reader_free(br);
                *why = "batch value mismatch"; return 1;
            }
        }
        carquet_row_batch_free(b);
        b = NULL;
    }
    carquet_batch_reader_free(br);
    if (row != DICT_ROWS) { *why = "batch row count"; return 1; }

    carquet_column_reader_t* cr = carquet_reader_get_column(r, 0, 0, &err);
    if (!cr) { *why = "column reader"; return 1; }
    int64_t vals[DICT_ROWS];
    int64_t skipped = carquet_column_skip(cr, 10);
    int64_t got = carquet_column_read_batch(cr, vals, DICT_ROWS, NULL, NULL);
    carquet_column_reader_free(cr);
    if (skipped != 10 || got != DICT_ROWS - 10) { *why = "skip/read count"; return 1; }
    for (int64_t i = 0; i < got; i++)
        if (vals[i] != (i + 10) % 7) { *why = "skip/read value mismatch"; return 1; }
    return 0;
}

static int run_dict_quirk(const char* name, footer_mutator_t mut) {
    uint8_t* orig = NULL;
    size_t orig_size = 0;
    if (write_dict_file(&orig, &orig_size) != 0) TEST_FAIL(name, "write");
    size_t sz = 0;
    uint8_t* f = rewrite_footer(orig, orig_size, mut, &sz);
    free(orig);
    if (!f) TEST_FAIL(name, "footer rewrite");

    const char* why = NULL;
    carquet_error_t err = CARQUET_ERROR_INIT;
    /* In-memory (mmap-style) page path. */
    carquet_reader_t* r = carquet_reader_open_buffer(f, sz, NULL, &err);
    int bad = !r ? (why = "open_buffer", 1) : check_dict_values(r, &why);
    carquet_reader_close(r);
    /* fread page path. */
    if (!bad) {
        char path[512];
        carquet_test_temp_path(path, sizeof(path), name);
        FILE* fp = fopen(path, "wb");
        if (fp) { fwrite(f, 1, sz, fp); fclose(fp); }
        err = (carquet_error_t)CARQUET_ERROR_INIT;
        r = carquet_reader_open(path, NULL, &err);
        bad = !r ? (why = "open(path)", 1) : check_dict_values(r, &why);
        carquet_reader_close(r);
        remove(path);
    }
    free(f);
    if (bad) TEST_FAIL(name, why);
    TEST_PASS(name);
    return 0;
}

/* parquet-testing alltypes_tiny_pages / hadoop_lz4_compressed: "Expected data
 * page" — the leading dictionary page was only loaded via
 * dictionary_page_offset. */
static int test_dictionary_without_offset(void) {
    return run_dict_quirk("dictionary_without_offset", drop_dictionary_offset);
}

/* parquet-testing dict-page-offset-zero: offset 0 sent the dictionary loader
 * to the "PAR1" magic. */
static int test_dictionary_offset_zero(void) {
    return run_dict_quirk("dictionary_offset_zero", zero_dictionary_offset);
}

/* parquet-testing dict-page-offset-zero (footer side): Dremio wrote
 * ColumnMetaData field 15 as list<struct>, the spec's i32
 * bloom_filter_length. Decoding it as i32 desynchronised the footer
 * ("Cannot skip STOP type"); a wire-type mismatch must be skipped. */
static int test_column_metadata_wire_type_mismatch(void) {
    const char* name = "column_metadata_wire_type_mismatch";
    carquet_buffer_t b;
    carquet_buffer_init(&b);
    thrift_encoder_t e;
    thrift_encoder_init(&e, &b);

    thrift_write_struct_begin(&e);                         /* FileMetaData */
    thrift_write_field_header(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 1);
    thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
    thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 2);
      thrift_write_struct_begin(&e);                       /* root */
      thrift_write_field_header(&e, THRIFT_TYPE_BINARY, 4); thrift_write_string(&e, "root");
      thrift_write_field_header(&e, THRIFT_TYPE_I32, 5); thrift_write_i32(&e, 1);
      thrift_write_struct_end(&e);
      thrift_write_struct_begin(&e);                       /* leaf a: INT32 */
      thrift_write_field_header(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 1);
      thrift_write_field_header(&e, THRIFT_TYPE_I32, 3); thrift_write_i32(&e, 0);
      thrift_write_field_header(&e, THRIFT_TYPE_BINARY, 4); thrift_write_string(&e, "a");
      thrift_write_struct_end(&e);
    thrift_write_field_header(&e, THRIFT_TYPE_I64, 3); thrift_write_i64(&e, 3);
    thrift_write_field_header(&e, THRIFT_TYPE_LIST, 4);
    thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 1);
      thrift_write_struct_begin(&e);                       /* RowGroup */
      thrift_write_field_header(&e, THRIFT_TYPE_LIST, 1);
      thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 1);
        thrift_write_struct_begin(&e);                     /* ColumnChunk */
        thrift_write_field_header(&e, THRIFT_TYPE_I64, 2); thrift_write_i64(&e, 4);
        thrift_write_field_header(&e, THRIFT_TYPE_STRUCT, 3);
          thrift_write_struct_begin(&e);                   /* ColumnMetaData */
          thrift_write_field_header(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 1);
          thrift_write_field_header(&e, THRIFT_TYPE_LIST, 2);
          thrift_write_list_begin(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 0);
          thrift_write_field_header(&e, THRIFT_TYPE_LIST, 3);
          thrift_write_list_begin(&e, THRIFT_TYPE_BINARY, 1); thrift_write_string(&e, "a");
          thrift_write_field_header(&e, THRIFT_TYPE_I32, 4); thrift_write_i32(&e, 0);
          thrift_write_field_header(&e, THRIFT_TYPE_I64, 5); thrift_write_i64(&e, 3);
          thrift_write_field_header(&e, THRIFT_TYPE_I64, 6); thrift_write_i64(&e, 40);
          thrift_write_field_header(&e, THRIFT_TYPE_I64, 7); thrift_write_i64(&e, 40);
          thrift_write_field_header(&e, THRIFT_TYPE_I64, 9); thrift_write_i64(&e, 4);
          /* Field 15 as list<struct{1: struct{..., 5: struct{...}}, 2: i64}>
           * (the shape Dremio wrote) instead of i32. */
          thrift_write_field_header(&e, THRIFT_TYPE_LIST, 15);
          thrift_write_list_begin(&e, THRIFT_TYPE_STRUCT, 1);
            thrift_write_struct_begin(&e);
            thrift_write_field_header(&e, THRIFT_TYPE_STRUCT, 1);
              thrift_write_struct_begin(&e);
              thrift_write_field_header(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 0);
              thrift_write_field_header(&e, THRIFT_TYPE_I32, 2); thrift_write_i32(&e, 162);
              thrift_write_field_header(&e, THRIFT_TYPE_I32, 3); thrift_write_i32(&e, 22);
              thrift_write_field_header(&e, THRIFT_TYPE_STRUCT, 5);
                thrift_write_struct_begin(&e);
                thrift_write_field_header(&e, THRIFT_TYPE_I32, 1); thrift_write_i32(&e, 39);
                thrift_write_field_header(&e, THRIFT_TYPE_I32, 2); thrift_write_i32(&e, 0);
                thrift_write_field_header(&e, THRIFT_TYPE_I32, 3); thrift_write_i32(&e, 3);
                thrift_write_field_header(&e, THRIFT_TYPE_I32, 4); thrift_write_i32(&e, 4);
                thrift_write_struct_end(&e);
              thrift_write_struct_end(&e);
            thrift_write_field_header(&e, THRIFT_TYPE_I64, 2); thrift_write_i64(&e, 22);
            thrift_write_struct_end(&e);
          thrift_write_struct_end(&e);
        thrift_write_struct_end(&e);
      thrift_write_field_header(&e, THRIFT_TYPE_I64, 2); thrift_write_i64(&e, 40);
      thrift_write_field_header(&e, THRIFT_TYPE_I64, 3); thrift_write_i64(&e, 3);
      thrift_write_struct_end(&e);
    thrift_write_field_header(&e, THRIFT_TYPE_BINARY, 6); thrift_write_string(&e, "w");
    thrift_write_struct_end(&e);                           /* writes STOP */

    carquet_arena_t arena;
    carquet_arena_init(&arena);
    parquet_file_metadata_t md;
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_status_t st = thrift_encoder_has_error(&e) ? CARQUET_ERROR_INTERNAL
        : parquet_parse_file_metadata(carquet_buffer_data(&b), carquet_buffer_size(&b),
                                      &arena, &md, &err);
    int ok = st == CARQUET_OK && md.num_row_groups == 1 &&
             md.row_groups[0].num_rows == 3 &&
             md.row_groups[0].columns[0].metadata.data_page_offset == 4 &&
             md.row_groups[0].columns[0].metadata.num_values == 3 &&
             md.created_by && strcmp(md.created_by, "w") == 0;
    if (!ok && st == CARQUET_OK) {
        printf("  rgs=%d rows=%lld dpo=%lld nv=%lld created_by=%s\n", md.num_row_groups,
               md.num_row_groups ? (long long)md.row_groups[0].num_rows : -1LL,
               md.num_row_groups ? (long long)md.row_groups[0].columns[0].metadata.data_page_offset : -1LL,
               md.num_row_groups ? (long long)md.row_groups[0].columns[0].metadata.num_values : -1LL,
               md.created_by ? md.created_by : "(null)");
    }
    if (st == CARQUET_OK) parquet_file_metadata_free(&md);
    carquet_arena_destroy(&arena);
    carquet_buffer_destroy(&b);
    if (!ok) {
        printf("  status=%d (%s)\n", st, err.message);
        TEST_FAIL(name, "footer with a mistyped field 15 did not parse intact");
    }
    TEST_PASS(name);
    return 0;
}

/* parquet-testing concatenated_gzip_members: a page holding two gzip
 * members decoded only the first. */
static int test_gzip_concatenated_members(void) {
    const char* name = "gzip_concatenated_members";
    uint8_t a[300], bb[200], comp[2048], out[600];
    for (int i = 0; i < 300; i++) a[i] = (uint8_t)(i * 7);
    for (int i = 0; i < 200; i++) bb[i] = (uint8_t)(255 - i);
    size_t ca = 0, cb = 0, n = 0;
    if (carquet_gzip_compress(a, sizeof(a), comp, sizeof(comp), &ca, 6) != CARQUET_OK ||
        carquet_gzip_compress(bb, sizeof(bb), comp + ca, sizeof(comp) - ca, &cb, 6) != CARQUET_OK)
        TEST_FAIL(name, "compress");
    if (carquet_gzip_decompress(comp, ca + cb, out, sizeof(out), &n) != CARQUET_OK ||
        n != sizeof(a) + sizeof(bb) || memcmp(out, a, sizeof(a)) != 0 ||
        memcmp(out + sizeof(a), bb, sizeof(bb)) != 0)
        TEST_FAIL(name, "second member not decoded");
    TEST_PASS(name);
    return 0;
}

/* parquet-testing non_hadoop_lz4_compressed: codec 5 (LZ4) holding a raw
 * LZ4 block without Hadoop framing must still decode. */
static int test_lz4_codec5_raw_block(void) {
    const char* name = "lz4_codec5_raw_block";
    uint8_t src[512], out[512];
    for (int i = 0; i < 512; i++) src[i] = (uint8_t)(i % 13);
    size_t cap = carquet_lz4_compress_bound(sizeof(src));
    uint8_t* comp = malloc(cap);
    size_t cn = 0, n = 0;
    int bad = !comp ||
              carquet_lz4_compress(src, sizeof(src), comp, cap, &cn) != CARQUET_OK ||
              carquet_decompress_page(CARQUET_COMPRESSION_LZ4, comp, cn, out,
                                      sizeof(out), &n) != CARQUET_OK ||
              n != sizeof(src) || memcmp(src, out, sizeof(src)) != 0;
    free(comp);
    if (bad) TEST_FAIL(name, "raw LZ4 block under codec 5 rejected");
    TEST_PASS(name);
    return 0;
}

/* Found by the CLI sweep over fuzz corpora: a leaf without a name made
 * carquet_schema_column_name() return NULL and the CLI strlen() it. */
static void drop_leaf_name(parquet_file_metadata_t* md) {
    md->schema[1].name = NULL;
}

static int test_column_name_never_null(void) {
    const char* name = "column_name_never_null";
    uint8_t* orig = NULL;
    size_t orig_size = 0, sz = 0;
    if (write_dict_file(&orig, &orig_size) != 0) TEST_FAIL(name, "write");
    uint8_t* f = rewrite_footer(orig, orig_size, drop_leaf_name, &sz);
    free(orig);
    if (!f) TEST_FAIL(name, "footer rewrite");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(f, sz, NULL, &err);
    int bad = 0;
    if (r) {   /* rejecting the file is fine; a NULL name is not */
        const char* n0 = carquet_schema_column_name(carquet_reader_schema(r), 0);
        bad = (n0 == NULL);
        carquet_reader_close(r);
    }
    free(f);
    if (bad) TEST_FAIL(name, "carquet_schema_column_name returned NULL");
    TEST_PASS(name);
    return 0;
}

int main(void) {
    (void)carquet_init();
    printf("=== Fuzz Regression Tests ===\n\n");
    int failures = 0;
    failures += test_nested_empty_row_group_terminates();
    failures += test_validate_rejects_missing_header_magic();
    failures += test_dictionary_without_offset();
    failures += test_dictionary_offset_zero();
    failures += test_column_metadata_wire_type_mismatch();
    failures += test_gzip_concatenated_members();
    failures += test_lz4_codec5_raw_block();
    failures += test_column_name_never_null();
    printf("\n%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
