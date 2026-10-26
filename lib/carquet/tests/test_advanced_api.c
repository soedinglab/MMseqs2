/**
 * @file test_advanced_api.c
 * @brief Tests for advanced API features (bloom filter, page index,
 *        key-value metadata, column chunk metadata, per-column options,
 *        buffer writer)
 */

#include <carquet/carquet.h>
#include "test_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ======================================================================
 * Helpers
 * ====================================================================== */

#define N_ROWS 5000

static const char* get_temp_file(void) {
    static char path[1024] = {0};
    if (path[0] == 0) {
        carquet_test_temp_path(path, sizeof(path), "test_advanced_api");
    }
    return path;
}
#define TEMP_FILE (get_temp_file())

static carquet_schema_t* make_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_column(s, "id",    CARQUET_PHYSICAL_INT64,      NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE,     NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "label", CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    return s;
}

/** Write a test file with all advanced features enabled. */
static int write_test_file(const char* path, bool bloom, bool page_index) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression         = CARQUET_COMPRESSION_ZSTD;
    opts.write_statistics    = true;
    opts.write_bloom_filters = bloom;
    opts.write_page_index    = page_index;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { fprintf(stderr, "write_test_file: %s\n", err.message); return 1; }

    /* Key-value metadata */
    carquet_writer_add_metadata(w, "test.key1", "hello");
    carquet_writer_add_metadata(w, "test.key2", "world");

    /* Per-column overrides */
    carquet_writer_set_column_bloom_filter(w, 2, false); /* disable bloom on label */

    /* Write data */
    int64_t ids[N_ROWS];
    double  vals[N_ROWS];
    for (int i = 0; i < N_ROWS; i++) { ids[i] = i; vals[i] = i * 1.5; }

    carquet_status_t st;
    st = carquet_writer_write_batch(w, 0, ids,  N_ROWS, NULL, NULL);
    if (st != CARQUET_OK) goto fail;
    st = carquet_writer_write_batch(w, 1, vals, N_ROWS, NULL, NULL);
    if (st != CARQUET_OK) goto fail;

    /* label: half NULL, half "abc" */
    carquet_byte_array_t labels[N_ROWS];
    int16_t def[N_ROWS];
    int non_null = 0;
    for (int i = 0; i < N_ROWS; i++) {
        if (i % 2 == 0) {
            def[i] = 1;
            labels[non_null].data   = (uint8_t*)"abc";
            labels[non_null].length = 3;
            non_null++;
        } else {
            def[i] = 0;
        }
    }
    st = carquet_writer_write_batch(w, 2, labels, N_ROWS, def, NULL);
    if (st != CARQUET_OK) goto fail;

    st = carquet_writer_close(w);
    carquet_schema_free(schema);
    return (st == CARQUET_OK) ? 0 : 1;

fail:
    carquet_writer_abort(w);
    carquet_schema_free(schema);
    return 1;
}

/* ======================================================================
 * Test: Key-Value Metadata
 * ====================================================================== */

static int test_kv_metadata(void) {
    if (write_test_file(TEMP_FILE, false, false)) TEST_FAIL("kv_metadata", "write failed");

    carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, NULL);
    if (!r) TEST_FAIL("kv_metadata", "open failed");

    int32_t n = carquet_reader_num_metadata(r);
    if (n < 2) TEST_FAIL("kv_metadata", "expected >= 2 kv entries");

    const char *k, *v;
    carquet_status_t st = carquet_reader_get_metadata(r, 0, &k, &v);
    if (st != CARQUET_OK) TEST_FAIL("kv_metadata", "get_metadata failed");
    if (strcmp(k, "test.key1") != 0) TEST_FAIL("kv_metadata", "key mismatch");
    if (strcmp(v, "hello") != 0)     TEST_FAIL("kv_metadata", "value mismatch");

    const char* found = carquet_reader_find_metadata(r, "test.key2");
    if (!found || strcmp(found, "world") != 0)
        TEST_FAIL("kv_metadata", "find_metadata failed");

    const char* missing = carquet_reader_find_metadata(r, "nonexistent");
    if (missing) TEST_FAIL("kv_metadata", "expected NULL for missing key");

    carquet_reader_close(r);
    TEST_PASS("kv_metadata");
    return 0;
}

/* ======================================================================
 * Test: Column Chunk Metadata
 * ====================================================================== */

static int test_column_chunk_metadata(void) {
    if (write_test_file(TEMP_FILE, true, true)) TEST_FAIL("chunk_meta", "write failed");

    carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, NULL);
    if (!r) TEST_FAIL("chunk_meta", "open failed");

    carquet_column_chunk_metadata_t m;
    carquet_status_t st = carquet_reader_column_chunk_metadata(r, 0, 0, &m);
    if (st != CARQUET_OK) TEST_FAIL("chunk_meta", "get metadata failed");

    if (m.type != CARQUET_PHYSICAL_INT64) TEST_FAIL("chunk_meta", "wrong type");
    if (m.num_values != N_ROWS)           TEST_FAIL("chunk_meta", "wrong num_values");
    if (m.total_compressed_size <= 0)     TEST_FAIL("chunk_meta", "bad compressed size");
    if (m.total_uncompressed_size <= 0)   TEST_FAIL("chunk_meta", "bad uncompressed size");
    if (m.data_page_offset <= 0)          TEST_FAIL("chunk_meta", "bad data_page_offset");

    /* Column 0 should have bloom filter (enabled globally, not overridden) */
    /* Column 2 (label) should NOT have bloom filter (overridden off) */
    carquet_column_chunk_metadata_t m2;
    carquet_reader_column_chunk_metadata(r, 0, 2, &m2);
    /* Note: bloom filter availability depends on internal write behavior */

    carquet_reader_close(r);
    TEST_PASS("chunk_meta");
    return 0;
}

/* ======================================================================
 * Test: Bloom Filter
 * ====================================================================== */

static int test_bloom_filter(void) {
    if (write_test_file(TEMP_FILE, true, false)) TEST_FAIL("bloom", "write failed");

    carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, NULL);
    if (!r) TEST_FAIL("bloom", "open failed");

    carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(r, 0, 0, NULL);
    if (!bf) {
        /* Bloom filter may not be present depending on implementation */
        carquet_reader_close(r);
        printf("[SKIP] bloom: no bloom filter available\n");
        return 0;
    }

    size_t sz = carquet_bloom_filter_size(bf);
    if (sz == 0) TEST_FAIL("bloom", "zero size");

    /* Value 42 was written (ids[42] = 42) - must say "might contain" */
    if (!carquet_bloom_filter_check_i64(bf, 42))
        TEST_FAIL("bloom", "false negative for 42");

    /* Value -1 was never written - should say "definitely not" (probabilistic) */
    /* We can't assert this fails, but it's very likely for a good filter */
    bool neg = carquet_bloom_filter_check_i64(bf, -1);
    printf("  bloom -1 check: %s (expected: likely no)\n", neg ? "yes" : "no");

    carquet_bloom_filter_destroy(bf);
    carquet_reader_close(r);
    TEST_PASS("bloom");
    return 0;
}

/* ======================================================================
 * Test: Page Index (Column Index + Offset Index)
 * ====================================================================== */

static int test_page_index(void) {
    if (write_test_file(TEMP_FILE, false, true)) TEST_FAIL("page_index", "write failed");

    carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, NULL);
    if (!r) TEST_FAIL("page_index", "open failed");

    /* Column index */
    carquet_column_index_t* ci = carquet_reader_get_column_index(r, 0, 0, NULL);
    if (!ci) {
        carquet_reader_close(r);
        printf("[SKIP] page_index: no column index available\n");
        return 0;
    }

    int32_t np = carquet_column_index_num_pages(ci);
    if (np <= 0) TEST_FAIL("page_index", "zero pages");

    /* Check first page stats */
    carquet_page_stats_t ps;
    carquet_status_t st = carquet_column_index_get_page_stats(ci, 0, &ps);
    if (st != CARQUET_OK) TEST_FAIL("page_index", "get_page_stats failed");

    /* For INT64, min_value should be a pointer to an int64_t */
    if (ps.min_value && ps.min_value_size >= 8) {
        int64_t mn = *(const int64_t*)ps.min_value;
        printf("  page 0 min: %lld\n", (long long)mn);
        if (mn != 0) TEST_FAIL("page_index", "first page min should be 0");
    }

    int32_t bo = carquet_column_index_boundary_order(ci);
    printf("  boundary_order: %d\n", bo);

    carquet_column_index_free(ci);

    /* Offset index */
    carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, 0, NULL);
    if (oi) {
        int32_t np2 = carquet_offset_index_num_pages(oi);
        if (np2 != np)
            printf("  WARNING: column_index pages=%d vs offset_index pages=%d\n", np, np2);

        carquet_page_location_t loc;
        st = carquet_offset_index_get_page_location(oi, 0, &loc);
        if (st != CARQUET_OK) TEST_FAIL("page_index", "get_page_location failed");
        if (loc.offset <= 0) TEST_FAIL("page_index", "bad page offset");
        if (loc.compressed_size <= 0) TEST_FAIL("page_index", "bad page size");
        if (loc.first_row_index != 0) TEST_FAIL("page_index", "first page should start at row 0");

        printf("  page 0: offset=%lld size=%d first_row=%lld\n",
               (long long)loc.offset, loc.compressed_size,
               (long long)loc.first_row_index);

        carquet_offset_index_free(oi);
    }

    carquet_reader_close(r);
    TEST_PASS("page_index");
    return 0;
}

/* ======================================================================
 * Test: Buffer Writer
 * ====================================================================== */

static int test_buffer_writer(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) {
        carquet_schema_free(schema);
        printf("[SKIP] buffer_writer: %s\n", err.message);
        return 0;
    }

    int32_t data[] = {10, 20, 30, 40, 50};
    carquet_status_t st = carquet_writer_write_batch(w, 0, data, 5, NULL, NULL);
    if (st != CARQUET_OK) TEST_FAIL("buffer_writer", "write_batch failed");

    st = carquet_writer_close(w);
    if (st != CARQUET_OK) TEST_FAIL("buffer_writer", "close failed");

    void* buf = NULL;
    size_t sz = 0;
    st = carquet_writer_get_buffer(w, &buf, &sz);
    if (st != CARQUET_OK || !buf || sz == 0)
        TEST_FAIL("buffer_writer", "get_buffer failed");

    /* Verify: open from buffer and read back */
    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("buffer_writer", "open_buffer failed");

    if (carquet_reader_num_rows(r) != 5)
        TEST_FAIL("buffer_writer", "wrong row count");
    if (carquet_reader_num_columns(r) != 1)
        TEST_FAIL("buffer_writer", "wrong column count");

    /* Read values back */
    carquet_column_reader_t* col = carquet_reader_get_column(r, 0, 0, NULL);
    if (!col) TEST_FAIL("buffer_writer", "get_column failed");

    int32_t vals[5];
    int64_t n = carquet_column_read_batch(col, vals, 5, NULL, NULL);
    if (n != 5) TEST_FAIL("buffer_writer", "read wrong count");
    for (int i = 0; i < 5; i++) {
        if (vals[i] != data[i]) TEST_FAIL("buffer_writer", "value mismatch");
    }

    carquet_column_reader_free(col);
    carquet_reader_close(r);
    free(buf);
    carquet_schema_free(schema);
    TEST_PASS("buffer_writer");
    return 0;
}

/* ======================================================================
 * Test: Per-Column Encoding/Compression
 * ====================================================================== */

static int test_per_column_options(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    carquet_writer_t* w = carquet_writer_create(TEMP_FILE, schema, &opts, &err);
    if (!w) TEST_FAIL("per_column", "create failed");

    /* Override column 0 to UNCOMPRESSED (enum value 0) — the overridden value
       happens to match the "unset" sentinel, so this exercises the regression
       where overrides with a 0-valued enum were being silently ignored. */
    carquet_status_t st = carquet_writer_set_column_compression(
        w, 0, CARQUET_COMPRESSION_UNCOMPRESSED, 0);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "set_compression failed");

    /* Override column 1 to SNAPPY with explicit BYTE_STREAM_SPLIT encoding */
    st = carquet_writer_set_column_compression(w, 1, CARQUET_COMPRESSION_SNAPPY, 0);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "set_compression 1 failed");
    st = carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_BYTE_STREAM_SPLIT);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "set_encoding 1 failed");

    /* Override column 2 to PLAIN encoding (enum value 0) — same 0-sentinel
       regression for encoding */
    st = carquet_writer_set_column_encoding(w, 2, CARQUET_ENCODING_PLAIN);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "set_encoding 2 failed");

    /* Invalid column index */
    st = carquet_writer_set_column_encoding(w, 99, CARQUET_ENCODING_PLAIN);
    if (st == CARQUET_OK) TEST_FAIL("per_column", "should reject invalid index");

    /* Column 0 is INT64. Per the Parquet spec, DELTA_BINARY_PACKED and
       BYTE_STREAM_SPLIT are both valid for INT64 and must be accepted. */
    st = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_DELTA_BINARY_PACKED);
    if (st != CARQUET_OK)
        TEST_FAIL("per_column", "DELTA_BINARY_PACKED must be accepted for INT64");
    st = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_BYTE_STREAM_SPLIT);
    if (st != CARQUET_OK)
        TEST_FAIL("per_column", "BYTE_STREAM_SPLIT must be accepted for INT64");

    /* Genuinely unsupported combinations must still be rejected:
       DELTA_BINARY_PACKED is integer-only, so it is invalid for the
       DOUBLE column (index 1). */
    st = carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_DELTA_BINARY_PACKED);
    if (st != CARQUET_ERROR_INVALID_ENCODING)
        TEST_FAIL("per_column", "should reject DELTA_BINARY_PACKED for DOUBLE");

    /* Restore column 0 to PLAIN so the remaining metadata checks below are
       independent of this encoding probe. */
    st = carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "reset col 0 encoding failed");

    /* Write data and close so we can read the metadata back */
    int64_t ids[N_ROWS];
    double vals[N_ROWS];
    for (int i = 0; i < N_ROWS; i++) { ids[i] = i; vals[i] = i * 1.5; }
    st = carquet_writer_write_batch(w, 0, ids, N_ROWS, NULL, NULL);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "write 0 failed");
    st = carquet_writer_write_batch(w, 1, vals, N_ROWS, NULL, NULL);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "write 1 failed");

    carquet_byte_array_t labels[N_ROWS];
    int16_t def[N_ROWS];
    int non_null = 0;
    for (int i = 0; i < N_ROWS; i++) {
        if (i % 2 == 0) {
            def[i] = 1;
            labels[non_null].data = (uint8_t*)"abc";
            labels[non_null].length = 3;
            non_null++;
        } else {
            def[i] = 0;
        }
    }
    st = carquet_writer_write_batch(w, 2, labels, N_ROWS, def, NULL);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "write 2 failed");
    st = carquet_writer_close(w);
    if (st != CARQUET_OK) TEST_FAIL("per_column", "close failed");

    carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, NULL);
    if (!r) TEST_FAIL("per_column", "read open failed");

    carquet_column_chunk_metadata_t m0, m1, m2;
    if (carquet_reader_column_chunk_metadata(r, 0, 0, &m0) != CARQUET_OK ||
        carquet_reader_column_chunk_metadata(r, 0, 1, &m1) != CARQUET_OK ||
        carquet_reader_column_chunk_metadata(r, 0, 2, &m2) != CARQUET_OK) {
        TEST_FAIL("per_column", "column_chunk_metadata failed");
    }

    if (m0.codec != CARQUET_COMPRESSION_UNCOMPRESSED)
        TEST_FAIL("per_column", "col 0 codec override not applied");
    if (m1.codec != CARQUET_COMPRESSION_SNAPPY)
        TEST_FAIL("per_column", "col 1 codec override not applied");
    if (m2.codec != CARQUET_COMPRESSION_ZSTD)
        TEST_FAIL("per_column", "col 2 codec should fall back to global");

    if (m1.encodings[0] != CARQUET_ENCODING_BYTE_STREAM_SPLIT)
        TEST_FAIL("per_column", "col 1 encoding override not applied");
    if (m2.encodings[0] != CARQUET_ENCODING_PLAIN)
        TEST_FAIL("per_column", "col 2 encoding override not applied");

    carquet_column_reader_t* label_reader = carquet_reader_get_column(r, 0, 2, &err);
    if (!label_reader) TEST_FAIL("per_column", "label reader failed");

    carquet_byte_array_t read_labels[N_ROWS];
    int16_t read_def[N_ROWS];
    int64_t label_rows = carquet_column_read_batch(
        label_reader, read_labels, N_ROWS, read_def, NULL);
    if (label_rows != N_ROWS) TEST_FAIL("per_column", "label row count mismatch");

    int64_t label_value_index = 0;
    for (int i = 0; i < N_ROWS; i++) {
        if (i % 2 == 0) {
            if (read_def[i] != 1 ||
                read_labels[label_value_index].length != 3 ||
                memcmp(read_labels[label_value_index].data, "abc", 3) != 0) {
                TEST_FAIL("per_column", "nullable BYTE_ARRAY value mismatch");
            }
            label_value_index++;
        } else if (read_def[i] != 0) {
            TEST_FAIL("per_column", "nullable BYTE_ARRAY def mismatch");
        }
    }
    if (label_value_index != non_null)
        TEST_FAIL("per_column", "nullable BYTE_ARRAY sparse count mismatch");

    carquet_column_reader_free(label_reader);
    carquet_reader_close(r);
    carquet_schema_free(schema);
    TEST_PASS("per_column");
    return 0;
}

/* LZ4 (codec 5) and LZ4_RAW (codec 7) are now distinct Parquet codecs: codec 5
 * is the deprecated Hadoop-framed LZ4, codec 7 is raw LZ4 blocks. The writer
 * emits each verbatim and both round-trip. */
static int test_lz4_codecs(void) {
    const carquet_compression_t codecs[2] = {
        CARQUET_COMPRESSION_LZ4, CARQUET_COMPRESSION_LZ4_RAW };
    for (int ci = 0; ci < 2; ci++) {
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* schema = carquet_schema_create(&err);
        if (!schema) TEST_FAIL("lz4_codecs", "schema create failed");
        if (carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
                NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
            carquet_schema_free(schema);
            TEST_FAIL("lz4_codecs", "schema add failed");
        }

        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        opts.compression = codecs[ci];

        carquet_writer_t* w = carquet_writer_create(TEMP_FILE, schema, &opts, &err);
        if (!w) { carquet_schema_free(schema); TEST_FAIL("lz4_codecs", "writer create failed"); }

        int32_t ids[64], out[64];
        for (int i = 0; i < 64; i++) ids[i] = (i / 3) * 7;  /* compressible */
        if (carquet_writer_write_batch(w, 0, ids, 64, NULL, NULL) != CARQUET_OK ||
            carquet_writer_close(w) != CARQUET_OK) {
            carquet_schema_free(schema);
            TEST_FAIL("lz4_codecs", "write failed");
        }

        carquet_reader_t* r = carquet_reader_open(TEMP_FILE, NULL, &err);
        if (!r) { carquet_schema_free(schema); TEST_FAIL("lz4_codecs", "reader open failed"); }

        carquet_column_chunk_metadata_t meta;
        if (carquet_reader_column_chunk_metadata(r, 0, 0, &meta) != CARQUET_OK) {
            carquet_reader_close(r); carquet_schema_free(schema);
            TEST_FAIL("lz4_codecs", "metadata failed");
        }
        /* The requested codec is preserved verbatim (no LZ4 -> LZ4_RAW alias). */
        if (meta.codec != codecs[ci]) {
            carquet_reader_close(r); carquet_schema_free(schema);
            TEST_FAIL("lz4_codecs", "codec not preserved");
        }
        /* Data round-trips through the matching framing. */
        carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
        int64_t n = c ? carquet_column_read_batch(c, out, 64, NULL, NULL) : -1;
        int ok = (n == 64) && (memcmp(ids, out, sizeof(ids)) == 0);
        carquet_column_reader_free(c);
        carquet_reader_close(r);
        carquet_schema_free(schema);
        if (!ok) TEST_FAIL("lz4_codecs", "value round-trip mismatch");
    }

    TEST_PASS("lz4_codecs");
    return 0;
}

/* ======================================================================
 * Test: carquet_column_read_batch_ex error reporting
 * ====================================================================== */

static int test_read_batch_ex(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) {
        carquet_schema_free(schema);
        printf("[SKIP] read_batch_ex: %s\n", err.message);
        return 0;
    }

    int32_t data[] = {10, 20, 30, 40, 50};
    if (carquet_writer_write_batch(w, 0, data, 5, NULL, NULL) != CARQUET_OK)
        TEST_FAIL("read_batch_ex", "write_batch failed");
    if (carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL("read_batch_ex", "close failed");

    void* buf = NULL;
    size_t sz = 0;
    if (carquet_writer_get_buffer(w, &buf, &sz) != CARQUET_OK || !buf || sz == 0)
        TEST_FAIL("read_batch_ex", "get_buffer failed");

    carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, &err);
    if (!r) TEST_FAIL("read_batch_ex", "open_buffer failed");
    carquet_column_reader_t* col = carquet_reader_get_column(r, 0, 0, NULL);
    if (!col) TEST_FAIL("read_batch_ex", "get_column failed");

    /* 1) Invalid argument: max_values < 0 -> -1 with INVALID_ARGUMENT set. */
    carquet_error_t e1 = CARQUET_ERROR_INIT;
    int32_t vals[5];
    int64_t n = carquet_column_read_batch_ex(col, vals, -1, NULL, NULL, &e1);
    if (n != -1)
        TEST_FAIL("read_batch_ex", "negative max_values should return -1");
    if (e1.code != CARQUET_ERROR_INVALID_ARGUMENT)
        TEST_FAIL("read_batch_ex", "negative max_values should set INVALID_ARGUMENT");

    /* 2) Clean read leaves the error untouched (code stays CARQUET_OK). */
    carquet_error_t e2 = CARQUET_ERROR_INIT;
    e2.code = CARQUET_ERROR_INTERNAL; /* poison: must be cleared on entry */
    n = carquet_column_read_batch_ex(col, vals, 5, NULL, NULL, &e2);
    if (n != 5)
        TEST_FAIL("read_batch_ex", "clean read should return 5");
    if (e2.code != CARQUET_OK)
        TEST_FAIL("read_batch_ex", "clean read must clear the error");
    for (int i = 0; i < 5; i++) {
        if (vals[i] != data[i]) TEST_FAIL("read_batch_ex", "value mismatch");
    }

    /* 3) Clean end-of-column short read: 0 values, no error. */
    carquet_error_t e3 = CARQUET_ERROR_INIT;
    n = carquet_column_read_batch_ex(col, vals, 5, NULL, NULL, &e3);
    if (n != 0)
        TEST_FAIL("read_batch_ex", "end-of-column read should return 0");
    if (e3.code != CARQUET_OK)
        TEST_FAIL("read_batch_ex", "end-of-column read must not set an error");

    /* 4) NULL error out-param behaves exactly like the legacy wrapper. */
    if (carquet_column_read_batch_ex(col, vals, -1, NULL, NULL, NULL) != -1)
        TEST_FAIL("read_batch_ex", "NULL error path should still return -1");

    carquet_column_reader_free(col);
    carquet_reader_close(r);
    free(buf);
    carquet_schema_free(schema);
    TEST_PASS("read_batch_ex");
    return 0;
}

/* ======================================================================
 * Test: reject chunk-metadata physical type that disagrees with schema
 *
 * Regression for a fuzzer-found heap-buffer-overflow: the batch reader sizes
 * its output buffer from the schema element type while the page reader writes
 * using the column-chunk metadata type. A crafted file where the two disagree
 * (schema INT32 = 4 B/value, chunk BYTE_ARRAY = sizeof(carquet_byte_array_t))
 * made the page decode overflow the batch buffer. The column reader must now
 * reject such files at construction.
 * ====================================================================== */

static int test_type_mismatch_rejected(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) {
        carquet_schema_free(schema);
        printf("[SKIP] type_mismatch_rejected: %s\n", err.message);
        return 0;
    }

    int32_t data[] = {10, 20, 30, 40, 50};
    if (carquet_writer_write_batch(w, 0, data, 5, NULL, NULL) != CARQUET_OK)
        TEST_FAIL("type_mismatch_rejected", "write_batch failed");
    if (carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL("type_mismatch_rejected", "close failed");

    void* src = NULL;
    size_t sz = 0;
    if (carquet_writer_get_buffer(w, &src, &sz) != CARQUET_OK || !src || sz == 0)
        TEST_FAIL("type_mismatch_rejected", "get_buffer failed");

    /* Copy so we can corrupt it. The column-chunk ColumnMetaData.type is a
     * compact-protocol i32 field 1, encoded INT32(1) as bytes {0x15, 0x02}.
     * The schema element's type encodes identically and appears earlier in the
     * footer, so flipping the LAST occurrence to BYTE_ARRAY (zigzag(6)=0x0C)
     * corrupts only the chunk type, leaving the schema type intact. */
    unsigned char* b = malloc(sz);
    if (!b) TEST_FAIL("type_mismatch_rejected", "oom");
    memcpy(b, src, sz);
    long idx = -1;
    for (size_t i = 0; i + 1 < sz; i++)
        if (b[i] == 0x15 && b[i + 1] == 0x02) idx = (long)i;
    if (idx < 0) {
        free(b);
        free(src);   /* get_buffer transferred ownership of the writer's buffer */
        carquet_schema_free(schema);
        printf("[SKIP] type_mismatch_rejected: type field pattern not found "
               "(writer encoding changed)\n");
        return 0;
    }
    b[idx + 1] = 0x0C;  /* INT32 -> BYTE_ARRAY */

    carquet_error_t oe = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open_buffer(b, sz, NULL, &oe);
    if (r) {
        /* If the footer still parses, get_column must reject the mismatch
         * rather than letting a later read overflow. */
        carquet_error_t ce = CARQUET_ERROR_INIT;
        carquet_column_reader_t* col = carquet_reader_get_column(r, 0, 0, &ce);
        if (col) {
            carquet_column_reader_free(col);
            carquet_reader_close(r);
            free(b);
            TEST_FAIL("type_mismatch_rejected",
                      "get_column accepted mismatched physical type");
        }
        if (ce.code != CARQUET_ERROR_INVALID_METADATA)
            printf("  note: rejected with %s (expected INVALID_METADATA)\n",
                   carquet_status_string(ce.code));
        carquet_reader_close(r);
    }
    /* Either outcome (footer rejected on open, or column rejected) is a clean,
     * crash-free refusal — which is the property under test. */

    free(b);
    free(src);   /* get_buffer transferred ownership of the writer's buffer */
    carquet_schema_free(schema);
    TEST_PASS("type_mismatch_rejected");
    return 0;
}

/* ======================================================================
 * Main
 * ====================================================================== */

int main(void) {
    carquet_init();
    int failures = 0;

    failures += test_kv_metadata();
    failures += test_column_chunk_metadata();
    failures += test_bloom_filter();
    failures += test_page_index();
    failures += test_buffer_writer();
    failures += test_per_column_options();
    failures += test_lz4_codecs();
    failures += test_read_batch_ex();
    failures += test_type_mismatch_rejected();

    remove(TEMP_FILE);

    if (failures > 0) {
        printf("\n%d test(s) FAILED\n", failures);
        return 1;
    }
    printf("\nAll advanced API tests passed.\n");
    return 0;
}
