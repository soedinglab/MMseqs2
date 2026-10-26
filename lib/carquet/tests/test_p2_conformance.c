/**
 * @file test_p2_conformance.c
 * @brief Tests for the P2 Parquet-conformance fixes:
 *   - FIXED_LEN_BYTE_ARRAY dictionary encoding
 *   - distinct_count in statistics (write, dictionary columns)
 *   - column_orders preserved on read
 *   - SizeStatistics (Parquet 2.9): unencoded_byte_array_data_bytes +
 *     definition/repetition level histograms
 *   - OffsetIndex field 2 = unencoded_byte_array_data_bytes (list<i64>)
 *   - RLE as a value encoding for BOOLEAN
 *   - Hadoop-framed LZ4 (codec 5) round-trip
 *   - FLOAT/DOUBLE min/max statistics skip NaN without losing real extremes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <carquet/carquet.h>
#include "reader/reader_internal.h"
#include "thrift/thrift_decode.h"
#include "test_helpers.h"

/* ---- FLBA dictionary + distinct_count ------------------------------------ */

static int test_flba_dictionary(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_flba_dict");
    const int L = 5, N = 300;
    unsigned char* vals = malloc((size_t)N * L);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < L; j++)
            vals[i * L + j] = (unsigned char)((i % 7) * 10 + j);  /* 7 distinct */

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (carquet_schema_add_column(s, "c", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
            NULL, CARQUET_REPETITION_REQUIRED, L, 0) != CARQUET_OK) {
        free(vals); carquet_schema_free(s);
        TEST_FAIL("flba_dictionary", "add_column failed");
    }
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY);
    if (carquet_writer_write_batch(w, 0, vals, N, NULL, NULL) != CARQUET_OK) {
        free(vals); carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL("flba_dictionary", "write_batch failed");
    }
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { free(vals); carquet_test_cleanup(path); TEST_FAIL("flba_dictionary", "open failed"); }

    /* Values round-trip. */
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    unsigned char* out = malloc((size_t)N * L);
    int64_t n = carquet_column_read_batch(c, out, N, NULL, NULL);
    int values_ok = (n == N) && (memcmp(vals, out, (size_t)N * L) == 0);

    /* A dictionary page was emitted (dict actually used, not PLAIN fallback). */
    const parquet_file_metadata_t* meta = &r->metadata;
    parquet_column_chunk_t* col = &meta->row_groups[0].columns[0];
    int dict_used = col->metadata.has_dictionary_page_offset;

    /* distinct_count computed from the dictionary. */
    carquet_column_statistics_t st;
    int have_stats = (carquet_reader_column_statistics(r, 0, 0, &st) == CARQUET_OK);
    int distinct_ok = have_stats && st.has_distinct_count && st.distinct_count == 7;

    free(vals); free(out);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);

    if (!values_ok) TEST_FAIL("flba_dictionary", "value mismatch");
    if (!dict_used) TEST_FAIL("flba_dictionary", "no dictionary page emitted");
    if (!distinct_ok) TEST_FAIL("flba_dictionary", "distinct_count != 7");
    TEST_PASS("flba_dictionary");
    return 0;
}

/* ---- column_orders preserved on read ------------------------------------- */

static int test_column_orders_preserved(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_colorders");
    int32_t a[16], b[16];
    for (int i = 0; i < 16; i++) { a[i] = i; b[i] = i * 2; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "a", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    carquet_writer_write_batch(w, 0, a, 16, NULL, NULL);
    carquet_writer_write_batch(w, 1, b, 16, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("column_orders_preserved", "open failed"); }
    const parquet_file_metadata_t* meta = &r->metadata;

    int ok = (meta->num_column_orders == 2) && (meta->column_order_types != NULL) &&
             (meta->column_order_types[0] == 1) && (meta->column_order_types[1] == 1);

    carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("column_orders_preserved", "column_order_types not preserved");
    TEST_PASS("column_orders_preserved");
    return 0;
}

/* ---- SizeStatistics (unencoded bytes + def-level histogram) --------------- */

static int test_size_statistics(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_sizestats");

    /* OPTIONAL byte-array column: rows [ "aa", NULL, "bbbb", "c", NULL ]. */
    carquet_byte_array_t vals[3];
    vals[0].data = (uint8_t*)"aa";   vals[0].length = 2;
    vals[1].data = (uint8_t*)"bbbb"; vals[1].length = 4;
    vals[2].data = (uint8_t*)"c";    vals[2].length = 1;
    int16_t def[5] = {1, 0, 1, 1, 0};  /* max_def_level == 1; 2 nulls, 3 present */
    int64_t total_value_bytes = 2 + 4 + 1;  /* 7 */

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "c", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    /* PLAIN so unencoded bytes come from the normal add_values path. */
    carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN);
    if (carquet_writer_write_batch(w, 0, vals, 5, def, NULL) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL("size_statistics", "write_batch failed");
    }
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("size_statistics", "open failed"); }
    const parquet_file_metadata_t* meta = &r->metadata;
    const parquet_column_metadata_t* cm = &meta->row_groups[0].columns[0].metadata;

    int ok = cm->has_size_statistics;
    const parquet_size_statistics_t* ss = &cm->size_statistics;
    ok = ok && ss->has_unencoded_byte_array_data_bytes &&
         ss->unencoded_byte_array_data_bytes == total_value_bytes;
    /* definition_level_histogram == [nulls=2, present=3]. */
    ok = ok && ss->definition_level_histogram_len == 2 &&
         ss->definition_level_histogram[0] == 2 &&
         ss->definition_level_histogram[1] == 3;
    /* A flat column has no meaningful repetition histogram: not emitted. */
    ok = ok && ss->repetition_level_histogram_len == 0;

    carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("size_statistics", "size statistics mismatch");
    TEST_PASS("size_statistics");
    return 0;
}

/* ---- OffsetIndex field 2 (unencoded_byte_array_data_bytes, list<i64>) ----- */

/* Read the raw OffsetIndex bytes and pull out field 2's total (sum over
 * pages). Field 2 is a list<i64>; we sum every element. Returns -1 on error. */
static int64_t read_offset_index_field2_sum(const char* path,
                                             int64_t off, int32_t len) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;
    uint8_t* buf = malloc((size_t)len);
    if (fseek(f, (long)off, SEEK_SET) != 0 ||
        fread(buf, 1, (size_t)len, f) != (size_t)len) {
        free(buf); fclose(f); return -1;
    }
    fclose(f);

    int64_t sum = 0;
    int found = 0;
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, buf, (size_t)len);
    thrift_read_struct_begin(&dec);
    thrift_type_t ft; int16_t fid;
    while (thrift_read_field_begin(&dec, &ft, &fid)) {
        if (fid == 2 && ft == THRIFT_TYPE_LIST) {
            thrift_type_t et; int32_t count;
            thrift_read_list_begin(&dec, &et, &count);
            found = 1;
            for (int32_t i = 0; i < count; i++) sum += thrift_read_i64(&dec);
        } else {
            thrift_skip(&dec, ft);
        }
    }
    thrift_read_struct_end(&dec);
    free(buf);
    if (thrift_decoder_has_error(&dec) || !found) return -1;
    return sum;
}

static int test_offset_index_field2(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_offidx");
    const int N = 400;
    carquet_byte_array_t* vals = malloc((size_t)N * sizeof(*vals));
    static char storage[400][8];
    int64_t expect_bytes = 0;
    for (int i = 0; i < N; i++) {
        int len = snprintf(storage[i], sizeof(storage[i]), "v%d", i);
        vals[i].data = (uint8_t*)storage[i];
        vals[i].length = len;
        expect_bytes += len;
    }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "c", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.write_page_index = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN);
    carquet_writer_write_batch(w, 0, vals, N, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);
    free(vals);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("offset_index_field2", "open failed"); }
    const parquet_file_metadata_t* meta = &r->metadata;
    parquet_column_chunk_t* col = &meta->row_groups[0].columns[0];
    int have = col->has_offset_index_offset && col->has_offset_index_length;
    int64_t off = col->offset_index_offset;
    int32_t olen = col->offset_index_length;
    carquet_reader_close(r);

    int64_t sum = have ? read_offset_index_field2_sum(path, off, olen) : -1;
    carquet_test_cleanup(path);

    if (!have) TEST_FAIL("offset_index_field2", "no offset index written");
    if (sum != expect_bytes) TEST_FAIL("offset_index_field2", "field 2 sum mismatch");
    TEST_PASS("offset_index_field2");
    return 0;
}

/* ---- RLE BOOLEAN value encoding round-trip ------------------------------- */

static int test_rle_boolean(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_rle_bool");
    const int N = 500;
    uint8_t* in = malloc((size_t)N);
    uint8_t* out = malloc((size_t)N);
    for (int i = 0; i < N; i++) in[i] = (uint8_t)(((i / 7) % 2) == 0);  /* runs */

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "b", CARQUET_PHYSICAL_BOOLEAN, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE) != CARQUET_OK) {
        free(in); free(out); carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL("rle_boolean", "set RLE encoding failed");
    }
    if (carquet_writer_write_batch(w, 0, in, N, NULL, NULL) != CARQUET_OK) {
        free(in); free(out); carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL("rle_boolean", "write_batch failed");
    }
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { free(in); free(out); carquet_test_cleanup(path); TEST_FAIL("rle_boolean", "open failed"); }
    /* The chunk must advertise RLE, not PLAIN, for the data encoding. */
    const parquet_file_metadata_t* meta = &r->metadata;
    const parquet_column_metadata_t* cm = &meta->row_groups[0].columns[0].metadata;
    int advertises_rle = 0;
    for (int32_t i = 0; i < cm->num_encodings; i++)
        if (cm->encodings[i] == CARQUET_ENCODING_RLE) advertises_rle = 1;

    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = carquet_column_read_batch(c, out, N, NULL, NULL);
    int ok = (n == N) && (memcmp(in, out, (size_t)N) == 0);

    free(in); free(out);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!advertises_rle) TEST_FAIL("rle_boolean", "chunk did not advertise RLE encoding");
    if (!ok) TEST_FAIL("rle_boolean", "boolean value mismatch");
    TEST_PASS("rle_boolean");
    return 0;
}

/* ---- Hadoop-framed LZ4 (codec 5) round-trip ------------------------------ */

static int test_hadoop_lz4(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "p2_lz4_hadoop");
    const int N = 2000;
    int32_t* in = malloc((size_t)N * sizeof(int32_t));
    int32_t* out = malloc((size_t)N * sizeof(int32_t));
    for (int i = 0; i < N; i++) in[i] = (i / 5) * 3;  /* compressible */

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT32, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_LZ4;  /* codec 5, Hadoop-framed */
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (carquet_writer_write_batch(w, 0, in, N, NULL, NULL) != CARQUET_OK) {
        free(in); free(out); carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL("hadoop_lz4", "write_batch failed");
    }
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { free(in); free(out); carquet_test_cleanup(path); TEST_FAIL("hadoop_lz4", "open failed"); }
    const parquet_file_metadata_t* meta = &r->metadata;
    int codec_ok = meta->row_groups[0].columns[0].metadata.codec == CARQUET_COMPRESSION_LZ4;
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = carquet_column_read_batch(c, out, N, NULL, NULL);
    int ok = (n == N) && (memcmp(in, out, (size_t)N * sizeof(int32_t)) == 0);

    free(in); free(out);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!codec_ok) TEST_FAIL("hadoop_lz4", "codec not LZ4");
    if (!ok) TEST_FAIL("hadoop_lz4", "lz4 value mismatch");
    TEST_PASS("hadoop_lz4");
    return 0;
}

/* ---- NaN-skipping FLOAT/DOUBLE min/max statistics ------------------------- */

/* NaN must never be a min/max, and must not hide the real extreme either. The
 * NEON double reduction let a NaN lane swallow the 9 below and then dropped
 * the NaN, recording max 4 (#35). PLAIN takes the copy+minmax path, the other
 * encodings the plain minmax path, so cover all three. */
static int test_nan_minmax_stats(void) {
    static const double dvals[8] = {1, 9, 2, NAN, 4, 3, 2, 1};
    static const float fvals[8] = {1, 9, 2, NAN, 4, 3, 2, 1};
    static const carquet_encoding_t encs[3] = {
        CARQUET_ENCODING_PLAIN, CARQUET_ENCODING_BYTE_STREAM_SPLIT,
        CARQUET_ENCODING_RLE_DICTIONARY,
    };

    for (int e = 0; e < 3; e++) {
        char path[512];
        carquet_test_temp_path(path, sizeof(path), "p2_nan_minmax");
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* s = carquet_schema_create(&err);
        carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_DOUBLE, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_schema_add_column(s, "f", CARQUET_PHYSICAL_FLOAT, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0);
        carquet_writer_options_t wo; carquet_writer_options_init(&wo);
        wo.write_statistics = true;
        carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
        if (!w) { carquet_schema_free(s); TEST_FAIL("nan_minmax_stats", "create failed"); }
        carquet_writer_set_column_encoding(w, 0, encs[e]);
        carquet_writer_set_column_encoding(w, 1, encs[e]);
        if (carquet_writer_write_batch(w, 0, dvals, 8, NULL, NULL) != CARQUET_OK ||
            carquet_writer_write_batch(w, 1, fvals, 8, NULL, NULL) != CARQUET_OK) {
            carquet_writer_close(w); carquet_schema_free(s);
            TEST_FAIL("nan_minmax_stats", "write_batch failed");
        }
        if (carquet_writer_close(w) != CARQUET_OK) {
            carquet_schema_free(s); TEST_FAIL("nan_minmax_stats", "close failed");
        }
        carquet_schema_free(s);

        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        if (!r) { carquet_test_cleanup(path); TEST_FAIL("nan_minmax_stats", "open failed"); }
        carquet_column_statistics_t ds, fs;
        int ok = carquet_reader_column_statistics(r, 0, 0, &ds) == CARQUET_OK &&
                 carquet_reader_column_statistics(r, 0, 1, &fs) == CARQUET_OK &&
                 ds.has_min_max && fs.has_min_max &&
                 ds.min_value_size == 8 && fs.min_value_size == 4;
        double dmin = 0, dmax = 0;
        float fmin = 0, fmax = 0;
        if (ok) {
            memcpy(&dmin, ds.min_value, 8); memcpy(&dmax, ds.max_value, 8);
            memcpy(&fmin, fs.min_value, 4); memcpy(&fmax, fs.max_value, 4);
        }
        carquet_reader_close(r); carquet_test_cleanup(path);
        if (!ok) TEST_FAIL("nan_minmax_stats", "no statistics");
        if (dmin != 1.0 || dmax != 9.0) {
            printf("  encoding %d: double min=%g max=%g\n", (int)encs[e], dmin, dmax);
            TEST_FAIL("nan_minmax_stats", "DOUBLE min/max wrong next to a NaN");
        }
        if (fmin != 1.0f || fmax != 9.0f) {
            printf("  encoding %d: float min=%g max=%g\n", (int)encs[e], fmin, fmax);
            TEST_FAIL("nan_minmax_stats", "FLOAT min/max wrong next to a NaN");
        }
    }
    TEST_PASS("nan_minmax_stats");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_flba_dictionary();
    failures += test_column_orders_preserved();
    failures += test_size_statistics();
    failures += test_offset_index_field2();
    failures += test_rle_boolean();
    failures += test_hadoop_lz4();
    failures += test_nan_minmax_stats();
    if (failures == 0) printf("\nAll P2 conformance tests passed\n");
    return failures ? 1 : 0;
}
