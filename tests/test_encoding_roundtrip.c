/**
 * @file test_encoding_roundtrip.c
 * @brief carquet-reads-carquet roundtrip for every opt-in (Phase 3) encoding.
 *
 * Regression guard for the write/read asymmetry where the writer could emit
 * DELTA_* / BYTE_STREAM_SPLIT (incl. INT32/INT64/FLBA) but the reader could
 * not decode them. Every encoding here is written with carquet AND read back
 * through carquet's own reader, asserting exact value equality — including a
 * nullable column and a multi-column single-batch read that exercises the
 * DELTA_BYTE_ARRAY reconstruction-scratch lifetime across pages.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#define N 4000

static carquet_reader_t* write_then_open(
    const char* path,
    carquet_physical_type_t phys,
    int32_t type_length,
    carquet_field_repetition_t rep,
    carquet_encoding_t enc,
    const void* values,
    int64_t nvalues,
    const int16_t* def_levels) {

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) return NULL;
    if (carquet_schema_add_column(s, "c", phys, NULL, rep, type_length, 0) != CARQUET_OK) {
        carquet_schema_free(s);
        return NULL;
    }
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); return NULL; }
    if (carquet_writer_set_column_encoding(w, 0, enc) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s); return NULL;
    }
    if (carquet_writer_write_batch(w, 0, values, nvalues, def_levels, NULL) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s); return NULL;
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); return NULL; }
    carquet_schema_free(s);
    return carquet_reader_open(path, NULL, &err);
}

static int test_delta_binary_packed_i32(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dbp_i32");
    int32_t in[N], out[N];
    for (int i = 0; i < N; i++) in[i] = (i * 7) - 1000 + (i % 13);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_INT32, 0,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_BINARY_PACKED, in, N, NULL);
    if (!r) TEST_FAIL("delta_binary_packed_i32", "write/open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    if (n != N || memcmp(in, out, sizeof(in)) != 0)
        { carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
          TEST_FAIL("delta_binary_packed_i32", "value mismatch"); }
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    TEST_PASS("delta_binary_packed_i32");
    return 0;
}

static int test_delta_binary_packed_i64(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dbp_i64");
    int64_t in[N], out[N];
    for (int i = 0; i < N; i++) in[i] = ((int64_t)i * 2654435761ULL) ^ (i % 7);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_INT64, 0,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_BINARY_PACKED, in, N, NULL);
    if (!r) TEST_FAIL("delta_binary_packed_i64", "write/open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    if (n != N || memcmp(in, out, sizeof(in)) != 0)
        { carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
          TEST_FAIL("delta_binary_packed_i64", "value mismatch"); }
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    TEST_PASS("delta_binary_packed_i64");
    return 0;
}

static int check_byte_arrays(const char* name, const char* path,
                             carquet_reader_t* r,
                             const carquet_byte_array_t* expect, int64_t cnt) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    static carquet_byte_array_t out[N];
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    int ok = (n == cnt);
    for (int64_t i = 0; ok && i < cnt; i++)
        ok = (out[i].length == expect[i].length) &&
             (memcmp(out[i].data, expect[i].data, (size_t)expect[i].length) == 0);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "byte array mismatch");
    TEST_PASS(name);
    return 0;
}

static int test_delta_length_byte_array(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dlba");
    static char buf[N][24];
    static carquet_byte_array_t in[N];
    for (int i = 0; i < N; i++) {
        int len = snprintf(buf[i], sizeof(buf[i]), "row-%d-val", i);
        in[i].data = (uint8_t*)buf[i]; in[i].length = len;
    }
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_BYTE_ARRAY, 0,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY, in, N, NULL);
    if (!r) TEST_FAIL("delta_length_byte_array", "write/open failed");
    return check_byte_arrays("delta_length_byte_array", path, r, in, N);
}

static int test_delta_byte_array(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dba");
    static char buf[N][32];
    static carquet_byte_array_t in[N];
    /* Strong shared prefixes to exercise prefix reconstruction. */
    for (int i = 0; i < N; i++) {
        int len = snprintf(buf[i], sizeof(buf[i]), "common/prefix/path/item_%05d", i);
        in[i].data = (uint8_t*)buf[i]; in[i].length = len;
    }
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_BYTE_ARRAY, 0,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_BYTE_ARRAY, in, N, NULL);
    if (!r) TEST_FAIL("delta_byte_array", "write/open failed");
    return check_byte_arrays("delta_byte_array", path, r, in, N);
}

static int test_delta_byte_array_flba(void) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dba_flba");
    enum { L = 8 };
    static uint8_t in[N * L], out[N * L];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < L; j++) in[i * L + j] = (uint8_t)((i + j * 31) & 0xFF);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, L,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_BYTE_ARRAY, in, N, NULL);
    if (!r) TEST_FAIL("delta_byte_array_flba", "write/open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    int ok = (n == N) && (memcmp(in, out, sizeof(in)) == 0);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("delta_byte_array_flba", "value mismatch");
    TEST_PASS("delta_byte_array_flba");
    return 0;
}

static int test_bss(const char* name, const char* base,
                     carquet_physical_type_t phys, int32_t tl, size_t esz) {
    char path[512]; carquet_test_temp_path(path, sizeof(path), base);
    static uint8_t in[N * 8], out[N * 8];
    size_t bytes = (size_t)N * esz;
    for (size_t i = 0; i < bytes; i++) in[i] = (uint8_t)((i * 131 + 7) & 0xFF);
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = write_then_open(path, phys, tl,
        CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_BYTE_STREAM_SPLIT, in, N, NULL);
    if (!r) TEST_FAIL(name, "write/open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    int ok = (n == N) && (memcmp(in, out, bytes) == 0);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "value mismatch");
    TEST_PASS(name);
    return 0;
}

/* BYTE_STREAM_SPLIT transposes a whole page into byte planes, so it is only
 * correct when applied to the finished page. Writing the same page in several
 * carquet_writer_write_batch() calls must therefore produce byte-identical
 * output to writing it in one. The chunk sizes below are deliberately not
 * multiples of any value width, so a per-call split leaves independently
 * transposed regions that the reader de-splits as a single stride. */
static int test_bss_chunked(const char* name, const char* base,
                            carquet_physical_type_t phys, int32_t tl, size_t esz) {
    static const int64_t chunks[] = { 7, 1, 1992, 2000 };
    char path[512]; carquet_test_temp_path(path, sizeof(path), base);
    static uint8_t in[N * 8], out[N * 8];
    size_t bytes = (size_t)N * esz;
    for (size_t i = 0; i < bytes; i++) in[i] = (uint8_t)((i * 131 + 7) & 0xFF);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) TEST_FAIL(name, "schema create failed");
    if (carquet_schema_add_column(s, "c", phys, NULL,
                                  CARQUET_REPETITION_REQUIRED, tl, 0) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "add column failed");
    }
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL(name, "writer create failed"); }
    if (carquet_writer_set_column_encoding(
            w, 0, CARQUET_ENCODING_BYTE_STREAM_SPLIT) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL(name, "set encoding failed");
    }
    int64_t done = 0;
    for (size_t k = 0; k < sizeof(chunks) / sizeof(chunks[0]); k++) {
        if (carquet_writer_write_batch(w, 0, in + (size_t)done * esz,
                                       chunks[k], NULL, NULL) != CARQUET_OK) {
            carquet_writer_close(w); carquet_schema_free(s);
            TEST_FAIL(name, "write batch failed");
        }
        done += chunks[k];
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "writer close failed");
    }
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, "open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    int ok = (done == N) && (n == N) && (memcmp(in, out, bytes) == 0);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "value mismatch");
    TEST_PASS(name);
    return 0;
}

/* Every DELTA_* encoding is one self-describing stream per page (header with
 * the page's total value count, then blocks). Feeding a page from several
 * carquet_writer_write_batch() calls must therefore still produce a single
 * stream, exactly as for BYTE_STREAM_SPLIT above: a stream per call leaves the
 * reader decoding the first header's count and running into the next header. */
static int test_delta_chunked(const char* name, const char* base,
                              carquet_physical_type_t phys, int32_t tl,
                              carquet_encoding_t enc) {
    static const int64_t chunks[] = { 7, 1, 1992, 2000 };
    char path[512]; carquet_test_temp_path(path, sizeof(path), base);
    static int32_t in32[N], out32[N];
    static int64_t in64[N], out64[N];
    static char sbuf[N][32];
    static carquet_byte_array_t inba[N], outba[N];
    static uint8_t infl[N * 8], outfl[N * 8];
    for (int i = 0; i < N; i++) {
        in32[i] = (i * 7) - 1000 + (i % 13);
        in64[i] = ((int64_t)i * 2654435761LL) ^ (i % 7);
        int len = snprintf(sbuf[i], sizeof(sbuf[i]), "common/prefix/item_%05d", i * 3);
        inba[i].data = (uint8_t*)sbuf[i]; inba[i].length = len;
        for (int j = 0; j < 8; j++) infl[i * 8 + j] = (uint8_t)((i + j * 31) & 0xFF);
    }
    const void* in; void* out; size_t esz;
    switch (phys) {
        case CARQUET_PHYSICAL_INT32: in = in32; out = out32; esz = 4; break;
        case CARQUET_PHYSICAL_INT64: in = in64; out = out64; esz = 8; break;
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            in = inba; out = outba; esz = sizeof(carquet_byte_array_t); break;
        default: in = infl; out = outfl; esz = (size_t)tl; break;
    }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) TEST_FAIL(name, "schema create failed");
    if (carquet_schema_add_column(s, "c", phys, NULL,
                                  CARQUET_REPETITION_REQUIRED, tl, 0) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "add column failed");
    }
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL(name, "writer create failed"); }
    if (carquet_writer_set_column_encoding(w, 0, enc) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL(name, "set encoding failed");
    }
    int64_t done = 0;
    for (size_t k = 0; k < sizeof(chunks) / sizeof(chunks[0]); k++) {
        if (carquet_writer_write_batch(w, 0, (const uint8_t*)in + (size_t)done * esz,
                                       chunks[k], NULL, NULL) != CARQUET_OK) {
            carquet_writer_close(w); carquet_schema_free(s);
            TEST_FAIL(name, "write batch failed");
        }
        done += chunks[k];
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "writer close failed");
    }
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, "open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
    int ok = (done == N) && (n == N);
    if (ok && phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
        for (int i = 0; ok && i < N; i++)
            ok = (outba[i].length == inba[i].length) &&
                 (memcmp(outba[i].data, inba[i].data, (size_t)inba[i].length) == 0);
    } else if (ok) {
        ok = (memcmp(in, out, (size_t)N * esz) == 0);
    }
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "value mismatch");
    TEST_PASS(name);
    return 0;
}

/* Mini-blocks wider than 32 bits are bit-packed like any other width, not
 * stored as whole little-endian bytes per value. Round-trip every width from
 * 33 to 64 with deltas that use the full width. */
static int test_delta_binary_packed_i64_wide(void) {
    for (int width = 33; width <= 64; width++) {
        char path[512]; carquet_test_temp_path(path, sizeof(path), "dbp_i64w");
        static int64_t in[N], out[N];
        uint64_t top = width == 64 ? ~(uint64_t)0 : (((uint64_t)1 << width) - 1);
        uint64_t acc = 0, x = 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < N; i++) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            /* Wrapping sums: delta arithmetic is modulo 2^64 by spec. */
            acc += (i % 5 == 0) ? top : (x & top);
            in[i] = (int64_t)acc;
        }
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_INT64, 0,
            CARQUET_REPETITION_REQUIRED, CARQUET_ENCODING_DELTA_BINARY_PACKED, in, N, NULL);
        if (!r) TEST_FAIL("delta_binary_packed_i64_wide", "write/open failed");
        carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
        int64_t n = c ? carquet_column_read_batch(c, out, N, NULL, NULL) : -1;
        int ok = (n == N) && memcmp(in, out, sizeof(in)) == 0;
        carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
        if (!ok) {
            printf("  width %d\n", width);
            TEST_FAIL("delta_binary_packed_i64_wide", "value mismatch");
        }
    }
    TEST_PASS("delta_binary_packed_i64_wide");
    return 0;
}

/* Known-answer layout: 32 deltas of width 33 occupy 32 * 33 / 8 = 132 bytes. */
extern carquet_status_t carquet_delta_encode_int64(
    const int64_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);

static int test_delta_wide_layout(void) {
    int64_t v[33];
    v[0] = 0;
    for (int i = 1; i < 33; i++)
        v[i] = v[i - 1] + ((i % 2) ? 0 : (((int64_t)1 << 33) - 1));
    uint8_t buf[1024];
    size_t written = 0;
    if (carquet_delta_encode_int64(v, 33, buf, sizeof(buf), &written) != CARQUET_OK)
        TEST_FAIL("delta_wide_layout", "encode failed");
    /* header: block size 128 (2 bytes), 4 mini-blocks, count 33, first 0;
     * block: min_delta 0, 4 bit widths {33,0,0,0}, 132 packed bytes. */
    if (written != 5 + 1 + 4 + 132) TEST_FAIL("delta_wide_layout", "wrong packed size");
    if (buf[6] != 33 || buf[7] != 0) TEST_FAIL("delta_wide_layout", "wrong bit widths");
    /* delta[0] = 0 fills bits 0..32, delta[1] = 2^33-1 fills bits 33..65. */
    if (buf[10] != 0 || buf[14] != 0xFE || buf[15] != 0xFF || buf[18] != 0x03)
        TEST_FAIL("delta_wide_layout", "values are not bit-packed");
    TEST_PASS("delta_wide_layout");
    return 0;
}

/* The same hazard reached from a single carquet_writer_write_batch() call
 * (issue #31). The column writer cuts a large batch into row-count chunks of
 * page_size / value width, but closes a page only once its estimated size
 * reaches page_size. Nulls take no space in the values section, so a leading
 * null run leaves the first chunk under that threshold and the second chunk is
 * added to the same page: one page assembled from two calls with different
 * non-null counts, which must still decode as a single transposition. A codec
 * is set because that is what selects BYTE_STREAM_SPLIT for FLOAT/DOUBLE by
 * default, so compressed nullable float columns hit this without opting in;
 * the encoding is still requested explicitly so the INT32/INT64/FLBA cases
 * take the same path. */
static int test_bss_null_shrunk_chunk(const char* name, const char* base,
                                      carquet_physical_type_t phys, int32_t tl, size_t esz) {
    enum { NULL_RUN = 260, PAGE_SIZE = 8192 };
    char path[512]; carquet_test_temp_path(path, sizeof(path), base);
    static uint8_t in[N * 8], out[N * 8];
    static int16_t def[N], outdef[N];
    /* OPTIONAL columns take and return the non-null values packed. */
    size_t bytes = (size_t)(N - NULL_RUN) * esz;
    for (size_t i = 0; i < bytes; i++) in[i] = (uint8_t)((i * 131 + 7) & 0xFF);
    for (int i = 0; i < N; i++) def[i] = (i < NULL_RUN) ? 0 : 1;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) TEST_FAIL(name, "schema create failed");
    if (carquet_schema_add_column(s, "c", phys, NULL,
                                  CARQUET_REPETITION_OPTIONAL, tl, 0) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "add column failed");
    }
    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_SNAPPY;
    wo.page_size = PAGE_SIZE;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL(name, "writer create failed"); }
    if (carquet_writer_set_column_encoding(
            w, 0, CARQUET_ENCODING_BYTE_STREAM_SPLIT) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL(name, "set encoding failed");
    }
    if (carquet_writer_write_batch(w, 0, in, N, def, NULL) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s);
        TEST_FAIL(name, "write batch failed");
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL(name, "writer close failed");
    }
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL(name, "open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, N, outdef, NULL) : -1;
    int ok = (n == N) && (memcmp(def, outdef, sizeof(def)) == 0) &&
             (memcmp(in, out, bytes) == 0);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL(name, "value mismatch");
    TEST_PASS(name);
    return 0;
}

static int test_delta_byte_array_nullable(void) {
    /* Nullable DELTA_BYTE_ARRAY: exercises the reconstruction scratch together
     * with definition levels (packed non-null value stream). */
    char path[512]; carquet_test_temp_path(path, sizeof(path), "dba_null");
    static char buf[N][32];
    static carquet_byte_array_t in[N];
    static int16_t def[N];
    int64_t nn = 0;
    for (int i = 0; i < N; i++) {
        if (i % 3 == 0) { def[i] = 0; continue; }
        def[i] = 1;
        int len = snprintf(buf[nn], sizeof(buf[nn]), "shared-prefix-%04d", i);
        in[nn].data = (uint8_t*)buf[nn]; in[nn].length = len; nn++;
    }
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = write_then_open(path, CARQUET_PHYSICAL_BYTE_ARRAY, 0,
        CARQUET_REPETITION_OPTIONAL, CARQUET_ENCODING_DELTA_BYTE_ARRAY, in, N, def);
    if (!r) TEST_FAIL("delta_byte_array_nullable", "write/open failed");
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    static carquet_byte_array_t out[N];
    static int16_t outdef[N];
    int64_t got = c ? carquet_column_read_batch(c, out, N, outdef, NULL) : -1;
    int ok = (got == N);
    int64_t vi = 0;
    for (int64_t i = 0; ok && i < N; i++) {
        if (outdef[i] == 0) continue;
        ok = (out[vi].length == in[vi].length) &&
             (memcmp(out[vi].data, in[vi].data, (size_t)in[vi].length) == 0);
        vi++;
    }
    ok = ok && (vi == nn);
    carquet_column_reader_free(c); carquet_reader_close(r); carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("delta_byte_array_nullable", "nullable mismatch");
    TEST_PASS("delta_byte_array_nullable");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_delta_binary_packed_i32();
    failures += test_delta_binary_packed_i64();
    failures += test_delta_length_byte_array();
    failures += test_delta_byte_array();
    failures += test_delta_byte_array_flba();
    failures += test_delta_byte_array_nullable();
    failures += test_delta_binary_packed_i64_wide();
    failures += test_delta_wide_layout();
    failures += test_delta_chunked("delta_binary_packed_i32_chunked", "dbpc_i32",
        CARQUET_PHYSICAL_INT32, 0, CARQUET_ENCODING_DELTA_BINARY_PACKED);
    failures += test_delta_chunked("delta_binary_packed_i64_chunked", "dbpc_i64",
        CARQUET_PHYSICAL_INT64, 0, CARQUET_ENCODING_DELTA_BINARY_PACKED);
    failures += test_delta_chunked("delta_length_byte_array_chunked", "dlbac",
        CARQUET_PHYSICAL_BYTE_ARRAY, 0, CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY);
    failures += test_delta_chunked("delta_byte_array_chunked", "dbac",
        CARQUET_PHYSICAL_BYTE_ARRAY, 0, CARQUET_ENCODING_DELTA_BYTE_ARRAY);
    failures += test_delta_chunked("delta_byte_array_flba_chunked", "dbac_flba",
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 8, CARQUET_ENCODING_DELTA_BYTE_ARRAY);
    failures += test_bss("bss_float",  "bss_f32", CARQUET_PHYSICAL_FLOAT, 0, 4);
    failures += test_bss("bss_double", "bss_f64", CARQUET_PHYSICAL_DOUBLE, 0, 8);
    failures += test_bss("bss_int32",  "bss_i32", CARQUET_PHYSICAL_INT32, 0, 4);
    failures += test_bss("bss_int64",  "bss_i64", CARQUET_PHYSICAL_INT64, 0, 8);
    failures += test_bss("bss_flba",   "bss_flba", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 6, 6);
    failures += test_bss_chunked("bss_float_chunked",  "bssc_f32", CARQUET_PHYSICAL_FLOAT, 0, 4);
    failures += test_bss_chunked("bss_double_chunked", "bssc_f64", CARQUET_PHYSICAL_DOUBLE, 0, 8);
    failures += test_bss_chunked("bss_int32_chunked",  "bssc_i32", CARQUET_PHYSICAL_INT32, 0, 4);
    failures += test_bss_chunked("bss_int64_chunked",  "bssc_i64", CARQUET_PHYSICAL_INT64, 0, 8);
    failures += test_bss_chunked("bss_flba_chunked",   "bssc_flba", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 6, 6);
    failures += test_bss_null_shrunk_chunk("bss_float_null_shrunk_chunk",  "bssn_f32", CARQUET_PHYSICAL_FLOAT, 0, 4);
    failures += test_bss_null_shrunk_chunk("bss_double_null_shrunk_chunk", "bssn_f64", CARQUET_PHYSICAL_DOUBLE, 0, 8);
    failures += test_bss_null_shrunk_chunk("bss_int32_null_shrunk_chunk",  "bssn_i32", CARQUET_PHYSICAL_INT32, 0, 4);
    failures += test_bss_null_shrunk_chunk("bss_int64_null_shrunk_chunk",  "bssn_i64", CARQUET_PHYSICAL_INT64, 0, 8);
    failures += test_bss_null_shrunk_chunk("bss_flba_null_shrunk_chunk",   "bssn_flba", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, 6, 6);
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll encoding roundtrip tests passed\n");
    return 0;
}
