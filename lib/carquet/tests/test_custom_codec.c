/**
 * @file test_custom_codec.c
 * @brief Tests for the pluggable custom-compression-codec registration API.
 *
 * carquet_register_codec() lets an application override or supply a codec for
 * a compression slot. These tests cover both halves of the contract:
 *   1. Argument validation — rejecting UNCOMPRESSED, out-of-range slots, and
 *      NULL function pointers (documented CARQUET_ERROR_INVALID_ARGUMENT).
 *   2. Real end-to-end use — registering a reversible XOR codec on the ZSTD
 *      slot, writing a file through the public writer with that compression,
 *      reading it back through the public reader, and asserting BOTH that the
 *      values roundtrip AND that our compress/decompress callbacks actually
 *      ran (a builtin-ZSTD fallback would corrupt XOR-obfuscated bytes, so a
 *      clean roundtrip is proof the custom path is wired end to end).
 *
 * This exercises the writer/reader compression dispatch a feature author is
 * most likely to break, and has no other coverage.
 */

#include <stdint.h>
#include <string.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

/* ---- a reversible XOR "codec" with call accounting --------------------- */

typedef struct {
    int compress_calls;
    int decompress_calls;
    int bound_calls;
    uint8_t key;
} xor_state_t;

static carquet_status_t xor_compress(const uint8_t* src, size_t n,
                                     uint8_t* dst, size_t cap, size_t* out,
                                     int32_t level, void* ud) {
    (void)level;
    xor_state_t* st = (xor_state_t*)ud;
    st->compress_calls++;
    if (cap < n) return CARQUET_ERROR_INVALID_ARGUMENT;
    for (size_t i = 0; i < n; i++) dst[i] = src[i] ^ st->key;
    *out = n;
    return CARQUET_OK;
}

static carquet_status_t xor_decompress(const uint8_t* src, size_t n,
                                       uint8_t* dst, size_t cap, size_t* out,
                                       void* ud) {
    xor_state_t* st = (xor_state_t*)ud;
    st->decompress_calls++;
    if (cap != n) return CARQUET_ERROR_INVALID_ARGUMENT;  /* exact declared size */
    for (size_t i = 0; i < n; i++) dst[i] = src[i] ^ st->key;
    *out = n;
    return CARQUET_OK;
}

static size_t xor_bound(size_t n, void* ud) {
    xor_state_t* st = (xor_state_t*)ud;
    st->bound_calls++;
    return n;  /* identity length */
}

/* ---- registration argument validation ---------------------------------- */

static int test_reject_uncompressed_slot(void) {
    xor_state_t st = {0, 0, 0, 0x5A};
    carquet_custom_codec_t codec = {xor_compress, xor_decompress, xor_bound, &st};
    carquet_status_t rc = carquet_register_codec(CARQUET_COMPRESSION_UNCOMPRESSED, &codec);
    if (rc != CARQUET_ERROR_INVALID_ARGUMENT)
        TEST_FAIL("reject_uncompressed_slot", "UNCOMPRESSED slot was not rejected");
    TEST_PASS("reject_uncompressed_slot");
    return 0;
}

static int test_reject_out_of_range(void) {
    xor_state_t st = {0, 0, 0, 0x5A};
    carquet_custom_codec_t codec = {xor_compress, xor_decompress, xor_bound, &st};
    carquet_status_t rc = carquet_register_codec((carquet_compression_t)9999, &codec);
    if (rc != CARQUET_ERROR_INVALID_ARGUMENT)
        TEST_FAIL("reject_out_of_range", "out-of-range slot was not rejected");
    TEST_PASS("reject_out_of_range");
    return 0;
}

static int test_reject_null_function(void) {
    xor_state_t st = {0, 0, 0, 0x5A};
    /* Each of the three required function pointers is individually mandatory. */
    carquet_custom_codec_t missing_compress = {NULL, xor_decompress, xor_bound, &st};
    carquet_custom_codec_t missing_decomp   = {xor_compress, NULL, xor_bound, &st};
    carquet_custom_codec_t missing_bound    = {xor_compress, xor_decompress, NULL, &st};
    if (carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &missing_compress) != CARQUET_ERROR_INVALID_ARGUMENT ||
        carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &missing_decomp) != CARQUET_ERROR_INVALID_ARGUMENT ||
        carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &missing_bound) != CARQUET_ERROR_INVALID_ARGUMENT)
        TEST_FAIL("reject_null_function", "NULL function pointer was not rejected");
    TEST_PASS("reject_null_function");
    return 0;
}

static int test_register_and_unregister(void) {
    xor_state_t st = {0, 0, 0, 0x5A};
    carquet_custom_codec_t codec = {xor_compress, xor_decompress, xor_bound, &st};
    if (carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &codec) != CARQUET_OK)
        TEST_FAIL("register_and_unregister", "valid registration failed");
    /* Passing NULL clears the slot and restores the builtin. */
    if (carquet_register_codec(CARQUET_COMPRESSION_ZSTD, NULL) != CARQUET_OK)
        TEST_FAIL("register_and_unregister", "unregister failed");
    TEST_PASS("register_and_unregister");
    return 0;
}

/* ---- end-to-end roundtrip through the public writer/reader ------------- */

#define ROWS 5000

static int roundtrip_with_codec(xor_state_t* st) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "custom_codec");
    carquet_error_t err = CARQUET_ERROR_INIT;
    int rc = 1;

    int64_t* in = malloc(sizeof(int64_t) * ROWS);
    int64_t* out = malloc(sizeof(int64_t) * ROWS);
    for (int i = 0; i < ROWS; i++) in[i] = ((int64_t)i * 6364136223846793005ULL) ^ (i << 3);

    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) goto done;
    if (carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT64, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) { carquet_schema_free(s); goto done; }

    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_ZSTD;  /* the slot we overrode */
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); goto done; }
    if (carquet_writer_write_batch(w, 0, in, ROWS, NULL, NULL) != CARQUET_OK) {
        carquet_writer_close(w); carquet_schema_free(s); goto done;
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); goto done; }
    carquet_schema_free(s);

    if (st->compress_calls == 0) { printf("  custom compress never invoked\n"); goto done; }

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) goto done;
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    int64_t n = c ? carquet_column_read_batch(c, out, ROWS, NULL, NULL) : -1;
    int ok = (n == ROWS) && (memcmp(in, out, sizeof(int64_t) * ROWS) == 0);
    carquet_column_reader_free(c);
    carquet_reader_close(r);
    if (!ok) { printf("  value mismatch (n=%lld)\n", (long long)n); goto done; }
    if (st->decompress_calls == 0) { printf("  custom decompress never invoked\n"); goto done; }
    rc = 0;

done:
    carquet_test_cleanup(path);
    free(in);
    free(out);
    return rc;
}

static int test_end_to_end_roundtrip(void) {
    xor_state_t st = {0, 0, 0, 0xA5};
    carquet_custom_codec_t codec = {xor_compress, xor_decompress, xor_bound, &st};
    if (carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &codec) != CARQUET_OK)
        TEST_FAIL("end_to_end_roundtrip", "registration failed");
    int failed = roundtrip_with_codec(&st);
    /* Always restore the builtin so a later test isn't affected. */
    carquet_register_codec(CARQUET_COMPRESSION_ZSTD, NULL);
    if (failed) TEST_FAIL("end_to_end_roundtrip", "custom codec roundtrip failed");
    printf("  (custom codec: %d compress, %d decompress calls)\n",
           st.compress_calls, st.decompress_calls);
    TEST_PASS("end_to_end_roundtrip");
    return 0;
}

static int test_builtin_restored_after_unregister(void) {
    /* After unregistering, a ZSTD write/read must succeed via the builtin,
     * with no stray calls into a stale custom codec. */
    xor_state_t st = {0, 0, 0, 0xA5};
    carquet_custom_codec_t codec = {xor_compress, xor_decompress, xor_bound, &st};
    carquet_register_codec(CARQUET_COMPRESSION_ZSTD, &codec);
    carquet_register_codec(CARQUET_COMPRESSION_ZSTD, NULL);  /* restore builtin */

    st.compress_calls = st.decompress_calls = 0;

    char path[512];
    carquet_test_temp_path(path, sizeof(path), "builtin_restored");
    carquet_error_t err = CARQUET_ERROR_INIT;
    int32_t in[256], out[256];
    for (int i = 0; i < 256; i++) in[i] = i * 3 - 100;

    carquet_schema_t* s = carquet_schema_create(&err);
    (void)carquet_schema_add_column(s, "v", CARQUET_PHYSICAL_INT32, NULL,
                                    CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.compression = CARQUET_COMPRESSION_ZSTD;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    int ok = w != NULL &&
             carquet_writer_write_batch(w, 0, in, 256, NULL, NULL) == CARQUET_OK &&
             carquet_writer_close(w) == CARQUET_OK;
    carquet_schema_free(s);
    if (ok) {
        carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
        carquet_column_reader_t* c = r ? carquet_reader_get_column(r, 0, 0, &err) : NULL;
        int64_t n = c ? carquet_column_read_batch(c, out, 256, NULL, NULL) : -1;
        ok = (n == 256) && (memcmp(in, out, sizeof(in)) == 0);
        carquet_column_reader_free(c);
        if (r) carquet_reader_close(r);
    }
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("builtin_restored_after_unregister", "builtin ZSTD roundtrip failed");
    if (st.compress_calls != 0 || st.decompress_calls != 0)
        TEST_FAIL("builtin_restored_after_unregister", "stale custom codec still invoked");
    TEST_PASS("builtin_restored_after_unregister");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_reject_uncompressed_slot();
    failures += test_reject_out_of_range();
    failures += test_reject_null_function();
    failures += test_register_and_unregister();
    failures += test_end_to_end_roundtrip();
    failures += test_builtin_restored_after_unregister();
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll custom codec tests passed\n");
    return 0;
}
