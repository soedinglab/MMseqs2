/**
 * @file test_encodings.c
 * @brief Tests for Parquet encodings
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "encoding/plain.h"
#include "encoding/rle.h"
#include "core/buffer.h"
#include <carquet/error.h>
#include <stdint.h>

/* Generic (FIXED_LEN_BYTE_ARRAY) BYTE_STREAM_SPLIT — no public header. */
extern carquet_status_t carquet_byte_stream_split_encode(
    const uint8_t* values, int64_t count, int32_t type_length,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);
extern carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data, size_t data_size, int32_t type_length,
    uint8_t* values, int64_t count);

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * PLAIN Encoding Tests
 * ============================================================================
 */

static int test_plain_int32(void) {
    int32_t input[] = {1, 2, 3, 4, 5, -100, 0, 2147483647};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Encode */
    assert(carquet_encode_plain_int32(input, count, &buf) == CARQUET_OK);
    assert(carquet_buffer_size(&buf) == (size_t)(count * 4));

    /* Decode */
    int32_t output[8];
    int64_t bytes = carquet_decode_plain_int32(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        output, count);
    (void)bytes;

    assert(bytes == count * 4);
    for (int i = 0; i < count; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_int32");
    return 0;
}

static int test_plain_int64(void) {
    int64_t input[] = {1, -1, 0x7FFFFFFFFFFFFFFFLL, (int64_t)0x8000000000000000LL};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_encode_plain_int64(input, count, &buf) == CARQUET_OK);

    int64_t output[4];
    int64_t bytes = carquet_decode_plain_int64(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        output, count);
    (void)bytes;

    assert(bytes == count * 8);
    for (int i = 0; i < count; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_int64");
    return 0;
}

static int test_plain_boolean(void) {
    uint8_t input[] = {1, 0, 1, 1, 0, 0, 1, 0, 1};  /* 9 booleans */
    int count = sizeof(input);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_encode_plain_boolean(input, count, &buf) == CARQUET_OK);
    assert(carquet_buffer_size(&buf) == 2);  /* 9 bits = 2 bytes */

    uint8_t output[9];
    int64_t bytes = carquet_decode_plain_boolean(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        output, count);
    (void)bytes;

    assert(bytes == 2);
    for (int i = 0; i < count; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_boolean");
    return 0;
}

static int test_plain_boolean_large(void) {
    uint8_t input[73];
    uint8_t output[73];

    for (int i = 0; i < 73; i++) {
        input[i] = (uint8_t)(((i * 5) + 3) & 1);
    }

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_encode_plain_boolean(input, 73, &buf) == CARQUET_OK);
    assert(carquet_buffer_size(&buf) == 10);  /* 73 bits = 10 bytes */

    int64_t bytes = carquet_decode_plain_boolean(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        output, 73);
    (void)bytes;

    assert(bytes == 10);
    for (int i = 0; i < 73; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_boolean_large");
    return 0;
}

static int test_plain_boolean_append(void) {
    /* A bit stream built in several appends must be byte-identical to the same
     * values packed in one go, wherever the cuts fall -- including cuts that
     * leave the stream mid-byte and appends that add nothing. The reference is
     * packed by hand here: value i is bit (i % 8) of byte (i / 8). */
    uint8_t input[73];
    for (int i = 0; i < 73; i++) {
        input[i] = (uint8_t)((((i * 5) + 3) % 3) == 0);
    }

    for (int n = 0; n <= 73; n++) {
        uint8_t expect[10] = {0};
        size_t expect_size = ((size_t)n + 7) / 8;
        for (int i = 0; i < n; i++) {
            if (input[i]) expect[i / 8] |= (uint8_t)(1u << (i % 8));
        }

        /* Two appends, cut at every position. */
        for (int k = 0; k <= n; k++) {
            carquet_buffer_t buf;
            carquet_buffer_init(&buf);
            assert(carquet_encode_plain_boolean_append(input, k, 0, &buf) == CARQUET_OK);
            assert(carquet_encode_plain_boolean_append(input + k, n - k, k, &buf) == CARQUET_OK);
            assert(carquet_buffer_size(&buf) == expect_size);
            assert(n == 0 || memcmp(carquet_buffer_data_const(&buf), expect, expect_size) == 0);
            carquet_buffer_destroy(&buf);
        }

        /* One value per append. */
        carquet_buffer_t buf;
        carquet_buffer_init(&buf);
        for (int i = 0; i < n; i++) {
            assert(carquet_encode_plain_boolean_append(input + i, 1, i, &buf) == CARQUET_OK);
        }
        assert(carquet_buffer_size(&buf) == expect_size);
        assert(n == 0 || memcmp(carquet_buffer_data_const(&buf), expect, expect_size) == 0);
        carquet_buffer_destroy(&buf);
    }

    /* A stream that claims to be mid-byte must have that byte. */
    carquet_buffer_t empty;
    carquet_buffer_init(&empty);
    assert(carquet_encode_plain_boolean_append(input, 4, 3, &empty) == CARQUET_ERROR_INVALID_ARGUMENT);
    assert(carquet_buffer_size(&empty) == 0);

    /* Encoding no values is not an error: an all-null batch does exactly that. */
    assert(carquet_encode_plain_boolean(input, 0, &empty) == CARQUET_OK);
    assert(carquet_buffer_size(&empty) == 0);
    carquet_buffer_destroy(&empty);

    TEST_PASS("plain_boolean_append");
    return 0;
}

static int test_plain_double(void) {
    double input[] = {0.0, 1.0, -1.0, 3.14159265359, 1e100};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_encode_plain_double(input, count, &buf) == CARQUET_OK);

    double output[5];
    int64_t bytes = carquet_decode_plain_double(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        output, count);
    (void)bytes;

    assert(bytes == count * 8);
    for (int i = 0; i < count; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_double");
    return 0;
}

/* ============================================================================
 * RLE Encoding Tests
 * ============================================================================
 */

static int test_rle_repeated_values(void) {
    /* 100 repeated zeros */
    uint32_t input[100];
    memset(input, 0, sizeof(input));

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_rle_encode_all(input, 100, 1, &buf) == CARQUET_OK);

    /* Should be very compact */
    assert(carquet_buffer_size(&buf) < 10);

    uint32_t output[100];
    int64_t count = carquet_rle_decode_all(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        1, output, 100);
    (void)count;

    assert(count == 100);
    for (int i = 0; i < 100; i++) {
        assert(output[i] == 0);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_repeated_values");
    return 0;
}

static int test_rle_alternating(void) {
    /* Alternating 0 and 1 */
    uint32_t input[16];
    for (int i = 0; i < 16; i++) {
        input[i] = i % 2;
    }

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_rle_encode_all(input, 16, 1, &buf);
    assert(status == CARQUET_OK);
    (void)status;
    (void)input[0];  /* Silence false -Wunused-but-set-variable warning */

    uint32_t output[16];
    int64_t count = carquet_rle_decode_all(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        1, output, 16);
    (void)count;

    assert(count == 16);
    for (int i = 0; i < 16; i++) {
        assert(output[i] == (uint32_t)(i % 2));
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_alternating");
    return 0;
}

static int test_rle_decoder_skip(void) {
    /* Create RLE data with repeated values */
    uint32_t input[100];
    for (int i = 0; i < 100; i++) {
        input[i] = i / 10;  /* 0,0,0,...,1,1,1,...,9,9,9,... */
    }

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_rle_encode_all(input, 100, 4, &buf);
    assert(status == CARQUET_OK);
    (void)status;
    (void)input[0];  /* Silence false -Wunused-but-set-variable warning */

    carquet_rle_decoder_t dec;
    carquet_rle_decoder_init(&dec,
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf), 4);

    /* Skip first 25 values */
    int64_t skipped = carquet_rle_decoder_skip(&dec, 25);
    (void)skipped;
    assert(skipped == 25);

    /* Read next 10 values (should be 2s and 3s) */
    uint32_t output[10];
    int64_t read = carquet_rle_decoder_get_batch(&dec, output, 10);
    (void)read;
    assert(read == 10);

    for (int i = 0; i < 5; i++) {
        assert(output[i] == 2);
    }
    for (int i = 5; i < 10; i++) {
        assert(output[i] == 3);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_decoder_skip");
    return 0;
}

static int test_rle_levels(void) {
    int16_t input[] = {0, 0, 1, 0, 1, 1, 0, 0, 1, 0};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    assert(carquet_rle_encode_levels(input, count, 1, &buf) == CARQUET_OK);

    int16_t output[10];
    int64_t decoded = carquet_rle_decode_levels(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        1, output, count);
    (void)decoded;

    assert(decoded == count);
    for (int i = 0; i < count; i++) {
        assert(output[i] == input[i]);
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_levels");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

/* ============================================================================
 * BYTE_STREAM_SPLIT FLBA alignment regression
 * ============================================================================
 * The type_length 4/8 fast paths reinterpret the value buffer as a float or
 * double pointer. For FIXED_LEN_BYTE_ARRAY the buffer carries no such alignment
 * guarantee, so the generic transpose must be used when the pointer is
 * unaligned. This test drives both an aligned and a deliberately misaligned
 * buffer and asserts they produce byte-identical, correct results (and do not
 * fault under UBSan).
 */
static int test_bss_flba_misaligned(void) {
    for (int32_t type_length = 4; type_length <= 8; type_length += 4) {
        const int64_t count = 17;  /* odd count, not a SIMD multiple */
        const size_t nbytes = (size_t)count * (size_t)type_length;

        /* Separate backing arrays so the aligned and misaligned inputs never
         * alias. Each holds the same logical bytes, one at a naturally aligned
         * base and one at a guaranteed-misaligned (+1) base. */
        uint8_t in_a_backing[8 * 17 + 16];
        uint8_t in_m_backing[8 * 17 + 16];
        uint8_t out_a_backing[8 * 17 + 16];
        uint8_t out_m_backing[8 * 17 + 16];

        uint8_t* in_aligned = in_a_backing;
        while (((uintptr_t)in_aligned & 7u) != 0) in_aligned++;
        uint8_t* in_mis = in_m_backing;
        while (((uintptr_t)in_mis & 7u) != 0) in_mis++;
        in_mis += 1;  /* now not 4/8-aligned */

        for (size_t i = 0; i < nbytes; i++) {
            uint8_t v = (uint8_t)(i * 31u + 7u);
            in_aligned[i] = v;
            in_mis[i] = v;
        }

        uint8_t enc_aligned[8 * 17];
        uint8_t enc_mis[8 * 17];
        size_t w1 = 0, w2 = 0;

        assert(carquet_byte_stream_split_encode(
            in_aligned, count, type_length, enc_aligned, sizeof(enc_aligned), &w1) == CARQUET_OK);
        assert(carquet_byte_stream_split_encode(
            in_mis, count, type_length, enc_mis, sizeof(enc_mis), &w2) == CARQUET_OK);
        /* The aligned fast path and the misaligned generic transpose must
         * produce identical, canonical BYTE_STREAM_SPLIT output. */
        assert(w1 == nbytes && w2 == nbytes);
        assert(memcmp(enc_aligned, enc_mis, nbytes) == 0);

        /* Decode back into an aligned and a misaligned buffer; both must
         * reconstruct the original exactly. */
        uint8_t* dec_aligned = out_a_backing;
        while (((uintptr_t)dec_aligned & 7u) != 0) dec_aligned++;
        uint8_t* dec_mis = out_m_backing;
        while (((uintptr_t)dec_mis & 7u) != 0) dec_mis++;
        dec_mis += 1;

        assert(carquet_byte_stream_split_decode(
            enc_aligned, nbytes, type_length, dec_aligned, count) == CARQUET_OK);
        assert(carquet_byte_stream_split_decode(
            enc_mis, nbytes, type_length, dec_mis, count) == CARQUET_OK);
        assert(memcmp(dec_aligned, in_aligned, nbytes) == 0);
        assert(memcmp(dec_mis, in_mis, nbytes) == 0);
    }

    TEST_PASS("bss_flba_misaligned");
    return 0;
}

int main(void) {
    int failures = 0;

    printf("=== Encoding Tests ===\n\n");

    /* PLAIN encoding tests */
    failures += test_plain_int32();
    failures += test_plain_int64();
    failures += test_plain_boolean();
    failures += test_plain_boolean_large();
    failures += test_plain_boolean_append();
    failures += test_plain_double();

    /* RLE encoding tests */
    failures += test_rle_repeated_values();
    failures += test_rle_alternating();
    failures += test_rle_decoder_skip();
    failures += test_rle_levels();

    /* BYTE_STREAM_SPLIT FLBA alignment regression */
    failures += test_bss_flba_misaligned();

    printf("\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
