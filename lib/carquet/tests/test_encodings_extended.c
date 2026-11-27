/**
 * @file test_encodings_extended.c
 * @brief Extended tests for Parquet encodings
 *
 * Tests for:
 * - DELTA_BINARY_PACKED encoding
 * - Dictionary encoding
 * - BYTE_STREAM_SPLIT encoding
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <carquet/error.h>
#include <carquet/types.h>
#include "core/buffer.h"
#include <stdbool.h>
#include <stdint.h>

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * Delta Encoding Function Declarations
 * ============================================================================
 */

carquet_status_t carquet_delta_encode_int32(
    const int32_t* values,
    int32_t num_values,
    uint8_t* data,
    size_t data_capacity,
    size_t* bytes_written);

carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data,
    size_t data_size,
    int32_t* values,
    int32_t num_values,
    size_t* bytes_consumed);

carquet_status_t carquet_delta_encode_int64(
    const int64_t* values,
    int32_t num_values,
    uint8_t* data,
    size_t data_capacity,
    size_t* bytes_written);

carquet_status_t carquet_delta_decode_int64(
    const uint8_t* data,
    size_t data_size,
    int64_t* values,
    int32_t num_values,
    size_t* bytes_consumed);

/* ============================================================================
 * Dictionary Encoding Function Declarations
 * ============================================================================
 */

carquet_status_t carquet_dictionary_encode_int32(
    const int32_t* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    int32_t* output,
    int64_t output_count);

carquet_status_t carquet_dictionary_encode_int64(
    const int64_t* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_decode_int64(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    int64_t* output,
    int64_t output_count);

carquet_status_t carquet_dictionary_encode_float(
    const float* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_decode_float(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    float* output,
    int64_t output_count);

carquet_status_t carquet_dictionary_encode_double(
    const double* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_decode_double(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    double* output,
    int64_t output_count);

/* ============================================================================
 * Byte Stream Split Function Declarations
 * ============================================================================
 */

carquet_status_t carquet_byte_stream_split_encode_float(
    const float* values,
    int64_t count,
    uint8_t* output,
    size_t output_capacity,
    size_t* bytes_written);

carquet_status_t carquet_byte_stream_split_decode_float(
    const uint8_t* data,
    size_t data_size,
    float* values,
    int64_t count);

carquet_status_t carquet_byte_stream_split_encode_double(
    const double* values,
    int64_t count,
    uint8_t* output,
    size_t output_capacity,
    size_t* bytes_written);

carquet_status_t carquet_byte_stream_split_decode_double(
    const uint8_t* data,
    size_t data_size,
    double* values,
    int64_t count);

carquet_status_t carquet_byte_stream_split_encode(
    const uint8_t* values,
    int64_t count,
    int32_t type_length,
    uint8_t* output,
    size_t output_capacity,
    size_t* bytes_written);

carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data,
    size_t data_size,
    int32_t type_length,
    uint8_t* values,
    int64_t count);

extern void carquet_dispatch_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_dispatch_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);

/* ============================================================================
 * Delta Encoding Tests
 * ============================================================================
 */

static int test_delta_int32_sequential(void) {
    int32_t input[100];
    for (int i = 0; i < 100; i++) {
        input[i] = i;
    }

    uint8_t encoded[4096];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int32(
        input, 100, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_sequential", "encode failed");
    }

    /* Sequential data should compress well */
    printf("  [DEBUG] Sequential: 100 int32s -> %zu bytes\n", bytes_written);

    int32_t output[100];
    size_t bytes_consumed;
    status = carquet_delta_decode_int32(
        encoded, bytes_written, output, 100, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_sequential", "decode failed");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            TEST_FAIL("delta_int32_sequential", "value mismatch");
        }
    }

    TEST_PASS("delta_int32_sequential");
    return 0;
}

static int test_delta_int32_negative(void) {
    int32_t input[] = {-100, -50, 0, 50, 100, -200, -150, -100};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t encoded[4096];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int32(
        input, count, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_negative", "encode failed");
    }

    int32_t output[8];
    size_t bytes_consumed;
    status = carquet_delta_decode_int32(
        encoded, bytes_written, output, count, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_negative", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            TEST_FAIL("delta_int32_negative", "value mismatch");
        }
    }

    TEST_PASS("delta_int32_negative");
    return 0;
}

static int test_delta_int32_large_jumps(void) {
    int32_t input[] = {0, 1000000, 2000000, 1500000, 0, -1000000};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t encoded[4096];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int32(
        input, count, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_large_jumps", "encode failed");
    }

    int32_t output[6];
    size_t bytes_consumed;
    status = carquet_delta_decode_int32(
        encoded, bytes_written, output, count, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int32_large_jumps", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("delta_int32_large_jumps", "value mismatch");
        }
    }

    TEST_PASS("delta_int32_large_jumps");
    return 0;
}

static int test_delta_int64_sequential(void) {
    int64_t input[100];
    for (int i = 0; i < 100; i++) {
        input[i] = (int64_t)i * 1000000000LL;
    }

    uint8_t encoded[8192];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int64(
        input, 100, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_sequential", "encode failed");
    }

    printf("  [DEBUG] Sequential int64: 100 values -> %zu bytes\n", bytes_written);

    int64_t output[100];
    size_t bytes_consumed;
    status = carquet_delta_decode_int64(
        encoded, bytes_written, output, 100, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_sequential", "decode failed");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %lld, got %lld\n",
                   i, (long long)input[i], (long long)output[i]);
            TEST_FAIL("delta_int64_sequential", "value mismatch");
        }
    }

    TEST_PASS("delta_int64_sequential");
    return 0;
}

static int test_delta_int64_timestamps(void) {
    /* Simulate millisecond timestamps with small variations */
    int64_t base = 1704067200000LL; /* 2024-01-01 00:00:00 UTC */
    int64_t input[50];
    for (int i = 0; i < 50; i++) {
        input[i] = base + i * 1000 + (i % 3) * 10;
    }

    uint8_t encoded[8192];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int64(
        input, 50, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_timestamps", "encode failed");
    }

    printf("  [DEBUG] Timestamps: 50 values -> %zu bytes\n", bytes_written);

    int64_t output[50];
    size_t bytes_consumed;
    status = carquet_delta_decode_int64(
        encoded, bytes_written, output, 50, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_timestamps", "decode failed");
    }

    for (int i = 0; i < 50; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("delta_int64_timestamps", "value mismatch");
        }
    }

    TEST_PASS("delta_int64_timestamps");
    return 0;
}

static int test_delta_single_value(void) {
    int32_t input[] = {42};

    uint8_t encoded[256];
    size_t bytes_written;

    carquet_status_t status = carquet_delta_encode_int32(
        input, 1, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_single_value", "encode failed");
    }

    int32_t output[1];
    size_t bytes_consumed;
    status = carquet_delta_decode_int32(
        encoded, bytes_written, output, 1, &bytes_consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_single_value", "decode failed");
    }

    if (output[0] != input[0]) {
        TEST_FAIL("delta_single_value", "value mismatch");
    }

    TEST_PASS("delta_single_value");
    return 0;
}

/* ============================================================================
 * Dictionary Encoding Tests
 * ============================================================================
 */

static int test_dictionary_int32_unique(void) {
    int32_t input[] = {100, 200, 300, 400, 500};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_int32(
        input, count, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int32_unique", "encode failed");
    }

    printf("  [DEBUG] Dict size: %zu, indices size: %zu\n",
           carquet_buffer_size(&dict_buf), carquet_buffer_size(&indices_buf));

    int32_t output[5];
    status = carquet_dictionary_decode_int32(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        count,
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int32_unique", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("dictionary_int32_unique", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_int32_unique");
    return 0;
}

static int test_dictionary_int32_repeated(void) {
    int32_t input[] = {1, 2, 1, 2, 1, 2, 3, 3, 3, 1};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_int32(
        input, count, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int32_repeated", "encode failed");
    }

    /* Should have only 3 unique values in dictionary */
    size_t expected_dict_size = 3 * sizeof(int32_t);
    if (carquet_buffer_size(&dict_buf) != expected_dict_size) {
        printf("  [DEBUG] Expected dict size %zu, got %zu\n",
               expected_dict_size, carquet_buffer_size(&dict_buf));
        TEST_FAIL("dictionary_int32_repeated", "dictionary size unexpected");
    }

    int32_t output[10];
    status = carquet_dictionary_decode_int32(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        3, /* dict_count */
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int32_repeated", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            TEST_FAIL("dictionary_int32_repeated", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_int32_repeated");
    return 0;
}

static int test_dictionary_int64(void) {
    int64_t input[] = {1000000000000LL, 2000000000000LL, 1000000000000LL, 3000000000000LL};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_int64(
        input, count, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int64", "encode failed");
    }

    int64_t output[4];
    status = carquet_dictionary_decode_int64(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        3, /* 3 unique values */
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_int64", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("dictionary_int64", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_int64");
    return 0;
}

static int test_dictionary_float(void) {
    float input[] = {1.0f, 2.0f, 1.0f, 3.0f, 2.0f, 1.0f};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_float(
        input, count, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_float", "encode failed");
    }

    float output[6];
    status = carquet_dictionary_decode_float(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        3, /* 3 unique values */
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_float", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("dictionary_float", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_float");
    return 0;
}

static int test_dictionary_double(void) {
    double input[] = {3.14159, 2.71828, 3.14159, 1.41421, 2.71828};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_double(
        input, count, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_double", "encode failed");
    }

    double output[5];
    status = carquet_dictionary_decode_double(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        3, /* 3 unique values */
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_double", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("dictionary_double", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_double");
    return 0;
}

static int test_dictionary_many_values(void) {
    /* Test with 1000 values but only 10 unique */
    int32_t input[1000];
    for (int i = 0; i < 1000; i++) {
        input[i] = (i * 17) % 10;
    }

    carquet_buffer_t dict_buf, indices_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&indices_buf);

    carquet_status_t status = carquet_dictionary_encode_int32(
        input, 1000, &dict_buf, &indices_buf);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_many_values", "encode failed");
    }

    printf("  [DEBUG] 1000 values, 10 unique: dict=%zu bytes, indices=%zu bytes\n",
           carquet_buffer_size(&dict_buf), carquet_buffer_size(&indices_buf));

    int32_t output[1000];
    status = carquet_dictionary_decode_int32(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        10, /* 10 unique values */
        carquet_buffer_data_const(&indices_buf), carquet_buffer_size(&indices_buf),
        output, 1000);
    if (status != CARQUET_OK) {
        TEST_FAIL("dictionary_many_values", "decode failed");
    }

    for (int i = 0; i < 1000; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("dictionary_many_values", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&indices_buf);
    TEST_PASS("dictionary_many_values");
    return 0;
}

/* ============================================================================
 * Byte Stream Split Tests
 * ============================================================================
 */

static int test_byte_stream_split_float(void) {
    float input[] = {1.5f, 2.5f, 3.5f, 4.5f, 5.5f};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t encoded[256];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode_float(
        input, count, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_float", "encode failed");
    }

    if (bytes_written != count * sizeof(float)) {
        TEST_FAIL("byte_stream_split_float", "unexpected output size");
    }

    float output[5];
    status = carquet_byte_stream_split_decode_float(
        encoded, bytes_written, output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_float", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %f, got %f\n", i, input[i], output[i]);
            TEST_FAIL("byte_stream_split_float", "value mismatch");
        }
    }

    TEST_PASS("byte_stream_split_float");
    return 0;
}

static int test_byte_stream_split_double(void) {
    double input[] = {1.123456789, 2.234567890, 3.345678901, 4.456789012};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t encoded[256];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode_double(
        input, count, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_double", "encode failed");
    }

    if (bytes_written != count * sizeof(double)) {
        TEST_FAIL("byte_stream_split_double", "unexpected output size");
    }

    double output[4];
    status = carquet_byte_stream_split_decode_double(
        encoded, bytes_written, output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_double", "decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("byte_stream_split_double", "value mismatch");
        }
    }

    TEST_PASS("byte_stream_split_double");
    return 0;
}

static int test_byte_stream_split_generic(void) {
    /* Test with 16-bit fixed length values */
    uint8_t input[] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
    int count = 4; /* 4 values of 2 bytes each */

    uint8_t encoded[256];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode(
        input, count, 2, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_generic", "encode failed");
    }

    /* Expected: all first bytes, then all second bytes */
    /* Input: [0x01,0x02], [0x03,0x04], [0x05,0x06], [0x07,0x08] */
    /* Expected: [0x01,0x03,0x05,0x07], [0x02,0x04,0x06,0x08] */
    uint8_t expected[] = {0x01, 0x03, 0x05, 0x07, 0x02, 0x04, 0x06, 0x08};
    if (memcmp(encoded, expected, 8) != 0) {
        printf("  [DEBUG] Encoded data mismatch\n");
        for (int i = 0; i < 8; i++) {
            printf("    [%d] expected 0x%02x, got 0x%02x\n", i, expected[i], encoded[i]);
        }
        TEST_FAIL("byte_stream_split_generic", "encoded data mismatch");
    }

    uint8_t output[8];
    status = carquet_byte_stream_split_decode(
        encoded, bytes_written, 2, output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_generic", "decode failed");
    }

    if (memcmp(output, input, 8) != 0) {
        TEST_FAIL("byte_stream_split_generic", "decoded data mismatch");
    }

    TEST_PASS("byte_stream_split_generic");
    return 0;
}

static int test_byte_stream_split_generic_32bit(void) {
    uint8_t input[] = {
        0x01, 0x02, 0x03, 0x04,
        0x11, 0x12, 0x13, 0x14,
        0x21, 0x22, 0x23, 0x24,
        0x31, 0x32, 0x33, 0x34
    };
    uint8_t expected[] = {
        0x01, 0x11, 0x21, 0x31,
        0x02, 0x12, 0x22, 0x32,
        0x03, 0x13, 0x23, 0x33,
        0x04, 0x14, 0x24, 0x34
    };
    uint8_t encoded[32];
    uint8_t output[sizeof(input)];
    size_t bytes_written = 0;

    carquet_status_t status = carquet_byte_stream_split_encode(
        input, 4, 4, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_generic_32bit", "encode failed");
    }
    if (memcmp(encoded, expected, sizeof(expected)) != 0) {
        TEST_FAIL("byte_stream_split_generic_32bit", "encoded data mismatch");
    }

    status = carquet_byte_stream_split_decode(encoded, bytes_written, 4, output, 4);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_generic_32bit", "decode failed");
    }
    if (memcmp(output, input, sizeof(input)) != 0) {
        TEST_FAIL("byte_stream_split_generic_32bit", "decoded data mismatch");
    }

    TEST_PASS("byte_stream_split_generic_32bit");
    return 0;
}

static int test_byte_stream_split_special_floats(void) {
    float input[] = {0.0f, -0.0f, INFINITY, -INFINITY, NAN};
    int count = 5;

    uint8_t encoded[256];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode_float(
        input, count, encoded, sizeof(encoded), &bytes_written);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_special_floats", "encode failed");
    }

    float output[5];
    status = carquet_byte_stream_split_decode_float(
        encoded, bytes_written, output, count);
    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_special_floats", "decode failed");
    }

    /* Check zero, negative zero, infinity, negative infinity */
    if (output[0] != 0.0f) TEST_FAIL("byte_stream_split_special_floats", "zero mismatch");
    if (output[2] != INFINITY) TEST_FAIL("byte_stream_split_special_floats", "inf mismatch");
    if (output[3] != -INFINITY) TEST_FAIL("byte_stream_split_special_floats", "-inf mismatch");
    if (!isnan(output[4])) TEST_FAIL("byte_stream_split_special_floats", "nan mismatch");

    TEST_PASS("byte_stream_split_special_floats");
    return 0;
}

/* ============================================================================
 * Batch RLE encoders vs the per-value encoder
 * ============================================================================
 *
 * carquet_rle_encode_levels / carquet_rle_encode_all locate long runs and
 * emit whole stretches at once (1-bit levels through a bit array); their
 * bytes must equal one carquet_rle_encoder_put() per value + flush(), which
 * the serial and page-parallel writer paths rely on. */

#include "encoding/rle.h"
#include "core/allocator.h"

static uint64_t rle_rng = 0x9E3779B97F4A7C15ULL;
static uint32_t rle_rand(void) {
    rle_rng ^= rle_rng << 13; rle_rng ^= rle_rng >> 7; rle_rng ^= rle_rng << 17;
    return (uint32_t)(rle_rng >> 11);
}

static void rle_put_all(const uint32_t* v, int64_t n, int bw, carquet_buffer_t* out) {
    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, out, bw);
    for (int64_t i = 0; i < n; i++) carquet_rle_encoder_put(&enc, v[i]);
    carquet_rle_encoder_flush(&enc);
}

/* Value stream shapes: random, runs of mean length ~1/p, 0/1 "nulls" with
 * density p, alternating, fixed-length runs, very long runs, values wider
 * than bw. */
static void rle_gen(uint32_t* v, int64_t n, int bw, int shape, uint32_t p1000) {
    uint32_t maxv = bw >= 32 ? 0xFFFFFFFFu : ((1u << bw) - 1u);
    uint32_t cur = 0;
    for (int64_t i = 0; i < n; i++) {
        switch (shape) {
        case 0: cur = rle_rand() & maxv; break;
        case 1: if (rle_rand() % 1000 < p1000) cur = rle_rand() & maxv; break;
        case 2: cur = (rle_rand() % 1000 < p1000) ? 0 : (1u & maxv); break;
        case 3: cur = (uint32_t)(i & 1) & maxv; break;
        case 4: cur = (uint32_t)((i / (1 + p1000 % 23)) & 1) & maxv; break;
        case 5: cur = (i % 4099 < 2000) ? (1u & maxv) : 0; break;
        default: cur = rle_rand(); break;
        }
        v[i] = cur;
    }
}

static int test_rle_batch_vs_put(void) {
    enum { MAXN = 70000 };
    uint32_t* v = malloc(MAXN * sizeof(uint32_t));
    int16_t* lv = malloc(MAXN * sizeof(int16_t));
    uint32_t* lu = malloc(MAXN * sizeof(uint32_t));
    if (!v || !lv || !lu) { free(v); free(lv); free(lu); TEST_FAIL("rle_batch_vs_put", "alloc"); }
    int bad = 0;
    for (int it = 0; it < 3000 && !bad; it++) {
        int64_t n;
        switch (rle_rand() % 4) {
        case 0: n = rle_rand() % 64; break;
        case 1: n = rle_rand() % 2000; break;
        case 2: n = 64 * (int64_t)(rle_rand() % 64) + (int64_t)(rle_rand() % 3) - 1; break;
        default: n = rle_rand() % MAXN; break;
        }
        if (n < 0) n = 0;
        int bw = (rle_rand() % 3) ? 1 : (int)(rle_rand() % 33);
        int shape = (int)(rle_rand() % 7);
        uint32_t p1000 = rle_rand() % 1000;
        p1000 = p1000 * p1000 / 1000;   /* bias towards sparse changes / nulls */
        rle_gen(v, n, bw, shape, p1000);
        for (int64_t i = 0; i < n; i++) lv[i] = (int16_t)v[i];
        if (bw == 1 && n > 0 && rle_rand() % 16 == 0) {
            lv[rle_rand() % n] = (rle_rand() & 1) ? 2 : -1;   /* not a 0/1 level */
        }
        for (int64_t i = 0; i < n; i++) lu[i] = (uint32_t)lv[i];

        carquet_buffer_t a, b;
        carquet_buffer_init(&a);
        carquet_buffer_init(&b);
        /* Levels: several streams back to back in one buffer, as the page
         * writer stores one stream per add_values call. */
        if (bw <= 16) {
            int64_t off = 0;
            do {
                int64_t chunk = (rle_rand() % 3 == 0) ? n - off
                                                      : (int64_t)(rle_rand() % (uint64_t)(n - off + 1));
                rle_put_all(lu + off, chunk, bw, &a);
                int64_t ones = -2;
                if (carquet_rle_encode_levels_counted(lv + off, chunk, bw, &b, &ones) != CARQUET_OK) bad = 1;
                if (ones >= 0) {
                    int64_t want = 0;
                    for (int64_t i = 0; i < chunk; i++) want += lv[off + i] == 1;
                    if (ones != want) bad = 2;
                } else if (ones != -1) {
                    bad = 3;
                }
                off += chunk;
            } while (off < n);
            if (a.size != b.size || (a.size && memcmp(a.data, b.data, a.size) != 0)) bad = 4;
        }
        carquet_buffer_clear(&a);
        carquet_buffer_clear(&b);
        rle_put_all(v, n, bw, &a);
        if (carquet_rle_encode_all(v, n, bw, &b) != CARQUET_OK) bad = 5;
        if (a.size != b.size || (a.size && memcmp(a.data, b.data, a.size) != 0)) bad = 6;
        carquet_buffer_destroy(&a);
        carquet_buffer_destroy(&b);
        if (bad) printf("  it=%d n=%lld bw=%d shape=%d p=%u code=%d\n", it, (long long)n, bw, shape, p1000, bad);
    }
    free(v); free(lv); free(lu);
    if (bad) TEST_FAIL("rle_batch_vs_put", "batch encoder bytes differ from per-value encoder");
    TEST_PASS("rle_batch_vs_put");
    return 0;
}

/* Long streams: runs and stretches far beyond one bit-packed run (63
 * groups), and an all-equal stream. */
static int test_rle_batch_long_streams(void) {
    const int64_t n = 1500007;
    uint32_t* v = malloc((size_t)n * sizeof(uint32_t));
    int16_t* lv = malloc((size_t)n * sizeof(int16_t));
    if (!v || !lv) { free(v); free(lv); TEST_FAIL("rle_batch_long_streams", "alloc"); }
    int bad = 0;
    for (int shape = 0; shape < 4 && !bad; shape++) {
        for (int64_t i = 0; i < n; i++) {
            uint32_t x;
            switch (shape) {
            case 0: x = 1; break;
            case 1: x = i < n / 2 ? (rle_rand() & 1) : 1; break;
            case 2: x = (i % 100000) < 99990; break;
            default: x = (rle_rand() % 10) != 0; break;
            }
            v[i] = x;
            lv[i] = (int16_t)x;
        }
        carquet_buffer_t a, b, c;
        carquet_buffer_init(&a); carquet_buffer_init(&b); carquet_buffer_init(&c);
        rle_put_all(v, n, 1, &a);
        if (carquet_rle_encode_levels(lv, n, 1, &b) != CARQUET_OK) bad = 1;
        if (carquet_rle_encode_all(v, n, 1, &c) != CARQUET_OK) bad = 1;
        if (a.size != b.size || memcmp(a.data, b.data, a.size) != 0) bad = 2;
        if (a.size != c.size || memcmp(a.data, c.data, a.size) != 0) bad = 3;
        carquet_buffer_destroy(&a); carquet_buffer_destroy(&b); carquet_buffer_destroy(&c);
    }
    free(v); free(lv);
    if (bad) TEST_FAIL("rle_batch_long_streams", "batch encoder bytes differ from per-value encoder");
    TEST_PASS("rle_batch_long_streams");
    return 0;
}

/* ============================================================================
 * Capped dictionary build vs a first-occurrence reference
 * ============================================================================ */

extern carquet_status_t carquet_dictionary_build_capped(
    carquet_physical_type_t type, int32_t type_length,
    const void* fixed_values, const carquet_byte_array_t* ba_values,
    int64_t count, size_t max_dict_bytes, carquet_buffer_t* dict_output,
    uint32_t** indices_out, int32_t* num_entries, bool* abandoned);

static int test_dictionary_build_capped_reference(void) {
    enum { N = 8000 };
    uint8_t* vals = malloc((size_t)N * 8);
    uint8_t* ref_dict = malloc((size_t)N * 8);
    uint32_t* ref_idx = malloc((size_t)N * sizeof(uint32_t));
    if (!vals || !ref_dict || !ref_idx) { free(vals); free(ref_dict); free(ref_idx); TEST_FAIL("dictionary_build_capped_reference", "alloc"); }
    const carquet_physical_type_t types[] = {
        CARQUET_PHYSICAL_INT32, CARQUET_PHYSICAL_INT64, CARQUET_PHYSICAL_FLOAT,
        CARQUET_PHYSICAL_DOUBLE, CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
        CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY };
    const size_t sizes[] = {4, 8, 4, 8, 4, 8, 3};
    int bad = 0;
    for (int it = 0; it < 140 && !bad; it++) {
        int k = it % 7;
        size_t size = sizes[k];
        int64_t n = 1 + rle_rand() % N;
        uint32_t card = 1 + rle_rand() % (it % 3 == 0 ? 20 : it % 3 == 1 ? 500 : 3000);
        /* Pool of distinct-ish bit patterns: small ints, doubles/floats
         * with integer values (only high bits vary), +-0.0, NaNs. */
        for (int64_t i = 0; i < n; i++) {
            uint32_t c = rle_rand() % card;
            uint8_t* dst = vals + (size_t)i * size;
            if (size == 8) {
                uint64_t x;
                switch (c % 4) {
                case 0: x = c; break;
                case 1: { double d = (double)c; memcpy(&x, &d, 8); break; }
                case 2: { double d = c / 3.0; memcpy(&x, &d, 8); break; }
                default: x = (c == 3) ? 0x8000000000000000ULL : (0x7FF8000000000000ULL | c); break;
                }
                memcpy(dst, &x, 8);
            } else if (size == 4) {
                uint32_t x;
                switch (c % 3) {
                case 0: x = c; break;
                case 1: { float f = (float)c; memcpy(&x, &f, 4); break; }
                default: x = 0x7FC00000u | c; break;
                }
                memcpy(dst, &x, 4);
            } else {
                dst[0] = (uint8_t)c; dst[1] = (uint8_t)(c >> 8); dst[2] = (uint8_t)(c >> 16);
            }
        }
        if (it % 5 == 0 && n > 3) {   /* runs */
            for (int64_t i = 1; i < n; i++) if (rle_rand() % 4) memcpy(vals + i * size, vals + (i - 1) * size, size);
        }
        /* Reference: first-occurrence numbering, PLAIN bytes in order. */
        int64_t nref = 0;
        for (int64_t i = 0; i < n; i++) {
            int64_t e = 0;
            while (e < nref && memcmp(ref_dict + (size_t)e * size, vals + (size_t)i * size, size) != 0) e++;
            if (e == nref) { memcpy(ref_dict + (size_t)nref * size, vals + (size_t)i * size, size); nref++; }
            ref_idx[i] = (uint32_t)e;
        }
        size_t cap = (it % 4 == 3) ? (size_t)(rle_rand() % (nref * size + 1)) : 0;
        carquet_buffer_t dict;
        carquet_buffer_init(&dict);
        carquet_buffer_append_byte(&dict, 0xEE);   /* the payload is appended */
        uint32_t* idx = NULL;
        int32_t ne = -1;
        bool ab = false;
        carquet_status_t st = carquet_dictionary_build_capped(
            types[k], (int32_t)size, vals, NULL, n, cap, &dict, &idx, &ne, &ab);
        bool want_abandon = cap && (size_t)nref * size > cap;
        if (st != CARQUET_OK) bad = 1;
        else if (ab != want_abandon) bad = 2;
        else if (ab) {
            if (idx || dict.size != 1) bad = 3;
        } else if (ne != nref || dict.size != 1 + (size_t)nref * size ||
                   memcmp(dict.data + 1, ref_dict, (size_t)nref * size) != 0 ||
                   memcmp(idx, ref_idx, (size_t)n * sizeof(uint32_t)) != 0) {
            bad = 4;
        }
        if (bad) printf("  it=%d type=%d size=%zu n=%lld card=%u cap=%zu code=%d\n", it, (int)types[k], size, (long long)n, card, cap, bad);
        carquet_mem_free(idx);
        carquet_buffer_destroy(&dict);
    }
    free(vals); free(ref_dict); free(ref_idx);
    if (bad) TEST_FAIL("dictionary_build_capped_reference", "dictionary or indices differ from reference");
    TEST_PASS("dictionary_build_capped_reference");
    return 0;
}

static int test_bool_pack_unpack_dispatch(void) {
    uint8_t packed[] = {0xA5, 0x3C, 0x81};
    uint8_t unpacked[24];
    uint8_t repacked[3];
    uint8_t expected[24];

    for (int i = 0; i < 24; i++) {
        expected[i] = (packed[i / 8] >> (i % 8)) & 1;
    }

    memset(unpacked, 0xFF, sizeof(unpacked));
    carquet_dispatch_unpack_bools(packed, unpacked, 24);
    if (memcmp(unpacked, expected, sizeof(expected)) != 0) {
        TEST_FAIL("bool_pack_unpack_dispatch", "unpack mismatch");
    }

    memset(repacked, 0, sizeof(repacked));
    carquet_dispatch_pack_bools(unpacked, repacked, 24);
    if (memcmp(repacked, packed, sizeof(packed)) != 0) {
        TEST_FAIL("bool_pack_unpack_dispatch", "pack mismatch");
    }

    TEST_PASS("bool_pack_unpack_dispatch");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    int failures = 0;

    printf("=== Extended Encoding Tests ===\n\n");

    printf("--- Delta Encoding Tests ---\n");
    failures += test_delta_int32_sequential();
    failures += test_delta_int32_negative();
    failures += test_delta_int32_large_jumps();
    failures += test_delta_int64_sequential();
    failures += test_delta_int64_timestamps();
    failures += test_delta_single_value();

    printf("\n--- Dictionary Encoding Tests ---\n");
    failures += test_dictionary_int32_unique();
    failures += test_dictionary_int32_repeated();
    failures += test_dictionary_int64();
    failures += test_dictionary_float();
    failures += test_dictionary_double();
    failures += test_dictionary_many_values();

    printf("\n--- Byte Stream Split Tests ---\n");
    failures += test_byte_stream_split_float();
    failures += test_byte_stream_split_double();
    failures += test_byte_stream_split_generic();
    failures += test_byte_stream_split_generic_32bit();
    failures += test_byte_stream_split_special_floats();
    failures += test_bool_pack_unpack_dispatch();

    printf("\n--- Batch RLE / Dictionary Build Tests ---\n");
    failures += test_rle_batch_vs_put();
    failures += test_rle_batch_long_streams();
    failures += test_dictionary_build_capped_reference();

    printf("\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
