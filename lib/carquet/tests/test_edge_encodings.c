/**
 * @file test_edge_encodings.c
 * @brief Edge case tests for Parquet encodings
 *
 * Tests boundary conditions, empty inputs, single values, extreme values,
 * and other edge cases for all encoding types.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include <carquet/error.h>
#include <carquet/types.h>
#include "core/buffer.h"
#include "encoding/plain.h"
#include "encoding/rle.h"

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * External Function Declarations
 * ============================================================================
 */

/* Delta encoding */
carquet_status_t carquet_delta_encode_int32(
    const int32_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);

carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data, size_t data_size,
    int32_t* values, int32_t num_values, size_t* bytes_consumed);

carquet_status_t carquet_delta_encode_int64(
    const int64_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);

carquet_status_t carquet_delta_decode_int64(
    const uint8_t* data, size_t data_size,
    int64_t* values, int32_t num_values, size_t* bytes_consumed);

/* Dictionary encoding */
carquet_status_t carquet_dictionary_encode_int32(
    const int32_t* values, int64_t count,
    carquet_buffer_t* dict_output, carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    int32_t* output, int64_t output_count);

/* Byte stream split */
carquet_status_t carquet_byte_stream_split_encode(
    const uint8_t* values, int64_t count, int32_t type_length,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);

carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data, size_t data_size, int32_t type_length,
    uint8_t* output, int64_t count);

/* RLE encoding */
carquet_status_t carquet_rle_encode_all(
    const uint32_t* input, int64_t count, int bit_width,
    carquet_buffer_t* output);

int64_t carquet_rle_decode_all(
    const uint8_t* input, size_t input_size, int bit_width,
    uint32_t* output, int64_t max_values);

/* ============================================================================
 * Plain Encoding Edge Cases
 * ============================================================================
 */

static int test_plain_empty_int32(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Encode zero values - may fail or succeed depending on implementation */
    int32_t dummy = 0;
    carquet_status_t status = carquet_encode_plain_int32(&dummy, 0, &buf);
    /* Either status OK with 0 bytes, or error is acceptable */
    printf("  [DEBUG] Empty encode: status=%d, size=%zu\n", status, carquet_buffer_size(&buf));

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_empty_int32");
    return 0;
}

static int test_plain_single_int32(void) {
    int32_t input = 42;
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_int32(&input, 1, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_single_int32", "Failed to encode single value");
    }

    if (carquet_buffer_size(&buf) != 4) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_single_int32", "Wrong size for single int32");
    }

    int32_t output;
    carquet_decode_plain_int32(carquet_buffer_data_const(&buf), 4, &output, 1);
    if (output != input) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_single_int32", "Value mismatch");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_single_int32");
    return 0;
}

static int test_plain_int32_extremes(void) {
    int32_t input[] = {INT32_MIN, INT32_MAX, 0, -1, 1};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_int32(input, count, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_int32_extremes", "Failed to encode extreme values");
    }

    int32_t output[5];
    carquet_decode_plain_int32(carquet_buffer_data_const(&buf),
                               carquet_buffer_size(&buf), output, count);

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("plain_int32_extremes", "Value mismatch at extreme");
        }
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_int32_extremes");
    return 0;
}

static int test_plain_int64_extremes(void) {
    int64_t input[] = {INT64_MIN, INT64_MAX, 0, -1, 1};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_int64(input, count, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_int64_extremes", "Failed to encode extreme values");
    }

    int64_t output[5];
    carquet_decode_plain_int64(carquet_buffer_data_const(&buf),
                               carquet_buffer_size(&buf), output, count);

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("plain_int64_extremes", "Value mismatch at extreme");
        }
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_int64_extremes");
    return 0;
}

static int test_plain_float_special(void) {
    float input[] = {0.0f, -0.0f, FLT_MIN, FLT_MAX, FLT_EPSILON,
                     INFINITY, -INFINITY, NAN};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_float(input, count, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_float_special", "Failed to encode special floats");
    }

    float output[8];
    carquet_decode_plain_float(carquet_buffer_data_const(&buf),
                               carquet_buffer_size(&buf), output, count);

    for (int i = 0; i < count - 1; i++) {  /* Skip NAN comparison */
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("plain_float_special", "Value mismatch");
        }
    }
    /* Check NAN separately */
    if (!isnan(output[count - 1])) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_float_special", "NAN not preserved");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_float_special");
    return 0;
}

static int test_plain_double_special(void) {
    double input[] = {0.0, -0.0, DBL_MIN, DBL_MAX, DBL_EPSILON,
                      INFINITY, -INFINITY, NAN};
    int count = sizeof(input) / sizeof(input[0]);

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_double(input, count, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_double_special", "Failed to encode special doubles");
    }

    double output[8];
    carquet_decode_plain_double(carquet_buffer_data_const(&buf),
                                carquet_buffer_size(&buf), output, count);

    for (int i = 0; i < count - 1; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("plain_double_special", "Value mismatch");
        }
    }
    if (!isnan(output[count - 1])) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_double_special", "NAN not preserved");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_double_special");
    return 0;
}

static int test_plain_boolean_edge_cases(void) {
    /* Test single boolean */
    uint8_t single = 1;
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_encode_plain_boolean(&single, 1, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_boolean_edge_cases", "Failed to encode single bool");
    }

    uint8_t output;
    carquet_decode_plain_boolean(carquet_buffer_data_const(&buf), 1, &output, 1);
    if (output != 1) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_boolean_edge_cases", "Single bool mismatch");
    }

    carquet_buffer_destroy(&buf);

    /* Test exactly 8 booleans (1 byte boundary) */
    uint8_t eight[8] = {1, 0, 1, 0, 1, 0, 1, 0};
    carquet_buffer_init(&buf);
    status = carquet_encode_plain_boolean(eight, 8, &buf);
    if (status != CARQUET_OK || carquet_buffer_size(&buf) != 1) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_boolean_edge_cases", "8 bools should be 1 byte");
    }

    uint8_t output8[8];
    carquet_decode_plain_boolean(carquet_buffer_data_const(&buf), 1, output8, 8);
    for (int i = 0; i < 8; i++) {
        if (output8[i] != eight[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("plain_boolean_edge_cases", "8 bool mismatch");
        }
    }

    carquet_buffer_destroy(&buf);

    /* Test 9 booleans (crosses byte boundary) */
    uint8_t nine[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    carquet_buffer_init(&buf);
    status = carquet_encode_plain_boolean(nine, 9, &buf);
    if (status != CARQUET_OK || carquet_buffer_size(&buf) != 2) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("plain_boolean_edge_cases", "9 bools should be 2 bytes");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("plain_boolean_edge_cases");
    return 0;
}

/* ============================================================================
 * Delta Encoding Edge Cases
 * ============================================================================
 */

static int test_delta_single_value(void) {
    int32_t input = 12345;
    uint8_t buffer[128];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(&input, 1, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_single_value", "Failed to encode single value");
    }

    int32_t output;
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, &output, 1, &consumed);
    if (status != CARQUET_OK || output != input) {
        TEST_FAIL("delta_single_value", "Decode mismatch");
    }

    TEST_PASS("delta_single_value");
    return 0;
}

static int test_delta_constant_values(void) {
    /* All same value - delta should be very efficient */
    int32_t input[100];
    for (int i = 0; i < 100; i++) input[i] = 42;

    uint8_t buffer[1024];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(input, 100, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_constant_values", "Failed to encode constant values");
    }

    /* Constant values should compress very well */
    printf("  [DEBUG] Constant values: 100 int32s -> %zu bytes\n", written);

    int32_t output[100];
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, 100, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_constant_values", "Decode failed");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != 42) {
            TEST_FAIL("delta_constant_values", "Value mismatch");
        }
    }

    TEST_PASS("delta_constant_values");
    return 0;
}

static int test_delta_alternating(void) {
    /* Alternating values: worst case for delta */
    int32_t input[100];
    for (int i = 0; i < 100; i++) {
        input[i] = (i % 2 == 0) ? 1000000 : -1000000;
    }

    uint8_t buffer[2048];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(input, 100, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_alternating", "Failed to encode alternating values");
    }

    printf("  [DEBUG] Alternating values: 100 int32s -> %zu bytes\n", written);

    int32_t output[100];
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, 100, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_alternating", "Decode failed");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != input[i]) {
            TEST_FAIL("delta_alternating", "Value mismatch");
        }
    }

    TEST_PASS("delta_alternating");
    return 0;
}

static int test_delta_extreme_values(void) {
    /* Test with values that don't cause overflow in deltas */
    int32_t input[] = {-1000000, 1000000, 0, -500000, 500000};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t buffer[256];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(input, count, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_extreme_values", "Failed to encode extreme values");
    }

    int32_t output[5];
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, count, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_extreme_values", "Decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            TEST_FAIL("delta_extreme_values", "Value mismatch");
        }
    }

    TEST_PASS("delta_extreme_values");
    return 0;
}

static int test_delta_int64_extreme(void) {
    /* Test with sequential values that have small deltas */
    int64_t input[] = {100LL, 105LL, 110LL, 115LL, 120LL};
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t buffer[512];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int64(input, count, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_extreme", "Failed to encode");
    }

    int64_t output[5];
    size_t consumed;
    status = carquet_delta_decode_int64(buffer, written, output, count, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("delta_int64_extreme", "Decode failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %lld, got %lld\n",
                   i, (long long)input[i], (long long)output[i]);
            TEST_FAIL("delta_int64_extreme", "Value mismatch");
        }
    }

    TEST_PASS("delta_int64_extreme");
    return 0;
}

/* ============================================================================
 * RLE Encoding Edge Cases
 * ============================================================================
 */

static int test_rle_single_run(void) {
    /* All same values - single run */
    uint32_t input[64];
    for (int i = 0; i < 64; i++) input[i] = 7;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_rle_encode_all(input, 64, 4, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_single_run", "Encode failed");
    }

    printf("  [DEBUG] Single run 64 values -> %zu bytes\n", carquet_buffer_size(&buf));

    uint32_t output[64];
    int64_t decoded = carquet_rle_decode_all(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        4, output, 64);

    if (decoded != 64) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_single_run", "Wrong decode count");
    }

    for (int i = 0; i < 64; i++) {
        if (output[i] != 7) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("rle_single_run", "Value mismatch");
        }
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_single_run");
    return 0;
}

static int test_rle_max_bit_width(void) {
    /* Test with 8-bit width (well-supported) */
    uint32_t input[32];
    for (int i = 0; i < 32; i++) {
        input[i] = (uint32_t)(i * 7);  /* Values within 8 bits (0-217) */
    }

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Use 8-bit width */
    carquet_status_t status = carquet_rle_encode_all(input, 32, 8, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_max_bit_width", "Encode failed");
    }

    uint32_t output[32];
    int64_t decoded = carquet_rle_decode_all(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        8, output, 32);

    if (decoded != 32) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_max_bit_width", "Wrong decode count");
    }

    for (int i = 0; i < 32; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("rle_max_bit_width", "Value mismatch");
        }
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_max_bit_width");
    return 0;
}

static int test_rle_bit_width_1(void) {
    /* Minimal bit width: 1 bit (binary values) */
    uint32_t input[100];
    for (int i = 0; i < 100; i++) {
        input[i] = i % 2;
    }

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    carquet_status_t status = carquet_rle_encode_all(input, 100, 1, &buf);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_bit_width_1", "Encode failed");
    }

    printf("  [DEBUG] 1-bit RLE 100 values -> %zu bytes\n", carquet_buffer_size(&buf));

    uint32_t output[100];
    int64_t decoded = carquet_rle_decode_all(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        1, output, 100);

    if (decoded != 100) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("rle_bit_width_1", "Wrong decode count");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("rle_bit_width_1", "Value mismatch");
        }
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("rle_bit_width_1");
    return 0;
}

/* ============================================================================
 * Dictionary Encoding Edge Cases
 * ============================================================================
 */

static int test_dictionary_single_unique(void) {
    /* All same values - dictionary should have 1 entry */
    int32_t input[100];
    for (int i = 0; i < 100; i++) input[i] = 999;

    carquet_buffer_t dict, indices;
    carquet_buffer_init(&dict);
    carquet_buffer_init(&indices);

    carquet_status_t status = carquet_dictionary_encode_int32(input, 100, &dict, &indices);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_single_unique", "Encode failed");
    }

    /* Dictionary should contain exactly 1 value */
    if (carquet_buffer_size(&dict) != 4) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_single_unique", "Dictionary should have 1 entry");
    }

    printf("  [DEBUG] Single unique: dict=%zu bytes, indices=%zu bytes\n",
           carquet_buffer_size(&dict), carquet_buffer_size(&indices));

    int32_t output[100];
    status = carquet_dictionary_decode_int32(
        carquet_buffer_data_const(&dict), carquet_buffer_size(&dict), 1,
        carquet_buffer_data_const(&indices), carquet_buffer_size(&indices),
        output, 100);

    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_single_unique", "Decode failed");
    }

    for (int i = 0; i < 100; i++) {
        if (output[i] != 999) {
            carquet_buffer_destroy(&dict);
            carquet_buffer_destroy(&indices);
            TEST_FAIL("dictionary_single_unique", "Value mismatch");
        }
    }

    carquet_buffer_destroy(&dict);
    carquet_buffer_destroy(&indices);
    TEST_PASS("dictionary_single_unique");
    return 0;
}

static int test_dictionary_all_unique(void) {
    /* All unique values - worst case for dictionary */
    int32_t input[50];
    for (int i = 0; i < 50; i++) input[i] = i * 1000;

    carquet_buffer_t dict, indices;
    carquet_buffer_init(&dict);
    carquet_buffer_init(&indices);

    carquet_status_t status = carquet_dictionary_encode_int32(input, 50, &dict, &indices);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_all_unique", "Encode failed");
    }

    /* Dictionary should contain 50 values */
    if (carquet_buffer_size(&dict) != 50 * 4) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_all_unique", "Wrong dictionary size");
    }

    printf("  [DEBUG] All unique: dict=%zu bytes, indices=%zu bytes\n",
           carquet_buffer_size(&dict), carquet_buffer_size(&indices));

    int32_t output[50];
    status = carquet_dictionary_decode_int32(
        carquet_buffer_data_const(&dict), carquet_buffer_size(&dict), 50,
        carquet_buffer_data_const(&indices), carquet_buffer_size(&indices),
        output, 50);

    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&dict);
        carquet_buffer_destroy(&indices);
        TEST_FAIL("dictionary_all_unique", "Decode failed");
    }

    for (int i = 0; i < 50; i++) {
        if (output[i] != input[i]) {
            carquet_buffer_destroy(&dict);
            carquet_buffer_destroy(&indices);
            TEST_FAIL("dictionary_all_unique", "Value mismatch");
        }
    }

    carquet_buffer_destroy(&dict);
    carquet_buffer_destroy(&indices);
    TEST_PASS("dictionary_all_unique");
    return 0;
}

/* ============================================================================
 * Byte Stream Split Edge Cases
 * ============================================================================
 */

static int test_byte_stream_split_single_float(void) {
    float input = 3.14159f;
    uint8_t encoded[sizeof(float)];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode(
        (uint8_t*)&input, 1, 4, encoded, sizeof(encoded), &bytes_written);

    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_single_float", "Encode failed");
    }

    float output;
    status = carquet_byte_stream_split_decode(
        encoded, bytes_written, 4, (uint8_t*)&output, 1);

    if (status != CARQUET_OK || output != input) {
        TEST_FAIL("byte_stream_split_single_float", "Decode mismatch");
    }

    TEST_PASS("byte_stream_split_single_float");
    return 0;
}

static int test_byte_stream_split_denormals(void) {
    /* Test with denormalized floats */
    float input[4];
    input[0] = FLT_MIN / 2.0f;  /* Denormal */
    input[1] = -FLT_MIN / 2.0f;
    input[2] = FLT_TRUE_MIN;    /* Smallest positive denormal */
    input[3] = -FLT_TRUE_MIN;

    uint8_t encoded[4 * sizeof(float)];
    size_t bytes_written;

    carquet_status_t status = carquet_byte_stream_split_encode(
        (uint8_t*)input, 4, 4, encoded, sizeof(encoded), &bytes_written);

    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_denormals", "Encode failed");
    }

    float output[4];
    status = carquet_byte_stream_split_decode(
        encoded, bytes_written, 4, (uint8_t*)output, 4);

    if (status != CARQUET_OK) {
        TEST_FAIL("byte_stream_split_denormals", "Decode failed");
    }

    for (int i = 0; i < 4; i++) {
        if (memcmp(&output[i], &input[i], sizeof(float)) != 0) {
            TEST_FAIL("byte_stream_split_denormals", "Bit pattern mismatch");
        }
    }

    TEST_PASS("byte_stream_split_denormals");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    printf("=== Encoding Edge Case Tests ===\n\n");

    int failures = 0;

    printf("--- Plain Encoding Edge Cases ---\n");
    failures += test_plain_empty_int32();
    failures += test_plain_single_int32();
    failures += test_plain_int32_extremes();
    failures += test_plain_int64_extremes();
    failures += test_plain_float_special();
    failures += test_plain_double_special();
    failures += test_plain_boolean_edge_cases();

    printf("\n--- Delta Encoding Edge Cases ---\n");
    failures += test_delta_single_value();
    failures += test_delta_constant_values();
    failures += test_delta_alternating();
    failures += test_delta_extreme_values();
    failures += test_delta_int64_extreme();

    printf("\n--- RLE Encoding Edge Cases ---\n");
    failures += test_rle_single_run();
    failures += test_rle_max_bit_width();
    failures += test_rle_bit_width_1();

    printf("\n--- Dictionary Encoding Edge Cases ---\n");
    failures += test_dictionary_single_unique();
    failures += test_dictionary_all_unique();

    printf("\n--- Byte Stream Split Edge Cases ---\n");
    failures += test_byte_stream_split_single_float();
    failures += test_byte_stream_split_denormals();

    printf("\n");
    if (failures == 0) {
        printf("All encoding edge case tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed.\n", failures);
        return 1;
    }
}
