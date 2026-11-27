/**
 * @file test_dict_decode.c
 * @brief Tests for the module-level dictionary decode API for the variable /
 *        fixed-width byte types: carquet_dictionary_decode_byte_array and
 *        carquet_dictionary_decode_fixed_len_byte_array.
 *
 * Encodes with the matching module encoder, decodes back, and asserts an exact
 * round-trip, plus error and boundary paths.
 */

#include <stdio.h>
#include <string.h>

#include <carquet/carquet.h>
#include "core/buffer.h"
#include "encoding/dictionary.h"
#include "test_helpers.h"

/* Encoders (declared inline at their call sites in the codebase). */
carquet_status_t carquet_dictionary_encode_byte_array(
    const carquet_byte_array_t* values, int64_t count,
    carquet_buffer_t* dict_output, carquet_buffer_t* indices_output);

carquet_status_t carquet_dictionary_encode_capped(
    carquet_physical_type_t type, int32_t type_length,
    const void* fixed_values, const carquet_byte_array_t* ba_values,
    int64_t count, size_t max_dict_bytes,
    carquet_buffer_t* dict_output, carquet_buffer_t* indices_output,
    bool* abandoned);

static carquet_byte_array_t ba(const char* s) {
    carquet_byte_array_t v;
    v.data = (uint8_t*)s;
    v.length = (int32_t)strlen(s);
    return v;
}

/* ---- BYTE_ARRAY round-trip (unique + repeated) -------------------------- */
static int test_byte_array_roundtrip(void) {
    /* Repeated values so the dictionary has fewer entries than rows. */
    carquet_byte_array_t input[] = {
        ba("apple"), ba("banana"), ba("apple"), ba("cherry"),
        ba("banana"), ba("apple"), ba(""), ba("date-fruit"),
    };
    int64_t count = (int64_t)(sizeof(input) / sizeof(input[0]));
    int32_t dict_count = 5;  /* apple, banana, cherry, "", date-fruit */

    carquet_buffer_t dict_buf, idx_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&idx_buf);

    if (carquet_dictionary_encode_byte_array(input, count, &dict_buf, &idx_buf)
        != CARQUET_OK) {
        TEST_FAIL("byte_array_roundtrip", "encode failed");
    }

    carquet_byte_array_t out[8];
    carquet_status_t st = carquet_dictionary_decode_byte_array(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        dict_count,
        carquet_buffer_data_const(&idx_buf), carquet_buffer_size(&idx_buf),
        out, count);
    if (st != CARQUET_OK) TEST_FAIL("byte_array_roundtrip", "decode failed");

    for (int64_t i = 0; i < count; i++) {
        if (out[i].length != input[i].length ||
            memcmp(out[i].data, input[i].data, (size_t)out[i].length) != 0) {
            TEST_FAIL("byte_array_roundtrip", "value mismatch");
        }
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&idx_buf);
    TEST_PASS("byte_array_roundtrip");
    return 0;
}

/* ---- FIXED_LEN_BYTE_ARRAY round-trip ------------------------------------ */
static int test_flba_roundtrip(void) {
    const int32_t tl = 4;
    /* 6 values of width 4, with repeats (3 unique). */
    const uint8_t vals[6 * 4] = {
        'a','a','a','a',  'b','b','b','b',  'a','a','a','a',
        'c','c','c','c',  'b','b','b','b',  'a','a','a','a',
    };
    int64_t count = 6;
    int32_t dict_count = 3;

    carquet_buffer_t dict_buf, idx_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&idx_buf);

    bool abandoned = false;
    if (carquet_dictionary_encode_capped(
            CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, tl, vals, NULL, count, 0,
            &dict_buf, &idx_buf, &abandoned) != CARQUET_OK || abandoned) {
        TEST_FAIL("flba_roundtrip", "encode failed");
    }

    uint8_t out[6 * 4];
    carquet_status_t st = carquet_dictionary_decode_fixed_len_byte_array(
        carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
        dict_count, tl,
        carquet_buffer_data_const(&idx_buf), carquet_buffer_size(&idx_buf),
        out, count);
    if (st != CARQUET_OK) TEST_FAIL("flba_roundtrip", "decode failed");

    if (memcmp(out, vals, sizeof(vals)) != 0) {
        TEST_FAIL("flba_roundtrip", "value mismatch");
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&idx_buf);
    TEST_PASS("flba_roundtrip");
    return 0;
}

/* ---- error / boundary paths -------------------------------------------- */
static int test_error_paths(void) {
    /* output_count == 0 is a no-op success even with NULL-ish inputs. */
    if (carquet_dictionary_decode_byte_array(NULL, 0, 0, NULL, 0, NULL, 0)
        != CARQUET_OK) {
        TEST_FAIL("error_paths", "zero output_count should be OK (ba)");
    }
    if (carquet_dictionary_decode_fixed_len_byte_array(NULL, 0, 0, 4, NULL, 0,
                                                       NULL, 0) != CARQUET_OK) {
        TEST_FAIL("error_paths", "zero output_count should be OK (flba)");
    }

    /* Encode a small BYTE_ARRAY dictionary to get a valid index stream. */
    carquet_byte_array_t input[] = { ba("x"), ba("yy"), ba("x") };
    carquet_buffer_t dict_buf, idx_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&idx_buf);
    if (carquet_dictionary_encode_byte_array(input, 3, &dict_buf, &idx_buf)
        != CARQUET_OK) {
        TEST_FAIL("error_paths", "setup encode failed");
    }

    carquet_byte_array_t out[3];

    /* dict_count <= 0 → DECODE. */
    if (carquet_dictionary_decode_byte_array(
            carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
            0, carquet_buffer_data_const(&idx_buf),
            carquet_buffer_size(&idx_buf), out, 3) != CARQUET_ERROR_DECODE) {
        TEST_FAIL("error_paths", "dict_count<=0 should be DECODE");
    }

    /* Truncated dictionary (claim more entries than the buffer holds) →
     * DECODE, no out-of-bounds read. */
    if (carquet_dictionary_decode_byte_array(
            carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
            99, carquet_buffer_data_const(&idx_buf),
            carquet_buffer_size(&idx_buf), out, 3) != CARQUET_ERROR_DECODE) {
        TEST_FAIL("error_paths", "over-claimed dict_count should be DECODE");
    }

    /* Empty index stream → DECODE. */
    if (carquet_dictionary_decode_byte_array(
            carquet_buffer_data_const(&dict_buf), carquet_buffer_size(&dict_buf),
            2, carquet_buffer_data_const(&idx_buf), 0, out, 3)
        != CARQUET_ERROR_DECODE) {
        TEST_FAIL("error_paths", "empty index stream should be DECODE");
    }

    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&idx_buf);

    /* FLBA: type_length <= 0 → INVALID_ARGUMENT. */
    uint8_t dummy_dict[4] = {0};
    uint8_t idx_stream[2] = {1 /*bit width*/, 0};
    uint8_t flba_out[4];
    if (carquet_dictionary_decode_fixed_len_byte_array(
            dummy_dict, sizeof(dummy_dict), 1, 0, idx_stream, sizeof(idx_stream),
            flba_out, 1) != CARQUET_ERROR_INVALID_ARGUMENT) {
        TEST_FAIL("error_paths", "flba type_length<=0 should be INVALID_ARGUMENT");
    }

    /* FLBA: dictionary too small for dict_count*type_length → DECODE. */
    if (carquet_dictionary_decode_fixed_len_byte_array(
            dummy_dict, sizeof(dummy_dict), 2 /*need 8 bytes*/, 4,
            idx_stream, sizeof(idx_stream), flba_out, 1)
        != CARQUET_ERROR_DECODE) {
        TEST_FAIL("error_paths", "flba short dict should be DECODE");
    }

    TEST_PASS("error_paths");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_byte_array_roundtrip();
    failures += test_flba_roundtrip();
    failures += test_error_paths();

    if (failures) {
        fprintf(stderr, "%d dictionary-decode test(s) failed\n", failures);
        return 1;
    }
    printf("All dictionary decode tests passed\n");
    return 0;
}
