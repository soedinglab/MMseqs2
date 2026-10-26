/**
 * @file delta_length.c
 * @brief DELTA_LENGTH_BYTE_ARRAY encoding implementation
 *
 * This encoding is used for variable-length byte arrays (strings).
 * It stores:
 * 1. The lengths of all byte arrays using DELTA_BINARY_PACKED encoding
 * 2. All the byte array data concatenated together
 *
 * Reference: https://parquet.apache.org/docs/file-format/data-pages/encodings/
 */

#include "core/allocator.h"
#include <carquet/error.h>
#include <carquet/types.h>
#include "core/buffer.h"
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Forward declaration from delta.c */
extern carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data,
    size_t data_size,
    int32_t* values,
    int32_t num_values,
    size_t* bytes_consumed);

extern carquet_status_t carquet_delta_encode_int32(
    const int32_t* values,
    int32_t num_values,
    uint8_t* data,
    size_t data_capacity,
    size_t* bytes_written);

/* ============================================================================
 * DELTA_LENGTH_BYTE_ARRAY Decoder
 * ============================================================================
 */

/**
 * Decode DELTA_LENGTH_BYTE_ARRAY encoded data.
 *
 * @param data Input buffer containing encoded data
 * @param data_size Size of input buffer
 * @param values Output array of byte arrays
 * @param num_values Number of values to decode
 * @param bytes_consumed Output: number of input bytes consumed
 * @return Status code
 */
carquet_status_t carquet_delta_length_decode(
    const uint8_t* data,
    size_t data_size,
    carquet_byte_array_t* values,
    int32_t num_values,
    size_t* bytes_consumed) {

    if (!data || num_values < 0 || (num_values > 0 && !values)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Empty (all-null) page: still consume the DELTA lengths header so the
     * caller's byte accounting stays correct, then yield zero values. */
    if (num_values == 0) {
        size_t hdr_consumed = 0;
        carquet_status_t s = carquet_delta_decode_int32(
            data, data_size, NULL, 0, &hdr_consumed);
        if (s != CARQUET_OK) return s;
        if (bytes_consumed) *bytes_consumed = hdr_consumed;
        return CARQUET_OK;
    }

    /* Allocate buffer for lengths. num_values comes from the (untrusted) page
     * header; guard the multiply so it cannot overflow size_t and yield an
     * undersized buffer (only reachable where size_t is 32-bit). */
    if ((size_t)num_values > SIZE_MAX / sizeof(int32_t)) {
        return CARQUET_ERROR_DECODE;
    }
    int32_t* lengths = carquet_mem_malloc((size_t)num_values * sizeof(int32_t));
    if (!lengths) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    /* Decode lengths using delta encoding */
    size_t lengths_consumed = 0;
    carquet_status_t status = carquet_delta_decode_int32(
        data, data_size, lengths, num_values, &lengths_consumed);

    if (status != CARQUET_OK) {
        carquet_mem_free(lengths);
        return status;
    }

    /* Calculate total data size and validate */
    size_t total_data_size = 0;
    for (int32_t i = 0; i < num_values; i++) {
        if (lengths[i] < 0) {
            carquet_mem_free(lengths);
            return CARQUET_ERROR_DECODE;
        }
        total_data_size += (size_t)lengths[i];
    }

    /* Check that we have enough data */
    if (lengths_consumed + total_data_size > data_size) {
        carquet_mem_free(lengths);
        return CARQUET_ERROR_DECODE;
    }

    /* Extract byte arrays from concatenated data */
    const uint8_t* byte_data = data + lengths_consumed;
    size_t offset = 0;

    for (int32_t i = 0; i < num_values; i++) {
        values[i].length = (uint32_t)lengths[i];
        /* Cast away const - the data is for reading only */
        values[i].data = (uint8_t*)(byte_data + offset);
        offset += lengths[i];
    }

    carquet_mem_free(lengths);

    if (bytes_consumed) {
        *bytes_consumed = lengths_consumed + total_data_size;
    }

    return CARQUET_OK;
}

/* ============================================================================
 * DELTA_LENGTH_BYTE_ARRAY Encoder
 * ============================================================================
 */

/**
 * Encode byte arrays using DELTA_LENGTH_BYTE_ARRAY encoding.
 *
 * @param values Input byte arrays to encode
 * @param num_values Number of values to encode
 * @param output Output buffer for encoded data
 * @return Status code
 */
carquet_status_t carquet_delta_length_encode(
    const carquet_byte_array_t* values,
    int32_t num_values,
    carquet_buffer_t* output) {

    if (!output || num_values < 0 || (num_values > 0 && !values)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* An all-null (zero value) page still needs the DELTA-encoded lengths
     * sub-stream header so the decoder doesn't hit EOF; emit it with no
     * trailing byte data. */
    if (num_values == 0) {
        uint8_t header[64];
        size_t written = 0;
        carquet_status_t s = carquet_delta_encode_int32(
            NULL, 0, header, sizeof(header), &written);
        if (s != CARQUET_OK) return s;
        return carquet_buffer_append(output, header, written);
    }

    /* Extract lengths */
    int32_t* lengths = carquet_mem_malloc(num_values * sizeof(int32_t));
    if (!lengths) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < num_values; i++) {
        if (values[i].length < 0 || (!values[i].data && values[i].length > 0)) {
            carquet_mem_free(lengths);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        lengths[i] = (int32_t)values[i].length;
    }

    /* Encode lengths using delta encoding */
    /* Estimate max size for delta encoding (generous estimate) */
    size_t lengths_capacity = (size_t)num_values * 10 + 100;
    uint8_t* lengths_buffer = carquet_mem_malloc(lengths_capacity);
    if (!lengths_buffer) {
        carquet_mem_free(lengths);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    size_t lengths_written = 0;
    carquet_status_t status = carquet_delta_encode_int32(
        lengths, num_values, lengths_buffer, lengths_capacity, &lengths_written);

    carquet_mem_free(lengths);

    if (status != CARQUET_OK) {
        carquet_mem_free(lengths_buffer);
        return status;
    }

    /* Write encoded lengths to output */
    status = carquet_buffer_append(output, lengths_buffer, lengths_written);
    carquet_mem_free(lengths_buffer);

    if (status != CARQUET_OK) {
        return status;
    }

    /* Write concatenated byte array data */
    for (int32_t i = 0; i < num_values; i++) {
        if (values[i].length > 0 && values[i].data) {
            status = carquet_buffer_append(output, values[i].data, values[i].length);
            if (status != CARQUET_OK) {
                return status;
            }
        }
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

/**
 * Estimate the maximum encoded size for DELTA_LENGTH_BYTE_ARRAY.
 *
 * @param values Input byte arrays
 * @param num_values Number of values
 * @return Estimated maximum encoded size
 */
size_t carquet_delta_length_max_encoded_size(
    const carquet_byte_array_t* values,
    int32_t num_values) {

    if (!values || num_values <= 0) {
        return 0;
    }

    /* Sum of all byte array lengths */
    size_t total_data_size = 0;
    for (int32_t i = 0; i < num_values; i++) {
        total_data_size += values[i].length;
    }

    /* Delta encoding overhead for lengths (very conservative estimate) */
    size_t lengths_overhead = (size_t)num_values * 5 + 100;

    return total_data_size + lengths_overhead;
}
