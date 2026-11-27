/**
 * @file delta_strings.c
 * @brief DELTA_BYTE_ARRAY encoding implementation
 *
 * This encoding uses incremental (prefix sharing) encoding for strings.
 * It stores:
 * 1. Prefix lengths (common prefix with previous string) using DELTA_BINARY_PACKED
 * 2. Suffix lengths using DELTA_BINARY_PACKED
 * 3. All suffix data concatenated
 *
 * This is particularly efficient for sorted string columns where
 * adjacent strings often share common prefixes.
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
 * Helper Functions
 * ============================================================================
 */

/**
 * Find the length of common prefix between two byte arrays.
 */
static int32_t common_prefix_length(
    const uint8_t* a, uint32_t a_len,
    const uint8_t* b, uint32_t b_len) {

    uint32_t min_len = a_len < b_len ? a_len : b_len;
    int32_t prefix_len = 0;

    for (uint32_t i = 0; i < min_len; i++) {
        if (a[i] != b[i]) break;
        prefix_len++;
    }

    return prefix_len;
}

/* ============================================================================
 * DELTA_BYTE_ARRAY Decoder
 * ============================================================================
 */

/**
 * Decode the two DELTA_BINARY_PACKED length sub-streams of a DELTA_BYTE_ARRAY
 * page once. Returns the prefix/suffix length arrays (heap, caller frees
 * both with carquet_mem_free), the number of header bytes they occupied and
 * the exact reconstruction size (sum of prefix + suffix over all values).
 * Every length is validated non-negative and the suffix payload is checked
 * to fit in `data_size`.
 */
carquet_status_t carquet_delta_strings_decode_lengths(
    const uint8_t* data,
    size_t data_size,
    int32_t num_values,
    int32_t** prefix_out,
    int32_t** suffix_out,
    size_t* header_bytes,
    size_t* total_size) {

    if (!data || !prefix_out || !suffix_out || !header_bytes || !total_size ||
        num_values <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    *prefix_out = NULL;
    *suffix_out = NULL;

    /* num_values comes from the (untrusted) page header; guard the multiply
     * so it cannot overflow size_t (only reachable where size_t is 32-bit). */
    if ((size_t)num_values > SIZE_MAX / sizeof(int32_t)) {
        return CARQUET_ERROR_DECODE;
    }
    int32_t* prefix_lengths = carquet_mem_malloc((size_t)num_values * sizeof(int32_t));
    int32_t* suffix_lengths = carquet_mem_malloc((size_t)num_values * sizeof(int32_t));
    if (!prefix_lengths || !suffix_lengths) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    size_t pos = 0;
    size_t consumed = 0;
    carquet_status_t status = carquet_delta_decode_int32(
        data, data_size, prefix_lengths, num_values, &consumed);
    if (status == CARQUET_OK) {
        pos += consumed;
        status = carquet_delta_decode_int32(
            data + pos, data_size - pos, suffix_lengths, num_values, &consumed);
    }
    if (status == CARQUET_OK) {
        pos += consumed;
        /* Validate as reconstruction will consume them: the first prefix is
         * 0, each prefix fits in the previous value, the suffixes fit in the
         * page. Then every value is at most total_suffix (<= data_size) bytes
         * long and the running total can be bounded without wrapping, even
         * with a 32-bit size_t. The total is capped at INT32_MAX: reachable
         * otherwise with tiny pages (bit-width-0 delta runs of maximal
         * prefixes), and the Arrow export cannot address more anyway. */
        size_t total = 0, total_suffix = 0, prev_len = 0;
        size_t suffix_avail = data_size - pos;
        for (int32_t i = 0; i < num_values; i++) {
            if (prefix_lengths[i] < 0 || suffix_lengths[i] < 0 ||
                (size_t)prefix_lengths[i] > prev_len ||
                (size_t)suffix_lengths[i] > suffix_avail - total_suffix) {
                status = CARQUET_ERROR_DECODE;
                break;
            }
            total_suffix += (size_t)suffix_lengths[i];
            prev_len = (size_t)prefix_lengths[i] + (size_t)suffix_lengths[i];
            if (prev_len > (size_t)INT32_MAX - total) {
                status = CARQUET_ERROR_DECODE;
                break;
            }
            total += prev_len;
        }
        if (status == CARQUET_OK) {
            *prefix_out = prefix_lengths;
            *suffix_out = suffix_lengths;
            *header_bytes = pos;
            *total_size = total;
            return CARQUET_OK;
        }
    }
    carquet_mem_free(prefix_lengths);
    carquet_mem_free(suffix_lengths);
    return status;
}

/**
 * Reconstruct the strings from lengths obtained with
 * carquet_delta_strings_decode_lengths(). `work_buffer` must hold at least
 * the `total_size` that call reported; every value points into it.
 */
carquet_status_t carquet_delta_strings_reconstruct(
    const uint8_t* data,
    size_t data_size,
    size_t header_bytes,
    const int32_t* prefix_lengths,
    const int32_t* suffix_lengths,
    carquet_byte_array_t* values,
    int32_t num_values,
    uint8_t* work_buffer,
    size_t work_buffer_size,
    size_t* bytes_consumed) {

    if (!data || !prefix_lengths || !suffix_lengths || num_values <= 0 || !values ||
        header_bytes > data_size || (!work_buffer && work_buffer_size > 0)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const uint8_t* suffix_data = data + header_bytes;
    size_t suffix_avail = data_size - header_bytes;
    size_t suffix_offset = 0;
    size_t work_offset = 0;
    uint8_t* prev_string = NULL;
    uint32_t prev_len = 0;

    for (int32_t i = 0; i < num_values; i++) {
        int32_t prefix_len = prefix_lengths[i];
        int32_t suffix_len = suffix_lengths[i];
        if (prefix_len < 0 || suffix_len < 0) return CARQUET_ERROR_DECODE;
        /* Sum in size_t so two near-INT32_MAX lengths cannot overflow. */
        size_t total_len = (size_t)prefix_len + (size_t)suffix_len;

        if (total_len > UINT32_MAX || total_len > work_buffer_size - work_offset) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if ((size_t)suffix_len > suffix_avail - suffix_offset) {
            return CARQUET_ERROR_DECODE;
        }

        uint8_t* dest = work_buffer + work_offset;

        /* Copy prefix from previous string */
        if (prefix_len > 0) {
            if (!prev_string || prefix_len > (int32_t)prev_len) {
                return CARQUET_ERROR_DECODE;
            }
            memcpy(dest, prev_string, (size_t)prefix_len);
        }

        /* Copy suffix from encoded data */
        if (suffix_len > 0) {
            memcpy(dest + prefix_len, suffix_data + suffix_offset, (size_t)suffix_len);
            suffix_offset += (size_t)suffix_len;
        }

        values[i].data = dest;
        values[i].length = (uint32_t)total_len;

        prev_string = dest;
        prev_len = (uint32_t)total_len;
        work_offset += total_len;
    }

    if (bytes_consumed) {
        *bytes_consumed = header_bytes + suffix_offset;
    }
    return CARQUET_OK;
}

/**
 * Decode DELTA_BYTE_ARRAY encoded data.
 *
 * @param data Input buffer containing encoded data
 * @param data_size Size of input buffer
 * @param values Output array of byte arrays (must be pre-allocated)
 * @param num_values Number of values to decode
 * @param work_buffer Work buffer for reconstructing strings
 * @param work_buffer_size Size of work buffer
 * @param bytes_consumed Output: number of input bytes consumed
 * @return Status code
 */
carquet_status_t carquet_delta_strings_decode(
    const uint8_t* data,
    size_t data_size,
    carquet_byte_array_t* values,
    int32_t num_values,
    uint8_t* work_buffer,
    size_t work_buffer_size,
    size_t* bytes_consumed) {

    if (!data || num_values < 0 || (num_values > 0 && !values)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Empty (all-null) page: consume both DELTA headers (prefix + suffix
     * length sub-streams) for correct byte accounting, then yield zero
     * values. */
    if (num_values == 0) {
        size_t c1 = 0, c2 = 0;
        carquet_status_t s = carquet_delta_decode_int32(
            data, data_size, NULL, 0, &c1);
        if (s != CARQUET_OK) return s;
        s = carquet_delta_decode_int32(
            data + c1, data_size - c1, NULL, 0, &c2);
        if (s != CARQUET_OK) return s;
        if (bytes_consumed) *bytes_consumed = c1 + c2;
        return CARQUET_OK;
    }

    int32_t* prefix_lengths = NULL;
    int32_t* suffix_lengths = NULL;
    size_t header_bytes = 0, total = 0;
    carquet_status_t status = carquet_delta_strings_decode_lengths(
        data, data_size, num_values, &prefix_lengths, &suffix_lengths,
        &header_bytes, &total);
    if (status != CARQUET_OK) return status;

    status = carquet_delta_strings_reconstruct(
        data, data_size, header_bytes, prefix_lengths, suffix_lengths,
        values, num_values, work_buffer, work_buffer_size, bytes_consumed);
    carquet_mem_free(prefix_lengths);
    carquet_mem_free(suffix_lengths);
    return status;
}

/**
 * Compute the exact work buffer size required to decode a DELTA_BYTE_ARRAY
 * page, without reconstructing the strings. Callers that go on to decode
 * should prefer carquet_delta_strings_decode_lengths() +
 * carquet_delta_strings_reconstruct(), which decode the length streams once.
 *
 * @param data Input buffer containing encoded data
 * @param data_size Size of input buffer
 * @param num_values Number of values to decode
 * @param required_size Output: exact work buffer size in bytes
 * @return Status code
 */
carquet_status_t carquet_delta_strings_decoded_size(
    const uint8_t* data,
    size_t data_size,
    int32_t num_values,
    size_t* required_size) {

    if (!data || !required_size || num_values <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int32_t* prefix_lengths = NULL;
    int32_t* suffix_lengths = NULL;
    size_t header_bytes = 0, total = 0;
    carquet_status_t status = carquet_delta_strings_decode_lengths(
        data, data_size, num_values, &prefix_lengths, &suffix_lengths,
        &header_bytes, &total);
    if (status != CARQUET_OK) return status;
    carquet_mem_free(prefix_lengths);
    carquet_mem_free(suffix_lengths);
    *required_size = total;
    return CARQUET_OK;
}

/* ============================================================================
 * DELTA_BYTE_ARRAY Encoder
 * ============================================================================
 */

/**
 * Encode byte arrays using DELTA_BYTE_ARRAY (incremental) encoding.
 *
 * @param values Input byte arrays to encode
 * @param num_values Number of values to encode
 * @param output Output buffer for encoded data
 * @return Status code
 */
carquet_status_t carquet_delta_strings_encode(
    const carquet_byte_array_t* values,
    int32_t num_values,
    carquet_buffer_t* output) {

    if (!output || num_values < 0 || (num_values > 0 && !values)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* An all-null (zero value) page still needs both DELTA headers
     * (prefix-length and suffix-length sub-streams) so the decoder doesn't
     * hit EOF; emit them with no trailing suffix bytes. */
    if (num_values == 0) {
        uint8_t header[64];
        size_t written = 0;
        carquet_status_t s = carquet_delta_encode_int32(
            NULL, 0, header, sizeof(header), &written);
        if (s != CARQUET_OK) return s;
        s = carquet_buffer_append(output, header, written);   /* prefix */
        if (s != CARQUET_OK) return s;
        return carquet_buffer_append(output, header, written); /* suffix */
    }

    /* Allocate arrays for prefix and suffix lengths */
    int32_t* prefix_lengths = carquet_mem_malloc(num_values * sizeof(int32_t));
    int32_t* suffix_lengths = carquet_mem_malloc(num_values * sizeof(int32_t));

    if (!prefix_lengths || !suffix_lengths) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    /* Calculate prefix and suffix lengths */
    const uint8_t* prev_data = NULL;
    uint32_t prev_len = 0;

    for (int32_t i = 0; i < num_values; i++) {
        if (i == 0) {
            prefix_lengths[i] = 0;
            suffix_lengths[i] = (int32_t)values[i].length;
        } else {
            int32_t prefix_len = common_prefix_length(
                prev_data, prev_len,
                values[i].data, values[i].length);
            prefix_lengths[i] = prefix_len;
            suffix_lengths[i] = (int32_t)(values[i].length - prefix_len);
        }

        prev_data = values[i].data;
        prev_len = values[i].length;
    }

    /* Encode prefix lengths */
    size_t delta_capacity = (size_t)num_values * 10 + 100;
    uint8_t* delta_buffer = carquet_mem_malloc(delta_capacity);
    if (!delta_buffer) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    size_t bytes_written = 0;
    carquet_status_t status = carquet_delta_encode_int32(
        prefix_lengths, num_values, delta_buffer, delta_capacity, &bytes_written);

    if (status != CARQUET_OK) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        carquet_mem_free(delta_buffer);
        return status;
    }

    status = carquet_buffer_append(output, delta_buffer, bytes_written);
    if (status != CARQUET_OK) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        carquet_mem_free(delta_buffer);
        return status;
    }

    /* Encode suffix lengths */
    status = carquet_delta_encode_int32(
        suffix_lengths, num_values, delta_buffer, delta_capacity, &bytes_written);

    if (status != CARQUET_OK) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        carquet_mem_free(delta_buffer);
        return status;
    }

    status = carquet_buffer_append(output, delta_buffer, bytes_written);
    carquet_mem_free(delta_buffer);

    if (status != CARQUET_OK) {
        carquet_mem_free(prefix_lengths);
        carquet_mem_free(suffix_lengths);
        return status;
    }

    /* Write suffix data (prefix lengths were computed above: no second
     * common-prefix scan per value). */
    prev_len = 0;
    for (int32_t i = 0; i < num_values; i++) {
        int32_t prefix_len = prefix_lengths[i];
        int32_t suffix_len = suffix_lengths[i];

        if (suffix_len > 0 && values[i].data) {
            status = carquet_buffer_append(output,
                values[i].data + prefix_len, suffix_len);
            if (status != CARQUET_OK) {
                carquet_mem_free(prefix_lengths);
                carquet_mem_free(suffix_lengths);
                return status;
            }
        }
    }

    carquet_mem_free(prefix_lengths);
                carquet_mem_free(suffix_lengths);
    return CARQUET_OK;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

/**
 * Estimate work buffer size needed for decoding.
 *
 * @param values Array of byte arrays (with only length information needed)
 * @param num_values Number of values
 * @return Required work buffer size
 */
size_t carquet_delta_strings_work_buffer_size(
    const carquet_byte_array_t* values,
    int32_t num_values) {

    if (!values || num_values <= 0) {
        return 0;
    }

    size_t total = 0;
    for (int32_t i = 0; i < num_values; i++) {
        total += values[i].length;
    }

    return total;
}

/**
 * Estimate maximum encoded size for DELTA_BYTE_ARRAY.
 *
 * @param values Input byte arrays
 * @param num_values Number of values
 * @return Estimated maximum encoded size
 */
size_t carquet_delta_strings_max_encoded_size(
    const carquet_byte_array_t* values,
    int32_t num_values) {

    if (!values || num_values <= 0) {
        return 0;
    }

    /* Sum of all string lengths */
    size_t total_size = 0;
    for (int32_t i = 0; i < num_values; i++) {
        total_size += values[i].length;
    }

    /* Overhead for two delta-encoded integer arrays (prefix and suffix lengths) */
    size_t overhead = (size_t)num_values * 10 + 200;

    return total_size + overhead;
}
