/**
 * @file plain.c
 * @brief PLAIN encoding implementation
 */

#include "plain.h"
#include "core/endian.h"
#include <string.h>

extern void carquet_dispatch_unpack_bools(const uint8_t* input, uint8_t* output, int64_t count);
extern void carquet_dispatch_pack_bools(const uint8_t* input, uint8_t* output, int64_t count);

/* ============================================================================
 * PLAIN Decoding
 * ============================================================================
 */

int64_t carquet_decode_plain_boolean(
    const uint8_t* input,
    size_t input_size,
    uint8_t* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    /* Booleans are packed 8 per byte */
    size_t bytes_needed = ((size_t)count + 7) / 8;
    if (input_size < bytes_needed) {
        return -1;
    }

    carquet_dispatch_unpack_bools(input, output, count);

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_int32(
    const uint8_t* input,
    size_t input_size,
    int32_t* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * 4;
    if (input_size < bytes_needed) {
        return -1;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    /* Fast path: direct memory copy on little-endian systems */
    memcpy(output, input, bytes_needed);
#else
    for (int64_t i = 0; i < count; i++) {
        output[i] = carquet_read_i32_le(input + i * 4);
    }
#endif

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_int64(
    const uint8_t* input,
    size_t input_size,
    int64_t* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * 8;
    if (input_size < bytes_needed) {
        return -1;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(output, input, bytes_needed);
#else
    for (int64_t i = 0; i < count; i++) {
        output[i] = carquet_read_i64_le(input + i * 8);
    }
#endif

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_int96(
    const uint8_t* input,
    size_t input_size,
    carquet_int96_t* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * 12;
    if (input_size < bytes_needed) {
        return -1;
    }

    for (int64_t i = 0; i < count; i++) {
        const uint8_t* p = input + i * 12;
        output[i].value[0] = carquet_read_u32_le(p);
        output[i].value[1] = carquet_read_u32_le(p + 4);
        output[i].value[2] = carquet_read_u32_le(p + 8);
    }

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_float(
    const uint8_t* input,
    size_t input_size,
    float* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * 4;
    if (input_size < bytes_needed) {
        return -1;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(output, input, bytes_needed);
#else
    for (int64_t i = 0; i < count; i++) {
        output[i] = carquet_read_f32_le(input + i * 4);
    }
#endif

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_double(
    const uint8_t* input,
    size_t input_size,
    double* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * 8;
    if (input_size < bytes_needed) {
        return -1;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(output, input, bytes_needed);
#else
    for (int64_t i = 0; i < count; i++) {
        output[i] = carquet_read_f64_le(input + i * 8);
    }
#endif

    return (int64_t)bytes_needed;
}

int64_t carquet_decode_plain_byte_array(
    const uint8_t* input,
    size_t input_size,
    carquet_byte_array_t* output,
    int64_t count) {

    if (!input || !output || count < 0) {
        return -1;
    }

    size_t pos = 0;

    for (int64_t i = 0; i < count; i++) {
        /* Read 4-byte length prefix */
        if (pos + 4 > input_size) {
            return -1;
        }

        int32_t len = carquet_read_i32_le(input + pos);
        pos += 4;

        if (len < 0 || pos + (size_t)len > input_size) {
            return -1;
        }

        output[i].data = (uint8_t*)(input + pos);
        output[i].length = len;
        pos += (size_t)len;
    }

    return (int64_t)pos;
}

int64_t carquet_decode_plain_fixed_byte_array(
    const uint8_t* input,
    size_t input_size,
    uint8_t* output,
    int64_t count,
    int32_t fixed_len) {

    if (!input || !output || count < 0 || fixed_len <= 0) {
        return -1;
    }

    size_t bytes_needed = (size_t)count * (size_t)fixed_len;
    if (input_size < bytes_needed) {
        return -1;
    }

    memcpy(output, input, bytes_needed);
    return (int64_t)bytes_needed;
}

/* ============================================================================
 * PLAIN Encoding
 * ============================================================================
 */

carquet_status_t carquet_encode_plain_boolean(
    const uint8_t* input,
    int64_t count,
    carquet_buffer_t* output) {

    /* A new stream is an append at bit 0. */
    return carquet_encode_plain_boolean_append(input, count, 0, output);
}

carquet_status_t carquet_encode_plain_boolean_append(
    const uint8_t* input,
    int64_t count,
    int64_t bit_count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0 || bit_count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* The stream's last byte has `used` bits taken when it does not end on a
     * byte boundary; the first `head` values go into the rest of it. */
    int used = (int)(bit_count & 7);
    int64_t head = 0;
    if (used != 0) {
        if (output->size == 0) {
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        head = 8 - used;
        if (head > count) head = count;
    }

    /* Whatever is left starts on a byte boundary. Reserve its bytes before
     * touching the partial byte, so a failed allocation leaves the stream
     * exactly as it was. */
    int64_t tail = count - head;
    size_t tail_bytes = ((size_t)tail + 7) / 8;
    size_t partial = (used != 0) ? output->size - 1 : 0;
    uint8_t* dest = NULL;
    if (tail_bytes > 0) {
        dest = carquet_buffer_advance(output, tail_bytes);
        if (!dest) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
    }

    for (int64_t i = 0; i < head; i++) {
        if (input[i]) {
            output->data[partial] |= (uint8_t)(1u << (used + (int)i));
        }
    }

    if (tail_bytes > 0) {
        memset(dest, 0, tail_bytes);
        carquet_dispatch_pack_bools(input + head, dest, tail);
    }

    return CARQUET_OK;
}

carquet_status_t carquet_encode_plain_int32(
    const int32_t* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    size_t bytes_needed = (size_t)count * 4;

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    return carquet_buffer_append(output, input, bytes_needed);
#else
    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_buffer_append_u32_le(output, (uint32_t)input[i]);
        if (status != CARQUET_OK) return status;
    }
    return CARQUET_OK;
#endif
}

carquet_status_t carquet_encode_plain_int64(
    const int64_t* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    return carquet_buffer_append(output, input, (size_t)count * 8);
#else
    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_buffer_append_u64_le(output, (uint64_t)input[i]);
        if (status != CARQUET_OK) return status;
    }
    return CARQUET_OK;
#endif
}

carquet_status_t carquet_encode_plain_int96(
    const carquet_int96_t* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status;
        status = carquet_buffer_append_u32_le(output, input[i].value[0]);
        if (status != CARQUET_OK) return status;
        status = carquet_buffer_append_u32_le(output, input[i].value[1]);
        if (status != CARQUET_OK) return status;
        status = carquet_buffer_append_u32_le(output, input[i].value[2]);
        if (status != CARQUET_OK) return status;
    }

    return CARQUET_OK;
}

carquet_status_t carquet_encode_plain_float(
    const float* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    return carquet_buffer_append(output, input, (size_t)count * 4);
#else
    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_buffer_append_f32_le(output, input[i]);
        if (status != CARQUET_OK) return status;
    }
    return CARQUET_OK;
#endif
}

carquet_status_t carquet_encode_plain_double(
    const double* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    return carquet_buffer_append(output, input, (size_t)count * 8);
#else
    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_buffer_append_f64_le(output, input[i]);
        if (status != CARQUET_OK) return status;
    }
    return CARQUET_OK;
#endif
}

carquet_status_t carquet_encode_plain_byte_array(
    const carquet_byte_array_t* input,
    int64_t count,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (count == 0) return CARQUET_OK;

    /* Sum the total encoded size (4-byte length prefix + payload per value) so
     * the output can be grown once instead of two capacity checks per value.
     * Lengths are validated here, before anything is reserved: a negative
     * int32 length would otherwise wrap to a tiny reservation and a ~4GB
     * copy. `total` is 64-bit and each term is < 2^32, so the sum cannot
     * realistically overflow, but guard the reservation cast. */
    uint64_t total = 0;
    for (int64_t i = 0; i < count; i++) {
        if (input[i].length < 0) return CARQUET_ERROR_INVALID_ARGUMENT;
        total += 4u + (uint64_t)input[i].length;
    }
    if (total > SIZE_MAX - output->size) return CARQUET_ERROR_OUT_OF_MEMORY;

    uint8_t* p = carquet_buffer_advance(output, (size_t)total);
    if (!p) return CARQUET_ERROR_OUT_OF_MEMORY;

    for (int64_t i = 0; i < count; i++) {
        uint32_t len = (uint32_t)input[i].length;
        p[0] = (uint8_t)(len);
        p[1] = (uint8_t)(len >> 8);
        p[2] = (uint8_t)(len >> 16);
        p[3] = (uint8_t)(len >> 24);
        p += 4;
        if (len > 0) {
            if (input[i].data) {
                memcpy(p, input[i].data, len);
            } else {
                /* Malformed input (length > 0, NULL data): zero-fill rather than
                 * leak the reserved-but-unwritten bytes into the output. The
                 * writer API rejects such entries upstream with
                 * CARQUET_ERROR_INVALID_ARGUMENT (carquet_column_writer_write_batch),
                 * so this only guards direct callers of the encoder. */
                memset(p, 0, len);
            }
        }
        p += len;
    }

    return CARQUET_OK;
}

carquet_status_t carquet_encode_plain_fixed_byte_array(
    const uint8_t* input,
    int64_t count,
    int32_t fixed_len,
    carquet_buffer_t* output) {

    if (!input || !output || count < 0 || fixed_len <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    return carquet_buffer_append(output, input, (size_t)count * (size_t)fixed_len);
}

/* ============================================================================
 * Generic PLAIN Function
 * ============================================================================
 */

int64_t carquet_decode_plain(
    const uint8_t* input,
    size_t input_size,
    carquet_physical_type_t type,
    int32_t type_length,
    void* output,
    int64_t count) {

    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            return carquet_decode_plain_boolean(input, input_size,
                (uint8_t*)output, count);

        case CARQUET_PHYSICAL_INT32:
            return carquet_decode_plain_int32(input, input_size,
                (int32_t*)output, count);

        case CARQUET_PHYSICAL_INT64:
            return carquet_decode_plain_int64(input, input_size,
                (int64_t*)output, count);

        case CARQUET_PHYSICAL_INT96:
            return carquet_decode_plain_int96(input, input_size,
                (carquet_int96_t*)output, count);

        case CARQUET_PHYSICAL_FLOAT:
            return carquet_decode_plain_float(input, input_size,
                (float*)output, count);

        case CARQUET_PHYSICAL_DOUBLE:
            return carquet_decode_plain_double(input, input_size,
                (double*)output, count);

        case CARQUET_PHYSICAL_BYTE_ARRAY:
            return carquet_decode_plain_byte_array(input, input_size,
                (carquet_byte_array_t*)output, count);

        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return carquet_decode_plain_fixed_byte_array(input, input_size,
                (uint8_t*)output, count, type_length);

        default:
            return -1;
    }
}
