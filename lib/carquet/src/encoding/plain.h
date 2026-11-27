/**
 * @file plain.h
 * @brief PLAIN encoding for Parquet
 *
 * PLAIN encoding stores values directly without any special encoding.
 * It's the simplest encoding and serves as a fallback.
 */

#ifndef CARQUET_ENCODING_PLAIN_H
#define CARQUET_ENCODING_PLAIN_H

#include <carquet/types.h>
#include <carquet/error.h>
#include "core/buffer.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * PLAIN Decoding
 * ============================================================================
 */

/**
 * Decode PLAIN encoded booleans.
 * Booleans are packed 8 per byte, LSB first.
 *
 * @param input Input data
 * @param input_size Size of input data
 * @param output Output boolean array (as uint8_t, 0 or 1)
 * @param count Number of values to decode
 * @return Number of bytes consumed, or -1 on error
 */
int64_t carquet_decode_plain_boolean(
    const uint8_t* input,
    size_t input_size,
    uint8_t* output,
    int64_t count);

/**
 * Decode PLAIN encoded 32-bit integers.
 *
 * @param input Input data
 * @param input_size Size of input data
 * @param output Output array
 * @param count Number of values to decode
 * @return Number of bytes consumed, or -1 on error
 */
int64_t carquet_decode_plain_int32(
    const uint8_t* input,
    size_t input_size,
    int32_t* output,
    int64_t count);

/**
 * Decode PLAIN encoded 64-bit integers.
 */
int64_t carquet_decode_plain_int64(
    const uint8_t* input,
    size_t input_size,
    int64_t* output,
    int64_t count);

/**
 * Decode PLAIN encoded INT96 values.
 */
int64_t carquet_decode_plain_int96(
    const uint8_t* input,
    size_t input_size,
    carquet_int96_t* output,
    int64_t count);

/**
 * Decode PLAIN encoded floats.
 */
int64_t carquet_decode_plain_float(
    const uint8_t* input,
    size_t input_size,
    float* output,
    int64_t count);

/**
 * Decode PLAIN encoded doubles.
 */
int64_t carquet_decode_plain_double(
    const uint8_t* input,
    size_t input_size,
    double* output,
    int64_t count);

/**
 * Decode PLAIN encoded byte arrays.
 * Each value is prefixed with a 4-byte little-endian length.
 *
 * @param input Input data
 * @param input_size Size of input data
 * @param output Output array of byte array structs
 * @param count Number of values to decode
 * @return Number of bytes consumed, or -1 on error
 */
int64_t carquet_decode_plain_byte_array(
    const uint8_t* input,
    size_t input_size,
    carquet_byte_array_t* output,
    int64_t count);

/**
 * Decode PLAIN encoded fixed-length byte arrays.
 *
 * @param input Input data
 * @param input_size Size of input data
 * @param output Output buffer (must be count * fixed_len bytes)
 * @param count Number of values to decode
 * @param fixed_len Length of each fixed array
 * @return Number of bytes consumed, or -1 on error
 */
int64_t carquet_decode_plain_fixed_byte_array(
    const uint8_t* input,
    size_t input_size,
    uint8_t* output,
    int64_t count,
    int32_t fixed_len);

/* ============================================================================
 * PLAIN Encoding
 * ============================================================================
 */

/**
 * Encode booleans using PLAIN encoding.
 *
 * @param input Input boolean array (0 or non-0)
 * @param count Number of values
 * @param output Output buffer
 * @return Status code
 */
carquet_status_t carquet_encode_plain_boolean(
    const uint8_t* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Append booleans to a PLAIN bit stream that is already under way.
 *
 * PLAIN packs a page's booleans LSB-first into one continuous stream, so
 * values added by a later call continue in the earlier call's trailing
 * partial byte instead of starting a new one. carquet_encode_plain_boolean()
 * always starts on a fresh byte, i.e. it is this with a bit_count of 0.
 *
 * @param input Input boolean array (0 or non-0)
 * @param count Number of values to append
 * @param bit_count Number of values the stream already holds; they must be
 *                  the last ceil(bit_count / 8) bytes of @p output
 * @param output Output buffer ending with the stream
 * @return Status code. On failure the stream is left unchanged.
 */
carquet_status_t carquet_encode_plain_boolean_append(
    const uint8_t* input,
    int64_t count,
    int64_t bit_count,
    carquet_buffer_t* output);

/**
 * Encode 32-bit integers using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_int32(
    const int32_t* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode 64-bit integers using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_int64(
    const int64_t* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode INT96 values using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_int96(
    const carquet_int96_t* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode floats using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_float(
    const float* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode doubles using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_double(
    const double* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode byte arrays using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_byte_array(
    const carquet_byte_array_t* input,
    int64_t count,
    carquet_buffer_t* output);

/**
 * Encode fixed-length byte arrays using PLAIN encoding.
 */
carquet_status_t carquet_encode_plain_fixed_byte_array(
    const uint8_t* input,
    int64_t count,
    int32_t fixed_len,
    carquet_buffer_t* output);

/* ============================================================================
 * Generic PLAIN Functions
 * ============================================================================
 */

/**
 * Decode values based on physical type.
 *
 * @param input Input data
 * @param input_size Size of input data
 * @param type Physical type
 * @param type_length Type length (for fixed arrays)
 * @param output Output buffer
 * @param count Number of values to decode
 * @return Number of bytes consumed, or -1 on error
 */
int64_t carquet_decode_plain(
    const uint8_t* input,
    size_t input_size,
    carquet_physical_type_t type,
    int32_t type_length,
    void* output,
    int64_t count);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_ENCODING_PLAIN_H */
