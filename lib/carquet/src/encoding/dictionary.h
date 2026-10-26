/**
 * @file dictionary.h
 * @brief Dictionary encoding module — decode API declarations.
 *
 * Declares the module-level dictionary *decode* entry points so callers do
 * not have to redeclare them as local externs. Each function takes a decoded
 * dictionary page (PLAIN-encoded values) plus the RLE-bit-packed index stream
 * (a leading bit-width byte followed by the hybrid RLE payload) and
 * materialises `output_count` values. Indices are validated against
 * `dict_count` and, for the variable/fixed-width byte types, payloads are
 * validated against `dict_size`.
 *
 * The encode side is declared inline at its call sites (writer/tests); this
 * header exists to make the decode side symmetric and discoverable.
 */

#ifndef CARQUET_ENCODING_DICTIONARY_H
#define CARQUET_ENCODING_DICTIONARY_H

#include <carquet/error.h>
#include <carquet/types.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    int32_t* output, int64_t output_count);

carquet_status_t carquet_dictionary_decode_int64(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    int64_t* output, int64_t output_count);

carquet_status_t carquet_dictionary_decode_float(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    float* output, int64_t output_count);

carquet_status_t carquet_dictionary_decode_double(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    double* output, int64_t output_count);

/**
 * Decode a BYTE_ARRAY dictionary. Each output value's `data` pointer aliases
 * into `dict_data`, so the dictionary buffer must outlive `output`.
 */
carquet_status_t carquet_dictionary_decode_byte_array(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    carquet_byte_array_t* output, int64_t output_count);

/**
 * Decode a FIXED_LEN_BYTE_ARRAY dictionary into a contiguous buffer of
 * `output_count * type_length` bytes.
 */
carquet_status_t carquet_dictionary_decode_fixed_len_byte_array(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    int32_t type_length,
    const uint8_t* indices_data, size_t indices_size,
    uint8_t* output, int64_t output_count);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_ENCODING_DICTIONARY_H */
