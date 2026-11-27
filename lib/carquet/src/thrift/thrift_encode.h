/**
 * @file thrift_encode.h
 * @brief Thrift Compact Protocol encoder
 *
 * Encodes data structures using the Thrift Compact Protocol for
 * writing Parquet file metadata.
 */

#ifndef CARQUET_THRIFT_ENCODE_H
#define CARQUET_THRIFT_ENCODE_H

#include <carquet/error.h>
#include "core/buffer.h"
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Thrift Encoder State
 * ============================================================================
 */

#define THRIFT_ENCODER_MAX_NESTING 32

typedef struct thrift_encoder {
    carquet_buffer_t* buffer;   /* Output buffer */

    /* Field ID tracking for delta encoding */
    int16_t last_field_id[THRIFT_ENCODER_MAX_NESTING];
    int nesting_level;

    /* Error state */
    carquet_status_t status;
} thrift_encoder_t;

/* ============================================================================
 * Encoder Lifecycle
 * ============================================================================
 */

/**
 * Initialize an encoder with an output buffer.
 */
void thrift_encoder_init(thrift_encoder_t* enc, carquet_buffer_t* buffer);

/**
 * Check if encoder is in error state.
 */
static inline bool thrift_encoder_has_error(const thrift_encoder_t* enc) {
    return enc->status != CARQUET_OK;
}

/**
 * Get bytes written so far.
 */
static inline size_t thrift_encoder_size(const thrift_encoder_t* enc) {
    return carquet_buffer_size(enc->buffer);
}

/* ============================================================================
 * Primitive Writing
 * ============================================================================
 */

/**
 * Write a single byte.
 */
void thrift_write_byte(thrift_encoder_t* enc, int8_t value);

/**
 * Write a 16-bit signed integer.
 */
void thrift_write_i16(thrift_encoder_t* enc, int16_t value);

/**
 * Write a 32-bit signed integer.
 */
void thrift_write_i32(thrift_encoder_t* enc, int32_t value);

/**
 * Write a 64-bit signed integer.
 */
void thrift_write_i64(thrift_encoder_t* enc, int64_t value);

/**
 * Write a double.
 */
void thrift_write_double(thrift_encoder_t* enc, double value);

/**
 * Write a boolean.
 */
void thrift_write_bool(thrift_encoder_t* enc, bool value);

/**
 * Write binary data.
 */
void thrift_write_binary(thrift_encoder_t* enc, const uint8_t* data, int32_t length);

/**
 * Write a null-terminated string.
 */
void thrift_write_string(thrift_encoder_t* enc, const char* str);

/**
 * Write a UUID.
 */
void thrift_write_uuid(thrift_encoder_t* enc, const uint8_t uuid[16]);

/* ============================================================================
 * Struct Writing
 * ============================================================================
 */

/**
 * Begin writing a struct.
 */
void thrift_write_struct_begin(thrift_encoder_t* enc);

/**
 * End writing a struct.
 */
void thrift_write_struct_end(thrift_encoder_t* enc);

/**
 * Write a field header.
 *
 * @param enc Encoder
 * @param type Field type
 * @param field_id Field ID
 */
void thrift_write_field_header(thrift_encoder_t* enc, int type, int16_t field_id);

/**
 * Write a field stop marker.
 */
void thrift_write_field_stop(thrift_encoder_t* enc);

/* Convenience macros for writing fields */
#define THRIFT_WRITE_FIELD_BYTE(enc, id, val) do { \
    thrift_write_field_header(enc, 3, id); \
    thrift_write_byte(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_I16(enc, id, val) do { \
    thrift_write_field_header(enc, 4, id); \
    thrift_write_i16(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_I32(enc, id, val) do { \
    thrift_write_field_header(enc, 5, id); \
    thrift_write_i32(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_I64(enc, id, val) do { \
    thrift_write_field_header(enc, 6, id); \
    thrift_write_i64(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_DOUBLE(enc, id, val) do { \
    thrift_write_field_header(enc, 7, id); \
    thrift_write_double(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_BOOL(enc, id, val) do { \
    thrift_write_field_header(enc, (val) ? 1 : 2, id); \
} while(0)

#define THRIFT_WRITE_FIELD_STRING(enc, id, val) do { \
    thrift_write_field_header(enc, 8, id); \
    thrift_write_string(enc, val); \
} while(0)

#define THRIFT_WRITE_FIELD_BINARY(enc, id, data, len) do { \
    thrift_write_field_header(enc, 8, id); \
    thrift_write_binary(enc, data, len); \
} while(0)

/* ============================================================================
 * Container Writing
 * ============================================================================
 */

/**
 * Begin writing a list.
 *
 * @param enc Encoder
 * @param elem_type Element type
 * @param count Number of elements
 */
void thrift_write_list_begin(thrift_encoder_t* enc, int elem_type, int32_t count);

/**
 * Begin writing a set.
 */
void thrift_write_set_begin(thrift_encoder_t* enc, int elem_type, int32_t count);

/**
 * Begin writing a map.
 *
 * @param enc Encoder
 * @param key_type Key type
 * @param value_type Value type
 * @param count Number of key-value pairs
 */
void thrift_write_map_begin(thrift_encoder_t* enc,
                             int key_type, int value_type, int32_t count);

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

/**
 * Write a varint (unsigned).
 */
void thrift_write_varint(thrift_encoder_t* enc, uint64_t value);

/**
 * Write a zigzag-encoded varint (signed).
 */
void thrift_write_zigzag(thrift_encoder_t* enc, int64_t value);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_THRIFT_ENCODE_H */
