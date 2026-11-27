/**
 * @file thrift_decode.h
 * @brief Thrift Compact Protocol decoder
 *
 * Parquet uses the Thrift Compact Protocol for metadata serialization.
 * This is a minimal implementation supporting only the features needed
 * for parsing Parquet files.
 *
 * Compact Protocol specification:
 * https://github.com/apache/thrift/blob/master/doc/specs/thrift-compact-protocol.md
 */

#ifndef CARQUET_THRIFT_DECODE_H
#define CARQUET_THRIFT_DECODE_H

#include <carquet/error.h>
#include "core/buffer.h"
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Thrift Type Constants
 * ============================================================================
 */

/**
 * Thrift wire types (compact protocol).
 * These are the types as they appear on the wire, not the Thrift type IDs.
 */
typedef enum thrift_type {
    THRIFT_TYPE_STOP = 0,       /* End of struct */
    THRIFT_TYPE_TRUE = 1,       /* Boolean true */
    THRIFT_TYPE_FALSE = 2,      /* Boolean false */
    THRIFT_TYPE_BYTE = 3,       /* Signed 8-bit integer */
    THRIFT_TYPE_I16 = 4,        /* Signed 16-bit integer */
    THRIFT_TYPE_I32 = 5,        /* Signed 32-bit integer */
    THRIFT_TYPE_I64 = 6,        /* Signed 64-bit integer */
    THRIFT_TYPE_DOUBLE = 7,     /* 64-bit floating point */
    THRIFT_TYPE_BINARY = 8,     /* Binary/string data */
    THRIFT_TYPE_LIST = 9,       /* List container */
    THRIFT_TYPE_SET = 10,       /* Set container */
    THRIFT_TYPE_MAP = 11,       /* Map container */
    THRIFT_TYPE_STRUCT = 12,    /* Struct/nested structure */
    THRIFT_TYPE_UUID = 13,      /* UUID (16 bytes) */
} thrift_type_t;

/* ============================================================================
 * Thrift Decoder State
 * ============================================================================
 */

/**
 * Maximum nesting depth for structs.
 */
#define THRIFT_MAX_NESTING 32

/**
 * Thrift decoder state.
 */
typedef struct thrift_decoder {
    carquet_buffer_reader_t reader;

    /* Field ID tracking for delta encoding */
    int16_t last_field_id[THRIFT_MAX_NESTING];
    int nesting_level;

    /* Boolean field tracking */
    bool bool_pending;
    bool bool_value;

    /* Error state */
    carquet_status_t status;
    char error_message[128];
} thrift_decoder_t;

/* ============================================================================
 * Decoder Lifecycle
 * ============================================================================
 */

/**
 * Initialize a decoder from a buffer.
 */
void thrift_decoder_init(thrift_decoder_t* dec, const uint8_t* data, size_t size);

/**
 * Initialize a decoder from a buffer reader.
 */
void thrift_decoder_init_reader(thrift_decoder_t* dec,
                                 const carquet_buffer_reader_t* reader);

/**
 * Check if decoder is in error state.
 */
static inline bool thrift_decoder_has_error(const thrift_decoder_t* dec) {
    return dec->status != CARQUET_OK;
}

/**
 * Get remaining bytes.
 */
static inline size_t thrift_decoder_remaining(const thrift_decoder_t* dec) {
    return carquet_buffer_reader_remaining(&dec->reader);
}

/* ============================================================================
 * Primitive Reading
 * ============================================================================
 */

/**
 * Read a single byte.
 */
int8_t thrift_read_byte(thrift_decoder_t* dec);

/**
 * Read a 16-bit signed integer (zigzag + varint).
 */
int16_t thrift_read_i16(thrift_decoder_t* dec);

/**
 * Read a 32-bit signed integer (zigzag + varint).
 */
int32_t thrift_read_i32(thrift_decoder_t* dec);

/**
 * Read a 64-bit signed integer (zigzag + varint).
 */
int64_t thrift_read_i64(thrift_decoder_t* dec);

/**
 * Read a double (8 bytes, IEEE 754).
 */
double thrift_read_double(thrift_decoder_t* dec);

/**
 * Read a boolean.
 */
bool thrift_read_bool(thrift_decoder_t* dec);

/**
 * Read a string/binary length and return pointer to data.
 * Does not copy the data - returns pointer into the buffer.
 *
 * @param dec Decoder
 * @param length Output: length of the binary data
 * @return Pointer to binary data, or NULL on error
 */
const uint8_t* thrift_read_binary(thrift_decoder_t* dec, int32_t* length);

/**
 * Read a string into a newly allocated buffer.
 * Caller must free the returned string.
 */
char* thrift_read_string_alloc(thrift_decoder_t* dec);

/**
 * Read a UUID (16 bytes).
 */
void thrift_read_uuid(thrift_decoder_t* dec, uint8_t uuid[16]);

/* ============================================================================
 * Struct Reading
 * ============================================================================
 */

/**
 * Begin reading a struct.
 * Must be paired with thrift_read_struct_end().
 */
void thrift_read_struct_begin(thrift_decoder_t* dec);

/**
 * End reading a struct.
 */
void thrift_read_struct_end(thrift_decoder_t* dec);

/**
 * Read a field header.
 *
 * @param dec Decoder
 * @param type Output: field type (or THRIFT_TYPE_STOP if no more fields)
 * @param field_id Output: field ID
 * @return true if a field was read, false if STOP encountered
 */
bool thrift_read_field_begin(thrift_decoder_t* dec,
                              thrift_type_t* type,
                              int16_t* field_id);

/**
 * Skip a field value based on its type.
 */
void thrift_skip_field(thrift_decoder_t* dec, thrift_type_t type);

/* ============================================================================
 * Container Reading
 * ============================================================================
 */

/**
 * Begin reading a list.
 *
 * @param dec Decoder
 * @param elem_type Output: element type
 * @param count Output: number of elements
 */
void thrift_read_list_begin(thrift_decoder_t* dec,
                             thrift_type_t* elem_type,
                             int32_t* count);

/**
 * Begin reading a set (same as list).
 */
void thrift_read_set_begin(thrift_decoder_t* dec,
                            thrift_type_t* elem_type,
                            int32_t* count);

/**
 * Begin reading a map.
 *
 * @param dec Decoder
 * @param key_type Output: key type
 * @param value_type Output: value type
 * @param count Output: number of key-value pairs
 */
void thrift_read_map_begin(thrift_decoder_t* dec,
                            thrift_type_t* key_type,
                            thrift_type_t* value_type,
                            int32_t* count);

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

/**
 * Read a varint (unsigned).
 */
uint64_t thrift_read_varint(thrift_decoder_t* dec);

/**
 * Read a zigzag-encoded varint (signed).
 */
int64_t thrift_read_zigzag(thrift_decoder_t* dec);

/**
 * Skip a value of the given type.
 */
void thrift_skip(thrift_decoder_t* dec, thrift_type_t type);

/**
 * Get type name for debugging.
 */
const char* thrift_type_name(thrift_type_t type);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_THRIFT_DECODE_H */
