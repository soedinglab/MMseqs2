/**
 * @file types.h
 * @brief Parquet physical and logical type definitions
 *
 * This header defines all Parquet data types according to the Apache Parquet
 * specification. Types are organized into physical types (storage format) and
 * logical types (semantic interpretation).
 */

#ifndef CARQUET_TYPES_H
#define CARQUET_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Physical Types
 * ============================================================================
 * Physical types represent how data is stored on disk. Parquet supports a
 * limited set of physical types to keep the format simple.
 */

typedef enum carquet_physical_type {
    CARQUET_PHYSICAL_BOOLEAN = 0,
    CARQUET_PHYSICAL_INT32 = 1,
    CARQUET_PHYSICAL_INT64 = 2,
    CARQUET_PHYSICAL_INT96 = 3,        /* Deprecated, used for timestamps */
    CARQUET_PHYSICAL_FLOAT = 4,
    CARQUET_PHYSICAL_DOUBLE = 5,
    CARQUET_PHYSICAL_BYTE_ARRAY = 6,
    CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY = 7,
} carquet_physical_type_t;

/* ============================================================================
 * Logical Types (ConvertedType - legacy)
 * ============================================================================
 * Legacy converted types for backwards compatibility.
 */

typedef enum carquet_converted_type {
    CARQUET_CONVERTED_NONE = -1,
    CARQUET_CONVERTED_UTF8 = 0,
    CARQUET_CONVERTED_MAP = 1,
    CARQUET_CONVERTED_MAP_KEY_VALUE = 2,
    CARQUET_CONVERTED_LIST = 3,
    CARQUET_CONVERTED_ENUM = 4,
    CARQUET_CONVERTED_DECIMAL = 5,
    CARQUET_CONVERTED_DATE = 6,
    CARQUET_CONVERTED_TIME_MILLIS = 7,
    CARQUET_CONVERTED_TIME_MICROS = 8,
    CARQUET_CONVERTED_TIMESTAMP_MILLIS = 9,
    CARQUET_CONVERTED_TIMESTAMP_MICROS = 10,
    CARQUET_CONVERTED_UINT_8 = 11,
    CARQUET_CONVERTED_UINT_16 = 12,
    CARQUET_CONVERTED_UINT_32 = 13,
    CARQUET_CONVERTED_UINT_64 = 14,
    CARQUET_CONVERTED_INT_8 = 15,
    CARQUET_CONVERTED_INT_16 = 16,
    CARQUET_CONVERTED_INT_32 = 17,
    CARQUET_CONVERTED_INT_64 = 18,
    CARQUET_CONVERTED_JSON = 19,
    CARQUET_CONVERTED_BSON = 20,
    CARQUET_CONVERTED_INTERVAL = 21,
} carquet_converted_type_t;

/* ============================================================================
 * Logical Types (Modern)
 * ============================================================================
 * Modern logical type system with more detailed type information.
 */

typedef enum carquet_logical_type_id {
    CARQUET_LOGICAL_UNKNOWN = 0,
    CARQUET_LOGICAL_STRING = 1,
    CARQUET_LOGICAL_MAP = 2,
    CARQUET_LOGICAL_LIST = 3,
    CARQUET_LOGICAL_ENUM = 4,
    CARQUET_LOGICAL_DECIMAL = 5,
    CARQUET_LOGICAL_DATE = 6,
    CARQUET_LOGICAL_TIME = 7,
    CARQUET_LOGICAL_TIMESTAMP = 8,
    CARQUET_LOGICAL_INTEGER = 9,
    CARQUET_LOGICAL_NULL = 10,
    CARQUET_LOGICAL_JSON = 11,
    CARQUET_LOGICAL_BSON = 12,
    CARQUET_LOGICAL_UUID = 13,
    CARQUET_LOGICAL_FLOAT16 = 14,
    CARQUET_LOGICAL_VARIANT = 15,
    CARQUET_LOGICAL_GEOMETRY = 16,
    CARQUET_LOGICAL_GEOGRAPHY = 17,
    /* INTERVAL has no modern LogicalType; it is ConvertedType-only and
       requires FIXED_LEN_BYTE_ARRAY with type_length == 12. */
    CARQUET_LOGICAL_INTERVAL = 18,
} carquet_logical_type_id_t;

/* Time unit for temporal types */
typedef enum carquet_time_unit {
    CARQUET_TIME_UNIT_MILLIS = 0,
    CARQUET_TIME_UNIT_MICROS = 1,
    CARQUET_TIME_UNIT_NANOS = 2,
} carquet_time_unit_t;

#define CARQUET_GEOSPATIAL_CRS_MAX 128

typedef enum carquet_geospatial_edge_algorithm {
    CARQUET_GEOSPATIAL_EDGE_SPHERICAL = 0,
    CARQUET_GEOSPATIAL_EDGE_VINCENTY = 1,
    CARQUET_GEOSPATIAL_EDGE_THOMAS = 2,
    CARQUET_GEOSPATIAL_EDGE_ANDOYER = 3,
    CARQUET_GEOSPATIAL_EDGE_KARNEY = 4,
} carquet_geospatial_edge_algorithm_t;

/* Logical type with parameters */
typedef struct carquet_logical_type {
    carquet_logical_type_id_t id;

    union {
        /* For DECIMAL */
        struct {
            int32_t precision;
            int32_t scale;
        } decimal;

        /* For INTEGER */
        struct {
            int8_t bit_width;   /* 8, 16, 32, or 64 */
            bool is_signed;
        } integer;

        /* For TIME */
        struct {
            carquet_time_unit_t unit;
            bool is_adjusted_to_utc;
        } time;

        /* For TIMESTAMP */
        struct {
            carquet_time_unit_t unit;
            bool is_adjusted_to_utc;
        } timestamp;

        /* For VARIANT */
        struct {
            int8_t specification_version;  /* 1 when unset/zero */
        } variant;

        /* For GEOMETRY */
        struct {
            char crs[CARQUET_GEOSPATIAL_CRS_MAX];  /* Optional, empty => OGC:CRS84 */
        } geometry;

        /* For GEOGRAPHY */
        struct {
            char crs[CARQUET_GEOSPATIAL_CRS_MAX];  /* Optional, empty => OGC:CRS84 */
            carquet_geospatial_edge_algorithm_t algorithm;
            bool has_algorithm;                    /* false => SPHERICAL */
        } geography;
    } params;
} carquet_logical_type_t;

/* ============================================================================
 * Field Repetition
 * ============================================================================
 */

typedef enum carquet_field_repetition {
    CARQUET_REPETITION_REQUIRED = 0,   /* Exactly one value */
    CARQUET_REPETITION_OPTIONAL = 1,   /* Zero or one value */
    CARQUET_REPETITION_REPEATED = 2,   /* Zero or more values */
} carquet_field_repetition_t;

/* ============================================================================
 * Encoding Types
 * ============================================================================
 */

typedef enum carquet_encoding {
    CARQUET_ENCODING_PLAIN = 0,
    CARQUET_ENCODING_PLAIN_DICTIONARY = 2,  /* Deprecated */
    CARQUET_ENCODING_RLE = 3,
    CARQUET_ENCODING_BIT_PACKED = 4,        /* Deprecated */
    CARQUET_ENCODING_DELTA_BINARY_PACKED = 5,
    CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY = 6,
    CARQUET_ENCODING_DELTA_BYTE_ARRAY = 7,
    CARQUET_ENCODING_RLE_DICTIONARY = 8,
    CARQUET_ENCODING_BYTE_STREAM_SPLIT = 9,
} carquet_encoding_t;

/* ============================================================================
 * Compression Codecs
 * ============================================================================
 */

typedef enum carquet_compression {
    CARQUET_COMPRESSION_UNCOMPRESSED = 0,
    CARQUET_COMPRESSION_SNAPPY = 1,
    CARQUET_COMPRESSION_GZIP = 2,
    CARQUET_COMPRESSION_LZO = 3,
    CARQUET_COMPRESSION_BROTLI = 4,
    CARQUET_COMPRESSION_LZ4 = 5,
    CARQUET_COMPRESSION_ZSTD = 6,
    CARQUET_COMPRESSION_LZ4_RAW = 7,
} carquet_compression_t;

/* ============================================================================
 * Page Types
 * ============================================================================
 */

typedef enum carquet_page_type {
    CARQUET_PAGE_DATA = 0,
    CARQUET_PAGE_INDEX = 1,
    CARQUET_PAGE_DICTIONARY = 2,
    CARQUET_PAGE_DATA_V2 = 3,
} carquet_page_type_t;

/* ============================================================================
 * Value Types for C API
 * ============================================================================
 */

/* Fixed-length byte array */
typedef struct carquet_fixed_byte_array {
    uint8_t* data;
    int32_t length;
} carquet_fixed_byte_array_t;

/* Variable-length byte array */
typedef struct carquet_byte_array {
    uint8_t* data;
    int32_t length;
} carquet_byte_array_t;

/* INT96 (deprecated, for legacy timestamp support) */
typedef struct carquet_int96 {
    uint32_t value[3];
} carquet_int96_t;

/* Decimal value (for high-precision decimals) */
typedef struct carquet_decimal128 {
    int64_t low;
    int64_t high;
} carquet_decimal128_t;

/* ============================================================================
 * Type Information Utilities
 * ============================================================================
 */

/**
 * Get the size in bytes of a physical type.
 * Returns -1 for variable-length types (BYTE_ARRAY).
 */
static inline int32_t carquet_physical_type_size(carquet_physical_type_t type) {
    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN: return 1;
        case CARQUET_PHYSICAL_INT32: return 4;
        case CARQUET_PHYSICAL_INT64: return 8;
        case CARQUET_PHYSICAL_INT96: return 12;
        case CARQUET_PHYSICAL_FLOAT: return 4;
        case CARQUET_PHYSICAL_DOUBLE: return 8;
        case CARQUET_PHYSICAL_BYTE_ARRAY: return -1;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return -1;
        default: return -1;
    }
}

/**
 * Get a human-readable name for a physical type.
 */
const char* carquet_physical_type_name(carquet_physical_type_t type);

/**
 * Get a human-readable name for a compression codec.
 */
const char* carquet_compression_name(carquet_compression_t codec);

/**
 * Get a human-readable name for an encoding.
 */
const char* carquet_encoding_name(carquet_encoding_t encoding);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_TYPES_H */
