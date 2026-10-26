/**
 * @file parquet_types.h
 * @brief Parquet Thrift structure definitions
 *
 * These structures match the Parquet Thrift specification.
 * They are parsed from the file footer metadata.
 */

#ifndef CARQUET_PARQUET_TYPES_H
#define CARQUET_PARQUET_TYPES_H

#include <carquet/types.h>
#include <carquet/error.h>
#include "core/arena.h"
#include "thrift_decode.h"
#include "thrift_encode.h"
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Forward Declarations
 * ============================================================================
 */

typedef struct parquet_schema_element parquet_schema_element_t;
typedef struct parquet_statistics parquet_statistics_t;
typedef struct parquet_page_encoding_stats parquet_page_encoding_stats_t;
typedef struct parquet_column_metadata parquet_column_metadata_t;
typedef struct parquet_column_chunk parquet_column_chunk_t;
typedef struct parquet_sorting_column parquet_sorting_column_t;
typedef struct parquet_row_group parquet_row_group_t;
typedef struct parquet_key_value parquet_key_value_t;
typedef struct parquet_file_metadata parquet_file_metadata_t;
typedef struct parquet_page_header parquet_page_header_t;
typedef struct parquet_data_page_header parquet_data_page_header_t;
typedef struct parquet_data_page_header_v2 parquet_data_page_header_v2_t;
typedef struct parquet_dictionary_page_header parquet_dictionary_page_header_t;

/* ============================================================================
 * Schema Element
 * ============================================================================
 */

struct parquet_schema_element {
    /* Field 1: type (optional for groups) */
    bool has_type;
    carquet_physical_type_t type;

    /* Field 2: type_length (for FIXED_LEN_BYTE_ARRAY) */
    int32_t type_length;

    /* Field 3: repetition_type */
    bool has_repetition;
    carquet_field_repetition_t repetition_type;

    /* Field 4: name */
    char* name;

    /* Field 5: num_children (for groups) */
    int32_t num_children;

    /* Field 6: converted_type (legacy logical type) */
    bool has_converted_type;
    carquet_converted_type_t converted_type;

    /* Field 7: scale (for DECIMAL) */
    int32_t scale;

    /* Field 8: precision (for DECIMAL) */
    int32_t precision;

    /* Field 9: field_id */
    bool has_field_id;
    int32_t field_id;

    /* Field 10: logicalType (modern logical type) */
    bool has_logical_type;
    carquet_logical_type_t logical_type;

    /* Carquet extension (NOT part of the Parquet SchemaElement wire format):
     * per-field key/value metadata mirroring Arrow's Field.custom_metadata,
     * used for variable labels/descriptions. Emitted only into the
     * "ARROW:schema" footer blob; never written to the Parquet schema itself.
     * Owned by whoever populates it (see carquet_schema_set_field_metadata /
     * store_schema_elements). NULL / 0 when the field carries no metadata. */
    int32_t num_field_metadata;
    parquet_key_value_t* field_metadata;

    /* Carquet extension: an Arrow type refinement recovered from the
     * "ARROW:schema" footer blob that the Parquet type system alone cannot
     * express (e.g. 64-bit-offset LargeUtf8 / LargeBinary / LargeList).
     * 0 (== CARQUET_ARROW_REFINE_NONE) when no refinement was recovered.
     * Values mirror carquet_arrow_type_refinement_t. Read-only: never written
     * to the Parquet schema. */
    int32_t arrow_type_refinement;
};

/* ============================================================================
 * Statistics
 * ============================================================================
 */

struct parquet_statistics {
    /* Field 1: max (deprecated, use max_value) */
    uint8_t* max_deprecated;
    int32_t max_deprecated_len;

    /* Field 2: min (deprecated, use min_value) */
    uint8_t* min_deprecated;
    int32_t min_deprecated_len;

    /* Field 3: null_count */
    bool has_null_count;
    int64_t null_count;

    /* Field 4: distinct_count */
    bool has_distinct_count;
    int64_t distinct_count;

    /* Field 5: max_value.
     * has_max_value records that the field was present, independently of its
     * length, so that a legitimately empty (zero-length) BYTE_ARRAY min/max is
     * not mistaken for absent. max_value may be NULL when the value is empty. */
    bool has_max_value;
    uint8_t* max_value;
    int32_t max_value_len;

    /* Field 6: min_value */
    bool has_min_value;
    uint8_t* min_value;
    int32_t min_value_len;

    /* Field 7: is_max_value_exact */
    bool has_is_max_value_exact;
    bool is_max_value_exact;

    /* Field 8: is_min_value_exact */
    bool has_is_min_value_exact;
    bool is_min_value_exact;
};

/* ============================================================================
 * Geospatial Statistics (for GEOMETRY / GEOGRAPHY logical types)
 * ============================================================================
 * Parquet GeospatialStatistics: a coordinate bounding box plus the set of
 * ISO-WKB geometry type codes present in the column chunk.
 */

#define CARQUET_GEO_MAX_TYPES 64

struct parquet_geospatial_statistics {
    bool valid;                 /* at least one finite coordinate accumulated */
    /* BoundingBox: x/y always present; z/m only if seen. */
    double xmin, xmax, ymin, ymax;
    bool has_z;
    double zmin, zmax;
    bool has_m;
    double mmin, mmax;
    /* Distinct ISO-WKB geometry type codes encountered (e.g. 1=Point XY,
     * 1001=Point XYZ). Empty list means "unknown". */
    int32_t types[CARQUET_GEO_MAX_TYPES];
    int32_t num_types;
};
typedef struct parquet_geospatial_statistics parquet_geospatial_statistics_t;

/* ============================================================================
 * Size Statistics (Parquet 2.9)
 * ============================================================================
 *
 * ColumnMetaData field 16. Histograms are heap/arena-allocated; a length of 0
 * means the corresponding optional list is absent.
 */

struct parquet_size_statistics {
    /* Field 1: total unencoded BYTE_ARRAY value bytes (length prefixes
     * excluded). Only meaningful for BYTE_ARRAY columns. */
    bool has_unencoded_byte_array_data_bytes;
    int64_t unencoded_byte_array_data_bytes;

    /* Field 2: repetition_level_histogram (length max_rep_level + 1). */
    int64_t* repetition_level_histogram;
    int32_t repetition_level_histogram_len;

    /* Field 3: definition_level_histogram (length max_def_level + 1). */
    int64_t* definition_level_histogram;
    int32_t definition_level_histogram_len;
};
typedef struct parquet_size_statistics parquet_size_statistics_t;

/* ============================================================================
 * Page Encoding Stats
 * ============================================================================
 */

struct parquet_page_encoding_stats {
    carquet_page_type_t page_type;
    carquet_encoding_t encoding;
    int32_t count;
};

/* ============================================================================
 * Column Metadata
 * ============================================================================
 */

struct parquet_column_metadata {
    /* Field 1: type */
    carquet_physical_type_t type;

    /* Field 2: encodings */
    carquet_encoding_t* encodings;
    int32_t num_encodings;

    /* Field 3: path_in_schema */
    char** path_in_schema;
    int32_t path_len;

    /* Field 4: codec */
    carquet_compression_t codec;

    /* Field 5: num_values */
    int64_t num_values;

    /* Field 6: total_uncompressed_size */
    int64_t total_uncompressed_size;

    /* Field 7: total_compressed_size */
    int64_t total_compressed_size;

    /* Field 8: key_value_metadata */
    parquet_key_value_t* key_value_metadata;
    int32_t num_key_value;

    /* Field 9: data_page_offset */
    int64_t data_page_offset;

    /* Field 10: index_page_offset */
    bool has_index_page_offset;
    int64_t index_page_offset;

    /* Field 11: dictionary_page_offset */
    bool has_dictionary_page_offset;
    int64_t dictionary_page_offset;

    /* Field 12: statistics */
    bool has_statistics;
    parquet_statistics_t statistics;

    /* Field 13: encoding_stats */
    parquet_page_encoding_stats_t* encoding_stats;
    int32_t num_encoding_stats;

    /* Field 14: bloom_filter_offset */
    bool has_bloom_filter_offset;
    int64_t bloom_filter_offset;

    /* Field 15: bloom_filter_length */
    bool has_bloom_filter_length;
    int32_t bloom_filter_length;

    /* Field 16: size_statistics (Parquet 2.9) */
    bool has_size_statistics;
    parquet_size_statistics_t size_statistics;

    /* Field 17: geospatial_statistics (GEOMETRY / GEOGRAPHY) */
    bool has_geospatial_statistics;
    parquet_geospatial_statistics_t geospatial_statistics;
};

/* ============================================================================
 * Column Chunk
 * ============================================================================
 */

struct parquet_column_chunk {
    /* Field 1: file_path (optional, for external files) */
    char* file_path;

    /* Field 2: file_offset */
    int64_t file_offset;

    /* Field 3: meta_data */
    bool has_metadata;
    parquet_column_metadata_t metadata;

    /* Field 4: offset_index_offset */
    bool has_offset_index_offset;
    int64_t offset_index_offset;

    /* Field 5: offset_index_length */
    bool has_offset_index_length;
    int32_t offset_index_length;

    /* Field 6: column_index_offset */
    bool has_column_index_offset;
    int64_t column_index_offset;

    /* Field 7: column_index_length */
    bool has_column_index_length;
    int32_t column_index_length;
};

/* ============================================================================
 * Row Group
 * ============================================================================
 */

/* Parquet Thrift: SortingColumn */
struct parquet_sorting_column {
    /* Field 1: column_idx (ordinal position of the column in this row group) */
    int32_t column_idx;
    /* Field 2: descending (true => sorted descending) */
    bool descending;
    /* Field 3: nulls_first (true => nulls before non-null values) */
    bool nulls_first;
};

struct parquet_row_group {
    /* Field 1: columns */
    parquet_column_chunk_t* columns;
    int32_t num_columns;

    /* Field 2: total_byte_size */
    int64_t total_byte_size;

    /* Field 3: num_rows */
    int64_t num_rows;

    /* Field 4: sorting_columns (optional) */
    parquet_sorting_column_t* sorting_columns;
    int32_t num_sorting_columns;

    /* Field 5: file_offset */
    bool has_file_offset;
    int64_t file_offset;

    /* Field 6: total_compressed_size */
    bool has_total_compressed_size;
    int64_t total_compressed_size;

    /* Field 7: ordinal */
    bool has_ordinal;
    int16_t ordinal;
};

/* ============================================================================
 * Key-Value Metadata
 * ============================================================================
 */

struct parquet_key_value {
    char* key;
    char* value;  /* Can be NULL */
};

/* ============================================================================
 * File Metadata
 * ============================================================================
 */

struct parquet_file_metadata {
    /* Field 1: version */
    int32_t version;

    /* Field 2: schema */
    parquet_schema_element_t* schema;
    int32_t num_schema_elements;

    /* Field 3: num_rows */
    int64_t num_rows;

    /* Field 4: row_groups */
    parquet_row_group_t* row_groups;
    int32_t num_row_groups;

    /* Field 5: key_value_metadata */
    parquet_key_value_t* key_value_metadata;
    int32_t num_key_value;

    /* Field 6: created_by */
    char* created_by;

    /* Field 7: column_orders (one ColumnOrder union per column).
     * column_order_types[i] holds the set union member's Thrift field id
     * (1 = TypeDefinedOrder, the only member the spec defines). NULL when the
     * footer omitted field 7; then only num_column_orders is meaningful. */
    int32_t num_column_orders;
    int16_t* column_order_types;

    /* Field 8: encryption_algorithm (we skip for now) */

    /* Field 9: footer_signing_key_metadata (we skip for now) */
};

/* ============================================================================
 * Page Headers
 * ============================================================================
 */

struct parquet_data_page_header {
    /* Field 1: num_values */
    int32_t num_values;

    /* Field 2: encoding */
    carquet_encoding_t encoding;

    /* Field 3: definition_level_encoding */
    carquet_encoding_t definition_level_encoding;

    /* Field 4: repetition_level_encoding */
    carquet_encoding_t repetition_level_encoding;

    /* Field 5: statistics */
    bool has_statistics;
    parquet_statistics_t statistics;
};

struct parquet_data_page_header_v2 {
    /* Field 1: num_values */
    int32_t num_values;

    /* Field 2: num_nulls */
    int32_t num_nulls;

    /* Field 3: num_rows */
    int32_t num_rows;

    /* Field 4: encoding */
    carquet_encoding_t encoding;

    /* Field 5: definition_levels_byte_length */
    int32_t definition_levels_byte_length;

    /* Field 6: repetition_levels_byte_length */
    int32_t repetition_levels_byte_length;

    /* Field 7: is_compressed */
    bool is_compressed;

    /* Field 8: statistics */
    bool has_statistics;
    parquet_statistics_t statistics;
};

struct parquet_dictionary_page_header {
    /* Field 1: num_values */
    int32_t num_values;

    /* Field 2: encoding */
    carquet_encoding_t encoding;

    /* Field 3: is_sorted */
    bool is_sorted;
};

struct parquet_page_header {
    /* Field 1: type */
    carquet_page_type_t type;

    /* Field 2: uncompressed_page_size */
    int32_t uncompressed_page_size;

    /* Field 3: compressed_page_size */
    int32_t compressed_page_size;

    /* Field 4: crc */
    bool has_crc;
    int32_t crc;

    /* One of these based on type */
    union {
        parquet_data_page_header_t data_page_header;
        parquet_data_page_header_v2_t data_page_header_v2;
        parquet_dictionary_page_header_t dictionary_page_header;
    };
};

/* ============================================================================
 * Parsing Functions
 * ============================================================================
 */

/**
 * Parse file metadata from Thrift data.
 *
 * @param data Thrift-encoded metadata
 * @param size Size of data
 * @param arena Arena for allocations
 * @param metadata Output metadata structure
 * @param error Error information
 * @return Status code
 */
carquet_status_t parquet_parse_file_metadata(
    const uint8_t* data,
    size_t size,
    carquet_arena_t* arena,
    parquet_file_metadata_t* metadata,
    carquet_error_t* error);

/**
 * Parse a page header from Thrift data.
 *
 * @param data Thrift-encoded page header
 * @param size Size of data
 * @param header Output page header
 * @param bytes_read Output: number of bytes consumed
 * @param error Error information
 * @return Status code
 */
carquet_status_t parquet_parse_page_header(
    const uint8_t* data,
    size_t size,
    parquet_page_header_t* header,
    size_t* bytes_read,
    carquet_error_t* error);

/**
 * Free file metadata (only frees non-arena allocations).
 */
void parquet_file_metadata_free(parquet_file_metadata_t* metadata);

/* ============================================================================
 * Writing Functions
 * ============================================================================
 */

/**
 * Write file metadata to a buffer.
 *
 * @param metadata Metadata to write
 * @param buffer Output buffer
 * @param error Error information
 * @return Status code
 */
carquet_status_t parquet_write_file_metadata(
    const parquet_file_metadata_t* metadata,
    carquet_buffer_t* buffer,
    carquet_error_t* error);

/**
 * Write a page header to a buffer.
 *
 * @param header Page header to write
 * @param buffer Output buffer
 * @param error Error information
 * @return Status code
 */
carquet_status_t parquet_write_page_header(
    const parquet_page_header_t* header,
    carquet_buffer_t* buffer,
    carquet_error_t* error);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_PARQUET_TYPES_H */
