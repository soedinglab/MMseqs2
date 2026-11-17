/**
 * @file arrow_leaf.h
 * @brief How a Parquet leaf column maps to an Arrow C Data value buffer.
 *
 * Shared by the batch exporter (arrow_c_export.c) and the nested assembler
 * (arrow_c_read.c), so both give a column the same Arrow type and the same
 * value bytes. The mapping matches pyarrow.parquet's reading of the column:
 * DECIMAL becomes decimal128 (precision <= 38) or decimal256 (<= 76) and the
 * deprecated INT96 becomes timestamp[ns]; a type Arrow has no lossless
 * equivalent for keeps its physical storage format (fixed-size binary,
 * binary, int).
 */
#ifndef CARQUET_READER_ARROW_LEAF_H
#define CARQUET_READER_ARROW_LEAF_H

#include <carquet/carquet.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
    CARQUET_ARROW_LEAF_COPY,      /* stored values, byte for byte */
    CARQUET_ARROW_LEAF_BOOL,      /* one byte per value -> one bit per value */
    CARQUET_ARROW_LEAF_BINARY,    /* BYTE_ARRAY -> int32 offsets + data */
    CARQUET_ARROW_LEAF_CONVERT,   /* per value through carquet_arrow_leaf_store() */
    CARQUET_ARROW_LEAF_NULL       /* Null logical type: Arrow null, no buffers */
} carquet_arrow_leaf_kind_t;

typedef struct {
    carquet_arrow_leaf_kind_t kind;
    carquet_physical_type_t physical;
    size_t src_width;   /* bytes per stored value (sizeof(carquet_byte_array_t) for BYTE_ARRAY) */
    size_t dst_width;   /* bytes per Arrow value; 0 for BOOL and BINARY */
    char format[32];    /* Arrow format string */
    /* CONVERT only: */
    uint8_t convert;    /* see arrow_leaf.c */
} carquet_arrow_leaf_t;

/**
 * Plan the export of a leaf of physical type @p pt (FLBA width @p type_length,
 * logical type @p lt, may be NULL).
 * @return CARQUET_OK, CARQUET_ERROR_INVALID_ARGUMENT for an FLBA without a
 *         width, CARQUET_ERROR_NOT_IMPLEMENTED for an unknown physical type.
 */
carquet_status_t carquet_arrow_leaf_plan(carquet_physical_type_t pt, int32_t type_length,
                                         const carquet_logical_type_t* lt,
                                         carquet_arrow_leaf_t* out);

/**
 * Store one value of a COPY or CONVERT leaf: @p src points at a stored value
 * (a carquet_byte_array_t for BYTE_ARRAY), @p dst receives dst_width bytes.
 * @return false when the value does not fit the Arrow type (a BYTE_ARRAY
 *         DECIMAL wider than its precision allows); @p dst is then zeroed.
 */
bool carquet_arrow_leaf_store(const carquet_arrow_leaf_t* plan, uint8_t* dst, const void* src);

#endif /* CARQUET_READER_ARROW_LEAF_H */
