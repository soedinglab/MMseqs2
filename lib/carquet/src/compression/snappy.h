/**
 * @file snappy.h
 * @brief Internal Snappy codec interface with detailed failure reporting.
 *
 * `carquet_snappy_decompress` collapses every rejection into
 * CARQUET_ERROR_INVALID_COMPRESSED_DATA, which makes a failure on a real file
 * impossible to triage without a debugger. `carquet_snappy_decompress_diag`
 * additionally reports *which* check fired and where, so a single failing run
 * on a multi-gigabyte file is enough to identify the cause.
 */

#ifndef CARQUET_COMPRESSION_SNAPPY_H
#define CARQUET_COMPRESSION_SNAPPY_H

#include <carquet/error.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Distinct reasons a Snappy block can be rejected. */
typedef enum carquet_snappy_reason {
    CARQUET_SNAPPY_OK = 0,
    CARQUET_SNAPPY_ERR_LENGTH_VARINT,     /**< Leading uncompressed-length varint malformed */
    CARQUET_SNAPPY_ERR_OUTPUT_TOO_SMALL,  /**< Declared length exceeds the caller's buffer */
    CARQUET_SNAPPY_ERR_LITERAL_LEN_TRUNC, /**< Long-literal length bytes run past input end */
    CARQUET_SNAPPY_ERR_LITERAL_INPUT,     /**< Literal payload runs past input end */
    CARQUET_SNAPPY_ERR_LITERAL_OUTPUT,    /**< Literal payload runs past output end */
    CARQUET_SNAPPY_ERR_COPY1_TRUNC,       /**< COPY_1 trailer byte missing */
    CARQUET_SNAPPY_ERR_COPY2_TRUNC,       /**< COPY_2 trailer bytes missing */
    CARQUET_SNAPPY_ERR_COPY4_TRUNC,       /**< COPY_4 trailer bytes missing */
    CARQUET_SNAPPY_ERR_COPY_OFFSET_ZERO,  /**< Copy offset of 0 is not representable */
    CARQUET_SNAPPY_ERR_COPY_OFFSET_RANGE, /**< Copy offset points before the output start */
    CARQUET_SNAPPY_ERR_COPY_OUTPUT,       /**< Copy runs past output end */
    CARQUET_SNAPPY_ERR_SHORT_OUTPUT,      /**< Input exhausted before the declared length */
    CARQUET_SNAPPY_ERR_TRAILING_INPUT     /**< Output complete but input bytes remain */
} carquet_snappy_reason_t;

/** @brief Where and why a Snappy block was rejected. */
typedef struct carquet_snappy_diag {
    carquet_snappy_reason_t reason;
    size_t input_pos;         /**< Offset in src of the tag being decoded */
    size_t output_pos;        /**< Bytes already produced when the check fired */
    uint64_t declared_length; /**< Uncompressed length read from the block header */
    int32_t tag;              /**< Failing tag byte, or -1 when not applicable */
    uint64_t operand;         /**< Literal length / copy offset, depending on reason */
} carquet_snappy_diag_t;

/** @brief Human-readable one-line description of a rejection reason. */
const char* carquet_snappy_reason_string(carquet_snappy_reason_t reason);

/**
 * @brief Format a diagnostic into `buf` as a single line.
 * @return `buf`, always NUL-terminated.
 */
const char* carquet_snappy_diag_format(const carquet_snappy_diag_t* diag,
                                       char* buf, size_t buf_size);

/** @brief Decompress a raw Snappy block, reporting the exact rejection cause. */
carquet_status_t carquet_snappy_decompress_diag(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size,
    carquet_snappy_diag_t* diag);

/** @brief Decompress a raw Snappy block. */
carquet_status_t carquet_snappy_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

size_t carquet_snappy_compress_bound(size_t src_size);

carquet_status_t carquet_snappy_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_COMPRESSION_SNAPPY_H */
