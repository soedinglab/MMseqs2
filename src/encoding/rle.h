/**
 * @file rle.h
 * @brief RLE/Bit-packing hybrid encoding for Parquet
 *
 * This encoding combines run-length encoding for repeated values with
 * bit-packing for sequences of distinct values. It's primarily used for
 * definition levels, repetition levels, and dictionary indices.
 *
 * Format:
 * - Each run starts with a header varint
 * - If (header & 1) == 0: RLE run, count = header >> 1, followed by value
 * - If (header & 1) == 1: Bit-packed run, count = (header >> 1) * 8, followed by packed values
 */

#ifndef CARQUET_ENCODING_RLE_H
#define CARQUET_ENCODING_RLE_H

#include <carquet/error.h>
#include "core/buffer.h"
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * RLE Decoder
 * ============================================================================
 */

/**
 * RLE decoder state.
 */
typedef struct carquet_rle_decoder {
    const uint8_t* data;
    size_t size;
    size_t pos;

    int bit_width;           /* Bits per value */
    uint32_t value_mask;     /* Mask for extracting values */

    /* Current run state */
    bool in_rle_run;
    int64_t run_remaining;   /* Values remaining in current run */
    uint32_t rle_value;      /* Value for RLE runs */

    /* Bit-pack buffer */
    uint32_t bitpack_buffer[8];
    int bitpack_pos;         /* Position within buffer */
    int bitpack_count;       /* Values in buffer */

    carquet_status_t status;
} carquet_rle_decoder_t;

/**
 * Initialize an RLE decoder.
 *
 * @param dec Decoder to initialize
 * @param data Input data
 * @param size Size of input data
 * @param bit_width Bits per value (0-32)
 */
void carquet_rle_decoder_init(
    carquet_rle_decoder_t* dec,
    const uint8_t* data,
    size_t size,
    int bit_width);

/**
 * Check if decoder has more values.
 */
bool carquet_rle_decoder_has_next(const carquet_rle_decoder_t* dec);

/**
 * Get a single value from the decoder.
 *
 * @param dec Decoder
 * @return Value, or 0 if no more values or error
 */
uint32_t carquet_rle_decoder_get(carquet_rle_decoder_t* dec);

/**
 * Get multiple values from the decoder.
 *
 * @param dec Decoder
 * @param output Output buffer
 * @param count Maximum values to get
 * @return Number of values actually read
 */
int64_t carquet_rle_decoder_get_batch(
    carquet_rle_decoder_t* dec,
    uint32_t* output,
    int64_t count);

/**
 * Skip values in the decoder.
 *
 * @param dec Decoder
 * @param count Number of values to skip
 * @return Number of values actually skipped
 */
int64_t carquet_rle_decoder_skip(
    carquet_rle_decoder_t* dec,
    int64_t count);

/**
 * Get decoder error status.
 */
static inline carquet_status_t carquet_rle_decoder_status(
    const carquet_rle_decoder_t* dec) {
    return dec->status;
}

/* ============================================================================
 * RLE Encoder
 * ============================================================================
 */

/**
 * RLE encoder state.
 */
typedef struct carquet_rle_encoder {
    carquet_buffer_t* buffer;
    int bit_width;

    /* Run detection */
    uint32_t prev_value;
    int64_t repeat_count;    /* Count of repeated values */
    bool has_prev;

    /* Bit-pack buffer: the next (incomplete) group of 8 values. */
    uint32_t bitpack_buffer[8];
    int bitpack_count;
    int64_t bitpack_total;   /* Total values in current bit-pack sequence */

    /* Complete packed groups waiting for a shared header. One bit-packed
     * run may hold up to 63 groups (header = (groups << 1) | 1, a 7-bit
     * varint), so a run-free stream costs one header byte per 504 values
     * instead of one per 8. */
    uint8_t staged[63 * 32];
    int staged_groups;

    carquet_status_t status;
} carquet_rle_encoder_t;

/**
 * Initialize an RLE encoder.
 *
 * @param enc Encoder to initialize
 * @param buffer Output buffer
 * @param bit_width Bits per value (0-32)
 */
void carquet_rle_encoder_init(
    carquet_rle_encoder_t* enc,
    carquet_buffer_t* buffer,
    int bit_width);

/**
 * Add a value to the encoder.
 *
 * @param enc Encoder
 * @param value Value to add
 * @return Status code
 */
carquet_status_t carquet_rle_encoder_put(
    carquet_rle_encoder_t* enc,
    uint32_t value);

/**
 * Add multiple identical values.
 *
 * @param enc Encoder
 * @param value Value to add
 * @param count Number of times to add
 * @return Status code
 */
carquet_status_t carquet_rle_encoder_put_repeat(
    carquet_rle_encoder_t* enc,
    uint32_t value,
    int64_t count);

/**
 * Flush any buffered data.
 * Must be called after all values have been added.
 *
 * @param enc Encoder
 * @return Status code
 */
carquet_status_t carquet_rle_encoder_flush(carquet_rle_encoder_t* enc);

/* ============================================================================
 * Convenience Functions
 * ============================================================================
 */

/**
 * Decode all RLE values into a buffer.
 *
 * @param input Input RLE data
 * @param input_size Size of input data
 * @param bit_width Bits per value
 * @param output Output buffer
 * @param max_values Maximum values to decode
 * @return Number of values decoded, or -1 on error
 */
int64_t carquet_rle_decode_all(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    uint32_t* output,
    int64_t max_values);

/**
 * Decode RLE values directly to int16 (for levels).
 */
int64_t carquet_rle_decode_levels(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    int16_t* output,
    int64_t max_values);

/**
 * Encode values using RLE.
 *
 * @param input Input values
 * @param count Number of values
 * @param bit_width Bits per value
 * @param output Output buffer
 * @return Status code
 */
carquet_status_t carquet_rle_encode_all(
    const uint32_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output);

/**
 * Encode levels (int16) using RLE.
 */
carquet_status_t carquet_rle_encode_levels(
    const int16_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output);

/* ============================================================================
 * Level Decoding with Prefix Length
 * ============================================================================
 */

/**
 * Decode levels that have a 4-byte length prefix.
 * This is the format used in Parquet data pages.
 *
 * @param input Input data (starts with 4-byte length)
 * @param input_size Size of input data
 * @param bit_width Bits per value
 * @param output Output buffer
 * @param max_values Maximum values to decode
 * @param bytes_consumed Output: total bytes consumed including length prefix
 * @return Number of values decoded, or -1 on error
 */
int64_t carquet_rle_decode_levels_prefixed(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    int16_t* output,
    int64_t max_values,
    size_t* bytes_consumed);

/**
 * Decode RLE-encoded 1-bit def levels directly into a null bitmap.
 *
 * Optimized for max_def_level == 1 (flat nullable columns). Decodes
 * RLE/bitpacked values directly into bitmap bits, skipping the
 * intermediate int16_t[] buffer and subsequent build_null_bitmap pass.
 *
 * Bitmap convention (matches build_null_bitmap):
 * bit set (1) = value is present (def_level == 1),
 * bit clear (0) = value is null (def_level == 0).
 *
 * @param input Input RLE data (no length prefix)
 * @param input_size Size of input data
 * @param bitmap Output bitmap (must be pre-allocated, (max_values+7)/8 bytes)
 * @param max_values Maximum values to decode
 * @param non_null_count Output: number of non-null values decoded
 * @return Number of values decoded, or -1 on error
 */
int64_t carquet_rle_decode_to_bitmap(
    const uint8_t* input,
    size_t input_size,
    uint8_t* bitmap,
    int64_t max_values,
    int64_t* non_null_count);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_ENCODING_RLE_H */
