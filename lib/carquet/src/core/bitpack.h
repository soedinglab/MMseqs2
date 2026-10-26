/**
 * @file bitpack.h
 * @brief Bit packing and unpacking utilities
 *
 * These functions handle packing and unpacking values at arbitrary bit widths,
 * which is essential for RLE/bit-packing hybrid encoding and delta encoding.
 */

#ifndef CARQUET_CORE_BITPACK_H
#define CARQUET_CORE_BITPACK_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Bit Manipulation Utilities
 * ============================================================================
 */

/**
 * Count leading zeros in a 32-bit integer.
 */
static inline int carquet_clz32(uint32_t v) {
    if (v == 0) return 32;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clz(v);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanReverse(&index, v);
    return 31 - (int)index;
#else
    int n = 0;
    if (v <= 0x0000FFFF) { n += 16; v <<= 16; }
    if (v <= 0x00FFFFFF) { n += 8;  v <<= 8; }
    if (v <= 0x0FFFFFFF) { n += 4;  v <<= 4; }
    if (v <= 0x3FFFFFFF) { n += 2;  v <<= 2; }
    if (v <= 0x7FFFFFFF) { n += 1; }
    return n;
#endif
}

/**
 * Count leading zeros in a 64-bit integer.
 */
static inline int carquet_clz64(uint64_t v) {
    if (v == 0) return 64;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(v);
#elif defined(_MSC_VER) && defined(_M_X64)
    unsigned long index;
    _BitScanReverse64(&index, v);
    return 63 - (int)index;
#else
    int n = 0;
    if (v <= 0x00000000FFFFFFFFULL) { n += 32; v <<= 32; }
    if (v <= 0x0000FFFFFFFFFFFFULL) { n += 16; v <<= 16; }
    if (v <= 0x00FFFFFFFFFFFFFFULL) { n += 8;  v <<= 8; }
    if (v <= 0x0FFFFFFFFFFFFFFFULL) { n += 4;  v <<= 4; }
    if (v <= 0x3FFFFFFFFFFFFFFFULL) { n += 2;  v <<= 2; }
    if (v <= 0x7FFFFFFFFFFFFFFFULL) { n += 1; }
    return n;
#endif
}

/**
 * Count trailing zeros in a 32-bit integer.
 */
static inline int carquet_ctz32(uint32_t v) {
    if (v == 0) return 32;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(v);
#elif defined(_MSC_VER)
    unsigned long index;
    _BitScanForward(&index, v);
    return (int)index;
#else
    int n = 31;
    if (v & 0x0000FFFF) { n -= 16; } else { v >>= 16; }
    if (v & 0x000000FF) { n -= 8; }  else { v >>= 8; }
    if (v & 0x0000000F) { n -= 4; }  else { v >>= 4; }
    if (v & 0x00000003) { n -= 2; }  else { v >>= 2; }
    if (v & 0x00000001) { n -= 1; }
    return n;
#endif
}

/**
 * Count population (number of set bits) in a 32-bit integer.
 */
static inline int carquet_popcount32(uint32_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(v);
#elif defined(_MSC_VER)
    return (int)__popcnt(v);
#else
    v = v - ((v >> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
    v = (v + (v >> 4)) & 0x0F0F0F0F;
    return (int)((v * 0x01010101) >> 24);
#endif
}

/**
 * Count population (number of set bits) in a 64-bit integer.
 */
static inline int carquet_popcount64(uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(v);
#elif defined(_MSC_VER) && defined(_M_X64)
    return (int)__popcnt64(v);
#else
    v = v - ((v >> 1) & 0x5555555555555555ULL);
    v = (v & 0x3333333333333333ULL) + ((v >> 2) & 0x3333333333333333ULL);
    v = (v + (v >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int)((v * 0x0101010101010101ULL) >> 56);
#endif
}

/**
 * Calculate bit width needed to represent a value.
 */
static inline int carquet_bit_width32(uint32_t v) {
    return v == 0 ? 0 : 32 - carquet_clz32(v);
}

/**
 * Calculate bit width needed to represent a value.
 */
static inline int carquet_bit_width64(uint64_t v) {
    return v == 0 ? 0 : 64 - carquet_clz64(v);
}

/* ============================================================================
 * Bit Packing (Scalar)
 * ============================================================================
 */

/**
 * Pack 8 values at the given bit width.
 *
 * This is the fundamental operation for bit-packing encoding.
 * Values are packed in LSB order within each byte.
 *
 * @param values Input values (8 values)
 * @param bit_width Bits per value (1-32)
 * @param output Output buffer (must have space for bit_width bytes)
 */
void carquet_bitpack8_32(const uint32_t* values, int bit_width, uint8_t* output);

/**
 * Unpack 8 values at the given bit width.
 *
 * @param input Input packed data (bit_width bytes)
 * @param bit_width Bits per value (1-32)
 * @param values Output values (8 values)
 */
void carquet_bitunpack8_32(const uint8_t* input, int bit_width, uint32_t* values);

/**
 * Pack N values at the given bit width.
 *
 * @param values Input values
 * @param count Number of values (should be multiple of 8 for efficiency)
 * @param bit_width Bits per value (1-32)
 * @param output Output buffer
 * @return Number of bytes written
 */
size_t carquet_bitpack_32(const uint32_t* values, size_t count,
                          int bit_width, uint8_t* output);

/**
 * Unpack N values at the given bit width.
 *
 * @param input Input packed data
 * @param count Number of values to unpack
 * @param bit_width Bits per value (1-32)
 * @param values Output values
 * @return Number of bytes consumed
 */
size_t carquet_bitunpack_32(const uint8_t* input, size_t count,
                            int bit_width, uint32_t* values);

/**
 * Calculate number of bytes needed to pack N values at given bit width.
 */
static inline size_t carquet_packed_size(size_t count, int bit_width) {
    return (count * (size_t)bit_width + 7) / 8;
}

/* ============================================================================
 * Specialized Unpack Functions (Performance Critical)
 * ============================================================================
 */

/**
 * Unpack values at specific bit widths (optimized versions).
 * These are called by the general unpack function but can be
 * called directly for known bit widths.
 */
void carquet_bitunpack8_1bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_2bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_3bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_4bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_5bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_6bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_7bit(const uint8_t* input, uint32_t* values);
void carquet_bitunpack8_8bit(const uint8_t* input, uint32_t* values);

/* ============================================================================
 * Function Pointer Type for SIMD Dispatch
 * ============================================================================
 */

typedef void (*carquet_bitunpack8_fn)(const uint8_t* input, uint32_t* values);
typedef void (*carquet_bitpack8_fn)(const uint32_t* values, uint8_t* output);

/**
 * Get the unpack function for a specific bit width.
 * Returns NULL for invalid bit widths.
 */
carquet_bitunpack8_fn carquet_get_bitunpack8_fn(int bit_width);

/**
 * Get the pack function for a specific bit width.
 * Returns NULL for invalid bit widths.
 */
carquet_bitpack8_fn carquet_get_bitpack8_fn(int bit_width);

/* ============================================================================
 * Bit Stream Reader/Writer
 * ============================================================================
 */

/**
 * Bit stream reader for arbitrary bit-level access.
 */
typedef struct carquet_bit_reader {
    const uint8_t* data;
    size_t size;
    size_t byte_pos;
    int bit_pos;  /* 0-7, bits remaining in current byte */
    uint64_t buffer;  /* Bit buffer for efficient reading */
    int buffer_bits;  /* Bits available in buffer */
} carquet_bit_reader_t;

/**
 * Initialize a bit reader.
 */
void carquet_bit_reader_init(carquet_bit_reader_t* reader,
                              const uint8_t* data, size_t size);

/**
 * Read a single bit.
 */
int carquet_bit_reader_read_bit(carquet_bit_reader_t* reader);

/**
 * Read up to 32 bits.
 */
uint32_t carquet_bit_reader_read_bits(carquet_bit_reader_t* reader, int num_bits);

/**
 * Read up to 64 bits.
 */
uint64_t carquet_bit_reader_read_bits64(carquet_bit_reader_t* reader, int num_bits);

/**
 * Check if reader has more data.
 */
bool carquet_bit_reader_has_more(const carquet_bit_reader_t* reader);

/**
 * Get remaining bits.
 */
size_t carquet_bit_reader_remaining_bits(const carquet_bit_reader_t* reader);

/**
 * Bit stream writer for arbitrary bit-level access.
 */
typedef struct carquet_bit_writer {
    uint8_t* data;
    size_t capacity;
    size_t byte_pos;
    int bit_pos;  /* 0-7, bits written in current byte */
    uint64_t buffer;  /* Bit buffer for efficient writing */
    int buffer_bits;  /* Bits in buffer */
} carquet_bit_writer_t;

/**
 * Initialize a bit writer.
 */
void carquet_bit_writer_init(carquet_bit_writer_t* writer,
                              uint8_t* data, size_t capacity);

/**
 * Write a single bit.
 */
void carquet_bit_writer_write_bit(carquet_bit_writer_t* writer, int bit);

/**
 * Write up to 32 bits.
 */
void carquet_bit_writer_write_bits(carquet_bit_writer_t* writer,
                                    uint32_t value, int num_bits);

/**
 * Write up to 64 bits.
 */
void carquet_bit_writer_write_bits64(carquet_bit_writer_t* writer,
                                      uint64_t value, int num_bits);

/**
 * Flush any remaining bits to output.
 */
void carquet_bit_writer_flush(carquet_bit_writer_t* writer);

/**
 * Get number of bytes written (after flush).
 */
size_t carquet_bit_writer_bytes_written(const carquet_bit_writer_t* writer);

/**
 * Decode the deprecated BIT_PACKED encoding (Parquet Encoding=4) used for
 * definition/repetition levels in legacy Data Page V1. Values are packed
 * MSB-first with no run headers and no length prefix; the byte length is
 * implied by ceil(count * bit_width / 8).
 *
 * @param data        Packed input.
 * @param data_size   Bytes available in @p data.
 * @param bit_width   Bits per value (0..16); 0 emits all-zero levels.
 * @param count       Number of level values to decode.
 * @param out         Output buffer for @p count int16 levels.
 * @param consumed    Set to the number of input bytes consumed.
 * @return 0 on success, -1 on bad arguments / truncated input.
 */
int carquet_decode_bitpacked_levels(const uint8_t* data, size_t data_size,
                                    int bit_width, int32_t count,
                                    int16_t* out, size_t* consumed);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_CORE_BITPACK_H */
