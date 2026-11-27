/**
 * @file buffer.h
 * @brief Growable byte buffer
 *
 * A simple growable buffer for building byte sequences.
 * Used for encoding and building output pages.
 */

#ifndef CARQUET_CORE_BUFFER_H
#define CARQUET_CORE_BUFFER_H

#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Constants
 * ============================================================================
 */

#define CARQUET_BUFFER_DEFAULT_CAPACITY 4096

/* ============================================================================
 * Types
 * ============================================================================
 */

/**
 * Growable byte buffer.
 */
typedef struct carquet_buffer {
    uint8_t* data;      /* Buffer data */
    size_t size;        /* Current size (bytes written) */
    size_t capacity;    /* Allocated capacity */
    bool owns_data;     /* Whether buffer owns the data (should free) */
} carquet_buffer_t;

/* ============================================================================
 * Buffer Operations
 * ============================================================================
 */

/**
 * Initialize an empty buffer.
 * @pre buf != NULL (asserts on violation)
 */
void carquet_buffer_init(carquet_buffer_t* buf);

/**
 * Initialize a buffer with a specific capacity.
 */
carquet_status_t carquet_buffer_init_capacity(carquet_buffer_t* buf, size_t capacity);

/**
 * Initialize a buffer wrapping existing data (non-owning).
 */
void carquet_buffer_init_wrap(carquet_buffer_t* buf, uint8_t* data, size_t size);

/**
 * Initialize a buffer with a copy of existing data.
 */
carquet_status_t carquet_buffer_init_copy(carquet_buffer_t* buf,
                                           const uint8_t* data, size_t size);

/**
 * Destroy a buffer and free memory if owned.
 */
void carquet_buffer_destroy(carquet_buffer_t* buf);

/**
 * Clear buffer contents without freeing memory.
 */
void carquet_buffer_clear(carquet_buffer_t* buf);

/**
 * Ensure buffer has at least the specified capacity.
 */
carquet_status_t carquet_buffer_reserve(carquet_buffer_t* buf, size_t capacity);

/**
 * Resize buffer to exact size, truncating or zero-filling.
 */
carquet_status_t carquet_buffer_resize(carquet_buffer_t* buf, size_t size);

/**
 * Shrink buffer capacity to match current size.
 */
carquet_status_t carquet_buffer_shrink_to_fit(carquet_buffer_t* buf);

/* ============================================================================
 * Write Operations
 * ============================================================================
 */

/**
 * Append bytes to the buffer.
 */
carquet_status_t carquet_buffer_append(carquet_buffer_t* buf,
                                        const void* data, size_t size);

/**
 * Append a single byte.
 */
carquet_status_t carquet_buffer_append_byte(carquet_buffer_t* buf, uint8_t byte);

/**
 * Append bytes, repeating a value.
 */
carquet_status_t carquet_buffer_append_fill(carquet_buffer_t* buf,
                                             uint8_t value, size_t count);

/**
 * Append a 16-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_append_u16_le(carquet_buffer_t* buf, uint16_t value);

/**
 * Append a 32-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_append_u32_le(carquet_buffer_t* buf, uint32_t value);

/**
 * Append a 64-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_append_u64_le(carquet_buffer_t* buf, uint64_t value);

/**
 * Append a 32-bit float (little-endian).
 */
carquet_status_t carquet_buffer_append_f32_le(carquet_buffer_t* buf, float value);

/**
 * Append a 64-bit double (little-endian).
 */
carquet_status_t carquet_buffer_append_f64_le(carquet_buffer_t* buf, double value);

/**
 * Reserve space and return pointer to write directly.
 * The buffer size is increased by `size`.
 */
uint8_t* carquet_buffer_advance(carquet_buffer_t* buf, size_t size);

/* ============================================================================
 * Read Operations (for cursor-based reading)
 * ============================================================================
 */

/**
 * Buffer reader cursor.
 */
typedef struct carquet_buffer_reader {
    const uint8_t* data;
    size_t size;
    size_t pos;
} carquet_buffer_reader_t;

/**
 * Initialize a reader from a buffer.
 */
void carquet_buffer_reader_init(carquet_buffer_reader_t* reader,
                                 const carquet_buffer_t* buf);

/**
 * Initialize a reader from raw data.
 */
void carquet_buffer_reader_init_data(carquet_buffer_reader_t* reader,
                                      const uint8_t* data, size_t size);

/**
 * Get remaining bytes in reader.
 */
static inline size_t carquet_buffer_reader_remaining(const carquet_buffer_reader_t* reader) {
    return reader->size - reader->pos;
}

/**
 * Check if reader has at least n bytes remaining.
 */
static inline bool carquet_buffer_reader_has(const carquet_buffer_reader_t* reader, size_t n) {
    return reader->pos + n <= reader->size;
}

/**
 * Get pointer to current position without advancing.
 */
static inline const uint8_t* carquet_buffer_reader_peek(const carquet_buffer_reader_t* reader) {
    return reader->data + reader->pos;
}

/**
 * Read bytes into a buffer.
 */
carquet_status_t carquet_buffer_reader_read(carquet_buffer_reader_t* reader,
                                             void* dest, size_t size);

/**
 * Skip bytes.
 */
carquet_status_t carquet_buffer_reader_skip(carquet_buffer_reader_t* reader, size_t size);

/**
 * Read a single byte.
 */
carquet_status_t carquet_buffer_reader_read_byte(carquet_buffer_reader_t* reader,
                                                  uint8_t* value);

/**
 * Read a 16-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_reader_read_u16_le(carquet_buffer_reader_t* reader,
                                                    uint16_t* value);

/**
 * Read a 32-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_reader_read_u32_le(carquet_buffer_reader_t* reader,
                                                    uint32_t* value);

/**
 * Read a 64-bit integer (little-endian).
 */
carquet_status_t carquet_buffer_reader_read_u64_le(carquet_buffer_reader_t* reader,
                                                    uint64_t* value);

/**
 * Read a 32-bit float (little-endian).
 */
carquet_status_t carquet_buffer_reader_read_f32_le(carquet_buffer_reader_t* reader,
                                                    float* value);

/**
 * Read a 64-bit double (little-endian).
 */
carquet_status_t carquet_buffer_reader_read_f64_le(carquet_buffer_reader_t* reader,
                                                    double* value);

/* ============================================================================
 * Accessors
 * ============================================================================
 */

/**
 * Get buffer data pointer.
 */
static inline uint8_t* carquet_buffer_data(carquet_buffer_t* buf) {
    return buf->data;
}

/**
 * Get buffer data pointer (const).
 */
static inline const uint8_t* carquet_buffer_data_const(const carquet_buffer_t* buf) {
    return buf->data;
}

/**
 * Get buffer size.
 */
static inline size_t carquet_buffer_size(const carquet_buffer_t* buf) {
    return buf->size;
}

/**
 * Get buffer capacity.
 */
static inline size_t carquet_buffer_capacity(const carquet_buffer_t* buf) {
    return buf->capacity;
}

/**
 * Check if buffer is empty.
 */
static inline bool carquet_buffer_empty(const carquet_buffer_t* buf) {
    return buf->size == 0;
}

/* ============================================================================
 * Utility Operations
 * ============================================================================
 */

/**
 * Detach buffer data (caller takes ownership).
 * Buffer is reset to empty state.
 */
uint8_t* carquet_buffer_detach(carquet_buffer_t* buf, size_t* size_out);

/**
 * Swap contents of two buffers.
 */
void carquet_buffer_swap(carquet_buffer_t* a, carquet_buffer_t* b);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_CORE_BUFFER_H */
