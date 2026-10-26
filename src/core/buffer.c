/**
 * @file buffer.c
 * @brief Growable byte buffer implementation
 */

#include "allocator.h"
#include "buffer.h"
#include "endian.h"
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static size_t next_power_of_two(size_t n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
#if SIZE_MAX > 0xFFFFFFFF
    n |= n >> 32;
#endif
    return n + 1;
}

static int add_overflows_size(size_t a, size_t b, size_t* out) {
    if (a > SIZE_MAX - b) {
        return 1;
    }
    *out = a + b;
    return 0;
}

static carquet_status_t ensure_capacity(carquet_buffer_t* buf, size_t needed) {
    if (needed <= buf->capacity) {
        return CARQUET_OK;
    }

    /* Don't grow non-owning buffers */
    if (!buf->owns_data && buf->data) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    size_t new_capacity = next_power_of_two(needed);
    if (new_capacity < needed) {
        new_capacity = needed;
    }
    if (new_capacity < CARQUET_BUFFER_DEFAULT_CAPACITY) {
        new_capacity = CARQUET_BUFFER_DEFAULT_CAPACITY;
    }

    uint8_t* new_data = (uint8_t*)carquet_mem_realloc(buf->data, new_capacity);
    if (!new_data) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    buf->data = new_data;
    buf->capacity = new_capacity;
    buf->owns_data = true;

    return CARQUET_OK;
}

/* ============================================================================
 * Buffer Operations
 * ============================================================================
 */

void carquet_buffer_init(carquet_buffer_t* buf) {
    assert(buf != NULL);

    buf->data = NULL;
    buf->size = 0;
    buf->capacity = 0;
    buf->owns_data = true;
}

carquet_status_t carquet_buffer_init_capacity(carquet_buffer_t* buf, size_t capacity) {
    assert(buf != NULL);
    carquet_buffer_init(buf);

    if (capacity > 0) {
        carquet_status_t status = carquet_buffer_reserve(buf, capacity);
        if (CARQUET_FAILED(status)) {
            return status;
        }
    }

    return CARQUET_OK;
}

void carquet_buffer_init_wrap(carquet_buffer_t* buf, uint8_t* data, size_t size) {
    assert(buf != NULL);

    buf->data = data;
    buf->size = size;
    buf->capacity = size;
    buf->owns_data = false;
}

carquet_status_t carquet_buffer_init_copy(carquet_buffer_t* buf,
                                           const uint8_t* data, size_t size) {
    carquet_status_t status = carquet_buffer_init_capacity(buf, size);
    if (CARQUET_FAILED(status)) {
        return status;
    }

    if (data && size > 0) {
        memcpy(buf->data, data, size);
        buf->size = size;
    }

    return CARQUET_OK;
}

void carquet_buffer_destroy(carquet_buffer_t* buf) {
    assert(buf != NULL);

    if (buf->owns_data && buf->data) {
        carquet_mem_free(buf->data);
    }

    buf->data = NULL;
    buf->size = 0;
    buf->capacity = 0;
    buf->owns_data = true;
}

void carquet_buffer_clear(carquet_buffer_t* buf) {
    assert(buf != NULL);
    buf->size = 0;
}

carquet_status_t carquet_buffer_reserve(carquet_buffer_t* buf, size_t capacity) {
    assert(buf != NULL);
    return ensure_capacity(buf, capacity);
}

carquet_status_t carquet_buffer_resize(carquet_buffer_t* buf, size_t size) {
    assert(buf != NULL);

    carquet_status_t status = ensure_capacity(buf, size);
    if (CARQUET_FAILED(status)) {
        return status;
    }

    /* Zero-fill if growing */
    if (size > buf->size) {
        memset(buf->data + buf->size, 0, size - buf->size);
    }

    buf->size = size;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_shrink_to_fit(carquet_buffer_t* buf) {
    assert(buf != NULL);
    assert(buf->owns_data);

    if (buf->size == 0) {
        carquet_mem_free(buf->data);
        buf->data = NULL;
        buf->capacity = 0;
        return CARQUET_OK;
    }

    if (buf->size < buf->capacity) {
        uint8_t* new_data = (uint8_t*)carquet_mem_realloc(buf->data, buf->size);
        if (new_data) {
            buf->data = new_data;
            buf->capacity = buf->size;
        }
        /* If realloc fails, keep the larger buffer */
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Write Operations
 * ============================================================================
 */

carquet_status_t carquet_buffer_append(carquet_buffer_t* buf,
                                        const void* data, size_t size) {
    assert(buf != NULL);
    if (size == 0) {
        return CARQUET_OK;
    }
    assert(data != NULL);

    size_t needed;
    if (add_overflows_size(buf->size, size, &needed)) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    carquet_status_t status = ensure_capacity(buf, needed);
    if (CARQUET_FAILED(status)) {
        return status;
    }

    memcpy(buf->data + buf->size, data, size);
    buf->size += size;

    return CARQUET_OK;
}

carquet_status_t carquet_buffer_append_byte(carquet_buffer_t* buf, uint8_t byte) {
    return carquet_buffer_append(buf, &byte, 1);
}

carquet_status_t carquet_buffer_append_fill(carquet_buffer_t* buf,
                                             uint8_t value, size_t count) {
    assert(buf != NULL);
    if (count == 0) {
        return CARQUET_OK;
    }

    size_t needed;
    if (add_overflows_size(buf->size, count, &needed)) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    carquet_status_t status = ensure_capacity(buf, needed);
    if (CARQUET_FAILED(status)) {
        return status;
    }

    memset(buf->data + buf->size, value, count);
    buf->size += count;

    return CARQUET_OK;
}

carquet_status_t carquet_buffer_append_u16_le(carquet_buffer_t* buf, uint16_t value) {
    uint8_t bytes[2];
    carquet_write_u16_le(bytes, value);
    return carquet_buffer_append(buf, bytes, 2);
}

carquet_status_t carquet_buffer_append_u32_le(carquet_buffer_t* buf, uint32_t value) {
    uint8_t bytes[4];
    carquet_write_u32_le(bytes, value);
    return carquet_buffer_append(buf, bytes, 4);
}

carquet_status_t carquet_buffer_append_u64_le(carquet_buffer_t* buf, uint64_t value) {
    uint8_t bytes[8];
    carquet_write_u64_le(bytes, value);
    return carquet_buffer_append(buf, bytes, 8);
}

carquet_status_t carquet_buffer_append_f32_le(carquet_buffer_t* buf, float value) {
    uint8_t bytes[4];
    carquet_write_f32_le(bytes, value);
    return carquet_buffer_append(buf, bytes, 4);
}

carquet_status_t carquet_buffer_append_f64_le(carquet_buffer_t* buf, double value) {
    uint8_t bytes[8];
    carquet_write_f64_le(bytes, value);
    return carquet_buffer_append(buf, bytes, 8);
}

uint8_t* carquet_buffer_advance(carquet_buffer_t* buf, size_t size) {
    assert(buf != NULL);
    if (size == 0) {
        return NULL;
    }

    size_t needed;
    if (add_overflows_size(buf->size, size, &needed)) {
        return NULL;
    }

    carquet_status_t status = ensure_capacity(buf, needed);
    if (CARQUET_FAILED(status)) {
        return NULL;
    }

    uint8_t* ptr = buf->data + buf->size;
    buf->size += size;
    return ptr;
}

/* ============================================================================
 * Reader Operations
 * ============================================================================
 */

void carquet_buffer_reader_init(carquet_buffer_reader_t* reader,
                                 const carquet_buffer_t* buf) {
    assert(reader != NULL);

    reader->data = buf ? buf->data : NULL;
    reader->size = buf ? buf->size : 0;
    reader->pos = 0;
}

void carquet_buffer_reader_init_data(carquet_buffer_reader_t* reader,
                                      const uint8_t* data, size_t size) {
    assert(reader != NULL);

    reader->data = data;
    reader->size = size;
    reader->pos = 0;
}

carquet_status_t carquet_buffer_reader_read(carquet_buffer_reader_t* reader,
                                             void* dest, size_t size) {
    assert(reader != NULL);
    assert(dest != NULL);
    if (!carquet_buffer_reader_has(reader, size)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }

    memcpy(dest, reader->data + reader->pos, size);
    reader->pos += size;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_skip(carquet_buffer_reader_t* reader, size_t size) {
    assert(reader != NULL);
    if (!carquet_buffer_reader_has(reader, size)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }

    reader->pos += size;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_byte(carquet_buffer_reader_t* reader,
                                                  uint8_t* value) {
    if (!carquet_buffer_reader_has(reader, 1)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = reader->data[reader->pos++];
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_u16_le(carquet_buffer_reader_t* reader,
                                                    uint16_t* value) {
    if (!carquet_buffer_reader_has(reader, 2)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = carquet_read_u16_le(reader->data + reader->pos);
    reader->pos += 2;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_u32_le(carquet_buffer_reader_t* reader,
                                                    uint32_t* value) {
    if (!carquet_buffer_reader_has(reader, 4)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = carquet_read_u32_le(reader->data + reader->pos);
    reader->pos += 4;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_u64_le(carquet_buffer_reader_t* reader,
                                                    uint64_t* value) {
    if (!carquet_buffer_reader_has(reader, 8)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = carquet_read_u64_le(reader->data + reader->pos);
    reader->pos += 8;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_f32_le(carquet_buffer_reader_t* reader,
                                                    float* value) {
    if (!carquet_buffer_reader_has(reader, 4)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = carquet_read_f32_le(reader->data + reader->pos);
    reader->pos += 4;
    return CARQUET_OK;
}

carquet_status_t carquet_buffer_reader_read_f64_le(carquet_buffer_reader_t* reader,
                                                    double* value) {
    if (!carquet_buffer_reader_has(reader, 8)) {
        return CARQUET_ERROR_FILE_TRUNCATED;
    }
    *value = carquet_read_f64_le(reader->data + reader->pos);
    reader->pos += 8;
    return CARQUET_OK;
}

/* ============================================================================
 * Utility Operations
 * ============================================================================
 */

uint8_t* carquet_buffer_detach(carquet_buffer_t* buf, size_t* size_out) {
    assert(buf != NULL);

    uint8_t* data = buf->data;
    if (size_out) {
        *size_out = buf->size;
    }

    buf->data = NULL;
    buf->size = 0;
    buf->capacity = 0;
    buf->owns_data = true;

    return data;
}

void carquet_buffer_swap(carquet_buffer_t* a, carquet_buffer_t* b) {
    assert(a != NULL);
    assert(b != NULL);

    carquet_buffer_t tmp = *a;
    *a = *b;
    *b = tmp;
}
