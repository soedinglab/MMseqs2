/**
 * @file bloom_filter.c
 * @brief Split Block Bloom Filter implementation for Parquet
 *
 * Parquet uses Split Block Bloom Filters (SBBF) for predicate pushdown.
 * The filter is divided into blocks of 256 bits (32 bytes), with each
 * block containing 8 32-bit words. Insertions set 8 bits using a
 * specific algorithm based on xxHash64.
 *
 * Reference: https://parquet.apache.org/docs/file-format/bloomfilter/
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include <carquet/error.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Constants
 * ============================================================================
 */

#define BLOOM_FILTER_BLOCK_SIZE 32     /* 256 bits = 32 bytes */
#define BLOOM_FILTER_WORDS_PER_BLOCK 8 /* 8 x 32-bit words */

/* Salt values used to generate bit positions within a block */
static const uint32_t SALT[8] = {
    0x47b6137bU, 0x44974d91U, 0x8824ad5bU, 0xa2b7289dU,
    0x705495c7U, 0x2df1424bU, 0x9efc4947U, 0x5c6bfb31U
};

/* xxHash64 function declaration (from xxhash.c) */
extern uint64_t carquet_xxhash64(const void* data, size_t length, uint64_t seed);

/* ============================================================================
 * Bloom Filter Structure
 * ============================================================================
 */

struct carquet_bloom_filter {
    uint8_t* data;           /* Filter bit array */
    size_t num_bytes;        /* Size of data in bytes */
    size_t num_blocks;       /* Number of 256-bit blocks */
    bool owns_data;          /* Whether we should free data */
};

/* ============================================================================
 * Core Bloom Filter Operations
 * ============================================================================
 */

/**
 * Generate block index from hash.
 */
static inline size_t bloom_filter_block_index(uint64_t hash, size_t num_blocks) {
    uint64_t top_bits = hash >> 32;
    return (size_t)((top_bits * (uint64_t)num_blocks) >> 32);
}

/**
 * Set bits in a block using the hash value.
 * Uses the SALT values to generate 8 different bit positions.
 */
static void bloom_filter_block_insert(uint32_t* block, uint64_t hash) {
    uint32_t key = (uint32_t)hash;

    for (int i = 0; i < 8; i++) {
        /* Compute mask from salt * key */
        uint32_t mask = SALT[i] * key;
        /* Use top 5 bits as bit position within the word */
        uint32_t bit_pos = mask >> 27;
        /* Set bit in the corresponding word */
        block[i] |= (1U << bit_pos);
    }
}

/**
 * Check if a value might be in the block.
 */
static bool bloom_filter_block_check(const uint32_t* block, uint64_t hash) {
    uint32_t key = (uint32_t)hash;

    for (int i = 0; i < 8; i++) {
        uint32_t mask = SALT[i] * key;
        uint32_t bit_pos = mask >> 27;
        if ((block[i] & (1U << bit_pos)) == 0) {
            return false;  /* Definitely not present */
        }
    }
    return true;  /* Might be present */
}

/* ============================================================================
 * Bloom Filter Creation and Destruction
 * ============================================================================
 */

carquet_bloom_filter_t* carquet_bloom_filter_create(size_t num_bytes) {
    /* Ensure size is a multiple of block size */
    if (num_bytes < BLOOM_FILTER_BLOCK_SIZE) {
        num_bytes = BLOOM_FILTER_BLOCK_SIZE;
    }
    num_bytes = (num_bytes + BLOOM_FILTER_BLOCK_SIZE - 1) /
                BLOOM_FILTER_BLOCK_SIZE * BLOOM_FILTER_BLOCK_SIZE;

    carquet_bloom_filter_t* filter = carquet_mem_malloc(sizeof(carquet_bloom_filter_t));
    if (!filter) {
        return NULL;
    }

    filter->data = carquet_mem_calloc(num_bytes, 1);
    if (!filter->data) {
        carquet_mem_free(filter);
        return NULL;
    }

    filter->num_bytes = num_bytes;
    filter->num_blocks = num_bytes / BLOOM_FILTER_BLOCK_SIZE;
    filter->owns_data = true;

    return filter;
}

carquet_bloom_filter_t* carquet_bloom_filter_create_with_ndv(
    int64_t ndv,
    double fpp) {

    if (ndv <= 0 || fpp <= 0.0 || fpp >= 1.0) {
        return NULL;
    }

    /* Calculate optimal size in bits:
     * m = -n * ln(p) / (ln(2)^2)
     * where n = number of distinct values, p = false positive probability
     */
    double ln2_squared = 0.4804530139182014246671025263266649717305529515945455;
    double bits = -(double)ndv * log(fpp) / ln2_squared;

    /* Convert to bytes, round up to block size */
    size_t num_bytes = (size_t)(bits / 8.0) + 1;

    return carquet_bloom_filter_create(num_bytes);
}

carquet_bloom_filter_t* carquet_bloom_filter_from_data(
    const uint8_t* data,
    size_t size) {

    if (!data || size < BLOOM_FILTER_BLOCK_SIZE) {
        return NULL;
    }

    /* Ensure size is valid */
    if (size % BLOOM_FILTER_BLOCK_SIZE != 0) {
        return NULL;
    }

    carquet_bloom_filter_t* filter = carquet_mem_malloc(sizeof(carquet_bloom_filter_t));
    if (!filter) {
        return NULL;
    }

    filter->data = carquet_mem_malloc(size);
    if (!filter->data) {
        carquet_mem_free(filter);
        return NULL;
    }

    memcpy(filter->data, data, size);
    filter->num_bytes = size;
    filter->num_blocks = size / BLOOM_FILTER_BLOCK_SIZE;
    filter->owns_data = true;

    return filter;
}

void carquet_bloom_filter_destroy(carquet_bloom_filter_t* filter) {
    if (filter) {
        if (filter->owns_data && filter->data) {
            carquet_mem_free(filter->data);
        }
        carquet_mem_free(filter);
    }
}

/* ============================================================================
 * Bloom Filter Insert Operations
 * ============================================================================
 */

void carquet_bloom_filter_insert_hash(carquet_bloom_filter_t* filter,
                                       uint64_t hash) {
    if (!filter || !filter->data) {
        return;
    }

    size_t block_idx = bloom_filter_block_index(hash, filter->num_blocks);
    uint32_t* block = (uint32_t*)(filter->data + block_idx * BLOOM_FILTER_BLOCK_SIZE);

    bloom_filter_block_insert(block, hash);
}

void carquet_bloom_filter_insert_i32(carquet_bloom_filter_t* filter,
                                      int32_t value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    carquet_bloom_filter_insert_hash(filter, hash);
}

void carquet_bloom_filter_insert_i64(carquet_bloom_filter_t* filter,
                                      int64_t value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    carquet_bloom_filter_insert_hash(filter, hash);
}

void carquet_bloom_filter_insert_float(carquet_bloom_filter_t* filter,
                                        float value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    carquet_bloom_filter_insert_hash(filter, hash);
}

void carquet_bloom_filter_insert_double(carquet_bloom_filter_t* filter,
                                         double value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    carquet_bloom_filter_insert_hash(filter, hash);
}

void carquet_bloom_filter_insert_bytes(carquet_bloom_filter_t* filter,
                                        const uint8_t* data,
                                        size_t len) {
    uint64_t hash = carquet_xxhash64(data, len, 0);
    carquet_bloom_filter_insert_hash(filter, hash);
}

/* ============================================================================
 * Bloom Filter Check Operations
 * ============================================================================
 */

bool carquet_bloom_filter_check_hash(const carquet_bloom_filter_t* filter,
                                      uint64_t hash) {
    if (!filter || !filter->data) {
        return true;  /* Assume present if no filter */
    }

    size_t block_idx = bloom_filter_block_index(hash, filter->num_blocks);
    const uint32_t* block = (const uint32_t*)(filter->data + block_idx * BLOOM_FILTER_BLOCK_SIZE);

    return bloom_filter_block_check(block, hash);
}

bool carquet_bloom_filter_check_i32(const carquet_bloom_filter_t* filter,
                                     int32_t value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    return carquet_bloom_filter_check_hash(filter, hash);
}

bool carquet_bloom_filter_check_i64(const carquet_bloom_filter_t* filter,
                                     int64_t value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    return carquet_bloom_filter_check_hash(filter, hash);
}

bool carquet_bloom_filter_check_float(const carquet_bloom_filter_t* filter,
                                       float value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    return carquet_bloom_filter_check_hash(filter, hash);
}

bool carquet_bloom_filter_check_double(const carquet_bloom_filter_t* filter,
                                        double value) {
    uint64_t hash = carquet_xxhash64(&value, sizeof(value), 0);
    return carquet_bloom_filter_check_hash(filter, hash);
}

bool carquet_bloom_filter_check_bytes(const carquet_bloom_filter_t* filter,
                                       const uint8_t* data,
                                       size_t len) {
    uint64_t hash = carquet_xxhash64(data, len, 0);
    return carquet_bloom_filter_check_hash(filter, hash);
}

/* ============================================================================
 * Bloom Filter Accessors
 * ============================================================================
 */

const uint8_t* carquet_bloom_filter_data(const carquet_bloom_filter_t* filter) {
    return filter ? filter->data : NULL;
}

size_t carquet_bloom_filter_size(const carquet_bloom_filter_t* filter) {
    /* filter is nonnull per API contract */
    return filter->num_bytes;
}

size_t carquet_bloom_filter_num_blocks(const carquet_bloom_filter_t* filter) {
    return filter ? filter->num_blocks : 0;
}

/* ============================================================================
 * Bloom Filter Serialization
 * ============================================================================
 */

carquet_status_t carquet_bloom_filter_write(
    const carquet_bloom_filter_t* filter,
    uint8_t* output,
    size_t output_capacity,
    size_t* bytes_written) {

    if (!filter || !output || !bytes_written) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (output_capacity < filter->num_bytes) {
        return CARQUET_ERROR_ENCODE;
    }

    memcpy(output, filter->data, filter->num_bytes);
    *bytes_written = filter->num_bytes;

    return CARQUET_OK;
}

carquet_status_t carquet_bloom_filter_read(
    carquet_bloom_filter_t** filter_out,
    const uint8_t* data,
    size_t data_size) {

    if (!filter_out || !data) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_bloom_filter_t* filter = carquet_bloom_filter_from_data(data, data_size);
    if (!filter) {
        /* A size that is below one block or not a multiple of 32 bytes is a
         * malformed filter, not an allocation failure. */
        return (data_size < 32 || data_size % 32 != 0)
            ? CARQUET_ERROR_INVALID_METADATA : CARQUET_ERROR_OUT_OF_MEMORY;
    }

    *filter_out = filter;
    return CARQUET_OK;
}

/* ============================================================================
 * Bloom Filter Merge
 * ============================================================================
 */

carquet_status_t carquet_bloom_filter_merge(
    carquet_bloom_filter_t* dest,
    const carquet_bloom_filter_t* src) {

    if (!dest || !src) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (dest->num_bytes != src->num_bytes) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* OR the bit arrays together */
    for (size_t i = 0; i < dest->num_bytes; i++) {
        dest->data[i] |= src->data[i];
    }

    return CARQUET_OK;
}
