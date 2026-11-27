/**
 * @file dictionary.c
 * @brief Dictionary encoding implementation
 *
 * Dictionary encoding stores unique values in a dictionary page,
 * and data pages contain RLE-encoded indices into the dictionary.
 */

#include "core/allocator.h"
#include <carquet/error.h>
#include <carquet/types.h>
#include "rle.h"
#include "core/buffer.h"
#include "core/endian.h"
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/* SIMD-dispatched gather functions for dictionary lookups */
extern void carquet_dispatch_gather_i32(const int32_t* dict, const uint32_t* indices, int64_t count, int32_t* output);
extern void carquet_dispatch_gather_i64(const int64_t* dict, const uint32_t* indices, int64_t count, int64_t* output);
extern void carquet_dispatch_gather_float(const float* dict, const uint32_t* indices, int64_t count, float* output);
extern void carquet_dispatch_gather_double(const double* dict, const uint32_t* indices, int64_t count, double* output);
extern bool carquet_dispatch_checked_gather_i32(const int32_t* dict, int32_t dict_count,
                                                 const uint32_t* indices, int64_t count, int32_t* output);
extern bool carquet_dispatch_checked_gather_i64(const int64_t* dict, int32_t dict_count,
                                                 const uint32_t* indices, int64_t count, int64_t* output);
extern bool carquet_dispatch_checked_gather_float(const float* dict, int32_t dict_count,
                                                   const uint32_t* indices, int64_t count, float* output);
extern bool carquet_dispatch_checked_gather_double(const double* dict, int32_t dict_count,
                                                    const uint32_t* indices, int64_t count, double* output);

/**
 * Ensure dict_data is aligned for type T before casting.
 * If misaligned, copies into a temporary aligned buffer.
 * Sets 'aligned_ptr' to the aligned pointer (type T*) and
 * 'aligned_buf' to the temp allocation (NULL if no copy needed).
 */
#define ENSURE_DICT_ALIGNED(dict_data, dict_bytes, T, aligned_ptr, aligned_buf) \
    do { \
        if (((uintptr_t)(dict_data)) % _Alignof(T) == 0) { \
            (aligned_ptr) = (const T*)(dict_data); \
            (aligned_buf) = NULL; \
        } else { \
            (aligned_buf) = carquet_mem_malloc(dict_bytes); \
            if (!(aligned_buf)) { \
                (aligned_ptr) = NULL; \
            } else { \
                memcpy((aligned_buf), (dict_data), (dict_bytes)); \
                (aligned_ptr) = (const T*)(aligned_buf); \
            } \
        } \
    } while (0)

/* ============================================================================
 * Dictionary Builder
 * ============================================================================
 */

typedef struct dict_entry {
    uint8_t* data;
    size_t size;
    uint32_t hash;
    uint32_t index;
    struct dict_entry* next;
} dict_entry_t;

typedef struct {
    dict_entry_t** buckets;
    size_t num_buckets;
    size_t count;

    carquet_buffer_t dict_buffer;  /* Stores dictionary values */
    uint32_t* indices;             /* Maps input index to dict index */
    size_t indices_capacity;
    size_t indices_count;

    size_t value_size;             /* For fixed-size types */
    bool is_variable_length;

    /* Early-abort cap. When max_dict_bytes is non-zero and the PLAIN
     * dictionary payload grows past it, the builder stops admitting new
     * entries and sets `abandoned` so the caller can bail out to PLAIN
     * without scanning the rest of the input or serializing indices. */
    size_t max_dict_bytes;
    bool abandoned;
} dict_builder_t;

#define DICT_BUILDER_INITIAL_BUCKETS 1024U  /* Must be power of 2 */
#define DICT_BUILDER_MAX_LOAD_NUM 3U
#define DICT_BUILDER_MAX_LOAD_DEN 4U

static carquet_status_t dict_builder_rehash(dict_builder_t* builder, size_t new_bucket_count) {
    dict_entry_t** new_buckets = carquet_mem_calloc(new_bucket_count, sizeof(dict_entry_t*));
    if (!new_buckets) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    /* Use bitmask instead of modulo (new_bucket_count is always power of 2) */
    size_t mask = new_bucket_count - 1;
    for (size_t i = 0; i < builder->num_buckets; i++) {
        dict_entry_t* entry = builder->buckets[i];
        while (entry) {
            dict_entry_t* next = entry->next;
            size_t bucket = entry->hash & mask;
            entry->next = new_buckets[bucket];
            new_buckets[bucket] = entry;
            entry = next;
        }
    }

    carquet_mem_free(builder->buckets);
    builder->buckets = new_buckets;
    builder->num_buckets = new_bucket_count;
    return CARQUET_OK;
}

/**
 * Fast hash for fixed-size values using murmur3-style finalizer.
 * Much faster than FNV-1a for 4/8 byte values because it avoids
 * the per-byte sequential dependency chain.
 */
static inline uint32_t dict_hash_fixed32(const uint8_t* data) {
    uint32_t h;
    memcpy(&h, data, 4);
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

static inline uint32_t dict_hash_fixed64(const uint8_t* data) {
    uint64_t k;
    memcpy(&k, data, 8);
    /* murmur3-style 64-to-32 mix */
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return (uint32_t)k;
}

static uint32_t dict_hash(const uint8_t* data, size_t size) {
    /* Fast paths for common fixed-size types */
    if (size == 4) return dict_hash_fixed32(data);
    if (size == 8) return dict_hash_fixed64(data);

    /* FNV-1a for variable-length data */
    uint32_t h = 0x811c9dc5;
    for (size_t i = 0; i < size; i++) {
        h ^= data[i];
        h *= 0x01000193;
    }
    return h;
}

static carquet_status_t dict_builder_init(dict_builder_t* builder,
                                           size_t expected_count,
                                           size_t value_size,
                                           bool is_variable_length) {
    memset(builder, 0, sizeof(*builder));

    builder->num_buckets = DICT_BUILDER_INITIAL_BUCKETS;
    builder->buckets = carquet_mem_calloc(builder->num_buckets, sizeof(dict_entry_t*));
    if (!builder->buckets) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    carquet_status_t status = carquet_buffer_init_capacity(&builder->dict_buffer, 4096);
    if (status != CARQUET_OK) {
        carquet_mem_free(builder->buckets);
        return status;
    }

    builder->indices_capacity = expected_count > 0 ? expected_count : 1024;
    builder->indices = carquet_mem_malloc(builder->indices_capacity * sizeof(uint32_t));
    if (!builder->indices) {
        carquet_buffer_destroy(&builder->dict_buffer);
        carquet_mem_free(builder->buckets);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    builder->value_size = value_size;
    builder->is_variable_length = is_variable_length;

    return CARQUET_OK;
}

static void dict_builder_destroy(dict_builder_t* builder) {
    if (builder->buckets) {
        for (size_t i = 0; i < builder->num_buckets; i++) {
            dict_entry_t* entry = builder->buckets[i];
            while (entry) {
                dict_entry_t* next = entry->next;
                carquet_mem_free(entry);
                entry = next;
            }
        }
        carquet_mem_free(builder->buckets);
    }
    carquet_buffer_destroy(&builder->dict_buffer);
    carquet_mem_free(builder->indices);
}

static carquet_status_t dict_builder_add(dict_builder_t* builder,
                                          const uint8_t* value,
                                          size_t value_size) {
    /* Ensure indices array has space */
    if (builder->indices_count >= builder->indices_capacity) {
        size_t new_cap = builder->indices_capacity * 2;
        uint32_t* new_indices = carquet_mem_realloc(builder->indices, new_cap * sizeof(uint32_t));
        if (!new_indices) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        builder->indices = new_indices;
        builder->indices_capacity = new_cap;
    }

    /* Look up in hash table (bitmask since num_buckets is power of 2) */
    uint32_t hash = dict_hash(value, value_size);
    size_t mask = builder->num_buckets - 1;
    size_t bucket = hash & mask;

    for (dict_entry_t* entry = builder->buckets[bucket]; entry; entry = entry->next) {
        if (entry->hash == hash &&
            entry->size == value_size &&
            memcmp(entry->data, value, value_size) == 0) {
            /* Found existing entry */
            builder->indices[builder->indices_count++] = entry->index;
            return CARQUET_OK;
        }
    }

    if ((builder->count + 1) * DICT_BUILDER_MAX_LOAD_DEN >
        builder->num_buckets * DICT_BUILDER_MAX_LOAD_NUM) {
        carquet_status_t status = dict_builder_rehash(builder, builder->num_buckets * 2);
        if (status != CARQUET_OK) {
            return status;
        }
        mask = builder->num_buckets - 1;
        bucket = hash & mask;
    }

    /* Add new entry */
    dict_entry_t* new_entry = carquet_mem_malloc(sizeof(dict_entry_t) + value_size);
    if (!new_entry) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    new_entry->data = (uint8_t*)(new_entry + 1);
    memcpy(new_entry->data, value, value_size);
    new_entry->size = value_size;
    new_entry->hash = hash;
    new_entry->index = (uint32_t)builder->count;
    new_entry->next = builder->buckets[bucket];
    builder->buckets[bucket] = new_entry;

    /* Add to dictionary buffer. An append failure must not leave an indexed
     * entry without its bytes (the indices would reference a missing value). */
    carquet_status_t append_status = CARQUET_OK;
    if (builder->is_variable_length) {
        /* Write length prefix */
        uint32_t len = (uint32_t)value_size;
        append_status = carquet_buffer_append_u32_le(&builder->dict_buffer, len);
    }
    if (append_status == CARQUET_OK) {
        append_status = carquet_buffer_append(&builder->dict_buffer, value, value_size);
    }
    if (append_status != CARQUET_OK) {
        builder->buckets[bucket] = new_entry->next;
        carquet_mem_free(new_entry);
        return append_status;
    }

    builder->indices[builder->indices_count++] = new_entry->index;
    builder->count++;

    /* Early-abort: the dictionary is no longer worthwhile once its PLAIN
     * payload exceeds the budget. Stop here; the caller detects `abandoned`
     * and falls back to PLAIN without touching the rest of the input. */
    if (builder->max_dict_bytes &&
        builder->dict_buffer.size > builder->max_dict_bytes) {
        builder->abandoned = true;
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Dictionary Encoding
 * ============================================================================
 */

static int bit_width_for_count(uint32_t count) {
    if (count == 0) return 0;
    count--; /* Max index */
    int width = 0;
    while (count > 0) {
        width++;
        count >>= 1;
    }
    return width > 0 ? width : 1;
}

carquet_status_t carquet_dictionary_encode_int32(
    const int32_t* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output) {

    dict_builder_t builder;
    carquet_status_t status = dict_builder_init(&builder, (size_t)count, sizeof(int32_t), false);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Build dictionary - convert each value to little-endian for Parquet format */
    for (int64_t i = 0; i < count; i++) {
        uint8_t le_bytes[sizeof(int32_t)];
        carquet_write_i32_le(le_bytes, values[i]);
        status = dict_builder_add(&builder, le_bytes, sizeof(int32_t));
        if (status != CARQUET_OK) {
            dict_builder_destroy(&builder);
            return status;
        }
    }

    /* Copy dictionary */
    status = carquet_buffer_append(dict_output, builder.dict_buffer.data, builder.dict_buffer.size);
    if (status != CARQUET_OK) {
        dict_builder_destroy(&builder);
        return status;
    }

    /* Encode indices with RLE */
    int bit_width = bit_width_for_count((uint32_t)builder.count);

    /* Write bit width byte */
    uint8_t bw = (uint8_t)bit_width;
    status = carquet_buffer_append_byte(indices_output, bw);
    if (status != CARQUET_OK) {
        dict_builder_destroy(&builder);
        return status;
    }

    /* RLE encode indices */
    status = carquet_rle_encode_all(builder.indices, count, bit_width, indices_output);

    dict_builder_destroy(&builder);
    return status;
}

carquet_status_t carquet_dictionary_encode_int64(
    const int64_t* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output) {

    dict_builder_t builder;
    carquet_status_t status = dict_builder_init(&builder, (size_t)count, sizeof(int64_t), false);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Convert each value to little-endian for Parquet format */
    for (int64_t i = 0; i < count; i++) {
        uint8_t le_bytes[sizeof(int64_t)];
        carquet_write_i64_le(le_bytes, values[i]);
        status = dict_builder_add(&builder, le_bytes, sizeof(int64_t));
        if (status != CARQUET_OK) {
            dict_builder_destroy(&builder);
            return status;
        }
    }

    status = carquet_buffer_append(dict_output, builder.dict_buffer.data, builder.dict_buffer.size);
    if (status == CARQUET_OK) {
        int bit_width = bit_width_for_count((uint32_t)builder.count);
        status = carquet_buffer_append_byte(indices_output, (uint8_t)bit_width);
        if (status == CARQUET_OK) {
            status = carquet_rle_encode_all(builder.indices, count, bit_width, indices_output);
        }
    }

    dict_builder_destroy(&builder);
    return status;
}

carquet_status_t carquet_dictionary_encode_float(
    const float* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output) {

    dict_builder_t builder;
    carquet_status_t status = dict_builder_init(&builder, (size_t)count, sizeof(float), false);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Convert each value to little-endian for Parquet format */
    for (int64_t i = 0; i < count; i++) {
        uint8_t le_bytes[sizeof(float)];
        carquet_write_f32_le(le_bytes, values[i]);
        status = dict_builder_add(&builder, le_bytes, sizeof(float));
        if (status != CARQUET_OK) {
            dict_builder_destroy(&builder);
            return status;
        }
    }

    status = carquet_buffer_append(dict_output, builder.dict_buffer.data, builder.dict_buffer.size);
    if (status == CARQUET_OK) {
        int bit_width = bit_width_for_count((uint32_t)builder.count);
        status = carquet_buffer_append_byte(indices_output, (uint8_t)bit_width);
        if (status == CARQUET_OK) {
            status = carquet_rle_encode_all(builder.indices, count, bit_width, indices_output);
        }
    }

    dict_builder_destroy(&builder);
    return status;
}

carquet_status_t carquet_dictionary_encode_double(
    const double* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output) {

    dict_builder_t builder;
    carquet_status_t status = dict_builder_init(&builder, (size_t)count, sizeof(double), false);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Convert each value to little-endian for Parquet format */
    for (int64_t i = 0; i < count; i++) {
        uint8_t le_bytes[sizeof(double)];
        carquet_write_f64_le(le_bytes, values[i]);
        status = dict_builder_add(&builder, le_bytes, sizeof(double));
        if (status != CARQUET_OK) {
            dict_builder_destroy(&builder);
            return status;
        }
    }

    status = carquet_buffer_append(dict_output, builder.dict_buffer.data, builder.dict_buffer.size);
    if (status == CARQUET_OK) {
        int bit_width = bit_width_for_count((uint32_t)builder.count);
        status = carquet_buffer_append_byte(indices_output, (uint8_t)bit_width);
        if (status == CARQUET_OK) {
            status = carquet_rle_encode_all(builder.indices, count, bit_width, indices_output);
        }
    }

    dict_builder_destroy(&builder);
    return status;
}

carquet_status_t carquet_dictionary_encode_byte_array(
    const carquet_byte_array_t* values,
    int64_t count,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output) {

    dict_builder_t builder;
    carquet_status_t status = dict_builder_init(&builder, (size_t)count, 0, true);
    if (status != CARQUET_OK) {
        return status;
    }

    for (int64_t i = 0; i < count; i++) {
        if (values[i].length < 0 || (!values[i].data && values[i].length > 0)) {
            dict_builder_destroy(&builder);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        status = dict_builder_add(&builder, values[i].data, (size_t)values[i].length);
        if (status != CARQUET_OK) {
            dict_builder_destroy(&builder);
            return status;
        }
    }

    status = carquet_buffer_append(dict_output, builder.dict_buffer.data, builder.dict_buffer.size);
    if (status == CARQUET_OK) {
        int bit_width = bit_width_for_count((uint32_t)builder.count);
        status = carquet_buffer_append_byte(indices_output, (uint8_t)bit_width);
        if (status == CARQUET_OK) {
            status = carquet_rle_encode_all(builder.indices, count, bit_width, indices_output);
        }
    }

    dict_builder_destroy(&builder);
    return status;
}

/* ============================================================================
 * Capped dictionary build (the column writer's path)
 * ============================================================================
 */

/* Fixed-width values (INT32/INT64/FLOAT/DOUBLE/FIXED_LEN_BYTE_ARRAY) are
 * deduplicated in an open-addressing table whose slots hold entry index + 1
 * (0 = empty). The entries themselves are the PLAIN dictionary payload being
 * built, so a lookup compares against payload bytes: no per-entry allocation,
 * no chain walk, no memcmp call for 4- and 8-byte values. Entries are
 * numbered in first-occurrence order, exactly like the chained builder. */
typedef struct {
    uint32_t* slots;
    size_t mask;
} fixed_table_t;

static inline uint32_t fixed_hash(const uint8_t* v, size_t size) {
    if (size == 4) return dict_hash_fixed32(v);
    if (size == 8) return dict_hash_fixed64(v);
    return dict_hash(v, size);
}

static inline bool fixed_equal(const uint8_t* a, const uint8_t* b, size_t size) {
    if (size == 4) {
        uint32_t x, y;
        memcpy(&x, a, 4); memcpy(&y, b, 4);
        return x == y;
    }
    if (size == 8) {
        uint64_t x, y;
        memcpy(&x, a, 8); memcpy(&y, b, 8);
        return x == y;
    }
    return memcmp(a, b, size) == 0;
}

static carquet_status_t fixed_table_grow(fixed_table_t* t, const uint8_t* entries,
                                         uint32_t count, size_t size) {
    size_t cap = t->slots ? (t->mask + 1) * 2 : 1024;
    uint32_t* slots = carquet_mem_calloc(cap, sizeof(uint32_t));
    if (!slots) return CARQUET_ERROR_OUT_OF_MEMORY;
    size_t mask = cap - 1;
    for (uint32_t e = 0; e < count; e++) {
        size_t b = fixed_hash(entries + (size_t)e * size, size) & mask;
        while (slots[b]) b = (b + 1) & mask;
        slots[b] = e + 1;
    }
    carquet_mem_free(t->slots);
    t->slots = slots;
    t->mask = mask;
    return CARQUET_OK;
}

/* Bring value i into PLAIN (little-endian) form. */
static inline const uint8_t* fixed_value_le(carquet_physical_type_t type,
                                            const void* values, int64_t i,
                                            size_t size, uint8_t* tmp) {
    switch (type) {
        case CARQUET_PHYSICAL_INT32:
            carquet_write_i32_le(tmp, ((const int32_t*)values)[i]); return tmp;
        case CARQUET_PHYSICAL_INT64:
            carquet_write_i64_le(tmp, ((const int64_t*)values)[i]); return tmp;
        case CARQUET_PHYSICAL_FLOAT:
            carquet_write_f32_le(tmp, ((const float*)values)[i]); return tmp;
        case CARQUET_PHYSICAL_DOUBLE:
            carquet_write_f64_le(tmp, ((const double*)values)[i]); return tmp;
        default:
            return (const uint8_t*)values + (size_t)i * size;
    }
}

#if CARQUET_LITTLE_ENDIAN
/* 4- and 8-byte values on a little-endian host, where a value's bytes in
 * memory already are its PLAIN form: each value is one machine word, hashed
 * with a multiplicative (Fibonacci) hash taking the top bits of the product,
 * so every input bit reaches the slot index. Entries live in a word array
 * (no unaligned payload loads), and a value equal to the previous one reuses
 * its index without a lookup. Same first-occurrence numbering and the same
 * early-abort rule as the generic loop below, so the same dictionary and
 * indices. */
static inline uint64_t word_load(const uint8_t* p, size_t size) {
    if (size == 4) {
        uint32_t v;
        memcpy(&v, p, 4);
        return v;
    }
    uint64_t v;
    memcpy(&v, p, 8);
    return v;
}

/* The xor-shift pre-mix brings the high bits (exponent and leading mantissa
 * bits of a float/double, which are often the only varying bits) down
 * before the multiply. Probe sequences must stay short: each extra probe is
 * a mispredicted branch after a dependent load, which costs more than the
 * hash. */
static inline size_t word_slot(uint64_t x, size_t size, unsigned bits) {
    if (size == 4) {
        uint32_t h = (uint32_t)x;
        h ^= h >> 16;
        return (size_t)((h * 0x9E3779B1u) >> (32 - bits));
    }
    x ^= x >> 29;
    x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 32;
    return (size_t)((x * 0x9E3779B97F4A7C15ULL) >> (64 - bits));
}

static inline carquet_status_t build_fixed_word(size_t size, const uint8_t* src,
                                                int64_t count, size_t max_dict_bytes,
                                                carquet_buffer_t* dict_output,
                                                uint32_t* indices, uint32_t* num_entries,
                                                bool* abandoned) {
    unsigned bits = 10;
    size_t mask = ((size_t)1 << bits) - 1;
    uint32_t* slots = carquet_mem_calloc(mask + 1, sizeof(uint32_t));
    size_t entries_cap = 256;
    uint64_t* entries = carquet_mem_malloc(entries_cap * sizeof(uint64_t));
    carquet_status_t status = CARQUET_OK;
    uint32_t n = 0;
    if (!slots || !entries) {
        status = CARQUET_ERROR_OUT_OF_MEMORY;
        count = 0;
    }
    uint64_t prev = 0;
    uint32_t prev_index = UINT32_MAX;   /* none yet */
    for (int64_t i = 0; i < count; i++) {
        uint64_t x = word_load(src + (size_t)i * size, size);
        if (x == prev && prev_index != UINT32_MAX) {
            indices[i] = prev_index;
            continue;
        }
        size_t b = word_slot(x, size, bits);
        uint32_t slot;
        while ((slot = slots[b]) != 0 && entries[slot - 1] != x) {
            b = (b + 1) & mask;
        }
        if (!slot) {
            if (n == UINT32_MAX) { status = CARQUET_ERROR_ENCODE; break; }
            if (n == entries_cap) {
                uint64_t* grown = (entries_cap <= SIZE_MAX / (2 * sizeof(uint64_t)))
                    ? carquet_mem_realloc(entries, entries_cap * 2 * sizeof(uint64_t))
                    : NULL;
                if (!grown) { status = CARQUET_ERROR_OUT_OF_MEMORY; break; }
                entries = grown;
                entries_cap *= 2;
            }
            entries[n] = x;
            slots[b] = n + 1;
            slot = ++n;
            /* Same budget rule as the generic loop. */
            if (max_dict_bytes && (size_t)n * size > max_dict_bytes) {
                *abandoned = true;
                break;
            }
            /* Load factor at most 1/8 while the table is small (up to 64K
             * slots, 256 KB), then at most 1/2. */
            if ((size_t)n * ((mask + 1) <= ((size_t)1 << 16) ? 8 : 2) > mask + 1) {
                unsigned nbits = bits + 1;
                size_t nmask = ((size_t)1 << nbits) - 1;
                uint32_t* nslots = carquet_mem_calloc(nmask + 1, sizeof(uint32_t));
                if (!nslots) { status = CARQUET_ERROR_OUT_OF_MEMORY; break; }
                for (uint32_t e = 0; e < n; e++) {
                    size_t nb = word_slot(entries[e], size, nbits);
                    while (nslots[nb]) nb = (nb + 1) & nmask;
                    nslots[nb] = e + 1;
                }
                carquet_mem_free(slots);
                slots = nslots;
                bits = nbits;
                mask = nmask;
            }
        }
        indices[i] = prev_index = slot - 1;
        prev = x;
    }
    if (status == CARQUET_OK && !*abandoned && n > 0) {
        size_t base = dict_output->size;
        status = carquet_buffer_resize(dict_output, base + (size_t)n * size);
        if (status == CARQUET_OK) {
            uint8_t* out = dict_output->data + base;
            if (size == 8) {
                memcpy(out, entries, (size_t)n * 8);
            } else {
                for (uint32_t e = 0; e < n; e++) {
                    uint32_t v = (uint32_t)entries[e];
                    memcpy(out + (size_t)e * 4, &v, 4);
                }
            }
        }
    }
    carquet_mem_free(slots);
    carquet_mem_free(entries);
    *num_entries = n;
    return status;
}
#endif

static carquet_status_t build_fixed(carquet_physical_type_t type, size_t size,
                                    const void* values, int64_t count,
                                    size_t max_dict_bytes,
                                    carquet_buffer_t* dict_output,
                                    uint32_t* indices, uint32_t* num_entries,
                                    bool* abandoned) {
#if CARQUET_LITTLE_ENDIAN
    if (size == 4) {
        return build_fixed_word(4, (const uint8_t*)values, count, max_dict_bytes,
                                dict_output, indices, num_entries, abandoned);
    }
    if (size == 8) {
        return build_fixed_word(8, (const uint8_t*)values, count, max_dict_bytes,
                                dict_output, indices, num_entries, abandoned);
    }
#endif
    fixed_table_t t = {NULL, 0};
    uint32_t n = 0;
    size_t base = dict_output->size;
    carquet_status_t status = fixed_table_grow(&t, NULL, 0, size);
    for (int64_t i = 0; status == CARQUET_OK && i < count; i++) {
        uint8_t tmp[8];
        const uint8_t* v = fixed_value_le(type, values, i, size, tmp);
        size_t b = fixed_hash(v, size) & t.mask;
        const uint8_t* entries = dict_output->data + base;
        uint32_t slot;
        while ((slot = t.slots[b]) != 0 &&
               !fixed_equal(entries + (size_t)(slot - 1) * size, v, size)) {
            b = (b + 1) & t.mask;
        }
        if (slot) {
            indices[i] = slot - 1;
            continue;
        }
        if (n == UINT32_MAX) { status = CARQUET_ERROR_ENCODE; break; }
        status = carquet_buffer_append(dict_output, v, size);
        if (status != CARQUET_OK) break;
        t.slots[b] = n + 1;
        indices[i] = n++;
        /* Same budget rule as the chained builder: stop as soon as the PLAIN
         * payload exceeds the cap. */
        if (max_dict_bytes && (size_t)n * size > max_dict_bytes) {
            *abandoned = true;
            break;
        }
        /* Keep the load factor at or below one half. */
        if ((size_t)n * 2 > t.mask + 1) {
            status = fixed_table_grow(&t, dict_output->data + base, n, size);
        }
    }
    carquet_mem_free(t.slots);
    *num_entries = n;
    return status;
}

/* Single-pass dictionary build with an early-abort budget, applying the same
 * PLAIN value marshaling as the per-type encoders above. On success the PLAIN
 * dictionary payload is appended to dict_output, *indices_out receives one
 * entry index per input value (count entries, carquet_mem_free to release)
 * and *num_entries the dictionary size. If the payload would exceed
 * max_dict_bytes, it stops immediately (without scanning the remaining input)
 * and reports *abandoned = true, with no indices, so the caller can fall back
 * to PLAIN. When max_dict_bytes is 0 the cap is disabled. */
carquet_status_t carquet_dictionary_build_capped(
    carquet_physical_type_t type,
    int32_t type_length,
    const void* fixed_values,
    const carquet_byte_array_t* ba_values,
    int64_t count,
    size_t max_dict_bytes,
    carquet_buffer_t* dict_output,
    uint32_t** indices_out,
    int32_t* num_entries,
    bool* abandoned) {

    if (abandoned) *abandoned = false;
    if (indices_out) *indices_out = NULL;
    if (num_entries) *num_entries = 0;
    if (count <= 0 || !dict_output || !indices_out || !num_entries ||
        (uint64_t)count > (uint64_t)(SIZE_MAX / sizeof(uint32_t))) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    size_t value_size;
    switch (type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:  value_size = 4; break;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE: value_size = 8; break;
        case CARQUET_PHYSICAL_BYTE_ARRAY: value_size = 0; break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (type_length <= 0) return CARQUET_ERROR_INVALID_ARGUMENT;
            value_size = (size_t)type_length;
            break;
        default: return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
    if (type == CARQUET_PHYSICAL_BYTE_ARRAY ? !ba_values : !fixed_values) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    bool dropped = false;
    carquet_status_t status;
    if (type != CARQUET_PHYSICAL_BYTE_ARRAY) {
        uint32_t* indices = carquet_mem_malloc((size_t)count * sizeof(uint32_t));
        if (!indices) return CARQUET_ERROR_OUT_OF_MEMORY;
        uint32_t n = 0;
        size_t base = dict_output->size;
        status = build_fixed(type, value_size, fixed_values, count, max_dict_bytes,
                             dict_output, indices, &n, &dropped);
        if (status != CARQUET_OK || dropped || n > INT32_MAX) {
            carquet_mem_free(indices);
            dict_output->size = base;
            if (status == CARQUET_OK && !dropped) status = CARQUET_ERROR_ENCODE;
            if (dropped && abandoned) *abandoned = true;
            return status;
        }
        *indices_out = indices;
        *num_entries = (int32_t)n;
        return CARQUET_OK;
    }

    dict_builder_t builder;
    status = dict_builder_init(&builder, (size_t)count, 0, true);
    if (status != CARQUET_OK) return status;
    builder.max_dict_bytes = max_dict_bytes;
    for (int64_t i = 0; i < count && !builder.abandoned; i++) {
        if (ba_values[i].length < 0 ||
            (!ba_values[i].data && ba_values[i].length > 0)) {
            status = CARQUET_ERROR_INVALID_ARGUMENT;
            break;
        }
        status = dict_builder_add(&builder, ba_values[i].data,
                                  (size_t)ba_values[i].length);
        if (status != CARQUET_OK) break;
    }
    if (status == CARQUET_OK && builder.abandoned) {
        if (abandoned) *abandoned = true;
    } else if (status == CARQUET_OK && builder.count > INT32_MAX) {
        status = CARQUET_ERROR_ENCODE;
    } else if (status == CARQUET_OK) {
        status = carquet_buffer_append(dict_output, builder.dict_buffer.data,
                                       builder.dict_buffer.size);
        if (status == CARQUET_OK) {
            /* Hand the index array over instead of copying it. */
            *indices_out = builder.indices;
            *num_entries = (int32_t)builder.count;
            builder.indices = NULL;
        }
    }
    dict_builder_destroy(&builder);
    return status;
}

/* Bit width of the indices of a dictionary of `count` entries. */
int carquet_dictionary_index_bit_width(int32_t count) {
    return bit_width_for_count(count > 0 ? (uint32_t)count : 0);
}

/* carquet_dictionary_build_capped, then the whole index stream encoded as one
 * RLE_DICTIONARY data page payload (bit-width byte + RLE/bit-packed runs). */
carquet_status_t carquet_dictionary_encode_capped(
    carquet_physical_type_t type,
    int32_t type_length,
    const void* fixed_values,
    const carquet_byte_array_t* ba_values,
    int64_t count,
    size_t max_dict_bytes,
    carquet_buffer_t* dict_output,
    carquet_buffer_t* indices_output,
    bool* abandoned) {

    if (abandoned) *abandoned = false;
    if (count <= 0) return CARQUET_ERROR_INVALID_ARGUMENT;
    uint32_t* indices = NULL;
    int32_t n = 0;
    bool dropped = false;
    carquet_status_t status = carquet_dictionary_build_capped(
        type, type_length, fixed_values, ba_values, count, max_dict_bytes,
        dict_output, &indices, &n, &dropped);
    if (status != CARQUET_OK) return status;
    if (dropped) {
        if (abandoned) *abandoned = true;
        return CARQUET_OK;
    }
    int bit_width = bit_width_for_count((uint32_t)n);
    status = carquet_buffer_append_byte(indices_output, (uint8_t)bit_width);
    if (status == CARQUET_OK) {
        status = carquet_rle_encode_all(indices, count, bit_width, indices_output);
    }
    carquet_mem_free(indices);
    return status;
}

/* output_count is an int64_t from the caller; on a 32-bit target
 * `output_count * sizeof(uint32_t)` can exceed SIZE_MAX and the converted
 * malloc size wraps to a small buffer that the RLE index decode then
 * overruns. Reject counts whose index array cannot be addressed. */
static inline bool indices_count_fits(int64_t output_count) {
    return output_count >= 0 &&
           (uint64_t)output_count <= (uint64_t)(SIZE_MAX / sizeof(uint32_t));
}

/* ============================================================================
 * Dictionary Decoding
 * ============================================================================
 */

carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    int32_t* output,
    int64_t output_count) {

    /* Early validation */
    if (output_count <= 0) {
        return CARQUET_OK;
    }

    if (dict_count <= 0 || dict_data == NULL) {
        return CARQUET_ERROR_DECODE;
    }

    if (dict_size / sizeof(int32_t) < (size_t)dict_count) {  /* no 32-bit wrap */
        return CARQUET_ERROR_DECODE;
    }

    /* Read bit width */
    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }
    int bit_width = indices_data[0];

    /* Decode RLE indices */
    if (!indices_count_fits(output_count)) {
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices = carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);

    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }

    /* Parquet stores dictionary values in little-endian format.
     * On little-endian systems (all x86, all modern ARM), we can cast
     * dict_data directly to int32_t* and use SIMD gather. */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    /* Big-endian: use scalar path with endian conversion */
    for (int64_t i = 0; i < decoded; i++) {
        if (indices[i] >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            return CARQUET_ERROR_DECODE;  /* index past the dictionary */
        }
        output[i] = carquet_read_i32_le(dict_data + indices[i] * sizeof(int32_t));
    }
#else
    const int32_t* aligned_dict;
    void* aligned_buf;
    size_t dict_bytes = (size_t)dict_count * sizeof(int32_t);
    ENSURE_DICT_ALIGNED(dict_data, dict_bytes, int32_t, aligned_dict, aligned_buf);
    if (!aligned_dict) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (!carquet_dispatch_checked_gather_i32(aligned_dict, dict_count,
                                             indices, decoded, output)) {
        carquet_mem_free(aligned_buf);
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }
    carquet_mem_free(aligned_buf);
#endif

    carquet_mem_free(indices);
    return CARQUET_OK;
}

carquet_status_t carquet_dictionary_decode_int64(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    int64_t* output,
    int64_t output_count) {

    /* Early validation */
    if (output_count <= 0) {
        return CARQUET_OK;
    }

    if (dict_count <= 0 || dict_data == NULL) {
        return CARQUET_ERROR_DECODE;
    }

    if (dict_size / sizeof(int64_t) < (size_t)dict_count) {  /* no 32-bit wrap */
        return CARQUET_ERROR_DECODE;
    }

    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }
    int bit_width = indices_data[0];

    if (!indices_count_fits(output_count)) {
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices = carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);

    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }

    /* Parquet stores dictionary values in little-endian format.
     * On little-endian systems, we can cast dict_data directly and use SIMD gather. */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    /* Big-endian: use scalar path with endian conversion */
    for (int64_t i = 0; i < decoded; i++) {
        if (indices[i] >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            return CARQUET_ERROR_DECODE;  /* index past the dictionary */
        }
        output[i] = carquet_read_i64_le(dict_data + indices[i] * sizeof(int64_t));
    }
#else
    const int64_t* aligned_dict;
    void* aligned_buf;
    size_t dict_bytes = (size_t)dict_count * sizeof(int64_t);
    ENSURE_DICT_ALIGNED(dict_data, dict_bytes, int64_t, aligned_dict, aligned_buf);
    if (!aligned_dict) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (!carquet_dispatch_checked_gather_i64(aligned_dict, dict_count,
                                             indices, decoded, output)) {
        carquet_mem_free(aligned_buf);
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }
    carquet_mem_free(aligned_buf);
#endif

    carquet_mem_free(indices);
    return CARQUET_OK;
}

carquet_status_t carquet_dictionary_decode_float(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    float* output,
    int64_t output_count) {

    /* Early validation */
    if (output_count <= 0) {
        return CARQUET_OK;
    }

    if (dict_count <= 0 || dict_data == NULL) {
        return CARQUET_ERROR_DECODE;
    }

    if (dict_size / sizeof(float) < (size_t)dict_count) {  /* no 32-bit wrap */
        return CARQUET_ERROR_DECODE;
    }

    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }
    int bit_width = indices_data[0];

    if (!indices_count_fits(output_count)) {
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices = carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);

    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }

    /* Parquet stores dictionary values in little-endian format.
     * On little-endian systems, we can cast dict_data directly and use SIMD gather. */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    /* Big-endian: use scalar path with endian conversion */
    for (int64_t i = 0; i < decoded; i++) {
        if (indices[i] >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            return CARQUET_ERROR_DECODE;  /* index past the dictionary */
        }
        output[i] = carquet_read_f32_le(dict_data + indices[i] * sizeof(float));
    }
#else
    const float* aligned_dict;
    void* aligned_buf;
    size_t dict_bytes = (size_t)dict_count * sizeof(float);
    ENSURE_DICT_ALIGNED(dict_data, dict_bytes, float, aligned_dict, aligned_buf);
    if (!aligned_dict) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (!carquet_dispatch_checked_gather_float(aligned_dict, dict_count,
                                               indices, decoded, output)) {
        carquet_mem_free(aligned_buf);
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }
    carquet_mem_free(aligned_buf);
#endif

    carquet_mem_free(indices);
    return CARQUET_OK;
}

carquet_status_t carquet_dictionary_decode_double(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    double* output,
    int64_t output_count) {

    /* Early validation */
    if (output_count <= 0) {
        return CARQUET_OK;
    }

    if (dict_count <= 0 || dict_data == NULL) {
        return CARQUET_ERROR_DECODE;
    }

    if (dict_size / sizeof(double) < (size_t)dict_count) {  /* no 32-bit wrap */
        return CARQUET_ERROR_DECODE;
    }

    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }
    int bit_width = indices_data[0];

    if (!indices_count_fits(output_count)) {
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices = carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);

    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }

    /* Parquet stores dictionary values in little-endian format.
     * On little-endian systems, we can cast dict_data directly and use SIMD gather. */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    /* Big-endian: use scalar path with endian conversion */
    for (int64_t i = 0; i < decoded; i++) {
        if (indices[i] >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            return CARQUET_ERROR_DECODE;  /* index past the dictionary */
        }
        output[i] = carquet_read_f64_le(dict_data + indices[i] * sizeof(double));
    }
#else
    const double* aligned_dict;
    void* aligned_buf;
    size_t dict_bytes = (size_t)dict_count * sizeof(double);
    ENSURE_DICT_ALIGNED(dict_data, dict_bytes, double, aligned_dict, aligned_buf);
    if (!aligned_dict) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (!carquet_dispatch_checked_gather_double(aligned_dict, dict_count,
                                                indices, decoded, output)) {
        carquet_mem_free(aligned_buf);
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }
    carquet_mem_free(aligned_buf);
#endif

    carquet_mem_free(indices);
    return CARQUET_OK;
}

/* Decode a BYTE_ARRAY dictionary page + RLE index stream into an array of
 * carquet_byte_array_t. Each output value points into dict_data (the caller
 * must keep the dictionary buffer alive for as long as the output is used).
 * The dictionary uses PLAIN BYTE_ARRAY layout: dict_count records of a 4-byte
 * little-endian length prefix followed by that many payload bytes. Every
 * length is validated against dict_size and every index against dict_count. */
carquet_status_t carquet_dictionary_decode_byte_array(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    const uint8_t* indices_data,
    size_t indices_size,
    carquet_byte_array_t* output,
    int64_t output_count) {

    if (output_count <= 0) {
        return CARQUET_OK;
    }
    if (dict_count <= 0 || dict_data == NULL || output == NULL) {
        return CARQUET_ERROR_DECODE;
    }
    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }

    /* Build an offset/length table over the length-prefixed dictionary,
     * validating that every entry stays within dict_size. Each entry needs
     * at least its 4-byte length prefix, so dict_count is bounded by
     * dict_size / 4 before anything is allocated from it. */
    if ((size_t)dict_count > dict_size / 4) {
        return CARQUET_ERROR_DECODE;
    }
    typedef struct { const uint8_t* data; int32_t length; } ba_entry_t;
    ba_entry_t* entries =
        carquet_mem_malloc((size_t)dict_count * sizeof(ba_entry_t));
    if (!entries) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    size_t pos = 0;
    for (int32_t i = 0; i < dict_count; i++) {
        if (pos + 4 > dict_size) {
            carquet_mem_free(entries);
            return CARQUET_ERROR_DECODE;
        }
        uint32_t len = carquet_read_u32_le(dict_data + pos);
        pos += 4;
        /* len bytes must remain; guard against overflow and truncation. */
        if (len > dict_size - pos || len > (uint32_t)INT32_MAX) {
            carquet_mem_free(entries);
            return CARQUET_ERROR_DECODE;
        }
        entries[i].data = dict_data + pos;
        entries[i].length = (int32_t)len;
        pos += len;
    }

    int bit_width = indices_data[0];
    if (!indices_count_fits(output_count)) {
        carquet_mem_free(entries);
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices =
        carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        carquet_mem_free(entries);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);
    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        carquet_mem_free(entries);
        return CARQUET_ERROR_DECODE;
    }

    for (int64_t i = 0; i < output_count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            carquet_mem_free(entries);
            return CARQUET_ERROR_DECODE;
        }
        /* Cast away const: output is a view into the caller's dictionary. */
        output[i].data = (uint8_t*)entries[idx].data;
        output[i].length = entries[idx].length;
    }

    carquet_mem_free(indices);
    carquet_mem_free(entries);
    return CARQUET_OK;
}

/* Decode a FIXED_LEN_BYTE_ARRAY dictionary page + RLE index stream into a
 * contiguous output buffer of output_count * type_length bytes (PLAIN FLBA
 * layout). The dictionary is dict_count fixed-width values concatenated. */
carquet_status_t carquet_dictionary_decode_fixed_len_byte_array(
    const uint8_t* dict_data,
    size_t dict_size,
    int32_t dict_count,
    int32_t type_length,
    const uint8_t* indices_data,
    size_t indices_size,
    uint8_t* output,
    int64_t output_count) {

    if (output_count <= 0) {
        return CARQUET_OK;
    }
    if (dict_count <= 0 || dict_data == NULL || output == NULL) {
        return CARQUET_ERROR_DECODE;
    }
    if (type_length <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (indices_size < 1) {
        return CARQUET_ERROR_DECODE;
    }

    size_t tl = (size_t)type_length;
    /* dict_size must hold dict_count fixed-width values (overflow-safe). */
    if ((size_t)dict_count > dict_size / tl) {
        return CARQUET_ERROR_DECODE;
    }

    int bit_width = indices_data[0];
    if (!indices_count_fits(output_count)) {
        return CARQUET_ERROR_DECODE;
    }
    uint32_t* indices =
        carquet_mem_malloc((size_t)output_count * sizeof(uint32_t));
    if (!indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t decoded = carquet_rle_decode_all(
        indices_data + 1, indices_size - 1, bit_width, indices, output_count);
    if (decoded < 0 || decoded < output_count) {
        carquet_mem_free(indices);
        return CARQUET_ERROR_DECODE;
    }

    for (int64_t i = 0; i < output_count; i++) {
        uint32_t idx = indices[i];
        if (idx >= (uint32_t)dict_count) {
            carquet_mem_free(indices);
            return CARQUET_ERROR_DECODE;
        }
        memcpy(output + (size_t)i * tl, dict_data + (size_t)idx * tl, tl);
    }

    carquet_mem_free(indices);
    return CARQUET_OK;
}
