/**
 * @file page_index.c
 * @brief Page index (ColumnIndex and OffsetIndex) implementation
 *
 * Page indexes enable predicate pushdown by storing per-page statistics.
 * - ColumnIndex: min/max values and null counts for each page
 * - OffsetIndex: file offset, compressed/uncompressed size for each page
 *
 * Reference: https://parquet.apache.org/docs/file-format/
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include <carquet/error.h>
#include "core/arena.h"
#include "core/buffer.h"
#include "thrift/thrift_encode.h"
#include "thrift/thrift_decode.h"
#include "thrift/parquet_types.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * ColumnIndex Structure
 * ============================================================================
 */

struct carquet_column_index {
    int32_t num_pages;

    /* Per-page null counts */
    int64_t* null_counts;
    int32_t num_null_counts;

    /* Per-page min/max values (packed binary) */
    uint8_t** min_values;
    int32_t* min_value_lens;
    int32_t num_min_values;
    uint8_t** max_values;
    int32_t* max_value_lens;
    int32_t num_max_values;

    /* Per-page null page flags */
    bool* null_pages;
    int32_t num_null_pages;

    /* Boundary order for efficient range queries */
    int32_t boundary_order;  /* 0=UNORDERED, 1=ASCENDING, 2=DESCENDING */
};

/* ============================================================================
 * OffsetIndex Structure
 * ============================================================================
 */

struct carquet_offset_index {
    int32_t num_pages;
    carquet_page_location_t* page_locations;
};

/* ============================================================================
 * Forward Declarations
 * ============================================================================
 */

typedef struct carquet_column_index_builder carquet_column_index_builder_t;
typedef struct carquet_offset_index_builder carquet_offset_index_builder_t;

void carquet_column_index_builder_destroy(carquet_column_index_builder_t* builder);
void carquet_offset_index_builder_destroy(carquet_offset_index_builder_t* builder);

/* ============================================================================
 * Column Index Builder
 * ============================================================================
 */

struct carquet_column_index_builder {
    carquet_physical_type_t type;
    carquet_logical_type_t logical_type;
    int32_t type_length;

    int32_t capacity;
    int32_t num_pages;

    int64_t* null_counts;
    uint8_t** min_values;
    int32_t* min_value_lens;
    uint8_t** max_values;
    int32_t* max_value_lens;
    bool* null_pages;

    int32_t boundary_order;

    /* Set once a non-null page without min/max was added (statistics
     * disabled, undefined sort order, all-NaN page). The spec has no way to
     * say "this page has no bounds", so the whole ColumnIndex is withheld,
     * as parquet-cpp does; the OffsetIndex is unaffected. */
    bool discarded;

    /* Per-page level histograms (Parquet 2.9), flattened page-major:
     * rep_level_histograms[page * rep_hist_len + level]. Allocated lazily on the
     * first add_page that supplies histograms; lengths come from the max
     * rep/def levels (max_level + 1). track_histograms gates emission. */
    bool track_histograms;
    int32_t rep_hist_len;
    int32_t def_hist_len;
    int64_t* rep_level_histograms;
    int64_t* def_level_histograms;
};

/**
 * Create a column index builder.
 */
carquet_column_index_builder_t* carquet_column_index_builder_create(
    carquet_physical_type_t type,
    const carquet_logical_type_t* logical_type,
    int32_t type_length) {

    carquet_column_index_builder_t* builder = carquet_mem_calloc(1, sizeof(*builder));
    if (!builder) return NULL;

    builder->type = type;
    if (logical_type) {
        builder->logical_type = *logical_type;
    }
    builder->type_length = type_length;
    builder->capacity = 16;
    builder->boundary_order = 0;  /* UNORDERED by default */

    builder->null_counts = carquet_mem_calloc(builder->capacity, sizeof(int64_t));
    builder->min_values = carquet_mem_calloc(builder->capacity, sizeof(uint8_t*));
    builder->min_value_lens = carquet_mem_calloc(builder->capacity, sizeof(int32_t));
    builder->max_values = carquet_mem_calloc(builder->capacity, sizeof(uint8_t*));
    builder->max_value_lens = carquet_mem_calloc(builder->capacity, sizeof(int32_t));
    builder->null_pages = carquet_mem_calloc(builder->capacity, sizeof(bool));

    if (!builder->null_counts || !builder->min_values || !builder->max_values ||
        !builder->min_value_lens || !builder->max_value_lens || !builder->null_pages) {
        carquet_column_index_builder_destroy(builder);
        return NULL;
    }

    return builder;
}

/**
 * Destroy a column index builder.
 */
void carquet_column_index_builder_destroy(carquet_column_index_builder_t* builder) {
    if (!builder) return;

    if (builder->min_values) {
        for (int32_t i = 0; i < builder->num_pages; i++) {
            carquet_mem_free(builder->min_values[i]);
        }
        carquet_mem_free(builder->min_values);
    }

    if (builder->max_values) {
        for (int32_t i = 0; i < builder->num_pages; i++) {
            carquet_mem_free(builder->max_values[i]);
        }
        carquet_mem_free(builder->max_values);
    }

    carquet_mem_free(builder->null_counts);
    carquet_mem_free(builder->min_value_lens);
    carquet_mem_free(builder->max_value_lens);
    carquet_mem_free(builder->null_pages);
    carquet_mem_free(builder->rep_level_histograms);
    carquet_mem_free(builder->def_level_histograms);
    carquet_mem_free(builder);
}

/**
 * Ensure capacity for more pages.
 */
/* Grow one builder array to `n` elements in place, returning OOM from the
 * enclosing function on failure with the array untouched. */
#define PAGE_INDEX_GROW(arr, n)                                              \
    do {                                                                     \
        void* grown_ = carquet_mem_realloc((arr), (size_t)(n) * sizeof(*(arr))); \
        if (!grown_) return CARQUET_ERROR_OUT_OF_MEMORY;                     \
        (arr) = grown_;                                                      \
    } while (0)

static carquet_status_t ensure_capacity(carquet_column_index_builder_t* builder) {
    if (builder->num_pages < builder->capacity) {
        return CARQUET_OK;
    }

    int32_t new_cap = builder->capacity * 2;

    /* Publish each array as soon as its realloc succeeds. If a later one
     * fails, every pointer in the builder is still live (and at least as large
     * as `capacity` says); holding the results back until all succeeded would
     * leave the builder pointing at blocks realloc had already freed. */
    PAGE_INDEX_GROW(builder->null_counts, new_cap);
    PAGE_INDEX_GROW(builder->min_values, new_cap);
    PAGE_INDEX_GROW(builder->min_value_lens, new_cap);
    PAGE_INDEX_GROW(builder->max_values, new_cap);
    PAGE_INDEX_GROW(builder->max_value_lens, new_cap);
    PAGE_INDEX_GROW(builder->null_pages, new_cap);

    /* Grow the flattened histogram arrays. The layout is page-major and
     * contiguous, so the existing num_pages*len prefix survives the realloc. */
    if (builder->track_histograms) {
        if (builder->rep_hist_len > 0) {
            int64_t* rh = carquet_mem_realloc(builder->rep_level_histograms,
                (size_t)new_cap * builder->rep_hist_len * sizeof(int64_t));
            if (!rh) return CARQUET_ERROR_OUT_OF_MEMORY;
            builder->rep_level_histograms = rh;
        }
        if (builder->def_hist_len > 0) {
            int64_t* dh = carquet_mem_realloc(builder->def_level_histograms,
                (size_t)new_cap * builder->def_hist_len * sizeof(int64_t));
            if (!dh) return CARQUET_ERROR_OUT_OF_MEMORY;
            builder->def_level_histograms = dh;
        }
    }

    builder->capacity = new_cap;

    /* Initialize new entries */
    for (int32_t i = builder->num_pages; i < new_cap; i++) {
        builder->null_counts[i] = 0;
        builder->min_values[i] = NULL;
        builder->min_value_lens[i] = 0;
        builder->max_values[i] = NULL;
        builder->max_value_lens[i] = 0;
        builder->null_pages[i] = false;
    }

    return CARQUET_OK;
}

/**
 * Add a page's statistics to the column index.
 */
carquet_status_t carquet_column_index_add_page(
    carquet_column_index_builder_t* builder,
    int64_t null_count,
    const void* min_value,
    int32_t min_value_len,
    const void* max_value,
    int32_t max_value_len,
    bool is_null_page,
    const int64_t* rep_level_hist,
    int32_t rep_level_hist_len,
    const int64_t* def_level_hist,
    int32_t def_level_hist_len) {

    if (!builder) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Histogram lengths are max_rep/def_level + 1 (Parquet nesting is shallow).
     * Guard against an out-of-range length so a caller bug can never drive a
     * runaway allocation below. */
    #define CARQUET_MAX_LEVEL_HIST_LEN 4096
    if (rep_level_hist_len < 0 || rep_level_hist_len > CARQUET_MAX_LEVEL_HIST_LEN) {
        rep_level_hist = NULL;
    }
    if (def_level_hist_len < 0 || def_level_hist_len > CARQUET_MAX_LEVEL_HIST_LEN) {
        def_level_hist = NULL;
    }
    #undef CARQUET_MAX_LEVEL_HIST_LEN

    /* Latch histogram tracking on the first page that supplies them. Lengths
     * are fixed for the whole column (derived from max rep/def levels). */
    if (!builder->track_histograms && (rep_level_hist || def_level_hist)) {
        builder->track_histograms = true;
        builder->rep_hist_len = rep_level_hist ? rep_level_hist_len : 0;
        builder->def_hist_len = def_level_hist ? def_level_hist_len : 0;
        if (builder->rep_hist_len > 0) {
            builder->rep_level_histograms = carquet_mem_calloc(
                (size_t)builder->capacity * builder->rep_hist_len, sizeof(int64_t));
            if (!builder->rep_level_histograms) return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (builder->def_hist_len > 0) {
            builder->def_level_histograms = carquet_mem_calloc(
                (size_t)builder->capacity * builder->def_hist_len, sizeof(int64_t));
            if (!builder->def_level_histograms) return CARQUET_ERROR_OUT_OF_MEMORY;
        }
    }

    carquet_status_t status = ensure_capacity(builder);
    if (status != CARQUET_OK) return status;

    int32_t idx = builder->num_pages;

    builder->null_counts[idx] = null_count;
    builder->null_pages[idx] = is_null_page;
    if (!is_null_page && (!min_value || !max_value)) {
        /* No statistics for a page holding values: see `discarded`. (A
         * zero-length bound with a non-NULL pointer is the empty string and
         * counts as present.) */
        builder->discarded = true;
    }

    if (builder->track_histograms) {
        if (builder->rep_hist_len > 0) {
            int64_t* dst = builder->rep_level_histograms +
                           (size_t)idx * builder->rep_hist_len;
            if (rep_level_hist && rep_level_hist_len == builder->rep_hist_len) {
                memcpy(dst, rep_level_hist,
                       (size_t)builder->rep_hist_len * sizeof(int64_t));
            } else {
                memset(dst, 0, (size_t)builder->rep_hist_len * sizeof(int64_t));
            }
        }
        if (builder->def_hist_len > 0) {
            int64_t* dst = builder->def_level_histograms +
                           (size_t)idx * builder->def_hist_len;
            if (def_level_hist && def_level_hist_len == builder->def_hist_len) {
                memcpy(dst, def_level_hist,
                       (size_t)builder->def_hist_len * sizeof(int64_t));
            } else {
                memset(dst, 0, (size_t)builder->def_hist_len * sizeof(int64_t));
            }
        }
    }

    /* Copy min value */
    if (min_value && min_value_len > 0) {
        builder->min_values[idx] = carquet_mem_malloc(min_value_len);
        if (!builder->min_values[idx]) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        memcpy(builder->min_values[idx], min_value, min_value_len);
        builder->min_value_lens[idx] = min_value_len;
    }

    /* Copy max value */
    if (max_value && max_value_len > 0) {
        builder->max_values[idx] = carquet_mem_malloc(max_value_len);
        if (!builder->max_values[idx]) {
            carquet_mem_free(builder->min_values[idx]);
            builder->min_values[idx] = NULL;
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        memcpy(builder->max_values[idx], max_value, max_value_len);
        builder->max_value_lens[idx] = max_value_len;
    }

    builder->num_pages++;
    return CARQUET_OK;
}

/* Withhold this chunk's ColumnIndex (see `discarded`). */
void carquet_column_index_builder_discard(carquet_column_index_builder_t* builder) {
    if (builder) builder->discarded = true;
}

/**
 * Set boundary order for the column index.
 */
void carquet_column_index_set_boundary_order(
    carquet_column_index_builder_t* builder,
    int32_t order) {
    if (builder) {
        builder->boundary_order = order;
    }
}

static bool index_logical_integer_is_unsigned(const carquet_logical_type_t* lt) {
    return lt &&
           lt->id == CARQUET_LOGICAL_INTEGER &&
           !lt->params.integer.is_signed;
}

static int compare_index_values(const carquet_column_index_builder_t* builder,
                                const uint8_t* a, int32_t alen,
                                const uint8_t* b, int32_t blen) {
    if (alen == blen) {
        switch (builder->type) {
            case CARQUET_PHYSICAL_FLOAT: {
                float av, bv;
                memcpy(&av, a, sizeof(av));
                memcpy(&bv, b, sizeof(bv));
                return (av < bv) ? -1 : (av > bv ? 1 : 0);
            }
            case CARQUET_PHYSICAL_DOUBLE: {
                double av, bv;
                memcpy(&av, a, sizeof(av));
                memcpy(&bv, b, sizeof(bv));
                return (av < bv) ? -1 : (av > bv ? 1 : 0);
            }
            case CARQUET_PHYSICAL_INT32:
                if (index_logical_integer_is_unsigned(&builder->logical_type)) {
                    uint32_t av, bv;
                    memcpy(&av, a, sizeof(av));
                    memcpy(&bv, b, sizeof(bv));
                    return (av < bv) ? -1 : (av > bv ? 1 : 0);
                } else {
                    int32_t av, bv;
                    memcpy(&av, a, sizeof(av));
                    memcpy(&bv, b, sizeof(bv));
                    return (av < bv) ? -1 : (av > bv ? 1 : 0);
                }
            case CARQUET_PHYSICAL_INT64:
                if (index_logical_integer_is_unsigned(&builder->logical_type)) {
                    uint64_t av, bv;
                    memcpy(&av, a, sizeof(av));
                    memcpy(&bv, b, sizeof(bv));
                    return (av < bv) ? -1 : (av > bv ? 1 : 0);
                } else {
                    int64_t av, bv;
                    memcpy(&av, a, sizeof(av));
                    memcpy(&bv, b, sizeof(bv));
                    return (av < bv) ? -1 : (av > bv ? 1 : 0);
                }
            default:
                break;
        }
    }

    int32_t n = alen < blen ? alen : blen;
    int c = n > 0 ? memcmp(a, b, (size_t)n) : 0;
    if (c != 0) return c;
    if (alen < blen) return -1;
    if (alen > blen) return 1;
    return 0;
}

/* ============================================================================
 * Offset Index Builder
 * ============================================================================
 */

struct carquet_offset_index_builder {
    int32_t capacity;
    int32_t num_pages;

    int64_t* offsets;
    int32_t* compressed_sizes;
    int64_t* first_row_indices;
    /* OffsetIndex field 2: unencoded_byte_array_data_bytes (Parquet 2.9),
     * list<i64>, one per page. Tracked only for BYTE_ARRAY columns. */
    int64_t* unencoded_bytes;
    bool track_unencoded;
};

/**
 * Create an offset index builder.
 */
carquet_offset_index_builder_t* carquet_offset_index_builder_create(
    bool track_unencoded) {

    carquet_offset_index_builder_t* builder = carquet_mem_calloc(1, sizeof(*builder));
    if (!builder) return NULL;

    builder->capacity = 16;
    builder->track_unencoded = track_unencoded;

    builder->offsets = carquet_mem_calloc(builder->capacity, sizeof(int64_t));
    builder->compressed_sizes = carquet_mem_calloc(builder->capacity, sizeof(int32_t));
    builder->first_row_indices = carquet_mem_calloc(builder->capacity, sizeof(int64_t));

    if (track_unencoded) {
        builder->unencoded_bytes = carquet_mem_calloc(builder->capacity, sizeof(int64_t));
    }

    if (!builder->offsets || !builder->compressed_sizes || !builder->first_row_indices ||
        (track_unencoded && !builder->unencoded_bytes)) {
        carquet_offset_index_builder_destroy(builder);
        return NULL;
    }

    return builder;
}

/**
 * Destroy an offset index builder.
 */
void carquet_offset_index_builder_destroy(carquet_offset_index_builder_t* builder) {
    if (!builder) return;

    carquet_mem_free(builder->offsets);
    carquet_mem_free(builder->compressed_sizes);
    carquet_mem_free(builder->first_row_indices);
    carquet_mem_free(builder->unencoded_bytes);
    carquet_mem_free(builder);
}

/**
 * Ensure capacity for more pages.
 */
static carquet_status_t offset_ensure_capacity(carquet_offset_index_builder_t* builder) {
    if (builder->num_pages < builder->capacity) {
        return CARQUET_OK;
    }

    int32_t new_cap = builder->capacity * 2;

    PAGE_INDEX_GROW(builder->offsets, new_cap);
    PAGE_INDEX_GROW(builder->compressed_sizes, new_cap);
    PAGE_INDEX_GROW(builder->first_row_indices, new_cap);

    if (builder->track_unencoded) {
        int64_t* new_unencoded = carquet_mem_realloc(builder->unencoded_bytes, new_cap * sizeof(int64_t));
        if (!new_unencoded) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        builder->unencoded_bytes = new_unencoded;
    }

    builder->capacity = new_cap;
    return CARQUET_OK;
}

/**
 * Shift every recorded page offset by `delta`. Used by the column writer
 * to convert per-column relative offsets (accumulated while values were
 * being flushed, before the column's absolute file offset was known) into
 * absolute file offsets at finalize time.
 */
void carquet_offset_index_builder_shift_offsets(
    carquet_offset_index_builder_t* builder, int64_t delta) {
    if (!builder || delta == 0) return;
    for (int32_t i = 0; i < builder->num_pages; i++) {
        builder->offsets[i] += delta;
    }
}

/**
 * Add a page's location to the offset index.
 */
carquet_status_t carquet_offset_index_add_page(
    carquet_offset_index_builder_t* builder,
    int64_t offset,
    int32_t compressed_size,
    int64_t first_row_index,
    int64_t unencoded_byte_array_bytes) {

    if (!builder) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = offset_ensure_capacity(builder);
    if (status != CARQUET_OK) return status;

    int32_t idx = builder->num_pages;

    builder->offsets[idx] = offset;
    builder->compressed_sizes[idx] = compressed_size;
    builder->first_row_indices[idx] = first_row_index;

    if (builder->track_unencoded) {
        builder->unencoded_bytes[idx] = unencoded_byte_array_bytes;
    }

    builder->num_pages++;
    return CARQUET_OK;
}

/* ============================================================================
 * Serialization to Thrift
 * ============================================================================
 */

/**
 * Serialize column index to buffer.
 */
carquet_status_t carquet_column_index_serialize(
    const carquet_column_index_builder_t* builder,
    carquet_buffer_t* output) {

    if (!builder || !output) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (builder->discarded) {
        return CARQUET_OK;  /* nothing serialized: no ColumnIndex for this chunk */
    }

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, output);

    thrift_write_struct_begin(&enc);

    /* Field 1: null_pages (list<bool>) */
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 1);
    thrift_write_list_begin(&enc, THRIFT_TYPE_TRUE, builder->num_pages);
    for (int32_t i = 0; i < builder->num_pages; i++) {
        thrift_write_bool(&enc, builder->null_pages[i]);
    }

    /* Field 2: min_values (list<binary>) */
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 2);
    thrift_write_list_begin(&enc, THRIFT_TYPE_BINARY, builder->num_pages);
    for (int32_t i = 0; i < builder->num_pages; i++) {
        if (builder->min_values[i]) {
            thrift_write_binary(&enc, builder->min_values[i], builder->min_value_lens[i]);
        } else {
            thrift_write_binary(&enc, NULL, 0);
        }
    }

    /* Field 3: max_values (list<binary>) */
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 3);
    thrift_write_list_begin(&enc, THRIFT_TYPE_BINARY, builder->num_pages);
    for (int32_t i = 0; i < builder->num_pages; i++) {
        if (builder->max_values[i]) {
            thrift_write_binary(&enc, builder->max_values[i], builder->max_value_lens[i]);
        } else {
            thrift_write_binary(&enc, NULL, 0);
        }
    }

    /* Field 4: boundary_order (i32) */
    thrift_write_field_header(&enc, THRIFT_TYPE_I32, 4);
    thrift_write_i32(&enc, builder->boundary_order);

    /* Field 5: null_counts (list<i64>) - optional */
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 5);
    thrift_write_list_begin(&enc, THRIFT_TYPE_I64, builder->num_pages);
    for (int32_t i = 0; i < builder->num_pages; i++) {
        thrift_write_i64(&enc, builder->null_counts[i]);
    }

    /* Field 6: repetition_level_histograms (list<i64>) - optional, Parquet 2.9.
     * Flattened page-major: for each page, (max_rep_level+1) buckets. Only
     * emitted for repeated columns (rep_hist_len > 1) since a flat column's
     * histogram is trivially [num_values] and carries no information. */
    if (builder->track_histograms && builder->rep_hist_len > 1 &&
        builder->rep_level_histograms) {
        int32_t total = builder->num_pages * builder->rep_hist_len;
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 6);
        thrift_write_list_begin(&enc, THRIFT_TYPE_I64, total);
        for (int32_t i = 0; i < total; i++) {
            thrift_write_i64(&enc, builder->rep_level_histograms[i]);
        }
    }

    /* Field 7: definition_level_histograms (list<i64>) - optional, Parquet 2.9.
     * Emitted when the column has definition levels (def_hist_len > 1), i.e.
     * it is nullable or nested; the histogram then encodes the null structure. */
    if (builder->track_histograms && builder->def_hist_len > 1 &&
        builder->def_level_histograms) {
        int32_t total = builder->num_pages * builder->def_hist_len;
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 7);
        thrift_write_list_begin(&enc, THRIFT_TYPE_I64, total);
        for (int32_t i = 0; i < total; i++) {
            thrift_write_i64(&enc, builder->def_level_histograms[i]);
        }
    }

    thrift_write_struct_end(&enc);
    return CARQUET_OK;
}

/**
 * Serialize offset index to buffer.
 */
carquet_status_t carquet_offset_index_serialize(
    const carquet_offset_index_builder_t* builder,
    carquet_buffer_t* output) {

    if (!builder || !output) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, output);

    thrift_write_struct_begin(&enc);

    /* Field 1: page_locations (list<PageLocation>) */
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 1);
    thrift_write_list_begin(&enc, THRIFT_TYPE_STRUCT, builder->num_pages);

    for (int32_t i = 0; i < builder->num_pages; i++) {
        thrift_write_struct_begin(&enc);

        /* PageLocation field 1: offset */
        thrift_write_field_header(&enc, THRIFT_TYPE_I64, 1);
        thrift_write_i64(&enc, builder->offsets[i]);

        /* PageLocation field 2: compressed_page_size */
        thrift_write_field_header(&enc, THRIFT_TYPE_I32, 2);
        thrift_write_i32(&enc, builder->compressed_sizes[i]);

        /* PageLocation field 3: first_row_index */
        thrift_write_field_header(&enc, THRIFT_TYPE_I64, 3);
        thrift_write_i64(&enc, builder->first_row_indices[i]);

        thrift_write_struct_end(&enc);
    }

    /* Field 2: unencoded_byte_array_data_bytes (list<i64>) - optional,
     * Parquet 2.9. Per-page total of BYTE_ARRAY value bytes assuming no
     * encoding (length prefixes excluded). Emitted only for BYTE_ARRAY. */
    if (builder->track_unencoded && builder->unencoded_bytes) {
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 2);
        thrift_write_list_begin(&enc, THRIFT_TYPE_I64, builder->num_pages);
        for (int32_t i = 0; i < builder->num_pages; i++) {
            thrift_write_i64(&enc, builder->unencoded_bytes[i]);
        }
    }

    thrift_write_struct_end(&enc);
    return CARQUET_OK;
}

/* ============================================================================
 * Page Filtering Using Column Index
 * ============================================================================
 */

/**
 * Check if a page might contain values in the given range.
 *
 * @param builder Column index builder
 * @param page_idx Page index
 * @param min_value Query min value (NULL for unbounded)
 * @param max_value Query max value (NULL for unbounded)
 * @param value_len Length of value for byte array types
 * @param might_match Output: true if page might contain matching values
 * @return Status code
 */
carquet_status_t carquet_column_index_page_might_match(
    const carquet_column_index_builder_t* builder,
    int32_t page_idx,
    const void* min_value,
    const void* max_value,
    int32_t value_len,
    bool* might_match) {

    if (!builder || !might_match || page_idx < 0 || page_idx >= builder->num_pages) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Null pages never match non-null predicates */
    if (builder->null_pages[page_idx]) {
        *might_match = false;
        return CARQUET_OK;
    }

    *might_match = true;  /* Assume match by default */

    /* If query max < page min, no match */
    if (max_value && builder->min_values[page_idx]) {
        int cmp = compare_index_values(
            builder, max_value, value_len,
            builder->min_values[page_idx], builder->min_value_lens[page_idx]);
        if (cmp < 0) {
            *might_match = false;
            return CARQUET_OK;
        }
    }

    /* If query min > page max, no match */
    if (min_value && builder->max_values[page_idx]) {
        int cmp = compare_index_values(
            builder, min_value, value_len,
            builder->max_values[page_idx], builder->max_value_lens[page_idx]);
        if (cmp > 0) {
            *might_match = false;
            return CARQUET_OK;
        }
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Deserialization from Thrift
 * ============================================================================
 */

/**
 * Free a parsed column index.
 */
void carquet_column_index_free(carquet_column_index_t* index) {
    if (!index) return;
    if (index->min_values) {
        for (int32_t i = 0; i < index->num_min_values; i++) carquet_mem_free(index->min_values[i]);
        carquet_mem_free(index->min_values);
    }
    if (index->max_values) {
        for (int32_t i = 0; i < index->num_max_values; i++) carquet_mem_free(index->max_values[i]);
        carquet_mem_free(index->max_values);
    }
    carquet_mem_free(index->min_value_lens);
    carquet_mem_free(index->max_value_lens);
    carquet_mem_free(index->null_counts);
    carquet_mem_free(index->null_pages);
    carquet_mem_free(index);
}

/**
 * Parse a Thrift-encoded ColumnIndex.
 *
 * Thrift schema:
 *   struct ColumnIndex {
 *     1: required list<bool> null_pages
 *     2: required list<binary> min_values
 *     3: required list<binary> max_values
 *     4: required BoundaryOrder boundary_order (i32 enum)
 *     5: optional list<i64> null_counts
 *   }
 *
 * @param data  Pointer to the Thrift-encoded data
 * @param size  Size of the data in bytes
 * @return Parsed column index, or NULL on failure. Caller must free with
 *         carquet_column_index_free().
 */
carquet_column_index_t* carquet_column_index_parse(const uint8_t* data, size_t size) {
    if (!data || size == 0) return NULL;

    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, size);

    struct carquet_column_index* ci = carquet_mem_calloc(1, sizeof(*ci));
    if (!ci) return NULL;

    thrift_read_struct_begin(&dec);

    thrift_type_t type;
    int16_t field_id;

    while (thrift_read_field_begin(&dec, &type, &field_id)) {
        switch (field_id) {
            case 1: { /* null_pages: list<bool> */
                if (ci->null_pages) { carquet_column_index_free(ci); return NULL; }
                thrift_type_t elem_type;
                int32_t count;
                thrift_read_list_begin(&dec, &elem_type, &count);
                if (count < 0 || count > 1000000) { carquet_column_index_free(ci); return NULL; }
                ci->num_pages = count;
                ci->null_pages = carquet_mem_calloc(count, sizeof(bool));
                if (!ci->null_pages) { carquet_column_index_free(ci); return NULL; }
                ci->num_null_pages = count;
                for (int32_t i = 0; i < count; i++) {
                    ci->null_pages[i] = thrift_read_bool(&dec);
                }
                break;
            }
            case 2: { /* min_values: list<binary> */
                if (ci->min_values) { carquet_column_index_free(ci); return NULL; }
                thrift_type_t elem_type;
                int32_t count;
                thrift_read_list_begin(&dec, &elem_type, &count);
                if (count < 0 || count > 1000000) { carquet_column_index_free(ci); return NULL; }
                ci->min_values = carquet_mem_calloc(count, sizeof(uint8_t*));
                ci->min_value_lens = carquet_mem_calloc(count, sizeof(int32_t));
                if (!ci->min_values || !ci->min_value_lens) {
                    carquet_column_index_free(ci);
                    return NULL;
                }
                ci->num_min_values = count;
                for (int32_t i = 0; i < count; i++) {
                    int32_t len;
                    const uint8_t* bin = thrift_read_binary(&dec, &len);
                    if (bin && len > 0) {
                        ci->min_values[i] = carquet_mem_malloc(len);
                        if (!ci->min_values[i]) {
                            carquet_column_index_free(ci);
                            return NULL;
                        }
                        memcpy(ci->min_values[i], bin, len);
                        ci->min_value_lens[i] = len;
                    }
                }
                break;
            }
            case 3: { /* max_values: list<binary> */
                if (ci->max_values) { carquet_column_index_free(ci); return NULL; }
                thrift_type_t elem_type;
                int32_t count;
                thrift_read_list_begin(&dec, &elem_type, &count);
                if (count < 0 || count > 1000000) { carquet_column_index_free(ci); return NULL; }
                ci->max_values = carquet_mem_calloc(count, sizeof(uint8_t*));
                ci->max_value_lens = carquet_mem_calloc(count, sizeof(int32_t));
                if (!ci->max_values || !ci->max_value_lens) {
                    carquet_column_index_free(ci);
                    return NULL;
                }
                ci->num_max_values = count;
                for (int32_t i = 0; i < count; i++) {
                    int32_t len;
                    const uint8_t* bin = thrift_read_binary(&dec, &len);
                    if (bin && len > 0) {
                        ci->max_values[i] = carquet_mem_malloc(len);
                        if (!ci->max_values[i]) {
                            carquet_column_index_free(ci);
                            return NULL;
                        }
                        memcpy(ci->max_values[i], bin, len);
                        ci->max_value_lens[i] = len;
                    }
                }
                break;
            }
            case 4: { /* boundary_order: i32 */
                ci->boundary_order = thrift_read_i32(&dec);
                break;
            }
            case 5: { /* null_counts: list<i64> */
                if (ci->null_counts) { carquet_column_index_free(ci); return NULL; }
                thrift_type_t elem_type;
                int32_t count;
                thrift_read_list_begin(&dec, &elem_type, &count);
                if (count < 0 || count > 1000000) { carquet_column_index_free(ci); return NULL; }
                ci->null_counts = carquet_mem_calloc(count, sizeof(int64_t));
                if (!ci->null_counts) { carquet_column_index_free(ci); return NULL; }
                ci->num_null_counts = count;
                for (int32_t i = 0; i < count; i++) {
                    ci->null_counts[i] = thrift_read_i64(&dec);
                }
                break;
            }
            default:
                thrift_skip_field(&dec, type);
                break;
        }
    }

    thrift_read_struct_end(&dec);

    if (thrift_decoder_has_error(&dec)) {
        carquet_column_index_free(ci);
        return NULL;
    }

    /* Clamp num_pages to the minimum of all parsed array sizes so that
     * accessors never read past any allocation — even with malformed
     * Thrift data where list counts disagree. */
    if (ci->null_pages && ci->num_null_pages < ci->num_pages)
        ci->num_pages = ci->num_null_pages;
    if (ci->min_values && ci->num_min_values < ci->num_pages)
        ci->num_pages = ci->num_min_values;
    if (ci->max_values && ci->num_max_values < ci->num_pages)
        ci->num_pages = ci->num_max_values;
    if (ci->null_counts && ci->num_null_counts < ci->num_pages)
        ci->num_pages = ci->num_null_counts;

    return ci;
}

/**
 * Parse a Thrift-encoded OffsetIndex.
 *
 * Thrift schema:
 *   struct OffsetIndex {
 *     1: required list<PageLocation> page_locations
 *   }
 *   struct PageLocation {
 *     1: required i64 offset
 *     2: required i32 compressed_page_size
 *     3: required i64 first_row_index
 *   }
 *
 * @param data  Pointer to the Thrift-encoded data
 * @param size  Size of the data in bytes
 * @return Parsed offset index, or NULL on failure. Caller must free with
 *         carquet_offset_index_free().
 */
carquet_offset_index_t* carquet_offset_index_parse(const uint8_t* data, size_t size) {
    if (!data || size == 0) return NULL;

    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, size);

    struct carquet_offset_index* oi = carquet_mem_calloc(1, sizeof(*oi));
    if (!oi) return NULL;

    thrift_read_struct_begin(&dec);

    thrift_type_t type;
    int16_t field_id;

    while (thrift_read_field_begin(&dec, &type, &field_id)) {
        switch (field_id) {
            case 1: { /* page_locations: list<PageLocation> */
                if (oi->page_locations) { carquet_offset_index_free(oi); return NULL; }
                thrift_type_t elem_type;
                int32_t count;
                thrift_read_list_begin(&dec, &elem_type, &count);
                if (count < 0 || count > 1000000) { carquet_mem_free(oi); return NULL; }
                oi->num_pages = count;
                oi->page_locations = carquet_mem_calloc(count, sizeof(carquet_page_location_t));
                if (!oi->page_locations) {
                    carquet_mem_free(oi);
                    return NULL;
                }
                for (int32_t i = 0; i < count; i++) {
                    thrift_read_struct_begin(&dec);

                    thrift_type_t ft;
                    int16_t fid;
                    while (thrift_read_field_begin(&dec, &ft, &fid)) {
                        switch (fid) {
                            case 1: /* offset: i64 */
                                oi->page_locations[i].offset = thrift_read_i64(&dec);
                                break;
                            case 2: /* compressed_page_size: i32 */
                                oi->page_locations[i].compressed_size = thrift_read_i32(&dec);
                                break;
                            case 3: /* first_row_index: i64 */
                                oi->page_locations[i].first_row_index = thrift_read_i64(&dec);
                                break;
                            default:
                                thrift_skip_field(&dec, ft);
                                break;
                        }
                    }

                    thrift_read_struct_end(&dec);
                }
                break;
            }
            default:
                thrift_skip_field(&dec, type);
                break;
        }
    }

    thrift_read_struct_end(&dec);

    if (thrift_decoder_has_error(&dec)) {
        carquet_offset_index_free(oi);
        return NULL;
    }

    return oi;
}

/* ============================================================================
 * Accessor Functions
 * ============================================================================
 */

/**
 * Get the number of pages in a column index.
 */
int32_t carquet_column_index_num_pages(const carquet_column_index_t* index) {
    /* index is nonnull per API contract */
    return index->num_pages;
}

/**
 * Get per-page statistics from a column index.
 */
carquet_status_t carquet_column_index_get_page_stats(
    const carquet_column_index_t* index,
    int32_t page_index,
    carquet_page_stats_t* stats) {

    /* index and stats are nonnull per API contract */
    if (page_index < 0 || page_index >= index->num_pages) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    stats->has_null_count = (index->null_counts && page_index < index->num_null_counts);
    stats->null_count = stats->has_null_count ? index->null_counts[page_index] : 0;
    stats->min_value = (index->min_values && page_index < index->num_min_values)
        ? index->min_values[page_index] : NULL;
    stats->min_value_size = (index->min_value_lens && page_index < index->num_min_values)
        ? index->min_value_lens[page_index] : 0;
    stats->max_value = (index->max_values && page_index < index->num_max_values)
        ? index->max_values[page_index] : NULL;
    stats->max_value_size = (index->max_value_lens && page_index < index->num_max_values)
        ? index->max_value_lens[page_index] : 0;
    stats->is_null_page = (index->null_pages && page_index < index->num_null_pages)
        ? index->null_pages[page_index] : false;
    return CARQUET_OK;
}

/**
 * Get boundary order of a column index.
 * @return 0=UNORDERED, 1=ASCENDING, 2=DESCENDING
 */
int32_t carquet_column_index_boundary_order(const carquet_column_index_t* index) {
    /* index is nonnull per API contract */
    return index->boundary_order;
}

/**
 * Get the number of pages in an offset index.
 */
int32_t carquet_offset_index_num_pages(const carquet_offset_index_t* index) {
    /* index is nonnull per API contract */
    return index->num_pages;
}

/**
 * Get page location from an offset index.
 */
carquet_status_t carquet_offset_index_get_page_location(
    const carquet_offset_index_t* index,
    int32_t page_index,
    carquet_page_location_t* location) {

    /* index and location are nonnull per API contract */
    if (page_index < 0 || page_index >= index->num_pages) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    *location = index->page_locations[page_index];
    return CARQUET_OK;
}

/**
 * Free a parsed offset index.
 */
void carquet_offset_index_free(carquet_offset_index_t* index) {
    if (!index) return;
    carquet_mem_free(index->page_locations);
    carquet_mem_free(index);
}
