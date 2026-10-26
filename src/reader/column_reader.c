/**
 * @file column_reader.c
 * @brief Column reading implementation
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include "reader_internal.h"
#include "thrift/parquet_types.h"
#include "encoding/plain.h"
#include "encoding/rle.h"
#include "core/endian.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Forward declaration */
extern carquet_status_t carquet_read_next_page(
    carquet_column_reader_t* reader,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels,
    int64_t* values_read,
    carquet_error_t* error);
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                 int16_t max_def_level);

/* ============================================================================
 * Batch Reading
 * ============================================================================
 */

static int64_t count_present_levels(
    const int16_t* def_levels,
    int64_t count,
    int16_t max_def_level) {
    return carquet_dispatch_count_non_nulls(def_levels, count, max_def_level);
}

int64_t carquet_column_read_batch_ex(
    carquet_column_reader_t* reader,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels,
    carquet_error_t* error) {

    /* Start from a clean slate so callers can rely on error->code == CARQUET_OK
     * meaning "no failure occurred", independent of the return value. */
    if (error) {
        carquet_error_clear(error);
    }

    /* max_values < 0 is invalid; max_values = 0 is a "peek" to trigger page loading */
    if (max_values < 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                          "max_values must be non-negative (got %lld)",
                          (long long)max_values);
        return -1;
    }
    if (max_values == 0) {
        /* Load page if needed, but don't read any values. A load failure here
         * is surfaced through `error` while preserving the historical return
         * value of 0 (no values read). */
        if (reader->values_remaining > 0 && !reader->page_loaded) {
            (void)carquet_column_ensure_page_loaded(reader, error);
        }
        return 0;
    }

    if (reader->values_remaining <= 0) {
        return 0;
    }

    carquet_error_t local_error = CARQUET_ERROR_INIT;
    int64_t total_read = 0;
    int64_t dense_values_read = 0;
    size_t value_size = 0;
    int16_t* scratch_def_levels = NULL;

    /* Determine value size for pointer arithmetic. Preserved dictionary pages
     * expose uint32_t indices instead of materialized physical values. */
    if (reader->preserve_dictionary) {
        value_size = sizeof(uint32_t);
    } else switch (reader->type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            value_size = 1;
            break;
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:
            value_size = 4;
            break;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE:
            value_size = 8;
            break;
        case CARQUET_PHYSICAL_INT96:
            value_size = 12;
            break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            value_size = reader->type_length;
            break;
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            /* Variable length - handled differently */
            value_size = sizeof(carquet_byte_array_t);
            break;
        default:
            CARQUET_SET_ERROR(error, CARQUET_ERROR_TYPE_MISMATCH,
                              "unknown physical type %d", (int)reader->type);
            return -1;
    }

    if (reader->max_def_level > 0 && !def_levels) {
        /* Nullable column read without caller-provided def_levels: decode the
         * levels into a per-reader scratch buffer that is grown on demand and
         * freed with the reader, instead of a malloc/free on every call. The
         * guard also rejects a byte count that would overflow size_t on a
         * 32-bit target ((size_t)max_values * 2 wraps for max_values >= 2^31). */
        if ((uint64_t)max_values > SIZE_MAX / sizeof(*scratch_def_levels)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                              "max_values %lld too large for scratch definition levels",
                              (long long)max_values);
            return -1;
        }
        if (reader->scratch_def_capacity < (size_t)max_values) {
            int16_t* grown = carquet_mem_realloc(
                reader->scratch_def_levels,
                (size_t)max_values * sizeof(*scratch_def_levels));
            if (!grown) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                                  "failed to allocate %lld scratch definition levels",
                                  (long long)max_values);
                return -1;
            }
            reader->scratch_def_levels = grown;
            reader->scratch_def_capacity = (size_t)max_values;
        }
        scratch_def_levels = reader->scratch_def_levels;
    }

    /* Read pages until we have enough values or run out */
    while (total_read < max_values && reader->values_remaining > 0) {
        int64_t values_read = 0;
        int64_t to_read = max_values - total_read;
        bool nullable = reader->max_def_level > 0;

        uint8_t* value_ptr = (uint8_t*)values +
            (size_t)(nullable ? dense_values_read : total_read) * value_size;
        int16_t* def_ptr = def_levels ? def_levels + total_read :
            (scratch_def_levels ? scratch_def_levels + total_read : NULL);
        int16_t* rep_ptr = rep_levels ? rep_levels + total_read : NULL;

        carquet_status_t status = carquet_read_next_page(
            reader, value_ptr, to_read, def_ptr, rep_ptr, &values_read, &local_error);

        if (status != CARQUET_OK) {
            /* Propagate the underlying failure verbatim. When some values were
             * already read we still return the salvaged count, but the error
             * out-parameter now lets the caller tell a truncated-by-error batch
             * apart from a clean end-of-column short read. */
            if (error) {
                carquet_error_copy(error, &local_error);
            }
            if (total_read > 0) {
                break;
            }
            return -1;
        }

        if (values_read == 0) {
            break;
        }

        if (nullable && def_ptr) {
            dense_values_read += count_present_levels(
                def_ptr, values_read, reader->max_def_level);
        } else {
            dense_values_read += values_read;
        }
        total_read += values_read;
    }

    /* scratch_def_levels (if used) is owned by the reader; nothing to free. */
    return total_read;
}

int64_t carquet_column_read_batch(
    carquet_column_reader_t* reader,
    void* values,
    int64_t max_values,
    int16_t* def_levels,
    int16_t* rep_levels) {
    /* Backward-compatible thin wrapper: identical behavior, error detail
     * discarded. New code that needs to distinguish failure modes should call
     * carquet_column_read_batch_ex() directly. */
    return carquet_column_read_batch_ex(
        reader, values, max_values, def_levels, rep_levels, NULL);
}

/* carquet_column_skip lives in page_reader.c: skipping is a page-state
 * operation that advances whole pages by parsing only their headers (no
 * decompression/decoding), so it belongs with the page-loading machinery. */
