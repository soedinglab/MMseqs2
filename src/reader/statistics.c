/**
 * @file statistics.c
 * @brief Row group statistics access and predicate pushdown
 *
 * Provides access to column statistics for intelligent row group filtering.
 * This enables predicate pushdown, allowing queries to skip entire row groups
 * that cannot contain matching data.
 */

#include <carquet/carquet.h>
#include "reader_internal.h"
#include "thrift/parquet_types.h"
#include "core/decimal_cmp.h"
#include "core/float16.h"
#include <string.h>

/* ============================================================================
 * Type-specific comparison
 * ============================================================================
 *
 * Statistics min/max are raw byte buffers from the file footer. They are not
 * guaranteed to be aligned, and their width must match the physical type
 * (1 byte for BOOLEAN, 2 for FLOAT16, 4/8 for ints/floats). All reads go
 * through memcpy to avoid unaligned access, and a width mismatch makes the
 * comparison "indeterminate" so the caller stays conservative.
 */

typedef enum {
    CMP_INT32_S,
    CMP_INT32_U,
    CMP_INT64_S,
    CMP_INT64_U,
    CMP_FLOAT,
    CMP_DOUBLE,
    CMP_BOOL,
    CMP_FLOAT16,
    CMP_DECIMAL_BE,  /* DECIMAL(FLBA/BYTE_ARRAY): signed big-endian two's comp. */
    CMP_BYTES,
    CMP_NONE         /* no defined order (INT96): never prune */
} cmp_kind_t;

#define CMP3(va, vb) (((va) > (vb)) - ((va) < (vb)))

/* Compare a predicate value against a statistics value of the given kind.
 * Returns true and stores the comparison (value <=> stat) in *out, or
 * returns false if the comparison cannot be performed safely (e.g. the
 * stat buffer is narrower than the type requires). */
static bool stat_compare(cmp_kind_t kind,
                         const void* val, size_t val_len,
                         const void* st, size_t st_len,
                         int* out) {
    switch (kind) {
        case CMP_BOOL: {
            if (val_len < 1 || st_len < 1) return false;
            uint8_t a = 0, b = 0;
            memcpy(&a, val, 1);
            memcpy(&b, st, 1);
            *out = CMP3(a ? 1 : 0, b ? 1 : 0);
            return true;
        }
        case CMP_INT32_S: {
            if (val_len < 4 || st_len < 4) return false;
            int32_t a, b;
            memcpy(&a, val, 4);
            memcpy(&b, st, 4);
            *out = CMP3(a, b);
            return true;
        }
        case CMP_INT32_U: {
            if (val_len < 4 || st_len < 4) return false;
            uint32_t a, b;
            memcpy(&a, val, 4);
            memcpy(&b, st, 4);
            *out = CMP3(a, b);
            return true;
        }
        case CMP_INT64_S: {
            if (val_len < 8 || st_len < 8) return false;
            int64_t a, b;
            memcpy(&a, val, 8);
            memcpy(&b, st, 8);
            *out = CMP3(a, b);
            return true;
        }
        case CMP_INT64_U: {
            if (val_len < 8 || st_len < 8) return false;
            uint64_t a, b;
            memcpy(&a, val, 8);
            memcpy(&b, st, 8);
            *out = CMP3(a, b);
            return true;
        }
        case CMP_FLOAT: {
            if (val_len < 4 || st_len < 4) return false;
            float a, b;
            memcpy(&a, val, 4);
            memcpy(&b, st, 4);
            if (a != a || b != b) return false;  /* NaN bound: indeterminate */
            *out = CMP3(a, b);
            return true;
        }
        case CMP_DOUBLE: {
            if (val_len < 8 || st_len < 8) return false;
            double a, b;
            memcpy(&a, val, 8);
            memcpy(&b, st, 8);
            if (a != a || b != b) return false;  /* NaN bound: indeterminate */
            *out = CMP3(a, b);
            return true;
        }
        case CMP_FLOAT16: {
            if (val_len < 2 || st_len < 2) return false;
            uint16_t ha, hb;
            memcpy(&ha, val, 2);
            memcpy(&hb, st, 2);
            float a = carquet_half_to_float(ha);
            float b = carquet_half_to_float(hb);
            if (a != a || b != b) return false;  /* NaN bound: indeterminate */
            *out = CMP3(a, b);
            return true;
        }
        case CMP_NONE:
            return false;
        case CMP_DECIMAL_BE: {
            *out = carquet_compare_decimal_be(val, val_len, st, st_len);
            return true;
        }
        case CMP_BYTES:
        default: {
            size_t min_len = val_len < st_len ? val_len : st_len;
            int cmp = (min_len > 0) ? memcmp(val, st, min_len) : 0;
            if (cmp != 0) {
                *out = cmp;
            } else {
                *out = (val_len > st_len) - (val_len < st_len);
            }
            return true;
        }
    }
}

/* Map a column's physical + logical type to a comparison kind that respects
 * signedness (UINT logical/converted types) and FLOAT16 numeric ordering. */
static cmp_kind_t get_cmp_kind(const parquet_schema_element_t* elem,
                               carquet_physical_type_t type) {
    bool is_unsigned = false;
    if (elem->has_logical_type &&
        elem->logical_type.id == CARQUET_LOGICAL_INTEGER) {
        is_unsigned = !elem->logical_type.params.integer.is_signed;
    } else if (elem->has_converted_type) {
        switch (elem->converted_type) {
            case CARQUET_CONVERTED_UINT_8:
            case CARQUET_CONVERTED_UINT_16:
            case CARQUET_CONVERTED_UINT_32:
            case CARQUET_CONVERTED_UINT_64:
                is_unsigned = true;
                break;
            default:
                break;
        }
    }

    /* DECIMAL backed by a byte physical type is ordered as signed big-endian
     * two's complement, not unsigned-lexicographic. INT32/INT64-backed DECIMAL
     * falls through to the native signed integer comparison. */
    bool is_decimal =
        (elem->has_logical_type && elem->logical_type.id == CARQUET_LOGICAL_DECIMAL) ||
        (elem->has_converted_type && elem->converted_type == CARQUET_CONVERTED_DECIMAL);

    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            return CMP_BOOL;
        case CARQUET_PHYSICAL_INT32:
            return is_unsigned ? CMP_INT32_U : CMP_INT32_S;
        case CARQUET_PHYSICAL_INT64:
            return is_unsigned ? CMP_INT64_U : CMP_INT64_S;
        case CARQUET_PHYSICAL_FLOAT:
            return CMP_FLOAT;
        case CARQUET_PHYSICAL_DOUBLE:
            return CMP_DOUBLE;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (elem->has_logical_type &&
                elem->logical_type.id == CARQUET_LOGICAL_FLOAT16) {
                return CMP_FLOAT16;
            }
            return is_decimal ? CMP_DECIMAL_BE : CMP_BYTES;
        case CARQUET_PHYSICAL_INT96:
            return CMP_NONE;  /* no defined sort order */
        default:
            return is_decimal ? CMP_DECIMAL_BE : CMP_BYTES;
    }
}

/* ============================================================================
 * Statistics Access
 * ============================================================================
 */

/* The deprecated `min` / `max` fields were computed with signed byte
 * comparison by parquet-mr before 1.10 (PARQUET-686), which is why they were
 * deprecated: for BYTE_ARRAY / FIXED_LEN_BYTE_ARRAY (UNSIGNED sort order) a
 * UTF-8 string starting with a byte >= 0x80 sorted *below* ASCII, so the
 * stored max can be smaller than real values and pruning on it drops rows.
 * Arrow / parquet-mr only trust the pair when the column's sort order is
 * SIGNED, where the two definitions agree. */
static bool deprecated_stats_usable(const carquet_reader_t* reader, int32_t column_index) {
    int32_t schema_idx = reader->schema->leaf_indices[column_index];
    const parquet_schema_element_t* elem = &reader->schema->elements[schema_idx];
    if (!elem->has_type) return false;
    switch (elem->type) {
        case CARQUET_PHYSICAL_BOOLEAN:
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE:
            return true;
        default:
            return false;  /* BYTE_ARRAY, FLBA, INT96 */
    }
}

carquet_status_t carquet_reader_column_statistics(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_column_statistics_t* stats) {

    /* reader and stats are nonnull per API contract */
    if (row_group_index < 0 || row_group_index >= reader->metadata.num_row_groups) {
        return CARQUET_ERROR_ROW_GROUP_NOT_FOUND;
    }

    if (column_index < 0 || column_index >= reader->schema->num_leaves) {
        return CARQUET_ERROR_COLUMN_NOT_FOUND;
    }

    memset(stats, 0, sizeof(*stats));

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];
    if (column_index >= rg->num_columns) {
        return CARQUET_ERROR_COLUMN_NOT_FOUND;
    }

    const parquet_column_chunk_t* chunk = &rg->columns[column_index];
    if (!chunk->has_metadata) {
        return CARQUET_OK;  /* No statistics available */
    }

    const parquet_column_metadata_t* meta = &chunk->metadata;
    stats->num_values = meta->num_values;

    if (!meta->has_statistics) {
        return CARQUET_OK;
    }

    const parquet_statistics_t* pstats = &meta->statistics;

    /* Null count */
    if (pstats->has_null_count) {
        stats->has_null_count = true;
        stats->null_count = pstats->null_count;
    }

    /* Distinct count */
    if (pstats->has_distinct_count) {
        stats->has_distinct_count = true;
        stats->distinct_count = pstats->distinct_count;
    }

    /* Min/max values - prefer new format, fall back to deprecated.
     * Presence is taken from the has_min_value/has_max_value flags rather than
     * from length, so a BYTE_ARRAY/STRING column whose true minimum is the
     * empty string still enables predicate pushdown instead of being treated
     * as having no stats. stat_compare handles a zero-length bound correctly. */
    if (pstats->has_min_value && pstats->has_max_value) {
        stats->has_min_max = true;
        stats->min_value = pstats->min_value;
        stats->min_value_size = pstats->min_value_len;
        stats->max_value = pstats->max_value;
        stats->max_value_size = pstats->max_value_len;
    } else if (pstats->min_deprecated && pstats->min_deprecated_len > 0 &&
               pstats->max_deprecated && pstats->max_deprecated_len > 0 &&
               deprecated_stats_usable(reader, column_index)) {
        stats->has_min_max = true;
        stats->min_value = pstats->min_deprecated;
        stats->min_value_size = pstats->min_deprecated_len;
        stats->max_value = pstats->max_deprecated;
        stats->max_value_size = pstats->max_deprecated_len;
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Predicate Pushdown
 * ============================================================================
 */

carquet_status_t carquet_reader_row_group_matches(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_compare_op_t op,
    const void* value,
    int32_t value_size,
    bool* might_match) {

    /* reader, value, might_match are nonnull per API contract */
    /* Default: might match (conservative) */
    *might_match = true;

    /* Get column statistics */
    carquet_column_statistics_t stats;
    carquet_status_t status = carquet_reader_column_statistics(
        reader, row_group_index, column_index, &stats);

    if (status != CARQUET_OK) {
        return status;
    }

    /* If no min/max stats, we can't filter */
    if (!stats.has_min_max) {
        return CARQUET_OK;
    }

    /* Get column type */
    int32_t schema_idx = reader->schema->leaf_indices[column_index];
    const parquet_schema_element_t* elem = &reader->schema->elements[schema_idx];
    carquet_physical_type_t type = elem->has_type ? elem->type : CARQUET_PHYSICAL_BYTE_ARRAY;

    cmp_kind_t kind = get_cmp_kind(elem, type);

    int cmp_min, cmp_max;
    if (!stat_compare(kind, value, (size_t)value_size,
                      stats.min_value, (size_t)stats.min_value_size, &cmp_min) ||
        !stat_compare(kind, value, (size_t)value_size,
                      stats.max_value, (size_t)stats.max_value_size, &cmp_max)) {
        /* Stats are not in the expected format for this type; cannot safely
         * prune. Stay conservative: the row group might match. */
        return CARQUET_OK;
    }

    /*
     * Determine if row group can be skipped based on comparison:
     *
     * For value comparison against [min, max] range:
     * - EQ: skip if value < min OR value > max
     * - NE: skip if min == max == value (all values are the same)
     * - LT: skip if min >= value (all values >= value)
     * - LE: skip if min > value
     * - GT: skip if max <= value
     * - GE: skip if max < value
     */

    switch (op) {
        case CARQUET_COMPARE_EQ:
            /* value == x: skip if value not in [min, max] */
            if (cmp_min < 0 || cmp_max > 0) {
                *might_match = false;
            }
            break;

        case CARQUET_COMPARE_NE:
            /* value != x: skip only if all values equal x */
            if (cmp_min == 0 && cmp_max == 0) {
                /* min == max == value, all values equal the search value */
                *might_match = false;
            }
            break;

        case CARQUET_COMPARE_LT:
            /* x < value: skip if min >= value */
            if (cmp_min <= 0) {
                *might_match = false;
            }
            break;

        case CARQUET_COMPARE_LE:
            /* x <= value: skip if min > value */
            if (cmp_min < 0) {
                *might_match = false;
            }
            break;

        case CARQUET_COMPARE_GT:
            /* x > value: skip if max <= value */
            if (cmp_max >= 0) {
                *might_match = false;
            }
            break;

        case CARQUET_COMPARE_GE:
            /* x >= value: skip if max < value */
            if (cmp_max > 0) {
                *might_match = false;
            }
            break;
    }

    return CARQUET_OK;
}

int32_t carquet_reader_filter_row_groups(
    const carquet_reader_t* reader,
    int32_t column_index,
    carquet_compare_op_t op,
    const void* value,
    int32_t value_size,
    int32_t* matching_indices,
    int32_t max_indices) {

    /* reader, value, matching_indices are nonnull per API contract */
    if (max_indices <= 0) {
        return -1;
    }

    int32_t num_row_groups = carquet_reader_num_row_groups(reader);
    int32_t num_matching = 0;

    for (int32_t i = 0; i < num_row_groups && num_matching < max_indices; i++) {
        bool might_match = true;

        carquet_status_t status = carquet_reader_row_group_matches(
            reader, i, column_index, op, value, value_size, &might_match);

        if (status != CARQUET_OK) {
            /* On error, include row group (conservative) */
            might_match = true;
        }

        if (might_match) {
            matching_indices[num_matching++] = i;
        }
    }

    return num_matching;
}
