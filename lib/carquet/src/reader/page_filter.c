/**
 * @file page_filter.c
 * @brief Page-level filter evaluation for the batch reader.
 *
 * Given a conjunction of clauses, builds the set of row ranges within a
 * row group that may contain matching values. Pages whose [min, max] stats
 * cannot overlap the predicate are pruned; the survivors are mapped to row
 * ranges via the offset index, and per-clause range lists are intersected
 * with a two-finger sweep.
 *
 * The evaluator handles every Parquet physical type whose sort order is
 * defined (BOOLEAN, INT32, INT64, FLOAT, DOUBLE, BYTE_ARRAY,
 * FIXED_LEN_BYTE_ARRAY), with logical-type adjustments for unsigned
 * integers (UINT8/16/32/64) and IEEE half-precision (FLOAT16). NaN
 * scalars in predicate values match nothing under ordered operators
 * (Arrow semantics).
 *
 * Bounds in the column index are conservative by construction — the
 * writer rounds truncated BYTE_ARRAY min downward (lex) and truncated
 * max upward — so we never need to treat them as exact.
 */

#include "page_filter.h"
#include "core/allocator.h"
#include "core/decimal_cmp.h"
#include "reader_internal.h"
#include "thrift/parquet_types.h"
#include <carquet/error.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Row range list
 * ============================================================================ */

void carquet_row_range_list_init(carquet_row_range_list_t* list) {
    list->ranges = NULL;
    list->count = 0;
    list->capacity = 0;
    list->total_rows = 0;
}

void carquet_row_range_list_destroy(carquet_row_range_list_t* list) {
    if (!list) return;
    carquet_mem_free(list->ranges);
    list->ranges = NULL;
    list->count = 0;
    list->capacity = 0;
    list->total_rows = 0;
}

void carquet_row_range_list_clear(carquet_row_range_list_t* list) {
    list->count = 0;
    list->total_rows = 0;
}

carquet_status_t carquet_row_range_list_append(
    carquet_row_range_list_t* list, int64_t first_row, int64_t num_rows) {

    if (num_rows <= 0) return CARQUET_OK;

    /* Coalesce with the last range if abutting. */
    if (list->count > 0) {
        carquet_row_range_t* last = &list->ranges[list->count - 1];
        if (last->first_row + last->num_rows == first_row) {
            last->num_rows += num_rows;
            list->total_rows += num_rows;
            return CARQUET_OK;
        }
    }

    if (list->count >= list->capacity) {
        int32_t new_cap = list->capacity > 0 ? list->capacity * 2 : 8;
        carquet_row_range_t* nr = carquet_mem_realloc(
            list->ranges, (size_t)new_cap * sizeof(carquet_row_range_t));
        if (!nr) return CARQUET_ERROR_OUT_OF_MEMORY;
        list->ranges = nr;
        list->capacity = new_cap;
    }

    list->ranges[list->count].first_row = first_row;
    list->ranges[list->count].num_rows = num_rows;
    list->count++;
    list->total_rows += num_rows;
    return CARQUET_OK;
}

/* ============================================================================
 * Schema lookup
 * ============================================================================ */

typedef struct {
    carquet_physical_type_t physical_type;
    int32_t type_length;
    bool is_unsigned;
    bool is_float16;
    bool is_decimal;  /* DECIMAL backed by FLBA/BYTE_ARRAY: signed big-endian order */
} column_type_info_t;

static carquet_status_t lookup_column_type(
    carquet_reader_t* file_reader, int32_t column_index,
    column_type_info_t* out, carquet_error_t* error) {

    const carquet_schema_t* schema = carquet_reader_schema(file_reader);
    if (!schema) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "No schema");
        return CARQUET_ERROR_INVALID_STATE;
    }
    if (column_index < 0 || column_index >= schema->num_leaves) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Filter column index %d out of range [0, %d)",
            column_index, schema->num_leaves);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int32_t schema_idx = schema->leaf_indices[column_index];
    const parquet_schema_element_t* elem = &schema->elements[schema_idx];

    out->physical_type = elem->has_type ? elem->type : CARQUET_PHYSICAL_BYTE_ARRAY;
    out->type_length = elem->type_length;
    out->is_unsigned = false;
    out->is_float16 = false;
    out->is_decimal = false;

    if (elem->has_logical_type) {
        if (elem->logical_type.id == CARQUET_LOGICAL_INTEGER &&
            !elem->logical_type.params.integer.is_signed) {
            out->is_unsigned = true;
        } else if (elem->logical_type.id == CARQUET_LOGICAL_FLOAT16) {
            out->is_float16 = true;
        } else if (elem->logical_type.id == CARQUET_LOGICAL_DECIMAL &&
                   (out->physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY ||
                    out->physical_type == CARQUET_PHYSICAL_BYTE_ARRAY)) {
            /* INT32/INT64-backed DECIMAL is ordered by the native signed
             * integer path; only the byte backings need signed big-endian. */
            out->is_decimal = true;
        }
    }
    if (elem->has_converted_type) {
        switch (elem->converted_type) {
            case CARQUET_CONVERTED_UINT_8:
            case CARQUET_CONVERTED_UINT_16:
            case CARQUET_CONVERTED_UINT_32:
            case CARQUET_CONVERTED_UINT_64:
                out->is_unsigned = true;
                break;
            case CARQUET_CONVERTED_DECIMAL:
                if (out->physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY ||
                    out->physical_type == CARQUET_PHYSICAL_BYTE_ARRAY) {
                    out->is_decimal = true;
                }
                break;
            default:
                break;
        }
    }

    return CARQUET_OK;
}

static size_t scalar_size(carquet_physical_type_t pt, int32_t type_length) {
    switch (pt) {
        case CARQUET_PHYSICAL_BOOLEAN: return 1;
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT: return 4;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE: return 8;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return (type_length > 0) ? (size_t)type_length : 0;
        default: return 0;
    }
}

/**
 * Whether a ColumnIndex min/max byte length is usable for this column type.
 *
 * compare_typed() loads a fixed native width for fixed-width physical types,
 * so a stat shorter (or longer) than that width must be rejected to avoid an
 * out-of-bounds read on a malformed file. BYTE_ARRAY is variable-width and
 * compared lexicographically, so any positive length is fine.
 */
static bool stat_len_ok(const column_type_info_t* ti, int32_t len) {
    if (len <= 0) return false;
    size_t expected = scalar_size(ti->physical_type, ti->type_length);
    if (expected == 0) {
        /* Variable-width (BYTE_ARRAY) — length-safe in compare_typed(). */
        return true;
    }
    return (size_t)len == expected;
}

/* ============================================================================
 * Typed comparator
 * ============================================================================
 *
 * Compares two values of the same column type. Returns < 0, 0, > 0.
 * For FLOAT/DOUBLE: NaN-bearing inputs are not expected here — predicate
 * scalars containing NaN are handled by the caller (matching nothing under
 * ordered ops); page min/max are never NaN per Parquet semantics.
 */

static float decode_float16(const uint8_t* b) {
    uint16_t raw = (uint16_t)b[0] | ((uint16_t)b[1] << 8);
    uint16_t sign = (raw >> 15) & 0x1;
    uint16_t exp = (raw >> 10) & 0x1F;
    uint16_t mant = raw & 0x3FF;
    uint32_t f;
    if (exp == 0) {
        if (mant == 0) {
            f = (uint32_t)sign << 31;
        } else {
            /* Subnormal: normalize. */
            int32_t e = -14;
            while ((mant & 0x400) == 0) { mant <<= 1; e--; }
            mant &= 0x3FF;
            f = ((uint32_t)sign << 31) |
                ((uint32_t)(e + 127) << 23) |
                ((uint32_t)mant << 13);
        }
    } else if (exp == 0x1F) {
        f = ((uint32_t)sign << 31) | (0xFFu << 23) | ((uint32_t)mant << 13);
    } else {
        f = ((uint32_t)sign << 31) |
            ((uint32_t)(exp - 15 + 127) << 23) |
            ((uint32_t)mant << 13);
    }
    float result;
    memcpy(&result, &f, sizeof(result));
    return result;
}

static int cmp_bytes_lex(const uint8_t* a, int32_t alen,
                         const uint8_t* b, int32_t blen) {
    int32_t n = alen < blen ? alen : blen;
    /* An empty-string predicate value (vlen == 0) reaches here with n == 0;
     * memcmp requires valid pointers even for a zero length, so skip the call
     * (same guard as page_index.c / reader/statistics.c). */
    int c = n > 0 ? memcmp(a, b, (size_t)n) : 0;
    if (c != 0) return c < 0 ? -1 : 1;
    if (alen == blen) return 0;
    return alen < blen ? -1 : 1;
}

static int compare_typed(const column_type_info_t* ti,
                         const uint8_t* a, int32_t alen,
                         const uint8_t* b, int32_t blen) {
    switch (ti->physical_type) {
        case CARQUET_PHYSICAL_BOOLEAN: {
            uint8_t av = a[0] ? 1 : 0;
            uint8_t bv = b[0] ? 1 : 0;
            return (av < bv) ? -1 : (av > bv ? 1 : 0);
        }
        case CARQUET_PHYSICAL_INT32: {
            if (ti->is_unsigned) {
                uint32_t av, bv;
                memcpy(&av, a, 4); memcpy(&bv, b, 4);
                return (av < bv) ? -1 : (av > bv ? 1 : 0);
            }
            int32_t av, bv;
            memcpy(&av, a, 4); memcpy(&bv, b, 4);
            return (av < bv) ? -1 : (av > bv ? 1 : 0);
        }
        case CARQUET_PHYSICAL_INT64: {
            if (ti->is_unsigned) {
                uint64_t av, bv;
                memcpy(&av, a, 8); memcpy(&bv, b, 8);
                return (av < bv) ? -1 : (av > bv ? 1 : 0);
            }
            int64_t av, bv;
            memcpy(&av, a, 8); memcpy(&bv, b, 8);
            return (av < bv) ? -1 : (av > bv ? 1 : 0);
        }
        case CARQUET_PHYSICAL_FLOAT: {
            float av, bv;
            memcpy(&av, a, 4); memcpy(&bv, b, 4);
            return (av < bv) ? -1 : (av > bv ? 1 : 0);
        }
        case CARQUET_PHYSICAL_DOUBLE: {
            double av, bv;
            memcpy(&av, a, 8); memcpy(&bv, b, 8);
            return (av < bv) ? -1 : (av > bv ? 1 : 0);
        }
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (ti->is_float16 && alen == 2 && blen == 2) {
                float av = decode_float16(a);
                float bv = decode_float16(b);
                return (av < bv) ? -1 : (av > bv ? 1 : 0);
            }
            if (ti->is_decimal) {
                return carquet_compare_decimal_be(a, (size_t)alen, b, (size_t)blen);
            }
            return cmp_bytes_lex(a, alen, b, blen);
        case CARQUET_PHYSICAL_BYTE_ARRAY:
            if (ti->is_decimal) {
                return carquet_compare_decimal_be(a, (size_t)alen, b, (size_t)blen);
            }
            return cmp_bytes_lex(a, alen, b, blen);
        default:
            return 0;
    }
}

static bool predicate_value_is_nan(const column_type_info_t* ti, const uint8_t* v) {
    if (!v) return false;
    if (ti->physical_type == CARQUET_PHYSICAL_FLOAT) {
        float f;
        memcpy(&f, v, 4);
        return isnan(f);
    }
    if (ti->physical_type == CARQUET_PHYSICAL_DOUBLE) {
        double d;
        memcpy(&d, v, 8);
        return isnan(d);
    }
    if (ti->is_float16 && ti->physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
        float f = decode_float16(v);
        return isnan(f);
    }
    return false;
}

/* ============================================================================
 * Clause validation
 * ============================================================================ */

static carquet_status_t validate_scalar_size(
    const column_type_info_t* ti, int32_t got_size,
    carquet_error_t* error) {

    size_t expected = scalar_size(ti->physical_type, ti->type_length);
    if (ti->physical_type == CARQUET_PHYSICAL_BYTE_ARRAY) {
        if (got_size < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Negative BYTE_ARRAY size in filter clause");
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        return CARQUET_OK;
    }
    if (expected == 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Unsupported physical type %d for filter",
            (int)ti->physical_type);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    /* For numeric types, value_size is informational; we use scalar_size. */
    if (ti->physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
        if (got_size != ti->type_length) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "FIXED_LEN_BYTE_ARRAY value_size %d != type_length %d",
                got_size, ti->type_length);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
    }
    return CARQUET_OK;
}

carquet_status_t carquet_page_filter_validate_clause(
    carquet_reader_t* file_reader,
    const carquet_filter_clause_t* clause,
    carquet_error_t* error) {

    column_type_info_t ti;
    carquet_status_t st = lookup_column_type(file_reader,
        clause->column_index, &ti, error);
    if (st != CARQUET_OK) return st;

    if (ti.physical_type == CARQUET_PHYSICAL_INT96) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "INT96 has no defined sort order; filter not supported");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    switch (clause->op) {
        case CARQUET_FILTER_EQ:
        case CARQUET_FILTER_NE:
        case CARQUET_FILTER_LT:
        case CARQUET_FILTER_LE:
        case CARQUET_FILTER_GT:
        case CARQUET_FILTER_GE:
            if (!clause->value) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Clause %d on column %d has NULL value",
                    (int)clause->op, clause->column_index);
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            return validate_scalar_size(&ti, clause->value_size, error);

        case CARQUET_FILTER_RANGE:
            if (!clause->has_lo && !clause->has_hi) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "RANGE clause has neither lower nor upper bound");
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            if (clause->has_lo) {
                if (!clause->lo) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                        "RANGE clause has_lo set but lo is NULL");
                    return CARQUET_ERROR_INVALID_ARGUMENT;
                }
                st = validate_scalar_size(&ti, clause->lo_size, error);
                if (st != CARQUET_OK) return st;
            }
            if (clause->has_hi) {
                if (!clause->hi) {
                    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                        "RANGE clause has_hi set but hi is NULL");
                    return CARQUET_ERROR_INVALID_ARGUMENT;
                }
                st = validate_scalar_size(&ti, clause->hi_size, error);
                if (st != CARQUET_OK) return st;
            }
            return CARQUET_OK;

        case CARQUET_FILTER_IN:
            if (!clause->values || clause->value_count <= 0) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "IN clause has empty value set");
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            if (clause->value_count > 256) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "IN clause value count %d exceeds 256",
                    clause->value_count);
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            return CARQUET_OK;

        case CARQUET_FILTER_IS_NULL:
        case CARQUET_FILTER_IS_NOT_NULL:
            return CARQUET_OK;
    }

    CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
        "Unknown filter op %d", (int)clause->op);
    return CARQUET_ERROR_INVALID_ARGUMENT;
}

/* ============================================================================
 * Per-page predicate evaluation
 * ============================================================================ */

/**
 * Extract the i-th IN value as (ptr, len).
 * For BYTE_ARRAY, values[] is carquet_byte_array_t entries.
 * For FIXED_LEN_BYTE_ARRAY, values[] is packed raw bytes of stride
 *  type_length.
 * For numeric, values[] is packed scalars of native stride.
 */
static void in_value_at(const column_type_info_t* ti,
                        const carquet_filter_clause_t* clause,
                        int32_t i,
                        const uint8_t** out_ptr, int32_t* out_len) {
    if (ti->physical_type == CARQUET_PHYSICAL_BYTE_ARRAY) {
        const carquet_byte_array_t* arr =
            (const carquet_byte_array_t*)clause->values;
        *out_ptr = arr[i].data;
        *out_len = arr[i].length;
        return;
    }
    if (ti->physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
        const uint8_t* base = (const uint8_t*)clause->values;
        *out_ptr = base + (size_t)i * (size_t)ti->type_length;
        *out_len = ti->type_length;
        return;
    }
    size_t sz = scalar_size(ti->physical_type, ti->type_length);
    const uint8_t* base = (const uint8_t*)clause->values;
    *out_ptr = base + (size_t)i * sz;
    *out_len = (int32_t)sz;
}

/**
 * Decide whether to keep this page given one clause.
 *
 * `keep = true` means: "the page may contain a row that satisfies the
 * clause." We err conservatively: a page we cannot prove empty is kept.
 */
static bool page_matches_clause(const column_type_info_t* ti,
                                const carquet_page_stats_t* stats,
                                const carquet_filter_clause_t* clause,
                                int64_t page_first_row,
                                int64_t page_num_rows) {
    (void)page_first_row;

    /* Null-presence ops: decide purely from the null count + null page
     * flag. */
    if (clause->op == CARQUET_FILTER_IS_NULL) {
        /* null_counts is optional in the ColumnIndex: without it a page
         * cannot be proven null-free. */
        if (!stats->has_null_count) return true;
        return stats->null_count > 0 || stats->is_null_page;
    }
    if (clause->op == CARQUET_FILTER_IS_NOT_NULL) {
        /* Page is all-nulls iff is_null_page is true. */
        if (stats->is_null_page) return false;
        /* If we can compare null_count to a known row count, do so for a
         * tighter decision; otherwise keep. */
        if (page_num_rows > 0 && stats->null_count >= page_num_rows) {
            return false;
        }
        return true;
    }

    /* Non-null ops: an all-null page can't match. */
    if (stats->is_null_page) return false;

    /* Page-num-rows known and null_count equals it ⇒ effectively all-null. */
    if (page_num_rows > 0 && stats->null_count >= page_num_rows) {
        return false;
    }

    /* Missing min/max ⇒ we have no bound information; keep conservatively.
     *
     * The min/max byte lengths come straight from the file's ColumnIndex and
     * are NOT guaranteed to match the column's physical width. compare_typed()
     * loads a fixed native width for numeric/BOOLEAN/FLOAT16 columns, so a
     * malformed (short) stat would read out of bounds. Treat any stat whose
     * length does not fit the column type as absent: we then keep the page
     * conservatively rather than prune on an unreliable bound. Variable-width
     * BYTE_ARRAY stats (and the lexicographic FLBA path) are length-safe in
     * compare_typed(), so they only require a positive length. */
    const uint8_t* pmin = (const uint8_t*)stats->min_value;
    const uint8_t* pmax = (const uint8_t*)stats->max_value;
    int32_t pmin_len = stats->min_value_size;
    int32_t pmax_len = stats->max_value_size;
    bool have_min = (pmin != NULL && stat_len_ok(ti, pmin_len));
    bool have_max = (pmax != NULL && stat_len_ok(ti, pmax_len));
    /* A NaN bound (written by older parquet-cpp and some other writers; the
     * spec only says writers *should* skip NaN) orders against nothing, so
     * it is no bound at all. */
    if (have_min && predicate_value_is_nan(ti, pmin)) have_min = false;
    if (have_max && predicate_value_is_nan(ti, pmax)) have_max = false;

    if (!have_min && !have_max) {
        /* No stats available — cannot prune. */
        return true;
    }

    switch (clause->op) {
        case CARQUET_FILTER_EQ: {
            const uint8_t* v = (const uint8_t*)clause->value;
            int32_t vlen = clause->value_size;
            if (predicate_value_is_nan(ti, v)) return false;
            if (have_min && compare_typed(ti, v, vlen, pmin, pmin_len) < 0) {
                return false;
            }
            if (have_max && compare_typed(ti, v, vlen, pmax, pmax_len) > 0) {
                return false;
            }
            return true;
        }
        case CARQUET_FILTER_NE: {
            /* Reject only if we can prove every value equals v.
             * That requires exact, equal min == max == v. With possibly
             * truncated bounds, equality is rare for BYTE_ARRAY but
             * always provable for numeric/FLBA-with-fixed-length. */
            if (have_min && have_max) {
                const uint8_t* v = (const uint8_t*)clause->value;
                int32_t vlen = clause->value_size;
                if (predicate_value_is_nan(ti, v)) return true;
                if (compare_typed(ti, pmin, pmin_len, pmax, pmax_len) == 0 &&
                    compare_typed(ti, v, vlen, pmin, pmin_len) == 0) {
                    return false;
                }
            }
            return true;
        }
        case CARQUET_FILTER_LT: {
            const uint8_t* v = (const uint8_t*)clause->value;
            int32_t vlen = clause->value_size;
            if (predicate_value_is_nan(ti, v)) return false;
            /* Keep iff min < v. */
            if (have_min && compare_typed(ti, pmin, pmin_len, v, vlen) >= 0) {
                return false;
            }
            return true;
        }
        case CARQUET_FILTER_LE: {
            const uint8_t* v = (const uint8_t*)clause->value;
            int32_t vlen = clause->value_size;
            if (predicate_value_is_nan(ti, v)) return false;
            if (have_min && compare_typed(ti, pmin, pmin_len, v, vlen) > 0) {
                return false;
            }
            return true;
        }
        case CARQUET_FILTER_GT: {
            const uint8_t* v = (const uint8_t*)clause->value;
            int32_t vlen = clause->value_size;
            if (predicate_value_is_nan(ti, v)) return false;
            if (have_max && compare_typed(ti, pmax, pmax_len, v, vlen) <= 0) {
                return false;
            }
            return true;
        }
        case CARQUET_FILTER_GE: {
            const uint8_t* v = (const uint8_t*)clause->value;
            int32_t vlen = clause->value_size;
            if (predicate_value_is_nan(ti, v)) return false;
            if (have_max && compare_typed(ti, pmax, pmax_len, v, vlen) < 0) {
                return false;
            }
            return true;
        }
        case CARQUET_FILTER_RANGE: {
            if (clause->has_lo) {
                const uint8_t* lo = (const uint8_t*)clause->lo;
                int32_t lo_len = clause->lo_size;
                if (predicate_value_is_nan(ti, lo)) return false;
                /* Need page max >= lo. */
                if (have_max && compare_typed(ti, pmax, pmax_len, lo, lo_len) < 0) {
                    return false;
                }
            }
            if (clause->has_hi) {
                const uint8_t* hi = (const uint8_t*)clause->hi;
                int32_t hi_len = clause->hi_size;
                if (predicate_value_is_nan(ti, hi)) return false;
                /* Need page min <= hi. */
                if (have_min && compare_typed(ti, pmin, pmin_len, hi, hi_len) > 0) {
                    return false;
                }
            }
            return true;
        }
        case CARQUET_FILTER_IN: {
            for (int32_t i = 0; i < clause->value_count; i++) {
                const uint8_t* v = NULL;
                int32_t vlen = 0;
                in_value_at(ti, clause, i, &v, &vlen);
                if (!v || predicate_value_is_nan(ti, v)) continue;
                bool below = have_min &&
                    compare_typed(ti, v, vlen, pmin, pmin_len) < 0;
                bool above = have_max &&
                    compare_typed(ti, v, vlen, pmax, pmax_len) > 0;
                if (!below && !above) return true;
            }
            return false;
        }
        case CARQUET_FILTER_IS_NULL:
        case CARQUET_FILTER_IS_NOT_NULL:
            /* Handled above. */
            return true;
    }
    return true;
}

/* ============================================================================
 * Row-group-level pruning (statistics + bloom filter)
 *
 * Before touching the page index, cheaply try to drop a whole row group
 * using the ColumnChunk-level statistics and, for equality-style clauses,
 * the column's bloom filter. This is purely additive to the page-level
 * range derivation below: it can only turn a provably-empty row group into
 * an empty range list (which the caller intersects to "skip"); it never
 * rescues a row group that the page path would otherwise reject, and it
 * works even for files that carry no page index. Every decision is
 * conservative — anything not provably empty is kept.
 * ============================================================================ */

/**
 * Decide whether a row group may match one clause given its ColumnChunk
 * statistics. Reuses page_matches_clause() with a synthesized page-stats
 * view over the whole row group. Returns true when the group may match
 * (including when the stats are missing or untrustworthy).
 */
static bool rg_stats_might_match(const column_type_info_t* ti,
                                 const carquet_filter_clause_t* clause,
                                 const carquet_column_statistics_t* s,
                                 int64_t rg_rows) {
    bool null_op = (clause->op == CARQUET_FILTER_IS_NULL ||
                    clause->op == CARQUET_FILTER_IS_NOT_NULL);
    /* Null-presence pruning is only sound when we actually know the null
     * count; a missing null_count must not be read as "zero nulls". */
    if (null_op && !s->has_null_count) return true;

    carquet_page_stats_t ps;
    ps.null_count = s->has_null_count ? s->null_count : 0;
    ps.min_value = s->has_min_max ? s->min_value : NULL;
    ps.min_value_size = s->has_min_max ? s->min_value_size : 0;
    ps.max_value = s->has_min_max ? s->max_value : NULL;
    ps.max_value_size = s->has_min_max ? s->max_value_size : 0;
    ps.is_null_page = false;

    /* page_matches_clause() uses page_num_rows only for its all-null
     * short-circuit (null_count >= page_num_rows). That is only meaningful
     * when the null count is trustworthy; pass 0 otherwise to disable it. */
    int64_t nrows = s->has_null_count ? rg_rows : 0;
    return page_matches_clause(ti, &ps, clause, 0, nrows);
}

/** Query one value against a bloom filter, dispatching on physical type. */
static bool bloom_check_value(const carquet_bloom_filter_t* bf,
                              const column_type_info_t* ti,
                              const void* v, int32_t vlen) {
    switch (ti->physical_type) {
        case CARQUET_PHYSICAL_INT32: {
            int32_t x; memcpy(&x, v, sizeof(x));
            return carquet_bloom_filter_check_i32(bf, x);
        }
        case CARQUET_PHYSICAL_INT64: {
            int64_t x; memcpy(&x, v, sizeof(x));
            return carquet_bloom_filter_check_i64(bf, x);
        }
        case CARQUET_PHYSICAL_FLOAT: {
            float x; memcpy(&x, v, sizeof(x));
            return carquet_bloom_filter_check_float(bf, x);
        }
        case CARQUET_PHYSICAL_DOUBLE: {
            double x; memcpy(&x, v, sizeof(x));
            return carquet_bloom_filter_check_double(bf, x);
        }
        case CARQUET_PHYSICAL_BYTE_ARRAY:
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return carquet_bloom_filter_check_bytes(bf, (const uint8_t*)v,
                                                    (size_t)vlen);
        default:
            /* INT96 / BOOLEAN: no bloom pruning. */
            return true;
    }
}

/**
 * For equality-style clauses (EQ / IN), consult the column's bloom filter.
 * Returns false only when the filter proves every candidate value absent.
 * A missing filter, an unsupported type, or any read error yields true.
 */
static bool rg_bloom_might_contain(carquet_reader_t* file_reader,
                                   int32_t row_group_index,
                                   const carquet_filter_clause_t* clause,
                                   const column_type_info_t* ti) {
    if (clause->op != CARQUET_FILTER_EQ && clause->op != CARQUET_FILTER_IN) {
        return true;
    }
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(
        file_reader, row_group_index, clause->column_index, &err);
    if (!bf) return true;  /* No bloom filter for this column ⇒ can't prune. */

    bool any = false;
    if (clause->op == CARQUET_FILTER_EQ) {
        any = bloom_check_value(bf, ti, clause->value, clause->value_size);
    } else {  /* CARQUET_FILTER_IN: keep if ANY listed value might be present. */
        for (int32_t i = 0; i < clause->value_count; i++) {
            const uint8_t* v = NULL;
            int32_t vlen = 0;
            in_value_at(ti, clause, i, &v, &vlen);
            if (v && bloom_check_value(bf, ti, v, vlen)) { any = true; break; }
        }
    }
    carquet_bloom_filter_destroy(bf);
    return any;
}

/* ============================================================================
 * Per-clause row-range derivation
 * ============================================================================ */

static carquet_status_t eval_clause_to_ranges(
    carquet_reader_t* file_reader,
    int32_t row_group_index,
    int64_t row_group_num_rows,
    const carquet_filter_clause_t* clause,
    carquet_row_range_list_t* out,
    carquet_error_t* error) {

    column_type_info_t ti;
    carquet_status_t st = lookup_column_type(file_reader,
        clause->column_index, &ti, error);
    if (st != CARQUET_OK) return st;

    /* Row-group-level pruning first: if the ColumnChunk statistics or the
     * bloom filter prove this row group cannot satisfy the clause, return an
     * empty range list. The caller intersects per-clause lists, so an empty
     * list skips the whole row group — with no page-index access at all. */
    carquet_column_statistics_t rg_stats;
    if (carquet_reader_column_statistics(file_reader, row_group_index,
            clause->column_index, &rg_stats) == CARQUET_OK) {
        if (!rg_stats_might_match(&ti, clause, &rg_stats, row_group_num_rows)) {
            carquet_row_range_list_clear(out);
            return CARQUET_OK;
        }
    }
    if (!rg_bloom_might_contain(file_reader, row_group_index, clause, &ti)) {
        carquet_row_range_list_clear(out);
        return CARQUET_OK;
    }

    carquet_column_index_t* ci = carquet_reader_get_column_index(
        file_reader, row_group_index, clause->column_index, error);
    carquet_offset_index_t* oi = carquet_reader_get_offset_index(
        file_reader, row_group_index, clause->column_index, error);

    if (!oi) {
        if (ci) carquet_column_index_free(ci);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_PAGE_INDEX_REQUIRED,
            "Page index missing for column %d in row group %d",
            clause->column_index, row_group_index);
        return CARQUET_ERROR_PAGE_INDEX_REQUIRED;
    }
    if (!ci) {
        /* OffsetIndex without ColumnIndex: the writer had no page bounds to
         * publish (statistics off, undefined sort order, ...). Nothing can be
         * proven absent, so every row of the group is a candidate. */
        carquet_offset_index_free(oi);
        carquet_row_range_list_clear(out);
        if (row_group_num_rows > 0) {
            return carquet_row_range_list_append(out, 0, row_group_num_rows);
        }
        return CARQUET_OK;
    }

    int32_t n_ci_pages = carquet_column_index_num_pages(ci);
    int32_t n_oi_pages = carquet_offset_index_num_pages(oi);
    int32_t n_pages = n_ci_pages < n_oi_pages ? n_ci_pages : n_oi_pages;

    carquet_row_range_list_clear(out);

    for (int32_t i = 0; i < n_pages; i++) {
        carquet_page_stats_t stats;
        if (carquet_column_index_get_page_stats(ci, i, &stats) != CARQUET_OK) {
            continue;
        }
        carquet_page_location_t loc;
        if (carquet_offset_index_get_page_location(oi, i, &loc) != CARQUET_OK) {
            continue;
        }

        int64_t first_row = loc.first_row_index;
        int64_t end_row;
        if (i + 1 < n_pages) {
            carquet_page_location_t next;
            (void)carquet_offset_index_get_page_location(oi, i + 1, &next);
            end_row = next.first_row_index;
        } else {
            end_row = row_group_num_rows;
        }
        int64_t page_rows = end_row - first_row;
        if (page_rows <= 0) continue;

        if (page_matches_clause(&ti, &stats, clause, first_row, page_rows)) {
            st = carquet_row_range_list_append(out, first_row, page_rows);
            if (st != CARQUET_OK) {
                carquet_column_index_free(ci);
                carquet_offset_index_free(oi);
                return st;
            }
        }
    }

    carquet_column_index_free(ci);
    carquet_offset_index_free(oi);
    return CARQUET_OK;
}

/* ============================================================================
 * Range intersection
 * ============================================================================ */

static carquet_status_t intersect_lists(
    const carquet_row_range_list_t* a,
    const carquet_row_range_list_t* b,
    carquet_row_range_list_t* out) {

    carquet_row_range_list_clear(out);
    int32_t i = 0, j = 0;
    while (i < a->count && j < b->count) {
        int64_t a_lo = a->ranges[i].first_row;
        int64_t a_hi = a_lo + a->ranges[i].num_rows;
        int64_t b_lo = b->ranges[j].first_row;
        int64_t b_hi = b_lo + b->ranges[j].num_rows;

        int64_t lo = a_lo > b_lo ? a_lo : b_lo;
        int64_t hi = a_hi < b_hi ? a_hi : b_hi;
        if (lo < hi) {
            carquet_status_t st = carquet_row_range_list_append(
                out, lo, hi - lo);
            if (st != CARQUET_OK) return st;
        }
        if (a_hi <= b_hi) i++;
        else j++;
    }
    return CARQUET_OK;
}

/* ============================================================================
 * Public entry point
 * ============================================================================ */

carquet_status_t carquet_page_filter_eval_row_group(
    carquet_reader_t* file_reader,
    int32_t row_group_index,
    const carquet_filter_clause_t* clauses,
    int32_t clause_count,
    carquet_row_range_list_t* out_ranges,
    carquet_error_t* error) {

    carquet_row_range_list_clear(out_ranges);
    if (clause_count <= 0) return CARQUET_OK;

    int32_t num_rg = carquet_reader_num_row_groups(file_reader);
    if (row_group_index < 0 || row_group_index >= num_rg) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_ROW_GROUP_NOT_FOUND,
            "Row group %d out of range", row_group_index);
        return CARQUET_ERROR_ROW_GROUP_NOT_FOUND;
    }

    int64_t row_group_num_rows =
        file_reader->metadata.row_groups[row_group_index].num_rows;
    if (row_group_num_rows <= 0) {
        return CARQUET_OK;  /* Empty row group ⇒ empty range list. */
    }

    /* Per-clause range lists, then iteratively intersect. */
    carquet_row_range_list_t accum;
    carquet_row_range_list_t scratch;
    carquet_row_range_list_t next;
    carquet_row_range_list_init(&accum);
    carquet_row_range_list_init(&scratch);
    carquet_row_range_list_init(&next);

    carquet_status_t st = eval_clause_to_ranges(file_reader, row_group_index,
        row_group_num_rows, &clauses[0], &accum, error);
    if (st != CARQUET_OK) goto done;

    for (int32_t c = 1; c < clause_count; c++) {
        if (accum.count == 0) break;  /* Already empty: short-circuit. */
        st = eval_clause_to_ranges(file_reader, row_group_index,
            row_group_num_rows, &clauses[c], &scratch, error);
        if (st != CARQUET_OK) goto done;

        st = intersect_lists(&accum, &scratch, &next);
        if (st != CARQUET_OK) goto done;

        /* Swap accum <- next, scratch keeps its capacity for reuse. */
        carquet_row_range_list_t tmp = accum;
        accum = next;
        next = tmp;
        carquet_row_range_list_clear(&next);
        carquet_row_range_list_clear(&scratch);
    }

    /* Move accum into out_ranges. */
    for (int32_t k = 0; k < accum.count; k++) {
        st = carquet_row_range_list_append(out_ranges,
            accum.ranges[k].first_row, accum.ranges[k].num_rows);
        if (st != CARQUET_OK) goto done;
    }

done:
    carquet_row_range_list_destroy(&accum);
    carquet_row_range_list_destroy(&scratch);
    carquet_row_range_list_destroy(&next);
    if (st != CARQUET_OK) {
        carquet_row_range_list_clear(out_ranges);
    }
    return st;
}
