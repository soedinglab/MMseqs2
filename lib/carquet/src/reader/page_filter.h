/**
 * @file page_filter.h
 * @brief Internal page-filter evaluation for the batch reader.
 *
 * Translates a conjunctive carquet_filter_clause_t[] into a sorted, non-
 * overlapping list of row ranges that must be materialized for a given row
 * group. Filtering reads only the column index + offset index of each
 * referenced column; predicate columns are never decompressed unless they
 * are also projected.
 */

#ifndef CARQUET_READER_PAGE_FILTER_H
#define CARQUET_READER_PAGE_FILTER_H

#include <carquet/carquet.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------------
 * Row-range list
 * ------------------------------------------------------------------------- */

typedef struct carquet_row_range {
    int64_t first_row;
    int64_t num_rows;
} carquet_row_range_t;

typedef struct carquet_row_range_list {
    carquet_row_range_t* ranges;
    int32_t count;
    int32_t capacity;
    int64_t total_rows;
} carquet_row_range_list_t;

void carquet_row_range_list_init(carquet_row_range_list_t* list);
void carquet_row_range_list_destroy(carquet_row_range_list_t* list);
void carquet_row_range_list_clear(carquet_row_range_list_t* list);
carquet_status_t carquet_row_range_list_append(
    carquet_row_range_list_t* list, int64_t first_row, int64_t num_rows);

/* ----------------------------------------------------------------------------
 * Clause validation (no row-group access required)
 * ------------------------------------------------------------------------- */

/**
 * Validate a filter clause against a reader's schema. Used by
 * set_page_filter() to surface synchronous errors. Returns
 * CARQUET_ERROR_INVALID_ARGUMENT for an out-of-range column, INT96, or a
 * size mismatch. Does NOT check column-index presence (that is checked
 * lazily during row-group evaluation).
 */
carquet_status_t carquet_page_filter_validate_clause(
    carquet_reader_t* file_reader,
    const carquet_filter_clause_t* clause,
    carquet_error_t* error);

/* ----------------------------------------------------------------------------
 * Row-group evaluation
 * ------------------------------------------------------------------------- */

/**
 * Compute matching row ranges for one row group, given an AND'd list of
 * clauses. On success out_ranges holds the sorted, non-overlapping ranges
 * (possibly empty if no rows match).
 *
 * Returns CARQUET_ERROR_PAGE_INDEX_REQUIRED if any referenced column lacks
 * a column index or offset index.
 */
carquet_status_t carquet_page_filter_eval_row_group(
    carquet_reader_t* file_reader,
    int32_t row_group_index,
    const carquet_filter_clause_t* clauses,
    int32_t clause_count,
    carquet_row_range_list_t* out_ranges,
    carquet_error_t* error);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_READER_PAGE_FILTER_H */
