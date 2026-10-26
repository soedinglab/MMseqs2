/**
 * @file page_tasks.h
 * @brief Page-granular parallel encode for the row-group finalize
 *
 * A deferred-encode column (see column_writer.c, defer_encode) holds its raw
 * input verbatim until the row group is finalized. For the common flat case
 * (REQUIRED, fixed-stride, PLAIN or BYTE_STREAM_SPLIT, no bloom filter or page
 * index) the page boundaries the eager encoder would pick depend only on the
 * value count, so the stash can be cut into pages up front and every page can
 * be encoded, CRC'd and compressed independently on any worker thread. The
 * row-group writer schedules one task per page across all columns, then each
 * column stitches its pages back together in order. Output is byte-identical
 * to the serial path: same values per page, same code path per page.
 */

#ifndef CARQUET_WRITER_PAGE_TASKS_H
#define CARQUET_WRITER_PAGE_TASKS_H

#include <carquet/error.h>
#include "core/buffer.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Header overhead added to a page's accumulated bytes when deciding whether
 * it has reached the target page size. Shared by the page writer's size
 * estimate and the page planner so both pick the same boundaries. */
#define CARQUET_PAGE_SIZE_ESTIMATE_OVERHEAD 64

typedef struct carquet_column_writer_internal carquet_column_writer_internal_t;
typedef struct carquet_page_writer carquet_page_writer_t;

typedef struct carquet_page_task {
    carquet_column_writer_internal_t* column;
    int64_t first_value;        /* Index into the column's deferred stash */
    int64_t num_values;
    carquet_buffer_t* out;      /* Encoded page: header + payload */

    /* Results, filled by carquet_column_writer_run_page_task */
    carquet_status_t status;
    size_t page_size;
    int32_t uncompressed_size;
    int32_t compressed_size;
    int64_t null_count;
    int64_t def_hist0;          /* Level histograms: a flat REQUIRED column has one bucket */
    int64_t rep_hist0;
    bool has_stats;
    /* Wall-clock seconds the encode of this page took, recorded by the
     * scheduler so the next row group can order its work (see
     * finalize_columns_parallel). */
    double elapsed;
    uint8_t* min_value;
    size_t min_size;
    size_t min_capacity;
    uint8_t* max_value;
    size_t max_size;
    size_t max_capacity;
} carquet_page_task_t;

/**
 * Cut the column's deferred stash into page tasks exactly where the eager
 * encoder would flush. Returns the number of pages, or 0 when the column is
 * not eligible (it then finalizes through the serial path). With @p tasks
 * non-NULL, fills up to @p capacity entries (column, first_value, num_values);
 * the return value is always the full count so a caller can size and retry.
 */
int32_t carquet_column_writer_plan_page_tasks(
    carquet_column_writer_internal_t* writer,
    carquet_page_task_t* tasks,
    int32_t capacity);

/**
 * Encode one page on the calling thread using @p scratch, a page writer that
 * belongs to this thread (see carquet_page_writer_create_scratch). Results
 * land in the task; task->out receives the page bytes.
 */
void carquet_column_writer_run_page_task(
    carquet_page_task_t* task,
    carquet_page_writer_t* scratch);

/**
 * Fold the column's page tasks into the column totals in order, exactly as
 * flushing each page would, and consume the deferred stash. With
 * @p copy_pages the page bytes are also appended to the column's chunk
 * buffer; without it they stay in the task buffers (the caller writes them
 * from there, see carquet_row_group_writer_column_pages) and only their total
 * size is reported in @p page_bytes. Must run on one thread per column.
 */
carquet_status_t carquet_column_writer_assemble_page_tasks(
    carquet_column_writer_internal_t* writer,
    carquet_page_task_t* tasks,
    int32_t count,
    bool copy_pages,
    size_t* page_bytes);

/** Release the min/max copies a task may hold. The task struct itself and
 *  task->out are owned by the caller. */
void carquet_page_task_release(carquet_page_task_t* task);

/** A page writer with no column binding, for use as per-thread scratch. */
carquet_page_writer_t* carquet_page_writer_create_scratch(void);

/**
 * Point @p dst at the same column configuration as @p src (type, encoding,
 * codec, level, CRC/statistics/V2 flags) and reset it. Both must be flat
 * (max definition and repetition level 0).
 */
carquet_status_t carquet_page_writer_adopt_config(
    carquet_page_writer_t* dst,
    const carquet_page_writer_t* src);

#endif /* CARQUET_WRITER_PAGE_TASKS_H */
