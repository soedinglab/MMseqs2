/**
 * @file file_writer.c
 * @brief Parquet file writing implementation
 *
 * Manages writing a complete Parquet file including:
 * - File header (PAR1 magic)
 * - Row groups via row_group_writer
 * - File metadata serialization
 * - Footer with metadata size and PAR1 magic
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include <carquet/error.h>
#include "core/buffer.h"
#include "core/arena.h"
#include "core/compat.h"
#include "reader/reader_internal.h"
#include "reader/worker_pool.h"
#include "thrift/thrift_encode.h"
#include "thrift/parquet_types.h"
#include "writer/arrow_schema.h"
#include "writer/page_tasks.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Parquet magic bytes */
static const uint8_t PARQUET_MAGIC[4] = {'P', 'A', 'R', '1'};
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                 int16_t max_def_level);

/* Forward declaration from row_group_writer.c */
typedef struct carquet_row_group_writer carquet_row_group_writer_t;

typedef struct column_chunk_info {
    int64_t file_offset;
    int64_t total_compressed_size;
    int64_t total_uncompressed_size;
    int64_t num_values;
    carquet_physical_type_t type;
    carquet_logical_type_t logical_type;
    carquet_encoding_t encoding;
    carquet_compression_t compression;
    int32_t type_length;
    char* path;
    /* Aggregated column statistics (mirrors row_group_writer.c). Min and max
     * are heap-allocated and may have different sizes (BYTE_ARRAY). */
    bool has_min_max;
    uint8_t* min_value;
    size_t min_value_size;
    uint8_t* max_value;
    size_t max_value_size;
    int64_t null_count;
    bool has_null_count;
    /* Dictionary page plumbing (mirrors row_group_writer.c). */
    bool has_dictionary_page;
    int64_t dictionary_page_size;
    /* GeospatialStatistics (mirrors row_group_writer.c). */
    bool has_geo_stats;
    parquet_geospatial_statistics_t geo_stats;
    /* Exact distinct non-null count (mirrors row_group_writer.c). */
    bool has_distinct_count;
    int64_t distinct_count;
    /* SizeStatistics (mirrors row_group_writer.c); histogram pointers alias the
     * column writer's buffers, copied out immediately below. */
    int64_t unencoded_ba_bytes;
    const int64_t* rep_level_hist;
    int32_t rep_hist_len;
    const int64_t* def_level_hist;
    int32_t def_hist_len;
} column_chunk_info_t;

extern carquet_row_group_writer_t* carquet_row_group_writer_create(
    const carquet_schema_t* schema,
    carquet_compression_t compression,
    size_t target_page_size,
    int64_t file_offset);

extern void carquet_row_group_writer_destroy(carquet_row_group_writer_t* writer);
extern void carquet_row_group_writer_reset(
    carquet_row_group_writer_t* writer,
    int64_t file_offset);

extern carquet_status_t carquet_row_group_writer_add_column(
    carquet_row_group_writer_t* writer,
    const char* name,
    carquet_physical_type_t type,
    const carquet_logical_type_t* logical_type,
    int16_t max_def_level,
    int16_t max_rep_level,
    int32_t type_length,
    carquet_encoding_t encoding,
    carquet_compression_t compression,
    int32_t compression_level);

extern void carquet_row_group_writer_configure_column_bloom(
    carquet_row_group_writer_t* writer,
    int column_index, bool enabled, int64_t ndv, double fpp);
extern void carquet_row_group_writer_set_column_max_rows_per_page(
    carquet_row_group_writer_t* writer,
    int column_index, int64_t max_rows);
extern void carquet_row_group_writer_set_column_write_batch_size(
    carquet_row_group_writer_t* writer,
    int column_index, int64_t batch_size);
extern void carquet_row_group_writer_set_column_page_size(
    carquet_row_group_writer_t* writer,
    int column_index, int64_t bytes);
extern void carquet_row_group_writer_set_column_data_page_v2(
    carquet_row_group_writer_t* writer,
    int column_index, bool enabled);

extern carquet_status_t carquet_row_group_writer_write_column(
    carquet_row_group_writer_t* writer,
    int column_index,
    const void* values,
    int64_t num_values,
    const int16_t* def_levels,
    const int16_t* rep_levels);

extern carquet_status_t carquet_row_group_writer_finalize(
    carquet_row_group_writer_t* writer,
    const uint8_t** data,
    size_t* size,
    int64_t num_rows);

extern carquet_status_t carquet_row_group_writer_finalize_columns(
    carquet_row_group_writer_t* writer,
    size_t* total_size,
    int64_t num_rows);
extern void carquet_row_group_writer_set_page_io(
    carquet_row_group_writer_t* writer, bool enabled);
extern void carquet_page_writer_destroy(carquet_page_writer_t* writer);
extern int32_t carquet_row_group_writer_begin_page_tasks(
    carquet_row_group_writer_t* writer);
extern carquet_page_task_t* carquet_row_group_writer_page_task(
    carquet_row_group_writer_t* writer, int32_t k);
extern carquet_status_t carquet_row_group_writer_complete_page_tasks(
    carquet_row_group_writer_t* writer, int64_t file_offset, int64_t num_rows,
    size_t* total_size);
extern const carquet_buffer_t* carquet_row_group_writer_column_pages(
    const carquet_row_group_writer_t* writer, int column_index, int32_t* count);
extern void carquet_row_group_writer_detach_column(
    carquet_row_group_writer_t* writer,
    int column_index,
    carquet_buffer_t* out,
    carquet_buffer_t* spare);
extern carquet_status_t carquet_row_group_writer_write_to_file(
    carquet_row_group_writer_t* writer,
    FILE* file,
    size_t* total_size,
    int64_t num_rows);

extern int carquet_row_group_writer_num_columns(const carquet_row_group_writer_t* writer);
extern int64_t carquet_row_group_writer_num_rows(const carquet_row_group_writer_t* writer);
extern int64_t carquet_row_group_writer_total_byte_size(const carquet_row_group_writer_t* writer);
extern const column_chunk_info_t* carquet_row_group_writer_get_column_info(
    const carquet_row_group_writer_t* writer, int index);

extern void carquet_row_group_writer_set_options(
    carquet_row_group_writer_t* writer,
    bool write_bloom_filters, bool write_page_index,
    bool write_statistics,
    bool write_crc,
    int32_t compression_level,
    int64_t dictionary_page_size);

/* Bloom filter and page index accessors */
typedef struct carquet_bloom_filter carquet_bloom_filter_t;
typedef struct carquet_column_index_builder carquet_column_index_builder_t;
typedef struct carquet_offset_index_builder carquet_offset_index_builder_t;

extern carquet_bloom_filter_t* carquet_row_group_writer_get_bloom_filter(
    const carquet_row_group_writer_t* writer, int index);
extern carquet_column_index_builder_t* carquet_row_group_writer_get_column_index(
    const carquet_row_group_writer_t* writer, int index);
extern carquet_offset_index_builder_t* carquet_row_group_writer_get_offset_index(
    const carquet_row_group_writer_t* writer, int index);

extern const uint8_t* carquet_bloom_filter_data(const carquet_bloom_filter_t* filter);
extern size_t carquet_bloom_filter_size(const carquet_bloom_filter_t* filter);
extern carquet_status_t carquet_column_index_serialize(
    const carquet_column_index_builder_t* builder, carquet_buffer_t* output);
extern carquet_status_t carquet_offset_index_serialize(
    const carquet_offset_index_builder_t* builder, carquet_buffer_t* output);

/* ============================================================================
 * Writer Schema Structure (for building)
 * ============================================================================
 */

typedef struct writer_column_def {
    char* name;
    carquet_physical_type_t physical_type;
    carquet_logical_type_t logical_type;
    carquet_field_repetition_t repetition;
    int32_t type_length;
    int16_t max_def_level;
    int16_t max_rep_level;
    bool statistics_sort_order_defined;
} writer_column_def_t;

/* ============================================================================
 * Row Group Metadata Storage
 * ============================================================================
 */

typedef struct row_group_column_info {
    int64_t file_offset;
    int64_t total_compressed_size;
    int64_t total_uncompressed_size;
    int64_t num_values;
    carquet_physical_type_t type;
    carquet_compression_t codec;
    int64_t data_page_offset;
    bool has_dictionary_page_offset;
    int64_t dictionary_page_offset;
    bool has_dictionary_page;
    /* Per-chunk ColumnMetaData.encodings, derived from whether this chunk
     * actually emitted a dictionary page (post-finalize), not from the static
     * configured-encoding cache. A dictionary chunk advertises
     * {PLAIN, RLE_DICTIONARY, RLE}; a plain/dict-fallback chunk advertises
     * {<data encoding>, RLE}. */
    carquet_encoding_t encodings[3];
    int32_t num_encodings;
    bool has_bloom_filter_offset;
    int64_t bloom_filter_offset;
    bool has_bloom_filter_length;
    int32_t bloom_filter_length;
    bool has_column_index_offset;
    int64_t column_index_offset;
    bool has_column_index_length;
    int32_t column_index_length;
    bool has_offset_index_offset;
    int64_t offset_index_offset;
    bool has_offset_index_length;
    int32_t offset_index_length;
    /* Aggregated column statistics */
    bool has_statistics;
    bool has_min_max;
    bool has_null_count;
    int64_t null_count;
    /* Heap-allocated; freed when the row_group_info_t is torn down. */
    uint8_t* min_value;
    int32_t min_value_size;
    uint8_t* max_value;
    int32_t max_value_size;
    /* GeospatialStatistics (GEOMETRY/GEOGRAPHY) */
    bool has_geo_stats;
    parquet_geospatial_statistics_t geo_stats;
    /* Exact distinct non-null count, when known (dictionary-encoded chunks). */
    bool has_distinct_count;
    int64_t distinct_count;
    /* SizeStatistics (Parquet 2.9); owned copies, freed with this info.
     * unencoded_ba_bytes is -1 for non-BYTE_ARRAY columns. */
    bool has_size_statistics;
    int64_t unencoded_ba_bytes;
    int64_t* rep_level_hist;
    int32_t rep_hist_len;
    int64_t* def_level_hist;
    int32_t def_hist_len;
} row_group_column_info_t;

typedef struct row_group_info {
    int64_t file_offset;
    int64_t num_rows;
    int64_t total_byte_size;
    int64_t total_compressed_size;
    int16_t ordinal;
    row_group_column_info_t* columns;
    int32_t num_columns;
} row_group_info_t;

/* One background write: a finalized column chunk headed for the file. */
typedef struct io_page {
    const uint8_t* data;
    size_t size;
} io_page_t;

/* One column chunk queued for the I/O thread: the chunk buffer (whole chunk
 * of a serially encoded column, usually empty otherwise) followed by the
 * column's pages, which the row-group writer keeps in its task buffers until
 * the job has been drained (page_io mode, see row_group_writer.c). The page
 * descriptors are copied here so the job never points into an array the
 * next row group may grow. */
typedef struct io_job {
    carquet_buffer_t buf;
    io_page_t* pages;
    int32_t num_pages;
    int32_t page_capacity;
    FILE* file;
    carquet_status_t status;
} io_job_t;

/* Pipelined encode (see "Pipelined encode" below). */
#ifdef _WIN32
typedef struct pipeline_sync {
    CRITICAL_SECTION mutex;
    CONDITION_VARIABLE cond;
} pipeline_sync_t;
#else
typedef struct pipeline_sync {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} pipeline_sync_t;
#endif

struct carquet_writer;
struct pipeline_rg;

typedef struct encode_task_ctx {
    struct carquet_writer* writer;
    struct pipeline_rg* rg;
    carquet_page_task_t* task;
} encode_task_ctx_t;

/* A row group whose pages are being encoded on the encode pool, or have
 * been and await their turn to be written. */
typedef struct pipeline_rg {
    carquet_row_group_writer_t* rg;   /* NULL: slot free */
    int64_t rows;
    int32_t remaining;                /* tasks still running; under pipe_sync */
    encode_task_ctx_t* ctx;
    int32_t ctx_capacity;
} pipeline_rg_t;

/* ============================================================================
 * Writer Structure
 * ============================================================================
 */

struct carquet_writer {
    FILE* file;
    bool owns_file;
    /* True when this writer was opened over an existing file via
     * carquet_writer_open_append(). Such files must never be deleted on abort:
     * remove() would destroy the user's pre-existing data. */
    bool is_append;
    char* path;

    /* Schema */
    writer_column_def_t* columns;
    int32_t num_columns;
    int32_t column_capacity;

    /* Full schema elements (including groups) for metadata serialization */
    parquet_schema_element_t* schema_elements;
    int32_t num_schema_elements;
    char*** column_paths;
    int32_t* column_path_lens;
    /* Per-column encodings list for ColumnMetaData.encodings. For dictionary
     * columns this is {PLAIN, RLE_DICTIONARY, RLE}; for plain columns
     * {<data encoding>, RLE}. column_num_encodings holds the live count. */
    carquet_encoding_t (*column_encodings)[3];
    int32_t* column_num_encodings;

    /* Options */
    carquet_writer_options_t options;

    /* Current row group */
    carquet_row_group_writer_t* current_row_group;
    int64_t current_row_group_rows;
    int64_t* column_values_written;  /* Values written per column in current row group */
    int64_t current_row_group_estimated_bytes;

    /* Completed row groups */
    row_group_info_t* row_groups;
    int32_t num_row_groups;
    int32_t row_groups_capacity;

    /* File state */
    int64_t file_offset;
    int64_t total_rows;
    bool header_written;

    /* Arena for metadata allocations */
    carquet_arena_t arena;

    /* Key-value metadata */
    parquet_key_value_t* kv_metadata;
    int32_t num_kv_metadata;
    int32_t kv_metadata_capacity;

    /* Max bytes for BYTE_ARRAY min/max in column statistics. The Parquet spec
       recommends truncating; carquet's historical default of 32 matches what
       Arrow uses. Configurable via carquet_writer_set_max_statistics_size. */
    int64_t max_statistics_size;

    /* Per-column overrides. _set flags distinguish "not overridden" from an
       overriding value of 0 (PLAIN / UNCOMPRESSED are both 0 in their enums). */
    carquet_encoding_t* column_encoding_overrides;
    bool* column_encoding_override_set;
    carquet_compression_t* column_compression_overrides;
    int32_t* column_compression_levels;
    bool* column_compression_override_set;
    bool* column_statistics_overrides;
    bool* column_bloom_filter_overrides;
    /* Per-column bloom NDV/FPP overrides (parallel arrays). _set distinguishes
       "not overridden" from an explicit value. */
    int64_t* column_bloom_ndv_overrides;
    double* column_bloom_fpp_overrides;
    bool* column_bloom_options_set;
    /* Set when the user explicitly called carquet_writer_set_column_bloom_filter
       for this column (as opposed to the value merely defaulting to the global
       flag). Lets the finalize path keep an explicit legacy enable even after
       the newer ndv/fpp options API takes per-column control. */
    bool* column_bloom_explicit;
    /* Per-column page-size override (target bytes). 0 / unset means use
       options.page_size. */
    int64_t* column_page_size_overrides;
    bool* column_page_size_override_set;
    bool column_overrides_allocated;

    /* Sorting columns metadata, applied to every row group */
    parquet_sorting_column_t* sorting_columns;
    int32_t num_sorting_columns;

    /* Buffer writer support */
    bool is_buffer_writer;
    uint8_t* output_buffer;
    size_t output_buffer_size;

    /* Background I/O (options.async_io). The pool has one thread and runs
     * io_jobs in submission order; jobs own their buffers until the main
     * thread reclaims them into io_spares after a drain. io_error is sticky:
     * once a background write fails, every later file operation fails too. */
    carquet_worker_pool_t* io_pool;
    io_job_t* io_jobs;            /* [num_columns] */
    int32_t io_jobs_in_flight;

    /* Pipelined encode. Two row-group writers alternate: while the user
     * fills current_row_group, the previous row group encodes on encode_pool
     * from its own writer (a pipe_slots entry); the idle twin waits in
     * spare_row_group. */
    carquet_worker_pool_t* encode_pool;
    int32_t encode_threads;
    carquet_row_group_writer_t* spare_row_group;
    pipeline_rg_t pipe_slots[2];
    int pipe_pending;             /* slot of the row group to complete next, -1 none */
    pipeline_sync_t pipe_sync;
    bool pipe_sync_ready;
    carquet_page_writer_t** scratch_stack;  /* idle per-thread page writers; under pipe_sync */
    int32_t scratch_count;
    int32_t scratch_capacity;
    bool pipeline_off;            /* could not set up: stay on the fork/join path */
    carquet_buffer_t* io_spares;  /* reclaimed chunk buffers, LIFO */
    int32_t io_spare_count;
    carquet_status_t io_error;
};

/* ============================================================================
 * Writer Options
 * ============================================================================
 */

void carquet_writer_options_init(carquet_writer_options_t* options) {
    /* options is nonnull per API contract */
    memset(options, 0, sizeof(*options));
    options->compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    options->compression_level = 0;
    options->row_group_size = 128 * 1024 * 1024;  /* 128 MB */
    options->page_size = 1024 * 1024;               /* 1 MB */
    options->write_statistics = true;
    options->write_crc = true;
    options->write_page_index = false;
    options->write_bloom_filters = false;
    options->dictionary_encoding = CARQUET_ENCODING_RLE_DICTIONARY;
    options->dictionary_page_size = 1024 * 1024;   /* 1 MB */
    options->created_by = "Carquet";
    options->max_rows_per_page = 0;                 /* unlimited */
    options->write_arrow_schema = false;
    options->data_page_version = 1;
    options->coerce_timestamps = false;
    options->coerce_timestamp_unit = CARQUET_TIME_UNIT_MICROS;
    options->allow_timestamp_truncation = false;
    options->write_batch_size = 0;
    options->file_format_version = 2;
    options->async_io = true;
}

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static carquet_status_t write_magic(FILE* file) {
    if (fwrite(PARQUET_MAGIC, 1, 4, file) != 4) {
        return CARQUET_ERROR_FILE_WRITE;
    }
    return CARQUET_OK;
}

static int64_t saturating_add_i64(int64_t lhs, int64_t rhs) {
    if (rhs <= 0 || lhs >= INT64_MAX - rhs) {
        return INT64_MAX;
    }
    return lhs + rhs;
}

static int64_t estimate_column_batch_bytes(
    const writer_column_def_t* column,
    const void* values,
    int64_t num_values,
    const int16_t* def_levels,
    const int16_t* rep_levels) {

    if (!column || !values || num_values <= 0) {
        return 0;
    }

    int64_t total = 0;

    switch (column->physical_type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            total = num_values * (int64_t)sizeof(uint8_t);
            break;
        case CARQUET_PHYSICAL_INT32:
            total = num_values * (int64_t)sizeof(int32_t);
            break;
        case CARQUET_PHYSICAL_INT64:
            total = num_values * (int64_t)sizeof(int64_t);
            break;
        case CARQUET_PHYSICAL_FLOAT:
            total = num_values * (int64_t)sizeof(float);
            break;
        case CARQUET_PHYSICAL_DOUBLE:
            total = num_values * (int64_t)sizeof(double);
            break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (column->type_length > 0) {
                total = num_values * (int64_t)column->type_length;
            }
            break;
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            const carquet_byte_array_t* arrays = (const carquet_byte_array_t*)values;
            /* For nullable columns `values` holds only the packed non-null
             * entries, so iterate the present count — not num_values (the
             * logical row count) — to avoid reading past the array. */
            int64_t value_count = num_values;
            if (def_levels && column->max_def_level > 0) {
                value_count = carquet_dispatch_count_non_nulls(
                    def_levels, num_values, column->max_def_level);
            }
            /* Accumulate in a local with a single saturation check at the end:
             * each length is <= UINT32_MAX, so value_count entries cannot overflow
             * uint64_t until ~2^32 values, which cannot be reached here. */
            uint64_t chunk_total = 0;
            for (int64_t i = 0; i < value_count; i++) {
                chunk_total += (uint64_t)sizeof(uint32_t) + arrays[i].length;
            }
            total = (chunk_total > (uint64_t)INT64_MAX)
                        ? INT64_MAX
                        : saturating_add_i64(total, (int64_t)chunk_total);
            break;
        }
        default:
            break;
    }

    if (column->max_def_level > 0 && def_levels) {
        total = saturating_add_i64(total, num_values * (int64_t)sizeof(*def_levels));
    }
    if (column->max_rep_level > 0 && rep_levels) {
        total = saturating_add_i64(total, num_values * (int64_t)sizeof(*rep_levels));
    }

    return total;
}

static bool writer_supports_aligned_auto_flush(const carquet_writer_t* writer) {
    if (!writer || writer->options.row_group_size <= 0 || writer->num_columns <= 0) {
        return false;
    }

    for (int32_t i = 0; i < writer->num_columns; i++) {
        if (writer->columns[i].max_rep_level > 0) {
            return false;
        }
    }

    return true;
}

static bool current_row_group_is_aligned(const carquet_writer_t* writer) {
    if (!writer || writer->current_row_group_rows <= 0) {
        return false;
    }

    int64_t expected_rows = writer->current_row_group_rows;
    for (int32_t i = 0; i < writer->num_columns; i++) {
        if (writer->column_values_written[i] != expected_rows) {
            return false;
        }
    }

    return true;
}

static carquet_status_t ensure_header_written(carquet_writer_t* writer) {
    if (writer->header_written) {
        return CARQUET_OK;
    }

    carquet_status_t status = write_magic(writer->file);
    if (status != CARQUET_OK) {
        return status;
    }

    writer->file_offset = 4;  /* PAR1 magic */
    writer->header_written = true;
    return CARQUET_OK;
}

static carquet_status_t add_column_internal(
    carquet_writer_t* writer,
    const char* name,
    carquet_physical_type_t physical_type,
    const carquet_logical_type_t* logical_type,
    carquet_field_repetition_t repetition,
    int32_t type_length,
    int16_t max_def_level,
    int16_t max_rep_level,
    bool statistics_sort_order_defined) {

    /* Expand capacity if needed */
    if (writer->num_columns >= writer->column_capacity) {
        int32_t new_cap = writer->column_capacity == 0 ? 8 : writer->column_capacity * 2;
        writer_column_def_t* new_cols = carquet_mem_realloc(writer->columns,
            new_cap * sizeof(writer_column_def_t));
        if (!new_cols) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        writer->columns = new_cols;

        int64_t* new_values = carquet_mem_realloc(writer->column_values_written,
            new_cap * sizeof(int64_t));
        if (!new_values) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        writer->column_values_written = new_values;

        writer->column_capacity = new_cap;
    }

    writer_column_def_t* col = &writer->columns[writer->num_columns];
    memset(col, 0, sizeof(*col));

    col->name = carquet_heap_strdup(name);
    if (!col->name) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    col->physical_type = physical_type;
    col->repetition = repetition;
    col->type_length = type_length;

    if (logical_type) {
        col->logical_type = *logical_type;
    }

    col->max_def_level = max_def_level;
    col->max_rep_level = max_rep_level;
    col->statistics_sort_order_defined = statistics_sort_order_defined;

    writer->column_values_written[writer->num_columns] = 0;
    writer->num_columns++;

    return CARQUET_OK;
}

/* Store the full schema elements (including groups) for metadata serialization */
/* Free the heap-owned per-field metadata deep-copied by store_schema_elements. */
static void free_schema_field_metadata(parquet_schema_element_t* elements,
                                       int32_t num_elements) {
    if (!elements) return;
    for (int32_t i = 0; i < num_elements; i++) {
        for (int32_t j = 0; j < elements[i].num_field_metadata; j++) {
            carquet_mem_free(elements[i].field_metadata[j].key);
            carquet_mem_free(elements[i].field_metadata[j].value);
        }
        carquet_mem_free(elements[i].field_metadata);
        elements[i].field_metadata = NULL;
        elements[i].num_field_metadata = 0;
    }
}

static carquet_status_t store_schema_elements(
    carquet_writer_t* writer,
    const carquet_schema_t* schema) {

    writer->num_schema_elements = schema->num_elements;
    writer->schema_elements = carquet_mem_calloc(schema->num_elements, sizeof(parquet_schema_element_t));
    if (!writer->schema_elements) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    bool coerce = writer->options.coerce_timestamps;
    carquet_time_unit_t tgt = writer->options.coerce_timestamp_unit;

    for (int32_t i = 0; i < schema->num_elements; i++) {
        writer->schema_elements[i] = schema->elements[i];
        if (schema->elements[i].name) {
            writer->schema_elements[i].name = carquet_heap_strdup(schema->elements[i].name);
            if (!writer->schema_elements[i].name) {
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
        }

        /* Deep-copy per-field metadata (Arrow custom_metadata) into heap the
         * writer owns, so it survives even if the caller frees the schema. The
         * copy is released in the writer teardown alongside name. */
        writer->schema_elements[i].field_metadata = NULL;
        writer->schema_elements[i].num_field_metadata = 0;
        int32_t nmeta = schema->elements[i].num_field_metadata;
        if (nmeta > 0 && schema->elements[i].field_metadata) {
            parquet_key_value_t* dst = carquet_mem_calloc(
                (size_t)nmeta, sizeof(parquet_key_value_t));
            if (!dst) return CARQUET_ERROR_OUT_OF_MEMORY;
            for (int32_t j = 0; j < nmeta; j++) {
                const parquet_key_value_t* src = &schema->elements[i].field_metadata[j];
                dst[j].key = src->key ? carquet_heap_strdup(src->key) : NULL;
                dst[j].value = src->value ? carquet_heap_strdup(src->value) : NULL;
                if ((src->key && !dst[j].key) || (src->value && !dst[j].value)) {
                    for (int32_t k = 0; k <= j; k++) {
                        carquet_mem_free(dst[k].key);
                        carquet_mem_free(dst[k].value);
                    }
                    carquet_mem_free(dst);
                    return CARQUET_ERROR_OUT_OF_MEMORY;
                }
            }
            writer->schema_elements[i].field_metadata = dst;
            writer->schema_elements[i].num_field_metadata = nmeta;
        }

        /* Coerce TIMESTAMP units in the emitted schema so file metadata
         * reflects the target unit (the values are rescaled on write). */
        if (coerce) {
            parquet_schema_element_t* e = &writer->schema_elements[i];
            if (e->has_logical_type &&
                e->logical_type.id == CARQUET_LOGICAL_TIMESTAMP) {
                e->logical_type.params.timestamp.unit = tgt;
                if (e->has_converted_type) {
                    if (tgt == CARQUET_TIME_UNIT_MILLIS) {
                        e->converted_type = CARQUET_CONVERTED_TIMESTAMP_MILLIS;
                    } else if (tgt == CARQUET_TIME_UNIT_MICROS) {
                        e->converted_type = CARQUET_CONVERTED_TIMESTAMP_MICROS;
                    } else {
                        /* NANOS has no legacy ConvertedType. */
                        e->has_converted_type = false;
                    }
                }
            }
        }
    }

    return CARQUET_OK;
}

static carquet_compression_t effective_column_compression(
    const carquet_writer_t* writer, int32_t column_index);
static int32_t effective_column_compression_level(
    const carquet_writer_t* writer, int32_t column_index);
static carquet_encoding_t effective_column_encoding(
    const carquet_writer_t* writer, int32_t column_index,
    carquet_physical_type_t type, carquet_compression_t compression);

static carquet_status_t build_column_metadata_cache(
    carquet_writer_t* writer,
    const carquet_schema_t* schema) {

    writer->column_paths = carquet_mem_calloc((size_t)schema->num_leaves, sizeof(char**));
    writer->column_path_lens = carquet_mem_calloc((size_t)schema->num_leaves, sizeof(int32_t));
    writer->column_encodings = carquet_mem_calloc((size_t)schema->num_leaves,
        sizeof(*writer->column_encodings));
    writer->column_num_encodings = carquet_mem_calloc((size_t)schema->num_leaves,
        sizeof(int32_t));

    if (!writer->column_paths || !writer->column_path_lens ||
        !writer->column_encodings || !writer->column_num_encodings) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < schema->num_leaves; i++) {
        int32_t elem_idx = schema->leaf_indices[i];
        int32_t depth = 0;

        for (int32_t cur = elem_idx; cur > 0; cur = schema->parent_indices[cur]) {
            depth++;
        }

        if (depth <= 0) {
            depth = 1;
        }

        writer->column_paths[i] = carquet_mem_calloc((size_t)depth, sizeof(char*));
        if (!writer->column_paths[i]) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        writer->column_path_lens[i] = depth;
        if (depth == 1) {
            writer->column_paths[i][0] = writer->schema_elements[elem_idx].name;
        } else {
            int32_t cur = elem_idx;
            for (int32_t pi = depth - 1; pi >= 0; pi--) {
                writer->column_paths[i][pi] = writer->schema_elements[cur].name;
                cur = schema->parent_indices[cur];
            }
        }
    }

    return CARQUET_OK;
}

static bool leaf_statistics_sort_order_defined(
    const carquet_schema_t* schema,
    int32_t leaf_elem_idx) {

    const parquet_schema_element_t* leaf = &schema->elements[leaf_elem_idx];
    if (leaf->has_logical_type) {
        switch (leaf->logical_type.id) {
            case CARQUET_LOGICAL_GEOMETRY:
            case CARQUET_LOGICAL_GEOGRAPHY:
            case CARQUET_LOGICAL_VARIANT:
                return false;
            default:
                break;
        }
    }

    for (int32_t cur = schema->parent_indices[leaf_elem_idx]; cur > 0;
         cur = schema->parent_indices[cur]) {
        const parquet_schema_element_t* elem = &schema->elements[cur];
        if (elem->has_logical_type &&
            elem->logical_type.id == CARQUET_LOGICAL_VARIANT) {
            return false;
        }
    }

    return true;
}

static carquet_status_t ensure_column_overrides(carquet_writer_t* writer) {
    if (writer->column_overrides_allocated) return CARQUET_OK;
    int32_t n = writer->num_columns;
    writer->column_encoding_overrides = carquet_mem_calloc(n, sizeof(carquet_encoding_t));
    writer->column_encoding_override_set = carquet_mem_calloc(n, sizeof(bool));
    writer->column_compression_overrides = carquet_mem_calloc(n, sizeof(carquet_compression_t));
    writer->column_compression_levels = carquet_mem_calloc(n, sizeof(int32_t));
    writer->column_compression_override_set = carquet_mem_calloc(n, sizeof(bool));
    writer->column_statistics_overrides = carquet_mem_calloc(n, sizeof(bool));
    writer->column_bloom_filter_overrides = carquet_mem_calloc(n, sizeof(bool));
    writer->column_bloom_ndv_overrides = carquet_mem_calloc(n, sizeof(int64_t));
    writer->column_bloom_fpp_overrides = carquet_mem_calloc(n, sizeof(double));
    writer->column_bloom_options_set = carquet_mem_calloc(n, sizeof(bool));
    writer->column_bloom_explicit = carquet_mem_calloc(n, sizeof(bool));
    writer->column_page_size_overrides = carquet_mem_calloc(n, sizeof(int64_t));
    writer->column_page_size_override_set = carquet_mem_calloc(n, sizeof(bool));
    if (!writer->column_encoding_overrides || !writer->column_encoding_override_set ||
        !writer->column_compression_overrides || !writer->column_compression_levels ||
        !writer->column_compression_override_set ||
        !writer->column_statistics_overrides || !writer->column_bloom_filter_overrides ||
        !writer->column_bloom_ndv_overrides || !writer->column_bloom_fpp_overrides ||
        !writer->column_bloom_options_set || !writer->column_bloom_explicit ||
        !writer->column_page_size_overrides || !writer->column_page_size_override_set) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    for (int32_t i = 0; i < n; i++) {
        writer->column_statistics_overrides[i] = writer->options.write_statistics;
        writer->column_bloom_filter_overrides[i] = writer->options.write_bloom_filters;
    }
    writer->column_overrides_allocated = true;
    return CARQUET_OK;
}

static void free_column_overrides(carquet_writer_t* writer) {
    carquet_mem_free(writer->column_encoding_overrides);
    carquet_mem_free(writer->column_encoding_override_set);
    carquet_mem_free(writer->column_compression_overrides);
    carquet_mem_free(writer->column_compression_levels);
    carquet_mem_free(writer->column_compression_override_set);
    carquet_mem_free(writer->column_statistics_overrides);
    carquet_mem_free(writer->column_bloom_filter_overrides);
    carquet_mem_free(writer->column_bloom_ndv_overrides);
    carquet_mem_free(writer->column_bloom_fpp_overrides);
    carquet_mem_free(writer->column_bloom_options_set);
    carquet_mem_free(writer->column_bloom_explicit);
    carquet_mem_free(writer->column_page_size_overrides);
    carquet_mem_free(writer->column_page_size_override_set);
    writer->column_encoding_overrides = NULL;
    writer->column_encoding_override_set = NULL;
    writer->column_compression_overrides = NULL;
    writer->column_compression_levels = NULL;
    writer->column_compression_override_set = NULL;
    writer->column_statistics_overrides = NULL;
    writer->column_bloom_filter_overrides = NULL;
    writer->column_bloom_ndv_overrides = NULL;
    writer->column_bloom_fpp_overrides = NULL;
    writer->column_bloom_options_set = NULL;
    writer->column_bloom_explicit = NULL;
    writer->column_page_size_overrides = NULL;
    writer->column_page_size_override_set = NULL;
    writer->column_overrides_allocated = false;
}

static bool any_column_bloom_options_set(const carquet_writer_t* writer) {
    if (!writer->column_overrides_allocated || !writer->column_bloom_options_set) {
        return false;
    }
    for (int32_t i = 0; i < writer->num_columns; i++) {
        if (writer->column_bloom_options_set[i]) return true;
    }
    return false;
}

static carquet_compression_t effective_column_compression(
    const carquet_writer_t* writer, int32_t column_index) {
    carquet_compression_t codec;
    if (writer->column_overrides_allocated &&
        writer->column_compression_override_set[column_index]) {
        codec = writer->column_compression_overrides[column_index];
    } else {
        codec = writer->options.compression;
    }

    /* LZ4 (codec 5) and LZ4_RAW (codec 7) are distinct Parquet codecs: codec 5
     * is the deprecated Hadoop-framed LZ4 (length-prefixed blocks), codec 7 is
     * raw LZ4 blocks. carquet honours the requested codec directly; the page
     * writer applies the matching framing. */
    return codec;
}

static int32_t effective_column_compression_level(
    const carquet_writer_t* writer, int32_t column_index) {
    if (writer->column_overrides_allocated &&
        writer->column_compression_override_set[column_index]) {
        return writer->column_compression_levels[column_index];
    }
    return writer->options.compression_level;
}

static carquet_encoding_t effective_column_encoding(
    const carquet_writer_t* writer,
    int32_t column_index,
    carquet_physical_type_t type,
    carquet_compression_t compression) {
    if (writer->column_overrides_allocated &&
        writer->column_encoding_override_set[column_index]) {
        return writer->column_encoding_overrides[column_index];
    }
    /* Default encoding policy (matches v0.4.4): PLAIN, with automatic
     * BYTE_STREAM_SPLIT for FLOAT/DOUBLE when a compression codec is set
     * (BSS makes the float byte planes far more compressible). Dictionary
     * encoding is a deliberate opt-in via carquet_writer_set_column_encoding()
     * — making it the default regressed read throughput badly (notably a
     * ~270x slowdown on the zero-copy uncompressed read path) for a size win
     * that is zero under zstd. The full dictionary writer remains available. */
    if (compression != CARQUET_COMPRESSION_UNCOMPRESSED &&
        (type == CARQUET_PHYSICAL_FLOAT || type == CARQUET_PHYSICAL_DOUBLE)) {
        return CARQUET_ENCODING_BYTE_STREAM_SPLIT;
    }
    return CARQUET_ENCODING_PLAIN;
}

/* Recompute the per-column ColumnMetaData.encodings cache. Must run after
 * per-column encoding/compression overrides are applied (the cache built at
 * writer-create time predates them), i.e. just before the first row group is
 * created. Dictionary columns advertise {PLAIN, RLE_DICTIONARY, RLE}; other
 * columns {<data encoding>, RLE}. This cache reflects the *configured*
 * encoding only; the actually-emitted ColumnMetaData.encodings is overridden
 * per chunk after finalize using row_group_column_info_t.encodings (see the
 * chunk-assembly loop), so a dict column that fell back to PLAIN advertises
 * only {PLAIN, RLE} and never the extra RLE_DICTIONARY entry. */
static void refresh_column_encodings_cache(carquet_writer_t* writer) {
    for (int32_t i = 0; i < writer->num_columns; i++) {
        carquet_physical_type_t phys = writer->columns[i].physical_type;
        carquet_compression_t col_comp = effective_column_compression(writer, i);
        carquet_encoding_t enc =
            effective_column_encoding(writer, i, phys, col_comp);
        if (enc == CARQUET_ENCODING_RLE_DICTIONARY ||
            enc == CARQUET_ENCODING_PLAIN_DICTIONARY) {
            writer->column_encodings[i][0] = CARQUET_ENCODING_PLAIN;
            writer->column_encodings[i][1] = CARQUET_ENCODING_RLE_DICTIONARY;
            writer->column_encodings[i][2] = CARQUET_ENCODING_RLE;
            writer->column_num_encodings[i] = 3;
        } else {
            writer->column_encodings[i][0] = enc;
            writer->column_encodings[i][1] = CARQUET_ENCODING_RLE;
            writer->column_num_encodings[i] = 2;
        }
    }
}

static bool writer_encoding_supported(
    carquet_encoding_t encoding,
    carquet_physical_type_t type) {

    switch (encoding) {
        case CARQUET_ENCODING_PLAIN:
            return true;
        case CARQUET_ENCODING_BYTE_STREAM_SPLIT:
            /* Parquet spec permits BYTE_STREAM_SPLIT for these physical
             * types (not just FLOAT/DOUBLE). */
            return type == CARQUET_PHYSICAL_FLOAT ||
                   type == CARQUET_PHYSICAL_DOUBLE ||
                   type == CARQUET_PHYSICAL_INT32 ||
                   type == CARQUET_PHYSICAL_INT64 ||
                   type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY;
        case CARQUET_ENCODING_DELTA_BINARY_PACKED:
            return type == CARQUET_PHYSICAL_INT32 ||
                   type == CARQUET_PHYSICAL_INT64;
        case CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY:
            return type == CARQUET_PHYSICAL_BYTE_ARRAY;
        case CARQUET_ENCODING_DELTA_BYTE_ARRAY:
            /* Spec permits DELTA_BYTE_ARRAY for BYTE_ARRAY and FLBA. */
            return type == CARQUET_PHYSICAL_BYTE_ARRAY ||
                   type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY;
        case CARQUET_ENCODING_RLE_DICTIONARY:
        case CARQUET_ENCODING_PLAIN_DICTIONARY:
            return type == CARQUET_PHYSICAL_INT32 ||
                   type == CARQUET_PHYSICAL_INT64 ||
                   type == CARQUET_PHYSICAL_FLOAT ||
                   type == CARQUET_PHYSICAL_DOUBLE ||
                   type == CARQUET_PHYSICAL_BYTE_ARRAY ||
                   type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY;
        case CARQUET_ENCODING_RLE:
            /* RLE as a value encoding is defined only for BOOLEAN. */
            return type == CARQUET_PHYSICAL_BOOLEAN;
        default:
            return false;
    }
}

static void free_row_groups(carquet_writer_t* writer) {
    if (!writer->row_groups) return;
    for (int32_t i = 0; i < writer->num_row_groups; i++) {
        row_group_info_t* rg = &writer->row_groups[i];
        if (rg->columns) {
            for (int32_t j = 0; j < rg->num_columns; j++) {
                carquet_mem_free(rg->columns[j].min_value);
                carquet_mem_free(rg->columns[j].max_value);
                carquet_mem_free(rg->columns[j].rep_level_hist);
                carquet_mem_free(rg->columns[j].def_level_hist);
            }
            carquet_mem_free(rg->columns);
        }
    }
    carquet_mem_free(writer->row_groups);
    writer->row_groups = NULL;
    writer->num_row_groups = 0;
}

static void free_kv_metadata(carquet_writer_t* writer) {
    if (writer->kv_metadata) {
        for (int32_t i = 0; i < writer->num_kv_metadata; i++) {
            carquet_mem_free(writer->kv_metadata[i].key);
            carquet_mem_free(writer->kv_metadata[i].value);
        }
        carquet_mem_free(writer->kv_metadata);
        writer->kv_metadata = NULL;
    }
}

/* A row-group writer configured for this file (columns, codecs, per-column
 * overrides). ensure_row_group makes the first; the pipelined encode makes a
 * second one the same way. */
static carquet_status_t create_row_group_writer(
    carquet_writer_t* writer, carquet_row_group_writer_t** out) {
    *out = NULL;

    /* Overrides may have been set after writer-create; refresh the encodings
     * cache so ColumnMetaData.encodings reflects the resolved encoding. */
    refresh_column_encodings_cache(writer);

    size_t target_page_size = (size_t)writer->options.page_size;

    carquet_row_group_writer_t* rg = carquet_row_group_writer_create(
        NULL,  /* Schema not used directly */
        writer->options.compression,
        target_page_size,
        writer->file_offset);

    if (!rg) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    /* Pass optional feature flags */
    carquet_row_group_writer_set_options(
        rg,
        writer->options.write_bloom_filters,
        writer->options.write_page_index,
        writer->options.write_statistics,
        writer->options.write_crc,
        writer->options.compression_level,
        writer->options.dictionary_page_size);

    /* Add all columns to the row group writer, resolving any per-column
       encoding/compression overrides into explicit values */
    for (int32_t i = 0; i < writer->num_columns; i++) {
        writer_column_def_t* col = &writer->columns[i];
        carquet_compression_t col_comp = effective_column_compression(writer, i);
        int32_t col_level = effective_column_compression_level(writer, i);
        carquet_encoding_t col_enc = effective_column_encoding(
            writer, i, col->physical_type, col_comp);
        /* When the leaf has no defined min/max sort order we pass a synthetic
         * logical type so the page writer suppresses min/max. GEOMETRY and
         * GEOGRAPHY are the exception: their real logical type must reach the
         * page writer so it accumulates GeospatialStatistics (min/max is
         * still suppressed for them by stats_order_defined_for_logical). */
        bool is_geo =
            (col->logical_type.id == CARQUET_LOGICAL_GEOMETRY ||
             col->logical_type.id == CARQUET_LOGICAL_GEOGRAPHY);
        carquet_logical_type_t no_stats_logical = { .id = CARQUET_LOGICAL_VARIANT };
        const carquet_logical_type_t* writer_logical_type =
            (col->statistics_sort_order_defined || is_geo)
                ? &col->logical_type : &no_stats_logical;

        carquet_status_t status = carquet_row_group_writer_add_column(
            rg,
            col->name,
            col->physical_type,
            writer_logical_type,
            col->max_def_level,
            col->max_rep_level,
            col->type_length,
            col_enc,
            col_comp,
            col_level);

        if (status != CARQUET_OK) {
            carquet_row_group_writer_destroy(rg);
            return status;
        }

        /* Apply the global max_rows_per_page knob (no-op when 0). */
        if (writer->options.max_rows_per_page > 0) {
            carquet_row_group_writer_set_column_max_rows_per_page(
                rg, i,
                writer->options.max_rows_per_page);
        }

        /* Apply the global write_batch_size knob (no-op when 0). */
        if (writer->options.write_batch_size > 0) {
            carquet_row_group_writer_set_column_write_batch_size(
                rg, i,
                writer->options.write_batch_size);
        }

        /* Apply per-column page-size override if set; otherwise the column
         * inherits the row-group default derived from options.page_size. */
        if (writer->column_overrides_allocated &&
            writer->column_page_size_override_set &&
            writer->column_page_size_override_set[i]) {
            carquet_row_group_writer_set_column_page_size(
                rg, i,
                writer->column_page_size_overrides[i]);
        }

        /* Opt-in Data Page V2 output. */
        if (writer->options.data_page_version == 2) {
            carquet_row_group_writer_set_column_data_page_v2(
                rg, i, true);
        }

        /* Apply per-column bloom NDV/FPP overrides. Once the new options API
         * has been used for ANY column, take full per-column control so a
         * column the user did not opt in does not silently get the default
         * bloom filter that the (now-on) global flag would create. When the
         * new API was never used this whole block is skipped, keeping the
         * legacy global-flag behavior byte-identical. */
        if (writer->column_overrides_allocated &&
            writer->column_bloom_options_set &&
            any_column_bloom_options_set(writer)) {
            /* A column keeps its bloom filter if it opted in through the newer
               ndv/fpp options API, OR was explicitly enabled/disabled through
               the legacy per-column setter. Columns that are true only because
               the global flag got flipped on (by the options API's side effect)
               are still suppressed, so they do not silently gain a default
               bloom the user never asked for. */
            bool col_enabled;
            if (writer->column_bloom_options_set[i] || writer->column_bloom_explicit[i]) {
                col_enabled = writer->column_bloom_filter_overrides[i];
            } else {
                col_enabled = false;
            }
            int64_t col_ndv = writer->column_bloom_options_set[i]
                ? writer->column_bloom_ndv_overrides[i] : 0;
            double col_fpp = writer->column_bloom_options_set[i]
                ? writer->column_bloom_fpp_overrides[i] : 0.0;
            carquet_row_group_writer_configure_column_bloom(
                rg, i,
                col_enabled, col_ndv, col_fpp);
        }
    }

    *out = rg;
    return CARQUET_OK;
}

/* Forget the rows of the row group just handed off; the next one is empty. */
static void reset_current_row_group_counters(carquet_writer_t* writer) {
    writer->current_row_group_rows = 0;
    writer->current_row_group_estimated_bytes = 0;
    for (int32_t i = 0; i < writer->num_columns; i++) {
        writer->column_values_written[i] = 0;
    }
}

static carquet_status_t ensure_row_group(carquet_writer_t* writer) {
    if (writer->current_row_group) {
        return CARQUET_OK;
    }
    carquet_status_t status = create_row_group_writer(writer, &writer->current_row_group);
    if (status != CARQUET_OK) {
        return status;
    }
    reset_current_row_group_counters(writer);
    return CARQUET_OK;
}

/* ============================================================================
 * Background I/O
 * ============================================================================
 */

static void io_write_task(void* arg) {
    io_job_t* job = (io_job_t*)arg;
    job->status = CARQUET_OK;
    if (job->buf.size > 0 &&
        fwrite(job->buf.data, 1, job->buf.size, job->file) != job->buf.size) {
        job->status = CARQUET_ERROR_FILE_WRITE;
        return;
    }
    for (int32_t k = 0; k < job->num_pages; k++) {
        const io_page_t* p = &job->pages[k];
        if (p->size > 0 && fwrite(p->data, 1, p->size, job->file) != p->size) {
            job->status = CARQUET_ERROR_FILE_WRITE;
            return;
        }
    }
}

/* Wait for every queued background write and take the chunk buffers back
 * for reuse. Returns the first write error, which also sticks on the writer.
 * Must precede any direct access to writer->file. */
static carquet_status_t writer_io_drain(carquet_writer_t* writer) {
    if (!writer->io_pool || writer->io_jobs_in_flight == 0) {
        return writer->io_error;
    }
    carquet_worker_pool_wait(writer->io_pool);
    for (int32_t i = 0; i < writer->io_jobs_in_flight; i++) {
        io_job_t* job = &writer->io_jobs[i];
        if (job->status != CARQUET_OK && writer->io_error == CARQUET_OK) {
            writer->io_error = job->status;
        }
        if (job->buf.owns_data && job->buf.data &&
            writer->io_spare_count < writer->num_columns) {
            writer->io_spares[writer->io_spare_count++] = job->buf;
        } else {
            carquet_buffer_destroy(&job->buf);
        }
        carquet_buffer_init(&job->buf);
        job->num_pages = 0;
    }
    writer->io_jobs_in_flight = 0;
    return writer->io_error;
}

static void writer_io_shutdown(carquet_writer_t* writer) {
    if (writer->io_pool) {
        /* Destroy runs the queue dry before joining, so nothing is lost. */
        (void)writer_io_drain(writer);
        carquet_worker_pool_destroy(writer->io_pool);
        writer->io_pool = NULL;
    }
    if (writer->io_jobs) {
        for (int32_t i = 0; i < writer->io_jobs_in_flight; i++) {
            carquet_buffer_destroy(&writer->io_jobs[i].buf);
        }
        for (int32_t i = 0; i < writer->num_columns; i++) {
            carquet_mem_free(writer->io_jobs[i].pages);
        }
        carquet_mem_free(writer->io_jobs);
        writer->io_jobs = NULL;
        writer->io_jobs_in_flight = 0;
    }
    if (writer->io_spares) {
        for (int32_t i = 0; i < writer->io_spare_count; i++) {
            carquet_buffer_destroy(&writer->io_spares[i]);
        }
        carquet_mem_free(writer->io_spares);
        writer->io_spares = NULL;
        writer->io_spare_count = 0;
    }
}

/* Start the I/O thread on first use. A failure here just means the row
 * group is written synchronously. */
static bool writer_io_ready(carquet_writer_t* writer) {
    if (!writer->options.async_io || writer->num_columns <= 0) return false;
    if (writer->io_pool) return true;
    writer->io_jobs = carquet_mem_calloc((size_t)writer->num_columns, sizeof(io_job_t));
    writer->io_spares = carquet_mem_calloc((size_t)writer->num_columns, sizeof(carquet_buffer_t));
    if (!writer->io_jobs || !writer->io_spares) {
        carquet_mem_free(writer->io_jobs);
        carquet_mem_free(writer->io_spares);
        writer->io_jobs = NULL;
        writer->io_spares = NULL;
        writer->options.async_io = false;
        return false;
    }
    for (int32_t i = 0; i < writer->num_columns; i++) {
        carquet_buffer_init(&writer->io_jobs[i].buf);
    }
    writer->io_pool = carquet_worker_pool_create(1);
    if (!writer->io_pool) {
        carquet_mem_free(writer->io_jobs);
        carquet_mem_free(writer->io_spares);
        writer->io_jobs = NULL;
        writer->io_spares = NULL;
        writer->options.async_io = false;
        return false;
    }
    return true;
}

/* Queue every column chunk of a finalized row group for the I/O thread. The
 * previous row group's jobs must have been drained (writer_io_drain). Pages
 * encoded in parallel are written straight from the row-group writer's page
 * buffers (page_io), which spares a copy of the whole row group. */
static carquet_status_t submit_row_group_io(carquet_writer_t* writer,
                                            carquet_row_group_writer_t* rg) {
    int num_cols = carquet_row_group_writer_num_columns(rg);
    for (int i = 0; i < num_cols && i < writer->num_columns; i++) {
        io_job_t* job = &writer->io_jobs[i];
        carquet_buffer_t spare;
        if (writer->io_spare_count > 0) {
            spare = writer->io_spares[--writer->io_spare_count];
        } else {
            carquet_buffer_init(&spare);
        }
        carquet_row_group_writer_detach_column(rg, i, &job->buf, &spare);
        int32_t num_pages = 0;
        const carquet_buffer_t* pages = carquet_row_group_writer_column_pages(
            rg, i, &num_pages);
        if (num_pages > job->page_capacity) {
            io_page_t* grown = carquet_mem_realloc(
                job->pages, (size_t)num_pages * sizeof(*grown));
            if (!grown) {
                /* Jobs 0..i-1 are already queued; the drain at the next
                 * flush or close reclaims them. */
                writer->io_error = CARQUET_ERROR_OUT_OF_MEMORY;
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            job->pages = grown;
            job->page_capacity = num_pages;
        }
        for (int32_t k = 0; k < num_pages; k++) {
            job->pages[k].data = pages[k].data;
            job->pages[k].size = pages[k].size;
        }
        job->num_pages = num_pages;
        job->file = writer->file;
        job->status = CARQUET_OK;
        writer->io_jobs_in_flight = i + 1;
        carquet_worker_pool_submit(writer->io_pool, io_write_task, job);
    }
    return CARQUET_OK;
}

/* Finalize the current row group and queue its column chunks for the I/O
 * thread. Encoding runs first, overlapping the previous row group's write;
 * only then do we wait for that write, reclaim its buffers and submit. */
static carquet_status_t write_row_group_async(carquet_writer_t* writer, size_t* size) {
    carquet_row_group_writer_set_page_io(writer->current_row_group, true);
    carquet_status_t status = carquet_row_group_writer_finalize_columns(
        writer->current_row_group, size, writer->current_row_group_rows);
    if (status != CARQUET_OK) {
        return status;
    }
    status = writer_io_drain(writer);
    if (status != CARQUET_OK) {
        return status;
    }
    return submit_row_group_io(writer, writer->current_row_group);
}

/* Append the metadata of a row group that has just been encoded and queued
 * (or written): its row_group_info_t with per-chunk offsets, sizes and
 * statistics, then any bloom filters and page indexes, which go to the file
 * right behind the row group's data. `size` is the row group's byte size. */
static carquet_status_t record_row_group(carquet_writer_t* writer,
                                         carquet_row_group_writer_t* rg,
                                         int64_t rows, size_t size) {
    carquet_status_t status;
    /* Store row group metadata */
    if (writer->num_row_groups >= writer->row_groups_capacity) {
        int32_t new_cap = writer->row_groups_capacity == 0 ? 4 : writer->row_groups_capacity * 2;
        row_group_info_t* new_rgs = carquet_mem_realloc(writer->row_groups,
            new_cap * sizeof(row_group_info_t));
        if (!new_rgs) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        writer->row_groups = new_rgs;
        writer->row_groups_capacity = new_cap;
    }

    row_group_info_t* rg_info = &writer->row_groups[writer->num_row_groups];
    memset(rg_info, 0, sizeof(*rg_info));

    rg_info->file_offset = writer->file_offset;
    rg_info->num_rows = rows;
    rg_info->total_byte_size = carquet_row_group_writer_total_byte_size(rg);
    rg_info->total_compressed_size = (int64_t)size;
    rg_info->ordinal = (int16_t)writer->num_row_groups;

    /* Build column chunks metadata */
    int num_cols = carquet_row_group_writer_num_columns(rg);
    rg_info->num_columns = num_cols;
    rg_info->columns = carquet_mem_calloc((size_t)num_cols, sizeof(*rg_info->columns));
    if (!rg_info->columns) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int i = 0; i < num_cols; i++) {
        const column_chunk_info_t* col_info = carquet_row_group_writer_get_column_info(
            rg, i);

        if (!col_info) continue;

        row_group_column_info_t* chunk = &rg_info->columns[i];
        chunk->file_offset = col_info->file_offset;
        chunk->type = col_info->type;
        chunk->codec = col_info->compression;
        chunk->num_values = col_info->num_values;
        chunk->total_compressed_size = col_info->total_compressed_size;
        chunk->total_uncompressed_size = col_info->total_uncompressed_size;
        if (col_info->has_dictionary_page) {
            /* The dictionary page is the first page of the chunk at
             * file_offset; data pages follow it. */
            chunk->has_dictionary_page = true;
            chunk->has_dictionary_page_offset = true;
            chunk->dictionary_page_offset = col_info->file_offset;
            chunk->data_page_offset =
                col_info->file_offset + col_info->dictionary_page_size;
        } else {
            chunk->data_page_offset = col_info->file_offset;
        }

        if (col_info->has_geo_stats) {
            chunk->has_geo_stats = true;
            chunk->geo_stats = col_info->geo_stats;
        }

        /* Copy SizeStatistics (Parquet 2.9) out of the column writer's aliased
         * histogram buffers into owned storage that lives with the file's row
         * group list. Independent of write_statistics (emitted below only when
         * it carries information). */
        chunk->unencoded_ba_bytes = col_info->unencoded_ba_bytes;
        if (col_info->rep_level_hist && col_info->rep_hist_len > 0) {
            chunk->rep_level_hist = carquet_mem_malloc(
                (size_t)col_info->rep_hist_len * sizeof(int64_t));
            if (chunk->rep_level_hist) {
                memcpy(chunk->rep_level_hist, col_info->rep_level_hist,
                       (size_t)col_info->rep_hist_len * sizeof(int64_t));
                chunk->rep_hist_len = col_info->rep_hist_len;
            }
        }
        if (col_info->def_level_hist && col_info->def_hist_len > 0) {
            chunk->def_level_hist = carquet_mem_malloc(
                (size_t)col_info->def_hist_len * sizeof(int64_t));
            if (chunk->def_level_hist) {
                memcpy(chunk->def_level_hist, col_info->def_level_hist,
                       (size_t)col_info->def_hist_len * sizeof(int64_t));
                chunk->def_hist_len = col_info->def_hist_len;
            }
        }

        /* Derive the per-chunk encodings list from whether this chunk actually
         * emitted a dictionary page. The static refresh_column_encodings_cache
         * is computed from the configured encoding before finalize knows about
         * dictionary fallback, so a dict column that fell back to PLAIN (dict
         * exceeded dictionary_page_size or was all-unique) must NOT advertise
         * RLE_DICTIONARY here. Per the Parquet spec, ColumnMetaData.encodings
         * is the set of encodings actually used in the chunk. */
        if (col_info->has_dictionary_page) {
            chunk->encodings[0] = CARQUET_ENCODING_PLAIN;
            chunk->encodings[1] = CARQUET_ENCODING_RLE_DICTIONARY;
            chunk->encodings[2] = CARQUET_ENCODING_RLE;
            chunk->num_encodings = 3;
        } else {
            /* Plain column, or dictionary column that fell back to PLAIN. For
             * the fallback case the data encoding is PLAIN; for an explicitly
             * non-dictionary column (e.g. BYTE_STREAM_SPLIT) it is the
             * configured encoding. A configured dictionary encoding that
             * produced no dictionary page implies a PLAIN fallback. */
            carquet_encoding_t data_enc = col_info->encoding;
            if (data_enc == CARQUET_ENCODING_RLE_DICTIONARY ||
                data_enc == CARQUET_ENCODING_PLAIN_DICTIONARY) {
                data_enc = CARQUET_ENCODING_PLAIN;
            }
            chunk->encodings[0] = data_enc;
            chunk->encodings[1] = CARQUET_ENCODING_RLE;
            chunk->num_encodings = 2;
        }

        /* Per-column statistics override may suppress emission for one column
         * even when global write_statistics is enabled. */
        bool emit_stats = writer->column_overrides_allocated && i < writer->num_columns
            ? writer->column_statistics_overrides[i]
            : writer->options.write_statistics;
        if (i < writer->num_columns && !writer->columns[i].statistics_sort_order_defined) {
            emit_stats = false;
        }

        /* SizeStatistics emission is gated on the write-statistics option but,
         * unlike min/max, NOT on sort-order definedness: the level histograms
         * are most valuable exactly for nested/repeated columns whose sort
         * order is undefined. Emit only when it carries information (BYTE_ARRAY
         * unencoded bytes, or a non-trivial rep/def histogram) so flat required
         * numeric columns stay byte-identical to before. */
        bool size_stats_enabled = writer->column_overrides_allocated && i < writer->num_columns
            ? writer->column_statistics_overrides[i]
            : writer->options.write_statistics;
        chunk->has_size_statistics = size_stats_enabled &&
            (chunk->type == CARQUET_PHYSICAL_BYTE_ARRAY ||
             chunk->rep_hist_len > 1 || chunk->def_hist_len > 1);

        if (emit_stats && (col_info->has_min_max || col_info->has_null_count)) {
            chunk->has_statistics = true;
            chunk->has_min_max = col_info->has_min_max;
            chunk->has_null_count = col_info->has_null_count;
            chunk->null_count = col_info->null_count;
            chunk->has_distinct_count = col_info->has_distinct_count;
            chunk->distinct_count = col_info->distinct_count;
            if (col_info->has_min_max) {
                /* zero-length = the empty BYTE_ARRAY, a real bound */
                chunk->min_value = carquet_mem_malloc(
                    col_info->min_value_size > 0 ? col_info->min_value_size : 1);
                chunk->max_value = carquet_mem_malloc(
                    col_info->max_value_size > 0 ? col_info->max_value_size : 1);
                if (chunk->min_value && chunk->max_value) {
                    if (col_info->min_value_size > 0)
                        memcpy(chunk->min_value, col_info->min_value,
                               col_info->min_value_size);
                    if (col_info->max_value_size > 0)
                        memcpy(chunk->max_value, col_info->max_value,
                               col_info->max_value_size);
                    chunk->min_value_size = (int32_t)col_info->min_value_size;
                    chunk->max_value_size = (int32_t)col_info->max_value_size;
                } else {
                    carquet_mem_free(chunk->min_value);
                    carquet_mem_free(chunk->max_value);
                    chunk->min_value = NULL;
                    chunk->max_value = NULL;
                    chunk->has_min_max = false;
                }
            }
        }
    }

    writer->num_row_groups++;
    writer->file_offset += (int64_t)size;
    writer->total_rows += rows;

    /* Write bloom filters for each column (after row group data) */
    if (writer->options.write_bloom_filters) {
        status = writer_io_drain(writer);
        if (status != CARQUET_OK) {
            return status;
        }
        for (int i = 0; i < num_cols; i++) {
            carquet_bloom_filter_t* bf = carquet_row_group_writer_get_bloom_filter(
                rg, i);
            if (!bf) continue;

            const uint8_t* bf_data = carquet_bloom_filter_data(bf);
            size_t bf_size = carquet_bloom_filter_size(bf);
            if (!bf_data || bf_size == 0) continue;

            /* Write Bloom Filter Header (Thrift):
             * numBytes: i32, algorithm: MURMUR3_X64_128, hash: XXHASH, compression: UNCOMPRESSED */
            carquet_buffer_t bf_header;
            carquet_buffer_init(&bf_header);
            {
                thrift_encoder_t enc;
                thrift_encoder_init(&enc, &bf_header);
                thrift_write_struct_begin(&enc);
                /* Field 1: numBytes (i32) */
                thrift_write_field_header(&enc, THRIFT_TYPE_I32, 1);
                thrift_write_i32(&enc, (int32_t)bf_size);
                /* Field 2: algorithm (BloomFilterAlgorithm struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 2);
                thrift_write_struct_begin(&enc);
                /* Field 1: SPLIT_BLOCK_BLOOM_FILTER (empty struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 1);
                thrift_write_struct_begin(&enc);
                thrift_write_struct_end(&enc);
                thrift_write_struct_end(&enc);
                /* Field 3: hash (BloomFilterHash struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 3);
                thrift_write_struct_begin(&enc);
                /* Field 1: XXHASH (empty struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 1);
                thrift_write_struct_begin(&enc);
                thrift_write_struct_end(&enc);
                thrift_write_struct_end(&enc);
                /* Field 4: compression (BloomFilterCompression struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 4);
                thrift_write_struct_begin(&enc);
                /* Field 1: UNCOMPRESSED (empty struct) */
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 1);
                thrift_write_struct_begin(&enc);
                thrift_write_struct_end(&enc);
                thrift_write_struct_end(&enc);
                thrift_write_struct_end(&enc);
            }

            /* Record offset in column metadata */
            row_group_column_info_t* chunk = &rg_info->columns[i];
            chunk->has_bloom_filter_offset = true;
            chunk->bloom_filter_offset = writer->file_offset;
            chunk->has_bloom_filter_length = true;
            chunk->bloom_filter_length = (int32_t)(bf_header.size + bf_size);

            /* Write header + data */
            if (fwrite(bf_header.data, 1, bf_header.size, writer->file) != bf_header.size) {
                carquet_buffer_destroy(&bf_header);
                return CARQUET_ERROR_FILE_WRITE;
            }
            writer->file_offset += (int64_t)bf_header.size;
            carquet_buffer_destroy(&bf_header);

            if (fwrite(bf_data, 1, bf_size, writer->file) != bf_size) {
                return CARQUET_ERROR_FILE_WRITE;
            }
            writer->file_offset += (int64_t)bf_size;
        }
    }

    /* Write column indexes and offset indexes (after bloom filters) */
    if (writer->options.write_page_index) {
        status = writer_io_drain(writer);
        if (status != CARQUET_OK) {
            return status;
        }
        for (int i = 0; i < num_cols; i++) {
            row_group_column_info_t* chunk = &rg_info->columns[i];

            /* Column index */
            carquet_column_index_builder_t* ci = carquet_row_group_writer_get_column_index(
                rg, i);
            if (ci) {
                carquet_buffer_t ci_buf;
                carquet_buffer_init(&ci_buf);
                carquet_column_index_serialize(ci, &ci_buf);
                if (ci_buf.size > 0) {
                    chunk->has_column_index_offset = true;
                    chunk->column_index_offset = writer->file_offset;
                    chunk->has_column_index_length = true;
                    chunk->column_index_length = (int32_t)ci_buf.size;
                    if (fwrite(ci_buf.data, 1, ci_buf.size, writer->file) != ci_buf.size) {
                        carquet_buffer_destroy(&ci_buf);
                        return CARQUET_ERROR_FILE_WRITE;
                    }
                    writer->file_offset += (int64_t)ci_buf.size;
                }
                carquet_buffer_destroy(&ci_buf);
            }

            /* Offset index */
            carquet_offset_index_builder_t* oi = carquet_row_group_writer_get_offset_index(
                rg, i);
            if (oi) {
                carquet_buffer_t oi_buf;
                carquet_buffer_init(&oi_buf);
                carquet_offset_index_serialize(oi, &oi_buf);
                if (oi_buf.size > 0) {
                    chunk->has_offset_index_offset = true;
                    chunk->offset_index_offset = writer->file_offset;
                    chunk->has_offset_index_length = true;
                    chunk->offset_index_length = (int32_t)oi_buf.size;
                    if (fwrite(oi_buf.data, 1, oi_buf.size, writer->file) != oi_buf.size) {
                        carquet_buffer_destroy(&oi_buf);
                        return CARQUET_ERROR_FILE_WRITE;
                    }
                    writer->file_offset += (int64_t)oi_buf.size;
                }
                carquet_buffer_destroy(&oi_buf);
            }
        }
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Pipelined encode
 * ============================================================================
 *
 * The fork/join finalize (finalize_columns_parallel) ends every row group
 * with a barrier: its wall time is the slowest page of the last wave, on a
 * loaded machine a page whose thread was descheduled, and the main thread
 * sits out the stash of the next row group and the I/O hand-off serially.
 * When every column is page-parallel, the row group is instead handed to a
 * persistent pool as soon as the user starts the next one: its pages encode
 * while the user fills the other row-group writer, and the stragglers of row
 * group n overlap the first pages of n+1. The row group is completed (pages
 * accounted, chunk infos filled, metadata recorded, I/O queued) one flush
 * later, or at close. Output is the same bytes in the same order; only
 * which thread encodes which page changes.
 */

#ifdef _OPENMP

/* Named PIPELINE_*, not PIPE_*: <winbase.h> (pulled in through
 * worker_pool.h) already defines PIPE_WAIT for named pipes, and the
 * redefinition is a C4005 on MSVC / -Wall warning on MinGW. */
#ifdef _WIN32
#define PIPELINE_LOCK(w)      EnterCriticalSection(&(w)->pipe_sync.mutex)
#define PIPELINE_UNLOCK(w)    LeaveCriticalSection(&(w)->pipe_sync.mutex)
#define PIPELINE_WAIT(w)      SleepConditionVariableCS(&(w)->pipe_sync.cond, &(w)->pipe_sync.mutex, INFINITE)
#define PIPELINE_BROADCAST(w) WakeAllConditionVariable(&(w)->pipe_sync.cond)
#else
#define PIPELINE_LOCK(w)      pthread_mutex_lock(&(w)->pipe_sync.mutex)
#define PIPELINE_UNLOCK(w)    pthread_mutex_unlock(&(w)->pipe_sync.mutex)
#define PIPELINE_WAIT(w)      pthread_cond_wait(&(w)->pipe_sync.cond, &(w)->pipe_sync.mutex)
#define PIPELINE_BROADCAST(w) pthread_cond_broadcast(&(w)->pipe_sync.cond)
#endif

static bool pipeline_sync_init(carquet_writer_t* writer) {
#ifdef _WIN32
    InitializeCriticalSection(&writer->pipe_sync.mutex);
    InitializeConditionVariable(&writer->pipe_sync.cond);
#else
    if (pthread_mutex_init(&writer->pipe_sync.mutex, NULL) != 0) return false;
    if (pthread_cond_init(&writer->pipe_sync.cond, NULL) != 0) {
        pthread_mutex_destroy(&writer->pipe_sync.mutex);
        return false;
    }
#endif
    writer->pipe_sync_ready = true;
    return true;
}

static void pipeline_sync_destroy(carquet_writer_t* writer) {
    if (!writer->pipe_sync_ready) return;
#ifdef _WIN32
    DeleteCriticalSection(&writer->pipe_sync.mutex);
#else
    pthread_mutex_destroy(&writer->pipe_sync.mutex);
    pthread_cond_destroy(&writer->pipe_sync.cond);
#endif
    writer->pipe_sync_ready = false;
}

/* Scratch page writers: one per pool thread at most, taken for the duration
 * of a task. Created on demand, kept across row groups so their buffers
 * stay warm. */
static carquet_page_writer_t* pipeline_take_scratch(carquet_writer_t* writer) {
    carquet_page_writer_t* scratch = NULL;
    PIPELINE_LOCK(writer);
    if (writer->scratch_count > 0) {
        scratch = writer->scratch_stack[--writer->scratch_count];
    }
    PIPELINE_UNLOCK(writer);
    if (!scratch) scratch = carquet_page_writer_create_scratch();
    return scratch;
}

static void pipeline_return_scratch(carquet_writer_t* writer,
                                    carquet_page_writer_t* scratch) {
    if (!scratch) return;
    PIPELINE_LOCK(writer);
    if (writer->scratch_count < writer->scratch_capacity) {
        writer->scratch_stack[writer->scratch_count++] = scratch;
        scratch = NULL;
    }
    PIPELINE_UNLOCK(writer);
    if (scratch) carquet_page_writer_destroy(scratch);
}

static void encode_task_fn(void* arg) {
    encode_task_ctx_t* ctx = (encode_task_ctx_t*)arg;
    carquet_writer_t* writer = ctx->writer;
    carquet_page_writer_t* scratch = pipeline_take_scratch(writer);
    double t0 = omp_get_wtime();
    if (scratch) {
        carquet_column_writer_run_page_task(ctx->task, scratch);
    } else {
        ctx->task->status = CARQUET_ERROR_OUT_OF_MEMORY;
    }
    ctx->task->elapsed = omp_get_wtime() - t0;
    pipeline_return_scratch(writer, scratch);
    PIPELINE_LOCK(writer);
    if (--ctx->rg->remaining == 0) {
        PIPELINE_BROADCAST(writer);
    }
    PIPELINE_UNLOCK(writer);
}

/* Block until every task of the slot has run. */
static void pipeline_wait(carquet_writer_t* writer, pipeline_rg_t* slot) {
    PIPELINE_LOCK(writer);
    while (slot->remaining > 0) {
        PIPELINE_WAIT(writer);
    }
    PIPELINE_UNLOCK(writer);
}

/* Set the pipeline up on first use. Needs the background I/O thread (the
 * completed row group is handed to it) and a machine worth threading. A
 * failure here just leaves the fork/join path in charge. */
static bool pipeline_ready(carquet_writer_t* writer) {
    if (writer->pipeline_off || writer->num_columns <= 1 ||
        writer->options.write_bloom_filters || writer->options.write_page_index) {
        return false;
    }
    if (!writer_io_ready(writer)) return false;
    if (writer->encode_pool) return true;

    int32_t threads = omp_get_max_threads();
    if (threads < 2) {
        writer->pipeline_off = true;
        return false;
    }
    if (!pipeline_sync_init(writer)) {
        writer->pipeline_off = true;
        return false;
    }
    writer->scratch_stack = carquet_mem_calloc((size_t)threads, sizeof(*writer->scratch_stack));
    writer->scratch_capacity = threads;
    writer->scratch_count = 0;
    writer->pipe_pending = -1;
    carquet_status_t status = writer->scratch_stack
        ? create_row_group_writer(writer, &writer->spare_row_group)
        : CARQUET_ERROR_OUT_OF_MEMORY;
    if (status == CARQUET_OK) {
        writer->encode_pool = carquet_worker_pool_create(threads);
        if (!writer->encode_pool) status = CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (status != CARQUET_OK) {
        if (writer->spare_row_group) {
            carquet_row_group_writer_destroy(writer->spare_row_group);
            writer->spare_row_group = NULL;
        }
        carquet_mem_free(writer->scratch_stack);
        writer->scratch_stack = NULL;
        writer->scratch_capacity = 0;
        pipeline_sync_destroy(writer);
        writer->pipeline_off = true;
        return false;
    }
    writer->encode_threads = threads;
    return true;
}

/* Start encoding the current row group on the pool, if every column allows
 * it. `*started` tells; when false nothing happened and the caller finalizes
 * the row group the ordinary way. */
static carquet_status_t pipeline_begin(carquet_writer_t* writer, bool* started,
                                       int* slot_index) {
    *started = false;
    int free_slot = (writer->pipe_pending == 0) ? 1 : 0;
    pipeline_rg_t* slot = &writer->pipe_slots[free_slot];
    if (slot->rg) return CARQUET_ERROR_INTERNAL;  /* both slots busy: cannot happen */

    carquet_row_group_writer_set_page_io(writer->current_row_group, true);
    int32_t n = carquet_row_group_writer_begin_page_tasks(writer->current_row_group);
    if (n <= 0) return CARQUET_OK;

    if (n > slot->ctx_capacity) {
        encode_task_ctx_t* grown = carquet_mem_realloc(
            slot->ctx, (size_t)n * sizeof(*grown));
        if (!grown) return CARQUET_ERROR_OUT_OF_MEMORY;
        slot->ctx = grown;
        slot->ctx_capacity = n;
    }
    for (int32_t k = 0; k < n; k++) {
        carquet_page_task_t* task = carquet_row_group_writer_page_task(
            writer->current_row_group, k);
        if (!task) return CARQUET_ERROR_INTERNAL;
        slot->ctx[k].writer = writer;
        slot->ctx[k].rg = slot;
        slot->ctx[k].task = task;
    }
    slot->rg = writer->current_row_group;
    slot->rows = writer->current_row_group_rows;
    PIPELINE_LOCK(writer);
    slot->remaining = n;
    PIPELINE_UNLOCK(writer);
    for (int32_t k = 0; k < n; k++) {
        carquet_worker_pool_submit(writer->encode_pool, encode_task_fn, &slot->ctx[k]);
    }
    *started = true;
    *slot_index = free_slot;
    return CARQUET_OK;
}

/* Finish the row group that is in the pipeline, if any: wait for its pages,
 * queue them for the I/O thread behind the previous row group, record its
 * metadata and return its writer to spare duty. On failure the writer keeps
 * the slot (its buffers may still be referenced) and the error sticks. */
static carquet_status_t pipeline_complete_pending(carquet_writer_t* writer) {
    if (writer->pipe_pending < 0) return CARQUET_OK;
    pipeline_rg_t* slot = &writer->pipe_slots[writer->pipe_pending];
    if (!slot->rg) {
        writer->pipe_pending = -1;
        return CARQUET_OK;
    }
    pipeline_wait(writer, slot);

    size_t size = 0;
    carquet_status_t status = writer_io_drain(writer);
    if (status == CARQUET_OK) {
        status = carquet_row_group_writer_complete_page_tasks(
            slot->rg, writer->file_offset, slot->rows, &size);
    }
    if (status == CARQUET_OK) {
        status = submit_row_group_io(writer, slot->rg);
    }
    if (status == CARQUET_OK) {
        status = record_row_group(writer, slot->rg, slot->rows, size);
    }
    if (status != CARQUET_OK) {
        if (writer->io_error == CARQUET_OK) writer->io_error = status;
        return status;
    }
    carquet_row_group_writer_reset(slot->rg, writer->file_offset);
    if (writer->spare_row_group) {
        carquet_row_group_writer_destroy(writer->spare_row_group);
    }
    writer->spare_row_group = slot->rg;
    slot->rg = NULL;
    writer->pipe_pending = -1;
    return CARQUET_OK;
}

/* Stop the encode pool: every queued task runs to the end first, so no
 * thread touches a row-group writer after this. The I/O thread may still be
 * reading page buffers; pipeline_free comes after writer_io_shutdown. */
static void pipeline_quiesce(carquet_writer_t* writer) {
    if (writer->encode_pool) {
        carquet_worker_pool_destroy(writer->encode_pool);
        writer->encode_pool = NULL;
    }
}

static void pipeline_free(carquet_writer_t* writer) {
    pipeline_quiesce(writer);
    for (int i = 0; i < 2; i++) {
        pipeline_rg_t* slot = &writer->pipe_slots[i];
        if (slot->rg) {
            carquet_row_group_writer_destroy(slot->rg);
            slot->rg = NULL;
        }
        carquet_mem_free(slot->ctx);
        slot->ctx = NULL;
        slot->ctx_capacity = 0;
    }
    writer->pipe_pending = -1;
    if (writer->spare_row_group) {
        carquet_row_group_writer_destroy(writer->spare_row_group);
        writer->spare_row_group = NULL;
    }
    if (writer->scratch_stack) {
        for (int32_t i = 0; i < writer->scratch_count; i++) {
            carquet_page_writer_destroy(writer->scratch_stack[i]);
        }
        carquet_mem_free(writer->scratch_stack);
        writer->scratch_stack = NULL;
        writer->scratch_count = 0;
        writer->scratch_capacity = 0;
    }
    pipeline_sync_destroy(writer);
}

#else  /* !_OPENMP: no page-parallel finalize, so nothing to pipeline */

static bool pipeline_ready(carquet_writer_t* writer) { (void)writer; return false; }
static carquet_status_t pipeline_begin(carquet_writer_t* writer, bool* started,
                                       int* slot_index) {
    (void)writer; (void)slot_index; *started = false; return CARQUET_OK;
}
static carquet_status_t pipeline_complete_pending(carquet_writer_t* writer) {
    (void)writer; return CARQUET_OK;
}
static void pipeline_quiesce(carquet_writer_t* writer) { (void)writer; }
static void pipeline_free(carquet_writer_t* writer) { (void)writer; }

#endif

static carquet_status_t flush_row_group(carquet_writer_t* writer) {
    if (!writer->current_row_group || writer->current_row_group_rows == 0) {
        return CARQUET_OK;
    }
    if (writer->io_error != CARQUET_OK) {
        return writer->io_error;
    }

    carquet_status_t status;
    if (pipeline_ready(writer)) {
        /* Start this row group's pages first, then settle the previous row
         * group: its last pages finish while these begin. */
        bool started = false;
        int slot = -1;
        status = pipeline_begin(writer, &started, &slot);
        if (status != CARQUET_OK) {
            return status;
        }
        status = pipeline_complete_pending(writer);
        if (status != CARQUET_OK) {
            return status;
        }
        if (started) {
            writer->pipe_pending = slot;
            writer->current_row_group = writer->spare_row_group;
            writer->spare_row_group = NULL;
            if (!writer->current_row_group) {
                status = create_row_group_writer(writer, &writer->current_row_group);
                if (status != CARQUET_OK) {
                    return status;
                }
            }
            reset_current_row_group_counters(writer);
            return CARQUET_OK;
        }
        /* Not page-parallel in full: the ordinary path below, in order. */
    } else {
        status = pipeline_complete_pending(writer);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    /* Finalize and write each column directly to file, avoiding
     * an intermediate copy of the entire row group into one buffer */
    size_t size;
    if (writer_io_ready(writer)) {
        status = write_row_group_async(writer, &size);
    } else {
        carquet_row_group_writer_set_page_io(writer->current_row_group, false);
        status = carquet_row_group_writer_write_to_file(
            writer->current_row_group, writer->file, &size,
            writer->current_row_group_rows);
    }
    if (status != CARQUET_OK) {
        return status;
    }

    status = record_row_group(writer, writer->current_row_group,
                              writer->current_row_group_rows, size);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Reuse the current row group writer for the next group. */
    carquet_row_group_writer_reset(writer->current_row_group, writer->file_offset);
    reset_current_row_group_counters(writer);
    return CARQUET_OK;
}

static bool deprecated_stats_are_compatible(const writer_column_def_t* column) {
    if (!column) return false;
    if (column->logical_type.id == CARQUET_LOGICAL_INTEGER &&
        !column->logical_type.params.integer.is_signed) {
        return false;
    }
    if (!column->statistics_sort_order_defined) {
        return false;
    }
    return true;
}

static carquet_status_t build_file_metadata(
    carquet_writer_t* writer,
    parquet_file_metadata_t* metadata) {

    memset(metadata, 0, sizeof(*metadata));

    metadata->version = (writer->options.file_format_version == 1) ? 1 : 2;
    metadata->num_rows = writer->total_rows;
    metadata->num_column_orders = writer->num_columns;
    metadata->created_by = carquet_arena_strdup(&writer->arena,
        writer->options.created_by ? writer->options.created_by : "Carquet");

    /* Optional "ARROW:schema" footer metadata. Skipped if the user already
     * supplied that key, if the schema is nested/unsupported, or on OOM. */
    char* arrow_schema_b64 = NULL;
    if (writer->options.write_arrow_schema) {
        bool user_set = false;
        for (int32_t i = 0; i < writer->num_kv_metadata; i++) {
            if (writer->kv_metadata[i].key &&
                strcmp(writer->kv_metadata[i].key, "ARROW:schema") == 0) {
                user_set = true;
                break;
            }
        }
        if (!user_set) {
            arrow_schema_b64 = carquet_build_arrow_schema_b64(
                writer->schema_elements, writer->num_schema_elements);
        }
    }
    int32_t extra_kv = arrow_schema_b64 ? 1 : 0;

    /* Key-value metadata */
    if (writer->num_kv_metadata + extra_kv > 0) {
        metadata->num_key_value = writer->num_kv_metadata + extra_kv;
        metadata->key_value_metadata = carquet_arena_calloc(&writer->arena,
            metadata->num_key_value, sizeof(parquet_key_value_t));
        if (!metadata->key_value_metadata) {
            carquet_mem_free(arrow_schema_b64);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        for (int32_t i = 0; i < writer->num_kv_metadata; i++) {
            metadata->key_value_metadata[i].key = carquet_arena_strdup(
                &writer->arena, writer->kv_metadata[i].key);
            metadata->key_value_metadata[i].value = writer->kv_metadata[i].value ?
                carquet_arena_strdup(&writer->arena, writer->kv_metadata[i].value) : NULL;
        }
        if (arrow_schema_b64) {
            parquet_key_value_t* kv =
                &metadata->key_value_metadata[writer->num_kv_metadata];
            kv->key = carquet_arena_strdup(&writer->arena, "ARROW:schema");
            kv->value = carquet_arena_strdup(&writer->arena, arrow_schema_b64);
            carquet_mem_free(arrow_schema_b64);
        }
    }

    /* Build schema from stored elements (includes groups for nested schemas) */
    metadata->num_schema_elements = writer->num_schema_elements;
    metadata->schema = carquet_arena_calloc(&writer->arena, writer->num_schema_elements,
        sizeof(parquet_schema_element_t));

    if (!metadata->schema) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < writer->num_schema_elements; i++) {
        metadata->schema[i] = writer->schema_elements[i];
        /* Duplicate strings into arena so they outlive the writer */
        if (writer->schema_elements[i].name) {
            metadata->schema[i].name = carquet_arena_strdup(
                &writer->arena, writer->schema_elements[i].name);
        }
    }

    /* Row groups */
    metadata->num_row_groups = writer->num_row_groups;
    metadata->row_groups = carquet_arena_calloc(&writer->arena, writer->num_row_groups,
        sizeof(parquet_row_group_t));

    if (!metadata->row_groups && writer->num_row_groups > 0) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < writer->num_row_groups; i++) {
        const row_group_info_t* src_rg = &writer->row_groups[i];
        parquet_row_group_t* dst_rg = &metadata->row_groups[i];

        dst_rg->num_rows = src_rg->num_rows;
        dst_rg->total_byte_size = src_rg->total_byte_size;
        dst_rg->has_file_offset = true;
        dst_rg->file_offset = src_rg->file_offset;
        dst_rg->has_total_compressed_size = true;
        dst_rg->total_compressed_size = src_rg->total_compressed_size;
        dst_rg->has_ordinal = true;
        dst_rg->ordinal = src_rg->ordinal;
        dst_rg->sorting_columns = writer->sorting_columns;
        dst_rg->num_sorting_columns = writer->num_sorting_columns;
        dst_rg->num_columns = src_rg->num_columns;
        dst_rg->columns = carquet_arena_calloc(&writer->arena, src_rg->num_columns,
            sizeof(parquet_column_chunk_t));
        if (!dst_rg->columns && src_rg->num_columns > 0) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        for (int32_t j = 0; j < src_rg->num_columns; j++) {
            const row_group_column_info_t* src_col = &src_rg->columns[j];
            parquet_column_chunk_t* dst_chunk = &dst_rg->columns[j];
            parquet_column_metadata_t* meta = &dst_chunk->metadata;

            dst_chunk->file_offset = src_col->file_offset;
            dst_chunk->has_metadata = true;

            meta->type = src_col->type;
            meta->codec = src_col->codec;
            meta->num_values = src_col->num_values;
            meta->total_compressed_size = src_col->total_compressed_size;
            meta->total_uncompressed_size = src_col->total_uncompressed_size;
            meta->data_page_offset = src_col->data_page_offset;
            if (src_col->has_dictionary_page_offset) {
                meta->has_dictionary_page_offset = true;
                meta->dictionary_page_offset = src_col->dictionary_page_offset;
            }
            /* Use the per-chunk encodings derived post-finalize from the
             * actual dictionary-page presence rather than the static
             * configured-encoding cache, so a dict-to-PLAIN fallback chunk
             * does not falsely advertise RLE_DICTIONARY. Fall back to the
             * static cache only if the per-chunk list was never populated. */
            if (src_col->num_encodings > 0) {
                meta->num_encodings = src_col->num_encodings;
                meta->encodings = (carquet_encoding_t*)src_col->encodings;
            } else {
                meta->num_encodings = writer->column_num_encodings[j] > 0
                    ? writer->column_num_encodings[j] : 2;
                meta->encodings = writer->column_encodings[j];
            }
            meta->path_len = writer->column_path_lens[j];
            meta->path_in_schema = writer->column_paths[j];

            meta->has_bloom_filter_offset = src_col->has_bloom_filter_offset;
            meta->bloom_filter_offset = src_col->bloom_filter_offset;
            meta->has_bloom_filter_length = src_col->has_bloom_filter_length;
            meta->bloom_filter_length = src_col->bloom_filter_length;

            dst_chunk->has_column_index_offset = src_col->has_column_index_offset;
            dst_chunk->column_index_offset = src_col->column_index_offset;
            dst_chunk->has_column_index_length = src_col->has_column_index_length;
            dst_chunk->column_index_length = src_col->column_index_length;
            dst_chunk->has_offset_index_offset = src_col->has_offset_index_offset;
            dst_chunk->offset_index_offset = src_col->offset_index_offset;
            dst_chunk->has_offset_index_length = src_col->has_offset_index_length;
            dst_chunk->offset_index_length = src_col->offset_index_length;

            if (src_col->has_geo_stats &&
                (src_col->geo_stats.valid || src_col->geo_stats.num_types > 0)) {
                meta->has_geospatial_statistics = true;
                meta->geospatial_statistics = src_col->geo_stats;
            }

            if (src_col->has_size_statistics) {
                meta->has_size_statistics = true;
                parquet_size_statistics_t* ss = &meta->size_statistics;
                memset(ss, 0, sizeof(*ss));
                if (src_col->type == CARQUET_PHYSICAL_BYTE_ARRAY &&
                    src_col->unencoded_ba_bytes >= 0) {
                    ss->has_unencoded_byte_array_data_bytes = true;
                    ss->unencoded_byte_array_data_bytes = src_col->unencoded_ba_bytes;
                }
                /* Emit histograms only when they carry information: rep for
                 * repeated columns, def for nullable/nested ones. */
                if (src_col->rep_hist_len > 1) {
                    ss->repetition_level_histogram = src_col->rep_level_hist;
                    ss->repetition_level_histogram_len = src_col->rep_hist_len;
                }
                if (src_col->def_hist_len > 1) {
                    ss->definition_level_histogram = src_col->def_level_hist;
                    ss->definition_level_histogram_len = src_col->def_hist_len;
                }
            }

            if (src_col->has_statistics) {
                meta->has_statistics = true;
                parquet_statistics_t* stats = &meta->statistics;
                memset(stats, 0, sizeof(*stats));

                if (src_col->has_null_count) {
                    stats->has_null_count = true;
                    stats->null_count = src_col->null_count;
                }

                if (src_col->has_distinct_count) {
                    stats->has_distinct_count = true;
                    stats->distinct_count = src_col->distinct_count;
                }

                if (src_col->has_min_max) {

                    /* Truncate variable-length min/max per the Parquet spec
                     * recommendation. Numeric / BOOLEAN / FLBA stats are
                     * already at their natural fixed width so pass through.
                     * The cap is configurable via
                     * carquet_writer_set_max_statistics_size. */
                    const size_t TRUNC = (size_t)writer->max_statistics_size;
                    bool variable_len =
                        (src_col->type == CARQUET_PHYSICAL_BYTE_ARRAY);

                    int32_t min_n = (int32_t)src_col->min_value_size;
                    int32_t max_n = (int32_t)src_col->max_value_size;
                    bool emit_max = true;
                    bool min_exact = true;
                    bool max_exact = true;

                    /* n == 0 is the empty BYTE_ARRAY: allocate one byte so
                     * the arena returns a usable (non-NULL) pointer. */
                    uint8_t* min_buf = carquet_arena_alloc(&writer->arena,
                        min_n > 0 ? (size_t)min_n : 1);
                    if (!min_buf) return CARQUET_ERROR_OUT_OF_MEMORY;
                    if (min_n > 0) memcpy(min_buf, src_col->min_value, (size_t)min_n);

                    if (variable_len && (size_t)min_n > TRUNC) {
                        /* Truncated prefix is lex <= original, valid as min. */
                        min_n = (int32_t)TRUNC;
                        min_exact = false;
                    }

                    uint8_t* max_buf = carquet_arena_alloc(&writer->arena,
                        max_n > 0 ? (size_t)max_n : 1);
                    if (!max_buf) return CARQUET_ERROR_OUT_OF_MEMORY;
                    if (max_n > 0) memcpy(max_buf, src_col->max_value, (size_t)max_n);

                    if (variable_len && (size_t)max_n > TRUNC) {
                        /* Truncate then increment to ensure result >= original.
                         * If all bytes in the prefix are 0xFF the increment
                         * wraps and we cannot emit a valid upper bound — drop
                         * max in that case. */
                        max_exact = false;
                        size_t new_len = TRUNC;
                        max_buf[new_len - 1]++;
                        while (new_len > 0 && max_buf[new_len - 1] == 0) {
                            new_len--;
                            if (new_len > 0) max_buf[new_len - 1]++;
                        }
                        if (new_len == 0) {
                            emit_max = false;
                        } else {
                            max_n = (int32_t)new_len;
                        }
                    }

                    stats->min_value = min_buf;
                    stats->min_value_len = min_n;
                    stats->has_is_min_value_exact = true;
                    stats->is_min_value_exact = min_exact;
                    if (emit_max) {
                        stats->max_value = max_buf;
                        stats->max_value_len = max_n;
                        stats->has_is_max_value_exact = true;
                        stats->is_max_value_exact = max_exact;
                    }

                    /* Mirror to deprecated min/max fields for older readers
                     * that ignore min_value/max_value (only for non-truncated
                     * stats to avoid leaking truncation semantics). */
                    bool deprecated_ok = j < writer->num_columns &&
                        deprecated_stats_are_compatible(&writer->columns[j]);
                    if (deprecated_ok &&
                        (!variable_len || (size_t)src_col->min_value_size <= TRUNC)) {
                        stats->min_deprecated = min_buf;
                        stats->min_deprecated_len = min_n;
                    }
                    if (deprecated_ok && emit_max &&
                        (!variable_len || (size_t)src_col->max_value_size <= TRUNC)) {
                        stats->max_deprecated = max_buf;
                        stats->max_deprecated_len = max_n;
                    }
                }
            }
        }
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Public API Implementation
 * ============================================================================
 */

carquet_writer_t* carquet_writer_create(
    const char* path,
    const carquet_schema_t* schema,
    const carquet_writer_options_t* options,
    carquet_error_t* error) {

    carquet_writer_t* writer = carquet_mem_calloc(1, sizeof(carquet_writer_t));
    if (!writer) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate writer");
        return NULL;
    }

    /* Initialize arena */
    if (carquet_arena_init_size(&writer->arena, 4096) != CARQUET_OK) {
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate arena");
        return NULL;
    }

    /* Open file */
    writer->file = fopen(path, "wb");
    if (!writer->file) {
        carquet_arena_destroy(&writer->arena);
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN, "Failed to open file for writing: %s", path);
        return NULL;
    }
    writer->owns_file = true;

    writer->path = carquet_heap_strdup(path);
    if (!writer->path) {
        fclose(writer->file);
        carquet_arena_destroy(&writer->arena);
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate path");
        return NULL;
    }

    /* Copy options */
    if (options) {
        writer->options = *options;
    } else {
        carquet_writer_options_init(&writer->options);
    }
    writer->max_statistics_size = 32;  /* Parquet spec recommendation; matches Arrow */

    /* Store full schema elements for metadata serialization */
    {
        carquet_status_t status = store_schema_elements(writer, schema);
        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to store schema elements");
            return NULL;
        }

        status = build_column_metadata_cache(writer, schema);
        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to build writer metadata cache");
            return NULL;
        }
    }

    /* Add leaf columns from schema (schema is nonnull per API contract) */
    for (int32_t i = 0; i < schema->num_leaves; i++) {
        int32_t elem_idx = schema->leaf_indices[i];
        parquet_schema_element_t* elem = &schema->elements[elem_idx];

        carquet_logical_type_t* lt = elem->has_logical_type ? &elem->logical_type : NULL;

        carquet_status_t status = add_column_internal(
            writer,
            elem->name,
            elem->type,
            lt,
            elem->repetition_type,
            elem->type_length,
            schema->max_def_levels[i],
            schema->max_rep_levels[i],
            leaf_statistics_sort_order_defined(schema, elem_idx));

        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to add column from schema");
            return NULL;
        }
    }

    return writer;
}

carquet_writer_t* carquet_writer_create_file(
    FILE* file,
    const carquet_schema_t* schema,
    const carquet_writer_options_t* options,
    carquet_error_t* error) {

    /* file and schema are nonnull per API contract */
    carquet_writer_t* writer = carquet_mem_calloc(1, sizeof(carquet_writer_t));
    if (!writer) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate writer");
        return NULL;
    }

    /* Initialize arena */
    if (carquet_arena_init_size(&writer->arena, 4096) != CARQUET_OK) {
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate arena");
        return NULL;
    }

    writer->file = file;
    writer->owns_file = false;

    /* Copy options */
    if (options) {
        writer->options = *options;
    } else {
        carquet_writer_options_init(&writer->options);
    }
    writer->max_statistics_size = 32;  /* Parquet spec recommendation; matches Arrow */

    /* Store full schema elements for metadata serialization */
    {
        carquet_status_t status = store_schema_elements(writer, schema);
        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to store schema elements");
            return NULL;
        }

        status = build_column_metadata_cache(writer, schema);
        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to build writer metadata cache");
            return NULL;
        }
    }

    /* Add leaf columns from schema (schema is nonnull per API contract) */
    for (int32_t i = 0; i < schema->num_leaves; i++) {
        int32_t elem_idx = schema->leaf_indices[i];
        parquet_schema_element_t* elem = &schema->elements[elem_idx];

        carquet_logical_type_t* lt = elem->has_logical_type ? &elem->logical_type : NULL;

        carquet_status_t status = add_column_internal(
            writer,
            elem->name,
            elem->type,
            lt,
            elem->repetition_type,
            elem->type_length,
            schema->max_def_levels[i],
            schema->max_rep_levels[i],
            leaf_statistics_sort_order_defined(schema, elem_idx));

        if (status != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, status, "Failed to add column from schema");
            return NULL;
        }
    }

    return writer;
}

/* ============================================================================
 * Append-mode helpers
 * ============================================================================
 *
 * Parse an existing file's footer + restore writer state so subsequent writes
 * are placed *after* the last byte of existing data (which sits just before
 * the existing footer's 8-byte tail). The existing footer is discarded; the
 * new one written by carquet_writer_close() lists existing + new row groups.
 */

/* Read PAR1 tail, footer-length, and footer bytes; parse into the writer's
 * arena. On success returns CARQUET_OK and fills *out_insert_offset with the
 * byte position where new data should start being written. */
static carquet_status_t append_parse_existing_footer(
    FILE* file,
    carquet_arena_t* arena,
    parquet_file_metadata_t* out_meta,
    int64_t* out_insert_offset,
    carquet_error_t* error) {

    if (carquet_fseek64(file, 0, SEEK_END) != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK,
            "append: failed to seek to end");
        return CARQUET_ERROR_FILE_SEEK;
    }
    int64_t file_size = carquet_ftell64(file);
    if (file_size < 12) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER,
            "append: file too small to contain a Parquet footer");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    /* Read the 8-byte tail: footer_length (uint32 LE) + PAR1 magic. */
    uint8_t tail[8];
    if (carquet_fseek64(file, file_size - 8, SEEK_SET) != 0 ||
        fread(tail, 1, 8, file) != 8) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ,
            "append: failed to read footer tail");
        return CARQUET_ERROR_FILE_READ;
    }
    if (memcmp(tail + 4, "PAR1", 4) != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC,
            "append: trailing PAR1 magic not found");
        return CARQUET_ERROR_INVALID_MAGIC;
    }
    uint32_t footer_len =
        (uint32_t)tail[0] |
        ((uint32_t)tail[1] << 8) |
        ((uint32_t)tail[2] << 16) |
        ((uint32_t)tail[3] << 24);
    if ((int64_t)footer_len + 8 > file_size) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER,
            "append: footer length exceeds file size");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    int64_t footer_offset = file_size - 8 - (int64_t)footer_len;
    uint8_t* footer_buf = carquet_mem_malloc(footer_len);
    if (!footer_buf) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "append: footer alloc");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (carquet_fseek64(file, footer_offset, SEEK_SET) != 0 ||
        fread(footer_buf, 1, footer_len, file) != footer_len) {
        carquet_mem_free(footer_buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ,
            "append: failed to read footer bytes");
        return CARQUET_ERROR_FILE_READ;
    }

    carquet_status_t status = parquet_parse_file_metadata(
        footer_buf, footer_len, arena, out_meta, error);
    carquet_mem_free(footer_buf);
    if (status != CARQUET_OK) return status;

    /* New row groups start where the existing footer used to start, so the
     * existing data pages and any bloom filters / page indexes are preserved
     * (those sit between the last row group's data and the old footer). */
    *out_insert_offset = footer_offset;
    return CARQUET_OK;
}

/* Verify the caller-supplied schema describes the same leaf columns as the
 * existing file. We only check the structural facts the writer relies on
 * (count, name, physical type, repetition) — strict enough that a successful
 * append produces a self-consistent file, loose enough that small writer-
 * option differences (compression, page size) don't trip it. */
static carquet_status_t append_validate_schema_matches(
    const carquet_schema_t* user_schema,
    const parquet_file_metadata_t* parsed,
    carquet_error_t* error) {

    /* Count parsed leaves. */
    int32_t parsed_leaves = 0;
    for (int32_t i = 1; i < parsed->num_schema_elements; i++) {
        if (parsed->schema[i].num_children == 0) parsed_leaves++;
    }
    if (parsed_leaves != user_schema->num_leaves) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_SCHEMA,
            "append: schema leaf count mismatch (existing=%d, supplied=%d)",
            parsed_leaves, user_schema->num_leaves);
        return CARQUET_ERROR_INVALID_SCHEMA;
    }

    /* Walk parsed leaves in document order and compare to user leaves. */
    int32_t leaf_idx = 0;
    for (int32_t i = 1; i < parsed->num_schema_elements; i++) {
        if (parsed->schema[i].num_children != 0) continue;
        int32_t user_elem_idx = user_schema->leaf_indices[leaf_idx];
        const parquet_schema_element_t* u = &user_schema->elements[user_elem_idx];
        const parquet_schema_element_t* p = &parsed->schema[i];

        if (u->type != p->type ||
            u->repetition_type != p->repetition_type ||
            (u->type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY &&
             u->type_length != p->type_length)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_SCHEMA,
                "append: column %d type/repetition mismatch", leaf_idx);
            return CARQUET_ERROR_INVALID_SCHEMA;
        }
        /* Logical type must match too: appending a BYTE_ARRAY annotated JSON
         * onto chunks annotated STRING (or any logical-type divergence) yields
         * a file whose row groups disagree on semantics, which stricter readers
         * (Java/Rust) reject. Physical type alone is not enough. */
        if (u->has_logical_type != p->has_logical_type ||
            (u->has_logical_type &&
             u->logical_type.id != p->logical_type.id)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_SCHEMA,
                "append: column %d logical type mismatch", leaf_idx);
            return CARQUET_ERROR_INVALID_SCHEMA;
        }
        if (!u->name || !p->name || strcmp(u->name, p->name) != 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_SCHEMA,
                "append: column %d name mismatch", leaf_idx);
            return CARQUET_ERROR_INVALID_SCHEMA;
        }
        leaf_idx++;
    }
    return CARQUET_OK;
}

/* Deep-copy one parsed column chunk into the writer's row_group_column_info_t
 * representation. Statistics min/max byte buffers are heap-allocated to match
 * the lifetime of writer->row_groups[i].columns[j], which free_row_groups
 * (called from writer_abort/close) frees individually. */
static carquet_status_t append_restore_column_chunk(
    const parquet_column_chunk_t* src,
    row_group_column_info_t* dst) {

    const parquet_column_metadata_t* m = &src->metadata;
    dst->file_offset = src->file_offset;
    dst->type = m->type;
    dst->codec = m->codec;
    dst->num_values = m->num_values;
    dst->total_compressed_size = m->total_compressed_size;
    dst->total_uncompressed_size = m->total_uncompressed_size;
    dst->data_page_offset = m->data_page_offset;
    if (m->has_dictionary_page_offset) {
        dst->has_dictionary_page_offset = true;
        dst->dictionary_page_offset = m->dictionary_page_offset;
        dst->has_dictionary_page = true;
    }
    /* Copy the existing encodings list (cap at the writer's storage of 3). */
    int32_t ne = m->num_encodings < 3 ? m->num_encodings : 3;
    for (int32_t k = 0; k < ne; k++) dst->encodings[k] = m->encodings[k];
    dst->num_encodings = ne;

    if (m->has_bloom_filter_offset) {
        dst->has_bloom_filter_offset = true;
        dst->bloom_filter_offset = m->bloom_filter_offset;
    }
    if (m->has_bloom_filter_length) {
        dst->has_bloom_filter_length = true;
        dst->bloom_filter_length = m->bloom_filter_length;
    }
    if (src->has_column_index_offset) {
        dst->has_column_index_offset = true;
        dst->column_index_offset = src->column_index_offset;
    }
    if (src->has_column_index_length) {
        dst->has_column_index_length = true;
        dst->column_index_length = src->column_index_length;
    }
    if (src->has_offset_index_offset) {
        dst->has_offset_index_offset = true;
        dst->offset_index_offset = src->offset_index_offset;
    }
    if (src->has_offset_index_length) {
        dst->has_offset_index_length = true;
        dst->offset_index_length = src->offset_index_length;
    }

    if (m->has_geospatial_statistics) {
        dst->has_geo_stats = true;
        dst->geo_stats = m->geospatial_statistics;
    }

    /* Preserve SizeStatistics (Parquet 2.9) across append by copying the parsed
     * histograms into owned storage. */
    if (m->has_size_statistics) {
        const parquet_size_statistics_t* ss = &m->size_statistics;
        dst->has_size_statistics = true;
        dst->unencoded_ba_bytes = ss->has_unencoded_byte_array_data_bytes
            ? ss->unencoded_byte_array_data_bytes : -1;
        if (ss->repetition_level_histogram && ss->repetition_level_histogram_len > 0) {
            dst->rep_level_hist = carquet_mem_malloc(
                (size_t)ss->repetition_level_histogram_len * sizeof(int64_t));
            if (dst->rep_level_hist) {
                memcpy(dst->rep_level_hist, ss->repetition_level_histogram,
                       (size_t)ss->repetition_level_histogram_len * sizeof(int64_t));
                dst->rep_hist_len = ss->repetition_level_histogram_len;
            }
        }
        if (ss->definition_level_histogram && ss->definition_level_histogram_len > 0) {
            dst->def_level_hist = carquet_mem_malloc(
                (size_t)ss->definition_level_histogram_len * sizeof(int64_t));
            if (dst->def_level_hist) {
                memcpy(dst->def_level_hist, ss->definition_level_histogram,
                       (size_t)ss->definition_level_histogram_len * sizeof(int64_t));
                dst->def_hist_len = ss->definition_level_histogram_len;
            }
        }
    }

    if (m->has_statistics) {
        dst->has_statistics = true;
        const parquet_statistics_t* s = &m->statistics;
        if (s->has_null_count) {
            dst->has_null_count = true;
            dst->null_count = s->null_count;
        }
        if (s->has_distinct_count) {
            dst->has_distinct_count = true;
            dst->distinct_count = s->distinct_count;
        }
        if (s->min_value && s->min_value_len > 0 &&
            s->max_value && s->max_value_len > 0) {
            dst->has_min_max = true;
            dst->min_value = carquet_mem_malloc((size_t)s->min_value_len);
            dst->max_value = carquet_mem_malloc((size_t)s->max_value_len);
            if (!dst->min_value || !dst->max_value) {
                carquet_mem_free(dst->min_value);
                carquet_mem_free(dst->max_value);
                dst->min_value = dst->max_value = NULL;
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            memcpy(dst->min_value, s->min_value, (size_t)s->min_value_len);
            memcpy(dst->max_value, s->max_value, (size_t)s->max_value_len);
            dst->min_value_size = s->min_value_len;
            dst->max_value_size = s->max_value_len;
        }
    }
    return CARQUET_OK;
}

static carquet_status_t append_restore_row_groups(
    carquet_writer_t* writer,
    const parquet_file_metadata_t* parsed) {

    if (parsed->num_row_groups == 0) return CARQUET_OK;

    writer->row_groups = carquet_mem_calloc(
        (size_t)parsed->num_row_groups, sizeof(row_group_info_t));
    if (!writer->row_groups) return CARQUET_ERROR_OUT_OF_MEMORY;
    writer->row_groups_capacity = parsed->num_row_groups;

    for (int32_t i = 0; i < parsed->num_row_groups; i++) {
        const parquet_row_group_t* src = &parsed->row_groups[i];
        row_group_info_t* dst = &writer->row_groups[i];

        /* The close-time footer rewrite indexes the writer's per-column
         * arrays (paths, encodings — all sized to the schema leaf count) by
         * each restored row group's column position. A malformed file can
         * declare a row group with a different num_columns than the schema
         * has leaves; restoring it would read those arrays out of bounds at
         * close. Reject such files (open_append fails, file untouched). */
        if (src->num_columns != writer->num_columns) {
            return CARQUET_ERROR_INVALID_SCHEMA;
        }

        dst->num_rows = src->num_rows;
        dst->total_byte_size = src->total_byte_size;
        dst->total_compressed_size = src->has_total_compressed_size
            ? src->total_compressed_size : 0;
        dst->file_offset = src->has_file_offset ? src->file_offset : 0;
        dst->ordinal = src->has_ordinal ? src->ordinal : (int16_t)i;
        dst->num_columns = src->num_columns;
        dst->columns = carquet_mem_calloc(
            (size_t)src->num_columns, sizeof(row_group_column_info_t));
        if (!dst->columns) return CARQUET_ERROR_OUT_OF_MEMORY;
        /* This row group now owns a columns allocation. Publish it to
         * writer->num_row_groups immediately so that if a later step fails
         * (column-chunk restore below, or a subsequent row group), the abort
         * path's free_row_groups() reclaims it instead of leaking. */
        writer->num_row_groups = i + 1;

        for (int32_t j = 0; j < src->num_columns; j++) {
            carquet_status_t s = append_restore_column_chunk(
                &src->columns[j], &dst->columns[j]);
            if (s != CARQUET_OK) return s;
        }

        writer->total_rows += src->num_rows;
    }
    writer->num_row_groups = parsed->num_row_groups;
    return CARQUET_OK;
}

static carquet_status_t append_restore_kv_metadata(
    carquet_writer_t* writer,
    const parquet_file_metadata_t* parsed) {

    if (parsed->num_key_value == 0) return CARQUET_OK;

    writer->kv_metadata = carquet_mem_calloc(
        (size_t)parsed->num_key_value, sizeof(parquet_key_value_t));
    if (!writer->kv_metadata) return CARQUET_ERROR_OUT_OF_MEMORY;
    writer->kv_metadata_capacity = parsed->num_key_value;

    for (int32_t i = 0; i < parsed->num_key_value; i++) {
        const parquet_key_value_t* src = &parsed->key_value_metadata[i];
        if (src->key) {
            writer->kv_metadata[i].key = carquet_heap_strdup(src->key);
            if (!writer->kv_metadata[i].key) return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (src->value) {
            writer->kv_metadata[i].value = carquet_heap_strdup(src->value);
            if (!writer->kv_metadata[i].value) return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        writer->num_kv_metadata++;
    }
    return CARQUET_OK;
}

carquet_writer_t* carquet_writer_open_append(
    const char* path,
    const carquet_schema_t* schema,
    const carquet_writer_options_t* options,
    carquet_error_t* error) {

    /* path and schema are nonnull per API contract */
    carquet_writer_t* writer = carquet_mem_calloc(1, sizeof(carquet_writer_t));
    if (!writer) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "append: writer alloc");
        return NULL;
    }
    if (carquet_arena_init_size(&writer->arena, 4096) != CARQUET_OK) {
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "append: arena init");
        return NULL;
    }

    /* "r+b": read+write, do not truncate, fail if missing. */
    writer->file = fopen(path, "r+b");
    if (!writer->file) {
        carquet_arena_destroy(&writer->arena);
        carquet_mem_free(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN,
            "append: cannot open %s for read+write", path);
        return NULL;
    }
    writer->owns_file = true;
    writer->is_append = true;
    writer->path = carquet_heap_strdup(path);

    if (options) writer->options = *options;
    else carquet_writer_options_init(&writer->options);
    writer->max_statistics_size = 32;

    /* Parse footer + validate it against the user's schema. */
    parquet_file_metadata_t parsed;
    int64_t insert_offset = 0;
    carquet_status_t st = append_parse_existing_footer(
        writer->file, &writer->arena, &parsed, &insert_offset, error);
    if (st != CARQUET_OK) { carquet_writer_abort(writer); return NULL; }

    st = append_validate_schema_matches(schema, &parsed, error);
    if (st != CARQUET_OK) { carquet_writer_abort(writer); return NULL; }

    /* From here on, set up the writer the same way create_file does. The
     * normal column metadata cache is built from the user's schema (which we
     * just proved matches), so encodings / paths used by the new row groups
     * are consistent with existing ones. */
    st = store_schema_elements(writer, schema);
    if (st != CARQUET_OK) {
        carquet_writer_abort(writer);
        CARQUET_SET_ERROR(error, st, "append: store schema elements");
        return NULL;
    }
    st = build_column_metadata_cache(writer, schema);
    if (st != CARQUET_OK) {
        carquet_writer_abort(writer);
        CARQUET_SET_ERROR(error, st, "append: build column metadata cache");
        return NULL;
    }
    for (int32_t i = 0; i < schema->num_leaves; i++) {
        int32_t elem_idx = schema->leaf_indices[i];
        parquet_schema_element_t* elem = &schema->elements[elem_idx];
        carquet_logical_type_t* lt = elem->has_logical_type ? &elem->logical_type : NULL;
        st = add_column_internal(writer, elem->name, elem->type, lt,
            elem->repetition_type, elem->type_length,
            schema->max_def_levels[i], schema->max_rep_levels[i],
            leaf_statistics_sort_order_defined(schema, elem_idx));
        if (st != CARQUET_OK) {
            carquet_writer_abort(writer);
            CARQUET_SET_ERROR(error, st, "append: add column from schema");
            return NULL;
        }
    }

    /* Restore the prior row groups + key-value metadata into the writer so
     * the close-time footer rewrite emits them ahead of the new ones. */
    st = append_restore_row_groups(writer, &parsed);
    if (st != CARQUET_OK) {
        carquet_writer_abort(writer);
        CARQUET_SET_ERROR(error, st, "append: restore row groups");
        return NULL;
    }
    st = append_restore_kv_metadata(writer, &parsed);
    if (st != CARQUET_OK) {
        carquet_writer_abort(writer);
        CARQUET_SET_ERROR(error, st, "append: restore key-value metadata");
        return NULL;
    }

    /* Position the file at the insertion point and mark the header as already
     * written so the close path does not re-emit PAR1. */
    if (carquet_fseek64(writer->file, insert_offset, SEEK_SET) != 0) {
        carquet_writer_abort(writer);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK,
            "append: seek to insertion offset");
        return NULL;
    }
    writer->file_offset = insert_offset;
    writer->header_written = true;
    return writer;
}

/* 1000^|rank diff| between time units (MILLIS=0, MICROS=1, NANOS=2). */
static int64_t timestamp_unit_factor(carquet_time_unit_t a,
                                     carquet_time_unit_t b) {
    int d = (int)a - (int)b;
    if (d < 0) d = -d;
    int64_t f = 1;
    for (int i = 0; i < d; i++) f *= 1000;
    return f;
}

/* Rescale `n` INT64 timestamps from `src` to `dst` unit in place into out.
 * Returns CARQUET_ERROR_INVALID_ARGUMENT on disallowed truncation/overflow. */
static carquet_status_t coerce_timestamp_values(
    const int64_t* in, int64_t* out, int64_t n,
    carquet_time_unit_t src, carquet_time_unit_t dst,
    bool allow_truncation) {

    int64_t factor = timestamp_unit_factor(src, dst);
    if ((int)dst > (int)src) {
        /* finer target: multiply, guard overflow */
        int64_t lim = INT64_MAX / factor;
        for (int64_t i = 0; i < n; i++) {
            int64_t v = in[i];
            if (v > lim || v < -lim) return CARQUET_ERROR_INVALID_ARGUMENT;
            out[i] = v * factor;
        }
    } else {
        /* coarser target: divide (truncates toward zero) */
        for (int64_t i = 0; i < n; i++) {
            int64_t v = in[i];
            if (!allow_truncation && (v % factor) != 0) {
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            out[i] = v / factor;
        }
    }
    return CARQUET_OK;
}

carquet_status_t carquet_writer_write_batch(
    carquet_writer_t* writer,
    int32_t column_index,
    const void* values,
    int64_t num_values,
    const int16_t* def_levels,
    const int16_t* rep_levels) {

    /* writer and values are nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Ensure header is written */
    carquet_status_t status = ensure_header_written(writer);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Ensure we have a row group */
    status = ensure_row_group(writer);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Optional TIMESTAMP coercion: rescale the (packed, non-null) INT64
     * values from the schema-declared unit to the target unit. */
    const void* write_values = values;
    int64_t* coerced = NULL;
    const writer_column_def_t* cdef = &writer->columns[column_index];
    if (writer->options.coerce_timestamps &&
        cdef->physical_type == CARQUET_PHYSICAL_INT64 &&
        cdef->logical_type.id == CARQUET_LOGICAL_TIMESTAMP) {
        carquet_time_unit_t src = cdef->logical_type.params.timestamp.unit;
        carquet_time_unit_t dst = writer->options.coerce_timestamp_unit;
        if (src != dst && num_values > 0) {
            int64_t present = num_values;
            if (def_levels && cdef->max_def_level > 0) {
                present = carquet_dispatch_count_non_nulls(
                    def_levels, num_values, cdef->max_def_level);
            }
            coerced = carquet_mem_malloc((size_t)present * sizeof(int64_t));
            if (!coerced) return CARQUET_ERROR_OUT_OF_MEMORY;
            status = coerce_timestamp_values(
                (const int64_t*)values, coerced, present, src, dst,
                writer->options.allow_timestamp_truncation);
            if (status != CARQUET_OK) {
                carquet_mem_free(coerced);
                return status;
            }
            write_values = coerced;
        }
    }

    /* Write to the row group */
    status = carquet_row_group_writer_write_column(
        writer->current_row_group,
        column_index,
        write_values,
        num_values,
        def_levels,
        rep_levels);

    carquet_mem_free(coerced);

    if (status != CARQUET_OK) {
        return status;
    }

    writer->column_values_written[column_index] += num_values;

    /* Track rows (use column 0 as reference).
     * For repeated columns (max_rep_level > 0), the number of logical rows
     * is the count of rep_level == 0 entries (new top-level records).
     * For non-repeated columns, num_values == num_rows. */
    if (column_index == 0) {
        if (rep_levels && writer->columns[0].max_rep_level > 0) {
            int64_t rows = 0;
            for (int64_t i = 0; i < num_values; i++) {
                if (rep_levels[i] == 0) rows++;
            }
            writer->current_row_group_rows += rows;
        } else {
            writer->current_row_group_rows += num_values;
        }
    }

    writer->current_row_group_estimated_bytes = saturating_add_i64(
        writer->current_row_group_estimated_bytes,
        estimate_column_batch_bytes(
            &writer->columns[column_index],
            values,
            num_values,
            def_levels,
            rep_levels));

    if (writer_supports_aligned_auto_flush(writer) &&
        writer->current_row_group_estimated_bytes > writer->options.row_group_size &&
        current_row_group_is_aligned(writer)) {
        status = flush_row_group(writer);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    return CARQUET_OK;
}

/* ============================================================================
 * Nested write helper (single-level LIST / MAP auto-shredding)
 * ============================================================================ */

/* Byte width of one packed value for a leaf column. */
static size_t leaf_value_size(const writer_column_def_t* c) {
    switch (c->physical_type) {
        case CARQUET_PHYSICAL_BOOLEAN:              return sizeof(uint8_t);
        case CARQUET_PHYSICAL_INT32:                return sizeof(int32_t);
        case CARQUET_PHYSICAL_INT64:                return sizeof(int64_t);
        case CARQUET_PHYSICAL_INT96:                return 12;
        case CARQUET_PHYSICAL_FLOAT:                return sizeof(float);
        case CARQUET_PHYSICAL_DOUBLE:               return sizeof(double);
        case CARQUET_PHYSICAL_BYTE_ARRAY:           return sizeof(carquet_byte_array_t);
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return (size_t)c->type_length;
        default:                                    return 0;
    }
}

/* Arrow-style (LSB-first) validity bit: 1 = valid/present. */
static inline bool validity_bit(const uint8_t* bm, int64_t i) {
    return (bm[i >> 3] >> (i & 7)) & 1u;
}

carquet_status_t carquet_writer_write_list_column(
    carquet_writer_t* writer,
    int32_t column_index,
    int64_t num_lists,
    const int32_t* offsets,
    const uint8_t* list_validity,
    const void* values,
    const uint8_t* value_validity,
    carquet_error_t* error) {

    if (column_index < 0 || column_index >= writer->num_columns) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "column_index %d out of range [0, %d)", column_index,
            writer->num_columns);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (num_lists < 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "num_lists must be non-negative");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const writer_column_def_t* col = &writer->columns[column_index];

    /* This front-end shreds the standard single-level LIST/MAP encoding only:
     * a repeated group with exactly one repetition level, an optional or
     * required container above it, and an optional or required leaf below it.
     * That covers every column produced by carquet_schema_add_list() and
     * carquet_schema_add_map(). Deeper nesting is out of scope. */
    if (col->max_rep_level != 1) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "column %d is not a single-level repeated column "
            "(max_rep_level=%d)", column_index, (int)col->max_rep_level);
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
    bool leaf_optional = (col->repetition == CARQUET_REPETITION_OPTIONAL);
    int16_t max_def = col->max_def_level;
    /* max_def = container_optional + 1 (repeated) + leaf_optional. */
    int container_optional = (int)max_def - 1 - (leaf_optional ? 1 : 0);
    if (container_optional != 0 && container_optional != 1) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "column %d has an unsupported nested shape "
            "(max_def_level=%d)", column_index, (int)max_def);
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
    int16_t def_present = max_def;
    int16_t def_null_elem = (int16_t)(container_optional + 1);
    int16_t def_empty = (int16_t)container_optional;

    if (num_lists > 0 && offsets == NULL) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "offsets must be non-NULL when num_lists > 0");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (num_lists == 0) {
        return CARQUET_OK;  /* Nothing to write. */
    }
    if (offsets[0] != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "offsets[0] must be 0 (sliced arrays are not supported)");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    /* Validate monotonic offsets and derive the child count. */
    for (int64_t i = 0; i < num_lists; i++) {
        if (offsets[i + 1] < offsets[i]) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "offsets must be non-decreasing (offsets[%lld]=%d < offsets[%lld]=%d)",
                (long long)(i + 1), offsets[i + 1], (long long)i, offsets[i]);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
    }
    int64_t total_children = offsets[num_lists];

    /* Upper bounds: one level per null/empty list plus one per child. */
    if (total_children > INT64_MAX - num_lists) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "list column too large");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int64_t level_cap = num_lists + total_children;

    size_t vsize = leaf_value_size(col);
    if (vsize == 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "unsupported leaf physical type for column %d", column_index);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int16_t* def_levels = carquet_mem_malloc((size_t)level_cap * sizeof(int16_t));
    int16_t* rep_levels = carquet_mem_malloc((size_t)level_cap * sizeof(int16_t));
    uint8_t* packed = (total_children > 0)
        ? carquet_mem_malloc((size_t)total_children * vsize) : NULL;
    if (!def_levels || !rep_levels || (total_children > 0 && !packed)) {
        carquet_mem_free(def_levels);
        carquet_mem_free(rep_levels);
        carquet_mem_free(packed);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Failed to allocate shredding buffers");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    const uint8_t* vbytes = (const uint8_t*)values;
    int64_t nlevels = 0;
    int64_t npacked = 0;
    carquet_status_t st = CARQUET_OK;

    for (int64_t i = 0; i < num_lists; i++) {
        bool present = list_validity ? validity_bit(list_validity, i) : true;
        if (!present) {
            if (!container_optional) {
                st = CARQUET_ERROR_INVALID_ARGUMENT;
                CARQUET_SET_ERROR(error, st,
                    "null list at row %lld but the container is REQUIRED",
                    (long long)i);
                goto done;
            }
            def_levels[nlevels] = 0;      /* null container */
            rep_levels[nlevels] = 0;
            nlevels++;
            continue;
        }
        int32_t start = offsets[i];
        int32_t end = offsets[i + 1];
        if (start == end) {
            def_levels[nlevels] = def_empty;  /* present but empty list */
            rep_levels[nlevels] = 0;
            nlevels++;
            continue;
        }
        for (int32_t j = start; j < end; j++) {
            rep_levels[nlevels] = (j == start) ? 0 : 1;
            bool vpresent = value_validity ? validity_bit(value_validity, j) : true;
            if (vpresent) {
                def_levels[nlevels] = def_present;
                memcpy(packed + (size_t)npacked * vsize,
                       vbytes + (size_t)j * vsize, vsize);
                npacked++;
            } else {
                if (!leaf_optional) {
                    st = CARQUET_ERROR_INVALID_ARGUMENT;
                    CARQUET_SET_ERROR(error, st,
                        "null element at child %d but the leaf is REQUIRED", j);
                    goto done;
                }
                def_levels[nlevels] = def_null_elem;
            }
            nlevels++;
        }
    }

    st = carquet_writer_write_batch(writer, column_index,
        packed ? (const void*)packed : (const void*)vbytes,
        nlevels, def_levels, rep_levels);

done:
    carquet_mem_free(def_levels);
    carquet_mem_free(rep_levels);
    carquet_mem_free(packed);
    return st;
}

carquet_status_t carquet_writer_new_row_group(carquet_writer_t* writer) {
    /* writer is nonnull per API contract */
    /* Ensure header is written */
    carquet_status_t status = ensure_header_written(writer);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Flush current row group if any */
    return flush_row_group(writer);
}

int32_t carquet_writer_num_columns(const carquet_writer_t* writer) {
    /* writer is nonnull per API contract */
    return writer->num_columns;
}

carquet_status_t carquet_writer_close(carquet_writer_t* writer) {
    /* writer is nonnull per API contract */
    carquet_status_t status = CARQUET_OK;

    /* Ensure header is written */
    status = ensure_header_written(writer);
    if (status != CARQUET_OK) {
        goto cleanup;
    }

    /* Flush any pending row group */
    status = flush_row_group(writer);
    if (status != CARQUET_OK) {
        goto cleanup;
    }
    status = pipeline_complete_pending(writer);
    if (status != CARQUET_OK) {
        goto cleanup;
    }

    /* Every row group must be on disk before the footer goes after it. */
    status = writer_io_drain(writer);
    if (status != CARQUET_OK) {
        goto cleanup;
    }

    /* Build file metadata */
    parquet_file_metadata_t metadata;
    status = build_file_metadata(writer, &metadata);
    if (status != CARQUET_OK) {
        goto cleanup;
    }

    /* Serialize metadata to buffer */
    carquet_buffer_t metadata_buffer;
    carquet_buffer_init(&metadata_buffer);

    status = parquet_write_file_metadata(&metadata, &metadata_buffer, NULL);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&metadata_buffer);
        goto cleanup;
    }

    /* Write metadata */
    if (fwrite(metadata_buffer.data, 1, metadata_buffer.size, writer->file) != metadata_buffer.size) {
        carquet_buffer_destroy(&metadata_buffer);
        status = CARQUET_ERROR_FILE_WRITE;
        goto cleanup;
    }

    /* Write metadata length (4 bytes, little-endian) */
    uint32_t metadata_len = (uint32_t)metadata_buffer.size;
    uint8_t len_bytes[4];
    len_bytes[0] = (uint8_t)(metadata_len & 0xFF);
    len_bytes[1] = (uint8_t)((metadata_len >> 8) & 0xFF);
    len_bytes[2] = (uint8_t)((metadata_len >> 16) & 0xFF);
    len_bytes[3] = (uint8_t)((metadata_len >> 24) & 0xFF);

    if (fwrite(len_bytes, 1, 4, writer->file) != 4) {
        carquet_buffer_destroy(&metadata_buffer);
        status = CARQUET_ERROR_FILE_WRITE;
        goto cleanup;
    }

    carquet_buffer_destroy(&metadata_buffer);

    /* Write footer magic */
    status = write_magic(writer->file);
    if (status != CARQUET_OK) {
        goto cleanup;
    }

    /* For buffer writers, read back the entire file into memory.
     * Use 64-bit seek/tell to handle files >2 GB on all platforms. */
    if (writer->is_buffer_writer && writer->file) {
        fflush(writer->file);
        int64_t file_size = -1;
        if (carquet_fseek64(writer->file, 0, SEEK_END) == 0)
            file_size = carquet_ftell64(writer->file);
        if (file_size > 0) {
            writer->output_buffer = carquet_mem_malloc((size_t)file_size);
            if (!writer->output_buffer) {
                status = CARQUET_ERROR_OUT_OF_MEMORY;
                goto cleanup;
            }
            rewind(writer->file);
            size_t nread = fread(writer->output_buffer, 1, (size_t)file_size, writer->file);
            if (nread != (size_t)file_size) {
                carquet_mem_free(writer->output_buffer);
                writer->output_buffer = NULL;
                status = CARQUET_ERROR_FILE_READ;
                goto cleanup;
            }
            writer->output_buffer_size = (size_t)file_size;
        }
    }

cleanup:
    /* Free resources. Encode tasks first (they fill page buffers), then the
     * I/O thread (it reads them), then the row-group writers that own them. */
    pipeline_quiesce(writer);
    writer_io_shutdown(writer);
    pipeline_free(writer);
    if (writer->current_row_group) {
        carquet_row_group_writer_destroy(writer->current_row_group);
        writer->current_row_group = NULL;
    }

    if (writer->owns_file && writer->file) {
        fclose(writer->file);
        writer->file = NULL;
    }

    /* Free column definitions. NULL every freed pointer: a buffer writer is
     * kept alive past close() (see the is_buffer_writer return below) and its
     * documented lifecycle is close() -> get_buffer() -> abort(). abort()
     * re-frees these same members, so leaving dangling pointers here causes a
     * double free / use-after-free. The free_* helpers below already NULL
     * their own state. */
    if (writer->columns) {
        for (int32_t i = 0; i < writer->num_columns; i++) {
            carquet_mem_free(writer->columns[i].name);
        }
        carquet_mem_free(writer->columns);
        writer->columns = NULL;
    }

    /* Free schema elements */
    if (writer->schema_elements) {
        free_schema_field_metadata(writer->schema_elements, writer->num_schema_elements);
        for (int32_t i = 0; i < writer->num_schema_elements; i++) {
            carquet_mem_free(writer->schema_elements[i].name);
        }
        carquet_mem_free(writer->schema_elements);
        writer->schema_elements = NULL;
    }
    if (writer->column_paths) {
        for (int32_t i = 0; i < writer->num_columns; i++) {
            carquet_mem_free(writer->column_paths[i]);
        }
        carquet_mem_free(writer->column_paths);
        writer->column_paths = NULL;
    }
    carquet_mem_free(writer->column_path_lens);
    writer->column_path_lens = NULL;
    carquet_mem_free(writer->column_encodings);
    writer->column_encodings = NULL;
    carquet_mem_free(writer->column_num_encodings);
    writer->column_num_encodings = NULL;

    carquet_mem_free(writer->column_values_written);
    writer->column_values_written = NULL;
    free_row_groups(writer);
    carquet_mem_free(writer->path);
    writer->path = NULL;

    /* Free key-value metadata */
    free_kv_metadata(writer);

    /* Free per-column overrides */
    free_column_overrides(writer);

    /* For buffer writers, keep the writer alive so get_buffer can be called.
     * The writer struct (and output_buffer) will be freed by get_buffer or
     * by a subsequent abort call. */
    if (writer->is_buffer_writer) {
        return status;
    }

    carquet_mem_free(writer->output_buffer);
    carquet_arena_destroy(&writer->arena);
    carquet_mem_free(writer);

    return status;
}

void carquet_writer_abort(carquet_writer_t* writer) {
    if (!writer) return;

    pipeline_quiesce(writer);
    writer_io_shutdown(writer);
    pipeline_free(writer);

    /* Cleanup row group */
    if (writer->current_row_group) {
        carquet_row_group_writer_destroy(writer->current_row_group);
        writer->current_row_group = NULL;
    }

    /* Close and delete file. For append writers the file existed before we
     * opened it, so deleting it on abort would destroy the user's data —
     * close it but leave it on disk untouched. */
    if (writer->owns_file && writer->file) {
        fclose(writer->file);
        writer->file = NULL;

        if (writer->path && !writer->is_append) {
            remove(writer->path);
        }
    }

    /* Free column definitions */
    if (writer->columns) {
        for (int32_t i = 0; i < writer->num_columns; i++) {
            carquet_mem_free(writer->columns[i].name);
        }
        carquet_mem_free(writer->columns);
    }

    /* Free schema elements */
    if (writer->schema_elements) {
        free_schema_field_metadata(writer->schema_elements, writer->num_schema_elements);
        for (int32_t i = 0; i < writer->num_schema_elements; i++) {
            carquet_mem_free(writer->schema_elements[i].name);
        }
        carquet_mem_free(writer->schema_elements);
    }
    if (writer->column_paths) {
        for (int32_t i = 0; i < writer->num_columns; i++) {
            carquet_mem_free(writer->column_paths[i]);
        }
        carquet_mem_free(writer->column_paths);
    }
    carquet_mem_free(writer->column_path_lens);
    carquet_mem_free(writer->column_encodings);
    carquet_mem_free(writer->column_num_encodings);

    carquet_mem_free(writer->column_values_written);
    free_row_groups(writer);
    carquet_mem_free(writer->path);

    /* Free key-value metadata */
    free_kv_metadata(writer);

    /* Free per-column overrides */
    free_column_overrides(writer);

    /* Free buffer writer output */
    carquet_mem_free(writer->output_buffer);

    carquet_arena_destroy(&writer->arena);
    carquet_mem_free(writer);
}

/* ============================================================================
 * Key-Value Metadata API
 * ============================================================================
 */

carquet_status_t carquet_writer_add_metadata(
    carquet_writer_t* writer,
    const char* key,
    const char* value) {

    /* writer and key are nonnull per API contract */
    if (writer->num_kv_metadata >= writer->kv_metadata_capacity) {
        int32_t new_cap = writer->kv_metadata_capacity == 0 ? 8 : writer->kv_metadata_capacity * 2;
        parquet_key_value_t* new_kv = carquet_mem_realloc(writer->kv_metadata,
            new_cap * sizeof(parquet_key_value_t));
        if (!new_kv) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        writer->kv_metadata = new_kv;
        writer->kv_metadata_capacity = new_cap;
    }

    parquet_key_value_t* entry = &writer->kv_metadata[writer->num_kv_metadata];
    entry->key = carquet_heap_strdup(key);
    if (!entry->key) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    entry->value = value ? carquet_heap_strdup(value) : NULL;
    if (value && !entry->value) {
        carquet_mem_free(entry->key);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    writer->num_kv_metadata++;
    return CARQUET_OK;
}

/* ============================================================================
 * Per-Column Writer Options API
 * ============================================================================
 */

carquet_status_t carquet_writer_set_column_encoding(
    carquet_writer_t* writer,
    int32_t column_index,
    carquet_encoding_t encoding) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (!writer_encoding_supported(
            encoding, writer->columns[column_index].physical_type)) {
        return CARQUET_ERROR_INVALID_ENCODING;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_encoding_overrides[column_index] = encoding;
    writer->column_encoding_override_set[column_index] = true;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_column_compression(
    carquet_writer_t* writer,
    int32_t column_index,
    carquet_compression_t codec,
    int32_t level) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_compression_overrides[column_index] = codec;
    writer->column_compression_levels[column_index] = level;
    writer->column_compression_override_set[column_index] = true;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_column_page_size(
    carquet_writer_t* writer,
    int32_t column_index,
    int64_t bytes) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns ||
        bytes <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_page_size_overrides[column_index] = bytes;
    writer->column_page_size_override_set[column_index] = true;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_max_statistics_size(
    carquet_writer_t* writer,
    int64_t bytes) {

    /* writer is nonnull per API contract */
    if (bytes <= 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    writer->max_statistics_size = bytes;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_column_statistics(
    carquet_writer_t* writer,
    int32_t column_index,
    bool enabled) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_statistics_overrides[column_index] = enabled;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_column_bloom_filter(
    carquet_writer_t* writer,
    int32_t column_index,
    bool enabled) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_bloom_filter_overrides[column_index] = enabled;
    writer->column_bloom_explicit[column_index] = true;
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_column_bloom_filter_options(
    carquet_writer_t* writer,
    int32_t column_index,
    bool enabled,
    int64_t ndv,
    double fpp) {

    /* writer is nonnull per API contract */
    if (column_index < 0 || column_index >= writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t status = ensure_column_overrides(writer);
    if (status != CARQUET_OK) return status;

    writer->column_bloom_filter_overrides[column_index] = enabled;
    writer->column_bloom_ndv_overrides[column_index] = ndv;
    writer->column_bloom_fpp_overrides[column_index] = fpp;
    writer->column_bloom_options_set[column_index] = true;
    /* Ensure the finalize path actually emits the per-column bloom filter
     * even when the global option was left off. Additive: only flips the
     * flag on; never disables a globally-enabled configuration. */
    if (enabled) {
        writer->options.write_bloom_filters = true;
    }
    return CARQUET_OK;
}

carquet_status_t carquet_writer_set_sorting_columns(
    carquet_writer_t* writer,
    const carquet_sorting_column_t* columns,
    int32_t count) {

    /* writer is nonnull per API contract */
    if (count < 0 || (count > 0 && !columns)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (count > writer->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (count == 0) {
        writer->sorting_columns = NULL;
        writer->num_sorting_columns = 0;
        return CARQUET_OK;
    }

    parquet_sorting_column_t* copy = carquet_arena_calloc(&writer->arena,
        count, sizeof(parquet_sorting_column_t));
    if (!copy) return CARQUET_ERROR_OUT_OF_MEMORY;

    for (int32_t i = 0; i < count; i++) {
        if (columns[i].column_index < 0 ||
            columns[i].column_index >= writer->num_columns) {
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        copy[i].column_idx = columns[i].column_index;
        copy[i].descending = columns[i].descending;
        copy[i].nulls_first = columns[i].nulls_first;
    }

    writer->sorting_columns = copy;
    writer->num_sorting_columns = count;
    return CARQUET_OK;
}

/* ============================================================================
 * Writer Buffer API
 * ============================================================================
 */

carquet_writer_t* carquet_writer_create_buffer(
    const carquet_schema_t* schema,
    const carquet_writer_options_t* options,
    carquet_error_t* error) {

    /* Create a temporary FILE* for writing.
     * tmpfile() can fail on Windows when the process lacks write access to
     * the root directory.  Fall back to a named temp file in that case. */
    FILE* tmp = tmpfile();
#ifdef _WIN32
    if (!tmp) {
        /* _tempnam uses %TMP%, %TEMP%, or the current directory */
        char* tpath = _tempnam(NULL, "cqt");   /* CRT malloc: free with the CRT free() */
        if (tpath) {
            tmp = fopen(tpath, "w+bTD");   /* T=short-lived, D=delete-on-close */
            free(tpath);   /* never carquet_mem_free(): a custom allocator must not free CRT memory */
        }
    }
#endif
    if (!tmp) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN,
            "Failed to create temporary file for buffer writer");
        return NULL;
    }

    /* Create a writer using the FILE* handle */
    carquet_writer_t* writer = carquet_writer_create_file(tmp, schema, options, error);
    if (!writer) {
        fclose(tmp);
        return NULL;
    }

    /* Mark as buffer writer and take ownership of the tmpfile */
    writer->is_buffer_writer = true;
    writer->owns_file = true;

    return writer;
}

carquet_status_t carquet_writer_get_buffer(
    carquet_writer_t* writer,
    void** buffer,
    size_t* size) {

    /* writer, buffer, size are nonnull per API contract */
    if (!writer->is_buffer_writer) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (!writer->output_buffer || writer->output_buffer_size == 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Transfer ownership of the buffer to the caller */
    *buffer = writer->output_buffer;
    *size = writer->output_buffer_size;
    writer->output_buffer = NULL;
    writer->output_buffer_size = 0;

    /* Free the writer struct (close already freed internal resources) */
    carquet_arena_destroy(&writer->arena);
    carquet_mem_free(writer);

    return CARQUET_OK;
}
