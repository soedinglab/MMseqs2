/**
 * @file batch_reader.c
 * @brief High-level batch reader with column projection and parallel I/O
 *
 * This provides a production-ready API for efficiently reading Parquet files
 * with support for:
 * - Column projection (only read needed columns)
 * - Parallel column reading
 * - Memory-mapped I/O
 * - Batched output
 * - Buffer pooling to minimize allocations
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include "reader_internal.h"
#include "worker_pool.h"
#include "page_filter.h"
#include "core/arena.h"
#include "core/bitpack.h"
#include "core/endian.h"
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if !defined(_WIN32)
#include <sys/mman.h>
#include <unistd.h>   /* sysconf(_SC_PAGESIZE) */
#endif

/* SIMD dispatch function for null bitmap construction */
extern void carquet_dispatch_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                                int16_t max_def_level, uint8_t* null_bitmap);
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                int16_t max_def_level);

#define CARQUET_MAX_PAGE_PAYLOAD_SIZE (256ULL * 1024 * 1024)

/* Window handed to parquet_parse_page_header() when scanning a chunk's pages
 * in the mmap fast paths. Page headers carry per-page statistics, so a wide
 * BYTE_ARRAY column (parquet-mr writes untruncated min/max there) routinely
 * needs more than a few hundred bytes; a window that is too small silently
 * demotes those columns to the slow path. The parser stops at the struct's
 * STOP byte, so a generous window costs nothing on normal headers. */
#define CARQUET_PAGE_HEADER_SCAN_WINDOW ((size_t)64 * 1024)

/* Upper bound for the pipeline's per-slot, per-column buffer pre-allocation.
 * Sized from attacker-controlled row_group.num_rows; above this we fall back
 * to lazy allocation in pipeline_fill instead of eagerly malloc'ing. */
#define CARQUET_MAX_PREALLOC_BYTES (1024ULL * 1024 * 1024)

/* ============================================================================
 * Internal Structures
 * ============================================================================
 */

/* ============================================================================
 * Shared slot buffers (zero-copy Arrow export)
 * ============================================================================
 * A pipeline slot's value buffer and validity bitmap of one column, owned by
 * reference count: one reference held by the slot while it keeps the
 * buffers, one per exported Arrow array that points into them. The Arrow
 * exporter takes a reference instead of copying a column the pipeline
 * served (carquet_row_batch_column_share). When the slot is refilled it
 * drops its reference: if that was the last one it keeps the buffers,
 * otherwise it leaves them to the arrays (the last release frees them) and
 * allocates new ones. Releases may come from any thread, hence the atomic
 * count; everything else happens on the thread driving the batch reader.
 */
#if defined(__GNUC__) || defined(__clang__)
#define BR_SHARE_ATOMICS 1
static long br_share_inc(long* p) { return __atomic_add_fetch(p, 1, __ATOMIC_RELAXED); }
static long br_share_dec(long* p) { return __atomic_sub_fetch(p, 1, __ATOMIC_ACQ_REL); }
#elif defined(_MSC_VER)
#define BR_SHARE_ATOMICS 1
static long br_share_inc(long* p) { return _InterlockedIncrement((volatile long*)p); }
static long br_share_dec(long* p) { return _InterlockedDecrement((volatile long*)p); }
#else
#define BR_SHARE_ATOMICS 0   /* no portable atomics: exports copy */
#endif

typedef struct br_share {
    long refs;
    void* values;
    uint8_t* bitmap;
} br_share_t;

static void br_share_free(br_share_t* sh) {
    carquet_mem_free(sh->values);
    carquet_mem_free(sh->bitmap);
    carquet_mem_free(sh);
}

typedef struct carquet_column_data {
    void* data;                 /* Column values (or uint32_t* indices if dict preserved) */
    uint8_t* null_bitmap;       /* Null bitmap (1 bit per value), NULL for REQUIRED */
    int64_t num_values;         /* Number of values */
    size_t data_capacity;       /* Allocated capacity for data */
    carquet_physical_type_t type;
    int32_t type_length;        /* For fixed-length types */
    carquet_data_ownership_t ownership;  /* OWNED or VIEW (for future zero-copy) */

    /* Dictionary preservation (when config.preserve_dictionaries == true) */
    bool is_dictionary;             /* True if this column has preserved dictionary */
    const uint8_t* dictionary_data; /* Pointer to dictionary bytes (view, not owned) */
    int32_t dictionary_count;       /* Number of dictionary entries */
    const uint32_t* dictionary_offsets; /* Offset table for BYTE_ARRAY (view) */

    /* Nested (single-level LIST/MAP-leaf) reconstruction. When list_offsets is
     * non-NULL this column is a list: `data`/`null_bitmap`/`num_values` describe
     * the flattened child (element) array (Arrow child layout — values with a
     * validity bitmap), and list_offsets[i]..list_offsets[i+1] delimit list i in
     * that child array. list_validity (may be NULL) is the list-level null
     * bitmap. num_lists is the logical row count. */
    /* Pipeline-served column whose storage an Arrow export may reference
     * instead of copying (see br_share_t), or NULL. bitmap_shared: the
     * null_bitmap is a view into that storage (not a shifted copy). */
    br_share_t* share;
    bool bitmap_shared;

    int32_t* list_offsets;      /* [num_lists + 1] Arrow list offsets, or NULL */
    uint8_t* list_validity;     /* list-level validity bitmap (LSB, present=1), or NULL */
    int64_t  num_lists;         /* number of logical rows (lists) */
    int16_t  max_rep_level;     /* > 0 marks a repeated (list) column */
} carquet_column_data_t;

/* Pre-allocated column buffer pool for reuse across batches */
typedef struct carquet_column_pool {
    void* data;                 /* Pre-allocated data buffer */
    size_t data_capacity;       /* Capacity in bytes */
    uint8_t* null_bitmap;       /* Pre-allocated null bitmap */
    size_t bitmap_capacity;     /* Capacity in bytes */
    int16_t* def_levels;        /* Pre-allocated def levels buffer */
    size_t def_levels_capacity; /* Capacity in elements */
    /* Nested (list) reconstruction scratch/output buffers */
    int16_t* rep_levels;        /* Pre-allocated rep levels buffer */
    size_t rep_levels_capacity; /* Capacity in elements */
    int32_t* list_offsets;      /* Pre-allocated list offsets buffer */
    size_t list_offsets_capacity; /* Capacity in elements */
    uint8_t* list_validity;     /* Pre-allocated list-level validity bitmap */
    size_t list_validity_capacity; /* Capacity in bytes */
} carquet_column_pool_t;

struct carquet_row_batch {
    carquet_column_data_t* columns;
    int32_t num_columns;
    int64_t num_rows;
    carquet_arena_t arena;
    bool pooled;                /* If true, data buffers are from batch_reader pool */
    const int32_t* projection;  /* File leaf index per column (owned by the batch
                                 * reader); NULL outside carquet_batch_reader_next */
    carquet_worker_pool_t* pool; /* owner's column pool for a parallel Arrow
                                 * export (non-OpenMP builds), else NULL */
};

/* Pipeline ring buffer slot: holds pre-read column data for one RG */
typedef struct {
    carquet_column_reader_t** col_readers;  /* [num_projected] readers, used for bulk read */
    int32_t rg_index;                       /* row group index, -1 = empty */
    bool ready;                             /* all columns fully read */
    int32_t pending_tasks;                  /* bulk-read tasks still running for this
                                             * slot; guarded by the reader's slot_mutex */

    /* Pre-read value buffers (entire column chunk per column) */
    void** col_values;       /* [num_projected] value buffers */
    size_t* col_buf_sizes;   /* [num_projected] buffer capacities in bytes */
    int64_t* col_num_values; /* [num_projected] values actually read */
    int64_t total_rows;      /* total rows in this slot (= range total when filter active) */
    int64_t rows_consumed;   /* rows already served to batch_reader_next */
    /* Columns read as several page-range segments: the task_args index of the
     * first segment and the segment count (<= 1 means a single whole-chunk
     * task, whose count lands in col_num_values directly). */
    int32_t* col_task_first; /* [num_projected] */
    int32_t* col_task_count; /* [num_projected] */
    /* Validity bitmaps (LSB first, 1 = present) of the columns that have
     * definition levels, one bit per slot row; NULL for REQUIRED columns. */
    uint8_t** col_bitmaps;   /* [num_projected] */
    size_t* col_bitmap_sizes; /* [num_projected] capacities in bytes */
    /* Reference-counted ownership of col_values[i] / col_bitmaps[i] (see
     * br_share_t); NULL when not allocated (exports then copy). */
    br_share_t** col_share;  /* [num_projected] */

    /* Page-filter row ranges for this slot (only populated when a page
     * filter is active). Per-projected-column offset indexes are cached
     * here so worker tasks can seek to matching pages without going
     * through the file reader again. */
    carquet_row_range_list_t ranges;
    bool filter_ranges_valid;
    carquet_offset_index_t** col_offset_indexes;  /* [num_projected], may be NULL entries */
} rg_slot_t;

/* Forward decl shared with filtered bulk-read task (defined later). */
static int32_t find_page_for_row(
    const carquet_offset_index_t* oi,
    int64_t row_group_num_rows,
    int64_t target_row,
    int64_t* page_first_row_out);

/* Forward declarations for coalesced read fast path.
 * data_base: pointer to file data (per-slot mmap or shared mmap).
 *            Byte at file offset N is at data_base[N]. */
static bool can_coalesce_column(const carquet_column_reader_t* cr);
static void coalesced_read_column_range(const carquet_column_reader_t* cr,
    const uint8_t* data_base, void* dest, int64_t max_values,
    int64_t start_offset, int64_t end_offset, int64_t* out_values_read);
static bool coalesced_read_segment(const carquet_column_reader_t* cr,
    const uint8_t* data_base, void* dest, int64_t max_values,
    int64_t start_offset, int64_t end_offset,
    struct carquet_batch_reader* scratch_owner, int64_t* out_values_read);
static void coalesced_read_column(const carquet_column_reader_t* cr,
    void* dest, int64_t max_values, int64_t* out_values_read);
static int32_t plan_page_splits(const carquet_column_reader_t* cr,
    const uint8_t* data_base, int64_t max_values, int32_t max_splits,
    int64_t* split_offsets, int64_t* split_values);

extern uint32_t carquet_crc32(const uint8_t* data, size_t length);
extern carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data, size_t data_size, int32_t type_length,
    uint8_t* values, int64_t count);
extern carquet_status_t carquet_byte_stream_split_decode_float(
    const uint8_t* data, size_t data_size, float* values, int64_t count);
extern carquet_status_t carquet_byte_stream_split_decode_double(
    const uint8_t* data, size_t data_size, double* values, int64_t count);

/* Task argument for parallel bulk column reading */
typedef struct {
    carquet_column_reader_t* col_reader;
    const uint8_t* data_base;  /* file data pointer (per-slot or shared mmap) */
    void* dest;
    int64_t max_values;
    int64_t* out_values_read;
    int64_t start_offset;      /* 0 = full chunk */
    int64_t end_offset;        /* 0 = full chunk */
    int64_t local_values_read; /* scratch for split tasks */
    bool is_split;             /* one page-range segment of a column chunk:
                                * never falls back inline (see bulk_read_task_body) */

    /* Page-filter mode: when ranges is non-NULL the task reads only the
     * matching pages, writing rows contiguously into dest. */
    const carquet_row_range_list_t* ranges;
    const carquet_offset_index_t* offset_index;
    int64_t rg_num_rows;
    size_t value_size;

    /* Row-aligned read (see row_aligned_read). Used for columns with
     * definition levels (flat OPTIONAL: max_def > 0, max_rep == 0) and for
     * segments decoded through their own column reader. The task reads
     * dense values plus levels into dest, spreads them to their rows in
     * place and writes validity bits for rows [row_start, row_start +
     * max_values) into bitmap. A segment shares its first and last bitmap
     * bytes with its neighbours, so it leaves those in head/tail for the
     * consumer to merge once every segment has finished. */
    int16_t max_def;
    uint8_t* bitmap;
    int64_t row_start;
    int64_t head_byte;         /* -1 = none */
    int64_t tail_byte;         /* -1 = none */
    uint8_t head_bits;
    uint8_t tail_bits;
    int16_t* levels;           /* task-owned scratch, kept across fills */
    size_t levels_cap;
    uint8_t* vbits;            /* a segment's own validity bits (bit 0 = its
                                * first row), task-owned like levels */
    size_t vbits_cap;
    /* Segment of a chunk the coalesced path cannot decode (dictionary or
     * other encodings, levels): read through seg_reader, a column reader
     * owned by this task entry and rebound on each fill, positioned at
     * start_offset with carquet_column_reader_seek_to_data_page. */
    bool use_seg_reader;
    carquet_column_reader_t* seg_reader;

    /* Completion accounting: the owning reader's slot_mutex/slot_cond and the
     * slot's pending_tasks counter, decremented when the task finishes so
     * the consumer can wait for one slot instead of draining the pool. */
    struct carquet_batch_reader* owner;
    int32_t* slot_pending;
} bulk_read_arg_t;

struct carquet_batch_reader {
    carquet_reader_t* reader;
    carquet_batch_reader_config_t config;

    /* Detail for the most recent carquet_batch_reader_next() failure, exposed
     * through carquet_batch_reader_last_error(). next() has no error
     * out-parameter, so without this the message built deeper in the reader
     * (e.g. which codec check rejected a page, and at what file offset) would
     * be discarded and the caller would see only a status code. */
    carquet_error_t last_error;

    /* Column projection */
    int32_t* projected_columns;  /* File column indices to read */
    int32_t num_projected;       /* Number of projected columns */
    carquet_physical_type_t* projected_types;
    int32_t* projected_type_lengths;
    int16_t* projected_max_defs;
    int16_t* projected_max_reps;
    size_t* projected_value_sizes;
    bool has_repeated;           /* true if any projected column has max_rep > 0 */

    /* Reading state */
    int32_t current_row_group;
    int64_t rows_read_in_group;
    int64_t total_rows_read;

    /* Column readers for current row group */
    carquet_column_reader_t** col_readers;

    /* Memory-mapped data */
    uint8_t* mmap_data;
    size_t mmap_size;

    /* Buffer pool for reuse across batches (one per projected column) */
    carquet_column_pool_t* col_pools;

    /* Cached batch struct to avoid repeated alloc/free */
    carquet_row_batch_t* cached_batch;

    /* Persistent worker pool for cross-RG parallel decompression */
    carquet_worker_pool_t* pool;
    bool pool_is_borrowed;  /* true when pool comes from config.thread_pool */

    /* Column-parallel decode for builds without OpenMP: the per-batch column
     * loops run on this pool (created on first use; unused with OpenMP). */
    carquet_worker_pool_t* col_pool;
    struct br_col_task* col_tasks;  /* [num_projected] task slots */
    void** col_task_args;           /* [num_projected] pointers into col_tasks */
    bool col_pool_unavailable;      /* serial: 1 thread, or alloc failed */

    /* Pipeline ring buffer for multi-RG parallel decompression.
     * Pre-decompresses pages for upcoming row groups so that by the time
     * batch_reader_next() needs data, it's already decompressed. */
    rg_slot_t* pipeline;          /* [pipeline_depth] ring buffer */
    /* Per-slot completion: lets next() wait for the head slot only, so the
     * other slots keep decompressing in the background (previously a global
     * pool wait stalled every row-group transition on the newest fill). */
#ifdef _WIN32
    CRITICAL_SECTION slot_mutex;
    CONDITION_VARIABLE slot_cond;
#else
    pthread_mutex_t slot_mutex;
    pthread_cond_t slot_cond;
#endif
    bool slot_sync_init;
    /* Intra-chunk splitting: a column chunk may be cut into up to max_splits
     * page-aligned segments, each its own task, so a few large row groups
     * still keep every pool thread busy. Segment tasks of one column run
     * concurrently and therefore cannot share the column reader's
     * decompress_buffer: BYTE_STREAM_SPLIT staging comes from this small
     * pool instead (one entry per pool thread; guarded by slot_mutex). */
    int32_t max_splits;
    uint8_t** scratch_bufs;       /* [scratch_cap]; entries [0, scratch_free) are free */
    size_t* scratch_sizes;        /* [scratch_cap] capacities (0 = not yet allocated) */
    int32_t scratch_free;
    int32_t scratch_cap;
    int32_t pipeline_depth;       /* window size */
    int32_t pipeline_head;        /* next slot to consume */
    int32_t pipeline_count;       /* slots in use */
    int32_t* rg_order;            /* pre-filtered list of RG indices */
    int32_t rg_order_len;         /* total filtered RGs */
    int32_t rg_order_next;        /* next RG to submit */
    bool pipeline_active;         /* multi-RG pipeline enabled */
    /* First failure raised while filling or reading pipeline slots. Slots that
     * were read successfully before it are still served; once the failure is
     * returned from next() the pipeline is dead and keeps returning it, so a
     * failed or short read can never be mistaken for END_OF_DATA. */
    carquet_error_t pipeline_error;
    bool pipeline_dead;

    /* Per-reader task args (replaces static global array) */
    bulk_read_arg_t* task_args;
    int32_t task_args_capacity;

    /* ====================================================================
     * Page filter state
     * ==================================================================== */
    /* Active filter clauses (caller-owned; not copied). NULL = no filter. */
    const carquet_filter_clause_t* filter_clauses;
    int32_t filter_clause_count;

    /* Row ranges that survive the conjunction for current_row_group.
     * Valid only when filter_rg_state_valid is true. */
    carquet_row_range_list_t current_rg_ranges;
    bool filter_rg_state_valid;
    int32_t current_range_index;
    int64_t current_range_rows_emitted;
    bool range_positioned;          /* Column readers seeked to current range start? */
    int64_t rows_skipped;           /* Diagnostic accumulator */

    /* Per-projected-column offset index cache for the current row group.
     * Loaded lazily on first range positioning, freed on RG transition. */
    carquet_offset_index_t** projected_offset_indexes;
    int32_t projected_oi_rg;        /* -1 when cache is empty */

    /* ====================================================================
     * Row-range read state ([offset, limit) independent of any predicate).
     * Set via carquet_batch_reader_set_row_range(). Mutually exclusive with
     * the page filter and with repeated (nested) columns. Drives the
     * dedicated sequential batch_reader_next_row_range() path, bypassing the
     * pipeline. ==================================================================== */
    bool     row_range_active;
    int64_t  row_range_start;       /* global first row, inclusive */
    int64_t  row_range_end;         /* global exclusive end (INT64_MAX = to EOF) */
    int64_t  rr_rg_global_start;    /* global row index of current_row_group's row 0 */
    int64_t  rr_local_lo;           /* RG-local first row of the current window */
    int64_t  rr_local_hi;           /* RG-local exclusive end of the current window */
    int64_t  rr_emitted;            /* rows already emitted from the current window */
    bool     rr_rg_open;            /* readers opened + positioned for current window */
};

/* ============================================================================
 * Configuration
 * ============================================================================
 */

void carquet_batch_reader_config_init(carquet_batch_reader_config_t* config) {
    /* config is nonnull per API contract */
    memset(config, 0, sizeof(*config));
    config->batch_size = 65536;  /* 64K rows per batch */
    config->num_threads = 0;     /* Auto-detect */
    config->use_mmap = false;
}

/* ============================================================================
 * Helper Functions
 * ============================================================================
 */

/* Maximum reasonable type_length for FIXED_LEN_BYTE_ARRAY (16 MB) */
#define CARQUET_MAX_TYPE_LENGTH (16 * 1024 * 1024)

static size_t get_type_size(carquet_physical_type_t type, int32_t type_length) {
    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN: return 1;
        case CARQUET_PHYSICAL_INT32: return 4;
        case CARQUET_PHYSICAL_INT64: return 8;
        case CARQUET_PHYSICAL_INT96: return 12;
        case CARQUET_PHYSICAL_FLOAT: return 4;
        case CARQUET_PHYSICAL_DOUBLE: return 8;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            /* Validate type_length to prevent overflow attacks */
            if (type_length <= 0 || type_length > CARQUET_MAX_TYPE_LENGTH) {
                return 0;  /* Invalid - will cause allocation to fail safely */
            }
            return (size_t)type_length;
        case CARQUET_PHYSICAL_BYTE_ARRAY: return sizeof(carquet_byte_array_t);
        default: return 0;
    }
}

/* Column readers of this file can load pages from several threads at once:
 * the file is mapped, or read with positional I/O (see positional_io). */
static bool br_parallel_io(const carquet_reader_t* reader) {
    return reader->mmap_data != NULL || reader->positional_io;
}

static int resolve_column_name(const carquet_reader_t* reader, const char* name) {
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    if (!schema) return -1;

    return carquet_schema_find_column(schema, name);
}

/* Ensure a pool buffer is at least 'needed' bytes, growing if necessary */
static void* pool_ensure_data(carquet_column_pool_t* pool, size_t needed) {
    if (needed <= pool->data_capacity) {
        return pool->data;
    }
    carquet_mem_free(pool->data);
    pool->data = carquet_mem_malloc(needed);
    pool->data_capacity = pool->data ? needed : 0;
    return pool->data;
}

static uint8_t* pool_ensure_bitmap(carquet_column_pool_t* pool, size_t needed) {
    if (needed <= pool->bitmap_capacity) {
        memset(pool->null_bitmap, 0, needed);
        return pool->null_bitmap;
    }
    carquet_mem_free(pool->null_bitmap);
    pool->null_bitmap = carquet_mem_calloc(1, needed);
    pool->bitmap_capacity = pool->null_bitmap ? needed : 0;
    return pool->null_bitmap;
}

static int16_t* pool_ensure_def_levels(carquet_column_pool_t* pool, size_t count) {
    if (count <= pool->def_levels_capacity) {
        return pool->def_levels;
    }
    carquet_mem_free(pool->def_levels);
    pool->def_levels = carquet_mem_malloc(sizeof(int16_t) * count);
    pool->def_levels_capacity = pool->def_levels ? count : 0;
    return pool->def_levels;
}

static int16_t* pool_ensure_rep_levels(carquet_column_pool_t* pool, size_t count) {
    if (count <= pool->rep_levels_capacity) {
        return pool->rep_levels;
    }
    carquet_mem_free(pool->rep_levels);
    pool->rep_levels = carquet_mem_malloc(sizeof(int16_t) * count);
    pool->rep_levels_capacity = pool->rep_levels ? count : 0;
    return pool->rep_levels;
}

static int32_t* pool_ensure_list_offsets(carquet_column_pool_t* pool, size_t count) {
    if (count <= pool->list_offsets_capacity) {
        return pool->list_offsets;
    }
    carquet_mem_free(pool->list_offsets);
    pool->list_offsets = carquet_mem_malloc(sizeof(int32_t) * count);
    pool->list_offsets_capacity = pool->list_offsets ? count : 0;
    return pool->list_offsets;
}

static uint8_t* pool_ensure_list_validity(carquet_column_pool_t* pool, size_t bytes) {
    if (bytes <= pool->list_validity_capacity) {
        memset(pool->list_validity, 0, bytes);
        return pool->list_validity;
    }
    carquet_mem_free(pool->list_validity);
    pool->list_validity = carquet_mem_calloc(1, bytes);
    pool->list_validity_capacity = pool->list_validity ? bytes : 0;
    return pool->list_validity;
}

static bool column_can_zero_copy_batch(
    const carquet_column_reader_t* col_reader,
    carquet_physical_type_t type,
    int16_t max_def,
    int64_t rows_to_read) {

    if (!col_reader || !col_reader->page_loaded ||
        col_reader->decoded_ownership != CARQUET_DATA_VIEW ||
        max_def != 0 || col_reader->max_rep_level > 0 ||
        type == CARQUET_PHYSICAL_BYTE_ARRAY) {
        return false;
    }

    int32_t page_available = col_reader->page_num_values - col_reader->page_values_read;
    return page_available > 0 && (int64_t)page_available >= rows_to_read;
}

static int64_t column_zero_copy_rows_available(
    const carquet_column_reader_t* col_reader,
    carquet_physical_type_t type,
    int16_t max_def) {

    if (!col_reader || !col_reader->page_loaded ||
        col_reader->decoded_ownership != CARQUET_DATA_VIEW ||
        max_def != 0 || col_reader->max_rep_level > 0 ||
        type == CARQUET_PHYSICAL_BYTE_ARRAY) {
        return 0;
    }

    int32_t page_available = col_reader->page_num_values - col_reader->page_values_read;
    return page_available > 0 ? page_available : 0;
}

static int64_t clamp_rows_to_zero_copy_window(
    const carquet_batch_reader_t* batch_reader,
    int64_t rows_to_read) {

    int64_t zero_copy_rows = rows_to_read;

    for (int32_t i = 0; i < batch_reader->num_projected; i++) {
        int64_t page_rows = column_zero_copy_rows_available(
            batch_reader->col_readers[i],
            batch_reader->projected_types[i],
            batch_reader->projected_max_defs[i]);

        if (page_rows <= 0) {
            return rows_to_read;
        }
        if (page_rows < zero_copy_rows) {
            zero_copy_rows = page_rows;
        }
    }

    return zero_copy_rows;
}

static bool column_is_zero_copy_candidate(
    const carquet_column_reader_t* col_reader,
    carquet_physical_type_t type,
    int16_t max_def) {

    if (!col_reader || !col_reader->file_reader ||
        col_reader->file_reader->mmap_data == NULL ||
        !col_reader->col_meta ||
        max_def != 0 || col_reader->max_rep_level > 0 ||
        col_reader->col_meta->codec != CARQUET_COMPRESSION_UNCOMPRESSED) {
        return false;
    }

    switch (type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_INT96:
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE:
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return true;
        case CARQUET_PHYSICAL_BOOLEAN:
        case CARQUET_PHYSICAL_BYTE_ARRAY:
        default:
            return false;
    }
}

/* Record why a column read failed on the column reader itself. Safe inside the
 * parallel column loop: each thread owns one column reader. */
#define BR_COLUMN_FAIL(col_reader, read_error, status_code, ...)         \
    do {                                                                \
        CARQUET_SET_ERROR(&(col_reader)->last_error, (status_code), __VA_ARGS__); \
        *(read_error) = true;                                           \
    } while (0)

/* Word k (64 rows) of a validity bitmap covering `rows` rows; bits past
 * `rows` read as 0 and no byte past the bitmap is touched. */
static inline uint64_t bitmap_word(const uint8_t* bm, int64_t k, int64_t rows) {
    int64_t first = k * 64;
    int64_t left = rows - first;
    uint64_t w = 0;
    if (left >= 64) {
        w = carquet_read_u64_le(bm + (size_t)k * 8);
    } else {
        size_t nbytes = (size_t)((left + 7) / 8);
        for (size_t b = 0; b < nbytes; b++) {
            w |= (uint64_t)bm[(size_t)k * 8 + b] << (8 * b);
        }
        w &= ((uint64_t)1 << left) - 1;   /* 0 < left < 64 */
    }
    return w;
}

/* Move `count` values from index `from` to index `to` (to >= from, ranges
 * may overlap). */
static inline void move_values_back(uint8_t* data, int64_t to, int64_t from,
                                    int64_t count, size_t value_size) {
    memmove(data + (size_t)to * value_size, data + (size_t)from * value_size,
            (size_t)count * value_size);
}

/* Spread the rows [base, base + h) of one mixed bitmap word `w` (bit j =
 * row base + j), back to front and branch-free: row r takes dense value s
 * (masked to zero when absent) and s steps down only past a present row.
 * s = -1 once every present value below has been placed; the load index is
 * then clamped to 0 and the value masked off. Reading row r's source never
 * aliases a row already written: a present row's source index is <= r and
 * an absent row's is < r, while every row written so far is > r. Groups of
 * four rows do their loads before their stores, which keeps that order. A
 * value is LANES words of type T (a 16-byte value is two uint64_t). */
#define SPREAD_WORD_BODY(T, LANES)                                               \
    do {                                                                         \
        uint8_t* vals = data;                                                    \
        int j = h - 1;                                                           \
        for (; j >= 3; j -= 4) {                                                 \
            int64_t b0 = (int64_t)((w >> j) & 1), b1 = (int64_t)((w >> (j - 1)) & 1); \
            int64_t b2 = (int64_t)((w >> (j - 2)) & 1), b3 = (int64_t)((w >> (j - 3)) & 1); \
            int64_t s0 = s, s1 = s0 - b0, s2 = s1 - b1, s3 = s2 - b2;            \
            int64_t i0 = s0 & ~(s0 >> 63), i1 = s1 & ~(s1 >> 63);               \
            int64_t i2 = s2 & ~(s2 >> 63), i3 = s3 & ~(s3 >> 63);               \
            T m0 = (T)0 - (T)b0, m1 = (T)0 - (T)b1, m2 = (T)0 - (T)b2, m3 = (T)0 - (T)b3; \
            T v0[LANES], v1[LANES], v2[LANES], v3[LANES];                        \
            memcpy(v0, vals + (size_t)i0 * sizeof(v0), sizeof(v0));               \
            memcpy(v1, vals + (size_t)i1 * sizeof(v1), sizeof(v1));               \
            memcpy(v2, vals + (size_t)i2 * sizeof(v2), sizeof(v2));               \
            memcpy(v3, vals + (size_t)i3 * sizeof(v3), sizeof(v3));               \
            for (int l = 0; l < LANES; l++) {                                    \
                v0[l] &= m0; v1[l] &= m1; v2[l] &= m2; v3[l] &= m3;              \
            }                                                                    \
            memcpy(vals + (size_t)(base + j) * sizeof(v0), v0, sizeof(v0));      \
            memcpy(vals + (size_t)(base + j - 1) * sizeof(v1), v1, sizeof(v1));  \
            memcpy(vals + (size_t)(base + j - 2) * sizeof(v2), v2, sizeof(v2));  \
            memcpy(vals + (size_t)(base + j - 3) * sizeof(v3), v3, sizeof(v3));  \
            s = s3 - b3;                                                         \
        }                                                                        \
        for (; j >= 0; j--) {                                                    \
            int64_t b = (int64_t)((w >> j) & 1);                                 \
            int64_t i = s & ~(s >> 63);                                          \
            T m = (T)0 - (T)b, v[LANES];                                         \
            memcpy(v, vals + (size_t)i * sizeof(v), sizeof(v));                  \
            for (int l = 0; l < LANES; l++) v[l] &= m;                           \
            memcpy(vals + (size_t)(base + j) * sizeof(v), v, sizeof(v));         \
            s -= b;                                                              \
        }                                                                        \
    } while (0)

/* The column reader returns dense non-null values (Parquet convention); the
 * batch reader's contract is row-aligned: value[i] belongs to row i, with
 * null slots zeroed. Spread the dense values in place to the rows whose bit
 * is set in `bitmap` (LSB first, `rows` bits), back to front, one 64-row
 * word at a time: a full word moves in one block, an empty one is zeroed,
 * a mixed one goes through SPREAD_WORD_BODY for 4-, 8- and 16-byte values
 * and moves runs of present rows (found with clz) otherwise. It stops as soon
 * as every remaining row is present, those are already in place. The
 * present count comes from the bitmap, so the levels are not scanned again.
 */
static void spread_dense_by_bitmap(uint8_t* data, const uint8_t* bitmap,
                                   int64_t rows, size_t value_size) {
    if (rows <= 0) return;
    int64_t nwords = (rows + 63) / 64;
    int64_t src = 0;   /* present values = end of the dense prefix */
    for (int64_t k = 0; k < nwords; k++) {
        src += carquet_popcount64(bitmap_word(bitmap, k, rows));
    }
    if (src >= rows) return;
    for (int64_t k = nwords - 1; k >= 0; k--) {
        int64_t base = k * 64;
        int h = (int)(rows - base < 64 ? rows - base : 64);
        if (src == base + h) return;   /* rows [0, base + h) all present */
        uint64_t w = bitmap_word(bitmap, k, rows);
        if (h == 64 && w == ~(uint64_t)0) {
            src -= 64;
            memmove(data + (size_t)base * value_size, data + (size_t)src * value_size,
                    64 * value_size);
            continue;
        }
        if (w == 0) {
            memset(data + (size_t)base * value_size, 0, (size_t)h * value_size);
            continue;
        }
        int64_t s = src - 1;   /* dense index of the next value to place */
        if (value_size == 8) {
            SPREAD_WORD_BODY(uint64_t, 1);
            src = s + 1;
            continue;
        }
        if (value_size == 4) {
            SPREAD_WORD_BODY(uint32_t, 1);
            src = s + 1;
            continue;
        }
        if (value_size == 16) {
            SPREAD_WORD_BODY(uint64_t, 2);
            src = s + 1;
            continue;
        }
        while (h > 0) {
            /* Bit h-1 of w is the top bit of x; the low 64-h bits of x are 0. */
            uint64_t x = w << (64 - h);
            int run;
            if (x >> 63) {
                run = carquet_clz64(~x);   /* ~x != 0: its low 64-h bits are 1 */
                if (run > h) run = h;
                src -= run;
                move_values_back(data, base + h - run, src, run, value_size);
            } else {
                run = carquet_clz64(x);    /* 64 when x == 0 */
                if (run > h) run = h;
                memset(data + (size_t)(base + h - run) * value_size, 0,
                       (size_t)run * value_size);
            }
            h -= run;
        }
    }
}
#undef SPREAD_WORD_BODY

/* Validity bitmap of `rows` levels into `bitmap`, then the dense values of
 * `data` spread to their rows (see spread_dense_by_bitmap). */
static void levels_to_rows(uint8_t* data, const int16_t* def_levels, int64_t rows,
                           int16_t max_def, size_t value_size, uint8_t* bitmap) {
    carquet_dispatch_build_null_bitmap(def_levels, rows, max_def, bitmap);
    spread_dense_by_bitmap(data, bitmap, rows, value_size);
}

static void read_projected_column(
    carquet_batch_reader_t* batch_reader,
    carquet_row_batch_t* new_batch,
    int32_t col_i,
    int64_t rows_to_read,
    bool allow_zero_copy,
    bool* read_error) {

    if (*read_error) {
        return;
    }

    carquet_column_reader_t* col_reader = batch_reader->col_readers[col_i];
    carquet_column_data_t* col_data = &new_batch->columns[col_i];
    carquet_column_pool_t* pool = &batch_reader->col_pools[col_i];
    size_t value_size = batch_reader->projected_value_sizes[col_i];
    int16_t max_def = batch_reader->projected_max_defs[col_i];

    col_data->type = batch_reader->projected_types[col_i];
    col_data->type_length = batch_reader->projected_type_lengths[col_i];
    col_data->is_dictionary = false;
    col_data->dictionary_data = NULL;
    col_data->dictionary_count = 0;
    col_data->dictionary_offsets = NULL;

    /* Dictionary preservation was decided in reset_column_reader_for_row_group
     * (before any page load), so the decode/copy width is consistent. Just
     * read the resolved flag here. */
    bool use_dict_preserve = col_reader->preserve_dictionary;

    /* When preserving dictionaries, value_size is sizeof(uint32_t) for indices */
    size_t effective_value_size = use_dict_preserve ? sizeof(uint32_t) : value_size;

    /* Check if direct page handoff is possible:
     * - Column is REQUIRED (no nulls, no definition levels)
     * - Page loader can expose a stable view (mmap or reusable page buffer)
     * - Entire page slice fits in this batch
     * - Not in dictionary-preserve mode (indices layout differs)
     */
    bool try_zero_copy = allow_zero_copy &&
                         (max_def == 0) &&
                         (!col_reader->page_loaded) &&
                         !use_dict_preserve;

    if (try_zero_copy) {
        /* Trigger page load to check if it's a zero-copy page */
        int64_t dummy_read = carquet_column_read_batch(
            col_reader, NULL, 0, NULL, NULL);
        (void)dummy_read;
    }

    bool use_zero_copy = allow_zero_copy && !use_dict_preserve &&
                         column_can_zero_copy_batch(
                             col_reader, col_data->type, max_def, rows_to_read);

    if (use_zero_copy) {
        /* ====== ZERO-COPY PATH ====== */
        /* Point directly to the currently loaded page slice. */
        size_t byte_offset = (size_t)col_reader->page_values_read * value_size;
        col_data->data = (uint8_t*)col_reader->decoded_values + byte_offset;
        col_data->data_capacity = 0;  /* Not our allocation */
        col_data->ownership = CARQUET_DATA_VIEW;
        col_data->num_values = rows_to_read;

        /* No nulls in REQUIRED columns - return NULL bitmap */
        col_data->null_bitmap = NULL;

        /* Mark page as consumed */
        col_reader->page_values_read += (int32_t)rows_to_read;
        col_reader->values_remaining -= rows_to_read;
        return;
    }

    /* ====== STANDARD PATH (with copy, using pooled buffers) ====== */

    /* Validate value_size and check for overflow */
    if (effective_value_size == 0 || rows_to_read <= 0) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Column %d: unusable batch geometry (value_size=%llu, rows=%lld)",
            col_i, (unsigned long long)effective_value_size, (long long)rows_to_read);
        return;
    }

    /* Check for multiplication overflow (max 1GB allocation) */
    #define CARQUET_MAX_BATCH_ALLOC (1024ULL * 1024 * 1024)
    if (effective_value_size > CARQUET_MAX_BATCH_ALLOC / (size_t)rows_to_read) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: batch of %lld values x %llu bytes exceeds the 1 GiB cap",
            col_i, (long long)rows_to_read, (unsigned long long)effective_value_size);
        return;
    }

    size_t data_size = effective_value_size * (size_t)rows_to_read;

    /* Use pooled data buffer (grows as needed, never shrinks) */
    col_data->data = pool_ensure_data(pool, data_size);
    if (!col_data->data) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: failed to allocate %llu-byte value buffer",
            col_i, (unsigned long long)data_size);
        return;
    }
    col_data->data_capacity = data_size;
    col_data->ownership = CARQUET_DATA_VIEW;  /* Pool owns the buffer */

    /* Only allocate null bitmap for OPTIONAL columns */
    if (max_def > 0) {
        size_t bitmap_size = ((size_t)rows_to_read + 7) / 8;
        col_data->null_bitmap = pool_ensure_bitmap(pool, bitmap_size);
    } else {
        col_data->null_bitmap = NULL;  /* REQUIRED columns have no nulls */
    }

    /* Read values (reuse pooled def_levels buffer) */
    int16_t* def_levels = NULL;
    if (max_def > 0) {
        def_levels = pool_ensure_def_levels(pool, (size_t)rows_to_read);
    }

    /* _ex, not the plain form: it clears col_reader->last_error on entry and
     * propagates the page/codec-level message verbatim on failure, which is
     * what carquet_batch_reader_last_error() then hands back to the caller. */
    int64_t values_read = carquet_column_read_batch_ex(
        col_reader, col_data->data, rows_to_read, def_levels, NULL,
        &col_reader->last_error);

    /* _ex can return the values decoded before a failure together with the
     * error; a short count alone must not pass for a complete read. */
    if (values_read < 0 || col_reader->last_error.code != CARQUET_OK) {
        if (col_reader->last_error.code == CARQUET_OK) {
            CARQUET_SET_ERROR(&col_reader->last_error, CARQUET_ERROR_DECODE,
                "Column %d: read failed with no detail reported", col_i);
        }
        *read_error = true;
        return;
    }

    col_data->num_values = values_read;

    /* Attach dictionary metadata for preserved dictionary columns */
    if (use_dict_preserve) {
        col_data->is_dictionary = true;
        col_data->dictionary_data = col_reader->dictionary_data;
        col_data->dictionary_count = col_reader->dictionary_count;
        col_data->dictionary_offsets = col_reader->dictionary_offsets;
    }

    /* Validity bitmap from the definition levels, then the dense values
     * spread to their rows. */
    if (def_levels && col_data->null_bitmap && values_read > 0) {
        levels_to_rows((uint8_t*)col_data->data, def_levels, values_read, max_def,
                       effective_value_size, col_data->null_bitmap);
    }
}

#ifndef _OPENMP
/* ============================================================================
 * Column-parallel fallback without OpenMP
 * ============================================================================
 * The prefetch and read loops of batch_reader_next_impl() run one column per
 * thread under OpenMP. Builds without OpenMP (e.g. the macOS / Windows Python
 * wheels) run the same per-column work on a persistent worker pool instead,
 * under the same conditions (mmap-backed data, more than one column).
 */
typedef struct br_col_task {
    carquet_batch_reader_t* br;
    carquet_row_batch_t* batch;
    int32_t col;
    int64_t rows;
    bool error;
} br_col_task_t;

static int32_t br_cpu_count(void) {
#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (int32_t)si.dwNumberOfProcessors;
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int32_t)n : 1;
#endif
}

/* The reader's column pool, or NULL to run the loop serially. Sized to the
 * thread budget rather than the column count: the Arrow exporter also cuts
 * large string columns into row ranges and runs them on this pool. */
static carquet_worker_pool_t* br_column_pool(carquet_batch_reader_t* br) {
    if (br->col_pool) return br->col_pool;
    if (br->col_pool_unavailable) return NULL;
    int32_t n = br->config.num_threads;
    if (n <= 0) n = br_cpu_count();
    if (n > 64) n = 64;
    if (n < 2) { br->col_pool_unavailable = true; return NULL; }
    br->col_tasks = carquet_mem_calloc((size_t)br->num_projected, sizeof(br_col_task_t));
    br->col_task_args = carquet_mem_calloc((size_t)br->num_projected, sizeof(void*));
    br->col_pool = (br->col_tasks && br->col_task_args) ? carquet_worker_pool_create(n) : NULL;
    if (!br->col_pool) {
        carquet_mem_free(br->col_tasks);
        carquet_mem_free(br->col_task_args);
        br->col_tasks = NULL;
        br->col_task_args = NULL;
        br->col_pool_unavailable = true;
        return NULL;
    }
    for (int32_t i = 0; i < br->num_projected; i++) br->col_task_args[i] = &br->col_tasks[i];
    return br->col_pool;
}

static void br_prefetch_task(void* arg) {
    br_col_task_t* t = (br_col_task_t*)arg;
    carquet_column_reader_t* col_reader = t->br->col_readers[t->col];
    if (col_reader && !col_reader->page_loaded && col_reader->values_remaining > 0) {
        (void)carquet_column_read_batch(col_reader, NULL, 0, NULL, NULL);
    }
}

static void br_read_task(void* arg) {
    br_col_task_t* t = (br_col_task_t*)arg;
    read_projected_column(t->br, t->batch, t->col, t->rows, true, &t->error);
}

/* Run `fn` once per projected column on the column pool. */
static void br_column_parallel_for(carquet_batch_reader_t* br, carquet_task_fn fn,
                                   carquet_row_batch_t* batch, int64_t rows) {
    for (int32_t i = 0; i < br->num_projected; i++) {
        br_col_task_t* t = &br->col_tasks[i];
        t->br = br;
        t->batch = batch;
        t->col = i;
        t->rows = rows;
        t->error = false;
    }
    carquet_worker_pool_parallel_for(br->col_pool, fn, br->col_task_args, br->num_projected);
}
#endif /* !_OPENMP */

/* ============================================================================
 * Nested (single-level LIST) reconstruction
 * ============================================================================
 * Reconstructs an Arrow List<T> layout for a repeated column (max_rep == 1) by
 * reading the entire column chunk's leaf slots (values + def/rep levels) and
 * folding them into: a flattened child (element) array, its validity bitmap,
 * a list-offsets buffer, and a list-level validity bitmap.
 *
 * This handles the shapes produced by carquet_schema_add_list (and each MAP
 * leaf), i.e. a single repeated ancestor. Deeper nesting (max_rep > 1) is not
 * supported and sets *read_error.
 *
 * Definition-level bands for a single-level list (max_def = D):
 *   - a slot is an *element* of its list when def >= elem_exists (D-1 if the
 *     element leaf is OPTIONAL, else D);
 *   - the element's value is *present* when def == D (else it is a null element);
 *   - a rep==0 slot with def == 0 is a *null list* (only reachable when an
 *     optional ancestor sits above the repeated group).
 */
static void read_nested_list_column(
    carquet_batch_reader_t* batch_reader,
    carquet_row_batch_t* new_batch,
    int32_t col_i,
    int64_t expected_rows,
    bool* read_error) {

    if (*read_error) {
        return;
    }

    carquet_column_reader_t* col_reader = batch_reader->col_readers[col_i];
    carquet_column_data_t* col_data = &new_batch->columns[col_i];
    carquet_column_pool_t* pool = &batch_reader->col_pools[col_i];
    size_t value_size = batch_reader->projected_value_sizes[col_i];
    int16_t max_def = batch_reader->projected_max_defs[col_i];
    int16_t max_rep = batch_reader->projected_max_reps[col_i];

    col_data->type = batch_reader->projected_types[col_i];
    col_data->type_length = batch_reader->projected_type_lengths[col_i];
    col_data->max_rep_level = max_rep;

    /* Only single-level lists are supported in this release. */
    if (max_rep != 1 || value_size == 0 || max_def < 1) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "Column %d: only single-level lists are supported "
            "(max_rep=%d, max_def=%d, value_size=%llu)",
            col_i, (int)max_rep, (int)max_def, (unsigned long long)value_size);
        return;
    }

    /* Element-optional flag from the leaf node's own repetition. */
    const carquet_schema_t* schema = batch_reader->reader->schema;
    int32_t file_col = batch_reader->projected_columns[col_i];
    const parquet_schema_element_t* leaf = &schema->elements[schema->leaf_indices[file_col]];
    bool elem_optional = (leaf->repetition_type == CARQUET_REPETITION_OPTIONAL);
    int16_t elem_exists = elem_optional ? (int16_t)(max_def - 1) : max_def;

    /* Definition level at which the list itself is present. It is the max
     * definition level of the list's container (the group above the repeated
     * node): OPTIONAL/REPEATED ancestors each add one. A REQUIRED top-level
     * list has container_def == 0, so def == 0 is an *empty* list, not a null
     * one; an OPTIONAL top-level list has container_def == 1. */
    int16_t container_def = 0;
    {
        int32_t leaf_idx = schema->leaf_indices[file_col];
        int32_t rep_node = leaf_idx;
        while (rep_node > 0 &&
               schema->elements[rep_node].repetition_type != CARQUET_REPETITION_REPEATED) {
            rep_node = schema->parent_indices[rep_node];
        }
        int32_t node = rep_node > 0 ? schema->parent_indices[rep_node] : 0;
        for (; node > 0; node = schema->parent_indices[node]) {
            if (schema->elements[node].repetition_type != CARQUET_REPETITION_REQUIRED) {
                container_def++;
            }
        }
    }

    /* Total leaf slots in this chunk = number of (def, rep) entries. */
    int64_t total_slots = carquet_column_remaining(col_reader);
    if (total_slots < 0) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_INVALID_STATE,
            "Column %d: negative remaining leaf-slot count", col_i);
        return;
    }

    /* Bound allocations: values plus one def and one rep level per slot
     * (total_slots comes from the untrusted chunk metadata). */
    if (total_slots > 0 &&
        value_size + 2 * sizeof(int16_t) >
            CARQUET_MAX_BATCH_ALLOC / (size_t)total_slots) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: %lld leaf slots x %llu bytes exceeds the 1 GiB cap",
            col_i, (long long)total_slots, (unsigned long long)value_size);
        return;
    }

    size_t slots_alloc = total_slots > 0 ? (size_t)total_slots : 1;
    void* data = pool_ensure_data(pool, value_size * slots_alloc);
    int16_t* def_levels = pool_ensure_def_levels(pool, slots_alloc);
    int16_t* rep_levels = pool_ensure_rep_levels(pool, slots_alloc);
    if (!data || !def_levels || !rep_levels) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: failed to allocate nested-list scratch for %llu slots",
            col_i, (unsigned long long)slots_alloc);
        return;
    }

    int64_t slots = carquet_column_read_batch_ex(
        col_reader, data, total_slots, def_levels, rep_levels,
        &col_reader->last_error);
    if (slots < 0 || col_reader->last_error.code != CARQUET_OK) {
        if (col_reader->last_error.code == CARQUET_OK) {
            CARQUET_SET_ERROR(&col_reader->last_error, CARQUET_ERROR_DECODE,
                "Column %d: nested read failed with no detail reported", col_i);
        }
        *read_error = true;
        return;
    }

    /* Pass 1: count lists (rep==0), child elements (def >= elem_exists). */
    int64_t num_lists = 0, child_count = 0;
    for (int64_t j = 0; j < slots; j++) {
        if (rep_levels[j] == 0) num_lists++;
        if (def_levels[j] >= elem_exists) child_count++;
    }
    if (num_lists > INT32_MAX || child_count > INT32_MAX) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "Column %d: %lld lists / %lld elements exceed the 32-bit Arrow "
            "list-offset limit", col_i, (long long)num_lists, (long long)child_count);
        return;
    }
    (void)expected_rows;  /* num_lists is authoritative; equals the RG row count */

    col_data->data = data;
    col_data->data_capacity = value_size * slots_alloc;
    col_data->ownership = CARQUET_DATA_VIEW;
    col_data->num_values = child_count;
    col_data->num_lists = num_lists;

    /* Offsets buffer (num_lists + 1). */
    int32_t* offsets = pool_ensure_list_offsets(pool, (size_t)num_lists + 1);
    if (!offsets) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: failed to allocate %lld list offsets",
            col_i, (long long)num_lists + 1);
        return;
    }
    col_data->list_offsets = offsets;

    /* List-level validity: only materialized if some list is null. */
    uint8_t* list_valid = pool_ensure_list_validity(pool, ((size_t)num_lists + 7) / 8 + 1);
    if (!list_valid) {
        BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Column %d: failed to allocate list validity bitmap", col_i);
        return;
    }

    /* Child-element validity: only needed when the element can be null. */
    uint8_t* child_valid = NULL;
    if (elem_optional) {
        child_valid = pool_ensure_bitmap(pool, ((size_t)child_count + 7) / 8 + 1);
        if (!child_valid) {
            BR_COLUMN_FAIL(col_reader, read_error, CARQUET_ERROR_OUT_OF_MEMORY,
                "Column %d: failed to allocate element validity bitmap", col_i);
            return;
        }
    }
    col_data->null_bitmap = child_valid;

    /* Pass 2: build offsets + list validity. */
    int64_t li = -1, cc = 0;
    bool any_list_null = false;
    for (int64_t j = 0; j < slots; j++) {
        if (rep_levels[j] == 0) {
            li++;
            offsets[li] = (int32_t)cc;
            if (def_levels[j] >= container_def) {
                list_valid[li >> 3] |= (uint8_t)(1u << (li & 7));
            } else {
                any_list_null = true;
            }
        }
        if (def_levels[j] >= elem_exists) cc++;
    }
    offsets[num_lists] = (int32_t)cc;
    col_data->list_validity = any_list_null ? list_valid : NULL;

    /* Pass 3: expand dense (present-only) values into child-slot positions,
     * back-to-front so it can run in place, and build child validity. The
     * reader wrote `child_present` dense values at the front of `data`. */
    if (child_count > 0) {
        int64_t child_present = 0;
        for (int64_t j = 0; j < slots; j++) {
            if (def_levels[j] == max_def) child_present++;
        }
        uint8_t* bytes = (uint8_t*)data;
        /* child index for each element slot, walked back-to-front */
        int64_t ci = child_count - 1;
        int64_t src = child_present - 1;
        for (int64_t j = slots - 1; j >= 0; j--) {
            if (def_levels[j] < elem_exists) continue;  /* not an element */
            bool present = (def_levels[j] == max_def);
            uint8_t* dp = bytes + (size_t)ci * value_size;
            if (present) {
                uint8_t* sp = bytes + (size_t)src * value_size;
                if (sp != dp) memmove(dp, sp, value_size);
                src--;
                if (child_valid) {
                    child_valid[ci >> 3] |= (uint8_t)(1u << (ci & 7));
                }
            } else {
                memset(dp, 0, value_size);
            }
            ci--;
        }
    }
}

/* ============================================================================
 * Batch Reader Implementation
 * ============================================================================
 */

/* Memory the pipeline may hold across all of its slots. Each slot keeps a
 * whole row group of every projected column: the row-aligned values and
 * validity bitmap, plus for BYTE_ARRAY the page bytes its values point into
 * (retained pages or the dictionary, bounded by the chunk's uncompressed
 * size). Wide or string-heavy row groups therefore get fewer slots; the
 * segment split (see max_splits) keeps the pool busy with fewer slots. */
#define CARQUET_PIPELINE_MEMORY_BUDGET (256LL * 1024 * 1024)

/* Most pipeline slots that fit CARQUET_PIPELINE_MEMORY_BUDGET for the
 * largest selected row group (at least 2). Metadata sizes are untrusted:
 * every sum saturates. */
static int32_t pipeline_memory_depth(const carquet_batch_reader_t* br) {
    const parquet_file_metadata_t* meta = &br->reader->metadata;
    int64_t worst = 0;
    for (int32_t r = 0; r < br->rg_order_len; r++) {
        const parquet_row_group_t* rg = &meta->row_groups[br->rg_order[r]];
        int64_t rows = rg->num_rows > 0 ? rg->num_rows : 0;
        int64_t bytes = 0;
        for (int32_t i = 0; i < br->num_projected; i++) {
            int64_t col = 0;
            size_t vs = br->projected_value_sizes[i];
            if (vs > 0 && rows <= INT64_MAX / (int64_t)vs) col = rows * (int64_t)vs;
            else col = INT64_MAX / 4;
            if (br->projected_max_defs[i] > 0) col += rows / 8 + 1;
            int32_t fc = br->projected_columns[i];
            if (br->projected_types[i] == CARQUET_PHYSICAL_BYTE_ARRAY &&
                rg->columns && fc < rg->num_columns && rg->columns[fc].has_metadata) {
                int64_t u = rg->columns[fc].metadata.total_uncompressed_size;
                if (u > 0) col += u < INT64_MAX / 4 ? u : INT64_MAX / 4;
            }
            bytes = (bytes > INT64_MAX / 2 - col) ? INT64_MAX / 2 : bytes + col;
        }
        if (bytes > worst) worst = bytes;
    }
    if (worst <= 0) return INT32_MAX;
    int64_t d = CARQUET_PIPELINE_MEMORY_BUDGET / worst;
    if (d < 2) d = 2;
    return d > INT32_MAX ? INT32_MAX : (int32_t)d;
}

carquet_batch_reader_t* carquet_batch_reader_create(
    carquet_reader_t* reader,
    const carquet_batch_reader_config_t* config,
    carquet_error_t* error) {

    /* reader is nonnull per API contract */
    carquet_batch_reader_t* batch_reader = carquet_mem_calloc(1, sizeof(carquet_batch_reader_t));
    if (!batch_reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate batch reader");
        return NULL;
    }

    batch_reader->reader = reader;

    /* Copy config or use defaults */
    if (config) {
        batch_reader->config = *config;
    } else {
        carquet_batch_reader_config_init(&batch_reader->config);
    }
    /* A zero-initialised config (not config_init) has batch_size 0, which
     * the sequential and pipeline paths would turn into endless empty
     * batches (and a negative size into a wrapped buffer offset). */
    if (batch_reader->config.batch_size <= 0) {
        batch_reader->config.batch_size = 65536;
    }

    /* Resolve column projection */
    int32_t total_columns = carquet_reader_num_columns(reader);

    if (batch_reader->config.column_indices && batch_reader->config.num_columns > 0) {
        /* Use provided column indices */
        batch_reader->num_projected = batch_reader->config.num_columns;
        batch_reader->projected_columns = carquet_mem_malloc(sizeof(int32_t) * batch_reader->num_projected);
        if (!batch_reader->projected_columns) {
            carquet_mem_free(batch_reader);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate projection");
            return NULL;
        }
        memcpy(batch_reader->projected_columns, batch_reader->config.column_indices,
               sizeof(int32_t) * batch_reader->num_projected);
        /* Caller-supplied indices go straight into schema->leaf_indices[]. */
        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            int32_t idx = batch_reader->projected_columns[i];
            if (idx < 0 || idx >= total_columns) {
                carquet_mem_free(batch_reader->projected_columns);
                carquet_mem_free(batch_reader);
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Projected column index %d out of range (file has %d columns)",
                    idx, total_columns);
                return NULL;
            }
        }
    } else if (batch_reader->config.column_names && batch_reader->config.num_column_names > 0) {
        /* Resolve column names to indices */
        batch_reader->num_projected = batch_reader->config.num_column_names;
        batch_reader->projected_columns = carquet_mem_malloc(sizeof(int32_t) * batch_reader->num_projected);
        if (!batch_reader->projected_columns) {
            carquet_mem_free(batch_reader);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate projection");
            return NULL;
        }

        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            const char* col_name = batch_reader->config.column_names[i];
            int32_t idx = resolve_column_name(reader, col_name);
            if (idx < 0) {
                carquet_mem_free(batch_reader->projected_columns);
                carquet_mem_free(batch_reader);
                CARQUET_SET_ERROR(error, CARQUET_ERROR_COLUMN_NOT_FOUND,
                    "Column not found: %s", col_name);
                return NULL;
            }
            batch_reader->projected_columns[i] = idx;
        }
    } else {
        /* Read all columns */
        batch_reader->num_projected = total_columns;
        batch_reader->projected_columns = carquet_mem_malloc(sizeof(int32_t) * total_columns);
        if (!batch_reader->projected_columns) {
            carquet_mem_free(batch_reader);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate projection");
            return NULL;
        }
        for (int32_t i = 0; i < total_columns; i++) {
            batch_reader->projected_columns[i] = i;
        }
    }

    /* Allocate column reader array */
    batch_reader->col_readers = carquet_mem_calloc(batch_reader->num_projected,
                                        sizeof(carquet_column_reader_t*));
    if (!batch_reader->col_readers) {
        carquet_mem_free(batch_reader->projected_columns);
        carquet_mem_free(batch_reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate column readers");
        return NULL;
    }

    batch_reader->projected_types = carquet_mem_malloc(sizeof(carquet_physical_type_t) *
                                           (size_t)batch_reader->num_projected);
    batch_reader->projected_type_lengths = carquet_mem_malloc(sizeof(int32_t) *
                                                  (size_t)batch_reader->num_projected);
    batch_reader->projected_max_defs = carquet_mem_malloc(sizeof(int16_t) *
                                              (size_t)batch_reader->num_projected);
    batch_reader->projected_max_reps = carquet_mem_malloc(sizeof(int16_t) *
                                              (size_t)batch_reader->num_projected);
    batch_reader->projected_value_sizes = carquet_mem_malloc(sizeof(size_t) *
                                                 (size_t)batch_reader->num_projected);
    if (!batch_reader->projected_types || !batch_reader->projected_type_lengths ||
        !batch_reader->projected_max_defs || !batch_reader->projected_max_reps ||
        !batch_reader->projected_value_sizes) {
        carquet_mem_free(batch_reader->projected_value_sizes);
        carquet_mem_free(batch_reader->projected_max_reps);
        carquet_mem_free(batch_reader->projected_max_defs);
        carquet_mem_free(batch_reader->projected_type_lengths);
        carquet_mem_free(batch_reader->projected_types);
        carquet_mem_free(batch_reader->col_readers);
        carquet_mem_free(batch_reader->projected_columns);
        carquet_mem_free(batch_reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate projection metadata");
        return NULL;
    }

    {
        const carquet_schema_t* schema = carquet_reader_schema(reader);
        batch_reader->has_repeated = false;
        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            int32_t file_col_idx = batch_reader->projected_columns[i];
            int32_t schema_idx = schema->leaf_indices[file_col_idx];
            const parquet_schema_element_t* elem = &schema->elements[schema_idx];

            batch_reader->projected_types[i] =
                elem->has_type ? elem->type : CARQUET_PHYSICAL_BYTE_ARRAY;
            batch_reader->projected_type_lengths[i] = elem->type_length;
            batch_reader->projected_max_defs[i] = schema->max_def_levels[file_col_idx];
            batch_reader->projected_max_reps[i] = schema->max_rep_levels[file_col_idx];
            if (batch_reader->projected_max_reps[i] > 0) {
                batch_reader->has_repeated = true;
            }
            batch_reader->projected_value_sizes[i] = get_type_size(
                batch_reader->projected_types[i],
                batch_reader->projected_type_lengths[i]);
        }
    }

    /* Allocate buffer pool (one per projected column) */
    batch_reader->col_pools = carquet_mem_calloc(batch_reader->num_projected,
                                      sizeof(carquet_column_pool_t));
    if (!batch_reader->col_pools) {
        carquet_mem_free(batch_reader->projected_value_sizes);
        carquet_mem_free(batch_reader->projected_max_reps);
        carquet_mem_free(batch_reader->projected_max_defs);
        carquet_mem_free(batch_reader->projected_type_lengths);
        carquet_mem_free(batch_reader->projected_types);
        carquet_mem_free(batch_reader->col_readers);
        carquet_mem_free(batch_reader->projected_columns);
        carquet_mem_free(batch_reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate buffer pool");
        return NULL;
    }

    batch_reader->current_row_group = -1;

    /* Filter state init */
    carquet_row_range_list_init(&batch_reader->current_rg_ranges);
    batch_reader->filter_rg_state_valid = false;
    batch_reader->current_range_index = 0;
    batch_reader->current_range_rows_emitted = 0;
    batch_reader->range_positioned = false;
    batch_reader->rows_skipped = 0;
    batch_reader->projected_offset_indexes = NULL;
    batch_reader->projected_oi_rg = -1;

    /* ====================================================================
     * Pre-compute filtered row group order
     * ==================================================================== */
    int32_t num_row_groups = carquet_reader_num_row_groups(reader);
    batch_reader->rg_order = carquet_mem_malloc(sizeof(int32_t) * (size_t)num_row_groups);
    if (!batch_reader->rg_order) {
        batch_reader->rg_order_len = 0;
    } else {
        int32_t count = 0;
        for (int32_t rg = 0; rg < num_row_groups; rg++) {
            if (batch_reader->config.row_group_filter) {
                bool should_read = batch_reader->config.row_group_filter(
                    reader, rg, batch_reader->config.row_group_filter_ctx);
                if (!should_read) continue;
            }
            batch_reader->rg_order[count++] = rg;
        }
        batch_reader->rg_order_len = count;
    }
    batch_reader->rg_order_next = 0;

    /* ====================================================================
     * Create worker pool + pipeline for compressed mmap multi-RG files
     * ==================================================================== */
    /* Enable the pipeline for multi-RG files, or single-RG files that are
     * large enough for the pipeline overhead to be amortized.  Small files
     * (< 500K rows) are faster with the simpler per-batch OMP path. */
    /* row_groups[].num_rows is parsed straight from attacker-controlled
     * metadata, so a crafted file can supply negative or absurd values.
     * Saturate the accumulation (this total only gates a >= 500000 check)
     * and ignore non-positive counts rather than overflowing int64_t. */
    int64_t total_pipeline_rows = 0;
    for (int32_t r = 0; r < batch_reader->rg_order_len; r++) {
        int64_t rg_rows = reader->metadata.row_groups[batch_reader->rg_order[r]].num_rows;
        if (rg_rows <= 0) continue;
        if (rg_rows > INT64_MAX - total_pipeline_rows) {
            total_pipeline_rows = INT64_MAX;
            break;
        }
        total_pipeline_rows += rg_rows;
    }
    bool pipeline_safe = true;
    /* The pipeline sizes one buffer per column from row_groups[].num_rows up
     * front. A row count beyond file_size*8 is legitimate (RLE runs, constant
     * data under any codec) but is also what a crafted file uses to drive a
     * huge allocation, so such files take the sequential path, whose buffers
     * are bounded by batch_size, instead of being rejected. */
    for (int32_t r = 0; r < batch_reader->rg_order_len; r++) {
        int64_t rg_rows = reader->metadata.row_groups[batch_reader->rg_order[r]].num_rows;
        if (rg_rows < 0 ||
            (uint64_t)rg_rows > (uint64_t)reader->file_size * 8u) {
            pipeline_safe = false;
            break;
        }
    }
    /* Repeated columns have their own driver. Columns with definition levels
     * (flat OPTIONAL) are pre-read row-aligned with a validity bitmap (see
     * row_aligned_read). BYTE_ARRAY values point into the dictionary or the
     * retained pages of the column reader that decoded them (the slot's own
     * reader, or a segment's seg_reader); both belong to the slot and are
     * only rebound when the slot is refilled, after its last batch has been
     * served. Slots hold physical values only, so a reader asked to
     * preserve dictionaries reads dictionary-encoded columns sequentially,
     * where the indices are exposed. */
    for (int32_t ci = 0; pipeline_safe && ci < batch_reader->num_projected; ci++) {
        int32_t file_col = batch_reader->projected_columns[ci];
        if (reader->schema->max_rep_levels[file_col] > 0) {
            pipeline_safe = false;
        } else if (reader->schema->max_def_levels[file_col] > 0 &&
                   batch_reader->projected_value_sizes[ci] == 0) {
            pipeline_safe = false;
        }
        for (int32_t r = 0; pipeline_safe && batch_reader->config.preserve_dictionaries &&
                            r < batch_reader->rg_order_len; r++) {
            const parquet_row_group_t* rgm = &reader->metadata.row_groups[batch_reader->rg_order[r]];
            if (rgm->columns && file_col < rgm->num_columns &&
                rgm->columns[file_col].has_metadata &&
                rgm->columns[file_col].metadata.has_dictionary_page_offset) {
                pipeline_safe = false;
            }
        }
    }

    if (pipeline_safe &&
        br_parallel_io(reader) &&
        (batch_reader->rg_order_len > 1 || total_pipeline_rows >= 500000)) {
        /* Worth pre-reading when a projected column needs real decoding:
         * a codec, definition levels or a dictionary. Uncompressed PLAIN
         * REQUIRED columns are served as page views by the sequential path. */
        bool needs_decode = false;
        const parquet_file_metadata_t* meta = &reader->metadata;
        if (meta->num_row_groups > 0 && meta->row_groups[0].columns) {
            for (int32_t ci = 0; ci < batch_reader->num_projected; ci++) {
                int32_t file_col = batch_reader->projected_columns[ci];
                if (batch_reader->projected_max_defs[ci] > 0) {
                    needs_decode = true;
                    break;
                }
                if (file_col < meta->row_groups[0].num_columns) {
                    const parquet_column_chunk_t* chunk = &meta->row_groups[0].columns[file_col];
                    if (chunk->has_metadata &&
                        (chunk->metadata.codec != CARQUET_COMPRESSION_UNCOMPRESSED ||
                         chunk->metadata.has_dictionary_page_offset)) {
                        needs_decode = true;
                        break;
                    }
                }
            }
        }

        if (needs_decode) {
            /* Optimization 6: borrow external pool if provided */
            if (batch_reader->config.thread_pool) {
                batch_reader->pool = (carquet_worker_pool_t*)batch_reader->config.thread_pool;
                batch_reader->pool_is_borrowed = true;
            } else {
                int32_t pt = batch_reader->config.num_threads;
#ifdef _OPENMP
                if (pt <= 0) pt = omp_get_max_threads();
#else
                /* As the OpenMP default (one thread per CPU): measured on an
                 * 8-CPU Mac, 8 pipeline threads read a 10M-row nullable file
                 * ~25% and a 10M-row strings file ~30% faster than 4. */
                if (pt <= 0) pt = br_cpu_count();
                if (pt > 64) pt = 64;
#endif
                if (pt < 2) pt = 2;

                batch_reader->pool = carquet_worker_pool_create(pt);
                batch_reader->pool_is_borrowed = false;
            }

            if (batch_reader->pool) {
                int32_t pt = batch_reader->pool->num_threads;
                int32_t depth = batch_reader->rg_order_len;
                if (depth > pt * 2) depth = pt * 2;
                int32_t mem_depth = pipeline_memory_depth(batch_reader);
                if (depth > mem_depth) depth = mem_depth;
                if (depth < 1) depth = 1;

                batch_reader->pipeline = carquet_mem_calloc(depth, sizeof(rg_slot_t));
                if (batch_reader->pipeline) {
                    batch_reader->pipeline_depth = depth;
#ifdef _WIN32
                    InitializeCriticalSection(&batch_reader->slot_mutex);
                    InitializeConditionVariable(&batch_reader->slot_cond);
                    batch_reader->slot_sync_init = true;
#else
                    batch_reader->slot_sync_init =
                        pthread_mutex_init(&batch_reader->slot_mutex, NULL) == 0;
                    if (batch_reader->slot_sync_init &&
                        pthread_cond_init(&batch_reader->slot_cond, NULL) != 0) {
                        pthread_mutex_destroy(&batch_reader->slot_mutex);
                        batch_reader->slot_sync_init = false;
                    }
#endif
                    batch_reader->pipeline_head = 0;
                    batch_reader->pipeline_count = 0;
                    batch_reader->pipeline_active = true;

                    /* Allocate per-reader task args. Each projected column of
                     * a slot is one task, or up to max_splits page-range
                     * segment tasks when the pipeline window would otherwise
                     * offer the pool too few tasks (few row groups and/or
                     * columns). Aim for ~32 tasks per thread per window. */
                    int32_t np = batch_reader->num_projected;
                    int32_t max_splits = 1;
                    if (batch_reader->slot_sync_init && np > 0) {
                        int64_t budget = 32 * (int64_t)pt;
                        int64_t base_tasks = (int64_t)depth * np;
                        int64_t ms = budget / base_tasks;
                        if (ms > 16) ms = 16;
                        if (ms < 1) ms = 1;
                        max_splits = (int32_t)ms;
                    }
                    batch_reader->max_splits = max_splits;
                    if (max_splits > 1) {
                        batch_reader->scratch_bufs = carquet_mem_calloc(pt, sizeof(uint8_t*));
                        batch_reader->scratch_sizes = carquet_mem_calloc(pt, sizeof(size_t));
                        if (batch_reader->scratch_bufs && batch_reader->scratch_sizes) {
                            batch_reader->scratch_cap = pt;
                            batch_reader->scratch_free = pt;
                        } else {
                            carquet_mem_free(batch_reader->scratch_bufs);
                            carquet_mem_free(batch_reader->scratch_sizes);
                            batch_reader->scratch_bufs = NULL;
                            batch_reader->scratch_sizes = NULL;
                            batch_reader->max_splits = max_splits = 1;
                        }
                    }
                    batch_reader->task_args_capacity = depth * np * (max_splits + 1);
                    batch_reader->task_args = carquet_mem_calloc(batch_reader->task_args_capacity,
                                                     sizeof(bulk_read_arg_t));
                    if (!batch_reader->task_args) {
                        carquet_mem_free(batch_reader->pipeline);
                        batch_reader->pipeline = NULL;
                        batch_reader->pipeline_active = false;
                        batch_reader->task_args_capacity = 0;
                    }

                    bool alloc_ok = batch_reader->pipeline_active;
                    for (int32_t s = 0; s < depth && alloc_ok; s++) {
                        batch_reader->pipeline[s].col_readers = carquet_mem_calloc(np, sizeof(carquet_column_reader_t*));
                        batch_reader->pipeline[s].col_values = carquet_mem_calloc(np, sizeof(void*));
                        batch_reader->pipeline[s].col_buf_sizes = carquet_mem_calloc(np, sizeof(size_t));
                        batch_reader->pipeline[s].col_num_values = carquet_mem_calloc(np, sizeof(int64_t));
                        batch_reader->pipeline[s].col_task_first = carquet_mem_calloc(np, sizeof(int32_t));
                        batch_reader->pipeline[s].col_task_count = carquet_mem_calloc(np, sizeof(int32_t));
                        batch_reader->pipeline[s].col_bitmaps = carquet_mem_calloc(np, sizeof(uint8_t*));
                        batch_reader->pipeline[s].col_bitmap_sizes = carquet_mem_calloc(np, sizeof(size_t));
                        batch_reader->pipeline[s].col_share = carquet_mem_calloc(np, sizeof(br_share_t*));
                        batch_reader->pipeline[s].rg_index = -1;
                        if (!batch_reader->pipeline[s].col_readers ||
                            !batch_reader->pipeline[s].col_values ||
                            !batch_reader->pipeline[s].col_buf_sizes ||
                            !batch_reader->pipeline[s].col_num_values ||
                            !batch_reader->pipeline[s].col_task_first ||
                            !batch_reader->pipeline[s].col_task_count ||
                            !batch_reader->pipeline[s].col_bitmaps ||
                            !batch_reader->pipeline[s].col_bitmap_sizes ||
                            !batch_reader->pipeline[s].col_share) {
                            /* Cleanup all slots on failure */
                            for (int32_t j = 0; j <= s; j++) {
                                carquet_mem_free(batch_reader->pipeline[j].col_readers);
                                carquet_mem_free(batch_reader->pipeline[j].col_values);
                                carquet_mem_free(batch_reader->pipeline[j].col_buf_sizes);
                                carquet_mem_free(batch_reader->pipeline[j].col_num_values);
                                carquet_mem_free(batch_reader->pipeline[j].col_task_first);
                                carquet_mem_free(batch_reader->pipeline[j].col_task_count);
                                carquet_mem_free(batch_reader->pipeline[j].col_bitmaps);
                                carquet_mem_free(batch_reader->pipeline[j].col_bitmap_sizes);
                                carquet_mem_free(batch_reader->pipeline[j].col_share);
                            }
                            carquet_mem_free(batch_reader->pipeline);
                            batch_reader->pipeline = NULL;
                            batch_reader->pipeline_active = false;
                            alloc_ok = false;
                        }
                    }

                    /* Optimization 3: pre-allocate value buffers based on
                     * max row group size so pipeline_fill avoids malloc
                     * storms on the hot path. */
                    if (alloc_ok) {
                        int64_t max_rg_rows = 0;
                        for (int32_t r = 0; r < batch_reader->rg_order_len; r++) {
                            int64_t rr = meta->row_groups[batch_reader->rg_order[r]].num_rows;
                            if (rr > max_rg_rows) max_rg_rows = rr;
                        }
                        /* This is a hot-path optimization only: a NULL/zero
                         * slot buffer is lazily (re)allocated by pipeline_fill.
                         * num_rows is attacker-controlled metadata, so skip the
                         * pre-allocation entirely when the size is non-positive,
                         * overflows size_t, or is implausibly large rather than
                         * trusting it into carquet_mem_malloc(). */
                        if (max_rg_rows > 0) {
                            for (int32_t s = 0; s < depth; s++) {
                                rg_slot_t* slot = &batch_reader->pipeline[s];
                                for (int32_t i = 0; i < np; i++) {
                                    size_t vsize = batch_reader->projected_value_sizes[i];
                                    if (vsize == 0 ||
                                        (uint64_t)max_rg_rows > SIZE_MAX / vsize) {
                                        continue;
                                    }
                                    size_t needed = (size_t)max_rg_rows * vsize;
                                    if (needed == 0 ||
                                        needed > CARQUET_MAX_PREALLOC_BYTES) {
                                        continue;
                                    }
                                    slot->col_values[i] = carquet_mem_malloc(needed);
                                    slot->col_buf_sizes[i] =
                                        slot->col_values[i] ? needed : 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return batch_reader;
}

static carquet_status_t open_row_group_readers(
    carquet_batch_reader_t* batch_reader,
    int32_t row_group_index,
    carquet_error_t* error) {

    /* Reuse existing readers if possible, otherwise create new ones */
    for (int32_t i = 0; i < batch_reader->num_projected; i++) {
        int32_t file_col_idx = batch_reader->projected_columns[i];

        if (batch_reader->col_readers[i]) {
            /* Reuse: reset state but keep allocated buffers */
            carquet_column_reader_rebind(
                batch_reader->col_readers[i],
                batch_reader->reader,
                row_group_index, file_col_idx,
                batch_reader->config.preserve_dictionaries);
        } else {
            /* First time: create new reader */
            batch_reader->col_readers[i] = carquet_reader_get_column(
                batch_reader->reader, row_group_index, file_col_idx, error);

            if (!batch_reader->col_readers[i]) {
                /* Close already opened readers */
                for (int32_t j = 0; j < i; j++) {
                    carquet_column_reader_free(batch_reader->col_readers[j]);
                    batch_reader->col_readers[j] = NULL;
                }
                return error ? error->code : CARQUET_ERROR_COLUMN_NOT_FOUND;
            }
        }

        /* Decide dictionary preservation up front for BOTH the reset and the
         * freshly-created reader, before any page is pre-loaded. The decode
         * buffer width depends on this (physical value vs. 4-byte index), and
         * pages are pre-loaded ahead of the read phase, so it must be fixed
         * now rather than at read time. */
        carquet_column_reader_t* cr = batch_reader->col_readers[i];
        cr->preserve_dictionary = batch_reader->config.preserve_dictionaries &&
            cr->col_meta && cr->col_meta->has_dictionary_page_offset &&
            cr->max_rep_level == 0;  /* nested columns can't expose preserved indices */
    }

    batch_reader->current_row_group = row_group_index;
    batch_reader->rows_read_in_group = 0;

    return CARQUET_OK;
}

/* ============================================================================
 * Pipeline Ring Buffer
 * ============================================================================
 *
 * Pre-decompresses pages for multiple row groups in parallel using the worker
 * pool. On the first call to batch_reader_next(), decompression tasks are
 * submitted for up to pipeline_depth row groups. As the user consumes
 * batches and exhausts a row group, that slot is retired and the next
 * uncovered RG is submitted. With pipeline_depth >= total_RGs (the common
 * benchmark case), ALL decompression happens upfront in parallel.
 */

/**
 * Bulk-read task: reads ALL values from a column reader into a pre-allocated
 * buffer. This forces decompression of ALL pages in the column chunk.
 */
static void bulk_read_task_body(bulk_read_arg_t* t);

static void slot_sync_lock(carquet_batch_reader_t* br) {
#ifdef _WIN32
    EnterCriticalSection(&br->slot_mutex);
#else
    pthread_mutex_lock(&br->slot_mutex);
#endif
}
static void slot_sync_unlock(carquet_batch_reader_t* br) {
#ifdef _WIN32
    LeaveCriticalSection(&br->slot_mutex);
#else
    pthread_mutex_unlock(&br->slot_mutex);
#endif
}

/* Block until every bulk-read task submitted for `slot` has finished. The
 * task results (col_num_values, last_error, the value buffers) are published
 * by the decrement under slot_mutex, so reading them after this returns is
 * race-free. */
static void pipeline_wait_slot(carquet_batch_reader_t* br, rg_slot_t* slot) {
    if (!br->slot_sync_init) {
        carquet_worker_pool_wait(br->pool);
        return;
    }
    slot_sync_lock(br);
    while (slot->pending_tasks > 0) {
#ifdef _WIN32
        SleepConditionVariableCS(&br->slot_cond, &br->slot_mutex, INFINITE);
#else
        pthread_cond_wait(&br->slot_cond, &br->slot_mutex);
#endif
    }
    slot_sync_unlock(br);
}

static void bulk_read_task(void* arg) {
    bulk_read_arg_t* t = (bulk_read_arg_t*)arg;
    bulk_read_task_body(t);
    if (t->owner && t->slot_pending && t->owner->slot_sync_init) {
        slot_sync_lock(t->owner);
        if (--(*t->slot_pending) == 0) {
#ifdef _WIN32
            WakeAllConditionVariable(&t->owner->slot_cond);
#else
            pthread_cond_broadcast(&t->owner->slot_cond);
#endif
        }
        slot_sync_unlock(t->owner);
    }
}

/* Copy a segment's validity bits (`bits`, bit 0 = row t->row_start, n
 * rows, readable one byte past its last) into the slot bitmap. The bytes the
 * segment covers alone are stored; its first byte when it starts mid-byte
 * and its last when it ends mid-byte are shared with the neighbouring
 * segments, so those go to head/tail for the consumer to merge
 * (merge_segment_bitmaps). */
static void write_segment_validity(bulk_read_arg_t* t, const uint8_t* bits, int64_t n) {
    uint8_t* bm = t->bitmap;
    int64_t start = t->row_start;
    int64_t end = start + n;
    int64_t first = start >> 3;
    int64_t last = (end - 1) >> 3;
    unsigned shift = (unsigned)(start & 7);
    for (int64_t ob = first; ob <= last; ob++) {
        /* Bits of absolute byte ob: segment bits [ob*8 - start, +8). */
        int64_t off = ob * 8 - start;
        unsigned v;
        if (off < 0) {
            v = (unsigned)bits[0] << shift;
        } else if (shift == 0) {
            v = bits[off >> 3];
        } else {
            int64_t q = off >> 3;          /* off = 8q + (8 - shift) */
            v = ((unsigned)bits[q] >> (8 - shift)) | ((unsigned)bits[q + 1] << shift);
        }
        if (ob == last && (end & 7)) v &= (1u << (end & 7)) - 1;
        if (ob == first && shift) {
            t->head_byte = ob;
            t->head_bits = (uint8_t)v;
        } else if (ob == last && (end & 7)) {
            t->tail_byte = ob;
            t->tail_bits = (uint8_t)v;
        } else {
            bm[ob] = (uint8_t)v;
        }
    }
}

/* Validity bits of `n` levels into `bits` at bit offset `at` (LSB first).
 * The bytes from at/8 on must start out zero; whole bytes are built with the
 * dispatched kernel, a partial first or last byte bit by bit. */
static void levels_to_bits_at(const int16_t* levels, int64_t n, int16_t max_def,
                              uint8_t* bits, int64_t at) {
    int64_t i = 0;
    for (; i < n && ((at + i) & 7); i++) {
        if (levels[i] == max_def) bits[(at + i) >> 3] |= (uint8_t)(1u << ((at + i) & 7));
    }
    int64_t whole = ((n - i) / 8) * 8;
    if (whole > 0) {
        carquet_dispatch_build_null_bitmap(levels + i, whole, max_def, bits + ((at + i) >> 3));
        i += whole;
    }
    for (; i < n; i++) {
        if (levels[i] == max_def) bits[(at + i) >> 3] |= (uint8_t)(1u << ((at + i) & 7));
    }
}

/* Scatter: row i of `dst` takes the next dense value of `src` when levels[i]
 * == max_def, else zero; `ndense` (> 0) values are available. Branch-free:
 * the source index is clamped to the last value and the value masked. */
#define SCATTER_BODY(T, LANES)                                                   \
    do {                                                                         \
        int64_t k = 0, last = ndense - 1;                                        \
        for (int64_t i = 0; i < n; i++) {                                        \
            int64_t b = levels[i] == max_def;                                    \
            int64_t idx = k < last ? k : last;                                   \
            T v[LANES], m = (T)0 - (T)b;                                         \
            memcpy(v, src + (size_t)idx * sizeof(v), sizeof(v));                 \
            for (int l = 0; l < LANES; l++) v[l] &= m;                           \
            memcpy(dst + (size_t)i * sizeof(v), v, sizeof(v));                   \
            k += b;                                                              \
        }                                                                        \
    } while (0)

static void scatter_dense_to_rows(uint8_t* dst, const uint8_t* src, int64_t ndense,
                                  const int16_t* levels, int64_t n, int16_t max_def,
                                  size_t value_size) {
    if (ndense <= 0) {
        memset(dst, 0, (size_t)n * value_size);
        return;
    }
    if (ndense == n) {
        memcpy(dst, src, (size_t)n * value_size);
        return;
    }
    switch (value_size) {
    case 4: SCATTER_BODY(uint32_t, 1); return;
    case 8: SCATTER_BODY(uint64_t, 1); return;
    case 16: SCATTER_BODY(uint64_t, 2); return;
    default: break;
    }
    int64_t k = 0;
    for (int64_t i = 0; i < n; i++) {
        uint8_t* d = dst + (size_t)i * value_size;
        if (levels[i] == max_def) {
            memcpy(d, src + (size_t)k * value_size, value_size);
            k++;
        } else {
            memset(d, 0, value_size);
        }
    }
}
#undef SCATTER_BODY

/* SCATTER_OFF[b][j]: for row j of a bitmap byte b, the dense offset of its
 * value from the byte's first value when bit j is set (the set bits below
 * j), else 0. A constant table, built by the preprocessor. */
#define SC_POP_BELOW(b, j) \
    ((((b) & ((1u << (j)) - 1)) >> 0 & 1) + (((b) & ((1u << (j)) - 1)) >> 1 & 1) + \
     (((b) & ((1u << (j)) - 1)) >> 2 & 1) + (((b) & ((1u << (j)) - 1)) >> 3 & 1) + \
     (((b) & ((1u << (j)) - 1)) >> 4 & 1) + (((b) & ((1u << (j)) - 1)) >> 5 & 1) + \
     (((b) & ((1u << (j)) - 1)) >> 6 & 1))
#define SC_OFF(b, j) (((b) >> (j) & 1) ? SC_POP_BELOW(b, j) : 0)
#define SC_ROW(b) { SC_OFF(b, 0), SC_OFF(b, 1), SC_OFF(b, 2), SC_OFF(b, 3), \
                    SC_OFF(b, 4), SC_OFF(b, 5), SC_OFF(b, 6), SC_OFF(b, 7) }
#define SC_R4(b) SC_ROW(b), SC_ROW((b) + 1), SC_ROW((b) + 2), SC_ROW((b) + 3)
#define SC_R16(b) SC_R4(b), SC_R4((b) + 4), SC_R4((b) + 8), SC_R4((b) + 12)
#define SC_R64(b) SC_R16(b), SC_R16((b) + 16), SC_R16((b) + 32), SC_R16((b) + 48)
static const uint8_t SCATTER_OFF[256][8] = {
    SC_R64(0u), SC_R64(64u), SC_R64(128u), SC_R64(192u)
};
#undef SC_R64
#undef SC_R16
#undef SC_R4
#undef SC_ROW
#undef SC_OFF
#undef SC_POP_BELOW

/* scatter_dense_to_rows driven by a validity bitmap: row i of `dst` takes
 * the next dense value of `src` when bit (at + i) of `bm` is set, else zero.
 * 64 rows at a time: a full word is one block copy, an empty one a memset.
 * A mixed word goes 8 rows (one bitmap byte) at a time with no dependency
 * between the rows: row j loads value k + SCATTER_OFF[byte][j] (k = values
 * placed before the byte) and masks it off when absent. An absent row loads
 * value k, which exists because a nonzero byte has a present row at or past
 * k; an all-zero byte is a memset. */
#define SCATTER_BITS_BODY(T, LANES)                                              \
    do {                                                                         \
        const size_t vs = sizeof(T) * LANES;                                     \
        int64_t k = 0;                                                           \
        for (int64_t i = 0; i < n; i += 64) {                                    \
            int h = n - i < 64 ? (int)(n - i) : 64;                              \
            uint64_t w = carquet_bitmap_load_bits(bm, at + i, h);                \
            uint8_t* d = dst + (size_t)i * vs;                                   \
            if (h == 64 && w == ~(uint64_t)0) {                                  \
                memcpy(d, src + (size_t)k * vs, 64 * vs);                        \
                k += 64;                                                         \
                continue;                                                        \
            }                                                                    \
            if (w == 0) {                                                        \
                memset(d, 0, (size_t)h * vs);                                    \
                continue;                                                        \
            }                                                                    \
            for (int g = 0; g < h; g += 8, w >>= 8, d += 8 * vs) {               \
                unsigned bb = (unsigned)(w & 0xFF);                              \
                int rows = h - g < 8 ? h - g : 8;                                \
                if (bb == 0) {                                                   \
                    memset(d, 0, (size_t)rows * vs);                             \
                    continue;                                                    \
                }                                                                \
                const uint8_t* off = SCATTER_OFF[bb];                            \
                const uint8_t* s0 = src + (size_t)k * vs;                        \
                if (rows == 8) {                                                 \
                    for (int j = 0; j < 8; j++) {                                \
                        T v[LANES], m = (T)0 - (T)((bb >> j) & 1);               \
                        memcpy(v, s0 + (size_t)off[j] * vs, vs);                 \
                        for (int l = 0; l < LANES; l++) v[l] &= m;               \
                        memcpy(d + (size_t)j * vs, v, vs);                       \
                    }                                                            \
                } else {                                                         \
                    for (int j = 0; j < rows; j++) {                             \
                        T v[LANES], m = (T)0 - (T)((bb >> j) & 1);               \
                        memcpy(v, s0 + (size_t)off[j] * vs, vs);                 \
                        for (int l = 0; l < LANES; l++) v[l] &= m;               \
                        memcpy(d + (size_t)j * vs, v, vs);                       \
                    }                                                            \
                }                                                                \
                k += carquet_popcount32(bb);                                     \
            }                                                                    \
        }                                                                        \
    } while (0)

static void scatter_dense_by_bits(uint8_t* dst, const uint8_t* src, int64_t ndense,
                                  const uint8_t* bm, int64_t at, int64_t n,
                                  size_t value_size) {
    if (ndense <= 0) {
        memset(dst, 0, (size_t)n * value_size);
        return;
    }
    if (ndense == n) {
        memcpy(dst, src, (size_t)n * value_size);
        return;
    }
    switch (value_size) {
    case 4: SCATTER_BITS_BODY(uint32_t, 1); return;
    case 8: SCATTER_BITS_BODY(uint64_t, 1); return;
    case 16: SCATTER_BITS_BODY(uint64_t, 2); return;
    default: break;
    }
    int64_t k = 0;
    for (int64_t i = 0; i < n; i++) {
        uint8_t* d = dst + (size_t)i * value_size;
        if ((bm[(at + i) >> 3] >> ((at + i) & 7)) & 1) {
            memcpy(d, src + (size_t)k * value_size, value_size);
            k++;
        } else {
            memset(d, 0, value_size);
        }
    }
}
#undef SCATTER_BITS_BODY

/* Read `want` rows of a flat column with definition levels straight out of
 * its decoded pages: values land at their rows (nulls zeroed) and validity
 * bits at bits[0, want) (zeroed here). This is carquet_read_next_page's
 * page walk and state update without its two copies (dense values, levels)
 * and the in-place spread that followed. Returns the rows read; a page
 * failure leaves the detail in cr->last_error. */
static int64_t read_rows_from_pages(carquet_column_reader_t* cr, uint8_t* dest,
                                    int64_t want, size_t value_size, int16_t max_def,
                                    uint8_t* bits) {
    memset(bits, 0, (size_t)((want + 7) / 8));
    cr->last_error.code = CARQUET_OK;
    cr->last_error.message[0] = '\0';
    int64_t done = 0;
    while (done < want && cr->values_remaining > 0) {
        /* A local error, as in carquet_column_read_batch_ex: a page load can
         * succeed after recording a recoverable failure (e.g. a page header
         * re-read with a larger window without mmap), which must not read
         * as an error. */
        carquet_error_t load_error = CARQUET_ERROR_INIT;
        if (carquet_column_ensure_page_loaded(cr, &load_error) != CARQUET_OK) {
            cr->last_error = load_error;
            if (cr->last_error.code == CARQUET_OK) {
                CARQUET_SET_ERROR(&cr->last_error, CARQUET_ERROR_DECODE,
                                  "Page load failed with no detail reported");
            }
            break;
        }
        int64_t take = (int64_t)cr->page_num_values - cr->page_values_read;
        if (take > want - done) take = want - done;
        if (take <= 0) break;
        uint8_t* dst = dest + (size_t)done * value_size;
        if (carquet_column_def_as_bits(cr) && cr->decoded_def_bits) {
            /* Levels decoded as a validity bitmap: counts are popcounts, the
             * spread is driven by the bits and they are copied as they are. */
            const uint8_t* pb = cr->decoded_def_bits;
            int64_t pvr = cr->page_values_read;
            int32_t dense_start;
            if (pvr == 0) {
                dense_start = 0;
            } else if (cr->dense_cursor_pos == cr->page_values_read) {
                dense_start = cr->dense_cursor_val;
            } else {
                dense_start = (int32_t)carquet_bitmap_count_range(pb, 0, pvr);
            }
            int64_t present = carquet_bitmap_count_range(pb, pvr, take);
            scatter_dense_by_bits(dst, (const uint8_t*)cr->decoded_values +
                                       (size_t)dense_start * value_size,
                                  present, pb, pvr, take, value_size);
            carquet_bitmap_copy_range(bits, done, pb, pvr, take);
            cr->dense_cursor_pos = cr->page_values_read + (int32_t)take;
            cr->dense_cursor_val = dense_start + (int32_t)present;
        } else if (!cr->decoded_def_levels) {
            /* No levels decoded for this page: every row reads as null, as
             * carquet_read_next_page reports it (zero levels). */
            memset(dst, 0, (size_t)take * value_size);
        } else {
            const int16_t* lv = cr->decoded_def_levels + cr->page_values_read;
            int32_t dense_start;
            if (cr->page_values_read == 0) {
                dense_start = 0;
            } else if (cr->dense_cursor_pos == cr->page_values_read) {
                dense_start = cr->dense_cursor_val;
            } else {
                dense_start = (int32_t)carquet_dispatch_count_non_nulls(
                    cr->decoded_def_levels, cr->page_values_read, max_def);
            }
            int64_t present = carquet_dispatch_count_non_nulls(lv, take, max_def);
            scatter_dense_to_rows(dst, (const uint8_t*)cr->decoded_values +
                                       (size_t)dense_start * value_size,
                                  present, lv, take, max_def, value_size);
            levels_to_bits_at(lv, take, max_def, bits, done);
            cr->dense_cursor_pos = cr->page_values_read + (int32_t)take;
            cr->dense_cursor_val = dense_start + (int32_t)present;
        }
        cr->page_values_read += (int32_t)take;
        cr->values_remaining -= take;
        done += take;
    }
    return done;
}

/* Rows per read-and-spread step of row_aligned_read (a multiple of 64). */
#define ROW_ALIGNED_STEP ((int64_t)16384)

/* Row-aligned read of t->max_values rows into t->dest (and, with levels,
 * t->bitmap) through the page-by-page column reader: the slot's own reader
 * for a whole chunk, or the task's seg_reader positioned at the segment's
 * first page. Any encoding the page reader supports works here. A failure
 * leaves *out_values_read short; the consumer then re-reads the chunk. */
static void row_aligned_read(bulk_read_arg_t* t) {
    *t->out_values_read = 0;
    t->head_byte = -1;
    t->tail_byte = -1;
    carquet_column_reader_t* cr = t->use_seg_reader ? t->seg_reader : t->col_reader;
    if (!cr) return;
    if (t->use_seg_reader &&
        carquet_column_reader_seek_to_data_page(cr, t->start_offset, 0,
                                                &cr->last_error) != CARQUET_OK) {
        return;
    }

    if (t->max_def <= 0) {
        int64_t n = carquet_column_read_batch_ex(cr, t->dest, t->max_values,
                                                 NULL, NULL, &cr->last_error);
        if (n != t->max_values || cr->last_error.code != CARQUET_OK) {
            *t->out_values_read = n > 0 && n < t->max_values ? n : 0;
            return;
        }
        *t->out_values_read = n;
        return;
    }

    /* With levels, read ROW_ALIGNED_STEP rows at a time and spread each
     * step right away, while its dense values are still in cache: spreading
     * a whole chunk afterwards is one more pass over memory. The step is a
     * multiple of 64, so every step's validity bits start on a byte (and
     * bitmap word) boundary. A whole chunk writes its bits straight into the
     * slot bitmap; a segment into its own vbits, copied (shifted) at the end. */
    const int64_t step = ROW_ALIGNED_STEP;
    int64_t level_count = t->max_values < step ? t->max_values : step;
    if (t->levels_cap < (size_t)level_count) {
        int16_t* grown = carquet_mem_realloc(t->levels,
                                             (size_t)level_count * sizeof(int16_t));
        if (!grown) {
            CARQUET_SET_ERROR(&cr->last_error, CARQUET_ERROR_OUT_OF_MEMORY,
                              "Cannot allocate %lld definition levels",
                              (long long)level_count);
            return;
        }
        t->levels = grown;
        t->levels_cap = (size_t)level_count;
    }
    uint8_t* bits = t->bitmap;
    size_t nbytes = (size_t)((t->max_values + 7) / 8) + 1;   /* +1: read past the last */
    if (t->is_split) {
        if (t->vbits_cap < nbytes) {
            uint8_t* grown = carquet_mem_realloc(t->vbits, nbytes);
            if (!grown) {
                CARQUET_SET_ERROR(&cr->last_error, CARQUET_ERROR_OUT_OF_MEMORY,
                                  "Cannot allocate a %llu-byte validity bitmap",
                                  (unsigned long long)nbytes);
                return;
            }
            t->vbits = grown;
            t->vbits_cap = nbytes;
        }
        bits = t->vbits;
        bits[nbytes - 1] = 0;
    }

    /* Flat columns read straight from the decoded pages (see
     * read_rows_from_pages); a reader that exposes dictionary indices
     * keeps the generic path. */
    bool direct = cr->max_rep_level == 0 && !cr->preserve_dictionary &&
                  cr->max_def_level == t->max_def;
    int64_t done = 0;
    while (done < t->max_values) {
        int64_t want = t->max_values - done < step ? t->max_values - done : step;
        uint8_t* dest = (uint8_t*)t->dest + (size_t)done * t->value_size;
        int64_t got;
        if (direct) {
            got = read_rows_from_pages(cr, dest, want, t->value_size, t->max_def,
                                       bits + (done >> 3));
        } else {
            got = carquet_column_read_batch_ex(cr, dest, want, t->levels, NULL,
                                               &cr->last_error);
        }
        if (got != want || cr->last_error.code != CARQUET_OK) {
            if (got > 0 && got < want) done += got;
            *t->out_values_read = done < t->max_values ? done : 0;
            return;
        }
        if (!direct) {
            levels_to_rows(dest, t->levels, got, t->max_def, t->value_size, bits + (done >> 3));
        }
        done += got;
    }
    if (t->is_split) write_segment_validity(t, bits, done);
    *t->out_values_read = done;
}

/* List `a` followed by list `b` (either may be NULL). */
static carquet_retained_page_t* retained_pages_concat(carquet_retained_page_t* a,
                                                      carquet_retained_page_t* b) {
    if (!a) return b;
    carquet_retained_page_t* p = a;
    while (p->next) p = p->next;
    p->next = b;
    return a;
}

static void bulk_read_task_body(bulk_read_arg_t* t) {
    if (!t->col_reader || !t->dest || t->max_values <= 0) {
        *t->out_values_read = 0;
        return;
    }

    /* ------------------------------------------------------------------
     * Filtered branch: read only matching pages for this column, writing
     * each range's rows contiguously into the slot buffer.
     *
     * When this column has an offset index we seek by file offset to the
     * page covering each range. When it doesn't (e.g. an externally
     * written file that supplied a page index for the predicate column
     * but not for this one), we degrade to monotonic read-and-discard:
     * skip the gap between the previous range end and the next range
     * start, then read the range's row count. This is the same fallback
     * the sequential filtered path uses (§6.5 of the design doc).
     * ------------------------------------------------------------------ */
    if (t->ranges && t->ranges->count > 0 && t->value_size > 0) {
        carquet_error_t err = CARQUET_ERROR_INIT;
        /* A seek flushes the reader's retained pages, which the BYTE_ARRAY
         * values of the ranges already read point into: set them aside
         * across each seek and hand them all back at the end. */
        carquet_retained_page_t* kept = NULL;
        size_t dest_offset_bytes = 0;
        int64_t total_read = 0;
        int64_t cursor_row = 0;   /* logical row position used by the
                                   * no-offset-index fallback. */
        for (int32_t r = 0; r < t->ranges->count; r++) {
            int64_t first = t->ranges->ranges[r].first_row;
            int64_t num = t->ranges->ranges[r].num_rows;
            if (num <= 0) continue;

            if (t->offset_index) {
                int64_t page_first_row = 0;
                int32_t page_idx = find_page_for_row(
                    t->offset_index, t->rg_num_rows, first, &page_first_row);
                if (page_idx < 0) break;

                carquet_page_location_t loc;
                if (carquet_offset_index_get_page_location(
                        t->offset_index, page_idx, &loc) != CARQUET_OK) {
                    break;
                }

                kept = retained_pages_concat(kept, t->col_reader->retained_pages);
                t->col_reader->retained_pages = NULL;
                if (carquet_column_reader_seek_to_data_page(
                        t->col_reader, loc.offset, 0, &err) != CARQUET_OK) {
                    break;
                }

                int64_t intra_skip = first - page_first_row;
                if (intra_skip > 0) {
                    int64_t skipped = carquet_column_skip(
                        t->col_reader, intra_skip);
                    if (skipped != intra_skip) break;
                }
            } else {
                /* Forward read-and-discard from cursor to range start. */
                int64_t gap = first - cursor_row;
                if (gap > 0) {
                    int64_t skipped = carquet_column_skip(t->col_reader, gap);
                    if (skipped != gap) break;
                }
            }

            uint8_t* dest_ptr = (uint8_t*)t->dest + dest_offset_bytes;
            int64_t got = carquet_column_read_batch_ex(
                t->col_reader, dest_ptr, num, NULL, NULL,
                &t->col_reader->last_error);
            if (got != num) break;
            dest_offset_bytes += (size_t)num * t->value_size;
            total_read += num;
            cursor_row = first + num;
        }
        t->col_reader->retained_pages =
            retained_pages_concat(t->col_reader->retained_pages, kept);
        *t->out_values_read = total_read;
        return;
    }

    /* ------------------------------------------------------------------
     * Unfiltered branch: existing fast paths.
     * ------------------------------------------------------------------ */
    if (t->max_def > 0 || t->use_seg_reader) {
        row_aligned_read(t);
    } else if (t->is_split) {
        /* One page-range segment of a chunk whose other segments are being
         * decoded concurrently by other tasks: the column reader is shared,
         * so no inline fallback through it. A bail leaves the segment's count
         * short; the consumer notices the column total is short after the
         * slot's wait and re-reads the chunk on its own thread. */
        (void)coalesced_read_segment(t->col_reader, t->data_base,
                                     t->dest, t->max_values,
                                     t->start_offset, t->end_offset,
                                     t->owner, t->out_values_read);
    } else if (can_coalesce_column(t->col_reader)) {
        coalesced_read_column(t->col_reader, t->dest, t->max_values,
                              t->out_values_read);
    } else {
        *t->out_values_read = carquet_column_read_batch_ex(
            t->col_reader, t->dest, t->max_values, NULL, NULL,
            &t->col_reader->last_error);
    }
}

/* ============================================================================
 * Coalesced Column Chunk Read — Fast Path for Pipeline Mode
 * ============================================================================
 * Instead of the page-by-page carquet_column_read_batch path, this scans all
 * page headers upfront and decompresses in a tight loop, writing directly to
 * the output buffer. Eliminates intermediate buffers and per-page state machine
 * overhead.
 *
 * Eligible: REQUIRED columns, fixed-width, PLAIN or BYTE_STREAM_SPLIT encoding,
 * no dictionary, compressed, mmap available.
 */

static bool can_coalesce_column(const carquet_column_reader_t* cr) {
    if (!cr || !cr->col_meta || !cr->file_reader || !cr->file_reader->mmap_data) return false;
    if (cr->max_def_level > 0 || cr->max_rep_level > 0) return false;
    if (cr->type == CARQUET_PHYSICAL_BOOLEAN || cr->type == CARQUET_PHYSICAL_BYTE_ARRAY) return false;
    if (cr->col_meta->has_dictionary_page_offset) return false;
    if (cr->col_meta->codec == CARQUET_COMPRESSION_UNCOMPRESSED) return false;
    return true;
}

/* Scratch pool for segment tasks (BYTE_STREAM_SPLIT staging). Entries
 * [0, scratch_free) are free; an entry's buffer may still be NULL (capacity 0)
 * until a task first grows it. Returns false when no entry is free, which
 * cannot happen while every task is run by a pool thread (one entry per
 * thread) but is handled anyway. */
static bool scratch_acquire(carquet_batch_reader_t* br, uint8_t** buf, size_t* cap) {
    bool ok = false;
    slot_sync_lock(br);
    if (br->scratch_free > 0) {
        br->scratch_free--;
        *buf = br->scratch_bufs[br->scratch_free];
        *cap = br->scratch_sizes[br->scratch_free];
        br->scratch_bufs[br->scratch_free] = NULL;
        br->scratch_sizes[br->scratch_free] = 0;
        ok = true;
    }
    slot_sync_unlock(br);
    return ok;
}

static void scratch_release(carquet_batch_reader_t* br, uint8_t* buf, size_t cap) {
    slot_sync_lock(br);
    if (br->scratch_free < br->scratch_cap) {
        br->scratch_bufs[br->scratch_free] = buf;
        br->scratch_sizes[br->scratch_free] = cap;
        br->scratch_free++;
        buf = NULL;
    }
    slot_sync_unlock(br);
    carquet_mem_free(buf);   /* pool full (cannot happen): don't leak */
}

/**
 * Decode every data page in [start_offset, end_offset) of a coalescable
 * column chunk straight into dest. Returns true when the whole range was
 * decoded; on false *out_values_read holds the number of values decoded
 * before the bail (whole pages only) and the caller recovers.
 *
 * scratch_owner != NULL selects the batch reader's scratch pool for
 * BYTE_STREAM_SPLIT staging instead of the column reader's decompress_buffer,
 * which is required when several segments of one chunk decode concurrently.
 */
static bool coalesced_read_segment(
    const carquet_column_reader_t* cr,
    const uint8_t* data_base,
    void* dest,
    int64_t max_values,
    int64_t start_offset,
    int64_t end_offset,
    carquet_batch_reader_t* scratch_owner,
    int64_t* out_values_read) {

    const carquet_reader_t* fr = cr->file_reader;
    carquet_column_reader_t* reader = (carquet_column_reader_t*)cr;
    const parquet_column_metadata_t* meta = cr->col_meta;
    size_t file_size = fr->file_size;

    uint8_t* out = (uint8_t*)dest;
    int64_t total_values = 0;

    /* BYTE_STREAM_SPLIT staging buffer: the column reader's own (single
     * task per chunk) or one borrowed from the scratch pool (segments). */
    uint8_t* stage = NULL;
    size_t stage_cap = 0;
    bool stage_from_pool = false;

    size_t value_size = 0;
    switch (cr->type) {
        case CARQUET_PHYSICAL_INT32: case CARQUET_PHYSICAL_FLOAT: value_size = 4; break;
        case CARQUET_PHYSICAL_INT64: case CARQUET_PHYSICAL_DOUBLE: value_size = 8; break;
        case CARQUET_PHYSICAL_INT96: value_size = 12; break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: value_size = (size_t)cr->type_length; break;
        default: goto bail;
    }
    if (value_size == 0) goto bail;

    int64_t offset = start_offset;
    int64_t chunk_end = end_offset;
    if (offset <= 0) offset = meta->data_page_offset;
    if (chunk_end <= 0) chunk_end = meta->data_page_offset + meta->total_compressed_size;
    if (offset < meta->data_page_offset || chunk_end > meta->data_page_offset + meta->total_compressed_size) {
        goto bail;
    }
    if (offset < 0 || chunk_end > (int64_t)file_size) goto bail;

    while (offset < chunk_end && total_values < max_values) {
        /* Parse page header */
        const uint8_t* ptr = data_base + offset;
        size_t remaining = (size_t)(chunk_end - offset);
        size_t max_hdr = remaining < CARQUET_PAGE_HEADER_SCAN_WINDOW
                       ? remaining : CARQUET_PAGE_HEADER_SCAN_WINDOW;

        parquet_page_header_t hdr;
        size_t hdr_size;
        if (parquet_parse_page_header(ptr, max_hdr, &hdr, &hdr_size, NULL) != CARQUET_OK)
            goto bail;

        if (hdr.type != CARQUET_PAGE_DATA)
            goto bail;

        int32_t num_values = hdr.data_page_header.num_values;
        carquet_encoding_t encoding = hdr.data_page_header.encoding;

        if (num_values <= 0 || total_values + num_values > max_values) goto bail;

        const uint8_t* compressed = ptr + hdr_size;
        size_t comp_size = (size_t)hdr.compressed_page_size;
        size_t uncomp_size = (size_t)hdr.uncompressed_page_size;
        if (hdr.compressed_page_size <= 0 ||
            hdr.uncompressed_page_size < 0 ||
            comp_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE ||
            uncomp_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE ||
            hdr_size > remaining ||
            comp_size > remaining - hdr_size) {
            goto bail;
        }

        /* This path decompresses straight into the caller's buffer, which has
         * room for exactly num_values * value_size bytes for this page. Only
         * proceed when the header's uncompressed size agrees; otherwise the
         * codec would be handed a capacity larger than the space that actually
         * exists. */
        if (uncomp_size != (size_t)num_values * value_size) {
            goto bail;
        }

        /* Same integrity check as the page-by-page reader, so this path does
         * not trade verification for speed. */
        if (hdr.has_crc && fr->options.verify_checksums &&
            carquet_crc32(compressed, comp_size) != (uint32_t)hdr.crc) {
            goto bail;
        }

        if (encoding == CARQUET_ENCODING_PLAIN) {
            /* Decompress directly to output — data IS the final values */
            size_t actual;
            if (carquet_decompress_page(meta->codec, compressed, comp_size,
                                         out, uncomp_size, &actual) != CARQUET_OK)
                goto bail;
            /* A short payload would leave the tail of this page's values
             * holding whatever the buffer held before. */
            if (actual != uncomp_size) goto bail;
        } else if (encoding == CARQUET_ENCODING_BYTE_STREAM_SPLIT) {
            /* Decompress to staging, then cache-tiled transpose to output */
            if (!stage && !stage_from_pool) {
                if (scratch_owner) {
                    if (!scratch_acquire(scratch_owner, &stage, &stage_cap)) goto bail;
                    stage_from_pool = true;
                } else {
                    stage = reader->decompress_buffer;
                    stage_cap = reader->decompress_capacity;
                }
            }
            if (uncomp_size > stage_cap) {
                uint8_t* new_buf = carquet_mem_realloc(stage, uncomp_size);
                if (!new_buf) {
                    goto bail;
                }
                stage = new_buf;
                stage_cap = uncomp_size;
                if (!stage_from_pool) {
                    reader->decompress_buffer = stage;
                    reader->decompress_capacity = stage_cap;
                }
            }
            size_t actual;
            if (carquet_decompress_page(meta->codec, compressed, comp_size,
                                         stage, uncomp_size, &actual) != CARQUET_OK)
                goto bail;
            if (actual != uncomp_size) goto bail;
            if (value_size == 4) {
                if (carquet_byte_stream_split_decode_float(
                        stage, actual, (float*)out, num_values) != CARQUET_OK) {
                    goto bail;
                }
            } else if (value_size == 8) {
                if (carquet_byte_stream_split_decode_double(
                        stage, actual, (double*)out, num_values) != CARQUET_OK) {
                    goto bail;
                }
            } else {
                /* Other widths (FIXED_LEN_BYTE_ARRAY) have value_size planes,
                 * which the 4- and 8-plane transposers above would misread. */
                if (carquet_byte_stream_split_decode(
                        stage, actual, (int32_t)value_size,
                        out, num_values) != CARQUET_OK) {
                    goto bail;
                }
            }
        } else {
            /* Unsupported encoding — bail */
            goto bail;
        }

        out += (size_t)num_values * value_size;
        total_values += num_values;
        offset += (int64_t)hdr_size + (int64_t)comp_size;
    }

    if (stage_from_pool) scratch_release(scratch_owner, stage, stage_cap);
    *out_values_read = total_values;
    return true;

bail:
    if (stage_from_pool) scratch_release(scratch_owner, stage, stage_cap);
    *out_values_read = total_values;
    return false;
}

static void coalesced_read_column_range(
    const carquet_column_reader_t* cr,
    const uint8_t* data_base,
    void* dest,
    int64_t max_values,
    int64_t start_offset,
    int64_t end_offset,
    int64_t* out_values_read) {

    if (coalesced_read_segment(cr, data_base, dest, max_values,
                               start_offset, end_offset, NULL,
                               out_values_read)) {
        return;
    }

    /* Fall back to standard page-by-page reader. It re-reads the chunk from
     * its start and applies every check itself, so whatever made the fast
     * path bail (including a CRC mismatch) is reported with its own detail.
     * Only valid for a whole-chunk read (segments never reach this). */
    carquet_column_reader_t* reader = (carquet_column_reader_t*)cr;
    *out_values_read = carquet_column_read_batch_ex(
        reader, dest, max_values, NULL, NULL, &reader->last_error);
}

static void coalesced_read_column(
    const carquet_column_reader_t* cr,
    void* dest,
    int64_t max_values,
    int64_t* out_values_read) {
    coalesced_read_column_range(cr, cr->file_reader->mmap_data,
                                dest, max_values, 0, 0, out_values_read);
}

/**
 * Plan N-way split of a column chunk for parallel decompression.
 * Walks page headers to find page-aligned boundaries that divide the chunk
 * into roughly equal pieces (by value count). Only V1 data pages are split,
 * and only flat columns (max_rep == 0), whose page value count is also its
 * row count.
 *
 * @param cr          Column reader of a flat column of a mapped file
 * @param max_values  Total values in the column chunk
 * @param max_splits  Maximum number of splits to produce (>= 2)
 * @param split_offsets  Output: [max_splits+1] byte offsets (start of each segment + end)
 * @param split_values   Output: [max_splits+1] cumulative value counts at each boundary
 * @return Number of segments (>= 2 on success, 0 if split not possible)
 */
static int32_t plan_page_splits(
    const carquet_column_reader_t* cr,
    const uint8_t* data_base,
    int64_t max_values,
    int32_t max_splits,
    int64_t* split_offsets,
    int64_t* split_values) {

    /* Every built-in codec (and any registered custom one) is already
     * entered concurrently by the per-column tasks, so segments add no new
     * thread-safety requirement. */
    if (!split_offsets || !split_values || !data_base || !cr || !cr->col_meta ||
        cr->max_rep_level > 0) return 0;
    if (max_splits < 2) return 0;
    const int64_t min_seg_values = 65536;
    if (max_values < 2 * min_seg_values) return 0;

    const carquet_reader_t* fr = cr->file_reader;
    const parquet_column_metadata_t* meta = cr->col_meta;
    int64_t offset = meta->data_page_offset;
    /* total_compressed_size spans the whole chunk, dictionary page included,
     * and that page comes first: measure the end from the chunk's start. A
     * dictionary offset at or past the data pages (legacy writers) is not
     * trusted as the start; such a chunk ends up unsplit below. */
    int64_t chunk_start = offset;
    if (meta->has_dictionary_page_offset && meta->dictionary_page_offset >= 0 &&
        meta->dictionary_page_offset < offset) {
        chunk_start = meta->dictionary_page_offset;
    }
    if (chunk_start < 0 || meta->total_compressed_size < 0 ||
        chunk_start > (int64_t)fr->file_size ||
        meta->total_compressed_size > (int64_t)fr->file_size - chunk_start) return 0;
    int64_t chunk_end = chunk_start + meta->total_compressed_size;
    if (offset >= chunk_end) return 0;

    /* First segment starts at the first data page */
    split_offsets[0] = offset;
    split_values[0] = 0;
    int32_t num_segments = 1;

    int64_t target_per_split = max_values / max_splits;
    /* Segments below ~64k values cost more in task dispatch than they gain */
    if (target_per_split < min_seg_values) {
        max_splits = (int32_t)(max_values / min_seg_values);
        if (max_splits < 2) return 0;
        target_per_split = max_values / max_splits;
    }
    int64_t next_target = target_per_split;
    int64_t values_so_far = 0;

    while (offset < chunk_end) {
        const uint8_t* ptr = data_base + offset;
        size_t remaining = (size_t)(chunk_end - offset);
        size_t max_hdr = remaining < CARQUET_PAGE_HEADER_SCAN_WINDOW
                       ? remaining : CARQUET_PAGE_HEADER_SCAN_WINDOW;

        parquet_page_header_t hdr;
        size_t hdr_size;
        if (parquet_parse_page_header(ptr, max_hdr, &hdr, &hdr_size, NULL) != CARQUET_OK)
            return 0;

        if (hdr.type != CARQUET_PAGE_DATA)
            return 0;

        if (hdr.compressed_page_size <= 0 ||
            hdr.uncompressed_page_size < 0 ||
            (size_t)hdr.compressed_page_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE ||
            (size_t)hdr.uncompressed_page_size > CARQUET_MAX_PAGE_PAYLOAD_SIZE ||
            hdr_size > remaining ||
            (size_t)hdr.compressed_page_size > remaining - hdr_size) {
            return 0;
        }

        int32_t num_values = hdr.data_page_header.num_values;
        if (num_values <= 0) return 0;

        values_so_far += num_values;
        /* Page headers are untrusted: each segment's destination is
         * col_values + split_values[k] * value_size inside a buffer of
         * max_values values, so the page counts must never exceed it. */
        if (values_so_far > max_values) return 0;
        offset += (int64_t)hdr_size + (int64_t)hdr.compressed_page_size;

        /* Place a split boundary when we've accumulated enough values,
         * but only if there's still data remaining for the next segment. */
        if (values_so_far >= next_target && offset < chunk_end &&
            num_segments < max_splits) {
            split_offsets[num_segments] = offset;
            split_values[num_segments] = values_so_far;
            num_segments++;
            next_target = values_so_far + target_per_split;
        }
    }

    /* A chunk whose pages do not add up to exactly max_values is not split
     * (the unsplit path checks every page against its own bound). */
    if (values_so_far != max_values) return 0;

    /* Close the final segment */
    split_offsets[num_segments] = chunk_end;
    split_values[num_segments] = values_so_far;

    return (num_segments >= 2) ? num_segments : 0;
}

/* Drop the slot's reference on column i's buffers (col_share[i] != NULL
 * exactly while the slot holds one). Last reference: the slot owns the
 * buffers outright again. Otherwise exported arrays still point into them:
 * they are left to those arrays (the last release frees them) and the slot
 * starts over with none. */
static void slot_reclaim_column(rg_slot_t* slot, int32_t i) {
    br_share_t* sh = slot->col_share[i];
    if (!sh) return;
    slot->col_share[i] = NULL;
#if BR_SHARE_ATOMICS
    if (br_share_dec(&sh->refs) == 0) {
        carquet_mem_free(sh);
        return;
    }
    slot->col_values[i] = NULL;
    slot->col_buf_sizes[i] = 0;
    slot->col_bitmaps[i] = NULL;
    slot->col_bitmap_sizes[i] = 0;
#else
    carquet_mem_free(sh);
#endif
}

/* Put column i's current buffers under a share record holding the slot's
 * reference. Failing to allocate one only means its exports copy. */
static void slot_share_column(rg_slot_t* slot, int32_t i) {
#if BR_SHARE_ATOMICS
    if (slot->col_share[i]) return;
    br_share_t* sh = carquet_mem_malloc(sizeof(*sh));
    if (!sh) return;
    sh->refs = 1;
    sh->values = slot->col_values[i];
    sh->bitmap = slot->col_bitmaps[i];
    slot->col_share[i] = sh;
#else
    (void)slot; (void)i;
#endif
}

static void slot_release_filter_state(rg_slot_t* slot, int32_t num_projected) {
    if (slot->col_offset_indexes) {
        for (int32_t i = 0; i < num_projected; i++) {
            if (slot->col_offset_indexes[i]) {
                carquet_offset_index_free(slot->col_offset_indexes[i]);
                slot->col_offset_indexes[i] = NULL;
            }
        }
    }
    if (slot->filter_ranges_valid) {
        carquet_row_range_list_destroy(&slot->ranges);
        slot->filter_ranges_valid = false;
    }
}

/**
 * Fill pipeline slots by reading entire column chunks in parallel.
 * Each task reads ALL values from one column in one row group.
 *
 * When a page filter is active, each slot's row-range list is computed
 * up front; the worker tasks read only matching pages, sized to the
 * range total rather than the whole row group. Row groups that match no
 * rows are skipped without consuming a pipeline slot (their rows are
 * still credited to rows_skipped).
 */
/* Record the pipeline's first failure; later ones are consequences of it. */
#define PIPELINE_FAIL(br, status_code, ...)                                  \
    do {                                                                     \
        if ((br)->pipeline_error.code == CARQUET_OK) {                       \
            CARQUET_SET_ERROR(&(br)->pipeline_error, (status_code), __VA_ARGS__); \
        }                                                                    \
    } while (0)

/* Most segments worth cutting a dictionary chunk into. Every segment's
 * column reader decodes the dictionary page again, so cap that repeated work
 * at 1/8 of the chunk's data. Returns 0 when the dictionary page header
 * cannot be read (no split). */
static int32_t dictionary_split_limit(const carquet_column_reader_t* cr,
                                      const uint8_t* data_base, int32_t max_splits) {
    const parquet_column_metadata_t* meta = cr->col_meta;
    if (!meta->has_dictionary_page_offset) return max_splits;
    int64_t off = meta->dictionary_page_offset;
    size_t file_size = cr->file_reader->file_size;
    if (off < 0 || (uint64_t)off >= (uint64_t)file_size) return 0;
    size_t remaining = file_size - (size_t)off;
    size_t window = remaining < CARQUET_PAGE_HEADER_SCAN_WINDOW
                  ? remaining : CARQUET_PAGE_HEADER_SCAN_WINDOW;
    parquet_page_header_t hdr;
    size_t hdr_size;
    if (parquet_parse_page_header(data_base + off, window, &hdr, &hdr_size, NULL) != CARQUET_OK ||
        hdr.type != CARQUET_PAGE_DICTIONARY || hdr.uncompressed_page_size < 0) {
        return 0;
    }
    int64_t dict_bytes = hdr.uncompressed_page_size;
    if (dict_bytes == 0) return max_splits;
    int64_t data_bytes = meta->total_uncompressed_size - dict_bytes;
    if (data_bytes <= 0) return 0;
    int64_t limit = 1 + data_bytes / (8 * dict_bytes);
    return limit < max_splits ? (int32_t)limit : max_splits;
}

/* Row-aligned layout of a task (see row_aligned_read); task entries are
 * reused across fills, so every submission sets all of it. */
static void task_set_layout(bulk_read_arg_t* t, int16_t max_def, uint8_t* bitmap,
                            int64_t row_start, bool use_seg_reader) {
    t->max_def = max_def;
    t->bitmap = max_def > 0 ? bitmap : NULL;
    t->row_start = row_start;
    t->head_byte = -1;
    t->tail_byte = -1;
    t->head_bits = 0;
    t->tail_bits = 0;
    t->use_seg_reader = use_seg_reader;
}

/* Give task entries [first, first + count) a column reader on (rg, col),
 * created on first use and rebound afterwards. Runs on the filling thread:
 * creating readers touches the file reader's shared state. */
static bool prepare_segment_readers(carquet_batch_reader_t* br, int32_t first,
                                    int32_t count, int32_t rg, int32_t file_col) {
    for (int32_t k = 0; k < count; k++) {
        bulk_read_arg_t* t = &br->task_args[first + k];
        if (t->seg_reader) {
            carquet_column_reader_rebind(t->seg_reader, br->reader, rg, file_col, false);
        } else {
            carquet_error_t err = CARQUET_ERROR_INIT;
            t->seg_reader = carquet_reader_get_column(br->reader, rg, file_col, &err);
            if (!t->seg_reader) return false;
        }
    }
    return true;
}

/* Register a task on its slot's pending count, then hand it to the pool. */
static void pipeline_submit(carquet_batch_reader_t* br, rg_slot_t* slot,
                            bulk_read_arg_t* t) {
    t->owner = br;
    t->slot_pending = &slot->pending_tasks;
    if (br->slot_sync_init) {
        slot_sync_lock(br);
        slot->pending_tasks++;
        slot_sync_unlock(br);
    }
    carquet_worker_pool_submit(br->pool, bulk_read_task, t);
}

static void pipeline_fill(carquet_batch_reader_t* br) {
    if (!br->pipeline_active || !br->pool) return;
    /* A failed fill is terminal: do not retry the same row group forever. */
    if (br->pipeline_error.code != CARQUET_OK) return;

    bool filter_active =
        br->filter_clauses != NULL && br->filter_clause_count > 0;

    while (br->pipeline_count < br->pipeline_depth &&
           br->rg_order_next < br->rg_order_len) {

        int32_t slot_idx = (br->pipeline_head + br->pipeline_count) % br->pipeline_depth;
        rg_slot_t* slot = &br->pipeline[slot_idx];
        int32_t target_rg = br->rg_order[br->rg_order_next];

        /* Get row count for this RG */
        const parquet_row_group_t* rg = &br->reader->metadata.row_groups[target_rg];
        int64_t rg_rows = rg->num_rows;

        /* Page filter evaluation: skip whole row groups that match nothing,
         * and clip per-slot allocations to the matching row count. */
        slot_release_filter_state(slot, br->num_projected);
        int64_t slot_rows = rg_rows;
        if (filter_active) {
            carquet_error_t feval_err = CARQUET_ERROR_INIT;
            carquet_row_range_list_init(&slot->ranges);
            carquet_status_t fst = carquet_page_filter_eval_row_group(
                br->reader, target_rg,
                br->filter_clauses, br->filter_clause_count,
                &slot->ranges, &feval_err);
            if (fst != CARQUET_OK) {
                /* Skipping the row group would silently drop its rows. */
                carquet_row_range_list_destroy(&slot->ranges);
                if (br->pipeline_error.code == CARQUET_OK) {
                    br->pipeline_error = feval_err;
                    if (br->pipeline_error.code == CARQUET_OK) {
                        PIPELINE_FAIL(br, fst, "Row group %d: page filter evaluation failed",
                                      target_rg);
                    }
                }
                return;
            }
            slot->filter_ranges_valid = true;

            br->rows_skipped += rg_rows - slot->ranges.total_rows;
            if (slot->ranges.count == 0) {
                /* No rows from this RG; do not occupy a pipeline slot. */
                carquet_row_range_list_destroy(&slot->ranges);
                slot->filter_ranges_valid = false;
                br->rg_order_next++;
                continue;
            }
            slot_rows = slot->ranges.total_rows;

            /* Load per-column offset indexes so worker tasks can seek
             * directly to matching pages. A NULL entry is allowed: that
             * column simply falls back to read-and-discard skip in the
             * worker (this can happen with externally-written files that
             * supplied a page index for the predicate column but not for
             * every projected column). */
            if (!slot->col_offset_indexes) {
                slot->col_offset_indexes = carquet_mem_calloc(
                    (size_t)br->num_projected,
                    sizeof(carquet_offset_index_t*));
                if (!slot->col_offset_indexes) {
                    PIPELINE_FAIL(br, CARQUET_ERROR_OUT_OF_MEMORY,
                                  "Row group %d: offset index table allocation failed", target_rg);
                    return;
                }
            }
            for (int32_t i = 0; i < br->num_projected; i++) {
                carquet_error_t oi_err = CARQUET_ERROR_INIT;
                slot->col_offset_indexes[i] = carquet_reader_get_offset_index(
                    br->reader, target_rg, br->projected_columns[i], &oi_err);
                /* NULL is fine — handled by the worker fallback. */
            }
        }

        /* Creation keeps row groups with more rows than file_size*8 off the
         * pipeline (see pipeline_safe); this only holds if that changes. */
        if (br->reader->file_size > 0 &&
            (uint64_t)slot_rows > (uint64_t)br->reader->file_size * 8u) {
            PIPELINE_FAIL(br, CARQUET_ERROR_INVALID_METADATA,
                          "Row group %d: %lld rows is too many to pre-read",
                          target_rg, (long long)slot_rows);
            return;
        }

        /* Ensure column readers exist and are reset for this slot */
        carquet_error_t err = CARQUET_ERROR_INIT;
        for (int32_t i = 0; i < br->num_projected; i++) {
            int32_t file_col_idx = br->projected_columns[i];
            if (slot->col_readers[i]) {
                carquet_column_reader_rebind(
                    slot->col_readers[i], br->reader,
                    target_rg, file_col_idx,
                    br->config.preserve_dictionaries);
            } else {
                slot->col_readers[i] = carquet_reader_get_column(
                    br->reader, target_rg, file_col_idx, &err);
                if (!slot->col_readers[i]) {
                    if (br->pipeline_error.code == CARQUET_OK) {
                        br->pipeline_error = err;
                        if (br->pipeline_error.code == CARQUET_OK) {
                            PIPELINE_FAIL(br, CARQUET_ERROR_INVALID_STATE,
                                          "Row group %d column %d: cannot open column",
                                          target_rg, file_col_idx);
                        }
                    }
                    return;
                }
            }

            /* Take the buffers back from exported Arrow arrays, or leave them
             * to those arrays and start over with fresh ones. */
            slot_reclaim_column(slot, i);

            /* Ensure value buffer is large enough (grow-only via realloc).
             * slot_rows derives from row-group metadata (num_rows), which is
             * attacker-controlled: guard against a negative count and against
             * size_t overflow in the multiply so a malformed file cannot drive
             * a wrapped-around (or absurd) allocation. */
            size_t vsz = br->projected_value_sizes[i];
            if (slot_rows < 0 ||
                (vsz != 0 && (uint64_t)slot_rows > (uint64_t)(SIZE_MAX / vsz))) {
                PIPELINE_FAIL(br, CARQUET_ERROR_INVALID_METADATA,
                              "Row group %d: invalid row count %lld",
                              target_rg, (long long)slot_rows);
                return;
            }
            size_t needed = (size_t)slot_rows * vsz;
            if (needed > slot->col_buf_sizes[i]) {
                void* new_buf = carquet_mem_realloc(slot->col_values[i], needed);
                if (!new_buf) {
                    /* %llu, not %zu: the MinGW/msvcrt printf has no z length modifier. */
                    PIPELINE_FAIL(br, CARQUET_ERROR_OUT_OF_MEMORY,
                                  "Row group %d column %d: cannot allocate %llu bytes",
                                  target_rg, file_col_idx, (unsigned long long)needed);
                    return;
                }
                slot->col_values[i] = new_buf;
                slot->col_buf_sizes[i] = needed;
            }

            /* Validity bitmap for a column with definition levels. */
            if (br->projected_max_defs[i] > 0) {
                size_t bm_needed = ((size_t)slot_rows + 7) / 8;
                if (bm_needed == 0) bm_needed = 1;
                if (bm_needed > slot->col_bitmap_sizes[i]) {
                    uint8_t* new_bm = carquet_mem_realloc(slot->col_bitmaps[i], bm_needed);
                    if (!new_bm) {
                        PIPELINE_FAIL(br, CARQUET_ERROR_OUT_OF_MEMORY,
                                      "Row group %d column %d: cannot allocate %llu-byte bitmap",
                                      target_rg, file_col_idx, (unsigned long long)bm_needed);
                        return;
                    }
                    slot->col_bitmaps[i] = new_bm;
                    slot->col_bitmap_sizes[i] = bm_needed;
                }
            }
            slot_share_column(slot, i);
        }

        slot->rg_index = target_rg;
        slot->ready = false;
        slot->total_rows = slot_rows;
        slot->rows_consumed = 0;
        slot->pending_tasks = 0;
        /* Until its task reports, a column has read nothing: a task that is
         * never submitted must not inherit the previous row group's count. */
        for (int32_t i = 0; i < br->num_projected; i++) {
            slot->col_num_values[i] = 0;
        }

        /* Tasks read the row group straight from the reader's file mapping.
         * Pre-fault its byte range from this (otherwise idle) thread: demand
         * faults taken by eight workers on the same file object serialise in
         * the kernel and cost ~2 us each, one WILLNEED pass maps the cached
         * pages at a fraction of that. A hint only: failure is fine. (A
         * private mmap per row group was measured slower: the extra
         * mmap/munmap pair per slot costs more than it saves.) */
        const uint8_t* slot_data = br->reader->mmap_data;
#if !defined(_WIN32)
        if (slot_data && br->reader->mmap_info && br->reader->mmap_info->fd >= 0) {
            /* Find byte range spanning all projected column chunks */
            int64_t range_lo = INT64_MAX, range_hi = 0;
            for (int32_t i = 0; i < br->num_projected; i++) {
                const parquet_column_metadata_t* cmeta = slot->col_readers[i]->col_meta;
                if (cmeta && cmeta->total_compressed_size > 0 && cmeta->data_page_offset >= 0) {
                    int64_t lo = cmeta->data_page_offset;
                    int64_t hi = lo + cmeta->total_compressed_size;
                    if (lo < range_lo) range_lo = lo;
                    if (hi > range_hi) range_hi = hi;
                }
            }
            if (range_lo < range_hi &&
                (uint64_t)range_hi <= (uint64_t)br->reader->mmap_info->size) {
                /* madvise wants a page-aligned start; the page is 16K on
                 * Apple Silicon and up to 64K on some Linux arm64/ppc64
                 * configs, so query it at runtime. */
                long ps = sysconf(_SC_PAGESIZE);
                int64_t page_mask = (ps > 0) ? (int64_t)ps - 1 : (int64_t)4095;
                int64_t page_lo = range_lo & ~page_mask;
                (void)madvise((void*)(uintptr_t)(slot_data + page_lo),
                              (size_t)(range_hi - page_lo), MADV_WILLNEED);
            }
        }
#endif

        /* Submit one task per compressed column, or up to max_splits
         * page-range segment tasks per column (see coalesced_read_segment). */
        int32_t max_splits_per_col = br->max_splits;
        if (max_splits_per_col > 512) max_splits_per_col = 512;
        for (int32_t i = 0; i < br->num_projected; i++) {
            slot->col_task_first[i] = 0;
            slot->col_task_count[i] = 0;
        }

        /* Bound by task_args space available for this pipeline slot. Indexed
         * by the slot, not by pipeline_count: other slots' tasks may still be
         * running when this one is filled, and their args must stay intact. */
        int32_t tasks_per_slot = br->task_args_capacity / (br->pipeline_depth > 0 ? br->pipeline_depth : 1);
        int32_t base = slot_idx * tasks_per_slot;
        int32_t task_offset = 0;

        int64_t split_offsets[513];
        int64_t split_values_arr[513];

        for (int32_t i = 0; i < br->num_projected; i++) {
            carquet_column_reader_t* cr = slot->col_readers[i];
            size_t value_size = br->projected_value_sizes[i];
            int16_t max_def = br->projected_max_defs[i];

            /* Filtered branch: one task per column, reads matching pages
             * only via the cached offset index. Splitting is not used —
             * range-skipping already constrains the work. (Columns with
             * definition levels never get here: setting a page filter turns
             * such a pipeline off, see carquet_batch_reader_set_page_filter.) */
            if (filter_active) {
                int32_t tidx = base + task_offset++;
                if (tidx >= br->task_args_capacity) break;
                bulk_read_arg_t* t = &br->task_args[tidx];
                task_set_layout(t, 0, NULL, 0, false);
                t->col_reader = cr;
                t->data_base = slot_data;
                t->dest = slot->col_values[i];
                t->max_values = slot_rows;
                t->out_values_read = &slot->col_num_values[i];
                t->start_offset = 0;
                t->end_offset = 0;
                t->local_values_read = 0;
                t->is_split = false;
                t->ranges = &slot->ranges;
                t->offset_index = slot->col_offset_indexes[i];
                t->rg_num_rows = rg_rows;
                t->value_size = value_size;
                slot->col_num_values[i] = slot_rows;
                pipeline_submit(br, slot, t);
                continue;
            }

            /* Page-range segments: decoded straight from the mapping when the
             * coalesced path can take the chunk, otherwise each through its
             * own column reader (dictionary or other encodings, levels). */
            bool coalesce = can_coalesce_column(cr);
            bool reader_segments = !coalesce && slot_data != NULL &&
                cr->max_rep_level == 0 && value_size > 0;
            int32_t nseg = 0;
            int32_t col_splits = max_splits_per_col;
            if (reader_segments) {
                col_splits = dictionary_split_limit(cr, slot_data, col_splits);
            }
            if ((coalesce || reader_segments) && col_splits >= 2) {
                nseg = plan_page_splits(cr, slot_data, rg_rows, col_splits,
                                        split_offsets, split_values_arr);
            }
            if (nseg >= 2 && base + task_offset + nseg > br->task_args_capacity) nseg = 0;
            if (nseg >= 2 && reader_segments &&
                !prepare_segment_readers(br, base + task_offset, nseg, target_rg,
                                         br->projected_columns[i])) {
                nseg = 0;
            }

            if (nseg >= 2) {
                /* The column total is folded from the segment counts once the
                 * slot is waited on (see batch_reader_next_impl). */
                slot->col_num_values[i] = 0;
                slot->col_task_first[i] = base + task_offset;
                slot->col_task_count[i] = nseg;

                for (int32_t s = 0; s < nseg; s++) {
                    int32_t tidx = base + task_offset++;
                    int64_t seg_start_val = split_values_arr[s];
                    int64_t seg_end_val = split_values_arr[s + 1];
                    bulk_read_arg_t* t = &br->task_args[tidx];

                    task_set_layout(t, max_def, slot->col_bitmaps[i], seg_start_val,
                                    reader_segments);
                    t->col_reader = cr;
                    t->data_base = slot_data;
                    t->dest = (uint8_t*)slot->col_values[i] +
                              (size_t)seg_start_val * value_size;
                    t->max_values = seg_end_val - seg_start_val;
                    t->out_values_read = &t->local_values_read;
                    t->start_offset = split_offsets[s];
                    t->end_offset = split_offsets[s + 1];
                    t->local_values_read = 0;
                    t->is_split = true;
                    t->ranges = NULL;
                    t->offset_index = NULL;
                    t->value_size = value_size;
                    pipeline_submit(br, slot, t);
                }
            } else {
                int32_t tidx = base + task_offset++;
                if (tidx >= br->task_args_capacity) break;
                bulk_read_arg_t* t = &br->task_args[tidx];
                task_set_layout(t, max_def, slot->col_bitmaps[i], 0, false);
                t->col_reader = cr;
                t->data_base = slot_data;
                t->dest = slot->col_values[i];
                t->max_values = rg_rows;
                t->out_values_read = &slot->col_num_values[i];
                t->start_offset = 0;
                t->end_offset = 0;
                t->local_values_read = 0;
                t->is_split = false;
                t->ranges = NULL;
                t->offset_index = NULL;
                t->value_size = value_size;
                pipeline_submit(br, slot, t);
            }
        }

        br->pipeline_count++;
        br->rg_order_next++;
    }
}

/* ============================================================================
 * Page Filter — internal helpers
 * ============================================================================ */

static void filter_release_offset_indexes(carquet_batch_reader_t* br) {
    if (!br->projected_offset_indexes) return;
    for (int32_t i = 0; i < br->num_projected; i++) {
        if (br->projected_offset_indexes[i]) {
            carquet_offset_index_free(br->projected_offset_indexes[i]);
            br->projected_offset_indexes[i] = NULL;
        }
    }
    br->projected_oi_rg = -1;
}

static carquet_status_t filter_load_offset_indexes(
    carquet_batch_reader_t* br, int32_t row_group_index,
    carquet_error_t* error) {

    if (br->projected_oi_rg == row_group_index &&
        br->projected_offset_indexes != NULL) {
        return CARQUET_OK;
    }
    filter_release_offset_indexes(br);
    if (!br->projected_offset_indexes) {
        br->projected_offset_indexes = carquet_mem_calloc(
            (size_t)br->num_projected, sizeof(carquet_offset_index_t*));
        if (!br->projected_offset_indexes) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
                "Failed to allocate offset index cache");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
    }
    for (int32_t i = 0; i < br->num_projected; i++) {
        carquet_error_t local = CARQUET_ERROR_INIT;
        br->projected_offset_indexes[i] = carquet_reader_get_offset_index(
            br->reader, row_group_index, br->projected_columns[i], &local);
        /* NULL is allowed (no offset index ⇒ fallback read-and-discard); the
         * skip path is only taken for that column. */
    }
    br->projected_oi_rg = row_group_index;
    return CARQUET_OK;
}

/* Returns the page index in oi that contains logical row `target_row`, or
 * -1 if not found. Sets *page_first_row to that page's first_row_index. */
static int32_t find_page_for_row(
    const carquet_offset_index_t* oi,
    int64_t row_group_num_rows,
    int64_t target_row,
    int64_t* page_first_row_out) {

    int32_t n = carquet_offset_index_num_pages(oi);
    /* Binary search: pages are sorted by first_row_index. */
    int32_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        int32_t mid = lo + (hi - lo) / 2;
        carquet_page_location_t loc;
        if (carquet_offset_index_get_page_location(oi, mid, &loc) != CARQUET_OK) {
            return -1;
        }
        int64_t end_row;
        if (mid + 1 < n) {
            carquet_page_location_t nxt;
            if (carquet_offset_index_get_page_location(oi, mid + 1, &nxt) !=
                CARQUET_OK) {
                return -1;
            }
            end_row = nxt.first_row_index;
        } else {
            end_row = row_group_num_rows;
        }
        if (target_row < loc.first_row_index) {
            hi = mid - 1;
        } else if (target_row >= end_row) {
            lo = mid + 1;
        } else {
            /* A negative first_row_index (malformed OffsetIndex) would make
             * the callers' `target_row - page_first_row` overflow. */
            if (loc.first_row_index < 0) return -1;
            *page_first_row_out = loc.first_row_index;
            return mid;
        }
    }
    return -1;
}

static carquet_status_t position_projected_column(
    carquet_batch_reader_t* br,
    int32_t pi,
    int64_t target_row,
    carquet_error_t* error) {

    carquet_column_reader_t* cr = br->col_readers[pi];
    int32_t file_col = br->projected_columns[pi];
    int64_t rg_num_rows = br->reader->metadata.row_groups[
        br->current_row_group].num_rows;

    carquet_offset_index_t* oi = br->projected_offset_indexes
        ? br->projected_offset_indexes[pi] : NULL;

    if (oi) {
        int64_t page_first_row = 0;
        int32_t page_idx = find_page_for_row(oi, rg_num_rows, target_row,
                                             &page_first_row);
        if (page_idx < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL,
                "Could not locate page for row %lld in column %d",
                (long long)target_row, file_col);
            return CARQUET_ERROR_INTERNAL;
        }
        carquet_page_location_t loc;
        if (carquet_offset_index_get_page_location(oi, page_idx, &loc) != CARQUET_OK) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL,
                "Could not read page location %d in column %d", (int)page_idx, file_col);
            return CARQUET_ERROR_INTERNAL;
        }

        carquet_status_t st = carquet_column_reader_seek_to_data_page(
            cr, loc.offset, 0, error);
        if (st != CARQUET_OK) return st;

        int64_t intra_skip = target_row - page_first_row;
        if (intra_skip > 0) {
            int64_t skipped = carquet_column_skip(cr, intra_skip);
            if (skipped != intra_skip) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL,
                    "Intra-page skip short (column %d): asked %lld, got %lld",
                    file_col, (long long)intra_skip, (long long)skipped);
                return CARQUET_ERROR_INTERNAL;
            }
        }
        return CARQUET_OK;
    }

    /* No offset index for this column: reset to chunk start and read-and-
     * discard up to target_row. Forward-only — backward seeks are handled
     * by the reset. */
    carquet_column_reader_rebind(cr, br->reader,
        br->current_row_group, file_col,
        br->config.preserve_dictionaries);
    if (target_row > 0) {
        int64_t skipped = carquet_column_skip(cr, target_row);
        if (skipped != target_row) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL,
                "Fallback skip short (column %d): asked %lld, got %lld",
                file_col, (long long)target_row, (long long)skipped);
            return CARQUET_ERROR_INTERNAL;
        }
    }
    return CARQUET_OK;
}

/**
 * Advance to a row group that survives both the user's row-group filter
 * and the active page filter. On success, the batch reader is positioned
 * at the first non-empty range of that row group, or returns
 * CARQUET_ERROR_END_OF_DATA when no more matching rows exist.
 */
static carquet_status_t filter_advance_to_active_row_group(
    carquet_batch_reader_t* br, carquet_error_t* error) {

    int32_t num_row_groups = carquet_reader_num_row_groups(br->reader);
    for (;;) {
        if (br->current_row_group < 0) {
            br->current_row_group = 0;
        } else if (!br->filter_rg_state_valid ||
                   br->current_range_index >= br->current_rg_ranges.count) {
            br->current_row_group++;
            br->filter_rg_state_valid = false;
        }
        if (br->current_row_group >= num_row_groups) {
            return CARQUET_ERROR_END_OF_DATA;
        }

        if (br->config.row_group_filter) {
            bool keep = br->config.row_group_filter(br->reader,
                br->current_row_group, br->config.row_group_filter_ctx);
            if (!keep) {
                /* Move on without counting rows toward rows_skipped (user-
                 * level RG filter, not page filter). */
                br->filter_rg_state_valid = false;
                br->current_range_index = br->current_rg_ranges.count;
                continue;
            }
        }

        if (!br->filter_rg_state_valid) {
            carquet_status_t st = carquet_page_filter_eval_row_group(
                br->reader, br->current_row_group,
                br->filter_clauses, br->filter_clause_count,
                &br->current_rg_ranges, error);
            if (st != CARQUET_OK) return st;
            br->filter_rg_state_valid = true;
            br->current_range_index = 0;
            br->current_range_rows_emitted = 0;
            br->range_positioned = false;

            int64_t rg_rows = br->reader->metadata.row_groups[
                br->current_row_group].num_rows;
            br->rows_skipped += rg_rows - br->current_rg_ranges.total_rows;

            /* Open the row group's column readers if we'll need them. */
            if (br->current_rg_ranges.count > 0) {
                carquet_status_t open_st = open_row_group_readers(
                    br, br->current_row_group, error);
                if (open_st != CARQUET_OK) return open_st;
                open_st = filter_load_offset_indexes(br,
                    br->current_row_group, error);
                if (open_st != CARQUET_OK) return open_st;
            } else {
                filter_release_offset_indexes(br);
            }
        }

        if (br->current_range_index < br->current_rg_ranges.count) {
            return CARQUET_OK;
        }
        /* Empty row group ⇒ try the next. */
    }
}

/**
 * Sequential next() path with an active page filter. Reads one range
 * (clipped to batch_size) per call, advancing range/row-group state.
 */
static carquet_status_t batch_reader_next_filtered(
    carquet_batch_reader_t* br, carquet_row_batch_t** batch) {

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_status_t st = filter_advance_to_active_row_group(br, &err);
    if (st != CARQUET_OK) {
        *batch = NULL;
        return st;
    }

    const carquet_row_range_t* range =
        &br->current_rg_ranges.ranges[br->current_range_index];
    int64_t range_remaining = range->num_rows - br->current_range_rows_emitted;
    int64_t batch_size = br->config.batch_size;
    if (batch_size <= 0) batch_size = 65536;
    int64_t rows_to_read = range_remaining < batch_size
        ? range_remaining : batch_size;

    if (!br->range_positioned) {
        int64_t target_row = range->first_row + br->current_range_rows_emitted;
        for (int32_t i = 0; i < br->num_projected; i++) {
            st = position_projected_column(br, i, target_row, &err);
            if (st != CARQUET_OK) return st;
        }
        br->range_positioned = true;
    }

    /* Reuse or allocate batch struct. */
    carquet_row_batch_t* new_batch = br->cached_batch;
    if (!new_batch) {
        new_batch = carquet_mem_calloc(1, sizeof(carquet_row_batch_t));
        if (!new_batch) return CARQUET_ERROR_OUT_OF_MEMORY;
        if (carquet_arena_init(&new_batch->arena) != CARQUET_OK) {
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->columns = carquet_arena_calloc(&new_batch->arena,
            br->num_projected, sizeof(carquet_column_data_t));
        if (!new_batch->columns) {
            carquet_arena_destroy(&new_batch->arena);
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->pooled = true;
        br->cached_batch = new_batch;
    }
    memset(new_batch->columns, 0,
           sizeof(carquet_column_data_t) * br->num_projected);
    new_batch->num_columns = br->num_projected;

    bool read_error = false;
    for (int32_t i = 0; i < br->num_projected; i++) {
        read_projected_column(br, new_batch, i, rows_to_read, true, &read_error);
    }
    if (read_error) {
        return CARQUET_ERROR_DECODE;
    }
    /* Every projected column must have produced the same row count: a
     * column whose chunk is shorter than column 0's would otherwise be
     * served over stale buffer bytes, and a column 0 that produced nothing
     * for a non-empty range would spin forever on empty batches. */
    for (int32_t i = 1; i < br->num_projected; i++) {
        if (new_batch->columns[i].num_values != new_batch->columns[0].num_values) {
            CARQUET_SET_ERROR(&br->last_error, CARQUET_ERROR_DECODE,
                "Column %d produced %lld values but column 0 produced %lld",
                i, (long long)new_batch->columns[i].num_values,
                (long long)new_batch->columns[0].num_values);
            return CARQUET_ERROR_DECODE;
        }
    }
    if (new_batch->columns[0].num_values <= 0 && rows_to_read > 0) {
        CARQUET_SET_ERROR(&br->last_error, CARQUET_ERROR_DECODE,
            "Column 0 produced no values for a range of %lld rows",
            (long long)rows_to_read);
        return CARQUET_ERROR_DECODE;
    }

    new_batch->num_rows = new_batch->columns[0].num_values;
    br->total_rows_read += new_batch->num_rows;

    br->current_range_rows_emitted += new_batch->num_rows;
    if (br->current_range_rows_emitted >= range->num_rows) {
        br->current_range_index++;
        br->current_range_rows_emitted = 0;
        br->range_positioned = false;
    }

    *batch = new_batch;
    return CARQUET_OK;
}

/* ============================================================================
 * Row-range read ([offset, limit))
 * ============================================================================
 * A predicate-independent slice of the logical row space. The global window
 * [row_range_start, row_range_end) is mapped onto each row group's local row
 * space; row groups entirely outside the window are skipped without I/O, and
 * the surviving prefix of a row group is skipped via the offset index (or a
 * read-and-discard fallback) before reading exactly the matching rows.
 */

/* Reuse or allocate the reader's cached batch struct. Returns NULL on OOM. */
static carquet_row_batch_t* ensure_cached_batch(carquet_batch_reader_t* br) {
    carquet_row_batch_t* b = br->cached_batch;
    if (!b) {
        b = carquet_mem_calloc(1, sizeof(carquet_row_batch_t));
        if (!b) return NULL;
        if (carquet_arena_init(&b->arena) != CARQUET_OK) {
            carquet_mem_free(b);
            return NULL;
        }
        b->columns = carquet_arena_calloc(&b->arena,
            br->num_projected, sizeof(carquet_column_data_t));
        if (!b->columns) {
            carquet_arena_destroy(&b->arena);
            carquet_mem_free(b);
            return NULL;
        }
        b->pooled = true;
        br->cached_batch = b;
    }
    memset(b->columns, 0, sizeof(carquet_column_data_t) * br->num_projected);
    b->num_columns = br->num_projected;
    return b;
}

/* Advance current_row_group to the next row group that intersects the active
 * [start,end) window; compute the RG-local window and open + position the
 * projected column readers at its first row. Returns CARQUET_ERROR_END_OF_DATA
 * when the window is exhausted. */
static carquet_status_t row_range_advance(
    carquet_batch_reader_t* br, carquet_error_t* error) {

    int32_t num_rg = carquet_reader_num_row_groups(br->reader);
    for (;;) {
        if (br->current_row_group < 0) {
            br->current_row_group = 0;
            br->rr_rg_global_start = 0;
        } else {
            int64_t prev_rows =
                br->reader->metadata.row_groups[br->current_row_group].num_rows;
            if (prev_rows > 0) br->rr_rg_global_start += prev_rows;
            br->current_row_group++;
        }
        if (br->current_row_group >= num_rg) {
            return CARQUET_ERROR_END_OF_DATA;
        }

        int64_t rg_rows =
            br->reader->metadata.row_groups[br->current_row_group].num_rows;
        if (rg_rows <= 0) continue;

        int64_t rg_gstart = br->rr_rg_global_start;
        int64_t rg_gend = rg_gstart + rg_rows;

        /* Past the window ⇒ no more rows can match. */
        if (rg_gstart >= br->row_range_end) {
            return CARQUET_ERROR_END_OF_DATA;
        }
        /* Entirely before the window ⇒ skip without any I/O. */
        if (rg_gend <= br->row_range_start) {
            continue;
        }
        /* Honour the user's row-group filter, if any. */
        if (br->config.row_group_filter &&
            !br->config.row_group_filter(br->reader, br->current_row_group,
                                         br->config.row_group_filter_ctx)) {
            continue;
        }

        int64_t lo = br->row_range_start - rg_gstart;
        if (lo < 0) lo = 0;
        int64_t hi = br->row_range_end - rg_gstart;
        if (hi > rg_rows) hi = rg_rows;
        if (lo >= hi) continue;  /* defensive: empty intersection */

        br->rr_local_lo = lo;
        br->rr_local_hi = hi;
        br->rr_emitted = 0;

        carquet_status_t st =
            open_row_group_readers(br, br->current_row_group, error);
        if (st != CARQUET_OK) return st;
        st = filter_load_offset_indexes(br, br->current_row_group, error);
        if (st != CARQUET_OK) return st;

        for (int32_t i = 0; i < br->num_projected; i++) {
            st = position_projected_column(br, i, lo, error);
            if (st != CARQUET_OK) return st;
        }
        br->rr_rg_open = true;
        return CARQUET_OK;
    }
}

/* next() driver for an active row range. Reads one window slice (clipped to
 * batch_size) per call, advancing across row groups as windows are consumed. */
static carquet_status_t batch_reader_next_row_range(
    carquet_batch_reader_t* br, carquet_row_batch_t** batch) {

    carquet_error_t err = CARQUET_ERROR_INIT;

    if (!br->rr_rg_open ||
        br->rr_emitted >= (br->rr_local_hi - br->rr_local_lo)) {
        carquet_status_t st = row_range_advance(br, &err);
        if (st != CARQUET_OK) { *batch = NULL; return st; }
    }

    int64_t window = br->rr_local_hi - br->rr_local_lo;
    int64_t remaining = window - br->rr_emitted;
    int64_t batch_size = br->config.batch_size;
    if (batch_size <= 0) batch_size = 65536;
    int64_t rows_to_read = remaining < batch_size ? remaining : batch_size;

    carquet_row_batch_t* new_batch = ensure_cached_batch(br);
    if (!new_batch) return CARQUET_ERROR_OUT_OF_MEMORY;

    bool read_error = false;
    for (int32_t i = 0; i < br->num_projected; i++) {
        read_projected_column(br, new_batch, i, rows_to_read, true, &read_error);
    }
    if (read_error) return CARQUET_ERROR_DECODE;

    new_batch->num_rows = new_batch->columns[0].num_values;
    if (new_batch->num_rows <= 0) {
        /* No forward progress: treat as end to avoid spinning. */
        *batch = NULL;
        return CARQUET_ERROR_END_OF_DATA;
    }
    br->total_rows_read += new_batch->num_rows;
    br->rr_emitted += new_batch->num_rows;
    if (br->rr_emitted >= window) {
        br->rr_rg_open = false;  /* advance to the next window on the next call */
    }

    *batch = new_batch;
    return CARQUET_OK;
}

carquet_status_t carquet_batch_reader_set_row_range(
    carquet_batch_reader_t* reader, int64_t offset, int64_t limit) {

    /* reader is nonnull per API contract. */

    /* offset==0 && limit<0 clears any active range (full sequential read). */
    if (offset == 0 && limit < 0) {
        reader->row_range_active = false;
        reader->rr_rg_open = false;
        reader->current_row_group = -1;
        reader->rows_read_in_group = 0;
        return CARQUET_OK;
    }

    if (offset < 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    /* Mutually exclusive with the page filter (untested interactions). */
    if (reader->filter_clauses && reader->filter_clause_count > 0) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int64_t end;
    if (limit < 0) {
        end = INT64_MAX;                 /* to end of file */
    } else if (offset > INT64_MAX - limit) {
        end = INT64_MAX;                 /* saturate on overflow */
    } else {
        end = offset + limit;
    }

    reader->row_range_active = true;
    reader->row_range_start = offset;
    reader->row_range_end = end;

    /* Restart iteration from the beginning: the window may precede where a
     * prior read already advanced to. */
    reader->rr_rg_open = false;
    reader->rr_emitted = 0;
    reader->rr_rg_global_start = 0;
    reader->current_row_group = -1;
    reader->rows_read_in_group = 0;
    return CARQUET_OK;
}

/* ============================================================================
 * Page Filter — public API
 * ============================================================================ */

carquet_status_t carquet_batch_reader_set_page_filter(
    carquet_batch_reader_t* reader,
    const carquet_filter_clause_t* clauses,
    int32_t count) {

    /* reader is nonnull per API contract. */

    if (clauses == NULL || count <= 0) {
        reader->filter_clauses = NULL;
        reader->filter_clause_count = 0;
        carquet_row_range_list_clear(&reader->current_rg_ranges);
        reader->filter_rg_state_valid = false;
        reader->current_range_index = 0;
        reader->current_range_rows_emitted = 0;
        reader->range_positioned = false;
        filter_release_offset_indexes(reader);
        return CARQUET_OK;
    }

    /* Mutually exclusive with an active row range (untested interactions). */
    if (reader->row_range_active) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Validate every clause up front. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    for (int32_t i = 0; i < count; i++) {
        carquet_status_t st = carquet_page_filter_validate_clause(
            reader->reader, &clauses[i], &err);
        if (st != CARQUET_OK) return st;
    }

    /* If a pipeline is currently active, drain in-flight tasks and drop
     * any pre-read slots whose contents predate the new filter state. */
    if (reader->pipeline_active && reader->pool) {
        carquet_worker_pool_wait(reader->pool);
        for (int32_t s = 0; s < reader->pipeline_depth; s++) {
            rg_slot_t* slot = &reader->pipeline[s];
            if (slot->rg_index >= 0) {
                slot_release_filter_state(slot, reader->num_projected);
                slot->rg_index = -1;
            }
        }
        reader->pipeline_head = 0;
        reader->pipeline_count = 0;
        reader->rg_order_next = 0;
        /* The pipeline's filtered pre-read is REQUIRED-only: with a column
         * that has definition levels, serve the filter from the sequential
         * range iterator instead. */
        for (int32_t i = 0; i < reader->num_projected; i++) {
            if (reader->projected_max_defs[i] > 0) {
                reader->pipeline_active = false;
                break;
            }
        }
    }

    reader->filter_clauses = clauses;
    reader->filter_clause_count = count;
    reader->filter_rg_state_valid = false;
    reader->current_range_index = 0;
    reader->current_range_rows_emitted = 0;
    reader->range_positioned = false;
    filter_release_offset_indexes(reader);

    /* A new filter restarts iteration from the beginning of the file:
     * the predicate may match row groups the previous filter (or
     * unfiltered read) already advanced past. */
    reader->current_row_group = -1;
    reader->rows_read_in_group = 0;
    return CARQUET_OK;
}

int64_t carquet_batch_reader_rows_skipped(
    const carquet_batch_reader_t* reader) {
    /* reader is nonnull per API contract. */
    return reader->rows_skipped;
}

/* ============================================================================
 * Nested batch driver
 * ============================================================================
 * Used when any projected column is repeated (max_rep > 0). Reads a whole row
 * group per batch (the natural granularity that avoids splitting a logical row
 * across batches) and reconstructs list columns via read_nested_list_column().
 * Flat columns in the same projection are read normally. Page filters are not
 * combined with nested reads in this release.
 */
static carquet_status_t batch_reader_next_nested(
    carquet_batch_reader_t* batch_reader,
    carquet_row_batch_t** batch) {

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Advance to the next row group when the current one is exhausted.
     *
     * Repeated columns are read a whole row group at a time, so a page filter
     * is composed at ROW-GROUP granularity: a row group whose statistics prove
     * no row can match is skipped entirely; a row group with any match is read
     * in full (sub-row-group page ranges are not applied to repeated leaves,
     * whose slot counts do not align with logical row ranges). The user-level
     * row_group_filter callback is honoured the same way. */
    bool have_page_filter = batch_reader->filter_clauses &&
                            batch_reader->filter_clause_count > 0;
    /* A row group is emitted as exactly one batch (rows_read_in_group is set
     * once it has been), so advance after that even if the first column
     * still reports values: a malformed chunk whose num_values exceeds what
     * the row group's rows consume would otherwise be re-read forever. */
    if (batch_reader->current_row_group < 0 ||
        batch_reader->rows_read_in_group > 0 ||
        !carquet_column_has_next(batch_reader->col_readers[0])) {

        int32_t num_row_groups = carquet_reader_num_row_groups(batch_reader->reader);
        for (;;) {
            batch_reader->current_row_group++;
            if (batch_reader->current_row_group >= num_row_groups) {
                *batch = NULL;
                return CARQUET_ERROR_END_OF_DATA;
            }
            /* Empty row groups carry no rows to emit (the flat paths skip
             * them too). Returning an empty batch here without consuming the
             * column readers spun forever when the chunk still had values. */
            if (batch_reader->reader->metadata.row_groups[
                    batch_reader->current_row_group].num_rows <= 0) {
                continue;
            }
            if (batch_reader->config.row_group_filter &&
                !batch_reader->config.row_group_filter(
                    batch_reader->reader, batch_reader->current_row_group,
                    batch_reader->config.row_group_filter_ctx)) {
                continue;
            }
            if (have_page_filter) {
                carquet_row_range_list_t ranges;
                carquet_row_range_list_init(&ranges);
                carquet_status_t fst = carquet_page_filter_eval_row_group(
                    batch_reader->reader, batch_reader->current_row_group,
                    batch_reader->filter_clauses, batch_reader->filter_clause_count,
                    &ranges, &err);
                if (fst != CARQUET_OK) {
                    carquet_row_range_list_destroy(&ranges);
                    return fst;
                }
                int64_t matched = ranges.total_rows;
                carquet_row_range_list_destroy(&ranges);
                if (matched == 0) {
                    continue;  /* statistics prove the row group has no match */
                }
            }
            break;
        }
        carquet_status_t status = open_row_group_readers(
            batch_reader, batch_reader->current_row_group, &err);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    /* Reuse or allocate batch struct. */
    carquet_row_batch_t* new_batch = batch_reader->cached_batch;
    if (!new_batch) {
        new_batch = carquet_mem_calloc(1, sizeof(carquet_row_batch_t));
        if (!new_batch) return CARQUET_ERROR_OUT_OF_MEMORY;
        if (carquet_arena_init(&new_batch->arena) != CARQUET_OK) {
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->columns = carquet_arena_calloc(&new_batch->arena,
            batch_reader->num_projected, sizeof(carquet_column_data_t));
        if (!new_batch->columns) {
            carquet_arena_destroy(&new_batch->arena);
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->pooled = true;
        batch_reader->cached_batch = new_batch;
    }

    memset(new_batch->columns, 0,
           sizeof(carquet_column_data_t) * batch_reader->num_projected);
    new_batch->num_columns = batch_reader->num_projected;

    int64_t rg_rows =
        batch_reader->reader->metadata.row_groups[batch_reader->current_row_group].num_rows;
    if (rg_rows <= 0) {
        new_batch->num_rows = 0;
        *batch = new_batch;
        return CARQUET_OK;
    }

    bool read_error = false;
    for (int32_t i = 0; i < batch_reader->num_projected; i++) {
        if (batch_reader->projected_max_reps[i] == 0) {
            read_projected_column(batch_reader, new_batch, i, rg_rows, false, &read_error);
        } else if (batch_reader->projected_max_reps[i] == 1) {
            read_nested_list_column(batch_reader, new_batch, i, rg_rows, &read_error);
        } else {
            CARQUET_SET_ERROR(&err, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Batch reader: nested column depth > 1 (max_rep=%d) not supported",
                batch_reader->projected_max_reps[i]);
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        if (read_error) {
            CARQUET_SET_ERROR(&err, CARQUET_ERROR_INTERNAL,
                "Batch reader: failed to read nested column %d", i);
            return CARQUET_ERROR_INTERNAL;
        }
    }

    new_batch->num_rows = rg_rows;
    batch_reader->total_rows_read += rg_rows;
    batch_reader->rows_read_in_group = rg_rows;
    *batch = new_batch;
    return CARQUET_OK;
}

/* Fold the bitmap bytes the segments of one column shared with their
 * neighbours (see write_segment_validity). Every byte a segment left behind belongs
 * to it and to at most one neighbour; clear them all, then OR the parts. */
static void merge_segment_bitmaps(uint8_t* bitmap, const bulk_read_arg_t* segs,
                                  int32_t nseg) {
    if (!bitmap) return;
    for (int32_t k = 0; k < nseg; k++) {
        if (segs[k].head_byte >= 0) bitmap[segs[k].head_byte] = 0;
        if (segs[k].tail_byte >= 0) bitmap[segs[k].tail_byte] = 0;
    }
    for (int32_t k = 0; k < nseg; k++) {
        if (segs[k].head_byte >= 0) bitmap[segs[k].head_byte] |= segs[k].head_bits;
        if (segs[k].tail_byte >= 0) bitmap[segs[k].tail_byte] |= segs[k].tail_bits;
    }
}

/* A segmented column came up short: re-read the whole chunk here, through
 * the slot's own (untouched) page-by-page reader, which applies every check
 * itself and leaves the detail of a failure in its last_error. Returns the
 * rows read; the slot's buffers are row-aligned again on success. */
static int64_t reread_chunk_row_aligned(carquet_batch_reader_t* br, rg_slot_t* slot,
                                        int32_t i) {
    carquet_column_reader_t* cr = slot->col_readers[i];
    int16_t max_def = br->projected_max_defs[i];
    if (max_def <= 0) {
        return carquet_column_read_batch_ex(cr, slot->col_values[i], slot->total_rows,
                                            NULL, NULL, &cr->last_error);
    }
    int16_t* levels = pool_ensure_def_levels(&br->col_pools[i], (size_t)slot->total_rows);
    if (!levels) {
        CARQUET_SET_ERROR(&cr->last_error, CARQUET_ERROR_OUT_OF_MEMORY,
                          "Cannot allocate %lld definition levels",
                          (long long)slot->total_rows);
        return 0;
    }
    int64_t n = carquet_column_read_batch_ex(cr, slot->col_values[i], slot->total_rows,
                                             levels, NULL, &cr->last_error);
    if (n == slot->total_rows && cr->last_error.code == CARQUET_OK) {
        levels_to_rows((uint8_t*)slot->col_values[i], levels, n, max_def,
                       br->projected_value_sizes[i], slot->col_bitmaps[i]);
    }
    return n;
}

/* Validity bitmap of slot rows [first, first + count) of column i. A view
 * into the slot's bitmap when `first` falls on a byte boundary (the usual
 * case: batch sizes are multiples of 8), else a shifted copy in the
 * column's pooled bitmap. NULL only when that copy cannot be allocated. */
static uint8_t* slot_bitmap_view(carquet_batch_reader_t* br, rg_slot_t* slot,
                                 int32_t i, int64_t first, int64_t count) {
    uint8_t* src = slot->col_bitmaps[i] + (first >> 3);
    unsigned shift = (unsigned)(first & 7);
    if (shift == 0) return src;

    size_t out_bytes = ((size_t)count + 7) / 8;
    size_t src_bytes = ((size_t)slot->total_rows + 7) / 8 - (size_t)(first >> 3);
    uint8_t* out = pool_ensure_bitmap(&br->col_pools[i], out_bytes);
    if (!out) return NULL;
    for (size_t j = 0; j < out_bytes; j++) {
        unsigned v = (unsigned)src[j] >> shift;
        if (j + 1 < src_bytes) v |= (unsigned)src[j + 1] << (8 - shift);
        out[j] = (uint8_t)v;
    }
    return out;
}

static carquet_status_t batch_reader_next_impl(
    carquet_batch_reader_t* batch_reader,
    carquet_row_batch_t** batch) {

    /* batch_reader and batch are nonnull per API contract */
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Row-range read ([offset,limit)) drives its own sequential path and is
     * not combined with repeated columns in this release. */
    if (batch_reader->row_range_active) {
        if (batch_reader->has_repeated) {
            CARQUET_SET_ERROR(&err, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Row-range read is not supported for repeated (nested) columns");
            *batch = NULL;
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        return batch_reader_next_row_range(batch_reader, batch);
    }

    /* Repeated (LIST/MAP-leaf) columns use the dedicated nested driver, which
     * reconstructs Arrow list layout a whole row group at a time. */
    if (batch_reader->has_repeated) {
        return batch_reader_next_nested(batch_reader, batch);
    }

    /* ====================================================================
     * FILTERED + SEQUENTIAL PATH (page filter active, pipeline disabled)
     *
     * The pipeline path's filtered variant (pipeline_fill below) handles
     * the case where pipeline_active is true: it pre-reads only matching
     * pages into the slot buffers and is then served by the pipeline
     * fast path further down. When the pipeline is not active (e.g.
     * single-row-group, uncompressed, OPTIONAL columns), we drive the
     * sequential range-iterator instead.
     * ==================================================================== */
    if (batch_reader->filter_clauses &&
        batch_reader->filter_clause_count > 0 &&
        !batch_reader->pipeline_active) {
        return batch_reader_next_filtered(batch_reader, batch);
    }

    /* ====================================================================
     * PIPELINE FAST PATH: serve pre-read data directly from ring buffer
     * ====================================================================
     * When pipeline is active, ALL column data has been bulk-read into
     * contiguous buffers by worker pool threads. We just memcpy batches
     * from those buffers. No column readers, no per-page overhead. */
    if (batch_reader->pipeline_active) {
        if (batch_reader->pipeline_dead) {
            *batch = NULL;
            return batch_reader->pipeline_error.code;
        }
        /* Check if we need to advance to the next pipeline slot */
        rg_slot_t* slot = NULL;
        if (batch_reader->pipeline_count > 0) {
            slot = &batch_reader->pipeline[batch_reader->pipeline_head];
            if (slot->rows_consumed >= slot->total_rows) {
                /* Current slot exhausted — retire it and advance */
                slot->rg_index = -1;
                slot_release_filter_state(slot, batch_reader->num_projected);
                batch_reader->pipeline_head = (batch_reader->pipeline_head + 1) % batch_reader->pipeline_depth;
                batch_reader->pipeline_count--;
                slot = NULL;
            }
        }

        if (!slot) {
            /* Fill and wait for new pipeline slots */
            pipeline_fill(batch_reader);
            if (batch_reader->pipeline_count == 0) {
                *batch = NULL;
                if (batch_reader->pipeline_error.code != CARQUET_OK) {
                    batch_reader->pipeline_dead = true;
                    return batch_reader->pipeline_error.code;
                }
                return CARQUET_ERROR_END_OF_DATA;
            }
            slot = &batch_reader->pipeline[batch_reader->pipeline_head];
            pipeline_wait_slot(batch_reader, slot);

            /* Columns read as page-range segments: fold the per-segment
             * counts into the column total. A segment that bailed leaves the
             * sum short; re-read the whole chunk here, on this thread, through
             * the page-by-page reader, which applies every check itself and
             * reports the detail of whatever stopped the fast path (the same
             * recovery a whole-chunk task performs inline). */
            for (int32_t i = 0; i < batch_reader->num_projected; i++) {
                if (slot->col_task_count[i] <= 1) continue;
                const bulk_read_arg_t* segs = &batch_reader->task_args[slot->col_task_first[i]];
                int32_t nseg = slot->col_task_count[i];
                int64_t sum = 0;
                for (int32_t k = 0; k < nseg; k++) {
                    sum += segs[k].local_values_read;
                }
                if (sum == slot->total_rows) {
                    slot->col_num_values[i] = sum;
                    merge_segment_bitmaps(slot->col_bitmaps[i], segs, nseg);
                    continue;
                }
                slot->col_num_values[i] = reread_chunk_row_aligned(batch_reader, slot, i);
            }

            /* Every column must have produced the slot's full row count before
             * any of it is exposed: a failed or truncated column leaves the
             * rest of its buffer uninitialized or holding an older row group. */
            for (int32_t i = 0; i < batch_reader->num_projected; i++) {
                if (slot->col_num_values[i] == slot->total_rows) continue;
                const carquet_column_reader_t* cr = slot->col_readers[i];
                if (cr && cr->last_error.code != CARQUET_OK) {
                    batch_reader->pipeline_error = cr->last_error;
                } else {
                    batch_reader->pipeline_error.code = CARQUET_OK;
                    CARQUET_SET_ERROR(&batch_reader->pipeline_error, CARQUET_ERROR_DECODE,
                        "Row group %d column %d: read %lld of %lld values",
                        slot->rg_index, batch_reader->projected_columns[i],
                        (long long)slot->col_num_values[i], (long long)slot->total_rows);
                }
                batch_reader->pipeline_dead = true;
                *batch = NULL;
                return batch_reader->pipeline_error.code;
            }

            /* Refill freed slots for next round */
            pipeline_fill(batch_reader);
        }

        /* Reuse or allocate batch struct */
        carquet_row_batch_t* new_batch = batch_reader->cached_batch;
        if (!new_batch) {
            new_batch = carquet_mem_calloc(1, sizeof(carquet_row_batch_t));
            if (!new_batch) return CARQUET_ERROR_OUT_OF_MEMORY;
            if (carquet_arena_init(&new_batch->arena) != CARQUET_OK) {
                carquet_mem_free(new_batch);
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            new_batch->columns = carquet_arena_calloc(&new_batch->arena,
                batch_reader->num_projected, sizeof(carquet_column_data_t));
            if (!new_batch->columns) {
                carquet_arena_destroy(&new_batch->arena);
                carquet_mem_free(new_batch);
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            new_batch->pooled = true;
            batch_reader->cached_batch = new_batch;
        }

        memset(new_batch->columns, 0, sizeof(carquet_column_data_t) * batch_reader->num_projected);
        new_batch->num_columns = batch_reader->num_projected;

        int64_t remaining = slot->total_rows - slot->rows_consumed;
        int64_t rows_to_read = remaining > batch_reader->config.batch_size
                             ? batch_reader->config.batch_size : remaining;

        /* Copy data from pre-read buffers into batch (zero-copy view) */
        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            carquet_column_data_t* col = &new_batch->columns[i];
            size_t vs = batch_reader->projected_value_sizes[i];
            size_t offset = (size_t)slot->rows_consumed * vs;

            col->data = (uint8_t*)slot->col_values[i] + offset;
            col->num_values = rows_to_read;
            col->type = batch_reader->projected_types[i];
            col->type_length = batch_reader->projected_type_lengths[i];
            col->ownership = CARQUET_DATA_VIEW;
            col->null_bitmap = NULL;  /* REQUIRED columns — no nulls */
            if (batch_reader->projected_max_defs[i] > 0) {
                col->null_bitmap = slot_bitmap_view(batch_reader, slot, i,
                                                    slot->rows_consumed, rows_to_read);
                if (!col->null_bitmap) return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            /* Offer the storage to an Arrow export when the batch covers a
             * good part of the row group: a small batch kept alive would pin
             * the whole row group's buffer, so those are copied. */
            col->share = rows_to_read * 4 >= slot->total_rows ? slot->col_share[i] : NULL;
            col->bitmap_shared = col->null_bitmap != NULL && (slot->rows_consumed & 7) == 0;
        }

        slot->rows_consumed += rows_to_read;
        new_batch->num_rows = rows_to_read;
        batch_reader->total_rows_read += rows_to_read;
        *batch = new_batch;
        return CARQUET_OK;
    }

    /* ====================================================================
     * SEQUENTIAL PATH (non-mmap, uncompressed, or single RG)
     * ==================================================================== */

    /* Check if we need to move to next row group */
    if (batch_reader->current_row_group < 0 ||
        !carquet_column_has_next(batch_reader->col_readers[0])) {

        int32_t num_row_groups = carquet_reader_num_row_groups(batch_reader->reader);
        batch_reader->current_row_group++;
        if (batch_reader->current_row_group >= num_row_groups) {
            *batch = NULL;
            return CARQUET_ERROR_END_OF_DATA;
        }

        /* Apply row group filter */
        while (batch_reader->config.row_group_filter) {
            bool should_read = batch_reader->config.row_group_filter(
                batch_reader->reader, batch_reader->current_row_group,
                batch_reader->config.row_group_filter_ctx);
            if (should_read) break;
            batch_reader->current_row_group++;
            if (batch_reader->current_row_group >= num_row_groups) {
                *batch = NULL;
                return CARQUET_ERROR_END_OF_DATA;
            }
        }

        carquet_status_t status = open_row_group_readers(
            batch_reader, batch_reader->current_row_group, &err);
        if (status != CARQUET_OK) {
            return status;
        }
    }

    /* Reuse or allocate batch struct */
    carquet_row_batch_t* new_batch = batch_reader->cached_batch;
    if (!new_batch) {
        new_batch = carquet_mem_calloc(1, sizeof(carquet_row_batch_t));
        if (!new_batch) {
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (carquet_arena_init(&new_batch->arena) != CARQUET_OK) {
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->columns = carquet_arena_calloc(&new_batch->arena,
            batch_reader->num_projected, sizeof(carquet_column_data_t));
        if (!new_batch->columns) {
            carquet_arena_destroy(&new_batch->arena);
            carquet_mem_free(new_batch);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        new_batch->pooled = true;
        batch_reader->cached_batch = new_batch;
    }

    /* Reset column data for this batch */
    memset(new_batch->columns, 0, sizeof(carquet_column_data_t) * batch_reader->num_projected);
    new_batch->num_columns = batch_reader->num_projected;

    int64_t batch_size = batch_reader->config.batch_size;
    int64_t rows_to_read = carquet_column_remaining(batch_reader->col_readers[0]);
    if (rows_to_read > batch_size) {
        rows_to_read = batch_size;
    }

    /* Handle empty row group - return empty batch, not an error */
    if (rows_to_read == 0) {
        new_batch->num_rows = 0;
        *batch = new_batch;
        return CARQUET_OK;
    }

    /* Read each column - potentially in parallel */
    bool read_error = false;

    /* Uncompressed fixed-width mmap columns can often be served entirely as
     * direct page views. Pre-load those pages serially before the parallel
     * decision so we can clamp the batch to the current page window and avoid
     * both the zero-byte peek copy path and the column-parallel barrier cost. */
    {
        bool zero_copy_candidates = true;
        for (int32_t zi = 0; zi < batch_reader->num_projected; zi++) {
            if (!column_is_zero_copy_candidate(
                    batch_reader->col_readers[zi],
                    batch_reader->projected_types[zi],
                    batch_reader->projected_max_defs[zi])) {
                zero_copy_candidates = false;
                break;
            }
        }

        if (zero_copy_candidates) {
            for (int32_t zi = 0; zi < batch_reader->num_projected; zi++) {
                carquet_column_reader_t* col_reader = batch_reader->col_readers[zi];
                if (col_reader && col_reader->values_remaining > 0) {
                    carquet_status_t status = carquet_column_ensure_page_loaded(col_reader, &err);
                    if (status != CARQUET_OK) {
                        /* Hand the page-level message to
                         * carquet_batch_reader_last_error(); this path returns
                         * before any read_projected_column() could record it. */
                        col_reader->last_error = err;
                        return status;
                    }
                }
            }
        }
    }

    /* ========================================================================
     * PARALLEL PAGE PREFETCH PHASE
     * ========================================================================
     * Pre-load pages for ALL columns in parallel BEFORE reading.
     * Uses persistent worker pool (no per-batch fork/join overhead) when
     * available, falls back to OpenMP, then serial.
     *
     * Only parallelize when column readers can load pages concurrently
     * (mmap or positional reads; a shared fseek+fread stream cannot) and
     * when there are columns needing decompression (uncompressed pages are
     * trivial). */
    bool parallel_io = br_parallel_io(batch_reader->reader);
    bool needs_decompression = false;
    for (int32_t pi = 0; pi < batch_reader->num_projected; pi++) {
        carquet_column_reader_t* cr = batch_reader->col_readers[pi];
        if (cr && cr->col_meta &&
            cr->col_meta->codec != CARQUET_COMPRESSION_UNCOMPRESSED) {
            needs_decompression = true;
            break;
        }
    }

#ifdef _OPENMP
    {
        int num_threads = batch_reader->config.num_threads;
        if (num_threads <= 0) num_threads = omp_get_max_threads();
        if (num_threads > batch_reader->num_projected) num_threads = batch_reader->num_projected;
        if (num_threads < 1) num_threads = 1;

        int32_t omp_i;
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1) if(parallel_io && needs_decompression && num_threads > 1)
        for (omp_i = 0; omp_i < batch_reader->num_projected; omp_i++) {
            carquet_column_reader_t* col_reader = batch_reader->col_readers[omp_i];
            if (col_reader && !col_reader->page_loaded && col_reader->values_remaining > 0) {
                (void)carquet_column_read_batch(col_reader, NULL, 0, NULL, NULL);
            }
        }
    }
#else
    bool big_batch = rows_to_read * (int64_t)batch_reader->num_projected >= CARQUET_PARALLEL_MIN_VALUES;
    if (parallel_io && needs_decompression && batch_reader->num_projected > 1 && big_batch &&
        br_column_pool(batch_reader)) {
        br_column_parallel_for(batch_reader, br_prefetch_task, NULL, 0);
    } else {
        for (int32_t pi = 0; pi < batch_reader->num_projected; pi++) {
            carquet_column_reader_t* col_reader = batch_reader->col_readers[pi];
            if (col_reader && !col_reader->page_loaded && col_reader->values_remaining > 0) {
                (void)carquet_column_read_batch(col_reader, NULL, 0, NULL, NULL);
            }
        }
    }
#endif

    /* If every projected column is backed by a direct page view, trim the
     * batch to the smallest currently available page slice. This avoids
     * copying across page boundaries and lets the main read phase stay on
     * the zero-copy path even when the requested batch size is larger than
     * an individual page. */
    {
        int64_t zero_copy_rows = clamp_rows_to_zero_copy_window(batch_reader, rows_to_read);
        if (zero_copy_rows > 0 && zero_copy_rows < rows_to_read) {
            rows_to_read = zero_copy_rows;
        }
    }

    /* ========================================================================
     * MAIN COLUMN READING PHASE
     * ========================================================================
     * Read from pre-loaded pages. Since pages are already decompressed,
     * this phase is mostly memory copies / zero-copy pointer setup.
     * Uses pooled buffers to avoid per-batch malloc/free.
     *
     * Worker pool is used for parallel column reading when pool is available
     * and columns need non-trivial work. Otherwise serial (which is often
     * optimal for zero-copy columns where read_projected_column is ~free).
     */
    bool all_zero_copy_ready = true;
    for (int32_t zi = 0; zi < batch_reader->num_projected; zi++) {
        if (!column_can_zero_copy_batch(
                batch_reader->col_readers[zi],
                batch_reader->projected_types[zi],
                batch_reader->projected_max_defs[zi],
                rows_to_read)) {
            all_zero_copy_ready = false;
            break;
        }
    }


    int32_t col_i;
#ifdef _OPENMP
    {
        int num_threads_read = batch_reader->config.num_threads;
        if (num_threads_read <= 0) num_threads_read = omp_get_max_threads();
        if (num_threads_read > batch_reader->num_projected) num_threads_read = batch_reader->num_projected;
        if (num_threads_read < 1) num_threads_read = 1;

        bool can_par = parallel_io && (num_threads_read > 1) &&
                       (batch_reader->num_projected > 1) && !all_zero_copy_ready;
        if (can_par) {
            /* Per-column error slots: each thread writes only its own slot, so
             * the failure flag is never a shared write across threads (no data
             * race). Reduce into read_error after the region. */
            bool* col_err = carquet_mem_calloc(
                (size_t)batch_reader->num_projected, sizeof(bool));
            if (!col_err) {
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            #pragma omp parallel for num_threads(num_threads_read) schedule(dynamic, 1)
            for (col_i = 0; col_i < batch_reader->num_projected; col_i++) {
                read_projected_column(batch_reader, new_batch, col_i, rows_to_read, true, &col_err[col_i]);
            }
            for (col_i = 0; col_i < batch_reader->num_projected; col_i++) {
                if (col_err[col_i]) read_error = true;
            }
            carquet_mem_free(col_err);
        } else {
            for (col_i = 0; col_i < batch_reader->num_projected; col_i++) {
                read_projected_column(batch_reader, new_batch, col_i, rows_to_read, true, &read_error);
            }
        }
    }
#else
    if (parallel_io && batch_reader->num_projected > 1 && !all_zero_copy_ready && big_batch &&
        br_column_pool(batch_reader)) {
        br_column_parallel_for(batch_reader, br_read_task, new_batch, rows_to_read);
        for (col_i = 0; col_i < batch_reader->num_projected; col_i++) {
            if (batch_reader->col_tasks[col_i].error) read_error = true;
        }
    } else {
        for (col_i = 0; col_i < batch_reader->num_projected; col_i++) {
            read_projected_column(batch_reader, new_batch, col_i, rows_to_read, true, &read_error);
        }
    }
#endif

    if (read_error) {
        /* Don't free cached_batch, just return error */
        return CARQUET_ERROR_DECODE;
    }
    /* See the filtered path: all columns must agree on the row count. */
    for (col_i = 1; col_i < batch_reader->num_projected; col_i++) {
        if (new_batch->columns[col_i].num_values != new_batch->columns[0].num_values) {
            CARQUET_SET_ERROR(&batch_reader->last_error, CARQUET_ERROR_DECODE,
                "Column %d produced %lld values but column 0 produced %lld",
                col_i, (long long)new_batch->columns[col_i].num_values,
                (long long)new_batch->columns[0].num_values);
            return CARQUET_ERROR_DECODE;
        }
    }

    new_batch->num_rows = new_batch->columns[0].num_values;
    batch_reader->total_rows_read += new_batch->num_rows;

    *batch = new_batch;
    return CARQUET_OK;
}

/* Clear the previous batch's failure detail. Columns not touched by this batch
 * would otherwise keep a stale message from an earlier one. */
static void batch_reader_clear_errors(carquet_batch_reader_t* batch_reader) {
    batch_reader->last_error.code = CARQUET_OK;
    batch_reader->last_error.message[0] = '\0';
    if (!batch_reader->col_readers) return;
    for (int32_t i = 0; i < batch_reader->num_projected; i++) {
        carquet_column_reader_t* cr = batch_reader->col_readers[i];
        if (cr) {
            cr->last_error.code = CARQUET_OK;
            cr->last_error.message[0] = '\0';
        }
    }
}

/* First column reader that recorded a failure during this batch, or NULL. */
static const carquet_error_t* batch_reader_first_column_error(
    const carquet_batch_reader_t* batch_reader) {
    if (!batch_reader->col_readers) return NULL;
    for (int32_t i = 0; i < batch_reader->num_projected; i++) {
        const carquet_column_reader_t* cr = batch_reader->col_readers[i];
        if (cr && cr->last_error.code != CARQUET_OK) {
            return &cr->last_error;
        }
    }
    return NULL;
}

carquet_status_t carquet_batch_reader_next(
    carquet_batch_reader_t* batch_reader,
    carquet_row_batch_t** batch) {

    batch_reader_clear_errors(batch_reader);

    carquet_status_t status = batch_reader_next_impl(batch_reader, batch);

    if (status == CARQUET_OK && *batch) {
        (*batch)->projection = batch_reader->projected_columns;
#ifdef _OPENMP
        (*batch)->pool = NULL;   /* the exporter parallelises with OpenMP */
#else
        bool big = (*batch)->num_rows * (int64_t)batch_reader->num_projected >= CARQUET_PARALLEL_MIN_VALUES;
        (*batch)->pool = big ? br_column_pool(batch_reader) : NULL;
#endif
    }

    /* END_OF_DATA is the normal terminator, not a failure. */
    if (status == CARQUET_OK || status == CARQUET_ERROR_END_OF_DATA) {
        return status;
    }

    /* Prefer the detailed message recorded by whichever column failed; the
     * generic fallback keeps last_error consistent with the returned status
     * for the paths that fail before or outside a column read. */
    const carquet_error_t* col_err = batch_reader_first_column_error(batch_reader);
    if (batch_reader->pipeline_dead && batch_reader->pipeline_error.code == status) {
        batch_reader->last_error = batch_reader->pipeline_error;
    } else if (col_err) {
        batch_reader->last_error = *col_err;
    } else {
        CARQUET_SET_ERROR(&batch_reader->last_error, status, "%s",
                          carquet_status_string(status));
    }
    return status;
}

const carquet_error_t* carquet_batch_reader_last_error(
    const carquet_batch_reader_t* batch_reader) {
    if (!batch_reader || batch_reader->last_error.code == CARQUET_OK) {
        return NULL;
    }
    return &batch_reader->last_error;
}

void carquet_batch_reader_free(carquet_batch_reader_t* batch_reader) {
    if (!batch_reader) return;

    /* Drain any in-flight pipeline tasks before freeing */
    if (batch_reader->pool) {
        carquet_worker_pool_wait(batch_reader->pool);
    }
    if (batch_reader->slot_sync_init) {
#ifdef _WIN32
        DeleteCriticalSection(&batch_reader->slot_mutex);
#else
        pthread_mutex_destroy(&batch_reader->slot_mutex);
        pthread_cond_destroy(&batch_reader->slot_cond);
#endif
        batch_reader->slot_sync_init = false;
    }

    /* Free pipeline slots */
    if (batch_reader->pipeline) {
        for (int32_t s = 0; s < batch_reader->pipeline_depth; s++) {
            rg_slot_t* slot = &batch_reader->pipeline[s];
            slot_release_filter_state(slot, batch_reader->num_projected);
            if (slot->col_offset_indexes) {
                carquet_mem_free(slot->col_offset_indexes);
                slot->col_offset_indexes = NULL;
            }
            if (slot->col_readers) {
                for (int32_t i = 0; i < batch_reader->num_projected; i++) {
                    if (slot->col_readers[i]) {
                        carquet_column_reader_free(slot->col_readers[i]);
                    }
                }
                carquet_mem_free(slot->col_readers);
            }
            /* Buffers still referenced by exported Arrow arrays stay alive;
             * the last release frees them. */
            for (int32_t i = 0; slot->col_share && i < batch_reader->num_projected; i++) {
                slot_reclaim_column(slot, i);
            }
            carquet_mem_free(slot->col_share);
            if (slot->col_values) {
                for (int32_t i = 0; i < batch_reader->num_projected; i++) {
                    carquet_mem_free(slot->col_values[i]);
                }
                carquet_mem_free(slot->col_values);
            }
            carquet_mem_free(slot->col_buf_sizes);
            carquet_mem_free(slot->col_num_values);
            carquet_mem_free(slot->col_task_first);
            carquet_mem_free(slot->col_task_count);
            if (slot->col_bitmaps) {
                for (int32_t i = 0; i < batch_reader->num_projected; i++) {
                    carquet_mem_free(slot->col_bitmaps[i]);
                }
                carquet_mem_free(slot->col_bitmaps);
            }
            carquet_mem_free(slot->col_bitmap_sizes);
        }
        carquet_mem_free(batch_reader->pipeline);
    }

    /* Free per-reader task args and the segment scratch pool (all tasks have
     * finished, so every scratch buffer is back in the pool). */
    if (batch_reader->task_args) {
        for (int32_t i = 0; i < batch_reader->task_args_capacity; i++) {
            carquet_mem_free(batch_reader->task_args[i].levels);
            carquet_mem_free(batch_reader->task_args[i].vbits);
            if (batch_reader->task_args[i].seg_reader) {
                carquet_column_reader_free(batch_reader->task_args[i].seg_reader);
            }
        }
    }
    carquet_mem_free(batch_reader->task_args);
    if (batch_reader->scratch_bufs) {
        for (int32_t i = 0; i < batch_reader->scratch_cap; i++) {
            carquet_mem_free(batch_reader->scratch_bufs[i]);
        }
        carquet_mem_free(batch_reader->scratch_bufs);
        carquet_mem_free(batch_reader->scratch_sizes);
    }

    /* Destroy worker pool (only if we created it) */
    if (!batch_reader->pool_is_borrowed) {
        carquet_worker_pool_destroy(batch_reader->pool);
    }
    if (batch_reader->col_pool) {
        carquet_worker_pool_destroy(batch_reader->col_pool);
    }
    carquet_mem_free(batch_reader->col_tasks);
    carquet_mem_free(batch_reader->col_task_args);

    /* Free column readers */
    if (batch_reader->col_readers) {
        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            if (batch_reader->col_readers[i]) {
                carquet_column_reader_free(batch_reader->col_readers[i]);
            }
        }
        carquet_mem_free(batch_reader->col_readers);
    }

    /* Free buffer pools */
    if (batch_reader->col_pools) {
        for (int32_t i = 0; i < batch_reader->num_projected; i++) {
            carquet_mem_free(batch_reader->col_pools[i].data);
            carquet_mem_free(batch_reader->col_pools[i].null_bitmap);
            carquet_mem_free(batch_reader->col_pools[i].def_levels);
            carquet_mem_free(batch_reader->col_pools[i].rep_levels);
            carquet_mem_free(batch_reader->col_pools[i].list_offsets);
            carquet_mem_free(batch_reader->col_pools[i].list_validity);
        }
        carquet_mem_free(batch_reader->col_pools);
    }

    /* Free cached batch struct (but NOT pool buffers - those are freed above) */
    if (batch_reader->cached_batch) {
        carquet_arena_destroy(&batch_reader->cached_batch->arena);
        carquet_mem_free(batch_reader->cached_batch);
    }

    /* Filter state cleanup */
    filter_release_offset_indexes(batch_reader);
    carquet_mem_free(batch_reader->projected_offset_indexes);
    carquet_row_range_list_destroy(&batch_reader->current_rg_ranges);

    carquet_mem_free(batch_reader->rg_order);
    carquet_mem_free(batch_reader->projected_value_sizes);
    carquet_mem_free(batch_reader->projected_max_reps);
    carquet_mem_free(batch_reader->projected_max_defs);
    carquet_mem_free(batch_reader->projected_type_lengths);
    carquet_mem_free(batch_reader->projected_types);
    carquet_mem_free(batch_reader->projected_columns);
    carquet_mem_free(batch_reader);
}

/* ============================================================================
 * Public Thread Pool API
 * ============================================================================
 */

carquet_thread_pool_t* carquet_thread_pool_create(int32_t num_threads) {
    if (num_threads <= 0) {
#ifdef _OPENMP
        num_threads = omp_get_max_threads();
#else
        num_threads = 4;
#endif
    }
    if (num_threads < 2) num_threads = 2;
    return (carquet_thread_pool_t*)carquet_worker_pool_create(num_threads);
}

void carquet_thread_pool_destroy(carquet_thread_pool_t* pool) {
    carquet_worker_pool_destroy((carquet_worker_pool_t*)pool);
}

/* ============================================================================
 * Row Batch Implementation
 * ============================================================================
 */

int64_t carquet_row_batch_num_rows(const carquet_row_batch_t* batch) {
    /* batch is nonnull per API contract */
    return batch->num_rows;
}

int32_t carquet_row_batch_num_columns(const carquet_row_batch_t* batch) {
    /* batch is nonnull per API contract */
    return batch->num_columns;
}

const int32_t* carquet_row_batch_projection(const carquet_row_batch_t* batch) {
    return batch ? batch->projection : NULL;
}

struct carquet_worker_pool* carquet_row_batch_pool(const carquet_row_batch_t* batch) {
    return batch ? batch->pool : NULL;
}

carquet_status_t carquet_row_batch_column(
    const carquet_row_batch_t* batch,
    int32_t column_index,
    const void** data,
    const uint8_t** null_bitmap,
    int64_t* num_values) {

    /* batch, data, null_bitmap, num_values are nonnull per API contract */
    if (column_index < 0 || column_index >= batch->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const carquet_column_data_t* col = &batch->columns[column_index];

    /* When preserve_dictionaries is enabled, col->data holds uint32_t indices,
     * not materialized values. Returning it through the value accessor would
     * hand the caller indices silently mis-cast as the column's physical type.
     * Force the caller to use carquet_row_batch_column_dictionary() instead. */
    if (col->is_dictionary) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* List (repeated) columns must be read via carquet_row_batch_column_list();
     * returning the flattened child through the flat accessor would silently
     * drop the list structure. */
    if (col->list_offsets) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    *data = col->data;
    *null_bitmap = col->null_bitmap;
    *num_values = col->num_values;

    return CARQUET_OK;
}

carquet_status_t carquet_row_batch_column_list(
    const carquet_row_batch_t* batch,
    int32_t column_index,
    const int32_t** offsets,
    int64_t* num_lists,
    const void** values,
    const uint8_t** value_validity,
    int64_t* num_values,
    const uint8_t** list_validity) {

    if (column_index < 0 || column_index >= batch->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    const carquet_column_data_t* col = &batch->columns[column_index];
    if (!col->list_offsets) {
        return CARQUET_ERROR_INVALID_ARGUMENT;  /* not a list column */
    }

    *offsets = col->list_offsets;
    *num_lists = col->num_lists;
    *values = col->data;
    *value_validity = col->null_bitmap;
    *num_values = col->num_values;
    if (list_validity) {
        *list_validity = col->list_validity;
    }
    return CARQUET_OK;
}

carquet_status_t carquet_row_batch_column_dictionary(
    const carquet_row_batch_t* batch,
    int32_t column_index,
    const uint32_t** indices,
    const uint8_t** null_bitmap,
    int64_t* num_values,
    const uint8_t** dictionary_data,
    int32_t* dictionary_count,
    const uint32_t** dictionary_offsets) {

    if (column_index < 0 || column_index >= batch->num_columns) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const carquet_column_data_t* col = &batch->columns[column_index];

    if (!col->is_dictionary) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    *indices = (const uint32_t*)col->data;
    *null_bitmap = col->null_bitmap;
    *num_values = col->num_values;
    *dictionary_data = col->dictionary_data;
    *dictionary_count = col->dictionary_count;
    if (dictionary_offsets) {
        *dictionary_offsets = col->dictionary_offsets;
    }

    return CARQUET_OK;
}

/* Internal, for the Arrow exporter (arrow_c_export.c): a new reference on
 * the storage behind batch column `column_index`, or NULL when the column
 * must be copied. *bitmap_shared tells whether its null bitmap lies in that
 * storage too. Release the reference with carquet_row_batch_share_release. */
void* carquet_row_batch_column_share(const carquet_row_batch_t* batch,
                                     int32_t column_index, bool* bitmap_shared) {
    *bitmap_shared = false;
#if BR_SHARE_ATOMICS
    if (!batch || column_index < 0 || column_index >= batch->num_columns) return NULL;
    const carquet_column_data_t* col = &batch->columns[column_index];
    if (!col->share || col->is_dictionary || col->list_offsets) return NULL;
    br_share_inc(&col->share->refs);
    *bitmap_shared = col->bitmap_shared;
    return col->share;
#else
    (void)batch; (void)column_index;
    return NULL;
#endif
}

void carquet_row_batch_share_release(void* share) {
#if BR_SHARE_ATOMICS
    br_share_t* sh = (br_share_t*)share;
    if (sh && br_share_dec(&sh->refs) == 0) br_share_free(sh);
#else
    (void)share;
#endif
}

void carquet_row_batch_free(carquet_row_batch_t* batch) {
    if (!batch) return;

    /* Pooled batches are owned by the batch_reader - don't free data */
    if (batch->pooled) {
        /* Data buffers belong to the batch_reader's pool.
         * The batch struct itself is cached and reused.
         * This is a no-op - the caller should just drop the pointer. */
        return;
    }

    /* Non-pooled batch: free column data (only if owned, not views into mmap) */
    for (int32_t i = 0; i < batch->num_columns; i++) {
        if (batch->columns[i].ownership == CARQUET_DATA_OWNED) {
            carquet_mem_free(batch->columns[i].data);
        }
        carquet_mem_free(batch->columns[i].null_bitmap);
    }

    carquet_arena_destroy(&batch->arena);
    carquet_mem_free(batch);
}
