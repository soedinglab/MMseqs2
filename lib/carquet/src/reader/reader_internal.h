/**
 * @file reader_internal.h
 * @brief Internal reader structures
 *
 * This header defines internal structures that are shared between
 * reader components but not exposed in the public API.
 */

#ifndef CARQUET_READER_INTERNAL_H
#define CARQUET_READER_INTERNAL_H

#include <carquet/carquet.h>
#include "thrift/parquet_types.h"
#include "core/arena.h"
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Memory Mapping Types
 * ============================================================================
 */

/**
 * Indicates whether data is owned (malloc'd) or a view (mmap pointer).
 */
typedef enum carquet_data_ownership {
    CARQUET_DATA_OWNED = 0,    /* Data is malloc'd, caller must free */
    CARQUET_DATA_VIEW = 1,     /* Data is view into mmap, do NOT free */
} carquet_data_ownership_t;

/**
 * Platform-specific memory mapping handle.
 */
typedef struct carquet_mmap_info {
    uint8_t* data;
    size_t size;
#ifdef _WIN32
    HANDLE file_handle;
    HANDLE mapping_handle;
#else
    int fd;
#endif
    bool is_valid;
} carquet_mmap_info_t;

/* ============================================================================
 * Internal Schema Structure
 * ============================================================================
 */

struct carquet_schema {
    carquet_arena_t arena;
    parquet_schema_element_t* elements;
    int32_t* parent_indices;    /* Parent element index for each element (-1 for root) */
    int32_t num_elements;
    int32_t capacity;           /* Capacity of elements/leaf arrays */

    /* Computed fields */
    int32_t* leaf_indices;      /* Map leaf index -> schema element index */
    int32_t num_leaves;         /* Number of leaf columns */
    int16_t* max_def_levels;    /* Max definition level per leaf */
    int16_t* max_rep_levels;    /* Max repetition level per leaf */
};

/* ============================================================================
 * Internal Reader Structure
 * ============================================================================
 */

/**
 * Pre-buffered I/O cache for coalesced reads.
 */
typedef struct carquet_prebuffer {
    uint8_t* data;          /* Coalesced read buffer (kept across calls) */
    size_t capacity;        /* Allocated bytes behind data */
    int64_t file_offset;    /* Start offset in file (-1 = no cached range) */
    size_t size;            /* Valid bytes cached from file_offset */
    int32_t row_group;      /* Row group this cache is for (-1 = none) */
} carquet_prebuffer_t;

struct carquet_reader {
    FILE* file;
    bool owns_file;
    /* `file` is a regular disk file this reader opened itself: page reads use
     * positional I/O (carquet_file_pread) instead of fseek+fread, so column
     * readers of one file can load pages from several threads at once. */
    bool positional_io;

    /* Memory-mapped data */
    const uint8_t* mmap_data;
    size_t file_size;
    carquet_mmap_info_t* mmap_info;  /* Platform-specific mmap handle, NULL if not using mmap */

    /* Metadata */
    carquet_arena_t arena;
    parquet_file_metadata_t metadata;
    carquet_schema_t* schema;

    /* Options */
    carquet_reader_options_t options;

    /* Pre-buffered I/O cache */
    carquet_prebuffer_t prebuffer;

    /* State */
    bool is_open;
};

/* ============================================================================
 * Internal Column Reader Structure
 * ============================================================================
 */

/**
 * Node in the BYTE_ARRAY page retention list.
 * Flexible array member holds the actual page bytes.
 */
typedef struct carquet_retained_page {
    struct carquet_retained_page* next;
    size_t size;
    uint8_t data[];
} carquet_retained_page_t;

struct carquet_column_reader {
    carquet_reader_t* file_reader;
    int32_t row_group_index;
    int32_t column_index;

    /* Column metadata */
    const parquet_column_chunk_t* chunk;
    const parquet_column_metadata_t* col_meta;

    /* Schema info */
    int16_t max_def_level;
    int16_t max_rep_level;
    carquet_physical_type_t type;
    int32_t type_length;

    /* Reading state */
    int64_t values_remaining;
    int64_t data_start_offset;    /* Actual offset of first data page in file */
    int64_t current_page;

    /* Page data */
    uint8_t* page_buffer;
    size_t page_buffer_size;
    size_t page_buffer_capacity;

    /* Dictionary */
    bool has_dictionary;
    uint8_t* dictionary_data;
    size_t dictionary_size;
    int32_t dictionary_count;
    uint32_t* dictionary_offsets;  /* Offset cache for O(1) BYTE_ARRAY lookup */
    carquet_data_ownership_t dictionary_ownership; /* OWNED or VIEW */

    /* Retained page data for BYTE_ARRAY value pointers.
     * When a BYTE_ARRAY PLAIN page is decoded, the resulting
     * carquet_byte_array_t.data pointers reference bytes inside the
     * (decompressed or file) page buffer. A single batch read may span
     * multiple pages, so we must keep EVERY page buffer alive until the
     * batch is consumed. We accumulate them in a singly linked list that
     * is flushed on row-group reset and column-reader free. */
    struct carquet_retained_page* retained_pages;

    /* Current page state for partial reads */
    bool page_loaded;           /* Is a page currently loaded? */
    int32_t page_num_values;    /* Total values in current page */
    int32_t page_values_read;   /* Values already read from current page */
    int32_t dense_cursor_pos;   /* page_values_read at which dense_cursor_val is valid; -1 after a page change */
    int32_t dense_cursor_val;   /* present values before dense_cursor_pos (OPTIONAL columns) */
    int32_t page_header_size;   /* Size of current page header */
    int32_t page_compressed_size; /* Size of current page compressed data */
    uint8_t* decoded_values;    /* Buffer for decoded values from current page */
    int16_t* decoded_def_levels; /* Buffer for decoded definition levels
                                  * (NULL for carquet_column_def_as_bits()
                                  * columns, see decoded_def_bits) */
    uint8_t* decoded_def_bits;   /* Flat max_def == 1 columns: the page's
                                  * definition levels as a validity bitmap,
                                  * bit i (LSB first) = level i == 1. The 1-bit
                                  * level stream already is that bitmap, so it
                                  * is decoded straight into it; int16 levels
                                  * are only expanded for a caller asking for
                                  * them (carquet_read_next_page). */
    size_t def_bits_capacity;    /* Bytes behind decoded_def_bits */
    int16_t* decoded_rep_levels; /* Buffer for decoded repetition levels */
    size_t decoded_capacity;    /* Capacity of decoded buffers, in values */
    size_t decoded_value_size;  /* Per-value byte width the values buffer was
                                 * allocated with; realloc when it changes (e.g.
                                 * preserve_dictionary flips index/value width) */
    carquet_data_ownership_t decoded_ownership; /* OWNED or VIEW (mmap) */

    /* Reusable buffers to reduce allocations */
    uint32_t* indices_buffer;   /* Reusable buffer for dictionary indices */
    size_t indices_capacity;    /* Capacity of indices buffer */
    uint8_t* decompress_buffer; /* Reusable decompression buffer */
    size_t decompress_capacity; /* Capacity of decompression buffer */
    int16_t* scratch_def_levels; /* Reusable def-level scratch for batch reads
                                  * on a nullable column when the caller passes
                                  * def_levels == NULL; grown on demand, freed
                                  * with the reader (was a malloc/free per call) */
    size_t scratch_def_capacity; /* Capacity of scratch_def_levels, in levels */

    /* Dictionary preservation mode */
    bool preserve_dictionary;   /* If true, skip materialization, keep indices */

    /* Last failure seen by this column reader, for callers that read through an
     * API with no error out-parameter (the batch reader). Written only by the
     * thread that owns this column reader — the batch reader's parallel column
     * loop gives each thread its own reader, so this is not a shared write.
     * Cleared by carquet_column_read_batch_ex() on entry. */
    carquet_error_t last_error;
};

/* ============================================================================
 * Internal Functions
 * ============================================================================
 */

/**
 * Whether the column reader decodes its pages' definition levels into a
 * validity bitmap (decoded_def_bits) instead of int16 levels: flat columns
 * whose only levels are 0/1 definition levels (max_def == 1, max_rep == 0),
 * i.e. top-level nullable columns.
 */
static inline bool carquet_column_def_as_bits(const carquet_column_reader_t* reader) {
    return reader->max_def_level == 1 && reader->max_rep_level == 0;
}

/**
 * Bits [at, at + n) of an LSB-first bitmap as a word (bit 0 = bit `at`),
 * 1 <= n <= 64. Reads only the bytes covering the range.
 */
static inline uint64_t carquet_bitmap_load_bits(const uint8_t* bm, int64_t at, int n) {
    const uint8_t* p = bm + (at >> 3);
    unsigned s = (unsigned)(at & 7);
    int nbytes = (int)((s + (unsigned)n + 7) >> 3);   /* 1..9 */
    uint64_t w = 0;
    if (nbytes >= 8) {
        for (int b = 0; b < 8; b++) w |= (uint64_t)p[b] << (8 * b);
        w >>= s;
        if (nbytes == 9) w |= (uint64_t)p[8] << (64 - s);   /* s >= 1 here */
    } else {
        for (int b = 0; b < nbytes; b++) w |= (uint64_t)p[b] << (8 * b);
        w >>= s;
    }
    if (n < 64) w &= ((uint64_t)1 << n) - 1;
    return w;
}

/** Set bits of bm[start, start + n). */
int64_t carquet_bitmap_count_range(const uint8_t* bm, int64_t start, int64_t n);

/**
 * Copy bits src[src_at, src_at + n) to dst[dst_at, ...). The bits of dst
 * from dst_at on must be zero; only the bytes covering the range are
 * written and only the bytes covering the source range are read.
 */
void carquet_bitmap_copy_range(uint8_t* dst, int64_t dst_at,
                               const uint8_t* src, int64_t src_at, int64_t n);

/** out[i] = bit (start + i) of bm, as a 0/1 level, for i < n. */
void carquet_bitmap_to_levels(const uint8_t* bm, int64_t start, int64_t n, int16_t* out);

/**
 * Build schema structure from parsed metadata.
 */
carquet_schema_t* build_schema(
    carquet_arena_t* arena,
    const parquet_file_metadata_t* metadata,
    carquet_error_t* error);

/**
 * Validate a row group index against a reader's loaded metadata.
 */
bool carquet_reader_row_group_index_valid(
    const carquet_reader_t* reader,
    int32_t row_group_index);

/**
 * Re-point a column reader at another row group of the same file, keeping its
 * page/decode/scratch buffers. Drops the dictionary, retained pages and all
 * per-chunk reading state. No validation: callers pass a valid row group and
 * the column index the reader was created for (see carquet_column_reader_reset
 * for the checked public variant).
 */
void carquet_column_reader_rebind(
    carquet_column_reader_t* col_reader,
    carquet_reader_t* file_reader,
    int32_t row_group_index,
    int32_t column_index,
    bool preserve_dictionaries);

/**
 * Open file with memory mapping.
 * Returns mmap_info on success, NULL on failure (fallback to fread).
 */
carquet_mmap_info_t* carquet_mmap_open(const char* path, carquet_error_t* error);

/**
 * Close memory mapping and release resources.
 */
void carquet_mmap_close(carquet_mmap_info_t* mmap_info);

/**
 * Check if a page is eligible for zero-copy reading.
 * Requires: uncompressed, PLAIN encoding, fixed-size type.
 */
bool carquet_page_is_zero_copy_eligible(
    carquet_compression_t codec,
    carquet_encoding_t encoding,
    carquet_physical_type_t type);

/**
 * Ensure the current page is loaded and ready for reading.
 * Advances to the next page when the current one has been fully consumed.
 */
carquet_status_t carquet_column_ensure_page_loaded(
    carquet_column_reader_t* reader,
    carquet_error_t* error);

/**
 * Retain a copy of a page data buffer on the column reader's retention list.
 * Returns the pointer to the retained bytes (stable until the retention list
 * is flushed) on success, or NULL on out-of-memory.
 * BYTE_ARRAY PLAIN decoding uses this so that carquet_byte_array_t.data
 * pointers remain valid across page boundaries within a batch.
 */
uint8_t* carquet_column_retain_page(
    carquet_column_reader_t* reader,
    const uint8_t* src,
    size_t size);

/**
 * Free and clear the retained page list on the column reader.
 */
void carquet_column_clear_retained_pages(carquet_column_reader_t* reader);

/**
 * Ensure the column's dictionary page (if any) is loaded. Idempotent.
 * After this returns CARQUET_OK, reader->data_start_offset points to the
 * first data page even when the writer recorded a stale data_page_offset.
 */
carquet_status_t carquet_column_ensure_dictionary_loaded(
    carquet_column_reader_t* reader,
    carquet_error_t* error);

/**
 * Seek the column reader to a specific data page so that the next page
 * load decodes that page. Dictionary state is preserved. The value count
 * for skipped pages is recovered by walking page headers between the
 * current position and the target, so the caller does not need to know
 * it; the `values_before_page` parameter is informational only.
 *
 * @param reader              Column reader
 * @param page_file_offset    Absolute file offset of the target data page
 * @param values_before_page  Reserved; currently unused (pass 0).
 */
carquet_status_t carquet_column_reader_seek_to_data_page(
    carquet_column_reader_t* reader,
    int64_t page_file_offset,
    int64_t values_before_page,
    carquet_error_t* error);

/**
 * @brief Pick the offset of the first data page of a dictionary-encoded chunk.
 *
 * @param data_page_offset ColumnMetaData.data_page_offset as recorded by the writer
 * @param dict_end         dictionary_page_offset + dict header + dict payload
 * @return @p data_page_offset unless it is provably wrong (at or before
 *         @p dict_end), in which case @p dict_end.
 */
int64_t carquet_resolve_data_start_offset(int64_t data_page_offset,
                                          int64_t dict_end);

/* ============================================================================
 * Page Decompression (shared between page_reader and batch_reader)
 * ============================================================================
 */

carquet_status_t carquet_decompress_page(
    carquet_compression_t codec,
    const uint8_t* compressed,
    size_t compressed_size,
    uint8_t* decompressed,
    size_t decompressed_capacity,
    size_t* decompressed_size);

/** @brief Size of the buffer to pass to carquet_decompress_page_diag(). */
#define CARQUET_DECOMPRESS_DIAG_SIZE 192

/**
 * @brief Decompress a page, additionally writing a codec-specific description
 *        of the failure into @p diag.
 *
 * On failure the caller can fold @p diag into its error message so that a
 * single failing run on a multi-gigabyte file identifies the exact check that
 * rejected the block. @p diag is always NUL-terminated; it is set to the empty
 * string when the codec has no detail to report.
 */
carquet_status_t carquet_decompress_page_diag(
    carquet_compression_t codec,
    const uint8_t* compressed,
    size_t compressed_size,
    uint8_t* decompressed,
    size_t decompressed_capacity,
    size_t* decompressed_size,
    char* diag,
    size_t diag_size);

/**
 * File leaf-column index behind each column of a batch produced by a batch
 * reader (length carquet_row_batch_num_columns()), or NULL for a batch that
 * did not come from carquet_batch_reader_next(). Lets the Arrow exporter map
 * a projected (or reordered) batch back onto the file schema.
 */
const int32_t* carquet_row_batch_projection(const carquet_row_batch_t* batch);

/* Below this many values (rows x columns) a batch is exported to Arrow on the
 * calling thread, and builds without OpenMP also decode it there: handing the
 * columns to other threads costs more than the work. */
#define CARQUET_PARALLEL_MIN_VALUES 65536

/**
 * Worker pool the Arrow exporter may use to export a batch's columns in
 * parallel (builds without OpenMP), or NULL to export serially.
 */
struct carquet_worker_pool* carquet_row_batch_pool(const carquet_row_batch_t* batch);

/* Select and reorder the children of a top-level struct ArrowSchema in place:
 * children keep[0..n-1] (indices into the current children) become the new
 * children, every other child is released. Implemented in arrow_c_read.c. */
carquet_status_t carquet_arrow_schema_select_children(
    struct ArrowSchema* schema, const int32_t* keep, int32_t n);

/* Row predicates for the Arrow stream (arrow_filter.c). `top_elems[0..ntop)`
 * are the schema element indices of the batch's top-level fields, in batch
 * order; every clause column must sit under one of them, with no LIST / MAP
 * ancestor. */
typedef struct carquet_compiled_predicate carquet_compiled_predicate_t;
carquet_status_t carquet_predicate_compile(
    const carquet_reader_t* reader, const carquet_predicate_t* pred,
    const int32_t* top_elems, int32_t ntop,
    carquet_compiled_predicate_t** out, carquet_error_t* error);
void carquet_predicate_free(carquet_compiled_predicate_t* p);
/* Statistics say some row of `rg` may match (`pred` is the uncompiled form). */
bool carquet_predicate_row_group_might_match(const carquet_compiled_predicate_t* p,
                                             const carquet_predicate_t* pred,
                                             const carquet_reader_t* r, int32_t rg);
/* Page-filter clauses (owned by `p`): the single conjunction's prunable ones;
 * *count is 0 for several conjunctions. */
const carquet_filter_clause_t* carquet_predicate_page_clauses(const carquet_compiled_predicate_t* p,
                                                              int32_t* count);
/* Distinct leaf columns the predicate reads (leaves needs room for every clause). */
void carquet_predicate_columns(const carquet_compiled_predicate_t* p, int32_t* leaves, int32_t* n);
/* keep[i] = row i of `batch` (a top-level struct array) matches; *kept = count. */
carquet_status_t carquet_predicate_select(const carquet_compiled_predicate_t* p,
                                          const struct ArrowArray* batch,
                                          uint8_t* keep, int64_t* kept);
/* "ARROW:schema" restore (arrow_restore.c). carquet_restore_plan() sets the
 * footer key/value metadata on `schema` (a top-level struct of exported
 * fields) and, with `restore`, rewrites its types to the stored ones; *plan
 * is NULL when the batches need no conversion. */
typedef struct carquet_restore_plan carquet_restore_plan_t;
carquet_status_t carquet_arrow_stored_schema(const char* b64, struct ArrowSchema* out);
carquet_status_t carquet_restore_plan(const carquet_reader_t* reader, bool restore,
                                      struct ArrowSchema* schema, carquet_restore_plan_t** plan);
carquet_status_t carquet_restore_apply(const carquet_restore_plan_t* plan, struct ArrowArray* batch);
void carquet_restore_free(carquet_restore_plan_t* plan);

/* An owned array of `len` nulls of `schema`'s type (arrow_filter.c). */
carquet_status_t carquet_arrow_null_array(const struct ArrowSchema* schema, int64_t len, struct ArrowArray* out);

/* New owned array with rows `rows[0..n)` of `in` (logical indices, ascending
 * not required). Formats: struct, list / large list, map, null, bool,
 * fixed-width primitives, utf8 / binary (+ large). */
carquet_status_t carquet_arrow_take(const struct ArrowArray* in, const struct ArrowSchema* schema,
                                    const int64_t* rows, int64_t n, struct ArrowArray* out);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_READER_INTERNAL_H */
