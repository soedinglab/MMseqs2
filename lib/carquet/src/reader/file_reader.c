/**
 * @file file_reader.c
 * @brief Parquet file reader implementation
 */

#include "core/allocator.h"
#include "core/compat.h"
#include <carquet/carquet.h>
#include "reader_internal.h"
#include "thrift/parquet_types.h"
#include "thrift/thrift_decode.h"
#include "core/arena.h"
#include "core/buffer.h"
#include "core/endian.h"
#include "encoding/plain.h"
#include "encoding/rle.h"
#include "arrow_schema_read.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

/* External functions from metadata modules */
extern carquet_status_t carquet_bloom_filter_read(carquet_bloom_filter_t** filter_out,
                                                   const uint8_t* data, size_t data_size);
extern carquet_column_index_t* carquet_column_index_parse(const uint8_t* data, size_t size);
extern carquet_offset_index_t* carquet_offset_index_parse(const uint8_t* data, size_t size);

/* ============================================================================
 * Constants
 * ============================================================================
 */

#define PARQUET_MAGIC "PAR1"
/* Trailing magic of a file whose footer is encrypted (Parquet modular
 * encryption). Recognised only to produce a clear error. */
#define PARQUET_ENCRYPTED_MAGIC "PARE"
#define PARQUET_MAGIC_LEN 4
#define PARQUET_FOOTER_SIZE_LEN 4

/* ============================================================================
 * Schema Building
 * ============================================================================
 */

static int32_t count_leaves(const parquet_schema_element_t* elements, int32_t count) {
    int32_t leaves = 0;
    for (int32_t i = 0; i < count; i++) {
        if (elements[i].num_children == 0) {
            leaves++;
        }
    }
    return leaves;
}

/**
 * Recursive schema traversal context for computing definition/repetition levels.
 */
typedef struct {
    const parquet_schema_element_t* elements;
    int32_t num_elements;
    int16_t* max_def;
    int16_t* max_rep;
    int32_t* leaf_indices;
    int32_t* parent_indices;
    int32_t leaf_idx;
    bool depth_exceeded;
} schema_traverse_ctx_t;

/* Cap schema nesting depth to keep the recursive traversal below the call-stack
 * limit even on threads with small stacks. A crafted file could otherwise nest
 * up to CARQUET_MAX_SCHEMA_ELEMENTS groups deep and overflow the stack. Real
 * Parquet schemas are only a handful of levels deep. */
#define CARQUET_MAX_SCHEMA_DEPTH 100

/**
 * Recursively traverse schema tree and compute definition/repetition levels.
 *
 * @param ctx Traversal context
 * @param element_idx Current element index in flat array
 * @param def_level Current definition level from ancestors
 * @param rep_level Current repetition level from ancestors
 * @return Next element index to process (after this subtree)
 */
static int32_t traverse_schema_recursive(
    schema_traverse_ctx_t* ctx,
    int32_t element_idx,
    int32_t parent_idx,
    int16_t def_level,
    int16_t rep_level,
    int32_t depth) {

    if (element_idx >= ctx->num_elements) {
        return element_idx;
    }

    if (depth > CARQUET_MAX_SCHEMA_DEPTH) {
        /* Refuse to recurse further; compute_levels reports this as an error. */
        ctx->depth_exceeded = true;
        return element_idx;
    }

    const parquet_schema_element_t* elem = &ctx->elements[element_idx];

    /* Record parent index */
    if (ctx->parent_indices) {
        ctx->parent_indices[element_idx] = parent_idx;
    }

    /* Calculate level contribution from this node's repetition type */
    int16_t this_def = def_level;
    int16_t this_rep = rep_level;

    if (elem->has_repetition) {
        switch (elem->repetition_type) {
            case CARQUET_REPETITION_OPTIONAL:
                /* Optional fields add 1 to definition level */
                this_def++;
                break;
            case CARQUET_REPETITION_REPEATED:
                /* Repeated fields add 1 to both definition and repetition levels */
                this_def++;
                this_rep++;
                break;
            case CARQUET_REPETITION_REQUIRED:
            default:
                /* Required fields don't add to levels */
                break;
        }
    }

    if (elem->num_children == 0) {
        /* Leaf node - record the accumulated levels */
        ctx->max_def[ctx->leaf_idx] = this_def;
        ctx->max_rep[ctx->leaf_idx] = this_rep;
        ctx->leaf_indices[ctx->leaf_idx] = element_idx;
        ctx->leaf_idx++;
        return element_idx + 1;
    }

    /* Group node - recursively process children */
    int32_t next_idx = element_idx + 1;
    for (int32_t child = 0; child < elem->num_children; child++) {
        /* Stop once the flat element array is exhausted. Without this, a crafted
         * footer declaring num_children up to INT32_MAX would spin billions of
         * no-op recursive calls (each returns immediately via the guard at the
         * top) — a CPU denial-of-service on an otherwise tiny file. */
        if (next_idx >= ctx->num_elements) {
            break;
        }
        next_idx = traverse_schema_recursive(ctx, next_idx, element_idx,
                                             this_def, this_rep, depth + 1);
    }

    return next_idx;
}

/**
 * Compute definition and repetition levels for all leaf columns.
 *
 * Parquet stores schema as a flat array in depth-first order. This function
 * recursively traverses the schema tree to compute the maximum definition
 * and repetition levels for each leaf column.
 *
 * Definition level: Number of optional/repeated ancestors + 1 if self is optional/repeated
 * Repetition level: Number of repeated ancestors + 1 if self is repeated
 *
 * Example schema:
 *   schema (root, required)
 *   ├── a (optional, int32)        -> def=1, rep=0
 *   ├── b (optional, group)
 *   │   ├── c (required, int32)    -> def=1, rep=0  (from parent b)
 *   │   └── d (optional, int32)    -> def=2, rep=0  (from b + self)
 *   └── e (repeated, group)
 *       ├── f (required, int32)    -> def=1, rep=1  (from parent e)
 *       └── g (optional, int32)    -> def=2, rep=1  (from e + self)
 */
static bool compute_levels(
    const parquet_schema_element_t* elements,
    int32_t num_elements,
    int16_t* max_def,
    int16_t* max_rep,
    int32_t* leaf_indices,
    int32_t* parent_indices) {

    if (num_elements <= 1) {
        return true;  /* Empty or root-only schema */
    }

    if (parent_indices) {
        parent_indices[0] = -1;  /* Root has no parent */
    }

    schema_traverse_ctx_t ctx = {
        .elements = elements,
        .num_elements = num_elements,
        .max_def = max_def,
        .max_rep = max_rep,
        .leaf_indices = leaf_indices,
        .parent_indices = parent_indices,
        .leaf_idx = 0,
        .depth_exceeded = false
    };

    /* Start traversal from root (index 0) with zero levels.
     * Root is required by definition, so it doesn't contribute to levels.
     * We process its children starting at index 1. */
    const parquet_schema_element_t* root = &elements[0];
    int32_t next_idx = 1;
    for (int32_t child = 0; child < root->num_children; child++) {
        /* Stop once the flat element array is exhausted. A crafted footer
         * declaring root->num_children up to INT32_MAX would otherwise spin
         * billions of no-op recursive calls (each returns immediately via the
         * guard at the top of traverse_schema_recursive) — a CPU denial-of-
         * service on a tiny file. Mirrors the guard in the recursive inner
         * loop above. */
        if (next_idx >= num_elements) {
            break;
        }
        next_idx = traverse_schema_recursive(&ctx, next_idx, 0, 0, 0, 1);
    }

    return !ctx.depth_exceeded;
}

carquet_schema_t* build_schema(
    carquet_arena_t* arena,
    const parquet_file_metadata_t* metadata,
    carquet_error_t* error) {

    carquet_schema_t* schema = carquet_arena_calloc(arena, 1, sizeof(carquet_schema_t));
    if (!schema) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema");
        return NULL;
    }

    schema->elements = metadata->schema;
    schema->num_elements = metadata->num_schema_elements;
    schema->capacity = metadata->num_schema_elements;  /* Fixed size from file */
    schema->num_leaves = count_leaves(metadata->schema, metadata->num_schema_elements);

    schema->parent_indices = carquet_arena_calloc(arena, schema->num_elements, sizeof(int32_t));
    schema->leaf_indices = carquet_arena_calloc(arena, schema->num_leaves, sizeof(int32_t));
    schema->max_def_levels = carquet_arena_calloc(arena, schema->num_leaves, sizeof(int16_t));
    schema->max_rep_levels = carquet_arena_calloc(arena, schema->num_leaves, sizeof(int16_t));

    if (!schema->parent_indices || !schema->leaf_indices ||
        !schema->max_def_levels || !schema->max_rep_levels) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema arrays");
        return NULL;
    }

    if (!compute_levels(schema->elements, schema->num_elements,
                        schema->max_def_levels, schema->max_rep_levels,
                        schema->leaf_indices, schema->parent_indices)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_SCHEMA,
                          "Schema nesting exceeds maximum depth of %d",
                          CARQUET_MAX_SCHEMA_DEPTH);
        return NULL;
    }

    /* Recover Arrow per-field custom_metadata (variable labels/descriptions)
     * from the "ARROW:schema" footer blob, if present. Best-effort: malformed
     * blobs are ignored and never fail the open. */
    for (int32_t i = 0; i < metadata->num_key_value; i++) {
        const parquet_key_value_t* kv = &metadata->key_value_metadata[i];
        if (kv->key && kv->value && strcmp(kv->key, "ARROW:schema") == 0) {
            carquet_apply_arrow_field_metadata(kv->value, schema->elements,
                                               schema->num_elements,
                                               schema->parent_indices, arena);
            break;
        }
    }

    return schema;
}

/* ============================================================================
 * File Reader Implementation
 * ============================================================================
 */

void carquet_reader_options_init(carquet_reader_options_t* options) {
    /* Parameter is nonnull per API contract */
    memset(options, 0, sizeof(*options));
    options->use_mmap = false;
    options->verify_checksums = true;
    options->buffer_size = 64 * 1024;
    options->num_threads = 0;
}

/**
 * Speculative footer read: read up to 64KB from end of file in a single I/O
 * call. Most Parquet footers fit within this, eliminating the second seek+read.
 * Falls back to a targeted read if the footer is larger than the initial read.
 */
#define CARQUET_FOOTER_SPECULATIVE_SIZE (64 * 1024)

static carquet_status_t read_footer(carquet_reader_t* reader, carquet_error_t* error) {
    /* Seek to end to get file size (64-bit aware for files >2 GiB on Windows) */
    if (carquet_fseek64(reader->file, 0, SEEK_END) != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK, "Failed to seek to end");
        return CARQUET_ERROR_FILE_SEEK;
    }

    int64_t file_size = carquet_ftell64(reader->file);
    if (file_size < 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to get file size");
        return CARQUET_ERROR_FILE_READ;
    }
    reader->file_size = (size_t)file_size;

    /* Check minimum size */
    if (reader->file_size < PARQUET_MAGIC_LEN * 2 + PARQUET_FOOTER_SIZE_LEN) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "File too small");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    /* Speculative read: grab min(file_size, 64KB) from end of file in one I/O.
     * This captures both the 8-byte tail (magic + footer length) and, for most
     * files, the entire Thrift-encoded footer in a single fread call. */
    size_t spec_size = reader->file_size < CARQUET_FOOTER_SPECULATIVE_SIZE
                     ? reader->file_size : CARQUET_FOOTER_SPECULATIVE_SIZE;
    uint8_t* spec_buf = carquet_mem_malloc(spec_size);
    if (!spec_buf) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate footer buffer");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    int64_t spec_offset = (int64_t)(reader->file_size - spec_size);
    if (carquet_fseek64(reader->file, spec_offset, SEEK_SET) != 0) {
        carquet_mem_free(spec_buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK, "Failed to seek to footer");
        return CARQUET_ERROR_FILE_SEEK;
    }

    if (fread(spec_buf, 1, spec_size, reader->file) != spec_size) {
        carquet_mem_free(spec_buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to read footer tail");
        return CARQUET_ERROR_FILE_READ;
    }

    /* Verify trailing magic (last 4 bytes of file) */
    if (memcmp(spec_buf + spec_size - 4, PARQUET_MAGIC, PARQUET_MAGIC_LEN) != 0) {
        bool encrypted = memcmp(spec_buf + spec_size - 4, PARQUET_ENCRYPTED_MAGIC,
                                PARQUET_MAGIC_LEN) == 0;
        carquet_mem_free(spec_buf);
        if (encrypted) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "File uses Parquet modular encryption (encrypted footer, magic "
                "\"PARE\"), which carquet does not support");
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC, "Invalid trailing magic");
        return CARQUET_ERROR_INVALID_MAGIC;
    }

    /* Get footer size (4 bytes before trailing magic) */
    uint32_t footer_size = carquet_read_u32_le(spec_buf + spec_size - 8);
    if (footer_size > reader->file_size - 8) {
        carquet_mem_free(spec_buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "Footer size too large");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    const uint8_t* footer_data;
    uint8_t* fallback_buf = NULL;

    if ((size_t)footer_size + 8 <= spec_size) {
        /* Fast path: footer fits within the speculative read - no second I/O.
         * The cast to size_t is required: `footer_size + 8` in 32-bit unsigned
         * arithmetic wraps for footer_size >= 0xFFFFFFF8 (reachable on files
         * larger than 4GB, which pass the range check above), which would take
         * this fast path and underflow the pointer computation below. */
        footer_data = spec_buf + spec_size - 8 - footer_size;
    } else {
        /* Slow path: footer is larger than speculative buffer, need second read */
        fallback_buf = carquet_mem_malloc(footer_size);
        if (!fallback_buf) {
            carquet_mem_free(spec_buf);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate footer buffer");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }

        int64_t footer_offset = (int64_t)(reader->file_size - 8 - footer_size);
        if (carquet_fseek64(reader->file, footer_offset, SEEK_SET) != 0) {
            carquet_mem_free(fallback_buf);
            carquet_mem_free(spec_buf);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK, "Failed to seek to footer data");
            return CARQUET_ERROR_FILE_SEEK;
        }

        if (fread(fallback_buf, 1, footer_size, reader->file) != footer_size) {
            carquet_mem_free(fallback_buf);
            carquet_mem_free(spec_buf);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to read footer data");
            return CARQUET_ERROR_FILE_READ;
        }

        footer_data = fallback_buf;
    }

    /* Parse metadata */
    carquet_status_t status = parquet_parse_file_metadata(
        footer_data, footer_size, &reader->arena, &reader->metadata, error);

    carquet_mem_free(fallback_buf);
    carquet_mem_free(spec_buf);

    if (status != CARQUET_OK) {
        return status;
    }

    /* Build schema */
    reader->schema = build_schema(&reader->arena, &reader->metadata, error);
    if (!reader->schema) {
        return CARQUET_ERROR_INVALID_SCHEMA;
    }

    return CARQUET_OK;
}

/**
 * Read footer from memory-mapped data.
 */
static carquet_status_t read_footer_mmap(carquet_reader_t* reader, carquet_error_t* error) {
    const uint8_t* data = reader->mmap_data;
    size_t file_size = reader->file_size;

    /* Check minimum size */
    if (file_size < PARQUET_MAGIC_LEN * 2 + PARQUET_FOOTER_SIZE_LEN) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "File too small");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    /* Verify magic bytes at start and end */
    if (memcmp(data, PARQUET_MAGIC, PARQUET_MAGIC_LEN) != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC, "Invalid header magic");
        return CARQUET_ERROR_INVALID_MAGIC;
    }

    const uint8_t* end = data + file_size;
    if (memcmp(end - 4, PARQUET_MAGIC, PARQUET_MAGIC_LEN) != 0) {
        if (memcmp(end - 4, PARQUET_ENCRYPTED_MAGIC, PARQUET_MAGIC_LEN) == 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "File uses Parquet modular encryption (encrypted footer, magic "
                "\"PARE\"), which carquet does not support");
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_MAGIC, "Invalid trailing magic");
        return CARQUET_ERROR_INVALID_MAGIC;
    }

    /* Get footer size */
    uint32_t footer_size = carquet_read_u32_le(end - 8);
    if (footer_size > file_size - 8) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_FOOTER, "Footer size too large");
        return CARQUET_ERROR_INVALID_FOOTER;
    }

    /* Parse metadata directly from mmap (zero-copy) */
    const uint8_t* footer_data = end - 8 - footer_size;
    carquet_status_t status = parquet_parse_file_metadata(
        footer_data, footer_size, &reader->arena, &reader->metadata, error);

    if (status != CARQUET_OK) {
        return status;
    }

    /* Build schema */
    reader->schema = build_schema(&reader->arena, &reader->metadata, error);
    if (!reader->schema) {
        return CARQUET_ERROR_INVALID_SCHEMA;
    }

    return CARQUET_OK;
}

carquet_reader_t* carquet_reader_open(
    const char* path,
    const carquet_reader_options_t* options,
    carquet_error_t* error) {

    carquet_reader_t* reader = carquet_mem_calloc(1, sizeof(carquet_reader_t));
    if (!reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate reader");
        return NULL;
    }

    if (options) {
        reader->options = *options;
    } else {
        carquet_reader_options_init(&reader->options);
    }

    /* Initialize arena */
    if (carquet_arena_init(&reader->arena) != CARQUET_OK) {
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to initialize arena");
        return NULL;
    }

    reader->prebuffer.row_group = -1;

    carquet_status_t status;

    /* Try mmap if requested */
    if (reader->options.use_mmap) {
        /* Use a scratch error for the mmap attempt: if mmap fails but the fread
         * fallback below succeeds, the caller's `error` must not be left holding
         * the (recovered-from) mmap failure. */
        carquet_error_t mmap_err = {0};
        carquet_mmap_info_t* mmap_info = carquet_mmap_open(path, &mmap_err);
        if (mmap_info) {
            reader->mmap_info = mmap_info;
            reader->mmap_data = mmap_info->data;
            reader->file_size = mmap_info->size;
            reader->owns_file = false;  /* mmap handles cleanup */

            /* Parse footer from mmap */
            status = read_footer_mmap(reader, error);
            if (status != CARQUET_OK) {
                carquet_mmap_close(reader->mmap_info);
                carquet_arena_destroy(&reader->arena);
                carquet_mem_free(reader);
                return NULL;
            }

            reader->is_open = true;
            return reader;
        }
        /* mmap failed, fall through to fread path */
    }

    /* Standard fread path */
    FILE* file = fopen(path, "rb");
    if (!file) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_OPEN, "Failed to open file: %s", path);
        return NULL;
    }

    reader->file = file;
    reader->owns_file = true;

    /* Read and parse footer */
    status = read_footer(reader, error);
    if (status != CARQUET_OK) {
        carquet_arena_destroy(&reader->arena);
        fclose(file);
        carquet_mem_free(reader);
        return NULL;
    }

    reader->is_open = true;
    return reader;
}

carquet_reader_t* carquet_reader_open_file(
    FILE* file,
    const carquet_reader_options_t* options,
    carquet_error_t* error) {

    /* file is nonnull per API contract (matches carquet_reader_open) */
    carquet_reader_t* reader = carquet_mem_calloc(1, sizeof(carquet_reader_t));
    if (!reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate reader");
        return NULL;
    }

    if (options) {
        reader->options = *options;
    } else {
        carquet_reader_options_init(&reader->options);
    }
    /* mmap is meaningless for a caller-provided stream. */
    reader->options.use_mmap = false;

    if (carquet_arena_init(&reader->arena) != CARQUET_OK) {
        carquet_mem_free(reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to initialize arena");
        return NULL;
    }

    reader->prebuffer.row_group = -1;
    reader->file = file;
    reader->owns_file = false;  /* Caller retains ownership of the FILE handle */

    carquet_status_t status = read_footer(reader, error);
    if (status != CARQUET_OK) {
        carquet_arena_destroy(&reader->arena);
        carquet_mem_free(reader);
        return NULL;
    }

    reader->is_open = true;
    return reader;
}

carquet_status_t carquet_get_file_info(
    const char* path,
    carquet_file_info_t* info,
    carquet_error_t* error) {

    /* path and info are nonnull per API contract */
    memset(info, 0, sizeof(*info));

    carquet_reader_t* reader = carquet_reader_open(path, NULL, error);
    if (!reader) {
        return error ? error->code : CARQUET_ERROR_FILE_OPEN;
    }

    info->file_size = (int64_t)reader->file_size;
    info->num_rows = reader->metadata.num_rows;
    info->num_row_groups = reader->metadata.num_row_groups;
    info->num_columns = reader->schema ? reader->schema->num_leaves : 0;
    info->version = reader->metadata.version;

    /* Copy created_by into the caller-owned inline buffer (truncating if
     * needed) so it stays valid after the reader is closed. */
    const char* cb = reader->metadata.created_by;
    if (cb) {
        size_t n = strlen(cb);
        if (n >= sizeof(info->created_by)) {
            n = sizeof(info->created_by) - 1;
        }
        memcpy(info->created_by, cb, n);
        info->created_by[n] = '\0';
    } else {
        info->created_by[0] = '\0';
    }

    carquet_reader_close(reader);
    return CARQUET_OK;
}

carquet_status_t carquet_validate_file(
    const char* path,
    carquet_error_t* error) {

    /* path is nonnull per API contract.
     *
     * Stage 1: carquet_reader_open performs structural validation - magic
     * bytes, footer size, Thrift metadata parse, and schema construction.
     *
     * Stage 2: stream every row group / column / page with verify_checksums
     * enabled. This forces each page to be read, CRC32-verified (where a CRC
     * is present), decompressed, and decoded - surfacing any corruption that
     * a footer-only check would miss. */
    carquet_error_t local = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.verify_checksums = true;

    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &local);
    if (!reader) {
        if (error) *error = local;
        return local.code;
    }

    carquet_status_t result = CARQUET_OK;

    /* Files with no columns or no rows have no pages to scan; structural
     * validation alone is conclusive for them. */
    if (reader->schema && reader->schema->num_leaves > 0 &&
        reader->metadata.num_rows > 0) {

        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);  /* NULL projection = all columns */

        carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &local);
        if (!br) {
            if (error) *error = local;
            result = local.code;
        } else {
            for (;;) {
                carquet_row_batch_t* batch = NULL;
                carquet_status_t st = carquet_batch_reader_next(br, &batch);
                if (st == CARQUET_ERROR_END_OF_DATA) {
                    break;
                }
                if (st != CARQUET_OK) {
                    CARQUET_SET_ERROR(&local, st,
                        "Page validation failed (checksum/decode error)");
                    if (error) *error = local;
                    result = st;
                    if (batch) carquet_row_batch_free(batch);
                    break;
                }
                if (batch) carquet_row_batch_free(batch);
            }
            carquet_batch_reader_free(br);
        }
    }

    carquet_reader_close(reader);
    return result;
}

void carquet_reader_close(carquet_reader_t* reader) {
    if (!reader) return;

    /* Release prebuffer cache */
    carquet_reader_release_prebuffer(reader);

    /* Close mmap if active */
    if (reader->mmap_info) {
        carquet_mmap_close(reader->mmap_info);
        reader->mmap_info = NULL;
        reader->mmap_data = NULL;
    }

    if (reader->owns_file && reader->file) {
        fclose(reader->file);
    }

    carquet_arena_destroy(&reader->arena);
    carquet_mem_free(reader);
}

const carquet_schema_t* carquet_reader_schema(const carquet_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->schema;
}

int64_t carquet_reader_num_rows(const carquet_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->metadata.num_rows;
}

int32_t carquet_reader_num_row_groups(const carquet_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->metadata.num_row_groups;
}

int32_t carquet_reader_num_columns(const carquet_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->schema->num_leaves;
}

carquet_status_t carquet_reader_row_group_metadata(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    carquet_row_group_metadata_t* metadata) {

    /* reader and metadata are nonnull per API contract */
    if (!carquet_reader_row_group_index_valid(reader, row_group_index)) {
        return CARQUET_ERROR_ROW_GROUP_NOT_FOUND;
    }

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];
    metadata->num_rows = rg->num_rows;
    metadata->total_byte_size = rg->total_byte_size;
    metadata->total_compressed_size = rg->has_total_compressed_size ?
        rg->total_compressed_size : rg->total_byte_size;

    return CARQUET_OK;
}

/* ============================================================================
 * I/O Coalescing (Pre-buffering)
 * ============================================================================
 */

/** Maximum gap between column ranges to coalesce (1 MB) */
#define CARQUET_COALESCE_HOLE_LIMIT (1024 * 1024)

carquet_status_t carquet_reader_prebuffer(
    carquet_reader_t* reader,
    int32_t row_group_index,
    const int32_t* column_indices,
    int32_t num_columns,
    carquet_error_t* error) {

    /* No-op for mmap readers (OS handles page coalescing) */
    if (reader->mmap_data) {
        return CARQUET_OK;
    }

    if (!reader->file) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "Reader has no file handle");
        return CARQUET_ERROR_INVALID_STATE;
    }

    if (!carquet_reader_row_group_index_valid(reader, row_group_index)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_ROW_GROUP_NOT_FOUND,
            "Row group %d not found", row_group_index);
        return CARQUET_ERROR_ROW_GROUP_NOT_FOUND;
    }

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];

    /* Determine which columns to pre-buffer */
    int32_t total_cols = rg->num_columns;
    bool all_columns = (column_indices == NULL || num_columns <= 0);
    int32_t cols_count = all_columns ? total_cols : num_columns;

    if (cols_count <= 0) {
        return CARQUET_OK;
    }

    /* Find the min and max byte offsets across all requested columns */
    int64_t min_offset = INT64_MAX;
    int64_t max_end = 0;

    for (int32_t i = 0; i < cols_count; i++) {
        int32_t ci = all_columns ? i : column_indices[i];
        if (ci < 0 || ci >= total_cols) continue;

        const parquet_column_chunk_t* chunk = &rg->columns[ci];
        if (!chunk->has_metadata) continue;

        const parquet_column_metadata_t* meta = &chunk->metadata;
        int64_t col_start = meta->data_page_offset;

        /* Include dictionary page if present */
        if (meta->has_dictionary_page_offset &&
            meta->dictionary_page_offset < col_start) {
            col_start = meta->dictionary_page_offset;
        }

        if (col_start < 0 || meta->total_compressed_size < 0 ||
            col_start > INT64_MAX - meta->total_compressed_size) {
            continue;
        }

        int64_t col_end = col_start + meta->total_compressed_size;
        /* total_compressed_size is untrusted: a chunk claiming more bytes
         * than the file holds must not drive a multi-GB allocation. */
        if (col_start >= (int64_t)reader->file_size) continue;
        if (col_end > (int64_t)reader->file_size) col_end = (int64_t)reader->file_size;

        if (col_start < min_offset) min_offset = col_start;
        if (col_end > max_end) max_end = col_end;
    }

    if (min_offset >= max_end || min_offset == INT64_MAX) {
        return CARQUET_OK;
    }

    size_t total_size = (size_t)(max_end - min_offset);

    /* Release previous prebuffer if any */
    carquet_reader_release_prebuffer(reader);

    /* Allocate and read the coalesced range */
    uint8_t* buf = carquet_mem_malloc(total_size);
    if (!buf) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
            /* %llu, not %zu: the MinGW/msvcrt printf has no z length modifier. */
            "Failed to allocate prebuffer (%llu bytes)",
            (unsigned long long)total_size);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    if (carquet_fseek64(reader->file, (int64_t)min_offset, SEEK_SET) != 0) {
        carquet_mem_free(buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK, "Failed to seek for prebuffer");
        return CARQUET_ERROR_FILE_SEEK;
    }

    if (fread(buf, 1, total_size, reader->file) != total_size) {
        carquet_mem_free(buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ, "Failed to read prebuffer data");
        return CARQUET_ERROR_FILE_READ;
    }

    reader->prebuffer.data = buf;
    reader->prebuffer.file_offset = min_offset;
    reader->prebuffer.size = total_size;
    reader->prebuffer.row_group = row_group_index;

    return CARQUET_OK;
}

void carquet_reader_release_prebuffer(carquet_reader_t* reader) {
    if (reader->prebuffer.data) {
        carquet_mem_free(reader->prebuffer.data);
        reader->prebuffer.data = NULL;
        reader->prebuffer.file_offset = 0;
        reader->prebuffer.size = 0;
        reader->prebuffer.row_group = -1;
    }
}

/* ============================================================================
 * Column Reader Implementation
 * ============================================================================
 */

carquet_column_reader_t* carquet_reader_get_column(
    carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_error_t* error) {

    /* reader is nonnull per API contract */
    if (!carquet_reader_row_group_index_valid(reader, row_group_index)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_ROW_GROUP_NOT_FOUND,
            "Row group %d not found", row_group_index);
        return NULL;
    }

    if (column_index < 0 || column_index >= reader->schema->num_leaves) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_COLUMN_NOT_FOUND,
            "Column %d not found", column_index);
        return NULL;
    }

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];

    if (column_index >= rg->num_columns) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_COLUMN_NOT_FOUND,
            "Column %d not in row group", column_index);
        return NULL;
    }

    carquet_column_reader_t* col_reader = carquet_mem_calloc(1, sizeof(carquet_column_reader_t));
    if (!col_reader) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Failed to allocate column reader");
        return NULL;
    }

    col_reader->file_reader = reader;
    col_reader->row_group_index = row_group_index;
    col_reader->column_index = column_index;
    col_reader->chunk = &rg->columns[column_index];

    if (col_reader->chunk->has_metadata) {
        col_reader->col_meta = &col_reader->chunk->metadata;
    } else {
        /* Metadata might be in separate file - not supported yet */
        carquet_mem_free(col_reader);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "External column metadata not supported");
        return NULL;
    }

    /* Get schema info */
    int32_t schema_idx = reader->schema->leaf_indices[column_index];
    const parquet_schema_element_t* schema_elem = &reader->schema->elements[schema_idx];

    /* Guard against malformed metadata whose column-chunk physical type
     * disagrees with the schema. The value width is derived from two different
     * sources on two different code paths: the batch reader sizes its output
     * buffer from the schema element type, while the page reader writes using
     * this column reader's type (taken from the chunk metadata below). If the
     * two disagree, the page decode writes past the batch buffer — a
     * heap-buffer-overflow reachable from untrusted input (e.g. schema INT32,
     * 4 B/value, vs chunk BYTE_ARRAY, sizeof(carquet_byte_array_t)/value).
     * Both sides must agree on the type, so reject the file when they do not.
     * The fallback mirrors the batch reader's handling of a typeless element. */
    carquet_physical_type_t schema_type =
        schema_elem->has_type ? schema_elem->type : CARQUET_PHYSICAL_BYTE_ARRAY;
    if (col_reader->col_meta->type != schema_type) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Column %d physical type in chunk metadata (%d) does not match "
            "schema type (%d)", column_index,
            (int)col_reader->col_meta->type, (int)schema_type);
        carquet_mem_free(col_reader);
        return NULL;
    }

    col_reader->max_def_level = reader->schema->max_def_levels[column_index];
    col_reader->max_rep_level = reader->schema->max_rep_levels[column_index];
    col_reader->type = col_reader->col_meta->type;
    col_reader->type_length = schema_elem->type_length;

    col_reader->values_remaining = col_reader->col_meta->num_values;
    col_reader->data_start_offset = col_reader->col_meta->data_page_offset;

    return col_reader;
}

void carquet_column_reader_free(carquet_column_reader_t* reader) {
    if (!reader) return;

    carquet_mem_free(reader->page_buffer);
    carquet_column_clear_retained_pages(reader);
    if (reader->dictionary_ownership == CARQUET_DATA_OWNED) {
        carquet_mem_free(reader->dictionary_data);
    }
    carquet_mem_free(reader->dictionary_offsets);

    /* Only free decoded_values if we own the memory (not a mmap view) */
    if (reader->decoded_ownership == CARQUET_DATA_OWNED) {
        carquet_mem_free(reader->decoded_values);
    }

    /* Levels are always owned (decoded from RLE) */
    carquet_mem_free(reader->decoded_def_levels);
    carquet_mem_free(reader->decoded_rep_levels);
    carquet_mem_free(reader->indices_buffer);
    carquet_mem_free(reader->decompress_buffer);
    carquet_mem_free(reader->scratch_def_levels);
    carquet_mem_free(reader);
}

bool carquet_column_has_next(const carquet_column_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->values_remaining > 0;
}

int64_t carquet_column_remaining(const carquet_column_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->values_remaining;
}

/* ============================================================================
 * Memory Mapping API
 * ============================================================================
 */

bool carquet_reader_is_mmap(const carquet_reader_t* reader) {
    /* reader is nonnull per API contract */
    return reader->mmap_info != NULL && reader->mmap_info->is_valid;
}

bool carquet_reader_can_zero_copy(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index) {

    /* reader is nonnull per API contract */

    /* Must have mmap enabled */
    if (!reader->mmap_info || !reader->mmap_info->is_valid) {
        return false;
    }

    /* Validate indices */
    if (!carquet_reader_row_group_index_valid(reader, row_group_index)) {
        return false;
    }
    if (column_index < 0 || column_index >= reader->schema->num_leaves) {
        return false;
    }

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];
    if (column_index >= rg->num_columns) {
        return false;
    }

    const parquet_column_chunk_t* chunk = &rg->columns[column_index];
    if (!chunk->has_metadata) {
        return false;
    }

    const parquet_column_metadata_t* col_meta = &chunk->metadata;

    /* Must be uncompressed */
    if (col_meta->codec != CARQUET_COMPRESSION_UNCOMPRESSED) {
        return false;
    }

    /* Check if column has definition levels (nullable) */
    int16_t max_def = reader->schema->max_def_levels[column_index];
    if (max_def > 0) {
        return false;  /* Nullable columns need level decoding */
    }

    /* Check physical type - must be fixed-size */
    carquet_physical_type_t type = col_meta->type;
    switch (type) {
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_FLOAT:
        case CARQUET_PHYSICAL_DOUBLE:
        case CARQUET_PHYSICAL_INT96:
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            return true;

        case CARQUET_PHYSICAL_BOOLEAN:
        case CARQUET_PHYSICAL_BYTE_ARRAY:
        default:
            return false;
    }
}

/* ============================================================================
 * Library Version
 * ============================================================================
 */

const char* carquet_version(void) {
    return CARQUET_VERSION_STRING;
}

void carquet_version_components(int* major, int* minor, int* patch) {
    if (major) *major = CARQUET_VERSION_MAJOR;
    if (minor) *minor = CARQUET_VERSION_MINOR;
    if (patch) *patch = CARQUET_VERSION_PATCH;
}

/* ============================================================================
 * Internal Helper: Read bytes from file at a given offset
 * ============================================================================
 */

/**
 * Read `size` bytes from the file at `offset` into `out_buf`.
 * Handles both mmap (direct pointer) and fread (seek + read) paths.
 *
 * For mmap readers, `out_buf` is set to point directly into the mapped region
 * and `*allocated` is set to false. For fread readers, a buffer is malloc'd,
 * `out_buf` points to it, and `*allocated` is set to true. The caller must
 * free the buffer when `*allocated` is true.
 */
static carquet_status_t reader_read_bytes(
    carquet_reader_t* reader,
    int64_t offset,
    int32_t size,
    const uint8_t** out_buf,
    bool* allocated,
    carquet_error_t* error) {

    *allocated = false;

    if (size <= 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "Invalid read size");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Bounds check */
    if (offset < 0 || offset > INT64_MAX - (int64_t)size ||
        (size_t)offset > reader->file_size ||
        (size_t)size > reader->file_size - (size_t)offset) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ,
            "Read at offset %lld size %d exceeds file size %lld",
            (long long)offset, size, (long long)reader->file_size);
        return CARQUET_ERROR_FILE_READ;
    }

    /* mmap path: direct pointer into mapped region */
    if (reader->mmap_data) {
        *out_buf = reader->mmap_data + offset;
        return CARQUET_OK;
    }

    /* fread path: allocate buffer and read */
    if (!reader->file) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_STATE, "Reader has no file handle");
        return CARQUET_ERROR_INVALID_STATE;
    }

    uint8_t* buf = carquet_mem_malloc((size_t)size);
    if (!buf) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY,
            "Failed to allocate %d bytes for read", size);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    if (carquet_fseek64(reader->file, (int64_t)offset, SEEK_SET) != 0) {
        carquet_mem_free(buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_SEEK,
            "Failed to seek to offset %lld", (long long)offset);
        return CARQUET_ERROR_FILE_SEEK;
    }

    if (fread(buf, 1, (size_t)size, reader->file) != (size_t)size) {
        carquet_mem_free(buf);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_FILE_READ,
            "Failed to read %d bytes at offset %lld", size, (long long)offset);
        return CARQUET_ERROR_FILE_READ;
    }

    *out_buf = buf;
    *allocated = true;
    return CARQUET_OK;
}

/**
 * Helper to validate row group and column indices and retrieve the column chunk.
 * Returns NULL on invalid indices and sets the error.
 */
static const parquet_column_chunk_t* reader_get_column_chunk(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_error_t* error) {

    if (!carquet_reader_row_group_index_valid(reader, row_group_index)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_ROW_GROUP_NOT_FOUND,
            "Row group %d not found", row_group_index);
        return NULL;
    }

    const parquet_row_group_t* rg = &reader->metadata.row_groups[row_group_index];

    if (column_index < 0 || column_index >= rg->num_columns) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_COLUMN_NOT_FOUND,
            "Column %d not found in row group %d", column_index, row_group_index);
        return NULL;
    }

    return &rg->columns[column_index];
}

/* ============================================================================
 * Bloom Filter API
 * ============================================================================
 */

carquet_bloom_filter_t* carquet_reader_get_bloom_filter(
    carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_error_t* error) {

    const parquet_column_chunk_t* chunk = reader_get_column_chunk(
        reader, row_group_index, column_index, error);
    if (!chunk) {
        return NULL;
    }

    if (!chunk->has_metadata) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Column chunk has no metadata");
        return NULL;
    }

    const parquet_column_metadata_t* col_meta = &chunk->metadata;

    /* Check if bloom filter is available */
    if (!col_meta->has_bloom_filter_offset || col_meta->bloom_filter_length <= 0) {
        /* No bloom filter -- not an error, just return NULL */
        return NULL;
    }

    /* Read bloom filter data from file */
    const uint8_t* data = NULL;
    bool allocated = false;
    carquet_status_t status = reader_read_bytes(
        reader, col_meta->bloom_filter_offset, col_meta->bloom_filter_length,
        &data, &allocated, error);
    if (status != CARQUET_OK) {
        return NULL;
    }

    /* The bloom filter region is a Thrift-compact BloomFilterHeader followed by
     * the raw filter bit array:
     *   struct BloomFilterHeader {
     *     1: required i32 numBytes;
     *     2: required BloomFilterAlgorithm algorithm;   // BLOCK
     *     3: required BloomFilterHash hash;             // XXHASH
     *     4: required BloomFilterCompression compression;
     *   }
     * We extract numBytes (field 1) and the compression union tag (field 4) so
     * that a header declaring anything other than UNCOMPRESSED is rejected
     * rather than misread as a raw filter. */
    size_t total_len = (size_t)col_meta->bloom_filter_length;
    int32_t num_bytes = 0;
    /* Default: field 4 absent => UNCOMPRESSED (tag 1). */
    int16_t compression_tag = 1;

    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, total_len);
    thrift_read_struct_begin(&dec);
    thrift_type_t ft;
    int16_t fid;
    while (thrift_read_field_begin(&dec, &ft, &fid)) {
        if (fid == 1 && ft == THRIFT_TYPE_I32) {
            num_bytes = thrift_read_i32(&dec);
        } else if (fid == 4 && ft == THRIFT_TYPE_STRUCT) {
            /* BloomFilterCompression union: exactly one field names the set
             * member (1 = UNCOMPRESSED, the only member the spec defines). */
            thrift_read_struct_begin(&dec);
            thrift_type_t ct;
            int16_t cfid;
            compression_tag = 0;
            while (thrift_read_field_begin(&dec, &ct, &cfid)) {
                if (compression_tag == 0) compression_tag = cfid;
                thrift_skip(&dec, ct);
            }
            thrift_read_struct_end(&dec);
        } else {
            thrift_skip(&dec, ft);
        }
    }
    thrift_read_struct_end(&dec);

    size_t header_size = total_len - thrift_decoder_remaining(&dec);

    if (thrift_decoder_has_error(&dec)) {
        if (allocated) carquet_mem_free((void*)data);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Malformed bloom filter header");
        return NULL;
    }

    /* Only UNCOMPRESSED (union tag 1) is defined by the Parquet spec. Reject any
     * other compression rather than feeding compressed bytes to the reader. */
    if (compression_tag != 1) {
        if (allocated) carquet_mem_free((void*)data);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_UNSUPPORTED_CODEC,
            "Unsupported bloom filter compression (tag %d)", (int)compression_tag);
        return NULL;
    }

    /* Validate and read the raw filter data after the header */
    if (num_bytes <= 0 || header_size + (size_t)num_bytes > total_len) {
        if (allocated) carquet_mem_free((void*)data);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Invalid bloom filter header (numBytes=%d, header=%llu, total=%llu)",
            num_bytes, (unsigned long long)header_size, (unsigned long long)total_len);
        return NULL;
    }

    carquet_bloom_filter_t* filter = NULL;
    status = carquet_bloom_filter_read(&filter, data + header_size, (size_t)num_bytes);

    if (allocated) {
        carquet_mem_free((void*)data);
    }

    if (status != CARQUET_OK) {
        CARQUET_SET_ERROR(error, status, "Failed to parse bloom filter");
        return NULL;
    }

    return filter;
}

/* ============================================================================
 * Key-Value Metadata API
 * ============================================================================
 */

int32_t carquet_reader_num_metadata(const carquet_reader_t* reader) {
    return reader->metadata.num_key_value;
}

carquet_status_t carquet_reader_get_metadata(
    const carquet_reader_t* reader,
    int32_t index,
    const char** key,
    const char** value) {

    if (index < 0 || index >= reader->metadata.num_key_value) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const parquet_key_value_t* kv = &reader->metadata.key_value_metadata[index];
    *key = kv->key;
    *value = kv->value;

    return CARQUET_OK;
}

const char* carquet_reader_find_metadata(
    const carquet_reader_t* reader,
    const char* key) {

    for (int32_t i = 0; i < reader->metadata.num_key_value; i++) {
        const parquet_key_value_t* kv = &reader->metadata.key_value_metadata[i];
        if (kv->key && strcmp(kv->key, key) == 0) {
            return kv->value;
        }
    }

    return NULL;
}

/* Resolve a leaf column index to its schema element, or NULL if out of range. */
static const parquet_schema_element_t* reader_column_element(
    const carquet_reader_t* reader, int32_t column_index) {
    if (!reader->schema || column_index < 0 ||
        column_index >= reader->schema->num_leaves) {
        return NULL;
    }
    int32_t elem_idx = reader->schema->leaf_indices[column_index];
    return &reader->schema->elements[elem_idx];
}

int32_t carquet_reader_column_num_metadata(
    const carquet_reader_t* reader,
    int32_t column_index) {
    const parquet_schema_element_t* e = reader_column_element(reader, column_index);
    return e ? e->num_field_metadata : 0;
}

carquet_status_t carquet_reader_column_get_metadata(
    const carquet_reader_t* reader,
    int32_t column_index,
    int32_t index,
    const char** key,
    const char** value) {
    const parquet_schema_element_t* e = reader_column_element(reader, column_index);
    if (!e || index < 0 || index >= e->num_field_metadata) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    *key = e->field_metadata[index].key;
    *value = e->field_metadata[index].value;
    return CARQUET_OK;
}

const char* carquet_reader_column_find_metadata(
    const carquet_reader_t* reader,
    int32_t column_index,
    const char* key) {
    const parquet_schema_element_t* e = reader_column_element(reader, column_index);
    if (!e) return NULL;
    for (int32_t i = 0; i < e->num_field_metadata; i++) {
        if (e->field_metadata[i].key &&
            strcmp(e->field_metadata[i].key, key) == 0) {
            return e->field_metadata[i].value;
        }
    }
    return NULL;
}

carquet_arrow_type_refinement_t carquet_reader_column_arrow_type_refinement(
    const carquet_reader_t* reader,
    int32_t column_index) {
    const parquet_schema_element_t* e = reader_column_element(reader, column_index);
    if (!e) return CARQUET_ARROW_REFINE_NONE;
    return (carquet_arrow_type_refinement_t)e->arrow_type_refinement;
}

/* ============================================================================
 * Column Chunk Metadata API
 * ============================================================================
 */

carquet_status_t carquet_reader_column_chunk_metadata(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_column_chunk_metadata_t* metadata) {

    const parquet_column_chunk_t* chunk = reader_get_column_chunk(
        reader, row_group_index, column_index, NULL);
    if (!chunk) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    if (!chunk->has_metadata) {
        return CARQUET_ERROR_INVALID_METADATA;
    }

    const parquet_column_metadata_t* col_meta = &chunk->metadata;

    memset(metadata, 0, sizeof(*metadata));

    metadata->type = col_meta->type;
    metadata->codec = col_meta->codec;
    metadata->num_values = col_meta->num_values;
    metadata->total_compressed_size = col_meta->total_compressed_size;
    metadata->total_uncompressed_size = col_meta->total_uncompressed_size;
    metadata->data_page_offset = col_meta->data_page_offset;

    metadata->has_dictionary_page = col_meta->has_dictionary_page_offset;
    metadata->dictionary_page_offset = col_meta->has_dictionary_page_offset
        ? col_meta->dictionary_page_offset : 0;

    /* Copy encodings (up to 4) */
    metadata->num_encodings = col_meta->num_encodings < 4
        ? col_meta->num_encodings : 4;
    for (int32_t i = 0; i < metadata->num_encodings; i++) {
        metadata->encodings[i] = col_meta->encodings[i];
    }

    /* Feature availability flags */
    metadata->has_bloom_filter = col_meta->has_bloom_filter_offset
        && col_meta->bloom_filter_length > 0;
    metadata->has_column_index = chunk->has_column_index_offset
        && chunk->has_column_index_length && chunk->column_index_length > 0;
    metadata->has_offset_index = chunk->has_offset_index_offset
        && chunk->has_offset_index_length && chunk->offset_index_length > 0;

    return CARQUET_OK;
}

carquet_status_t carquet_reader_geospatial_statistics(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_geospatial_statistics_t* stats) {

    const parquet_column_chunk_t* chunk = reader_get_column_chunk(
        reader, row_group_index, column_index, NULL);
    if (!chunk) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (!chunk->has_metadata ||
        !chunk->metadata.has_geospatial_statistics) {
        return CARQUET_ERROR_INVALID_METADATA;
    }

    const parquet_geospatial_statistics_t* g =
        &chunk->metadata.geospatial_statistics;

    memset(stats, 0, sizeof(*stats));
    stats->has_bbox = g->valid;
    stats->xmin = g->xmin; stats->xmax = g->xmax;
    stats->ymin = g->ymin; stats->ymax = g->ymax;
    stats->has_z = g->has_z; stats->zmin = g->zmin; stats->zmax = g->zmax;
    stats->has_m = g->has_m; stats->mmin = g->mmin; stats->mmax = g->mmax;

    int32_t n = g->num_types;
    if (n > CARQUET_MAX_GEOSPATIAL_TYPES) n = CARQUET_MAX_GEOSPATIAL_TYPES;
    stats->num_geometry_types = n;
    for (int32_t i = 0; i < n; i++) {
        stats->geometry_types[i] = g->types[i];
    }
    return CARQUET_OK;
}

/* ============================================================================
 * Page Index API (Column Index + Offset Index)
 * ============================================================================
 */

carquet_column_index_t* carquet_reader_get_column_index(
    carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_error_t* error) {

    const parquet_column_chunk_t* chunk = reader_get_column_chunk(
        reader, row_group_index, column_index, error);
    if (!chunk) {
        return NULL;
    }

    /* Check if column index is available */
    if (!chunk->has_column_index_offset || !chunk->has_column_index_length
        || chunk->column_index_length <= 0) {
        /* No column index -- not an error */
        return NULL;
    }

    /* Read column index data from file */
    const uint8_t* data = NULL;
    bool allocated = false;
    carquet_status_t status = reader_read_bytes(
        reader, chunk->column_index_offset, chunk->column_index_length,
        &data, &allocated, error);
    if (status != CARQUET_OK) {
        return NULL;
    }

    /* Parse column index */
    carquet_column_index_t* ci = carquet_column_index_parse(
        data, (size_t)chunk->column_index_length);

    if (allocated) {
        carquet_mem_free((void*)data);
    }

    if (!ci) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Failed to parse column index");
    }

    return ci;
}

carquet_offset_index_t* carquet_reader_get_offset_index(
    carquet_reader_t* reader,
    int32_t row_group_index,
    int32_t column_index,
    carquet_error_t* error) {

    const parquet_column_chunk_t* chunk = reader_get_column_chunk(
        reader, row_group_index, column_index, error);
    if (!chunk) {
        return NULL;
    }

    /* Check if offset index is available */
    if (!chunk->has_offset_index_offset || !chunk->has_offset_index_length
        || chunk->offset_index_length <= 0) {
        /* No offset index -- not an error */
        return NULL;
    }

    /* Read offset index data from file */
    const uint8_t* data = NULL;
    bool allocated = false;
    carquet_status_t status = reader_read_bytes(
        reader, chunk->offset_index_offset, chunk->offset_index_length,
        &data, &allocated, error);
    if (status != CARQUET_OK) {
        return NULL;
    }

    /* Parse offset index */
    carquet_offset_index_t* oi = carquet_offset_index_parse(
        data, (size_t)chunk->offset_index_length);

    if (allocated) {
        carquet_mem_free((void*)data);
    }

    if (!oi) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_METADATA,
            "Failed to parse offset index");
    }

    return oi;
}
