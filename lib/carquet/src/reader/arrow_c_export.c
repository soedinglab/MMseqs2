/**
 * @file arrow_c_export.c
 * @brief Export Carquet schema/row batches to the Arrow C Data Interface.
 *
 * Produces standard `ArrowSchema` / `ArrowArray` structs (see carquet.h) whose
 * memory is owned by the produced struct and released through its `release`
 * callback, so the exported structs stay valid after the source row batch is
 * freed and can be released by any Arrow consumer, from any thread. Buffers
 * are copies allocated with the C standard allocator (malloc/free), except
 * for fixed-width columns served by the batch reader's pipeline: those arrays
 * point into the reader's row-group buffers and hold a reference on them
 * (see export_shared_leaf), which outlives the batch and the reader.
 *
 * Scope: flat columns plus single-level LIST<T> and MAP<K,V> — the nesting a
 * carquet_row_batch_t can represent (the batch reader materializes repeated
 * `max_rep_level == 1` leaves; deeper nesting is served by
 * carquet_reader_read_arrow). STRUCT and deeper-than-single-level nesting
 * return CARQUET_ERROR_NOT_IMPLEMENTED here.
 */

#include <carquet/carquet.h>

#include <stdlib.h>
#include <string.h>

#include "reader_internal.h"  /* struct carquet_schema (leaf_indices) */
#include "arrow_leaf.h"
#include "core/bitpack.h"
#include "worker_pool.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Release callbacks (ownership: producer allocates, consumer releases)
 * ============================================================================
 */

static void arrow_schema_release(struct ArrowSchema* schema) {
    if (!schema || !schema->release) {
        return;
    }
    free((void*)schema->format);
    free((void*)schema->name);
    free((void*)schema->metadata);
    for (int64_t i = 0; i < schema->n_children; i++) {
        struct ArrowSchema* child = schema->children[i];
        if (child) {
            if (child->release) {
                child->release(child);
            }
            free(child);
        }
    }
    free(schema->children);
    if (schema->dictionary) {
        if (schema->dictionary->release) {
            schema->dictionary->release(schema->dictionary);
        }
        free(schema->dictionary);
    }
    schema->release = NULL;
    schema->private_data = NULL;
}

static void arrow_array_release(struct ArrowArray* array) {
    if (!array || !array->release) {
        return;
    }
    /* Every non-NULL entry in buffers[] is an owned malloc (or NULL). */
    if (array->buffers) {
        for (int64_t i = 0; i < array->n_buffers; i++) {
            free((void*)array->buffers[i]);
        }
        free(array->buffers);
    }
    for (int64_t i = 0; i < array->n_children; i++) {
        struct ArrowArray* child = array->children[i];
        if (child) {
            if (child->release) {
                child->release(child);
            }
            free(child);
        }
    }
    free(array->children);
    if (array->dictionary) {
        if (array->dictionary->release) {
            array->dictionary->release(array->dictionary);
        }
        free(array->dictionary);
    }
    array->release = NULL;
    array->private_data = NULL;
}

/* Implemented in batch_reader.c: a reference on the storage behind a batch
 * column the batch reader's pipeline served (NULL when it must be copied),
 * and its release. */
extern void* carquet_row_batch_column_share(const carquet_row_batch_t* batch,
                                            int32_t column_index, bool* bitmap_shared);
extern void carquet_row_batch_share_release(void* share);

/* A leaf array whose buffers point into batch-reader storage: buffers whose
 * bit is set in `owned` are this array's mallocs, the others are views kept
 * alive by the reference `share`. */
typedef struct {
    void* share;
    unsigned owned;
} ex_shared_leaf_t;

static void arrow_array_release_shared(struct ArrowArray* array) {
    if (!array || !array->release) {
        return;
    }
    ex_shared_leaf_t* priv = (ex_shared_leaf_t*)array->private_data;
    if (array->buffers) {
        for (int64_t i = 0; i < array->n_buffers; i++) {
            if (priv && (priv->owned & (1u << i))) free((void*)array->buffers[i]);
        }
        free(array->buffers);
    }
    if (priv) {
        carquet_row_batch_share_release(priv->share);
        free(priv);
    }
    array->release = NULL;
    array->private_data = NULL;
}

/* strdup via the C allocator; NULL input yields NULL. */
static char* c_strdup(const char* s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char* out = (char*)malloc(n);
    if (out) memcpy(out, s, n);
    return out;
}

/* Populate one leaf ArrowSchema child in place. */
static carquet_status_t export_child_schema(
    struct ArrowSchema* child,
    const char* name,
    carquet_physical_type_t pt,
    int32_t type_length,
    const carquet_logical_type_t* lt,
    bool nullable) {

    carquet_arrow_leaf_t plan;
    carquet_status_t st = carquet_arrow_leaf_plan(pt, type_length, lt, &plan);
    if (st != CARQUET_OK) return st;

    child->format = c_strdup(plan.format);
    child->name = c_strdup(name ? name : "");
    child->metadata = NULL;
    child->flags = nullable ? ARROW_FLAG_NULLABLE : 0;
    child->n_children = 0;
    child->children = NULL;
    child->dictionary = NULL;
    child->release = arrow_schema_release;
    child->private_data = NULL;

    if (!child->format || !child->name) {
        arrow_schema_release(child);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    return CARQUET_OK;
}

/* Build a struct ("+s") ArrowSchema from a flat carquet schema. */
static carquet_status_t build_struct_schema(
    const carquet_schema_t* schema,
    struct ArrowSchema* out,
    carquet_error_t* error) {

    int32_t ncols = carquet_schema_num_columns(schema);

    memset(out, 0, sizeof(*out));
    out->format = c_strdup("+s");
    out->name = NULL;
    out->metadata = NULL;
    out->flags = 0;
    out->n_children = ncols;
    out->children = ncols > 0 ? (struct ArrowSchema**)calloc((size_t)ncols,
                                    sizeof(struct ArrowSchema*)) : NULL;
    out->dictionary = NULL;
    out->release = arrow_schema_release;
    out->private_data = NULL;

    if (!out->format || (ncols > 0 && !out->children)) {
        arrow_schema_release(out);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow schema alloc failed");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < ncols; i++) {
        if (carquet_schema_max_rep_level(schema, i) > 0) {
            arrow_schema_release(out);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Arrow export: nested/repeated column %d not supported", i);
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }

        int32_t elem_idx = schema->leaf_indices[i];
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, elem_idx);
        const carquet_logical_type_t* lt = node ? carquet_schema_node_logical_type(node) : NULL;
        int32_t type_length = node ? carquet_schema_node_type_length(node) : 0;
        carquet_field_repetition_t rep =
            node ? carquet_schema_node_repetition(node) : CARQUET_REPETITION_OPTIONAL;
        bool nullable = (rep != CARQUET_REPETITION_REQUIRED);

        struct ArrowSchema* child = (struct ArrowSchema*)calloc(1, sizeof(struct ArrowSchema));
        if (!child) {
            arrow_schema_release(out);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow child alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        out->children[i] = child;

        carquet_status_t st = export_child_schema(
            child, carquet_schema_column_name(schema, i),
            carquet_schema_column_type(schema, i), type_length, lt, nullable);
        if (st != CARQUET_OK) {
            arrow_schema_release(out);
            CARQUET_SET_ERROR(error, st, "Arrow export: unsupported type for column %d", i);
            return st;
        }
    }
    return CARQUET_OK;
}

/* Implemented in arrow_c_read.c: builds a full nested ArrowSchema tree
 * (struct / list / map at any depth). */
extern carquet_status_t carquet_arrow_build_schema_tree(
    const carquet_schema_t* schema, struct ArrowSchema* out, carquet_error_t* error);

carquet_status_t carquet_arrow_export_schema(
    const carquet_schema_t* schema,
    struct ArrowSchema* out,
    carquet_error_t* error) {

    if (!schema || !out) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL schema or out");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    return carquet_arrow_build_schema_tree(schema, out, error);
}

/* ============================================================================
 * Array export
 * ============================================================================
 */

/* Count nulls in the first `n` bits of a validity bitmap (present bit = 1). */
static int64_t count_nulls(const uint8_t* validity, int64_t n) {
    if (!validity || n <= 0) return 0;
    /* Popcount whole 64-bit words, then the tail bytes, then the last
     * partial byte (bits past n are ignored). A per-bit loop here cost
     * ~1 ns per row, which dominated the export of nullable columns. */
    int64_t present = 0;
    int64_t full_bytes = n >> 3;
    int64_t i = 0;
    for (; i + 8 <= full_bytes; i += 8) {
        uint64_t w;
        memcpy(&w, validity + i, sizeof(w));   /* bit order is irrelevant to a count */
        present += carquet_popcount64(w);
    }
    for (; i < full_bytes; i++) {
        present += carquet_popcount32(validity[i]);
    }
    for (int64_t j = full_bytes << 3; j < n; j++) {
        if (validity[j >> 3] & (uint8_t)(1u << (j & 7))) present++;
    }
    return n - present;
}

/* Copy an Arrow validity buffer from a carquet null bitmap (identical layout:
 * LSB-first, present = 1). Returns NULL when the source is NULL (all valid). */
static uint8_t* copy_validity(const uint8_t* src, int64_t n) {
    if (!src) return NULL;
    size_t bytes = (size_t)((n + 7) / 8);
    if (bytes == 0) bytes = 1;
    uint8_t* out = (uint8_t*)malloc(bytes);
    if (out) memcpy(out, src, (size_t)((n + 7) / 8));
    return out;
}

/* Copy n bytes between non-overlapping buffers. Strings are mostly short, and
 * a memcpy call per value dominated the BYTE_ARRAY export: lengths up to 64
 * are copied as two fixed-size, possibly overlapping, head/tail moves that
 * the compiler inlines. Every load and store stays inside [s, s + n) and
 * [d, d + n), so this never reads past the end of a value. */
static inline void ex_copy_bytes(uint8_t* d, const uint8_t* s, size_t n) {
    if (n <= 16) {
        if (n >= 8) {
            uint64_t a, b;
            memcpy(&a, s, 8); memcpy(&b, s + n - 8, 8);
            memcpy(d, &a, 8); memcpy(d + n - 8, &b, 8);
        } else if (n >= 4) {
            uint32_t a, b;
            memcpy(&a, s, 4); memcpy(&b, s + n - 4, 4);
            memcpy(d, &a, 4); memcpy(d + n - 4, &b, 4);
        } else if (n > 0) {
            uint8_t a = s[0], b = s[n >> 1], c = s[n - 1];
            d[0] = a; d[n >> 1] = b; d[n - 1] = c;
        }
    } else if (n <= 32) {
        uint8_t a[16], b[16];
        memcpy(a, s, 16); memcpy(b, s + n - 16, 16);
        memcpy(d, a, 16); memcpy(d + n - 16, b, 16);
    } else if (n <= 64) {
        uint8_t a[32], b[32];
        memcpy(a, s, 32); memcpy(b, s + n - 32, 32);
        memcpy(d, a, 32); memcpy(d + n - 32, b, 32);
    } else {
        memcpy(d, s, n);
    }
}

/* Copy the bytes of values [begin, end) of a BYTE_ARRAY column to their
 * place in an Arrow data buffer laid out by `offsets`. */
static void ex_copy_strings(uint8_t* bytes, const int32_t* offsets,
                            const carquet_byte_array_t* ba, int64_t begin, int64_t end) {
    for (int64_t i = begin; i < end; i++) {
        int32_t len = offsets[i + 1] - offsets[i];
        if (len > 0) ex_copy_bytes(bytes + offsets[i], ba[i].data, (size_t)len);
    }
}

/* A BYTE_ARRAY data copy left to the caller (see export_child_array): the
 * offsets and the data buffer are in place, the bytes still have to be
 * copied, which carquet_arrow_export_batch spreads over threads in row
 * ranges. */
typedef struct {
    bool active;
    const carquet_byte_array_t* src;
    const int32_t* offsets;
    uint8_t* bytes;
    int64_t n;
} ex_defer_t;

/* Rows per string-copy task; columns with at least two of these defer. */
#define EX_COPY_CHUNK_ROWS 65536

/* Populate one leaf ArrowArray child from a carquet batch column. With a
 * non-NULL `defer`, a large BYTE_ARRAY column only gets its offsets: the
 * byte copy is described in *defer for the caller to run. */
static carquet_status_t export_child_array(
    struct ArrowArray* child,
    const void* data,
    const uint8_t* validity,
    int64_t n,
    const carquet_arrow_leaf_t* plan,
    ex_defer_t* defer) {

    memset(child, 0, sizeof(*child));
    child->length = n;
    child->null_count = count_nulls(validity, n);
    child->offset = 0;
    child->n_children = 0;
    child->children = NULL;
    child->dictionary = NULL;
    child->release = arrow_array_release;
    child->private_data = NULL;

    if (plan->kind == CARQUET_ARROW_LEAF_NULL) {
        child->null_count = n;
        child->n_buffers = 0;
        return CARQUET_OK;
    }

    uint8_t* val = copy_validity(validity, n);
    if (validity && !val) return CARQUET_ERROR_OUT_OF_MEMORY;

    if (plan->kind == CARQUET_ARROW_LEAF_BINARY) {
        /* [validity, offsets(int32, n+1), data] */
        /* One pass builds the offsets (a value with no bytes, or a NULL data
         * pointer, has length 0), a second copies the bytes. */
        const carquet_byte_array_t* ba = (const carquet_byte_array_t*)data;
        int32_t* offsets = (int32_t*)malloc((size_t)(n + 1) * sizeof(int32_t));
        const void** buffers = (const void**)malloc(3 * sizeof(void*));
        if (!offsets || !buffers) {
            free(val); free(offsets); free(buffers);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        int64_t pos = 0;
        for (int64_t i = 0; i < n; i++) {
            offsets[i] = (int32_t)pos;
            if (ba[i].length > 0 && ba[i].data) {
                pos += ba[i].length;
                if (pos > INT32_MAX) {
                    free(val); free(offsets); free(buffers);
                    return CARQUET_ERROR_NOT_IMPLEMENTED;  /* needs large-utf8/binary */
                }
            }
        }
        offsets[n] = (int32_t)pos;
        uint8_t* bytes = (uint8_t*)malloc(pos > 0 ? (size_t)pos : 1);
        if (!bytes) {
            free(val); free(offsets); free(buffers);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (defer && n >= 2 * EX_COPY_CHUNK_ROWS) {
            defer->active = true;
            defer->src = ba;
            defer->offsets = offsets;
            defer->bytes = bytes;
            defer->n = n;
        } else {
            ex_copy_strings(bytes, offsets, ba, 0, n);
        }
        buffers[0] = val;
        buffers[1] = offsets;
        buffers[2] = bytes;
        child->n_buffers = 3;
        child->buffers = buffers;
        return CARQUET_OK;
    }

    /* Fixed-width (incl. BOOLEAN, FIXED_LEN_BYTE_ARRAY and the converted
     * DECIMAL / INT96 / narrow INTEGER values): [validity, data]. */
    uint8_t* out_data = NULL;
    if (plan->kind == CARQUET_ARROW_LEAF_BOOL) {
        /* carquet stores 1 byte/value; Arrow wants bit-packed LSB-first. */
        size_t bytes = (size_t)((n + 7) / 8);
        out_data = (uint8_t*)calloc(bytes > 0 ? bytes : 1, 1);
        if (!out_data) { free(val); return CARQUET_ERROR_OUT_OF_MEMORY; }
        const uint8_t* src = (const uint8_t*)data;
        for (int64_t i = 0; i < n; i++) {
            if (src[i]) out_data[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
    } else {
        if (plan->src_width == 0 || plan->dst_width == 0) {
            free(val);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        size_t total = (size_t)n * plan->dst_width;
        out_data = (uint8_t*)malloc(total > 0 ? total : 1);
        if (!out_data) { free(val); return CARQUET_ERROR_OUT_OF_MEMORY; }
        if (total > 0 && data && plan->kind == CARQUET_ARROW_LEAF_COPY) {
            memcpy(out_data, data, total);
        } else if (total > 0 && data) {
            /* Null slots hold no value to convert: zero them. */
            const uint8_t* src = (const uint8_t*)data;
            for (int64_t i = 0; i < n; i++) {
                uint8_t* dst = out_data + (size_t)i * plan->dst_width;
                if (validity && !(validity[i >> 3] & (1u << (i & 7)))) {
                    memset(dst, 0, plan->dst_width);
                } else if (!carquet_arrow_leaf_store(plan, dst, src + (size_t)i * plan->src_width)) {
                    free(val); free(out_data);
                    return CARQUET_ERROR_DECODE;   /* DECIMAL wider than its precision */
                }
            }
        }
    }

    const void** buffers = (const void**)malloc(2 * sizeof(void*));
    if (!buffers) { free(val); free(out_data); return CARQUET_ERROR_OUT_OF_MEMORY; }
    buffers[0] = val;
    buffers[1] = out_data;
    child->n_buffers = 2;
    child->buffers = buffers;
    return CARQUET_OK;
}

/* Export a top-level leaf without copying its values: the child points into
 * the batch reader's storage and holds a reference on it, so it stays valid
 * after the batch is freed or the reader moves on (see br_share_t in
 * batch_reader.c). Only for values Arrow takes as stored (a COPY plan: no
 * BOOLEAN bit packing, no converted value, no BYTE_ARRAY offsets) and buffers aligned to
 * 8 bytes; the validity bitmap is shared when it is such a view too, else
 * copied. Returns false (nothing allocated, `out` untouched) to fall back to
 * the copying export. */
static bool export_shared_leaf(struct ArrowArray* out, const carquet_row_batch_t* batch,
                               int32_t col, const void* data, const uint8_t* validity,
                               int64_t n, const carquet_arrow_leaf_t* plan) {
    if (plan->kind != CARQUET_ARROW_LEAF_COPY || !data || n <= 0 ||
        ((uintptr_t)data & 7) != 0) {
        return false;
    }
    bool bitmap_shared = false;
    void* share = carquet_row_batch_column_share(batch, col, &bitmap_shared);
    if (!share) return false;
    const void** buffers = (const void**)malloc(2 * sizeof(void*));
    ex_shared_leaf_t* priv = (ex_shared_leaf_t*)malloc(sizeof(*priv));
    uint8_t* val = NULL;
    bool view_bitmap = validity && bitmap_shared && ((uintptr_t)validity & 7) == 0;
    if (validity && !view_bitmap) val = copy_validity(validity, n);
    if (!buffers || !priv || (validity && !view_bitmap && !val)) {
        free(buffers); free(priv); free(val);
        carquet_row_batch_share_release(share);
        return false;
    }
    priv->share = share;
    priv->owned = val ? 1u : 0u;   /* buffer 0 when it is a copy */
    buffers[0] = view_bitmap ? (const void*)validity : (const void*)val;
    buffers[1] = data;
    memset(out, 0, sizeof(*out));
    out->length = n;
    out->null_count = count_nulls(validity, n);
    out->offset = 0;
    out->n_buffers = 2;
    out->buffers = buffers;
    out->n_children = 0;
    out->children = NULL;
    out->dictionary = NULL;
    out->release = arrow_array_release_shared;
    out->private_data = priv;
    return true;
}

/* ============================================================================
 * Nested field navigation (mirrors the classification in arrow_c_read.c)
 * ============================================================================
 */
typedef enum { EK_LEAF, EK_STRUCT, EK_LIST, EK_MAP } ex_kind_t;

static ex_kind_t ex_elem_kind(const carquet_schema_t* cs, int32_t e) {
    const parquet_schema_element_t* el = &cs->elements[e];
    if (el->has_type) return EK_LEAF;
    if (el->has_logical_type && el->logical_type.id == CARQUET_LOGICAL_LIST) return EK_LIST;
    if (el->has_logical_type && el->logical_type.id == CARQUET_LOGICAL_MAP) return EK_MAP;
    return EK_STRUCT;
}

/* Fill `out` (capacity cap) with the child element indices of group `g`, in
 * creation order. Returns the count (which may exceed cap). */
static int32_t ex_elem_children(const carquet_schema_t* cs, int32_t g,
                                int32_t* out, int32_t cap) {
    int32_t n = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == g) {
            if (n < cap) out[n] = i;
            n++;
        }
    }
    return n;
}

/* Map a leaf element index to its leaf column ordinal (-1 if not a leaf). */
static int32_t ex_elem_to_leaf(const carquet_schema_t* cs, int32_t e) {
    for (int32_t l = 0; l < cs->num_leaves; l++) {
        if (cs->leaf_indices[l] == e) return l;
    }
    return -1;
}

/* True iff every top-level field is a primitive leaf (the pre-0.7 flat case). */
static bool ex_schema_is_flat(const carquet_schema_t* cs) {
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == 0 && ex_elem_kind(cs, i) != EK_LEAF) return false;
    }
    return true;
}

/* Copy `count` int32 offsets into a fresh owned buffer. */
static int32_t* copy_offsets(const int32_t* src, int64_t count) {
    int32_t* out = (int32_t*)malloc((size_t)count * sizeof(int32_t));
    if (out && src) memcpy(out, src, (size_t)count * sizeof(int32_t));
    return out;
}

/* The Arrow mapping of leaf element `e`. */
static carquet_status_t ex_leaf_plan(const carquet_schema_t* cs, int32_t e,
                                     carquet_arrow_leaf_t* plan) {
    const parquet_schema_element_t* el = &cs->elements[e];
    return carquet_arrow_leaf_plan(el->type, el->type_length,
                                   el->has_logical_type ? &el->logical_type : NULL, plan);
}

/* Build one flattened leaf child array for element `e` from raw values. */
static carquet_status_t build_leaf_child(struct ArrowArray* child,
                                         const carquet_schema_t* cs, int32_t e,
                                         const void* values, const uint8_t* validity,
                                         int64_t n, ex_defer_t* defer) {
    carquet_arrow_leaf_t plan;
    carquet_status_t st = ex_leaf_plan(cs, e, &plan);
    if (st != CARQUET_OK) return st;
    return export_child_array(child, values, validity, n, &plan, defer);
}

/* Export a single top-level field `e` (leaf / single-level LIST / MAP) into a
 * fully owned ArrowArray. On error `out` is left releasable-empty (release
 * NULL) and every temporary allocation is freed, so the caller's top-level
 * arrow_array_release stays leak-free. */
static carquet_status_t export_field_array(
    const carquet_schema_t* cs, int32_t e,
    const carquet_row_batch_t* batch, int64_t num_rows,
    const int32_t* leaf_to_col,
    struct ArrowArray* out, ex_defer_t* defer, carquet_error_t* error) {

/* Batch column holding leaf `leaf` (identity for an unprojected batch). */
#define EX_BATCH_COL(leaf) ((leaf) < 0 ? -1 : leaf_to_col ? leaf_to_col[(leaf)] : (leaf))

    switch (ex_elem_kind(cs, e)) {
    case EK_LEAF: {
        int32_t lc = EX_BATCH_COL(ex_elem_to_leaf(cs, e));
        const void* data = NULL; const uint8_t* validity = NULL; int64_t n = 0;
        carquet_status_t st = carquet_row_batch_column(batch, lc, &data, &validity, &n);
        if (st != CARQUET_OK) {
            CARQUET_SET_ERROR(error, st, "Arrow export: leaf column %d not accessible", lc);
            return st;
        }
        if (n != num_rows) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow export: column %d has %lld values, expected %lld rows",
                lc, (long long)n, (long long)num_rows);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        carquet_arrow_leaf_t plan;
        if (ex_leaf_plan(cs, e, &plan) == CARQUET_OK &&
            export_shared_leaf(out, batch, lc, data, validity, n, &plan)) {
            return CARQUET_OK;
        }
        st = build_leaf_child(out, cs, e, data, validity, n, defer);
        if (st != CARQUET_OK)
            CARQUET_SET_ERROR(error, st, "Arrow export: unsupported type for column %d", lc);
        return st;
    }

    case EK_LIST: {
        int32_t rep[4]; int32_t nr = ex_elem_children(cs, e, rep, 4);
        int32_t elem[4]; int32_t ne = (nr == 1) ? ex_elem_children(cs, rep[0], elem, 4) : 0;
        if (nr != 1 || ne != 1) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow export: malformed LIST group at element %d", e);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        if (ex_elem_kind(cs, elem[0]) != EK_LEAF) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Arrow export: nested LIST element (use carquet_reader_read_arrow)");
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        int32_t lc = EX_BATCH_COL(ex_elem_to_leaf(cs, elem[0]));
        const int32_t* offsets = NULL; int64_t num_lists = 0;
        const void* values = NULL; const uint8_t* value_validity = NULL;
        int64_t num_values = 0; const uint8_t* list_validity = NULL;
        carquet_status_t st = carquet_row_batch_column_list(
            batch, lc, &offsets, &num_lists, &values, &value_validity,
            &num_values, &list_validity);
        if (st != CARQUET_OK) {
            CARQUET_SET_ERROR(error, st, "Arrow export: LIST column %d not accessible", lc);
            return st;
        }
        if (num_lists != num_rows) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow export: LIST column %d has %lld rows, expected %lld",
                lc, (long long)num_lists, (long long)num_rows);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }

        struct ArrowArray* elem_child = (struct ArrowArray*)calloc(1, sizeof(*elem_child));
        const void** buffers = (const void**)calloc(2, sizeof(void*));
        struct ArrowArray** children = (struct ArrowArray**)calloc(1, sizeof(void*));
        int32_t* off_copy = copy_offsets(offsets, num_lists + 1);
        uint8_t* vld = copy_validity(list_validity, num_lists);
        if (!elem_child || !buffers || !children || !off_copy ||
            (list_validity && !vld)) {
            free(elem_child); free(buffers); free(children); free(off_copy); free(vld);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow LIST alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        st = build_leaf_child(elem_child, cs, elem[0], values, value_validity, num_values, NULL);
        if (st != CARQUET_OK) {
            free(elem_child); free(buffers); free(children); free(off_copy); free(vld);
            CARQUET_SET_ERROR(error, st, "Arrow export: LIST element type unsupported");
            return st;
        }

        memset(out, 0, sizeof(*out));
        out->length = num_lists;
        out->null_count = count_nulls(list_validity, num_lists);
        out->offset = 0;
        buffers[0] = vld;        /* validity (may be NULL) */
        buffers[1] = off_copy;   /* int32 offsets, num_lists + 1 entries */
        out->n_buffers = 2;
        out->buffers = buffers;
        children[0] = elem_child;
        out->n_children = 1;
        out->children = children;
        out->release = arrow_array_release;
        return CARQUET_OK;
    }

    case EK_MAP: {
        int32_t kv[4]; int32_t nkv = ex_elem_children(cs, e, kv, 4);
        int32_t pair[4]; int32_t np = (nkv == 1) ? ex_elem_children(cs, kv[0], pair, 4) : 0;
        if (nkv != 1 || np != 2) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow export: malformed MAP group at element %d", e);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        if (ex_elem_kind(cs, pair[0]) != EK_LEAF || ex_elem_kind(cs, pair[1]) != EK_LEAF) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Arrow export: nested MAP key/value (use carquet_reader_read_arrow)");
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        int32_t kc = EX_BATCH_COL(ex_elem_to_leaf(cs, pair[0]));
        int32_t vc = EX_BATCH_COL(ex_elem_to_leaf(cs, pair[1]));
        const int32_t* off_k = NULL; int64_t nl_k = 0; const void* keys = NULL;
        const uint8_t* key_validity = NULL; int64_t nk = 0; const uint8_t* map_null_k = NULL;
        const int32_t* off_v = NULL; int64_t nl_v = 0; const void* vals = NULL;
        const uint8_t* val_validity = NULL; int64_t nv = 0; const uint8_t* map_null_v = NULL;
        carquet_status_t st = carquet_row_batch_column_list(
            batch, kc, &off_k, &nl_k, &keys, &key_validity, &nk, &map_null_k);
        if (st == CARQUET_OK)
            st = carquet_row_batch_column_list(
                batch, vc, &off_v, &nl_v, &vals, &val_validity, &nv, &map_null_v);
        if (st != CARQUET_OK) {
            CARQUET_SET_ERROR(error, st, "Arrow export: MAP columns %d/%d not accessible", kc, vc);
            return st;
        }
        if (nl_k != num_rows || nl_k != nl_v || nk != nv) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow export: MAP key/value shape mismatch");
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }

        struct ArrowArray* key_child = (struct ArrowArray*)calloc(1, sizeof(*key_child));
        struct ArrowArray* val_child = (struct ArrowArray*)calloc(1, sizeof(*val_child));
        struct ArrowArray* entries = (struct ArrowArray*)calloc(1, sizeof(*entries));
        struct ArrowArray** ent_children = (struct ArrowArray**)calloc(2, sizeof(void*));
        const void** ent_buffers = (const void**)calloc(1, sizeof(void*));
        struct ArrowArray** map_children = (struct ArrowArray**)calloc(1, sizeof(void*));
        const void** map_buffers = (const void**)calloc(2, sizeof(void*));
        int32_t* off_copy = copy_offsets(off_k, nl_k + 1);
        uint8_t* map_vld = copy_validity(map_null_k, nl_k);
        if (!key_child || !val_child || !entries || !ent_children || !ent_buffers ||
            !map_children || !map_buffers || !off_copy || (map_null_k && !map_vld)) {
            free(key_child); free(val_child); free(entries); free(ent_children);
            free(ent_buffers); free(map_children); free(map_buffers);
            free(off_copy); free(map_vld);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow MAP alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        st = build_leaf_child(key_child, cs, pair[0], keys, key_validity, nk, NULL);
        if (st == CARQUET_OK)
            st = build_leaf_child(val_child, cs, pair[1], vals, val_validity, nv, NULL);
        if (st != CARQUET_OK) {
            /* key_child may already own buffers; release it before discarding. */
            if (key_child->release) arrow_array_release(key_child);
            free(key_child); free(val_child); free(entries); free(ent_children);
            free(ent_buffers); free(map_children); free(map_buffers);
            free(off_copy); free(map_vld);
            CARQUET_SET_ERROR(error, st, "Arrow export: MAP key/value type unsupported");
            return st;
        }

        /* entries: non-nullable struct { key, value }, one row per map entry. */
        memset(entries, 0, sizeof(*entries));
        entries->length = nk;
        entries->null_count = 0;
        entries->offset = 0;
        ent_children[0] = key_child;
        ent_children[1] = val_child;
        entries->n_children = 2;
        entries->children = ent_children;
        ent_buffers[0] = NULL;        /* struct validity absent (non-null) */
        entries->n_buffers = 1;
        entries->buffers = ent_buffers;
        entries->release = arrow_array_release;

        memset(out, 0, sizeof(*out));
        out->length = nl_k;
        out->null_count = count_nulls(map_null_k, nl_k);
        out->offset = 0;
        map_buffers[0] = map_vld;    /* map-level validity (may be NULL) */
        map_buffers[1] = off_copy;   /* int32 offsets, nl_k + 1 entries */
        out->n_buffers = 2;
        out->buffers = map_buffers;
        map_children[0] = entries;
        out->n_children = 1;
        out->children = map_children;
        out->release = arrow_array_release;
        return CARQUET_OK;
    }

    case EK_STRUCT:
    default:
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "Arrow export: STRUCT columns not supported by the batch bridge "
            "(use carquet_reader_read_arrow)");
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
#undef EX_BATCH_COL
}

/* Leaf column ordinals [first, first + count) covered by element `e`. */
static void ex_leaf_span(const carquet_schema_t* cs, int32_t e, int32_t* first, int32_t* count) {
    *first = -1; *count = 0;
    for (int32_t l = 0; l < cs->num_leaves; l++) {
        int32_t x = cs->leaf_indices[l];
        while (x > 0 && x != e) x = cs->parent_indices[x];
        if (x == e) {
            if (*first < 0) *first = l;
            (*count)++;
        }
    }
}

typedef struct {
    const carquet_schema_t* cs;
    int32_t e;
    const carquet_row_batch_t* batch;
    int64_t num_rows;
    const int32_t* leaf_to_col;
    struct ArrowArray* out;
    ex_defer_t defer;
    carquet_status_t status;
    carquet_error_t error;
} ex_field_task_t;

static void ex_field_task(void* arg) {
    ex_field_task_t* t = (ex_field_task_t*)arg;
    carquet_error_t err = CARQUET_ERROR_INIT;
    t->status = export_field_array(t->cs, t->e, t->batch, t->num_rows, t->leaf_to_col,
                                   t->out, &t->defer, &err);
    t->error = err;
}

/* Rows [begin, end) of a deferred string copy. */
typedef struct {
    const ex_defer_t* d;
    int64_t begin;
    int64_t end;
} ex_copy_task_t;

static void ex_copy_task(void* arg) {
    const ex_copy_task_t* c = (const ex_copy_task_t*)arg;
    ex_copy_strings(c->d->bytes, c->d->offsets, c->d->src, c->begin, c->end);
}

/* Run the string copies deferred by the field tasks, cut into row ranges and
 * spread over the threads (OpenMP, or the batch reader's column pool). The
 * ranges write disjoint byte spans, so they need no synchronisation. When
 * the range list cannot be allocated the copies run serially. */
static void ex_run_deferred_copies(const ex_field_task_t* tasks, int32_t nfields,
                                   const carquet_row_batch_t* batch) {
    int64_t count = 0;
    for (int32_t f = 0; f < nfields; f++) {
        if (tasks[f].status == CARQUET_OK && tasks[f].defer.active) {
            count += (tasks[f].defer.n + EX_COPY_CHUNK_ROWS - 1) / EX_COPY_CHUNK_ROWS;
        }
    }
    if (count == 0) return;
    ex_copy_task_t* jobs = count <= INT32_MAX
        ? (ex_copy_task_t*)malloc((size_t)count * sizeof(*jobs)) : NULL;
    if (!jobs) {
        for (int32_t f = 0; f < nfields; f++) {
            const ex_defer_t* d = &tasks[f].defer;
            if (tasks[f].status == CARQUET_OK && d->active)
                ex_copy_strings(d->bytes, d->offsets, d->src, 0, d->n);
        }
        return;
    }
    int32_t k = 0;
    for (int32_t f = 0; f < nfields; f++) {
        const ex_defer_t* d = &tasks[f].defer;
        if (tasks[f].status != CARQUET_OK || !d->active) continue;
        for (int64_t b = 0; b < d->n; b += EX_COPY_CHUNK_ROWS) {
            jobs[k].d = d;
            jobs[k].begin = b;
            jobs[k].end = b + EX_COPY_CHUNK_ROWS < d->n ? b + EX_COPY_CHUNK_ROWS : d->n;
            k++;
        }
    }
#ifdef _OPENMP
    (void)batch;
    int32_t j;
    #pragma omp parallel for schedule(dynamic, 1) if(k > 1)
    for (j = 0; j < k; j++) ex_copy_task(&jobs[j]);
#else
    carquet_worker_pool_t* pool = k > 1 ? carquet_row_batch_pool(batch) : NULL;
    void** args = pool ? (void**)malloc((size_t)k * sizeof(void*)) : NULL;
    if (args) {
        for (int32_t j = 0; j < k; j++) args[j] = &jobs[j];
        carquet_worker_pool_parallel_for(pool, ex_copy_task, args, k);
        free(args);
    } else {
        for (int32_t j = 0; j < k; j++) ex_copy_task(&jobs[j]);
    }
#endif
    free(jobs);
}

carquet_status_t carquet_arrow_export_batch(
    const carquet_row_batch_t* batch,
    const carquet_schema_t* schema,
    struct ArrowSchema* out_schema,
    struct ArrowArray* out_array,
    carquet_error_t* error) {

    if (!batch || !schema || !out_array) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL argument");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int32_t ncols = carquet_row_batch_num_columns(batch);
    int32_t nleaves = carquet_schema_num_columns(schema);
    const int32_t* projection = carquet_row_batch_projection(batch);
    if (!projection && ncols != nleaves) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "batch has %d columns but schema has %d leaves (projection not supported)",
            ncols, nleaves);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int64_t num_rows = carquet_row_batch_num_rows(batch);

    /* Walk the top-level schema fields. For a flat schema each field is a leaf
     * and this is 1:1 with the batch columns (output byte-identical to before);
     * a LIST / MAP field consumes one / two leaf columns and expands into an
     * Arrow list / map child. */
    int32_t top[1024];
    int32_t ntop = ex_elem_children(schema, 0, top, 1024);
    if (ntop > 1024) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow export: too many top-level fields (%d)", ntop);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* A batch from carquet_batch_reader_next() may be projected or reordered:
     * map each file leaf to its batch column, then export the top-level fields
     * whose leaves are all present, ordered by where the projection first
     * names them. Without a projection every field is exported in schema order. */
    int32_t sel[1024];          /* selected top-level field ordinals, in output order */
    int32_t nfields = ntop;
    int32_t* leaf_to_col = NULL;
    if (projection) {
        leaf_to_col = (int32_t*)malloc((size_t)(nleaves > 0 ? nleaves : 1) * sizeof(int32_t));
        if (!leaf_to_col) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow export: alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        for (int32_t l = 0; l < nleaves; l++) leaf_to_col[l] = -1;
        for (int32_t c = ncols - 1; c >= 0; c--) {   /* first occurrence wins */
            int32_t l = projection[c];
            if (l < 0 || l >= nleaves) {
                free(leaf_to_col);
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Arrow export: projected leaf %d out of range", l);
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            leaf_to_col[l] = c;
        }
        int32_t first_col[1024];
        nfields = 0;
        for (int32_t f = 0; f < ntop; f++) {
            int32_t first = 0, count = 0, present = 0, min_col = INT32_MAX;
            ex_leaf_span(schema, top[f], &first, &count);
            for (int32_t l = first; count > 0 && l < first + count; l++) {
                if (leaf_to_col[l] >= 0) {
                    present++;
                    if (leaf_to_col[l] < min_col) min_col = leaf_to_col[l];
                }
            }
            if (present == 0) continue;
            if (present != count) {
                free(leaf_to_col);
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Arrow export: projection covers %d of the %d leaves of field %d; "
                    "project all of them or none", present, count, f);
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
            /* insertion sort by first projected column */
            int32_t k = nfields++;
            while (k > 0 && first_col[k - 1] > min_col) {
                first_col[k] = first_col[k - 1]; sel[k] = sel[k - 1]; k--;
            }
            first_col[k] = min_col; sel[k] = f;
        }
    } else {
        for (int32_t f = 0; f < ntop; f++) sel[f] = f;
    }

    memset(out_array, 0, sizeof(*out_array));
    out_array->length = num_rows;
    out_array->null_count = 0;
    out_array->offset = 0;
    out_array->n_children = nfields;
    out_array->children = nfields > 0 ? (struct ArrowArray**)calloc((size_t)nfields,
                                            sizeof(struct ArrowArray*)) : NULL;
    out_array->dictionary = NULL;
    out_array->release = arrow_array_release;
    out_array->private_data = NULL;
    /* struct array carries a single (absent) validity buffer */
    out_array->n_buffers = 1;
    out_array->buffers = (const void**)calloc(1, sizeof(void*));

    if (!out_array->buffers || (nfields > 0 && !out_array->children)) {
        arrow_array_release(out_array);
        free(leaf_to_col);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow array alloc failed");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t f = 0; f < nfields; f++) {
        out_array->children[f] = (struct ArrowArray*)calloc(1, sizeof(struct ArrowArray));
        if (!out_array->children[f]) {
            arrow_array_release(out_array);
            free(leaf_to_col);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow child alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
    }

    /* Fields are independent copies: export them in parallel when the batch
     * is large enough for the fan-out to pay (OpenMP, or the batch reader's
     * column pool in builds without it). */
    ex_field_task_t* tasks = (ex_field_task_t*)calloc((size_t)(nfields > 0 ? nfields : 1), sizeof(*tasks));
    if (!tasks) {
        arrow_array_release(out_array);
        free(leaf_to_col);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow export: alloc failed");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    for (int32_t f = 0; f < nfields; f++) {
        ex_field_task_t* t = &tasks[f];
        t->cs = schema;
        t->e = top[sel[f]];
        t->batch = batch;
        t->num_rows = num_rows;
        t->leaf_to_col = leaf_to_col;
        t->out = out_array->children[f];
        t->status = CARQUET_OK;
    }
    bool parallel = nfields > 1 && num_rows * (int64_t)nfields >= CARQUET_PARALLEL_MIN_VALUES;
#ifdef _OPENMP
    if (parallel) {
        int32_t f;
        #pragma omp parallel for schedule(dynamic, 1)
        for (f = 0; f < nfields; f++) ex_field_task(&tasks[f]);
    } else {
        for (int32_t f = 0; f < nfields; f++) ex_field_task(&tasks[f]);
    }
#else
    carquet_worker_pool_t* pool = parallel ? carquet_row_batch_pool(batch) : NULL;
    void** args = pool ? (void**)malloc((size_t)nfields * sizeof(void*)) : NULL;
    if (args) {
        for (int32_t f = 0; f < nfields; f++) args[f] = &tasks[f];
        carquet_worker_pool_parallel_for(pool, ex_field_task, args, nfields);
        free(args);
    } else {
        for (int32_t f = 0; f < nfields; f++) ex_field_task(&tasks[f]);
    }
#endif
    bool all_ok = true;
    for (int32_t f = 0; f < nfields; f++) {
        if (tasks[f].status != CARQUET_OK) all_ok = false;
    }
    if (all_ok) ex_run_deferred_copies(tasks, nfields, batch);
    for (int32_t f = 0; f < nfields; f++) {
        if (tasks[f].status != CARQUET_OK) {
            carquet_status_t st = tasks[f].status;
            if (error) *error = tasks[f].error;
            free(tasks);
            arrow_array_release(out_array);
            free(leaf_to_col);
            return st;
        }
    }
    free(tasks);
    free(leaf_to_col);

    /* Schema last, so an earlier array failure leaves nothing for the caller to
     * release. A flat schema keeps the byte-identical flat builder; a nested one
     * uses the recursive tree builder (shared with carquet_arrow_export_schema). */
    if (out_schema) {
        carquet_status_t st = ex_schema_is_flat(schema)
            ? build_struct_schema(schema, out_schema, error)
            : carquet_arrow_build_schema_tree(schema, out_schema, error);
        if (st == CARQUET_OK && projection) {
            st = carquet_arrow_schema_select_children(out_schema, sel, nfields);
            if (st != CARQUET_OK) {
                out_schema->release(out_schema);
                CARQUET_SET_ERROR(error, st, "Arrow export: schema projection failed");
            }
        }
        if (st != CARQUET_OK) {
            arrow_array_release(out_array);
            return st;
        }
    }

    return CARQUET_OK;
}
