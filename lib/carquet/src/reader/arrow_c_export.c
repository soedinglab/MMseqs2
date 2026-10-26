/**
 * @file arrow_c_export.c
 * @brief Export Carquet schema/row batches to the Arrow C Data Interface.
 *
 * Produces standard `ArrowSchema` / `ArrowArray` structs (see carquet.h) whose
 * memory is fully owned by the produced struct and released through its
 * `release` callback. All buffers are independent copies allocated with the
 * C standard allocator (malloc/free) so the exported structs stay valid after
 * the source row batch is freed and can be released by any Arrow consumer
 * without knowledge of Carquet's internal allocator.
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

/* ============================================================================
 * Carquet type -> Arrow format string
 * ============================================================================
 * Only lossless, well-defined mappings refine the base physical format with a
 * logical annotation. Types that Arrow can only express via decimal128 /
 * extension metadata (DECIMAL, UUID, INTERVAL, GEOMETRY, ...) fall back to the
 * underlying physical storage format (e.g. "w:16"), which is a correct — if
 * un-annotated — representation of the bytes.
 */
static carquet_status_t arrow_format_string(
    carquet_physical_type_t pt,
    int32_t type_length,
    const carquet_logical_type_t* lt,
    char* buf,
    size_t buf_size) {

    carquet_logical_type_id_t lid = lt ? lt->id : CARQUET_LOGICAL_UNKNOWN;

    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN:
        snprintf(buf, buf_size, "b");
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT32:
        if (lid == CARQUET_LOGICAL_DATE) {
            snprintf(buf, buf_size, "tdD");   /* date32[days] */
        } else if (lid == CARQUET_LOGICAL_TIME) {
            snprintf(buf, buf_size, "ttm");   /* time32[ms] */
        } else if (lid == CARQUET_LOGICAL_INTEGER) {
            int bw = lt->params.integer.bit_width;
            bool s = lt->params.integer.is_signed;
            if (bw == 8)       snprintf(buf, buf_size, s ? "c" : "C");
            else if (bw == 16) snprintf(buf, buf_size, s ? "s" : "S");
            else               snprintf(buf, buf_size, s ? "i" : "I");
        } else {
            snprintf(buf, buf_size, "i");
        }
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT64:
        if (lid == CARQUET_LOGICAL_TIMESTAMP) {
            char u = lt->params.timestamp.unit == CARQUET_TIME_UNIT_MILLIS ? 'm'
                   : lt->params.timestamp.unit == CARQUET_TIME_UNIT_MICROS ? 'u'
                   : 'n';
            snprintf(buf, buf_size, "ts%c:%s", u,
                     lt->params.timestamp.is_adjusted_to_utc ? "UTC" : "");
        } else if (lid == CARQUET_LOGICAL_TIME) {
            char u = lt->params.time.unit == CARQUET_TIME_UNIT_MICROS ? 'u' : 'n';
            snprintf(buf, buf_size, "tt%c", u);   /* time64[us|ns] */
        } else if (lid == CARQUET_LOGICAL_INTEGER) {
            snprintf(buf, buf_size, lt->params.integer.is_signed ? "l" : "L");
        } else {
            snprintf(buf, buf_size, "l");
        }
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT96:
        snprintf(buf, buf_size, "w:12");   /* fixed_size_binary[12] */
        return CARQUET_OK;

    case CARQUET_PHYSICAL_FLOAT:
        snprintf(buf, buf_size, "f");
        return CARQUET_OK;

    case CARQUET_PHYSICAL_DOUBLE:
        snprintf(buf, buf_size, "g");
        return CARQUET_OK;

    case CARQUET_PHYSICAL_BYTE_ARRAY:
        if (lid == CARQUET_LOGICAL_STRING || lid == CARQUET_LOGICAL_ENUM ||
            lid == CARQUET_LOGICAL_JSON) {
            snprintf(buf, buf_size, "u");   /* utf8 */
        } else {
            snprintf(buf, buf_size, "z");   /* binary */
        }
        return CARQUET_OK;

    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
        if (lid == CARQUET_LOGICAL_FLOAT16) {
            snprintf(buf, buf_size, "e");   /* halffloat */
        } else {
            if (type_length <= 0) return CARQUET_ERROR_INVALID_ARGUMENT;
            snprintf(buf, buf_size, "w:%d", type_length);
        }
        return CARQUET_OK;

    default:
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
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

    char fmt[32];
    carquet_status_t st = arrow_format_string(pt, type_length, lt, fmt, sizeof(fmt));
    if (st != CARQUET_OK) return st;

    child->format = c_strdup(fmt);
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
    if (!validity) return 0;
    int64_t nulls = 0;
    for (int64_t i = 0; i < n; i++) {
        if (!(validity[i >> 3] & (uint8_t)(1u << (i & 7)))) nulls++;
    }
    return nulls;
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

/* Arrow value width of an INT32 column annotated INTEGER(8|16): the schema
 * exports it as int8/uint8/int16/uint16, so the values buffer must hold 1- or
 * 2-byte elements rather than carquet's 4-byte physical values. 0 otherwise. */
static size_t arrow_narrow_int_width(carquet_physical_type_t pt,
                                     const carquet_logical_type_t* lt) {
    if (pt != CARQUET_PHYSICAL_INT32 || !lt || lt->id != CARQUET_LOGICAL_INTEGER) return 0;
    if (lt->params.integer.bit_width == 8) return 1;
    if (lt->params.integer.bit_width == 16) return 2;
    return 0;
}

/* Truncate `n` INT32 values to `width`-byte elements (two's complement, so the
 * same code serves the signed and unsigned Arrow types). */
static void arrow_narrow_int_store(uint8_t* dst, const uint8_t* src_i32, size_t width) {
    int32_t v;
    memcpy(&v, src_i32, sizeof(v));
    if (width == 1) {
        dst[0] = (uint8_t)v;
    } else {
        uint16_t t = (uint16_t)v;
        memcpy(dst, &t, sizeof(t));
    }
}

/* Populate one leaf ArrowArray child from a carquet batch column. */
static carquet_status_t export_child_array(
    struct ArrowArray* child,
    const void* data,
    const uint8_t* validity,
    int64_t n,
    carquet_physical_type_t pt,
    int32_t type_length,
    size_t narrow) {

    memset(child, 0, sizeof(*child));
    child->length = n;
    child->null_count = count_nulls(validity, n);
    child->offset = 0;
    child->n_children = 0;
    child->children = NULL;
    child->dictionary = NULL;
    child->release = arrow_array_release;
    child->private_data = NULL;

    uint8_t* val = copy_validity(validity, n);
    if (validity && !val) return CARQUET_ERROR_OUT_OF_MEMORY;

    if (pt == CARQUET_PHYSICAL_BYTE_ARRAY) {
        /* [validity, offsets(int32, n+1), data] */
        const carquet_byte_array_t* ba = (const carquet_byte_array_t*)data;
        int64_t total = 0;
        for (int64_t i = 0; i < n; i++) {
            if (ba[i].length > 0) total += ba[i].length;
        }
        if (total > INT32_MAX) {
            free(val);
            return CARQUET_ERROR_NOT_IMPLEMENTED;  /* needs large-utf8/binary */
        }
        int32_t* offsets = (int32_t*)malloc((size_t)(n + 1) * sizeof(int32_t));
        uint8_t* bytes = (uint8_t*)malloc(total > 0 ? (size_t)total : 1);
        const void** buffers = (const void**)malloc(3 * sizeof(void*));
        if (!offsets || !bytes || !buffers) {
            free(val); free(offsets); free(bytes); free(buffers);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        int32_t pos = 0;
        for (int64_t i = 0; i < n; i++) {
            offsets[i] = pos;
            if (ba[i].length > 0 && ba[i].data) {
                memcpy(bytes + pos, ba[i].data, (size_t)ba[i].length);
                pos += ba[i].length;
            }
        }
        offsets[n] = pos;
        buffers[0] = val;
        buffers[1] = offsets;
        buffers[2] = bytes;
        child->n_buffers = 3;
        child->buffers = buffers;
        return CARQUET_OK;
    }

    /* Fixed-width primitive (incl. BOOLEAN and FIXED_LEN_BYTE_ARRAY):
     * [validity, data]. */
    uint8_t* out_data = NULL;
    if (pt == CARQUET_PHYSICAL_BOOLEAN) {
        /* carquet stores 1 byte/value; Arrow wants bit-packed LSB-first. */
        size_t bytes = (size_t)((n + 7) / 8);
        out_data = (uint8_t*)calloc(bytes > 0 ? bytes : 1, 1);
        if (!out_data) { free(val); return CARQUET_ERROR_OUT_OF_MEMORY; }
        const uint8_t* src = (const uint8_t*)data;
        for (int64_t i = 0; i < n; i++) {
            if (src[i]) out_data[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
    } else {
        size_t stride;
        switch (pt) {
        case CARQUET_PHYSICAL_INT32:  stride = 4; break;
        case CARQUET_PHYSICAL_INT64:  stride = 8; break;
        case CARQUET_PHYSICAL_INT96:  stride = 12; break;
        case CARQUET_PHYSICAL_FLOAT:  stride = 4; break;
        case CARQUET_PHYSICAL_DOUBLE: stride = 8; break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            if (type_length <= 0) { free(val); return CARQUET_ERROR_INVALID_ARGUMENT; }
            stride = (size_t)type_length; break;
        default:
            free(val);
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        size_t total = (size_t)n * (narrow ? narrow : stride);
        out_data = (uint8_t*)malloc(total > 0 ? total : 1);
        if (!out_data) { free(val); return CARQUET_ERROR_OUT_OF_MEMORY; }
        if (total > 0 && data && narrow) {
            for (int64_t i = 0; i < n; i++) {
                arrow_narrow_int_store(out_data + (size_t)i * narrow,
                                       (const uint8_t*)data + (size_t)i * stride, narrow);
            }
        } else if (total > 0 && data) {
            memcpy(out_data, data, total);
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

/* Build one flattened leaf child array for element `e` from raw values. */
static carquet_status_t build_leaf_child(struct ArrowArray* child,
                                         const carquet_schema_t* cs, int32_t e,
                                         const void* values, const uint8_t* validity,
                                         int64_t n) {
    const parquet_schema_element_t* el = &cs->elements[e];
    return export_child_array(child, values, validity, n, el->type, el->type_length,
        arrow_narrow_int_width(el->type, el->has_logical_type ? &el->logical_type : NULL));
}

/* Export a single top-level field `e` (leaf / single-level LIST / MAP) into a
 * fully owned ArrowArray. On error `out` is left releasable-empty (release
 * NULL) and every temporary allocation is freed, so the caller's top-level
 * arrow_array_release stays leak-free. */
static carquet_status_t export_field_array(
    const carquet_schema_t* cs, int32_t e,
    const carquet_row_batch_t* batch, int64_t num_rows,
    struct ArrowArray* out, carquet_error_t* error) {

    switch (ex_elem_kind(cs, e)) {
    case EK_LEAF: {
        int32_t lc = ex_elem_to_leaf(cs, e);
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
        st = build_leaf_child(out, cs, e, data, validity, n);
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
        int32_t lc = ex_elem_to_leaf(cs, elem[0]);
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
        st = build_leaf_child(elem_child, cs, elem[0], values, value_validity, num_values);
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
        int32_t kc = ex_elem_to_leaf(cs, pair[0]);
        int32_t vc = ex_elem_to_leaf(cs, pair[1]);
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
        st = build_leaf_child(key_child, cs, pair[0], keys, key_validity, nk);
        if (st == CARQUET_OK)
            st = build_leaf_child(val_child, cs, pair[1], vals, val_validity, nv);
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
    if (ncols != carquet_schema_num_columns(schema)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "batch has %d columns but schema has %d leaves (projection not supported)",
            ncols, carquet_schema_num_columns(schema));
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    int64_t num_rows = carquet_row_batch_num_rows(batch);

    /* Walk the top-level schema fields. For a flat schema each field is a leaf
     * and this is 1:1 with the batch columns (output byte-identical to before);
     * a LIST / MAP field consumes one / two leaf columns and expands into an
     * Arrow list / map child. */
    int32_t top[1024];
    int32_t nfields = ex_elem_children(schema, 0, top, 1024);
    if (nfields > 1024) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow export: too many top-level fields (%d)", nfields);
        return CARQUET_ERROR_INVALID_ARGUMENT;
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
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow array alloc failed");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }

    for (int32_t f = 0; f < nfields; f++) {
        struct ArrowArray* child = (struct ArrowArray*)calloc(1, sizeof(struct ArrowArray));
        if (!child) {
            arrow_array_release(out_array);
            CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow child alloc failed");
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        out_array->children[f] = child;

        carquet_status_t st = export_field_array(schema, top[f], batch, num_rows, child, error);
        if (st != CARQUET_OK) {
            arrow_array_release(out_array);
            return st;
        }
    }

    /* Schema last, so an earlier array failure leaves nothing for the caller to
     * release. A flat schema keeps the byte-identical flat builder; a nested one
     * uses the recursive tree builder (shared with carquet_arrow_export_schema). */
    if (out_schema) {
        carquet_status_t st = ex_schema_is_flat(schema)
            ? build_struct_schema(schema, out_schema, error)
            : carquet_arrow_build_schema_tree(schema, out_schema, error);
        if (st != CARQUET_OK) {
            arrow_array_release(out_array);
            return st;
        }
    }

    return CARQUET_OK;
}
