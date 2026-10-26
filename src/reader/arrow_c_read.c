/**
 * @file arrow_c_read.c
 * @brief Read a Parquet row group directly into a nested Arrow C Data array.
 *
 * This is the read-side counterpart to the generic Dremel shredder in
 * src/writer/arrow_c_import.c. For a row group it reads every leaf column's
 * raw (repetition, definition, value) stream via the public column-reader API,
 * then reassembles the original nested structure — struct, list, large-list
 * and map, composed to any depth — as a standard `ArrowArray` tree.
 *
 * Reassembly is driven by the Carquet schema element tree. Each node is built
 * from its subtree's leaf streams using two threaded quantities that mirror the
 * shredder exactly:
 *   - exist_def : the definition level at or above which a slot for this node
 *                 materialises (a slot below this level belongs to an empty or
 *                 null ancestor list and is not an element of this node).
 *   - rd        : the repetition depth (rep level of the innermost repeated
 *                 ancestor); a new element of the node begins at any slot with
 *                 rep <= rd.
 * A node instance is present iff def >= exist_def + (nullable ? 1 : 0).
 *
 * All buffers are independent malloc copies owned by the produced structs and
 * released through their `release` callbacks, so the result outlives the
 * reader.
 */

#include <carquet/carquet.h>

#include <stdlib.h>
#include <string.h>

#include "reader_internal.h"  /* struct carquet_schema, struct carquet_reader */

/* ============================================================================
 * Release callbacks (producer allocates with malloc, consumer releases)
 * ============================================================================
 */
static void release_schema(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format); free((void*)s->name); free((void*)s->metadata);
    for (int64_t i = 0; i < s->n_children; i++) {
        if (s->children[i]) { if (s->children[i]->release) s->children[i]->release(s->children[i]); free(s->children[i]); }
    }
    free(s->children);
    s->release = NULL; s->private_data = NULL;
}
static void release_array(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) { for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]); free(a->buffers); }
    for (int64_t i = 0; i < a->n_children; i++) {
        if (a->children[i]) { if (a->children[i]->release) a->children[i]->release(a->children[i]); free(a->children[i]); }
    }
    free(a->children);
    a->release = NULL; a->private_data = NULL;
}

static char* c_strdup(const char* s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char* o = (char*)malloc(n);
    if (o) memcpy(o, s, n);
    return o;
}

/* ============================================================================
 * Carquet schema element tree navigation
 * ============================================================================
 */
typedef enum { K_LEAF, K_STRUCT, K_LIST, K_MAP } elem_kind_t;

static elem_kind_t elem_kind(const carquet_schema_t* cs, int32_t e) {
    const parquet_schema_element_t* el = &cs->elements[e];
    if (el->has_type) return K_LEAF;
    if (el->has_logical_type && el->logical_type.id == CARQUET_LOGICAL_LIST) return K_LIST;
    if (el->has_logical_type && el->logical_type.id == CARQUET_LOGICAL_MAP) return K_MAP;
    return K_STRUCT;
}

static bool elem_nullable(const carquet_schema_t* cs, int32_t e) {
    return cs->elements[e].repetition_type != CARQUET_REPETITION_REQUIRED;
}

/* Fill `out` (capacity cap) with the child element indices of group `g`, in
 * creation order. Returns the count. */
static int32_t elem_children(const carquet_schema_t* cs, int32_t g, int32_t* out, int32_t cap) {
    int32_t n = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == g) {
            if (n < cap) out[n] = i;
            n++;
        }
    }
    return n;
}

/* Number of leaf columns under element `e`. */
static int32_t count_leaves(const carquet_schema_t* cs, int32_t e) {
    if (cs->elements[e].has_type) return 1;
    int32_t total = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == e) total += count_leaves(cs, i);
    }
    return total;
}

/* Map a leaf element index to its leaf column ordinal (-1 if not a leaf). */
static int32_t elem_to_leaf(const carquet_schema_t* cs, int32_t e) {
    for (int32_t l = 0; l < cs->num_leaves; l++) {
        if (cs->leaf_indices[l] == e) return l;
    }
    return -1;
}

/* ============================================================================
 * Carquet type -> Arrow format string
 * ============================================================================
 */
static carquet_status_t arrow_format_string(
    carquet_physical_type_t pt, int32_t type_length,
    const carquet_logical_type_t* lt, char* buf, size_t buf_size) {

    carquet_logical_type_id_t lid = lt ? lt->id : CARQUET_LOGICAL_UNKNOWN;
    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN: snprintf(buf, buf_size, "b"); return CARQUET_OK;
    case CARQUET_PHYSICAL_INT32:
        if (lid == CARQUET_LOGICAL_DATE) snprintf(buf, buf_size, "tdD");
        else if (lid == CARQUET_LOGICAL_TIME) snprintf(buf, buf_size, "ttm");
        else if (lid == CARQUET_LOGICAL_INTEGER) {
            int bw = lt->params.integer.bit_width; bool s = lt->params.integer.is_signed;
            if (bw == 8) snprintf(buf, buf_size, s ? "c" : "C");
            else if (bw == 16) snprintf(buf, buf_size, s ? "s" : "S");
            else snprintf(buf, buf_size, s ? "i" : "I");
        } else snprintf(buf, buf_size, "i");
        return CARQUET_OK;
    case CARQUET_PHYSICAL_INT64:
        if (lid == CARQUET_LOGICAL_TIMESTAMP) {
            char u = lt->params.timestamp.unit == CARQUET_TIME_UNIT_MILLIS ? 'm'
                   : lt->params.timestamp.unit == CARQUET_TIME_UNIT_MICROS ? 'u' : 'n';
            snprintf(buf, buf_size, "ts%c:%s", u, lt->params.timestamp.is_adjusted_to_utc ? "UTC" : "");
        } else if (lid == CARQUET_LOGICAL_TIME) {
            char u = lt->params.time.unit == CARQUET_TIME_UNIT_MICROS ? 'u' : 'n';
            snprintf(buf, buf_size, "tt%c", u);
        } else if (lid == CARQUET_LOGICAL_INTEGER) {
            snprintf(buf, buf_size, lt->params.integer.is_signed ? "l" : "L");
        } else snprintf(buf, buf_size, "l");
        return CARQUET_OK;
    case CARQUET_PHYSICAL_INT96: snprintf(buf, buf_size, "w:12"); return CARQUET_OK;
    case CARQUET_PHYSICAL_FLOAT: snprintf(buf, buf_size, "f"); return CARQUET_OK;
    case CARQUET_PHYSICAL_DOUBLE: snprintf(buf, buf_size, "g"); return CARQUET_OK;
    case CARQUET_PHYSICAL_BYTE_ARRAY:
        if (lid == CARQUET_LOGICAL_STRING || lid == CARQUET_LOGICAL_ENUM || lid == CARQUET_LOGICAL_JSON)
            snprintf(buf, buf_size, "u");
        else snprintf(buf, buf_size, "z");
        return CARQUET_OK;
    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
        if (lid == CARQUET_LOGICAL_FLOAT16) snprintf(buf, buf_size, "e");
        else { if (type_length <= 0) return CARQUET_ERROR_INVALID_ARGUMENT;
               snprintf(buf, buf_size, "w:%d", type_length); }
        return CARQUET_OK;
    default: return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
}

/* ============================================================================
 * Nested ArrowSchema builder (shared with carquet_arrow_export_schema)
 * ============================================================================
 */
static carquet_status_t schema_node(const carquet_schema_t* cs, int32_t e,
                                    const char* name_override,
                                    bool force_nonnull, struct ArrowSchema* out);

static struct ArrowSchema* new_schema_child(void) {
    return (struct ArrowSchema*)calloc(1, sizeof(struct ArrowSchema));
}

static carquet_status_t schema_alloc_children(struct ArrowSchema* out, int32_t n) {
    out->n_children = n;
    out->children = n ? (struct ArrowSchema**)calloc((size_t)n, sizeof(void*)) : NULL;
    if (n && !out->children) return CARQUET_ERROR_OUT_OF_MEMORY;
    return CARQUET_OK;
}

static carquet_status_t schema_node(const carquet_schema_t* cs, int32_t e,
                                    const char* name_override,
                                    bool force_nonnull, struct ArrowSchema* out) {
    const parquet_schema_element_t* el = &cs->elements[e];
    const char* name = name_override ? name_override : (el->name ? el->name : "");
    bool nullable = force_nonnull ? false : elem_nullable(cs, e);

    memset(out, 0, sizeof(*out));
    out->name = c_strdup(name);
    out->flags = nullable ? ARROW_FLAG_NULLABLE : 0;
    out->release = release_schema;
    if (!out->name) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }

    switch (elem_kind(cs, e)) {
    case K_LEAF: {
        char fmt[32];
        const carquet_logical_type_t* lt = el->has_logical_type ? &el->logical_type : NULL;
        carquet_status_t st = arrow_format_string(el->type, el->type_length, lt, fmt, sizeof(fmt));
        if (st != CARQUET_OK) { release_schema(out); return st; }
        out->format = c_strdup(fmt);
        if (!out->format) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        return CARQUET_OK;
    }
    case K_STRUCT: {
        out->format = c_strdup("+s");
        int32_t kids[256]; int32_t nk = elem_children(cs, e, kids, 256);
        if (!out->format || nk > 256) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        if (schema_alloc_children(out, nk) != CARQUET_OK) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        for (int32_t i = 0; i < nk; i++) {
            out->children[i] = new_schema_child();
            if (!out->children[i]) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
            carquet_status_t st = schema_node(cs, kids[i], NULL, false, out->children[i]);
            if (st != CARQUET_OK) { release_schema(out); return st; }
        }
        return CARQUET_OK;
    }
    case K_LIST: {
        /* list group -> repeated "list" group -> element */
        out->format = c_strdup("+l");
        if (!out->format) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        int32_t rep[4]; int32_t nr = elem_children(cs, e, rep, 4);
        if (nr != 1) { release_schema(out); return CARQUET_ERROR_INVALID_ARGUMENT; }
        int32_t elem[4]; int32_t ne = elem_children(cs, rep[0], elem, 4);
        if (ne != 1) { release_schema(out); return CARQUET_ERROR_INVALID_ARGUMENT; }
        if (schema_alloc_children(out, 1) != CARQUET_OK) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->children[0] = new_schema_child();
        if (!out->children[0]) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        carquet_status_t st = schema_node(cs, elem[0], "element", false, out->children[0]);
        if (st != CARQUET_OK) { release_schema(out); return st; }
        return CARQUET_OK;
    }
    case K_MAP: {
        /* map group -> repeated "key_value" -> {key, value} : Arrow "+m" with a
         * non-nullable struct "entries" child holding [key, value]. */
        out->format = c_strdup("+m");
        if (!out->format) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        int32_t kv[4]; int32_t nkv = elem_children(cs, e, kv, 4);
        if (nkv != 1) { release_schema(out); return CARQUET_ERROR_INVALID_ARGUMENT; }
        int32_t pair[4]; int32_t np = elem_children(cs, kv[0], pair, 4);
        if (np != 2) { release_schema(out); return CARQUET_ERROR_INVALID_ARGUMENT; }
        if (schema_alloc_children(out, 1) != CARQUET_OK) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        struct ArrowSchema* entries = new_schema_child();
        out->children[0] = entries;
        if (!entries) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        memset(entries, 0, sizeof(*entries));
        entries->format = c_strdup("+s");
        entries->name = c_strdup("entries");
        entries->flags = 0;               /* entries struct is non-nullable */
        entries->release = release_schema;
        if (!entries->format || !entries->name) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        if (schema_alloc_children(entries, 2) != CARQUET_OK) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        entries->children[0] = new_schema_child();
        entries->children[1] = new_schema_child();
        if (!entries->children[0] || !entries->children[1]) { release_schema(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        carquet_status_t st = schema_node(cs, pair[0], "key", true, entries->children[0]);
        if (st != CARQUET_OK) { release_schema(out); return st; }
        st = schema_node(cs, pair[1], "value", false, entries->children[1]);
        if (st != CARQUET_OK) { release_schema(out); return st; }
        return CARQUET_OK;
    }
    }
    release_schema(out);
    return CARQUET_ERROR_INTERNAL;
}

carquet_status_t carquet_arrow_build_schema_tree(
    const carquet_schema_t* cs, struct ArrowSchema* out, carquet_error_t* error) {

    memset(out, 0, sizeof(*out));
    out->format = c_strdup("+s");
    out->name = NULL;
    out->release = release_schema;
    if (!out->format) { CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "alloc"); return CARQUET_ERROR_OUT_OF_MEMORY; }

    int32_t top[1024]; int32_t nt = elem_children(cs, 0, top, 1024);
    if (nt > 1024) { release_schema(out); CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "too many fields"); return CARQUET_ERROR_INVALID_ARGUMENT; }
    if (schema_alloc_children(out, nt) != CARQUET_OK) { release_schema(out); CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "alloc"); return CARQUET_ERROR_OUT_OF_MEMORY; }
    for (int32_t i = 0; i < nt; i++) {
        out->children[i] = new_schema_child();
        if (!out->children[i]) { release_schema(out); CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "alloc"); return CARQUET_ERROR_OUT_OF_MEMORY; }
        carquet_status_t st = schema_node(cs, top[i], NULL, false, out->children[i]);
        if (st != CARQUET_OK) { release_schema(out); CARQUET_SET_ERROR(error, st, "Arrow export: field %d", i); return st; }
    }
    return CARQUET_OK;
}

/* ============================================================================
 * Leaf reading
 * ============================================================================
 */
typedef struct {
    int16_t* rep;
    int16_t* def;
    int64_t  nslots;
    int16_t  max_def;
    int16_t  max_rep;
    carquet_physical_type_t pt;
    int32_t  type_length;
    size_t   stride;        /* fixed-width byte stride; 0 for byte array */
    size_t   narrow;        /* Arrow element width for INTEGER(8|16), else 0 */
    bool     is_bool;
    bool     is_bytearray;
    void*    values;        /* dense present values */
    int64_t  present;       /* count(def == max_def) */
    carquet_column_reader_t* cr;
} rleaf_t;

/* ============================================================================
 * Array assembly
 * ============================================================================
 */

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
typedef struct {
    const carquet_schema_t* cs;
    rleaf_t* leaves;        /* [num_leaves] */
    int32_t  num_leaves;
    carquet_error_t* error;
} rctx_t;

static struct ArrowArray* new_array(void) {
    struct ArrowArray* a = (struct ArrowArray*)calloc(1, sizeof(struct ArrowArray));
    if (a) a->release = release_array;
    return a;
}

/* Build a bit-packed (LSB-first) validity buffer; returns NULL if all present
 * (null_count 0). *null_count receives the number of unset bits. */
static uint8_t* build_validity(const bool* present, int64_t n, int64_t* null_count) {
    int64_t nulls = 0;
    for (int64_t i = 0; i < n; i++) if (!present[i]) nulls++;
    *null_count = nulls;
    if (nulls == 0) return NULL;
    size_t bytes = (size_t)((n + 7) / 8); if (bytes == 0) bytes = 1;
    uint8_t* v = (uint8_t*)calloc(bytes, 1);
    if (!v) return NULL;
    for (int64_t i = 0; i < n; i++) if (present[i]) v[i >> 3] |= (uint8_t)(1u << (i & 7));
    return v;
}

/* Build the ArrowArray for a primitive leaf. exist_def = slot-exists level. */
static carquet_status_t build_leaf(rctx_t* ctx, int32_t leaf_col, int16_t exist_def,
                                   struct ArrowArray** out) {
    rleaf_t* L = &ctx->leaves[leaf_col];
    /* Count element slots and gather present flags. */
    int64_t len = 0;
    for (int64_t s = 0; s < L->nslots; s++) if (L->def[s] >= exist_def) len++;

    bool* present = (bool*)malloc((size_t)(len > 0 ? len : 1) * sizeof(bool));
    if (!present) return CARQUET_ERROR_OUT_OF_MEMORY;
    int64_t k = 0;
    for (int64_t s = 0; s < L->nslots; s++) {
        if (L->def[s] >= exist_def) present[k++] = (L->def[s] == L->max_def);
    }

    struct ArrowArray* a = new_array();
    if (!a) { free(present); return CARQUET_ERROR_OUT_OF_MEMORY; }
    a->length = len; a->offset = 0;

    int64_t null_count = 0;
    uint8_t* validity = build_validity(present, len, &null_count);
    a->null_count = null_count;

    carquet_status_t rc = CARQUET_OK;
    if (L->is_bytearray) {
        const carquet_byte_array_t* src = (const carquet_byte_array_t*)L->values;
        int64_t total = 0;
        for (int64_t i = 0; i < L->present; i++) total += src[i].length;
        if (total > INT32_MAX) { rc = CARQUET_ERROR_NOT_IMPLEMENTED; goto fail; } /* needs large binary */
        int32_t* offs = (int32_t*)malloc((size_t)(len + 1) * sizeof(int32_t));
        uint8_t* data = (uint8_t*)malloc((size_t)(total > 0 ? total : 1));
        const void** bufs = (const void**)malloc(3 * sizeof(void*));
        if (!offs || !data || !bufs) { free(offs); free(data); free(bufs); rc = CARQUET_ERROR_OUT_OF_MEMORY; goto fail; }
        int32_t pos = 0; int64_t vc = 0;
        for (int64_t i = 0; i < len; i++) {
            offs[i] = pos;
            if (present[i]) {
                const carquet_byte_array_t* b = &src[vc++];
                if (b->length && b->data) { memcpy(data + pos, b->data, b->length); pos += (int32_t)b->length; }
            }
        }
        offs[len] = pos;
        bufs[0] = validity; bufs[1] = offs; bufs[2] = data;
        a->n_buffers = 3; a->buffers = bufs;
    } else if (L->is_bool) {
        size_t bytes = (size_t)((len + 7) / 8); if (bytes == 0) bytes = 1;
        uint8_t* data = (uint8_t*)calloc(bytes, 1);
        const void** bufs = (const void**)malloc(2 * sizeof(void*));
        if (!data || !bufs) { free(data); free(bufs); rc = CARQUET_ERROR_OUT_OF_MEMORY; goto fail; }
        const uint8_t* src = (const uint8_t*)L->values;
        int64_t vc = 0;
        for (int64_t i = 0; i < len; i++) {
            if (present[i]) { if (src[vc++]) data[i >> 3] |= (uint8_t)(1u << (i & 7)); }
        }
        bufs[0] = validity; bufs[1] = data;
        a->n_buffers = 2; a->buffers = bufs;
    } else {
        size_t stride = L->stride;
        size_t ostride = L->narrow ? L->narrow : stride;
        uint8_t* data = (uint8_t*)calloc((size_t)(len > 0 ? len : 1) * ostride, 1);
        const void** bufs = (const void**)malloc(2 * sizeof(void*));
        if (!data || !bufs) { free(data); free(bufs); rc = CARQUET_ERROR_OUT_OF_MEMORY; goto fail; }
        const uint8_t* src = (const uint8_t*)L->values;
        int64_t vc = 0;
        for (int64_t i = 0; i < len; i++) {
            if (!present[i]) continue;
            if (L->narrow) {
                arrow_narrow_int_store(data + (size_t)i * ostride,
                                       src + (size_t)vc * stride, ostride);
            } else {
                memcpy(data + (size_t)i * stride, src + (size_t)vc * stride, stride);
            }
            vc++;
        }
        bufs[0] = validity; bufs[1] = data;
        a->n_buffers = 2; a->buffers = bufs;
    }
    free(present);
    *out = a;
    return CARQUET_OK;
fail:
    free(present); free(validity);
    a->release = NULL; free(a);
    return rc;
}

/*
 * Two threaded quantities (see file header). `def_in` mirrors the write
 * shredder (a present struct hands children Dpres(struct); a list hands its
 * element Dpres(list)+1) and drives presence/band tests. `exist` is the def at
 * which the node's *slot* materialises as an Arrow element — a struct passes it
 * to children unchanged (a null struct still yields a child slot), while a list
 * raises it to Dpres+1 (an empty/null list yields no element slot). They differ
 * only across an OPTIONAL struct.
 */
static carquet_status_t build_node(rctx_t* ctx, int32_t e, int16_t def_in,
                                   int16_t exist, int16_t rd, int32_t base,
                                   struct ArrowArray** out);

/* Build a struct-shaped node with an explicit child element list. Used for
 * plain structs and for a map's synthetic "entries" struct (children = the
 * key/value elements of the REPEATED key_value group). */
static carquet_status_t build_struct_like(rctx_t* ctx, const int32_t* kids, int32_t nk,
                                          bool nullable, int16_t def_in, int16_t exist,
                                          int16_t rd, int32_t base, struct ArrowArray** out) {
    rleaf_t* rep_leaf = &ctx->leaves[base];
    int16_t dpres = (int16_t)(def_in + (nullable ? 1 : 0));

    /* Length + per-instance presence from the representative (leftmost) leaf. */
    int64_t len = 0;
    for (int64_t s = 0; s < rep_leaf->nslots; s++)
        if (rep_leaf->rep[s] <= rd && rep_leaf->def[s] >= exist) len++;

    bool* present = (bool*)malloc((size_t)(len > 0 ? len : 1) * sizeof(bool));
    if (!present) return CARQUET_ERROR_OUT_OF_MEMORY;
    int64_t k = 0;
    for (int64_t s = 0; s < rep_leaf->nslots; s++)
        if (rep_leaf->rep[s] <= rd && rep_leaf->def[s] >= exist)
            present[k++] = (rep_leaf->def[s] >= dpres);

    struct ArrowArray* a = new_array();
    if (!a) { free(present); return CARQUET_ERROR_OUT_OF_MEMORY; }
    a->length = len; a->offset = 0;
    int64_t null_count = 0;
    uint8_t* validity = nullable ? build_validity(present, len, &null_count) : NULL;
    a->null_count = nullable ? null_count : 0;
    free(present);

    const void** bufs = (const void**)malloc(1 * sizeof(void*));
    if (!bufs) { free(validity); a->release = NULL; free(a); return CARQUET_ERROR_OUT_OF_MEMORY; }
    bufs[0] = validity;
    a->n_buffers = 1; a->buffers = bufs;

    a->n_children = nk;
    a->children = nk ? (struct ArrowArray**)calloc((size_t)nk, sizeof(void*)) : NULL;
    if (nk && !a->children) { a->release(a); free(a); return CARQUET_ERROR_OUT_OF_MEMORY; }

    int32_t child_base = base;
    for (int32_t i = 0; i < nk; i++) {
        /* children: def_in = Dpres(struct); exist unchanged (null struct still
         * yields a child slot). */
        carquet_status_t st = build_node(ctx, kids[i], dpres, exist, rd, child_base, &a->children[i]);
        if (st != CARQUET_OK) { release_array(a); free(a); return st; }
        if (a->children[i]->length != len) {
            /* struct children must align 1:1 with the struct's elements */
            release_array(a); free(a);
            return CARQUET_ERROR_INTERNAL;
        }
        child_base += count_leaves(ctx->cs, kids[i]);
    }
    *out = a;
    return CARQUET_OK;
}

/* Build a list/map node. For a list, `child_e` is the element element index;
 * for a map, `map_pair` holds the {key, value} element indices. */
static carquet_status_t build_list_like(rctx_t* ctx, bool nullable, int16_t def_in,
                                        int16_t exist, int16_t rd, int32_t base,
                                        int32_t child_e, const int32_t* map_pair,
                                        struct ArrowArray** out) {
    rleaf_t* rep_leaf = &ctx->leaves[base];
    int16_t dpres = (int16_t)(def_in + (nullable ? 1 : 0));
    int16_t def_in_child = (int16_t)(dpres + 1);   /* repeated group present */
    int16_t rd_child = (int16_t)(rd + 1);

    /* First pass: count list instances (one per parent element) and child
     * elements. A slot is a new list boundary when rep <= rd and the slot
     * belongs to this list level (def >= exist); a child element is a slot with
     * rep <= rd_child and def >= def_in_child. */
    int64_t num_lists = 0, child_total = 0;
    for (int64_t s = 0; s < rep_leaf->nslots; s++) {
        if (rep_leaf->rep[s] <= rd && rep_leaf->def[s] >= exist) num_lists++;
        if (rep_leaf->def[s] >= def_in_child && rep_leaf->rep[s] <= rd_child) child_total++;
    }
    if (num_lists > INT32_MAX || child_total > INT32_MAX) return CARQUET_ERROR_INVALID_ARGUMENT;

    int32_t* offs = (int32_t*)malloc((size_t)(num_lists + 1) * sizeof(int32_t));
    bool* present = (bool*)malloc((size_t)(num_lists > 0 ? num_lists : 1) * sizeof(bool));
    if (!offs || !present) { free(offs); free(present); return CARQUET_ERROR_OUT_OF_MEMORY; }

    int64_t li = -1, cc = 0;
    for (int64_t s = 0; s < rep_leaf->nslots; s++) {
        if (rep_leaf->rep[s] <= rd && rep_leaf->def[s] >= exist) {
            li++;
            offs[li] = (int32_t)cc;
            present[li] = (rep_leaf->def[s] >= dpres);
        }
        if (rep_leaf->def[s] >= def_in_child && rep_leaf->rep[s] <= rd_child) cc++;
    }
    offs[num_lists] = (int32_t)cc;

    struct ArrowArray* a = new_array();
    if (!a) { free(offs); free(present); return CARQUET_ERROR_OUT_OF_MEMORY; }
    a->length = num_lists; a->offset = 0;
    int64_t null_count = 0;
    uint8_t* validity = nullable ? build_validity(present, num_lists, &null_count) : NULL;
    a->null_count = nullable ? null_count : 0;
    free(present);

    const void** bufs = (const void**)malloc(2 * sizeof(void*));
    if (!bufs) { free(validity); free(offs); a->release = NULL; free(a); return CARQUET_ERROR_OUT_OF_MEMORY; }
    bufs[0] = validity; bufs[1] = offs;
    a->n_buffers = 2; a->buffers = bufs;
    a->n_children = 1;
    a->children = (struct ArrowArray**)calloc(1, sizeof(void*));
    if (!a->children) { a->release(a); free(a); return CARQUET_ERROR_OUT_OF_MEMORY; }

    /* child: both def_in and exist become def_in_child (a list raises exist). */
    carquet_status_t st;
    if (map_pair) {
        st = build_struct_like(ctx, map_pair, 2, /*nullable=*/false,
                               def_in_child, def_in_child, rd_child, base, &a->children[0]);
    } else {
        st = build_node(ctx, child_e, def_in_child, def_in_child, rd_child, base, &a->children[0]);
    }
    if (st != CARQUET_OK) { release_array(a); free(a); return st; }
    if (a->children[0]->length != cc) { release_array(a); free(a); return CARQUET_ERROR_INTERNAL; }
    *out = a;
    return CARQUET_OK;
}

static carquet_status_t build_node(rctx_t* ctx, int32_t e, int16_t def_in,
                                   int16_t exist, int16_t rd, int32_t base,
                                   struct ArrowArray** out) {
    const carquet_schema_t* cs = ctx->cs;
    (void)def_in;  /* leaves need only `exist`; groups thread both */
    /* `base` (the leftmost leaf column of this node) must index a real leaf.
     * A malformed / inconsistent schema tree can drive it out of range; every
     * node reads its representative leaf ctx->leaves[base], so guard here. */
    if (base < 0 || base >= ctx->num_leaves) return CARQUET_ERROR_INVALID_ARGUMENT;
    switch (elem_kind(cs, e)) {
    case K_LEAF: {
        int32_t leaf = elem_to_leaf(cs, e);
        if (leaf < 0 || leaf != base) return CARQUET_ERROR_INTERNAL;
        return build_leaf(ctx, leaf, exist, out);
    }
    case K_STRUCT: {
        int32_t kids[256]; int32_t nk = elem_children(cs, e, kids, 256);
        if (nk > 256) return CARQUET_ERROR_INVALID_ARGUMENT;
        return build_struct_like(ctx, kids, nk, elem_nullable(cs, e), def_in, exist, rd, base, out);
    }
    case K_LIST: {
        int32_t rep[4]; if (elem_children(cs, e, rep, 4) != 1) return CARQUET_ERROR_INVALID_ARGUMENT;
        int32_t elem[4]; if (elem_children(cs, rep[0], elem, 4) != 1) return CARQUET_ERROR_INVALID_ARGUMENT;
        return build_list_like(ctx, elem_nullable(cs, e), def_in, exist, rd, base, elem[0], NULL, out);
    }
    case K_MAP: {
        int32_t kv[4]; if (elem_children(cs, e, kv, 4) != 1) return CARQUET_ERROR_INVALID_ARGUMENT;
        int32_t pair[4]; if (elem_children(cs, kv[0], pair, 4) != 2) return CARQUET_ERROR_INVALID_ARGUMENT;
        return build_list_like(ctx, elem_nullable(cs, e), def_in, exist, rd, base, -1, pair, out);
    }
    }
    return CARQUET_ERROR_INTERNAL;
}

/* ============================================================================
 * Public entry point
 * ============================================================================
 */
carquet_status_t carquet_reader_read_arrow(
    carquet_reader_t* reader,
    int32_t row_group_index,
    struct ArrowSchema* out_schema,
    struct ArrowArray* out_array,
    carquet_error_t* error) {

    if (!reader || !out_array) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader or out_array");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int32_t nrg = carquet_reader_num_row_groups(reader);
    if (row_group_index < 0 || row_group_index >= nrg) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "row_group_index out of range");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    const carquet_schema_t* cs = reader->schema;
    int32_t num_leaves = cs->num_leaves;

    rctx_t ctx = { cs, NULL, num_leaves, error };
    ctx.leaves = (rleaf_t*)calloc((size_t)(num_leaves > 0 ? num_leaves : 1), sizeof(rleaf_t));
    if (!ctx.leaves) { CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "alloc"); return CARQUET_ERROR_OUT_OF_MEMORY; }

    carquet_status_t rc = CARQUET_OK;

    /* Read every leaf column's full (rep, def, value) stream. */
    for (int32_t l = 0; l < num_leaves; l++) {
        rleaf_t* L = &ctx.leaves[l];
        int32_t elem = cs->leaf_indices[l];
        const parquet_schema_element_t* el = &cs->elements[elem];
        L->pt = el->type;
        L->type_length = el->type_length;
        L->max_def = cs->max_def_levels[l];
        L->max_rep = cs->max_rep_levels[l];
        L->is_bool = (el->type == CARQUET_PHYSICAL_BOOLEAN);
        L->is_bytearray = (el->type == CARQUET_PHYSICAL_BYTE_ARRAY);
        L->narrow = arrow_narrow_int_width(
            el->type, el->has_logical_type ? &el->logical_type : NULL);
        switch (el->type) {
        case CARQUET_PHYSICAL_INT32: case CARQUET_PHYSICAL_FLOAT: L->stride = 4; break;
        case CARQUET_PHYSICAL_INT64: case CARQUET_PHYSICAL_DOUBLE: L->stride = 8; break;
        case CARQUET_PHYSICAL_INT96: L->stride = 12; break;
        case CARQUET_PHYSICAL_BOOLEAN: L->stride = 1; break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: L->stride = (size_t)el->type_length; break;
        default: L->stride = 0; break;  /* byte array */
        }

        L->cr = carquet_reader_get_column(reader, row_group_index, l, error);
        if (!L->cr) { rc = CARQUET_ERROR_INTERNAL; goto cleanup; }
        int64_t total = carquet_column_remaining(L->cr);
        if (total < 0) { rc = CARQUET_ERROR_INTERNAL; goto cleanup; }
        int64_t alloc = total > 0 ? total : 1;

        /* Buffer stride must cover what carquet_read_next_page will write, which
         * uses the column reader's *own* physical type/length. On a malformed
         * file that can differ from the schema element type we shred against, so
         * size the value buffer to the larger of the two to stay in bounds
         * regardless (data for such files is undefined, but memory-safe). */
        size_t schema_vstride = L->is_bytearray ? sizeof(carquet_byte_array_t)
                              : (L->stride ? L->stride : 1);
        size_t reader_vstride;
        switch (L->cr->type) {
        case CARQUET_PHYSICAL_BOOLEAN: reader_vstride = 1; break;
        case CARQUET_PHYSICAL_INT32: case CARQUET_PHYSICAL_FLOAT: reader_vstride = 4; break;
        case CARQUET_PHYSICAL_INT64: case CARQUET_PHYSICAL_DOUBLE: reader_vstride = 8; break;
        case CARQUET_PHYSICAL_INT96: reader_vstride = 12; break;
        case CARQUET_PHYSICAL_BYTE_ARRAY: reader_vstride = sizeof(carquet_byte_array_t); break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            reader_vstride = L->cr->type_length > 0 ? (size_t)L->cr->type_length : 1; break;
        default: reader_vstride = sizeof(carquet_byte_array_t); break;
        }
        size_t vstride = schema_vstride > reader_vstride ? schema_vstride : reader_vstride;

        /* `total` is the column chunk's claimed num_values, taken from
         * attacker-controllable metadata. Bound every per-slot allocation
         * (def/rep levels and values) to a sane maximum so a malformed file
         * can't request a multi-terabyte buffer. Matches the batch reader's
         * CARQUET_MAX_BATCH_ALLOC guard. */
        #define CARQUET_ARROW_MAX_ALLOC (1024ULL * 1024 * 1024)
        size_t max_stride = vstride > sizeof(int16_t) ? vstride : sizeof(int16_t);
        if ((uint64_t)alloc > (uint64_t)(CARQUET_ARROW_MAX_ALLOC / max_stride)) {
            rc = CARQUET_ERROR_INVALID_ARGUMENT; goto cleanup;
        }
        #undef CARQUET_ARROW_MAX_ALLOC

        L->def = (int16_t*)malloc((size_t)alloc * sizeof(int16_t));
        L->rep = (int16_t*)malloc((size_t)alloc * sizeof(int16_t));
        L->values = malloc((size_t)alloc * vstride);
        if (!L->def || !L->rep || !L->values) { rc = CARQUET_ERROR_OUT_OF_MEMORY; goto cleanup; }

        int64_t ns = carquet_column_read_batch(L->cr, L->values, total, L->def, L->rep);
        if (ns < 0) { rc = CARQUET_ERROR_INTERNAL; goto cleanup; }
        L->nslots = ns;
        int64_t pres = 0;
        for (int64_t s = 0; s < ns; s++) if (L->def[s] == L->max_def) pres++;
        L->present = pres;
    }

    /* Assemble the top-level struct array. */
    {
        int32_t top[1024]; int32_t nt = elem_children(cs, 0, top, 1024);
        if (nt > 1024) { rc = CARQUET_ERROR_INVALID_ARGUMENT; goto cleanup; }
        int64_t num_rows = reader->metadata.row_groups[row_group_index].num_rows;

        memset(out_array, 0, sizeof(*out_array));
        out_array->length = num_rows;
        out_array->null_count = 0;
        out_array->offset = 0;
        out_array->n_buffers = 1;
        out_array->buffers = (const void**)calloc(1, sizeof(void*));  /* struct validity (absent) */
        out_array->n_children = nt;
        out_array->children = nt ? (struct ArrowArray**)calloc((size_t)nt, sizeof(void*)) : NULL;
        out_array->release = release_array;
        if (!out_array->buffers || (nt && !out_array->children)) {
            release_array(out_array); rc = CARQUET_ERROR_OUT_OF_MEMORY;
            CARQUET_SET_ERROR(error, rc, "alloc"); goto cleanup;
        }

        int32_t base = 0;
        for (int32_t i = 0; i < nt; i++) {
            rc = build_node(&ctx, top[i], 0, 0, 0, base, &out_array->children[i]);
            if (rc != CARQUET_OK) {
                release_array(out_array);
                CARQUET_SET_ERROR(error, rc, "Arrow read: failed to assemble field %d", i);
                goto cleanup;
            }
            int64_t got_len = out_array->children[i]->length;
            if (got_len != num_rows) {
                release_array(out_array);
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL,
                    "Arrow read: field %d length %lld != %lld rows",
                    i, (long long)got_len, (long long)num_rows);
                rc = CARQUET_ERROR_INTERNAL; goto cleanup;
            }
            base += count_leaves(cs, top[i]);
        }
    }

    if (out_schema) {
        rc = carquet_arrow_build_schema_tree(cs, out_schema, error);
        if (rc != CARQUET_OK) { release_array(out_array); goto cleanup; }
    }

cleanup:
    for (int32_t l = 0; l < num_leaves; l++) {
        free(ctx.leaves[l].def);
        free(ctx.leaves[l].rep);
        free(ctx.leaves[l].values);
        if (ctx.leaves[l].cr) carquet_column_reader_free(ctx.leaves[l].cr);
    }
    free(ctx.leaves);
    return rc;
}
