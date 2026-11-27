/**
 * @file arrow_schema.c
 * @brief Minimal backward FlatBuffer builder + Arrow IPC Schema message.
 *
 * Produces the base64-encoded encapsulated Arrow IPC Schema message stored
 * under the Parquet footer key "ARROW:schema". The FlatBuffer builder is a
 * faithful, compact port of the reference (Python) flatbuffers builder
 * algorithm: it writes back-to-front, references point forward, and tables
 * carry a soffset to a trailing vtable.
 */

#include <carquet/carquet.h>

#include "core/allocator.h"
#include "arrow_schema.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Arrow flatbuf constants (format/Message.fbs, Schema.fbs, Type.fbs) -- */
enum { ARROW_METADATA_V5 = 4 };
enum { MSG_HEADER_SCHEMA = 1 };
/* Type union tags (format/Type.fbs) */
enum {
    AT_NULL = 1, AT_INT = 2, AT_FLOAT = 3, AT_BINARY = 4, AT_UTF8 = 5,
    AT_BOOL = 6, AT_DECIMAL = 7, AT_DATE = 8, AT_TIME = 9, AT_TIMESTAMP = 10,
    AT_LIST = 12, AT_STRUCT = 13, AT_FIXEDSIZEBINARY = 15, AT_MAP = 17
};
/* Deepest schema nesting the builder will emit (guards runaway recursion). */
enum { ARROW_MAX_NEST_DEPTH = 100 };
/* TimeUnit: SEC=0 MILLI=1 MICRO=2 NANO=3 ; DateUnit DAY=0 ; Precision DOUBLE=2 SINGLE=1 */

/* ============================================================================
 * Backward FlatBuffer builder
 * ============================================================================
 */

typedef struct {
    uint8_t* buf;       /* buffer; live data is buf[head .. cap) */
    size_t cap;
    size_t head;        /* write cursor; prepending decreases head */
    size_t minalign;
    int32_t* vtable;    /* current object's field->Offset() slots */
    int vtable_len;
    size_t object_end;  /* Offset() at StartObject time */
    int oom;
} fbb;

static int fbb_init(fbb* b) {
    b->cap = 1024;
    b->buf = (uint8_t*)carquet_mem_malloc(b->cap);
    b->head = b->cap;
    b->minalign = 1;
    b->vtable = NULL;
    b->vtable_len = 0;
    b->object_end = 0;
    b->oom = (b->buf == NULL);
    return !b->oom;
}

static void fbb_free(fbb* b) {
    carquet_mem_free(b->buf);
    carquet_mem_free(b->vtable);
}

static size_t fbb_offset(const fbb* b) { return b->cap - b->head; }

/* Grow the buffer, keeping live bytes anchored at the high end. */
static int fbb_grow(fbb* b) {
    size_t old_cap = b->cap;
    size_t live = old_cap - b->head;
    size_t new_cap = old_cap * 2;
    uint8_t* nb = (uint8_t*)carquet_mem_malloc(new_cap);
    if (!nb) { b->oom = 1; return 0; }
    memcpy(nb + (new_cap - live), b->buf + b->head, live);
    carquet_mem_free(b->buf);
    b->buf = nb;
    b->cap = new_cap;
    b->head = new_cap - live;
    return 1;
}

static void fbb_pad(fbb* b, size_t n) {
    while (n--) {
        if (b->head == 0 && !fbb_grow(b)) return;
        b->buf[--b->head] = 0;
    }
}

static void fbb_prep(fbb* b, size_t size, size_t additional) {
    if (size > b->minalign) b->minalign = size;
    size_t buf_used = b->cap - b->head + additional;
    size_t align_size = ((~buf_used) + 1) & (size - 1);
    while (b->head < align_size + size + additional) {
        if (!fbb_grow(b)) return;
    }
    fbb_pad(b, align_size);
}

static void fbb_put(fbb* b, const void* data, size_t n) {
    while (b->head < n) { if (!fbb_grow(b)) return; }
    b->head -= n;
    memcpy(b->buf + b->head, data, n);
}

static void fbb_place_u8(fbb* b, uint8_t v)  { fbb_put(b, &v, 1); }
static void fbb_place_u16(fbb* b, uint16_t v) {
    uint8_t t[2] = { (uint8_t)(v & 0xFF), (uint8_t)(v >> 8) };
    fbb_put(b, t, 2);
}
static void fbb_place_u32(fbb* b, uint32_t v) {
    uint8_t t[4] = { (uint8_t)v, (uint8_t)(v >> 8),
                     (uint8_t)(v >> 16), (uint8_t)(v >> 24) };
    fbb_put(b, t, 4);
}

/* Prepend a forward uoffset that refers to a previously built object. */
static void fbb_prepend_uoffset(fbb* b, size_t off) {
    fbb_prep(b, 4, 0);
    uint32_t v = (uint32_t)(fbb_offset(b) - off + 4);
    fbb_place_u32(b, v);
}

/* ---- strings & vectors ---- */
static size_t fbb_create_string(fbb* b, const char* s) {
    size_t len = strlen(s);
    fbb_prep(b, 4, len + 1);
    fbb_place_u8(b, 0);                 /* NUL terminator (not counted) */
    fbb_put(b, s, len);
    fbb_place_u32(b, (uint32_t)len);    /* vector length prefix */
    return fbb_offset(b);
}

static void fbb_start_vector(fbb* b, size_t elem_size, size_t n, size_t align) {
    fbb_prep(b, 4, elem_size * n);
    fbb_prep(b, align, elem_size * n);
}
static size_t fbb_end_vector(fbb* b, size_t n) {
    fbb_place_u32(b, (uint32_t)n);
    return fbb_offset(b);
}

/* ---- tables ---- */
static int fbb_start_table(fbb* b, int numfields) {
    carquet_mem_free(b->vtable);
    /* At least one slot: a 0-byte allocation may return NULL (empty type tables). */
    b->vtable = (int32_t*)carquet_mem_calloc((size_t)(numfields > 0 ? numfields : 1), sizeof(int32_t));
    if (!b->vtable) { b->oom = 1; return 0; }
    b->vtable_len = numfields;
    b->object_end = fbb_offset(b);
    return 1;
}
static void fbb_slot(fbb* b, int slot) { b->vtable[slot] = (int32_t)fbb_offset(b); }

static void fbb_add_uoffset(fbb* b, int slot, size_t off) {
    fbb_prepend_uoffset(b, off);
    fbb_slot(b, slot);
}
static void fbb_add_u8(fbb* b, int slot, uint8_t v) {
    fbb_prep(b, 1, 0); fbb_place_u8(b, v); fbb_slot(b, slot);
}
static void fbb_add_i16(fbb* b, int slot, int16_t v) {
    fbb_prep(b, 2, 0); fbb_place_u16(b, (uint16_t)v); fbb_slot(b, slot);
}
static void fbb_add_i32(fbb* b, int slot, int32_t v) {
    fbb_prep(b, 4, 0); fbb_place_u32(b, (uint32_t)v); fbb_slot(b, slot);
}

static size_t fbb_end_table(fbb* b) {
    /* placeholder soffset (patched below) */
    fbb_prep(b, 4, 0);
    fbb_place_u32(b, 0);
    size_t object_offset = fbb_offset(b);

    int trimmed = b->vtable_len;
    while (trimmed > 0 && b->vtable[trimmed - 1] == 0) trimmed--;

    for (int i = trimmed - 1; i >= 0; i--) {
        uint16_t off = b->vtable[i]
            ? (uint16_t)(object_offset - (size_t)b->vtable[i]) : 0;
        fbb_prep(b, 2, 0); fbb_place_u16(b, off);
    }
    fbb_prep(b, 2, 0);
    fbb_place_u16(b, (uint16_t)(object_offset - b->object_end));   /* object size */
    fbb_prep(b, 2, 0);
    fbb_place_u16(b, (uint16_t)((trimmed + 2) * 2));               /* vtable size */

    /* patch soffset at the table start */
    int32_t soffset = (int32_t)(fbb_offset(b) - object_offset);
    size_t pos = b->cap - object_offset;
    b->buf[pos + 0] = (uint8_t)soffset;
    b->buf[pos + 1] = (uint8_t)(soffset >> 8);
    b->buf[pos + 2] = (uint8_t)(soffset >> 16);
    b->buf[pos + 3] = (uint8_t)(soffset >> 24);

    carquet_mem_free(b->vtable);
    b->vtable = NULL;
    b->vtable_len = 0;
    return object_offset;
}

static void fbb_finish(fbb* b, size_t root) {
    fbb_prep(b, b->minalign, 4);
    fbb_prepend_uoffset(b, root);
}

/* ============================================================================
 * Arrow type construction
 * ============================================================================
 */

/* Arrow Decimal type table. The precision/scale live in the DECIMAL logical
 * type; the legacy SchemaElement.precision/scale fields are only a fallback
 * (they are unset for schemas built from logical types, and a precision of 0
 * makes the whole ARROW:schema blob unreadable for Arrow C++ / PyArrow). */
static int build_decimal_type(fbb* b, const parquet_schema_element_t* e,
                              size_t* type_off, uint8_t* type_tag) {
    int32_t precision = e->logical_type.params.decimal.precision;
    int32_t scale = e->logical_type.params.decimal.scale;
    if (precision <= 0) { precision = e->precision; scale = e->scale; }
    if (precision <= 0 || precision > 76) return 0;  /* not Arrow-representable */
    fbb_start_table(b, 3);
    fbb_add_i32(b, 0, precision);
    fbb_add_i32(b, 1, scale);
    fbb_add_i32(b, 2, precision > 38 ? 256 : 128);
    *type_off = fbb_end_table(b); *type_tag = AT_DECIMAL; return 1;
}

/* Returns 0 on unsupported type (caller aborts ARROW:schema emission). */
static int build_arrow_type(fbb* b, const parquet_schema_element_t* e,
                            uint8_t* type_tag, size_t* type_off) {
    carquet_physical_type_t pt = e->type;
    int has_lt = e->has_logical_type;
    carquet_logical_type_id_t lt = has_lt ? e->logical_type.id : CARQUET_LOGICAL_UNKNOWN;
    carquet_converted_type_t ct = e->has_converted_type ? e->converted_type
                                                        : CARQUET_CONVERTED_NONE;

    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN:
        fbb_start_table(b, 0); *type_off = fbb_end_table(b);
        *type_tag = AT_BOOL; return 1;

    case CARQUET_PHYSICAL_FLOAT:
        fbb_start_table(b, 1); fbb_add_i16(b, 0, 1 /*SINGLE*/);
        *type_off = fbb_end_table(b); *type_tag = AT_FLOAT; return 1;

    case CARQUET_PHYSICAL_DOUBLE:
        fbb_start_table(b, 1); fbb_add_i16(b, 0, 2 /*DOUBLE*/);
        *type_off = fbb_end_table(b); *type_tag = AT_FLOAT; return 1;

    case CARQUET_PHYSICAL_INT32:
    case CARQUET_PHYSICAL_INT64: {
        int is64 = (pt == CARQUET_PHYSICAL_INT64);
        if (has_lt && lt == CARQUET_LOGICAL_TIMESTAMP && is64) {
            int16_t unit = (e->logical_type.params.timestamp.unit ==
                            CARQUET_TIME_UNIT_MILLIS) ? 1 :
                           (e->logical_type.params.timestamp.unit ==
                            CARQUET_TIME_UNIT_MICROS) ? 2 : 3;
            size_t tz = 0;
            int adj = e->logical_type.params.timestamp.is_adjusted_to_utc;
            if (adj) tz = fbb_create_string(b, "UTC");
            fbb_start_table(b, 2);
            fbb_add_i16(b, 0, unit);
            if (adj) fbb_add_uoffset(b, 1, tz);
            *type_off = fbb_end_table(b); *type_tag = AT_TIMESTAMP; return 1;
        }
        if (has_lt && lt == CARQUET_LOGICAL_DATE && !is64) {
            fbb_start_table(b, 1); /* DateUnit DAY=0 default */
            *type_off = fbb_end_table(b); *type_tag = AT_DATE; return 1;
        }
        if (has_lt && lt == CARQUET_LOGICAL_TIME) {
            int16_t unit = (e->logical_type.params.time.unit ==
                            CARQUET_TIME_UNIT_MILLIS) ? 1 :
                           (e->logical_type.params.time.unit ==
                            CARQUET_TIME_UNIT_MICROS) ? 2 : 3;
            fbb_start_table(b, 2);
            fbb_add_i16(b, 0, unit);
            fbb_add_i32(b, 1, is64 ? 64 : 32);
            *type_off = fbb_end_table(b); *type_tag = AT_TIME; return 1;
        }
        if (has_lt && lt == CARQUET_LOGICAL_DECIMAL) {
            return build_decimal_type(b, e, type_off, type_tag);
        }
        {
            int32_t bw = is64 ? 64 : 32;
            int is_signed = 1;
            if (has_lt && lt == CARQUET_LOGICAL_INTEGER) {
                bw = e->logical_type.params.integer.bit_width;
                is_signed = e->logical_type.params.integer.is_signed;
            } else if (ct >= CARQUET_CONVERTED_UINT_8 &&
                       ct <= CARQUET_CONVERTED_INT_64) {
                static const int8_t cbw[] = { 8, 16, 32, 64, 8, 16, 32, 64 };
                bw = cbw[ct - CARQUET_CONVERTED_UINT_8];
                is_signed = (ct >= CARQUET_CONVERTED_INT_8);
            }
            fbb_start_table(b, 2);
            fbb_add_i32(b, 0, bw);
            if (is_signed) fbb_add_u8(b, 1, 1);
            *type_off = fbb_end_table(b); *type_tag = AT_INT; return 1;
        }
    }

    case CARQUET_PHYSICAL_INT96:
        /* Arrow maps the deprecated INT96 to timestamp[ns] (no tz). */
        fbb_start_table(b, 2);
        fbb_add_i16(b, 0, 3 /*NANO*/);
        *type_off = fbb_end_table(b); *type_tag = AT_TIMESTAMP; return 1;

    case CARQUET_PHYSICAL_BYTE_ARRAY:
        if ((has_lt && (lt == CARQUET_LOGICAL_STRING ||
                        lt == CARQUET_LOGICAL_JSON ||
                        lt == CARQUET_LOGICAL_ENUM)) ||
            ct == CARQUET_CONVERTED_UTF8 || ct == CARQUET_CONVERTED_JSON ||
            ct == CARQUET_CONVERTED_ENUM) {
            fbb_start_table(b, 0); *type_off = fbb_end_table(b);
            *type_tag = AT_UTF8; return 1;
        }
        if (has_lt && lt == CARQUET_LOGICAL_DECIMAL) return 0; /* not Arrow-representable */
        fbb_start_table(b, 0); *type_off = fbb_end_table(b);
        *type_tag = AT_BINARY; return 1;

    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
        if (has_lt && lt == CARQUET_LOGICAL_DECIMAL) {
            return build_decimal_type(b, e, type_off, type_tag);
        }
        fbb_start_table(b, 1);
        fbb_add_i32(b, 0, e->type_length);
        *type_off = fbb_end_table(b); *type_tag = AT_FIXEDSIZEBINARY; return 1;

    default:
        return 0;
    }
}

/* Build the Arrow Field.custom_metadata [KeyValue] vector for one field.
 * Returns the vector offset, or 0 when the field has no metadata (the caller
 * then leaves slot 6 unset). Sets b->oom on allocation failure. */
static size_t build_field_metadata(fbb* b, const parquet_key_value_t* kv,
                                   int32_t n) {
    if (n <= 0 || !kv) return 0;
    size_t* kv_offs = (size_t*)carquet_mem_malloc((size_t)n * sizeof(size_t));
    if (!kv_offs) { b->oom = 1; return 0; }
    for (int32_t i = 0; i < n; i++) {
        /* Strings/tables are built before the vector that references them. */
        size_t key_off = fbb_create_string(b, kv[i].key ? kv[i].key : "");
        size_t val_off = kv[i].value ? fbb_create_string(b, kv[i].value) : 0;
        fbb_start_table(b, 2);              /* KeyValue { key, value } */
        fbb_add_uoffset(b, 0, key_off);
        if (val_off) fbb_add_uoffset(b, 1, val_off);
        kv_offs[i] = fbb_end_table(b);
    }
    fbb_start_vector(b, 4, (size_t)n, 4);
    for (int32_t i = n - 1; i >= 0; i--) fbb_prepend_uoffset(b, kv_offs[i]);
    size_t vec = fbb_end_vector(b, (size_t)n);
    carquet_mem_free(kv_offs);
    return vec;
}

/* ============================================================================
 * Recursive Field emission (nested LIST / MAP / STRUCT)
 * ============================================================================
 *
 * The Parquet schema arrives flattened in depth-first pre-order; the tree is
 * reconstructed from each element's num_children (via a precomputed parent[]
 * map). We walk it and emit Arrow Fields, collapsing the intermediate levels
 * that Arrow does not model: a Parquet 3-level LIST becomes Arrow List<child>,
 * and a Parquet MAP becomes Arrow Map<entries:struct<key,value>>.
 */

/* Fill out[] with the direct children of `node` (ascending index == Parquet
 * declaration order == Arrow field order). Returns the child count. */
static int32_t collect_children(const int32_t* parent, int32_t num_elements,
                                int32_t node, int32_t* out, int32_t max_out) {
    int32_t n = 0;
    for (int32_t j = 0; j < num_elements; j++) {
        if (parent[j] == node) {
            if (n < max_out) out[n] = j;
            n++;
        }
    }
    return n;
}

/* Emit the trailing slots common to every group Field: the (empty or single-
 * slot) type table plus name and custom_metadata, then the Field table with
 * its children vector. Assumes the children vector `kids_vec` is already
 * built. `type_ntab_slots` is the number of slots in the union type table. */
static size_t emit_group_field(fbb* b, const parquet_schema_element_t* e,
                               uint8_t type_tag, int type_tab_slots,
                               size_t kids_vec, int nullable, int* ok) {
    fbb_start_table(b, type_tab_slots);
    size_t type_off = fbb_end_table(b);  /* List/Struct/Map: no set fields */
    size_t name_off = e->name ? fbb_create_string(b, e->name) : 0;
    size_t meta_off = build_field_metadata(b, e->field_metadata,
                                           e->num_field_metadata);
    if (b->oom) { *ok = 0; return 0; }

    fbb_start_table(b, 7);
    if (name_off) fbb_add_uoffset(b, 0, name_off);
    if (nullable) fbb_add_u8(b, 1, 1);
    fbb_add_u8(b, 2, type_tag);
    fbb_add_uoffset(b, 3, type_off);
    fbb_add_uoffset(b, 5, kids_vec);          /* Field.children */
    if (meta_off) fbb_add_uoffset(b, 6, meta_off);
    return fbb_end_table(b);
}

/* Build one Arrow Field for schema element `idx` and return its offset.
 * `force_nonnull` overrides nullability (Arrow map keys / entries structs must
 * be non-nullable). Sets *ok = 0 and returns 0 on any unsupported shape. */
static size_t build_field(fbb* b, const parquet_schema_element_t* schema,
                          int32_t num_elements, const int32_t* parent,
                          int32_t idx, int force_nonnull, int depth, int* ok) {
    if (depth > ARROW_MAX_NEST_DEPTH) { *ok = 0; return 0; }
    const parquet_schema_element_t* e = &schema[idx];

    /* collect_children() returns the full child count but stores at most
     * max_out indexes: take the first child here (lists and maps need only
     * that one) and the full list in the struct branch, sized from nk. */
    int32_t first_kid = -1;
    int32_t nk = collect_children(parent, num_elements, idx, &first_kid, 1);
    int nullable = force_nonnull ? 0
        : !(e->has_repetition && e->repetition_type == CARQUET_REPETITION_REQUIRED);

    if (nk == 0) {
        /* Primitive leaf — identical to the historical flat emission. */
        uint8_t tag = 0; size_t toff = 0;
        if (!build_arrow_type(b, e, &tag, &toff)) { *ok = 0; return 0; }
        size_t name_off = e->name ? fbb_create_string(b, e->name) : 0;
        size_t meta_off = build_field_metadata(b, e->field_metadata,
                                               e->num_field_metadata);
        if (b->oom) { *ok = 0; return 0; }
        fbb_start_table(b, 7);
        if (name_off) fbb_add_uoffset(b, 0, name_off);
        if (nullable) fbb_add_u8(b, 1, 1);
        fbb_add_u8(b, 2, tag);
        fbb_add_uoffset(b, 3, toff);
        if (meta_off) fbb_add_uoffset(b, 6, meta_off);
        return fbb_end_table(b);
    }

    int is_list = (e->has_logical_type && e->logical_type.id == CARQUET_LOGICAL_LIST) ||
                  (e->has_converted_type && e->converted_type == CARQUET_CONVERTED_LIST);
    int is_map = (e->has_logical_type && e->logical_type.id == CARQUET_LOGICAL_MAP) ||
                 (e->has_converted_type &&
                  (e->converted_type == CARQUET_CONVERTED_MAP ||
                   e->converted_type == CARQUET_CONVERTED_MAP_KEY_VALUE));

    if (is_list) {
        /* group(LIST) → repeated group → element. Arrow List<element>. */
        if (nk != 1) { *ok = 0; return 0; }
        int32_t elem_kids[4];
        int32_t nek = collect_children(parent, num_elements, first_kid, elem_kids, 4);
        if (nek != 1) { *ok = 0; return 0; }
        size_t elem_off = build_field(b, schema, num_elements, parent,
                                      elem_kids[0], 0, depth + 1, ok);
        if (!*ok) return 0;
        fbb_start_vector(b, 4, 1, 4);
        fbb_prepend_uoffset(b, elem_off);
        size_t kids_vec = fbb_end_vector(b, 1);
        return emit_group_field(b, e, AT_LIST, 0, kids_vec, nullable, ok);
    }

    if (is_map) {
        /* group(MAP) → repeated key_value → {key, value?}.
         * Arrow Map<entries: struct<key, value>> (entries non-nullable). */
        if (nk != 1) { *ok = 0; return 0; }
        int32_t kv = first_kid;
        int32_t kv_kids[4];
        int32_t nkv = collect_children(parent, num_elements, kv, kv_kids, 4);
        if (nkv < 1 || nkv > 2) { *ok = 0; return 0; }
        size_t key_off = build_field(b, schema, num_elements, parent,
                                     kv_kids[0], 1, depth + 1, ok);
        if (!*ok) return 0;
        size_t val_off = 0;
        if (nkv == 2) {
            val_off = build_field(b, schema, num_elements, parent,
                                  kv_kids[1], 0, depth + 1, ok);
            if (!*ok) return 0;
        }
        /* entries struct children vector [key(, value)] */
        fbb_start_vector(b, 4, (size_t)nkv, 4);
        if (nkv == 2) fbb_prepend_uoffset(b, val_off);
        fbb_prepend_uoffset(b, key_off);
        size_t struct_kids = fbb_end_vector(b, (size_t)nkv);
        /* entries struct field (non-nullable), named after the key_value group */
        const parquet_schema_element_t* kve = &schema[kv];
        fbb_start_table(b, 0);
        size_t struct_type = fbb_end_table(b);
        size_t entries_name = fbb_create_string(b, kve->name ? kve->name : "key_value");
        if (b->oom) { *ok = 0; return 0; }
        fbb_start_table(b, 7);
        fbb_add_uoffset(b, 0, entries_name);
        fbb_add_u8(b, 2, AT_STRUCT);
        fbb_add_uoffset(b, 3, struct_type);
        fbb_add_uoffset(b, 5, struct_kids);
        size_t entries_off = fbb_end_table(b);
        /* Map field children vector [entries]. */
        fbb_start_vector(b, 4, 1, 4);
        fbb_prepend_uoffset(b, entries_off);
        size_t map_kids = fbb_end_vector(b, 1);
        /* Map type table carries a keysSorted bool (slot 0, default false). */
        return emit_group_field(b, e, AT_MAP, 1, map_kids, nullable, ok);
    }

    /* Plain STRUCT group. */
    {
        size_t* coffs = (size_t*)carquet_mem_malloc((size_t)nk * sizeof(size_t));
        int32_t* kids = (int32_t*)carquet_mem_malloc((size_t)nk * sizeof(int32_t));
        if (!coffs || !kids) {
            carquet_mem_free(coffs);
            carquet_mem_free(kids);
            b->oom = 1; *ok = 0; return 0;
        }
        collect_children(parent, num_elements, idx, kids, nk);
        for (int32_t i = 0; i < nk; i++) {
            coffs[i] = build_field(b, schema, num_elements, parent,
                                   kids[i], 0, depth + 1, ok);
            if (!*ok) { carquet_mem_free(coffs); carquet_mem_free(kids); return 0; }
        }
        carquet_mem_free(kids);
        fbb_start_vector(b, 4, (size_t)nk, 4);
        for (int32_t i = nk - 1; i >= 0; i--) fbb_prepend_uoffset(b, coffs[i]);
        size_t kids_vec = fbb_end_vector(b, (size_t)nk);
        carquet_mem_free(coffs);
        return emit_group_field(b, e, AT_STRUCT, 0, kids_vec, nullable, ok);
    }
}

/* Compute parent[i] for every element from the flattened num_children layout
 * (depth-first pre-order). Returns 0 on allocation failure. */
static int compute_parents(const parquet_schema_element_t* schema,
                           int32_t num_elements, int32_t* parent) {
    int32_t* st_idx = (int32_t*)carquet_mem_malloc((size_t)num_elements * sizeof(int32_t));
    int32_t* st_rem = (int32_t*)carquet_mem_malloc((size_t)num_elements * sizeof(int32_t));
    if (!st_idx || !st_rem) { carquet_mem_free(st_idx); carquet_mem_free(st_rem); return 0; }
    int sp = 0;
    for (int32_t i = 0; i < num_elements; i++) {
        while (sp > 0 && st_rem[sp - 1] == 0) sp--;
        if (sp > 0) { parent[i] = st_idx[sp - 1]; st_rem[sp - 1]--; }
        else parent[i] = -1;
        if (schema[i].num_children > 0) {
            st_idx[sp] = i; st_rem[sp] = schema[i].num_children; sp++;
        }
    }
    carquet_mem_free(st_idx);
    carquet_mem_free(st_rem);
    return 1;
}

/* ============================================================================
 * base64
 * ============================================================================
 */

static char* base64(const uint8_t* in, size_t n) {
    static const char T[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    size_t out_len = ((n + 2) / 3) * 4;
    char* out = (char*)carquet_mem_malloc(out_len + 1);
    if (!out) return NULL;
    size_t o = 0;
    for (size_t i = 0; i < n; i += 3) {
        uint32_t v = (uint32_t)in[i] << 16;
        if (i + 1 < n) v |= (uint32_t)in[i + 1] << 8;
        if (i + 2 < n) v |= in[i + 2];
        out[o++] = T[(v >> 18) & 0x3F];
        out[o++] = T[(v >> 12) & 0x3F];
        out[o++] = (i + 1 < n) ? T[(v >> 6) & 0x3F] : '=';
        out[o++] = (i + 2 < n) ? T[v & 0x3F] : '=';
    }
    out[o] = '\0';
    return out;
}

/* Finish the Schema message around `nfields` built Field offsets (and an
 * optional Schema.custom_metadata vector), encapsulate it as an IPC message
 * and base64 it. Frees the builder. */
static char* finish_schema_message(fbb* b, const size_t* field_offs, int32_t nfields, size_t meta_vec) {
    fbb_start_vector(b, 4, (size_t)nfields, 4);
    for (int32_t i = nfields - 1; i >= 0; i--) fbb_prepend_uoffset(b, field_offs[i]);
    size_t fields_vec = fbb_end_vector(b, (size_t)nfields);

    /* Schema table: only set fields (endianness Little=0 is the default). */
    fbb_start_table(b, 4);
    fbb_add_uoffset(b, 1, fields_vec);
    if (meta_vec) fbb_add_uoffset(b, 2, meta_vec);
    size_t schema_off = fbb_end_table(b);

    /* Message table */
    fbb_start_table(b, 5);
    fbb_add_i16(b, 0, ARROW_METADATA_V5);
    fbb_add_u8(b, 1, MSG_HEADER_SCHEMA);
    fbb_add_uoffset(b, 2, schema_off);
    size_t msg_off = fbb_end_table(b);

    fbb_finish(b, msg_off);
    if (b->oom) { fbb_free(b); return NULL; }

    size_t fb_len = fbb_offset(b);
    const uint8_t* fb = b->buf + b->head;

    /* Encapsulated IPC message: 0xFFFFFFFF | int32 metadata_len | flatbuffer,
     * with the flatbuffer padded so metadata_len is a multiple of 8. */
    size_t padded = (fb_len + 7) & ~(size_t)7;
    size_t enc_len = 8 + padded;
    uint8_t* enc = (uint8_t*)carquet_mem_malloc(enc_len);
    if (!enc) { fbb_free(b); return NULL; }
    enc[0] = enc[1] = enc[2] = enc[3] = 0xFF;
    uint32_t m = (uint32_t)padded;
    enc[4] = (uint8_t)m; enc[5] = (uint8_t)(m >> 8);
    enc[6] = (uint8_t)(m >> 16); enc[7] = (uint8_t)(m >> 24);
    memcpy(enc + 8, fb, fb_len);
    memset(enc + 8 + fb_len, 0, padded - fb_len);
    fbb_free(b);

    char* b64 = base64(enc, enc_len);
    carquet_mem_free(enc);
    return b64;
}

/* ============================================================================
 * Public entry
 * ============================================================================
 */

char* carquet_build_arrow_schema_b64(
    const parquet_schema_element_t* schema,
    int32_t num_elements) {

    if (!schema || num_elements < 2) return NULL;

    /* Reconstruct the tree (parent per element) from the flattened layout. */
    int32_t* parent = (int32_t*)carquet_mem_malloc((size_t)num_elements * sizeof(int32_t));
    if (!parent) return NULL;
    if (!compute_parents(schema, num_elements, parent)) {
        carquet_mem_free(parent);
        return NULL;
    }

    /* Top-level fields are the root's direct children. */
    int32_t* top = (int32_t*)carquet_mem_malloc((size_t)num_elements * sizeof(int32_t));
    if (!top) { carquet_mem_free(parent); return NULL; }
    int32_t nfields = collect_children(parent, num_elements, 0, top, num_elements);
    if (nfields < 1) { carquet_mem_free(top); carquet_mem_free(parent); return NULL; }

    fbb b;
    if (!fbb_init(&b)) { carquet_mem_free(top); carquet_mem_free(parent); return NULL; }

    /* Build each top-level Field (recursively descending into nested groups).
     * Any element we cannot faithfully map to Arrow aborts the whole emission
     * (return NULL) so we never write a schema that disagrees with Parquet. */
    size_t* field_offs = (size_t*)carquet_mem_malloc((size_t)nfields * sizeof(size_t));
    if (!field_offs) { fbb_free(&b); carquet_mem_free(top); carquet_mem_free(parent); return NULL; }

    int ok = 1;
    for (int32_t i = 0; i < nfields && ok; i++) {
        field_offs[i] = build_field(&b, schema, num_elements, parent,
                                    top[i], 0, 0, &ok);
    }
    carquet_mem_free(top);
    carquet_mem_free(parent);
    if (!ok || b.oom) { carquet_mem_free(field_offs); fbb_free(&b); return NULL; }

    char* out = finish_schema_message(&b, field_offs, nfields, 0);
    carquet_mem_free(field_offs);
    return out;
}

/* ============================================================================
 * From an Arrow schema (the original types of data written via write_arrow)
 * ============================================================================
 */

enum {
    AT_INTERVAL = 11, AT_FIXEDSIZELIST = 16, AT_DURATION = 18, AT_LARGEBINARY_ = 19,
    AT_LARGEUTF8_ = 20, AT_LARGELIST_ = 21, AT_BINARYVIEW = 23, AT_UTF8VIEW = 24
};

static int16_t unit_code(char c) {
    return c == 's' ? 0 : c == 'm' ? 1 : c == 'u' ? 2 : c == 'n' ? 3 : -1;
}

/* Arrow metadata blob -> [KeyValue] vector offset (0 when none). */
static size_t metadata_vector(fbb* b, const char* md) {
    if (!md) return 0;
    int32_t n;
    memcpy(&n, md, 4);
    if (n <= 0 || n > 1 << 16) return 0;
    size_t* offs = (size_t*)carquet_mem_malloc((size_t)n * sizeof(size_t));
    if (!offs) { b->oom = 1; return 0; }
    size_t pos = 4;
    for (int32_t i = 0; i < n; i++) {
        int32_t kl, vl;
        memcpy(&kl, md + pos, 4);
        char* k = (char*)carquet_mem_malloc((size_t)kl + 1);
        if (!k) { b->oom = 1; carquet_mem_free(offs); return 0; }
        memcpy(k, md + pos + 4, (size_t)kl); k[kl] = '\0';
        pos += 4 + (size_t)kl;
        memcpy(&vl, md + pos, 4);
        char* v = (char*)carquet_mem_malloc((size_t)vl + 1);
        if (!v) { b->oom = 1; carquet_mem_free(k); carquet_mem_free(offs); return 0; }
        memcpy(v, md + pos + 4, (size_t)vl); v[vl] = '\0';
        pos += 4 + (size_t)vl;
        size_t ko = fbb_create_string(b, k);
        size_t vo = fbb_create_string(b, v);
        carquet_mem_free(k);
        carquet_mem_free(v);
        fbb_start_table(b, 2);
        fbb_add_uoffset(b, 0, ko);
        fbb_add_uoffset(b, 1, vo);
        offs[i] = fbb_end_table(b);
    }
    fbb_start_vector(b, 4, (size_t)n, 4);
    for (int32_t i = n - 1; i >= 0; i--) fbb_prepend_uoffset(b, offs[i]);
    size_t vec = fbb_end_vector(b, (size_t)n);
    carquet_mem_free(offs);
    return vec;
}

static int int_bits(const char* f, int32_t* bits, int* is_signed) {
    if (f[1] != '\0') return 0;
    switch (f[0]) {
    case 'c': *bits = 8;  *is_signed = 1; return 1;
    case 'C': *bits = 8;  *is_signed = 0; return 1;
    case 's': *bits = 16; *is_signed = 1; return 1;
    case 'S': *bits = 16; *is_signed = 0; return 1;
    case 'i': *bits = 32; *is_signed = 1; return 1;
    case 'I': *bits = 32; *is_signed = 0; return 1;
    case 'l': *bits = 64; *is_signed = 1; return 1;
    case 'L': *bits = 64; *is_signed = 0; return 1;
    default: return 0;
    }
}

/* Type table for an Arrow format string; 0 when it has no IPC equivalent. */
static int arrow_format_type(fbb* b, const char* f, int keys_sorted, uint8_t* tag, size_t* off) {
    int32_t bits; int sg;
    if (int_bits(f, &bits, &sg)) {
        fbb_start_table(b, 2);
        fbb_add_i32(b, 0, bits);
        if (sg) fbb_add_u8(b, 1, 1);
        *off = fbb_end_table(b); *tag = AT_INT; return 1;
    }
    if (f[1] == '\0') {
        int16_t prec = -1;
        switch (f[0]) {
        case 'n': *tag = AT_NULL; break;
        case 'b': *tag = AT_BOOL; break;
        case 'z': *tag = AT_BINARY; break;
        case 'u': *tag = AT_UTF8; break;
        case 'Z': *tag = AT_LARGEBINARY_; break;
        case 'U': *tag = AT_LARGEUTF8_; break;
        case 'e': prec = 0; break;
        case 'f': prec = 1; break;
        case 'g': prec = 2; break;
        default: return 0;
        }
        if (prec >= 0) {
            fbb_start_table(b, 1);
            fbb_add_i16(b, 0, prec);
            *off = fbb_end_table(b); *tag = AT_FLOAT; return 1;
        }
        fbb_start_table(b, 0);
        *off = fbb_end_table(b); return 1;
    }
    if (strcmp(f, "vz") == 0 || strcmp(f, "vu") == 0) {
        fbb_start_table(b, 0);
        *off = fbb_end_table(b); *tag = f[1] == 'z' ? AT_BINARYVIEW : AT_UTF8VIEW; return 1;
    }
    if (f[0] == 'd' && f[1] == ':') {
        int p = 0, sc = 0, bw = 128;
        int n = sscanf(f + 2, "%d,%d,%d", &p, &sc, &bw);
        if (n < 2 || p <= 0) return 0;
        fbb_start_table(b, 3);
        fbb_add_i32(b, 0, p);
        fbb_add_i32(b, 1, sc);
        fbb_add_i32(b, 2, bw);
        *off = fbb_end_table(b); *tag = AT_DECIMAL; return 1;
    }
    if (f[0] == 'w' && f[1] == ':') {
        int w = atoi(f + 2);
        if (w <= 0) return 0;
        fbb_start_table(b, 1);
        fbb_add_i32(b, 0, w);
        *off = fbb_end_table(b); *tag = AT_FIXEDSIZEBINARY; return 1;
    }
    if (f[0] == 't') {
        if (f[1] == 'd' && (f[2] == 'D' || f[2] == 'm') && f[3] == '\0') {
            fbb_start_table(b, 1);
            fbb_add_i16(b, 0, f[2] == 'D' ? 0 : 1);
            *off = fbb_end_table(b); *tag = AT_DATE; return 1;
        }
        if (f[1] == 't' && unit_code(f[2]) >= 0 && f[3] == '\0') {
            fbb_start_table(b, 2);
            fbb_add_i16(b, 0, unit_code(f[2]));
            fbb_add_i32(b, 1, (f[2] == 's' || f[2] == 'm') ? 32 : 64);
            *off = fbb_end_table(b); *tag = AT_TIME; return 1;
        }
        if (f[1] == 's' && unit_code(f[2]) >= 0 && f[3] == ':') {
            size_t tz = f[4] ? fbb_create_string(b, f + 4) : 0;
            fbb_start_table(b, 2);
            fbb_add_i16(b, 0, unit_code(f[2]));
            if (tz) fbb_add_uoffset(b, 1, tz);
            *off = fbb_end_table(b); *tag = AT_TIMESTAMP; return 1;
        }
        if (f[1] == 'D' && unit_code(f[2]) >= 0 && f[3] == '\0') {
            fbb_start_table(b, 1);
            fbb_add_i16(b, 0, unit_code(f[2]));
            *off = fbb_end_table(b); *tag = AT_DURATION; return 1;
        }
        if (f[1] == 'i' && (f[2] == 'M' || f[2] == 'D' || f[2] == 'n') && f[3] == '\0') {
            fbb_start_table(b, 1);
            fbb_add_i16(b, 0, f[2] == 'M' ? 0 : f[2] == 'D' ? 1 : 2);
            *off = fbb_end_table(b); *tag = AT_INTERVAL; return 1;
        }
        return 0;
    }
    if (f[0] == '+') {
        if (strcmp(f, "+l") == 0 || strcmp(f, "+L") == 0 || strcmp(f, "+s") == 0) {
            fbb_start_table(b, 0);
            *off = fbb_end_table(b);
            *tag = f[1] == 'l' ? AT_LIST : f[1] == 'L' ? AT_LARGELIST_ : AT_STRUCT;
            return 1;
        }
        if (strcmp(f, "+m") == 0) {
            fbb_start_table(b, 1);
            if (keys_sorted) fbb_add_u8(b, 0, 1);
            *off = fbb_end_table(b); *tag = AT_MAP; return 1;
        }
        if (f[1] == 'w' && f[2] == ':' && atoi(f + 3) > 0) {
            fbb_start_table(b, 1);
            fbb_add_i32(b, 0, atoi(f + 3));
            *off = fbb_end_table(b); *tag = AT_FIXEDSIZELIST; return 1;
        }
    }
    return 0;
}

static size_t arrow_field(fbb* b, const struct ArrowSchema* s, int depth, int64_t* dict_id, int* ok) {
    if (depth > ARROW_MAX_NEST_DEPTH || !s->format || s->n_children < 0) { *ok = 0; return 0; }
    /* A dictionary-encoded field has the value type; its own format is the index type. */
    const struct ArrowSchema* vt = s->dictionary ? s->dictionary : s;
    size_t* kids = NULL;
    if (vt->n_children > 0) {
        kids = (size_t*)carquet_mem_malloc((size_t)vt->n_children * sizeof(size_t));
        if (!kids) { b->oom = 1; *ok = 0; return 0; }
        for (int64_t i = 0; i < vt->n_children; i++) {
            kids[i] = arrow_field(b, vt->children[i], depth + 1, dict_id, ok);
            if (!*ok) { carquet_mem_free(kids); return 0; }
        }
    }
    size_t kids_vec = 0;
    if (vt->n_children > 0) {
        fbb_start_vector(b, 4, (size_t)vt->n_children, 4);
        for (int64_t i = vt->n_children - 1; i >= 0; i--) fbb_prepend_uoffset(b, kids[i]);
        kids_vec = fbb_end_vector(b, (size_t)vt->n_children);
    }
    carquet_mem_free(kids);

    uint8_t tag = 0;
    size_t type_off = 0;
    if (!arrow_format_type(b, vt->format, (s->flags & ARROW_FLAG_MAP_KEYS_SORTED) != 0, &tag, &type_off)) {
        *ok = 0;
        return 0;
    }
    size_t dict_off = 0;
    if (s->dictionary) {
        int32_t bits; int sg;
        if (!int_bits(s->format, &bits, &sg)) { *ok = 0; return 0; }
        fbb_start_table(b, 2);
        fbb_add_i32(b, 0, bits);
        if (sg) fbb_add_u8(b, 1, 1);
        size_t idx = fbb_end_table(b);
        fbb_start_table(b, 4);
        int64_t id = (*dict_id)++;
        /* DictionaryEncoding.id is a long: two little-endian words. */
        fbb_prep(b, 8, 0);
        fbb_place_u32(b, (uint32_t)((uint64_t)id >> 32));
        fbb_place_u32(b, (uint32_t)id);
        fbb_slot(b, 0);
        fbb_add_uoffset(b, 1, idx);
        if (s->flags & ARROW_FLAG_DICTIONARY_ORDERED) fbb_add_u8(b, 2, 1);
        dict_off = fbb_end_table(b);
    }
    size_t name_off = fbb_create_string(b, s->name ? s->name : "");
    size_t meta_off = metadata_vector(b, s->metadata);
    if (b->oom) { *ok = 0; return 0; }
    fbb_start_table(b, 7);
    fbb_add_uoffset(b, 0, name_off);
    if (s->flags & ARROW_FLAG_NULLABLE) fbb_add_u8(b, 1, 1);
    fbb_add_u8(b, 2, tag);
    fbb_add_uoffset(b, 3, type_off);
    if (dict_off) fbb_add_uoffset(b, 4, dict_off);
    if (kids_vec) fbb_add_uoffset(b, 5, kids_vec);
    if (meta_off) fbb_add_uoffset(b, 6, meta_off);
    return fbb_end_table(b);
}

char* carquet_build_arrow_schema_b64_from_arrow(const struct ArrowSchema* schema) {
    if (!schema || !schema->format || strcmp(schema->format, "+s") != 0 || schema->n_children < 1) return NULL;
    fbb b;
    if (!fbb_init(&b)) return NULL;
    int32_t n = (int32_t)schema->n_children;
    size_t* offs = (size_t*)carquet_mem_malloc((size_t)n * sizeof(size_t));
    if (!offs) { fbb_free(&b); return NULL; }
    int ok = 1;
    int64_t dict_id = 0;
    for (int32_t i = 0; i < n && ok; i++) offs[i] = arrow_field(&b, schema->children[i], 0, &dict_id, &ok);
    size_t meta = ok ? metadata_vector(&b, schema->metadata) : 0;
    if (!ok || b.oom) { carquet_mem_free(offs); fbb_free(&b); return NULL; }
    char* out = finish_schema_message(&b, offs, n, meta);
    carquet_mem_free(offs);
    return out;
}
