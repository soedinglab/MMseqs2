/**
 * @file arrow_schema_read.c
 * @brief Minimal bounds-checked FlatBuffer reader for the "ARROW:schema" blob.
 *
 * Navigates Message -> Schema -> [Field] -> Field.custom_metadata -> [KeyValue]
 * and attaches each field's metadata to the matching Parquet schema element.
 * The input is untrusted (it comes from the file), so every offset, length and
 * vector span is validated before use; anything inconsistent aborts the parse
 * of that sub-tree without touching the rest.
 */

#include "arrow_schema_read.h"
#include "core/allocator.h"
#include <stdlib.h>
#include <string.h>

/* Arrow FlatBuffer field ids (format/Message.fbs, Schema.fbs). A union in the
 * schema occupies two vtable slots (type, then value), which is why the Schema
 * table sits at Message slot 2. */
enum { MSG_HEADER_SLOT = 2 };      /* Message.header (union value) */
enum { SCHEMA_FIELDS_SLOT = 1 };   /* Schema.fields */
enum { FIELD_NAME_SLOT = 0 };      /* Field.name */
enum { FIELD_TYPETYPE_SLOT = 2 };  /* Field.type_type (union tag, u8) */
enum { FIELD_META_SLOT = 6 };      /* Field.custom_metadata */
enum { KV_KEY_SLOT = 0, KV_VALUE_SLOT = 1 };

/* Arrow Type union tags for the 64-bit-offset variants Parquet can't express
 * (format/Type.fbs). Kept in sync with carquet_arrow_type_refinement_t. */
enum { AT_LARGEBINARY = 19, AT_LARGEUTF8 = 20, AT_LARGELIST = 21 };

/* Sanity cap on vector element counts from crafted input. */
enum { MAX_VECTOR_ELEMS = 1 << 20 };

/* ---- FlatBuffer reader over a bounded byte range ---- */

typedef struct { const uint8_t* buf; size_t len; } fbr;

static int fbr_u8(const fbr* r, size_t pos, uint8_t* out) {
    if (pos + 1 > r->len) return 0;
    *out = r->buf[pos];
    return 1;
}
static int fbr_u16(const fbr* r, size_t pos, uint16_t* out) {
    if (pos + 2 > r->len) return 0;
    *out = (uint16_t)(r->buf[pos] | ((uint16_t)r->buf[pos + 1] << 8));
    return 1;
}
static int fbr_u32(const fbr* r, size_t pos, uint32_t* out) {
    if (pos + 4 > r->len) return 0;
    *out = (uint32_t)r->buf[pos] | ((uint32_t)r->buf[pos + 1] << 8) |
           ((uint32_t)r->buf[pos + 2] << 16) | ((uint32_t)r->buf[pos + 3] << 24);
    return 1;
}
static int fbr_i32(const fbr* r, size_t pos, int32_t* out) {
    uint32_t u;
    if (!fbr_u32(r, pos, &u)) return 0;
    *out = (int32_t)u;
    return 1;
}

/* Follow the forward uoffset stored at `pos`; sets *out to the target pos. */
static int fbr_indirect(const fbr* r, size_t pos, size_t* out) {
    uint32_t off;
    if (!fbr_u32(r, pos, &off) || off == 0) return 0;
    size_t target = pos + off;
    if (target < pos || target > r->len) return 0;  /* overflow / OOB */
    *out = target;
    return 1;
}

/* Locate field `field_id` of the table at `table`. On success sets *out to the
 * field's data position, or 0 when the field is absent. Returns 0 on a
 * structurally invalid table/vtable. */
static int fbr_field(const fbr* r, size_t table, int field_id, size_t* out) {
    int32_t soffset;
    if (!fbr_i32(r, table, &soffset)) return 0;
    int64_t vt_signed = (int64_t)table - (int64_t)soffset;
    if (vt_signed < 0 || (uint64_t)vt_signed >= r->len) return 0;
    size_t vt = (size_t)vt_signed;

    uint16_t vt_size;
    if (!fbr_u16(r, vt, &vt_size)) return 0;
    size_t slot = 4 + (size_t)field_id * 2;
    if (slot + 2 > vt_size) { *out = 0; return 1; }  /* field not in vtable */

    uint16_t voff;
    if (!fbr_u16(r, vt + slot, &voff)) return 0;
    if (voff == 0) { *out = 0; return 1; }           /* field absent */
    size_t fpos = table + voff;
    if (fpos < table || fpos > r->len) return 0;
    *out = fpos;
    return 1;
}

/* Read a FlatBuffer string located at `str_pos` into a NUL-terminated arena
 * copy. Returns NULL on OOB or OOM. */
static char* fbr_string_at(const fbr* r, size_t str_pos, carquet_arena_t* arena) {
    uint32_t slen;
    if (!fbr_u32(r, str_pos, &slen)) return NULL;
    if (str_pos + 4 + slen < str_pos || str_pos + 4 + slen > r->len) return NULL;
    char* s = (char*)carquet_arena_alloc(arena, (size_t)slen + 1);
    if (!s) return NULL;
    memcpy(s, r->buf + str_pos + 4, slen);
    s[slen] = '\0';
    return s;
}

/* Read a string-typed field of a table. Returns NULL when absent/invalid. */
static char* fbr_field_string(const fbr* r, size_t table, int field_id,
                              carquet_arena_t* arena) {
    size_t f;
    if (!fbr_field(r, table, field_id, &f) || f == 0) return NULL;
    size_t sp;
    if (!fbr_indirect(r, f, &sp)) return NULL;
    return fbr_string_at(r, sp, arena);
}

/* Resolve a vector field of a table: sets *vec_data to the position of the
 * first element (past the count prefix) and *count to a validated element
 * count. Returns 0 when the field is absent or the span is inconsistent. */
static int fbr_field_vector(const fbr* r, size_t table, int field_id,
                            size_t elem_size, size_t* vec_data, uint32_t* count) {
    size_t f;
    if (!fbr_field(r, table, field_id, &f) || f == 0) return 0;
    size_t vec;
    if (!fbr_indirect(r, f, &vec)) return 0;
    uint32_t n;
    if (!fbr_u32(r, vec, &n)) return 0;
    if (n > MAX_VECTOR_ELEMS) return 0;
    size_t data = vec + 4;
    size_t span = (size_t)n * elem_size;
    if (data + span < data || data + span > r->len) return 0;
    *vec_data = data;
    *count = n;
    return 1;
}

/* ---- base64 decode (tolerant: skips whitespace/newlines, stops at pad) ---- */

static int b64_val(unsigned char c) {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
}

static uint8_t* base64_decode(const char* in, size_t* out_len) {
    size_t in_len = strlen(in);
    uint8_t* out = (uint8_t*)carquet_mem_malloc(in_len / 4 * 3 + 4);
    if (!out) return NULL;
    size_t o = 0;
    int quad[4];
    int qn = 0;
    for (size_t i = 0; i < in_len; i++) {
        int v = b64_val((unsigned char)in[i]);
        if (v < 0) continue;  /* skip '=', newlines, stray bytes */
        quad[qn++] = v;
        if (qn == 4) {
            out[o++] = (uint8_t)((quad[0] << 2) | (quad[1] >> 4));
            out[o++] = (uint8_t)((quad[1] << 4) | (quad[2] >> 2));
            out[o++] = (uint8_t)((quad[2] << 6) | quad[3]);
            qn = 0;
        }
    }
    if (qn >= 2) {
        out[o++] = (uint8_t)((quad[0] << 2) | (quad[1] >> 4));
        if (qn >= 3) out[o++] = (uint8_t)((quad[1] << 4) | (quad[2] >> 2));
    }
    *out_len = o;
    return out;
}

/* Map a Field's inline type_type (union tag, u8) to a carquet refinement.
 * Returns 0 (CARQUET_ARROW_REFINE_NONE) when absent or not a 64-bit variant. */
static int32_t field_type_refinement(const fbr* r, size_t field_tbl) {
    size_t f;
    if (!fbr_field(r, field_tbl, FIELD_TYPETYPE_SLOT, &f) || f == 0) return 0;
    uint8_t tag;
    if (!fbr_u8(r, f, &tag)) return 0;
    switch (tag) {
        case AT_LARGEUTF8:   return 1;  /* CARQUET_ARROW_REFINE_LARGE_UTF8   */
        case AT_LARGEBINARY: return 2;  /* CARQUET_ARROW_REFINE_LARGE_BINARY */
        case AT_LARGELIST:   return 3;  /* CARQUET_ARROW_REFINE_LARGE_LIST   */
        default:             return 0;
    }
}

/* ---- attach metadata to the matching top-level schema element ---- */

static int32_t match_element(const char* name,
                             parquet_schema_element_t* elements,
                             int32_t num_elements,
                             const int32_t* parent_indices,
                             const uint8_t* used) {
    if (!name) return -1;
    /* Match by name among direct children of the root that have no metadata
     * yet (used[] guards against duplicate names being over-written). */
    for (int32_t i = 1; i < num_elements; i++) {
        int32_t parent = parent_indices ? parent_indices[i] : 0;
        if (parent != 0) continue;
        if (used[i]) continue;
        if (elements[i].name && strcmp(elements[i].name, name) == 0) return i;
    }
    return -1;
}

void carquet_apply_arrow_field_metadata(
    const char* b64_value,
    parquet_schema_element_t* elements,
    int32_t num_elements,
    const int32_t* parent_indices,
    carquet_arena_t* arena) {

    if (!b64_value || !elements || num_elements < 2 || !arena) return;

    size_t raw_len = 0;
    uint8_t* raw = base64_decode(b64_value, &raw_len);
    if (!raw) return;

    /* Skip the Arrow IPC encapsulation prefix: either the modern
     * 0xFFFFFFFF continuation marker + u32 length, or the legacy bare u32. */
    size_t fb_start;
    if (raw_len >= 8 && raw[0] == 0xFF && raw[1] == 0xFF &&
        raw[2] == 0xFF && raw[3] == 0xFF) {
        fb_start = 8;
    } else if (raw_len >= 4) {
        fb_start = 4;
    } else {
        carquet_mem_free(raw);
        return;
    }

    fbr r = { raw + fb_start, raw_len - fb_start };

    uint8_t* used = (uint8_t*)carquet_mem_calloc((size_t)num_elements, 1);
    if (!used) { carquet_mem_free(raw); return; }

    size_t msg, hdr_field, schema_tbl, fields_data;
    uint32_t nfields;
    if (!fbr_indirect(&r, 0, &msg)) goto done;                         /* root Message */
    if (!fbr_field(&r, msg, MSG_HEADER_SLOT, &hdr_field) || hdr_field == 0) goto done;
    if (!fbr_indirect(&r, hdr_field, &schema_tbl)) goto done;          /* Schema */
    if (!fbr_field_vector(&r, schema_tbl, SCHEMA_FIELDS_SLOT, 4,
                          &fields_data, &nfields)) goto done;

    for (uint32_t j = 0; j < nfields; j++) {
        size_t field_tbl;
        if (!fbr_indirect(&r, fields_data + (size_t)j * 4, &field_tbl)) continue;

        char* name = fbr_field_string(&r, field_tbl, FIELD_NAME_SLOT, arena);

        /* Match the Arrow field to its Parquet element up front so both the
         * type refinement and the custom_metadata can be attached, even when
         * the field carries no metadata. */
        int32_t idx = match_element(name, elements, num_elements,
                                    parent_indices, used);
        if (idx < 0) continue;
        used[idx] = 1;

        /* Type refinement: a 64-bit-offset Arrow type Parquet can't express. */
        int32_t refine = field_type_refinement(&r, field_tbl);
        if (refine != 0) elements[idx].arrow_type_refinement = refine;

        /* custom_metadata (variable labels/descriptions). */
        size_t meta_data;
        uint32_t nmeta;
        if (!fbr_field_vector(&r, field_tbl, FIELD_META_SLOT, 4,
                              &meta_data, &nmeta) || nmeta == 0) continue;

        parquet_key_value_t* kvs = (parquet_key_value_t*)carquet_arena_calloc(
            arena, nmeta, sizeof(parquet_key_value_t));
        if (!kvs) continue;

        int32_t got = 0;
        for (uint32_t k = 0; k < nmeta; k++) {
            size_t kv_tbl;
            if (!fbr_indirect(&r, meta_data + (size_t)k * 4, &kv_tbl)) continue;
            char* key = fbr_field_string(&r, kv_tbl, KV_KEY_SLOT, arena);
            if (!key) continue;  /* KeyValue.key is required to be useful */
            kvs[got].key = key;
            kvs[got].value = fbr_field_string(&r, kv_tbl, KV_VALUE_SLOT, arena);
            got++;
        }
        if (got == 0) continue;

        elements[idx].field_metadata = kvs;
        elements[idx].num_field_metadata = got;
    }

done:
    carquet_mem_free(used);
    carquet_mem_free(raw);
}
