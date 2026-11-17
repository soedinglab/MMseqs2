/**
 * @file fb_read.h
 * @brief Bounds-checked FlatBuffer reading and base64 decoding for the
 *        "ARROW:schema" footer blob (untrusted input). Shared by
 *        arrow_schema_read.c and arrow_restore.c.
 */
#ifndef CARQUET_READER_FB_READ_H
#define CARQUET_READER_FB_READ_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "core/allocator.h"
#include "core/arena.h"

/* Sanity cap on vector element counts from crafted input. */
enum { MAX_VECTOR_ELEMS = 1 << 20 };

/* ---- FlatBuffer reader over a bounded byte range ---- */

typedef struct { const uint8_t* buf; size_t len; } fbr;

static inline int fbr_u8(const fbr* r, size_t pos, uint8_t* out) {
    if (pos + 1 > r->len) return 0;
    *out = r->buf[pos];
    return 1;
}
static inline int fbr_u16(const fbr* r, size_t pos, uint16_t* out) {
    if (pos + 2 > r->len) return 0;
    *out = (uint16_t)(r->buf[pos] | ((uint16_t)r->buf[pos + 1] << 8));
    return 1;
}
static inline int fbr_u32(const fbr* r, size_t pos, uint32_t* out) {
    if (pos + 4 > r->len) return 0;
    *out = (uint32_t)r->buf[pos] | ((uint32_t)r->buf[pos + 1] << 8) |
           ((uint32_t)r->buf[pos + 2] << 16) | ((uint32_t)r->buf[pos + 3] << 24);
    return 1;
}
static inline int fbr_i32(const fbr* r, size_t pos, int32_t* out) {
    uint32_t u;
    if (!fbr_u32(r, pos, &u)) return 0;
    *out = (int32_t)u;
    return 1;
}

/* Follow the forward uoffset stored at `pos`; sets *out to the target pos. */
static inline int fbr_indirect(const fbr* r, size_t pos, size_t* out) {
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
static inline int fbr_field(const fbr* r, size_t table, int field_id, size_t* out) {
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
 * copy. Returns NULL on OOB, OOM, or when the copy would exceed *budget.
 * The budget bounds the total bytes copied out of the blob: FlatBuffer
 * offsets let many table entries share one large string, so without it a
 * small ARROW:schema value could make the reader copy N x S bytes. */
static inline char* fbr_string_at(const fbr* r, size_t str_pos, carquet_arena_t* arena,
                           size_t* budget) {
    uint32_t slen;
    if (!fbr_u32(r, str_pos, &slen)) return NULL;
    if (str_pos + 4 + slen < str_pos || str_pos + 4 + slen > r->len) return NULL;
    if ((size_t)slen + 1 > *budget) return NULL;
    *budget -= (size_t)slen + 1;
    char* s = (char*)carquet_arena_alloc(arena, (size_t)slen + 1);
    if (!s) return NULL;
    memcpy(s, r->buf + str_pos + 4, slen);
    s[slen] = '\0';
    return s;
}

/* Read a string-typed field of a table. Returns NULL when absent/invalid. */
static inline char* fbr_field_string(const fbr* r, size_t table, int field_id,
                              carquet_arena_t* arena, size_t* budget) {
    size_t f;
    if (!fbr_field(r, table, field_id, &f) || f == 0) return NULL;
    size_t sp;
    if (!fbr_indirect(r, f, &sp)) return NULL;
    return fbr_string_at(r, sp, arena, budget);
}

/* View a string-typed field of a table in place (no copy). Returns 0 when
 * absent/invalid. */
static inline int fbr_field_string_view(const fbr* r, size_t table, int field_id,
                                 const char** str, size_t* len) {
    size_t f;
    if (!fbr_field(r, table, field_id, &f) || f == 0) return 0;
    size_t sp;
    if (!fbr_indirect(r, f, &sp)) return 0;
    uint32_t slen;
    if (!fbr_u32(r, sp, &slen)) return 0;
    if (sp + 4 + slen < sp || sp + 4 + slen > r->len) return 0;
    *str = (const char*)(r->buf + sp + 4);
    *len = slen;
    return 1;
}

/* Resolve a vector field of a table: sets *vec_data to the position of the
 * first element (past the count prefix) and *count to a validated element
 * count. Returns 0 when the field is absent or the span is inconsistent. */
static inline int fbr_field_vector(const fbr* r, size_t table, int field_id,
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

static inline int b64_val(unsigned char c) {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
}

static inline uint8_t* base64_decode(const char* in, size_t* out_len) {
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


#endif /* CARQUET_READER_FB_READ_H */
