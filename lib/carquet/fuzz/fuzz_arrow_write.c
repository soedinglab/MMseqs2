/**
 * @file fuzz_arrow_write.c
 * @brief Differential roundtrip fuzzer for the Arrow C Data import/write path.
 *
 * Builds a *valid* Arrow struct array from fuzz bytes — flat primitives
 * (incl. the coerced int8/uint16/uint64, date32, timestamp, decimal128,
 * fixed-size binary, large_utf8 and dictionary-encoded columns) plus LIST,
 * LARGE_LIST, STRUCT and MAP children, with fuzz-chosen validity bitmaps and
 * slice offsets at every level — then:
 *
 *   carquet_arrow_import_schema -> carquet_writer_create_buffer (fuzz-chosen
 *   codec / page version / page size / per-leaf encodings / bloom / page
 *   index) -> carquet_writer_write_arrow -> close -> carquet_reader_open_buffer
 *   -> carquet_reader_read_arrow (every row group)
 *
 * and compares every value of the read-back arrays against the source arrays
 * with a format-driven comparator (so coerced types compare by value, not by
 * Arrow format). Every call must succeed on these well-formed inputs, every
 * value must survive, and every non-null value of a flat INT32 / INT64 /
 * STRING column with a bloom filter must test positive in that filter.
 *
 * fuzz_arrow_c_data only covers the schema format-string parser and
 * fuzz_nested_write the single-level LIST/MAP helper; neither checks values.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <carquet/carquet.h>

#define MAX_COLS 6
#define MAX_ROWS 300
#define NUM_KINDS 20

typedef struct { const uint8_t* p; size_t n; } cur_t;

static uint8_t cu8(cur_t* c) {
    if (!c->n) return 0;
    c->n--;
    return *c->p++;
}
static void cbytes(cur_t* c, uint8_t* out, size_t n) {
    for (size_t i = 0; i < n; i++) out[i] = cu8(c);
}

/* ── Arrow struct construction (everything heap-owned, freed by release) ── */

static void rel_schema(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    for (int64_t i = 0; i < s->n_children; i++) {
        struct ArrowSchema* k = s->children[i];
        if (k) { if (k->release) k->release(k); free(k); }
    }
    free(s->children);
    if (s->dictionary) {
        if (s->dictionary->release) s->dictionary->release(s->dictionary);
        free(s->dictionary);
    }
    s->release = NULL;
}

static void rel_array(struct ArrowArray* a) {
    if (!a || !a->release) return;
    for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]);
    free(a->buffers);
    for (int64_t i = 0; i < a->n_children; i++) {
        struct ArrowArray* k = a->children[i];
        if (k) { if (k->release) k->release(k); free(k); }
    }
    free(a->children);
    if (a->dictionary) {
        if (a->dictionary->release) a->dictionary->release(a->dictionary);
        free(a->dictionary);
    }
    a->release = NULL;
}

static struct ArrowSchema* sch(const char* fmt, const char* name, bool nullable,
                               int64_t nch) {
    struct ArrowSchema* s = calloc(1, sizeof(*s));
    s->format = fmt;
    s->name = name;
    s->flags = nullable ? ARROW_FLAG_NULLABLE : 0;
    s->n_children = nch;
    s->children = nch ? calloc((size_t)nch, sizeof(*s->children)) : NULL;
    s->release = rel_schema;
    return s;
}

static struct ArrowArray* arr(int64_t phys_len, int64_t off, int64_t nbuf,
                              int64_t nch) {
    struct ArrowArray* a = calloc(1, sizeof(*a));
    a->offset = off;
    a->length = phys_len - off;
    a->n_buffers = nbuf;
    a->buffers = calloc((size_t)nbuf, sizeof(*a->buffers));
    a->n_children = nch;
    a->children = nch ? calloc((size_t)nch, sizeof(*a->children)) : NULL;
    a->release = rel_array;
    return a;
}

static bool bit(const uint8_t* bm, int64_t i) { return (bm[i >> 3] >> (i & 7)) & 1; }

/* Validity bitmap over [0, phys_len); sets null_count over the logical window. */
static void set_validity(cur_t* c, struct ArrowArray* a, int64_t phys_len,
                         bool nullable) {
    a->null_count = 0;
    if (!nullable || (cu8(c) & 3) == 0) return;          /* no bitmap */
    uint8_t* bm = calloc((size_t)(phys_len / 8 + 1), 1);
    for (int64_t i = 0; i < phys_len; i++)
        if ((cu8(c) & 3) != 0) bm[i >> 3] |= (uint8_t)(1u << (i & 7));
    for (int64_t i = a->offset; i < phys_len; i++)
        if (!bit(bm, i)) a->null_count++;
    a->buffers[0] = bm;
}

static int32_t* make_offsets32(cur_t* c, int64_t n, int max_len) {
    int32_t* o = malloc(sizeof(int32_t) * (size_t)(n + 1));
    o[0] = 0;
    for (int64_t i = 0; i < n; i++) o[i + 1] = o[i] + (int32_t)(cu8(c) % (max_len + 1));
    return o;
}

static const char* kind_format(int k) {
    static const char* f[NUM_KINDS] = {
        "i", "l", "f", "g", "b", "u", "U", "z", "c", "S",
        "tdD", "tsu:UTC", "w:5", "+l", "+s", "i" /* dict idx */, "+L", "+m",
        "d:9,2", "L"
    };
    return f[k];
}

static int fixed_width(int k) {
    switch (k) {
        case 0: case 10: case 15: return 4;
        case 1: case 11: case 19: return 8;
        case 2: return 4;
        case 3: return 8;
        case 8: return 1;
        case 9: return 2;
        case 12: return 5;
        case 18: return 16;
        default: return 0;
    }
}

static int kind_leaves(int k) { return (k == 14 || k == 17) ? 2 : 1; }

static const char* g_names[MAX_COLS] = {"c0", "c1", "c2", "c3", "c4", "c5"};

static void build_bytes_leaf(cur_t* c, struct ArrowArray* a, int64_t phys_len,
                             bool large) {
    int32_t* o32 = NULL;
    int64_t* o64 = NULL;
    int64_t total = 0;
    if (large) o64 = malloc(sizeof(int64_t) * (size_t)(phys_len + 1));
    else o32 = malloc(sizeof(int32_t) * (size_t)(phys_len + 1));
    if (large) o64[0] = 0; else o32[0] = 0;
    for (int64_t i = 0; i < phys_len; i++) {
        total += cu8(c) % 9;
        if (large) o64[i + 1] = total; else o32[i + 1] = (int32_t)total;
    }
    uint8_t* d = malloc((size_t)total + 1);
    /* Small alphabet so dictionary / DELTA_BYTE_ARRAY see repeats and prefixes. */
    for (int64_t i = 0; i < total; i++) d[i] = (uint8_t)('a' + (cu8(c) & 3));
    a->buffers[1] = large ? (void*)o64 : (void*)o32;
    a->buffers[2] = d;
}

/* Build one column (schema + array) of kind k with physical length phys_len
 * and own slice offset off. */
static void build_col(cur_t* c, int k, const char* name, bool nullable,
                      int64_t phys_len, int64_t off,
                      struct ArrowSchema** out_s, struct ArrowArray** out_a) {
    struct ArrowSchema* s = NULL;
    struct ArrowArray* a = NULL;
    int w = fixed_width(k);
    switch (k) {
    case 4: {   /* bool: bit-packed data */
        s = sch("b", name, nullable, 0);
        a = arr(phys_len, off, 2, 0);
        uint8_t* bm = calloc((size_t)(phys_len / 8 + 1), 1);
        for (int64_t i = 0; i < phys_len; i++)
            if (cu8(c) & 1) bm[i >> 3] |= (uint8_t)(1u << (i & 7));
        a->buffers[1] = bm;
        break;
    }
    case 5: case 6: case 7:
        s = sch(kind_format(k), name, nullable, 0);
        a = arr(phys_len, off, 3, 0);
        build_bytes_leaf(c, a, phys_len, k == 6);
        break;
    case 13: case 16: {   /* list<int32> / large_list<double> */
        bool large = (k == 16);
        s = sch(large ? "+L" : "+l", name, nullable, 1);
        a = arr(phys_len, off, 2, 1);
        int32_t* o32 = make_offsets32(c, phys_len, 4);
        int64_t child_len = o32[phys_len];
        if (large) {
            int64_t* o64 = malloc(sizeof(int64_t) * (size_t)(phys_len + 1));
            for (int64_t i = 0; i <= phys_len; i++) o64[i] = o32[i];
            free(o32);
            a->buffers[1] = o64;
        } else {
            a->buffers[1] = o32;
        }
        build_col(c, large ? 3 : 0, "element", (cu8(c) & 1) != 0, child_len, 0,
                  &s->children[0], &a->children[0]);
        break;
    }
    case 14: {  /* struct<a:int32, b:utf8> */
        s = sch("+s", name, nullable, 2);
        a = arr(phys_len, off, 1, 2);
        build_col(c, 0, "a", true, phys_len, 0, &s->children[0], &a->children[0]);
        build_col(c, 5, "b", true, phys_len, 0, &s->children[1], &a->children[1]);
        break;
    }
    case 15: {  /* dictionary<int32 -> utf8> */
        s = sch("i", name, nullable, 0);
        a = arr(phys_len, off, 2, 0);
        int64_t dict_len = 1 + cu8(c) % 5;
        int32_t* idx = malloc(sizeof(int32_t) * (size_t)(phys_len + 1));
        for (int64_t i = 0; i < phys_len; i++) idx[i] = (int32_t)(cu8(c) % dict_len);
        a->buffers[1] = idx;
        build_col(c, 5, NULL, false, dict_len, 0, &s->dictionary, &a->dictionary);
        break;
    }
    case 17: {  /* map<utf8 (non-null), int64> */
        s = sch("+m", name, nullable, 1);
        a = arr(phys_len, off, 2, 1);
        int32_t* o32 = make_offsets32(c, phys_len, 3);
        a->buffers[1] = o32;
        int64_t n = o32[phys_len];
        struct ArrowSchema* es = sch("+s", "entries", false, 2);
        struct ArrowArray* ea = arr(n, 0, 1, 2);
        build_col(c, 5, "key", false, n, 0, &es->children[0], &ea->children[0]);
        build_col(c, 1, "value", true, n, 0, &es->children[1], &ea->children[1]);
        s->children[0] = es;
        a->children[0] = ea;
        break;
    }
    default: {  /* fixed-width primitives */
        s = sch(kind_format(k), name, nullable, 0);
        a = arr(phys_len, off, 2, 0);
        uint8_t* d = calloc((size_t)phys_len + 1, (size_t)w);
        for (int64_t i = 0; i < phys_len; i++) {
            uint8_t* v = d + (size_t)i * (size_t)w;
            if (k == 18) {
                /* decimal128(9,2): keep |unscaled| < 10^9 so it fits. */
                int32_t x;
                cbytes(c, (uint8_t*)&x, 4);
                int64_t u = (int64_t)(x % 1000000000);
                memcpy(v, &u, 8);
                memset(v + 8, u < 0 ? 0xFF : 0x00, 8);
            } else if (k == 2 || k == 3) {
                /* Mostly small integral values (dictionary-friendly), but
                 * also raw bit patterns: NaN payloads, -0.0, infinities. */
                if (cu8(c) & 1) cbytes(c, v, (size_t)w);
                else if (w == 4) { float f = (float)(int8_t)cu8(c); memcpy(v, &f, 4); }
                else { double f = (double)(int8_t)cu8(c); memcpy(v, &f, 8); }
            } else if (cu8(c) & 1) {
                cbytes(c, v, (size_t)w);
            } else {
                /* Small signed value, sign-extended (delta / dictionary friendly). */
                int64_t sv = (int8_t)cu8(c);
                memcpy(v, &sv, (size_t)(w < 8 ? w : 8));
                if (w > 8) memset(v + 8, sv < 0 ? 0xFF : 0, (size_t)w - 8);
            }
        }
        a->buffers[1] = d;
        break;
    }
    }
    set_validity(c, a, phys_len, nullable);
    *out_s = s;
    *out_a = a;
}

typedef struct {
    int ncols;
    int kinds[MAX_COLS];
    bool nullable[MAX_COLS];
    int64_t child_off[MAX_COLS];
    int64_t nrows;
    int64_t top_off;
} spec_t;

static void build_top(const spec_t* sp, cur_t c, struct ArrowSchema** os,
                      struct ArrowArray** oa) {
    struct ArrowSchema* s = sch("+s", "schema", false, sp->ncols);
    struct ArrowArray* a = arr(sp->top_off + sp->nrows, sp->top_off, 1, sp->ncols);
    for (int i = 0; i < sp->ncols; i++) {
        int64_t phys = sp->child_off[i] + sp->top_off + sp->nrows;
        build_col(&c, sp->kinds[i], g_names[i], sp->nullable[i], phys,
                  sp->child_off[i], &s->children[i], &a->children[i]);
    }
    *os = s;
    *oa = a;
}

/* ── Format-driven comparator ─────────────────────────────────────────── */

typedef enum { C_INT, C_FLOAT, C_BOOL, C_BYTES, C_DEC, C_LIST, C_LLIST,
               C_STRUCT, C_MAP, C_OTHER } cat_t;

static cat_t classify(const char* f, int* width, bool* is_unsigned) {
    *width = 0;
    *is_unsigned = false;
    switch (f[0]) {
        case 'c': *width = 1; return C_INT;
        case 'C': *width = 1; *is_unsigned = true; return C_INT;
        case 's': *width = 2; return C_INT;
        case 'S': *width = 2; *is_unsigned = true; return C_INT;
        case 'i': *width = 4; return C_INT;
        case 'I': *width = 4; *is_unsigned = true; return C_INT;
        case 'l': *width = 8; return C_INT;
        case 'L': *width = 8; *is_unsigned = true; return C_INT;
        case 'f': *width = 4; return C_FLOAT;
        case 'g': *width = 8; return C_FLOAT;
        case 'b': return C_BOOL;
        case 'u': case 'z': *width = 4; return C_BYTES;
        case 'U': case 'Z': *width = 8; return C_BYTES;
        case 'w': *width = atoi(f + 2); return C_BYTES;
        case 'd': {
            const char* p = strchr(f, ',');
            p = p ? strchr(p + 1, ',') : NULL;
            *width = p ? atoi(p + 1) / 8 : 16;
            return C_DEC;
        }
        case 't':
            if (f[1] == 'd') { *width = f[2] == 'D' ? 4 : 8; return C_INT; }
            if (f[1] == 't') { *width = (f[2] == 's' || f[2] == 'm') ? 4 : 8; return C_INT; }
            if (f[1] == 's') { *width = 8; return C_INT; }
            return C_OTHER;
        case '+':
            if (f[1] == 'l') return C_LIST;
            if (f[1] == 'L') return C_LLIST;
            if (f[1] == 's') return C_STRUCT;
            if (f[1] == 'm') return C_MAP;
            return C_OTHER;
        default: return C_OTHER;
    }
}

static bool valid_at(const struct ArrowArray* a, int64_t i) {
    const uint8_t* bm = (const uint8_t*)a->buffers[0];
    if (!bm || a->null_count == 0) return true;
    return bit(bm, a->offset + i);
}

static int64_t int_at(const struct ArrowArray* a, int width, bool uns, int64_t i) {
    const uint8_t* d = (const uint8_t*)a->buffers[1] + (size_t)(a->offset + i) * (size_t)width;
    switch (width) {
        case 1: return uns ? (int64_t)*d : (int64_t)(int8_t)*d;
        case 2: { uint16_t v; memcpy(&v, d, 2); return uns ? (int64_t)v : (int64_t)(int16_t)v; }
        case 4: { uint32_t v; memcpy(&v, d, 4); return uns ? (int64_t)v : (int64_t)(int32_t)v; }
        default: { int64_t v; memcpy(&v, d, 8); return v; }
    }
}

static void bytes_at(const struct ArrowArray* a, int width, const char* f,
                     int64_t i, const uint8_t** p, int64_t* len) {
    int64_t phys = a->offset + i;
    if (f[0] == 'w') {
        *p = (const uint8_t*)a->buffers[1] + phys * width;
        *len = width;
    } else if (width == 4) {
        const int32_t* o = (const int32_t*)a->buffers[1];
        *p = (const uint8_t*)a->buffers[2] + o[phys];
        *len = o[phys + 1] - o[phys];
    } else {
        const int64_t* o = (const int64_t*)a->buffers[1];
        *p = (const uint8_t*)a->buffers[2] + o[phys];
        *len = o[phys + 1] - o[phys];
    }
}

static void list_range(const struct ArrowArray* a, bool large, int64_t i,
                       int64_t* start, int64_t* end) {
    int64_t phys = a->offset + i;
    if (large) {
        const int64_t* o = (const int64_t*)a->buffers[1];
        *start = o[phys]; *end = o[phys + 1];
    } else {
        const int32_t* o = (const int32_t*)a->buffers[1];
        *start = o[phys]; *end = o[phys + 1];
    }
}

static const char* g_fail_path;

#define MISMATCH(what) do { \
    fprintf(stderr, "fuzz_arrow_write: mismatch at %s: %s (src fmt %s, read fmt %s, src idx %lld, read idx %lld)\n", \
            g_fail_path, what, ss->format, rs->format, (long long)si, (long long)ri); \
    __builtin_trap(); \
} while (0)

static void cmp_value(const struct ArrowSchema* ss, const struct ArrowArray* sa, int64_t si,
                      const struct ArrowSchema* rs, const struct ArrowArray* ra, int64_t ri) {
    bool sv = valid_at(sa, si);
    bool rv = valid_at(ra, ri);
    if (sv != rv) MISMATCH("validity");
    if (!sv) return;

    /* Resolve a source dictionary column to its value. */
    if (ss->dictionary) {
        int64_t idx = int_at(sa, 4, false, si);
        cmp_value(ss->dictionary, sa->dictionary, idx, rs, ra, ri);
        return;
    }

    int sw, rw;
    bool su, ru;
    cat_t sc = classify(ss->format, &sw, &su);
    cat_t rc = classify(rs->format, &rw, &ru);

    if (sc == C_DEC || rc == C_DEC) {
        /* Compare unscaled values (fit in int64 by construction). */
        int64_t a = 0, b = 0;
        if (sc == C_DEC) {
            memcpy(&a, (const uint8_t*)sa->buffers[1] + (size_t)(sa->offset + si) * (size_t)sw, 8);
        } else if (sc == C_INT) {
            a = int_at(sa, sw, su, si);
        } else MISMATCH("decimal category");
        if (rc == C_DEC) {
            const uint8_t* d = (const uint8_t*)ra->buffers[1] + (size_t)(ra->offset + ri) * (size_t)rw;
            if (rw >= 8) memcpy(&b, d, 8);
            else { int32_t t; memcpy(&t, d, 4); b = t; }
        } else if (rc == C_INT) {
            b = int_at(ra, rw, ru, ri);
        } else MISMATCH("decimal category");
        if (a != b) MISMATCH("decimal value");
        return;
    }
    if (sc != rc && !(sc == C_LIST && rc == C_LLIST) && !(sc == C_LLIST && rc == C_LIST))
        MISMATCH("category");

    switch (sc) {
    case C_INT: {
        int64_t a = int_at(sa, sw, su, si), b = int_at(ra, rw, ru, ri);
        if (a != b) MISMATCH("int value");
        break;
    }
    case C_FLOAT: {
        double a, b;
        if (sw == 4) { float f; memcpy(&f, (const uint8_t*)sa->buffers[1] + (sa->offset + si) * 4, 4); a = f; }
        else memcpy(&a, (const uint8_t*)sa->buffers[1] + (sa->offset + si) * 8, 8);
        if (rw == 4) { float f; memcpy(&f, (const uint8_t*)ra->buffers[1] + (ra->offset + ri) * 4, 4); b = f; }
        else memcpy(&b, (const uint8_t*)ra->buffers[1] + (ra->offset + ri) * 8, 8);
        if (isnan(a) && isnan(b)) break;
        if (sw == rw && sw == 8 ? memcmp(&a, &b, 8) != 0
                                : (a != b || signbit(a) != signbit(b)))
            MISMATCH("float value");
        break;
    }
    case C_BOOL: {
        bool a = bit((const uint8_t*)sa->buffers[1], sa->offset + si);
        bool b = bit((const uint8_t*)ra->buffers[1], ra->offset + ri);
        if (a != b) MISMATCH("bool value");
        break;
    }
    case C_BYTES: {
        const uint8_t *pa, *pb;
        int64_t la, lb;
        bytes_at(sa, sw, ss->format, si, &pa, &la);
        bytes_at(ra, rw, rs->format, ri, &pb, &lb);
        if (la != lb || (la && memcmp(pa, pb, (size_t)la) != 0)) MISMATCH("bytes value");
        break;
    }
    case C_LIST: case C_LLIST: case C_MAP: {
        int64_t s0, s1, r0, r1;
        list_range(sa, sc == C_LLIST, si, &s0, &s1);
        list_range(ra, rc == C_LLIST, ri, &r0, &r1);
        if (s1 - s0 != r1 - r0) MISMATCH("list length");
        for (int64_t j = 0; j < s1 - s0; j++)
            cmp_value(ss->children[0], sa->children[0], s0 + j,
                      rs->children[0], ra->children[0], r0 + j);
        break;
    }
    case C_STRUCT: {
        if (ss->n_children != rs->n_children) MISMATCH("struct arity");
        for (int64_t j = 0; j < ss->n_children; j++)
            cmp_value(ss->children[j], sa->children[j], sa->offset + si,
                      rs->children[j], ra->children[j], ra->offset + ri);
        break;
    }
    default:
        MISMATCH("unhandled category");
    }
}

/* ── Bloom no-false-negative check for flat INT32 / INT64 / STRING ─────── */

static void check_bloom(carquet_reader_t* r, int32_t rg, int32_t leaf, int kind,
                        const struct ArrowSchema* ss, const struct ArrowArray* sa,
                        int64_t first, int64_t n) {
    carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(r, rg, leaf, NULL);
    if (!bf) return;
    for (int64_t i = first; i < first + n; i++) {
        if (!valid_at(sa, i)) continue;
        bool hit = true;
        if (kind == 0) hit = carquet_bloom_filter_check_i32(bf, (int32_t)int_at(sa, 4, false, i));
        else if (kind == 1) hit = carquet_bloom_filter_check_i64(bf, int_at(sa, 8, false, i));
        else if (kind == 5) {
            const uint8_t* p; int64_t len;
            bytes_at(sa, 4, ss->format, i, &p, &len);
            hit = carquet_bloom_filter_check_bytes(bf, p, (size_t)len);
        }
        if (!hit) {
            fprintf(stderr, "fuzz_arrow_write: bloom false negative (leaf %d kind %d row %lld)\n",
                    leaf, kind, (long long)i);
            __builtin_trap();
        }
    }
    carquet_bloom_filter_destroy(bf);
}

#define REQUIRE(cond, what) do { if (!(cond)) { \
    fprintf(stderr, "fuzz_arrow_write: %s failed: %s\n", what, err.message); \
    __builtin_trap(); } } while (0)

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 8 || size > 65536) return 0;
    (void)carquet_init();
    cur_t c = {data, size};

    spec_t sp;
    memset(&sp, 0, sizeof(sp));
    uint8_t h0 = cu8(&c), h1 = cu8(&c), h2 = cu8(&c), h3 = cu8(&c);
    sp.ncols = 1 + h0 % MAX_COLS;
    sp.nrows = (int64_t)(cu8(&c) | ((uint32_t)(h1 & 1) << 8)) % (MAX_ROWS + 1);
    sp.top_off = (h1 >> 1) & 3;
    for (int i = 0; i < sp.ncols; i++) {
        uint8_t b = cu8(&c);
        sp.kinds[i] = b % NUM_KINDS;
        sp.nullable[i] = (b & 0x80) != 0;
        sp.child_off[i] = (b >> 5) & 3;
    }
    uint8_t enc_sel[MAX_COLS * 2];
    for (int i = 0; i < MAX_COLS * 2; i++) enc_sel[i] = cu8(&c);

    /* The rest of the input is value bytes; every build re-reads them from
     * the same position, so the three builds below are identical. */
    cur_t vals = c;

    static const carquet_compression_t codecs[] = {
        CARQUET_COMPRESSION_UNCOMPRESSED, CARQUET_COMPRESSION_SNAPPY,
        CARQUET_COMPRESSION_GZIP, CARQUET_COMPRESSION_ZSTD,
        CARQUET_COMPRESSION_LZ4_RAW, CARQUET_COMPRESSION_LZ4,
        CARQUET_COMPRESSION_BROTLI
    };
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = codecs[h2 % (sizeof(codecs) / sizeof(codecs[0]))];
    if (opts.compression == CARQUET_COMPRESSION_BROTLI) opts.compression_level = 1;
    opts.data_page_version = (h2 & 0x08) ? 2 : 1;
    opts.write_page_index = (h2 & 0x10) != 0;
    opts.write_arrow_schema = (h2 & 0x20) != 0;
    opts.async_io = (h2 & 0x40) != 0;
    if (h3 & 1) opts.page_size = 64 + (h3 >> 1) * 16;
    if (h3 & 2) opts.max_rows_per_page = 1 + (h3 >> 2);
    if ((h3 & 0x0C) == 0x0C) opts.row_group_size = 256 + h3 * 8;   /* many row groups */

    carquet_error_t err = CARQUET_ERROR_INIT;
    struct ArrowSchema* imp_s;
    struct ArrowArray* imp_a;
    build_top(&sp, vals, &imp_s, &imp_a);
    imp_a->release(imp_a);
    free(imp_a);
    carquet_schema_t* schema = NULL;
    carquet_status_t st = carquet_arrow_import_schema(imp_s, &schema, &err);
    free(imp_s);
    REQUIRE(st == CARQUET_OK && schema, "import_schema");

    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    REQUIRE(w, "create_buffer");

    int32_t nleaves = carquet_schema_num_columns(schema);
    bool bloom_on[MAX_COLS * 2] = {0};
    for (int32_t l = 0; l < nleaves && l < MAX_COLS * 2; l++) {
        carquet_physical_type_t t = carquet_schema_column_type(schema, l);
        uint8_t e = enc_sel[l];
        carquet_encoding_t enc = CARQUET_ENCODING_PLAIN;
        switch (e & 7) {
            case 1: enc = CARQUET_ENCODING_RLE_DICTIONARY; break;
            case 2:
                if (t == CARQUET_PHYSICAL_INT32 || t == CARQUET_PHYSICAL_INT64)
                    enc = CARQUET_ENCODING_DELTA_BINARY_PACKED;
                else if (t == CARQUET_PHYSICAL_BYTE_ARRAY)
                    enc = CARQUET_ENCODING_DELTA_LENGTH_BYTE_ARRAY;
                break;
            case 3:
                if (t == CARQUET_PHYSICAL_BYTE_ARRAY || t == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
                    enc = CARQUET_ENCODING_DELTA_BYTE_ARRAY;
                break;
            case 4:
                if (t == CARQUET_PHYSICAL_FLOAT || t == CARQUET_PHYSICAL_DOUBLE ||
                    t == CARQUET_PHYSICAL_INT32 || t == CARQUET_PHYSICAL_INT64 ||
                    t == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
                    enc = CARQUET_ENCODING_BYTE_STREAM_SPLIT;
                break;
            case 5:
                if (t == CARQUET_PHYSICAL_BOOLEAN) enc = CARQUET_ENCODING_RLE;
                break;
            default: break;
        }
        if (enc != CARQUET_ENCODING_PLAIN)
            (void)carquet_writer_set_column_encoding(w, l, enc);
        if (e & 0x10) {
            bloom_on[l] = carquet_writer_set_column_bloom_filter(w, l, true) == CARQUET_OK;
        }
    }

    struct ArrowSchema* wr_s;
    struct ArrowArray* wr_a;
    build_top(&sp, vals, &wr_s, &wr_a);
    st = carquet_writer_write_arrow(w, wr_a, wr_s, &err);
    free(wr_a);
    free(wr_s);
    REQUIRE(st == CARQUET_OK, "write_arrow");
    st = carquet_writer_close(w);
    REQUIRE(st == CARQUET_OK, "writer_close");
    void* buf = NULL;
    size_t bsz = 0;
    st = carquet_writer_get_buffer(w, &buf, &bsz);
    REQUIRE(st == CARQUET_OK && buf, "get_buffer");

    /* Source copy for comparison. */
    struct ArrowSchema* src_s;
    struct ArrowArray* src_a;
    build_top(&sp, vals, &src_s, &src_a);

    carquet_reader_t* r = carquet_reader_open_buffer(buf, bsz, NULL, &err);
    REQUIRE(r, "reader_open_buffer");
    REQUIRE(carquet_reader_num_rows(r) == sp.nrows, "row count");

    int64_t row = 0;
    int32_t nrg = carquet_reader_num_row_groups(r);
    for (int32_t rg = 0; rg < nrg; rg++) {
        struct ArrowSchema rs;
        struct ArrowArray ra;
        memset(&rs, 0, sizeof(rs));
        memset(&ra, 0, sizeof(ra));
        st = carquet_reader_read_arrow(r, rg, &rs, &ra, &err);
        REQUIRE(st == CARQUET_OK, "read_arrow");
        REQUIRE(rs.n_children == sp.ncols && ra.n_children == sp.ncols, "read arity");
        REQUIRE(row + ra.length <= sp.nrows, "row overflow");
        for (int i = 0; i < sp.ncols; i++) {
            g_fail_path = g_names[i];
            for (int64_t j = 0; j < ra.length; j++) {
                cmp_value(src_s->children[i], src_a->children[i], src_a->offset + row + j,
                          rs.children[i], ra.children[i], ra.offset + j);
            }
        }
        /* Bloom: leaf index of column i is the sum of leaves before it. */
        int32_t leaf = 0;
        for (int i = 0; i < sp.ncols; i++) {
            int k = sp.kinds[i];
            if ((k == 0 || k == 1 || k == 5) && leaf < MAX_COLS * 2 && bloom_on[leaf])
                check_bloom(r, rg, leaf, k, src_s->children[i], src_a->children[i],
                            src_a->offset + row, ra.length);
            leaf += kind_leaves(k);
        }
        row += ra.length;
        ra.release(&ra);
        rs.release(&rs);
    }
    REQUIRE(row == sp.nrows, "total rows read");

    carquet_reader_close(r);
    src_a->release(src_a);
    free(src_a);
    src_s->release(src_s);
    free(src_s);
    free(buf);
    carquet_schema_free(schema);
    return 0;
}

#ifdef AFL_MAIN
#include <sys/stat.h>
int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <input_file>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 1; }
    struct stat st; fstat(fileno(f), &st);
    uint8_t* d = malloc((size_t)st.st_size);
    if (!d) { fclose(f); return 1; }
    fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, (size_t)st.st_size);
    free(d); return r;
}
#endif
