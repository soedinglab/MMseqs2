/**
 * @file arrow_restore.c
 * @brief Restore the Arrow types recorded in a file's "ARROW:schema" entry.
 *
 * Writers that start from Arrow data (PyArrow, Arrow C++, carquet's own
 * write_arrow) store the original Arrow schema in the footer, base64 of an
 * IPC Schema message, because Parquet cannot express every Arrow type: time
 * zones, large / view / dictionary types, durations, date64, second units,
 * the null type, fixed-size lists, decimal widths, child field names and
 * field metadata. pyarrow.parquet reads columns back with those types; this
 * module does the same for the Arrow stream:
 *
 *   1. carquet_arrow_stored_schema() decodes the blob into an ArrowSchema
 *      (bounds-checked: the blob comes from the file);
 *   2. carquet_restore_plan() compares it, field by field, with the schema
 *      the export produces and builds the restored schema plus a plan of the
 *      value conversions it needs; a node whose stored type is not reachable
 *      from the exported one keeps the exported type;
 *   3. carquet_restore_apply() converts each batch.
 *
 * The footer's key/value metadata (minus ARROW:schema, plus the stored
 * schema's own) becomes the schema metadata, as with pyarrow.
 */

#include <carquet/carquet.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fb_read.h"
#include "reader_internal.h"

/* ============================================================================
 * Small helpers
 * ============================================================================
 */

static char* dupn(const char* s, size_t n) {
    char* o = (char*)malloc(n + 1);
    if (!o) return NULL;
    if (n) memcpy(o, s, n);
    o[n] = '\0';
    return o;
}

static void schema_release(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format);
    free((void*)s->name);
    free((void*)s->metadata);
    for (int64_t i = 0; i < s->n_children; i++) {
        if (s->children[i]) {
            if (s->children[i]->release) s->children[i]->release(s->children[i]);
            free(s->children[i]);
        }
    }
    free(s->children);
    if (s->dictionary) {
        if (s->dictionary->release) s->dictionary->release(s->dictionary);
        free(s->dictionary);
    }
    s->release = NULL;
}

/* Arrow metadata blob builder: int32 count, then int32 length + bytes for
 * each key and value. */
typedef struct {
    char* buf;
    size_t len, cap;
    int32_t n;
} md_t;

static int md_put(md_t* m, const void* p, size_t n) {
    if (m->len + n > m->cap) {
        size_t cap = m->cap ? m->cap * 2 : 64;
        while (cap < m->len + n) cap *= 2;
        char* b = (char*)realloc(m->buf, cap);
        if (!b) return 0;
        m->buf = b;
        m->cap = cap;
    }
    memcpy(m->buf + m->len, p, n);
    m->len += n;
    return 1;
}

static int md_add(md_t* m, const char* k, size_t kn, const char* v, size_t vn) {
    if (m->len == 0) {
        int32_t zero = 0;
        if (!md_put(m, &zero, 4)) return 0;
    }
    int32_t a = (int32_t)kn, b = (int32_t)vn;
    if (!md_put(m, &a, 4) || !md_put(m, k, kn) || !md_put(m, &b, 4) || !md_put(m, v, vn)) return 0;
    m->n++;
    memcpy(m->buf, &m->n, 4);
    return 1;
}

/* Find key `k` in a metadata blob. */
static int md_has(const char* md, const char* k, size_t kn) {
    if (!md) return 0;
    int32_t n;
    memcpy(&n, md, 4);
    size_t pos = 4;
    for (int32_t i = 0; i < n; i++) {
        int32_t a, b;
        memcpy(&a, md + pos, 4);
        int hit = (size_t)a == kn && memcmp(md + pos + 4, k, kn) == 0;
        pos += 4 + (size_t)a;
        memcpy(&b, md + pos, 4);
        pos += 4 + (size_t)b;
        if (hit) return 1;
    }
    return 0;
}

/* Append every pair of `src` whose key `dst` lacks. */
static int md_merge(md_t* dst, const char* src) {
    if (!src) return 1;
    int32_t n;
    memcpy(&n, src, 4);
    size_t pos = 4;
    for (int32_t i = 0; i < n; i++) {
        int32_t a, b;
        memcpy(&a, src + pos, 4);
        const char* k = src + pos + 4;
        pos += 4 + (size_t)a;
        memcpy(&b, src + pos, 4);
        const char* v = src + pos + 4;
        pos += 4 + (size_t)b;
        if (!md_has(dst->buf, k, (size_t)a) && !md_add(dst, k, (size_t)a, v, (size_t)b)) return 0;
    }
    return 1;
}

static size_t md_size(const char* md) {
    if (!md) return 0;
    int32_t n;
    memcpy(&n, md, 4);
    size_t pos = 4;
    for (int32_t i = 0; i < 2 * n; i++) {
        int32_t a;
        memcpy(&a, md + pos, 4);
        pos += 4 + (size_t)a;
    }
    return pos;
}

static char* md_dup(const char* md) {
    if (!md) return NULL;
    size_t n = md_size(md);
    char* o = (char*)malloc(n);
    if (o) memcpy(o, md, n);
    return o;
}

/* ============================================================================
 * IPC Schema -> ArrowSchema
 * ============================================================================
 */

enum {
    T_NULL = 1, T_INT = 2, T_FLOAT = 3, T_BINARY = 4, T_UTF8 = 5, T_BOOL = 6, T_DECIMAL = 7,
    T_DATE = 8, T_TIME = 9, T_TIMESTAMP = 10, T_INTERVAL = 11, T_LIST = 12, T_STRUCT = 13,
    T_UNION = 14, T_FSB = 15, T_FSL = 16, T_MAP = 17, T_DURATION = 18, T_LARGEBINARY = 19,
    T_LARGEUTF8 = 20, T_LARGELIST = 21, T_REE = 22, T_BINARYVIEW = 23, T_UTF8VIEW = 24,
    T_LISTVIEW = 25, T_LARGELISTVIEW = 26
};

typedef struct {
    fbr r;
    size_t budget;   /* bytes of strings / metadata we still accept to copy */
} dec_t;

static int fb_scalar_i32(const fbr* r, size_t tbl, int slot, int32_t def, int32_t* out) {
    size_t f;
    if (!fbr_field(r, tbl, slot, &f)) return 0;
    if (f == 0) { *out = def; return 1; }
    return fbr_i32(r, f, out);
}

static int fb_scalar_i16(const fbr* r, size_t tbl, int slot, int32_t def, int32_t* out) {
    size_t f;
    if (!fbr_field(r, tbl, slot, &f)) return 0;
    if (f == 0) { *out = def; return 1; }
    uint16_t v;
    if (!fbr_u16(r, f, &v)) return 0;
    *out = (int16_t)v;
    return 1;
}

static int fb_scalar_bool(const fbr* r, size_t tbl, int slot, int* out) {
    size_t f;
    if (!fbr_field(r, tbl, slot, &f)) return 0;
    if (f == 0) { *out = 0; return 1; }
    uint8_t v;
    if (!fbr_u8(r, f, &v)) return 0;
    *out = v != 0;
    return 1;
}

static const char UNIT_CHAR[4] = {'s', 'm', 'u', 'n'};

static const char* int_format(int32_t bits, int is_signed) {
    switch (bits) {
    case 8: return is_signed ? "c" : "C";
    case 16: return is_signed ? "s" : "S";
    case 32: return is_signed ? "i" : "I";
    case 64: return is_signed ? "l" : "L";
    default: return NULL;
    }
}

/* Format string of an IPC Type table; "?" for types this module does not
 * restore (unions, run-end encoding...). */
static int type_format(dec_t* d, uint8_t tag, size_t t, char* fmt, size_t fs, int* keys_sorted) {
    const fbr* r = &d->r;
    int32_t a = 0, b = 0, c = 0;
    *keys_sorted = 0;
    switch (tag) {
    case T_NULL: snprintf(fmt, fs, "n"); return 1;
    case T_INT: {
        int sg;
        if (!t || !fb_scalar_i32(r, t, 0, 0, &a) || !fb_scalar_bool(r, t, 1, &sg)) return 0;
        const char* f = int_format(a, sg);
        if (!f) return 0;
        snprintf(fmt, fs, "%s", f);
        return 1;
    }
    case T_FLOAT:
        if (!t || !fb_scalar_i16(r, t, 0, 0, &a) || a < 0 || a > 2) return 0;
        snprintf(fmt, fs, "%s", a == 0 ? "e" : a == 1 ? "f" : "g");
        return 1;
    case T_BINARY: snprintf(fmt, fs, "z"); return 1;
    case T_UTF8: snprintf(fmt, fs, "u"); return 1;
    case T_BOOL: snprintf(fmt, fs, "b"); return 1;
    case T_DECIMAL:
        if (!t || !fb_scalar_i32(r, t, 0, 0, &a) || !fb_scalar_i32(r, t, 1, 0, &b) ||
            !fb_scalar_i32(r, t, 2, 128, &c)) return 0;
        if (c == 128) snprintf(fmt, fs, "d:%d,%d", (int)a, (int)b);
        else snprintf(fmt, fs, "d:%d,%d,%d", (int)a, (int)b, (int)c);
        return 1;
    case T_DATE:
        if (!fb_scalar_i16(r, t, 0, 1, &a) || a < 0 || a > 1) return 0;
        snprintf(fmt, fs, "%s", a == 0 ? "tdD" : "tdm");
        return 1;
    case T_TIME:
        if (!fb_scalar_i16(r, t, 0, 1, &a) || !fb_scalar_i32(r, t, 1, 32, &b) || a < 0 || a > 3) return 0;
        snprintf(fmt, fs, "tt%c", UNIT_CHAR[a]);
        return 1;
    case T_TIMESTAMP: {
        if (!fb_scalar_i16(r, t, 0, 0, &a) || a < 0 || a > 3) return 0;
        const char* tz = NULL;
        size_t tzn = 0;
        if (t) fbr_field_string_view(r, t, 1, &tz, &tzn);
        if (tzn > 200) return 0;
        snprintf(fmt, fs, "ts%c:%.*s", UNIT_CHAR[a], (int)tzn, tz ? tz : "");
        return 1;
    }
    case T_INTERVAL:
        if (!fb_scalar_i16(r, t, 0, 0, &a) || a < 0 || a > 2) return 0;
        snprintf(fmt, fs, "%s", a == 0 ? "tiM" : a == 1 ? "tiD" : "tin");
        return 1;
    case T_LIST: snprintf(fmt, fs, "+l"); return 1;
    case T_STRUCT: snprintf(fmt, fs, "+s"); return 1;
    case T_FSB:
        if (!t || !fb_scalar_i32(r, t, 0, 0, &a) || a <= 0) return 0;
        snprintf(fmt, fs, "w:%d", (int)a);
        return 1;
    case T_FSL:
        if (!t || !fb_scalar_i32(r, t, 0, 0, &a) || a <= 0) return 0;
        snprintf(fmt, fs, "+w:%d", (int)a);
        return 1;
    case T_MAP: {
        int ks = 0;
        if (t && !fb_scalar_bool(r, t, 0, &ks)) return 0;
        *keys_sorted = ks;
        snprintf(fmt, fs, "+m");
        return 1;
    }
    case T_DURATION:
        if (!fb_scalar_i16(r, t, 0, 1, &a) || a < 0 || a > 3) return 0;
        snprintf(fmt, fs, "tD%c", UNIT_CHAR[a]);
        return 1;
    case T_LARGEBINARY: snprintf(fmt, fs, "Z"); return 1;
    case T_LARGEUTF8: snprintf(fmt, fs, "U"); return 1;
    case T_LARGELIST: snprintf(fmt, fs, "+L"); return 1;
    case T_BINARYVIEW: snprintf(fmt, fs, "vz"); return 1;
    case T_UTF8VIEW: snprintf(fmt, fs, "vu"); return 1;
    default:
        snprintf(fmt, fs, "?");
        return 1;
    }
}

/* KeyValue vector -> metadata blob (NULL when empty). */
static int decode_kv(dec_t* d, size_t tbl, int slot, char** out) {
    *out = NULL;
    size_t data;
    uint32_t n;
    if (!fbr_field_vector(&d->r, tbl, slot, 4, &data, &n) || n == 0) return 1;
    md_t m = {0};
    for (uint32_t i = 0; i < n; i++) {
        size_t kv;
        const char* k = NULL;
        const char* v = "";
        size_t kn = 0, vn = 0;
        if (!fbr_indirect(&d->r, data + (size_t)i * 4, &kv)) continue;
        if (!fbr_field_string_view(&d->r, kv, 0, &k, &kn)) continue;
        fbr_field_string_view(&d->r, kv, 1, &v, &vn);
        if (kn + vn + 8 > d->budget) { free(m.buf); return 0; }
        d->budget -= kn + vn + 8;
        if (!md_add(&m, k, kn, v, vn)) { free(m.buf); return 0; }
    }
    *out = m.buf;
    return 1;
}

static int decode_field(dec_t* d, size_t fld, struct ArrowSchema* out, int depth) {
    memset(out, 0, sizeof(*out));
    out->release = schema_release;
    if (depth > 64) return 0;
    const fbr* r = &d->r;

    const char* name = NULL;
    size_t nl = 0;
    fbr_field_string_view(r, fld, 0, &name, &nl);
    if (nl > d->budget) return 0;
    d->budget -= nl;
    out->name = dupn(name ? name : "", nl);
    int nullable = 0;
    if (!out->name || !fb_scalar_bool(r, fld, 1, &nullable)) return 0;
    if (nullable) out->flags |= ARROW_FLAG_NULLABLE;

    size_t f;
    uint8_t tag = 0;
    if (!fbr_field(r, fld, 2, &f)) return 0;
    if (f && !fbr_u8(r, f, &tag)) return 0;
    size_t ttbl = 0;
    if (!fbr_field(r, fld, 3, &f)) return 0;
    if (f && !fbr_indirect(r, f, &ttbl)) return 0;
    char fmt[256];
    int keys_sorted = 0;
    if (!type_format(d, tag, ttbl, fmt, sizeof(fmt), &keys_sorted)) return 0;
    if (keys_sorted) out->flags |= ARROW_FLAG_MAP_KEYS_SORTED;

    /* Children. */
    size_t cdata;
    uint32_t nc = 0;
    if (fbr_field_vector(r, fld, 5, 4, &cdata, &nc) && nc > 0) {
        if ((size_t)nc * sizeof(struct ArrowSchema) > d->budget) return 0;
        d->budget -= (size_t)nc * sizeof(struct ArrowSchema);
        out->children = (struct ArrowSchema**)calloc(nc, sizeof(void*));
        if (!out->children) return 0;
        out->n_children = nc;
        for (uint32_t i = 0; i < nc; i++) {
            size_t ct;
            out->children[i] = (struct ArrowSchema*)calloc(1, sizeof(struct ArrowSchema));
            if (!out->children[i] || !fbr_indirect(r, cdata + (size_t)i * 4, &ct)) return 0;
            if (!decode_field(d, ct, out->children[i], depth + 1)) return 0;
        }
    }
    if (!decode_kv(d, fld, 6, (char**)&out->metadata)) return 0;

    /* Dictionary encoding: the node becomes the index type, the dictionary
     * child carries the value type. */
    size_t dict = 0;
    if (!fbr_field(r, fld, 4, &f)) return 0;
    if (f && !fbr_indirect(r, f, &dict)) return 0;
    if (dict) {
        int32_t bits = 32;
        int sg = 1, ordered = 0;
        size_t it = 0;
        if (!fbr_field(r, dict, 1, &f)) return 0;
        if (f && !fbr_indirect(r, f, &it)) return 0;
        if (it && (!fb_scalar_i32(r, it, 0, 32, &bits) || !fb_scalar_bool(r, it, 1, &sg))) return 0;
        if (!it) sg = 1;
        if (!fb_scalar_bool(r, dict, 2, &ordered)) return 0;
        const char* ifmt = int_format(bits, sg);
        if (!ifmt) return 0;
        out->dictionary = (struct ArrowSchema*)calloc(1, sizeof(struct ArrowSchema));
        if (!out->dictionary) return 0;
        out->dictionary->release = schema_release;
        out->dictionary->format = dupn(fmt, strlen(fmt));
        out->dictionary->name = dupn("", 0);
        if (!out->dictionary->format || !out->dictionary->name) return 0;
        out->dictionary->flags = ARROW_FLAG_NULLABLE;
        /* Children belong to the value type. */
        out->dictionary->children = out->children;
        out->dictionary->n_children = out->n_children;
        out->children = NULL;
        out->n_children = 0;
        if (ordered) out->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
        out->format = dupn(ifmt, strlen(ifmt));
    } else {
        out->format = dupn(fmt, strlen(fmt));
    }
    return out->format != NULL;
}

carquet_status_t carquet_arrow_stored_schema(const char* b64, struct ArrowSchema* out) {
    memset(out, 0, sizeof(*out));
    if (!b64) return CARQUET_ERROR_INVALID_ARGUMENT;
    size_t raw_len = 0;
    uint8_t* raw = base64_decode(b64, &raw_len);
    if (!raw) return CARQUET_ERROR_OUT_OF_MEMORY;
    size_t start;
    if (raw_len >= 8 && raw[0] == 0xFF && raw[1] == 0xFF && raw[2] == 0xFF && raw[3] == 0xFF) start = 8;
    else if (raw_len >= 4) start = 4;
    else { carquet_mem_free(raw); return CARQUET_ERROR_INVALID_METADATA; }

    dec_t d;
    d.r.buf = raw + start;
    d.r.len = raw_len - start;
    d.budget = d.r.len * 8 + 65536;

    carquet_status_t st = CARQUET_ERROR_INVALID_METADATA;
    size_t msg, hdr, schema, fields;
    uint32_t nf = 0;
    out->release = schema_release;
    out->format = dupn("+s", 2);
    if (!out->format) { st = CARQUET_ERROR_OUT_OF_MEMORY; goto done; }
    if (!fbr_indirect(&d.r, 0, &msg) || !fbr_field(&d.r, msg, 2, &hdr) || hdr == 0 ||
        !fbr_indirect(&d.r, hdr, &schema)) goto done;
    if (fbr_field_vector(&d.r, schema, 1, 4, &fields, &nf) && nf > 0) {
        if (nf > 1u << 16) goto done;
        out->children = (struct ArrowSchema**)calloc(nf, sizeof(void*));
        if (!out->children) { st = CARQUET_ERROR_OUT_OF_MEMORY; goto done; }
        out->n_children = nf;
        for (uint32_t i = 0; i < nf; i++) {
            size_t ft;
            out->children[i] = (struct ArrowSchema*)calloc(1, sizeof(struct ArrowSchema));
            if (!out->children[i] || !fbr_indirect(&d.r, fields + (size_t)i * 4, &ft) ||
                !decode_field(&d, ft, out->children[i], 0)) goto done;
        }
    }
    if (!decode_kv(&d, schema, 2, (char**)&out->metadata)) goto done;
    st = CARQUET_OK;
done:
    carquet_mem_free(raw);
    if (st != CARQUET_OK) schema_release(out);
    return st;
}

/* ============================================================================
 * Plan
 * ============================================================================
 */

enum {
    OP_KEEP = 0,      /* same values; children may have their own plans */
    OP_OFFSETS64,     /* utf8/binary/list -> large: int32 -> int64 offsets */
    OP_VIEW,          /* utf8/binary -> string/binary view */
    OP_DATE64,        /* int32 days -> int64 milliseconds */
    OP_DIV1000_32,    /* time32[ms] -> time32[s] */
    OP_DIV1000_64,    /* timestamp[ms] -> timestamp[s] */
    OP_NULL,          /* all-null column -> null type */
    OP_FSL,           /* list -> fixed-size list */
    OP_DICT,          /* values -> dictionary of the index type */
    OP_DECIMAL        /* decimal width change */
};

struct carquet_restore_plan {
    int op;
    int32_t arg;          /* FSL size, DICT index width (bytes), decimal target width */
    int32_t arg2;         /* DICT: index signed; DECIMAL: source width */
    char value_kind;      /* DICT: 'v' variable (utf8/binary), 'V' large, 'f' fixed */
    int32_t value_width;  /* DICT fixed: bytes */
    int32_t nkids;
    struct carquet_restore_plan* kids;   /* per child (struct / list / map / fsl) */
    struct ArrowSchema raw_child;        /* FSL: the child's exported schema, to compact it */
};

static int is_var(const char* f) { return (f[0] == 'u' || f[0] == 'z') && f[1] == '\0'; }

static int fixed_bytes(const char* f) {
    if (f[1] == '\0') {
        switch (f[0]) {
        case 'c': case 'C': return 1;
        case 's': case 'S': case 'e': return 2;
        case 'i': case 'I': case 'f': return 4;
        case 'l': case 'L': case 'g': return 8;
        default: return 0;
        }
    }
    if (f[0] == 'w' && f[1] == ':') return atoi(f + 2);
    if (f[0] == 'd' && f[1] == ':') {
        int commas = 0;
        for (const char* p = f; *p; p++) commas += *p == ',';
        if (commas < 2) return 16;
        int bits = atoi(strrchr(f, ',') + 1);
        return bits / 8;
    }
    if (f[0] == 't') {
        if (f[1] == 'd') return f[2] == 'D' ? 4 : 8;
        if (f[1] == 't') return (f[2] == 's' || f[2] == 'm') ? 4 : 8;
        if (f[1] == 's' || f[1] == 'D') return 8;
    }
    return 0;
}

/* "d:p,s" (any width) -> p,s */
static int decimal_ps(const char* f, int* p, int* s) {
    if (f[0] != 'd' || f[1] != ':') return 0;
    return sscanf(f + 2, "%d,%d", p, s) == 2;
}

static int set_str(const char** slot, const char* v) {
    char* c = dupn(v, strlen(v));
    if (!c) return 0;
    free((void*)*slot);
    *slot = c;
    return 1;
}

/* Plan the conversion of node `raw` (in `out`, a copy being edited) to
 * `stored`. Returns 0 only when out of memory. */
static int plan_node(struct ArrowSchema* out, const struct ArrowSchema* stored,
                     carquet_restore_plan_t* p, int depth) {
    memset(p, 0, sizeof(*p));
    if (depth > 64) return 1;
    /* A copy: set_str() below replaces (frees) out->format. */
    char rf[256];
    snprintf(rf, sizeof(rf), "%s", out->format ? out->format : "");
    const char* sf = stored->format;

    /* Metadata and flags first: they apply whatever the conversion. */
    if (stored->metadata) {
        free((void*)out->metadata);
        out->metadata = md_dup(stored->metadata);
        if (!out->metadata) return 0;
    }
    out->flags = (out->flags & ARROW_FLAG_NULLABLE) |
                 (stored->flags & (ARROW_FLAG_MAP_KEYS_SORTED | ARROW_FLAG_DICTIONARY_ORDERED));

    if (stored->dictionary) {
        const char* vf = stored->dictionary->format;
        int ib = fixed_bytes(sf);
        int fixed = fixed_bytes(rf);
        bool var_ok = (is_var(rf) && (strcmp(vf, rf) == 0 || (vf[0] == rf[0] - 32 && vf[1] == '\0')));
        bool fixed_ok = fixed > 0 && strcmp(vf, rf) == 0 && rf[0] != '+';
        if (ib >= 1 && ib <= 8 && (var_ok || fixed_ok)) {
            p->op = OP_DICT;
            p->arg = ib;
            p->arg2 = sf[0] >= 'a';   /* lower-case formats are signed */
            p->value_kind = var_ok ? ((vf[0] == 'U' || vf[0] == 'Z') ? 'V' : 'v') : 'f';
            p->value_width = fixed;
            struct ArrowSchema* d = (struct ArrowSchema*)calloc(1, sizeof(*d));
            if (!d) return 0;
            d->release = schema_release;
            d->format = dupn(vf, strlen(vf));
            d->name = dupn("", 0);
            d->flags = ARROW_FLAG_NULLABLE;
            if (!d->format || !d->name || !set_str(&out->format, sf)) { schema_release(d); free(d); return 0; }
            out->dictionary = d;
            return 1;
        }
        return 1;   /* keep the plain values */
    }

    if (strcmp(rf, sf) == 0) {
        /* Same type: recurse into nested children. */
        if (rf[0] == '+' && out->n_children == stored->n_children && out->n_children > 0) {
            p->kids = (carquet_restore_plan_t*)calloc((size_t)out->n_children, sizeof(*p->kids));
            if (!p->kids) return 0;
            p->nkids = (int32_t)out->n_children;
            /* Struct children are named by the stored schema; list / map
             * children keep their Parquet names (element, key_value...), as
             * pyarrow reads them. */
            bool rename = strcmp(rf, "+s") == 0;
            for (int64_t i = 0; i < out->n_children; i++) {
                if (rename && stored->children[i]->name &&
                    !set_str(&out->children[i]->name, stored->children[i]->name)) return 0;
                if (!plan_node(out->children[i], stored->children[i], &p->kids[i], depth + 1)) return 0;
            }
        }
        return 1;
    }

    int op = -1;
    if (is_var(rf) && sf[1] == '\0' && sf[0] == rf[0] - 32) op = OP_OFFSETS64;          /* u->U, z->Z */
    else if (is_var(rf) && sf[0] == 'v' && sf[1] == rf[0] && sf[2] == '\0') op = OP_VIEW;  /* u->vu, z->vz */
    else if (strncmp(rf, "ts", 2) == 0 && strncmp(sf, "ts", 2) == 0 && rf[3] == ':' && sf[3] == ':') {
        /* The stored time zone, on the unit Parquet holds: pyarrow does not
         * restore timestamp[s] (nor date64 / time32[s]), it reads them in the
         * coarsest unit Parquet has. */
        char f[256];
        snprintf(f, sizeof(f), "ts%c:%s", rf[2], sf + 4);
        if (!set_str(&out->format, f)) return 0;
        return 1;
    }
    else if (strcmp(rf, "l") == 0 && strncmp(sf, "tD", 2) == 0 && sf[3] == '\0') op = OP_KEEP;          /* duration */
    else if (strcmp(sf, "n") == 0 && rf[0] != '+') op = OP_NULL;
    else if (strcmp(rf, "+l") == 0 && strcmp(sf, "+L") == 0) op = OP_OFFSETS64;
    else if (strcmp(rf, "+l") == 0 && strncmp(sf, "+w:", 3) == 0 && atoi(sf + 3) > 0) op = OP_FSL;
    else {
        int p1, s1, p2, s2;
        if (decimal_ps(rf, &p1, &s1) && decimal_ps(sf, &p2, &s2) && p1 == p2 && s1 == s2) {
            /* The stored width must be one its precision can have (the blob
             * is untrusted): 9 / 18 / 38 / 76 digits. */
            int w = fixed_bytes(sf);
            int max_p = w == 4 ? 9 : w == 8 ? 18 : w == 16 ? 38 : w == 32 ? 76 : 0;
            if (max_p && p2 >= 1 && p2 <= max_p) {
                op = OP_DECIMAL;
                p->arg = w;
                p->arg2 = fixed_bytes(rf);
            }
        }
    }
    if (op < 0) return 1;   /* not reachable: keep the exported type */
    p->op = op;
    if (op == OP_FSL) {
        p->arg = atoi(sf + 3);
        if (out->n_children != 1 || carquet_arrow_schema_copy(out->children[0], &p->raw_child) != CARQUET_OK) {
            p->op = OP_KEEP;
            return 1;
        }
    }
    if (!set_str(&out->format, sf)) return 0;
    if (op == OP_NULL) out->flags |= ARROW_FLAG_NULLABLE;
    /* Lists keep restoring below. */
    if ((op == OP_OFFSETS64 || op == OP_FSL) && rf[0] == '+' && out->n_children == 1 && stored->n_children == 1) {
        p->kids = (carquet_restore_plan_t*)calloc(1, sizeof(*p->kids));
        if (!p->kids) return 0;
        p->nkids = 1;
        if (!plan_node(out->children[0], stored->children[0], &p->kids[0], depth + 1)) return 0;
    }
    return 1;
}

static void plan_free_kids(carquet_restore_plan_t* p) {
    for (int32_t i = 0; i < p->nkids; i++) plan_free_kids(&p->kids[i]);
    free(p->kids);
    if (p->raw_child.release) p->raw_child.release(&p->raw_child);
}

void carquet_restore_free(carquet_restore_plan_t* p) {
    if (!p) return;
    plan_free_kids(p);
    free(p);
}

static bool plan_is_noop(const carquet_restore_plan_t* p) {
    if (p->op != OP_KEEP) return false;
    for (int32_t i = 0; i < p->nkids; i++) if (!plan_is_noop(&p->kids[i])) return false;
    return true;
}

carquet_status_t carquet_restore_plan(const carquet_reader_t* reader, bool restore,
                                      struct ArrowSchema* schema, carquet_restore_plan_t** plan) {
    *plan = NULL;
    /* Footer key/value metadata -> schema metadata. */
    md_t m = {0};
    const char* b64 = NULL;
    int32_t nkv = carquet_reader_num_metadata(reader);
    for (int32_t i = 0; i < nkv; i++) {
        const char* k = NULL;
        const char* v = NULL;
        if (carquet_reader_get_metadata(reader, i, &k, &v) != CARQUET_OK || !k) continue;
        if (strcmp(k, "ARROW:schema") == 0) { b64 = v; continue; }
        if (!md_has(m.buf, k, strlen(k)) && !md_add(&m, k, strlen(k), v ? v : "", v ? strlen(v) : 0)) {
            free(m.buf);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
    }

    struct ArrowSchema stored;
    memset(&stored, 0, sizeof(stored));
    bool have_stored = restore && b64 && carquet_arrow_stored_schema(b64, &stored) == CARQUET_OK;
    if (have_stored && !md_merge(&m, stored.metadata)) {
        free(m.buf);
        schema_release(&stored);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    free((void*)schema->metadata);
    schema->metadata = m.buf;
    if (!have_stored) return CARQUET_OK;

    carquet_restore_plan_t* p = (carquet_restore_plan_t*)calloc(1, sizeof(*p));
    if (!p) { schema_release(&stored); return CARQUET_ERROR_OUT_OF_MEMORY; }
    if (schema->n_children > 0) {
        p->kids = (carquet_restore_plan_t*)calloc((size_t)schema->n_children, sizeof(*p->kids));
        if (!p->kids) { free(p); schema_release(&stored); return CARQUET_ERROR_OUT_OF_MEMORY; }
        p->nkids = (int32_t)schema->n_children;
    }
    bool* used = (bool*)calloc((size_t)(stored.n_children > 0 ? stored.n_children : 1), sizeof(bool));
    carquet_status_t st = used ? CARQUET_OK : CARQUET_ERROR_OUT_OF_MEMORY;
    for (int64_t i = 0; st == CARQUET_OK && i < schema->n_children; i++) {
        const char* name = schema->children[i]->name ? schema->children[i]->name : "";
        for (int64_t j = 0; j < stored.n_children; j++) {
            const char* sn = stored.children[j]->name ? stored.children[j]->name : "";
            if (!used[j] && strcmp(sn, name) == 0) {
                used[j] = true;
                if (!plan_node(schema->children[i], stored.children[j], &p->kids[i], 0)) st = CARQUET_ERROR_OUT_OF_MEMORY;
                break;
            }
        }
    }
    free(used);
    schema_release(&stored);
    if (st != CARQUET_OK || plan_is_noop(p)) {
        carquet_restore_free(p);
        return st;
    }
    *plan = p;
    return CARQUET_OK;
}

/* ============================================================================
 * Apply
 * ============================================================================
 */

static void owned_release(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) {
        for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]);
        free(a->buffers);
    }
    for (int64_t i = 0; i < a->n_children; i++) {
        if (a->children[i]) {
            if (a->children[i]->release) a->children[i]->release(a->children[i]);
            free(a->children[i]);
        }
    }
    free(a->children);
    if (a->dictionary) {
        if (a->dictionary->release) a->dictionary->release(a->dictionary);
        free(a->dictionary);
    }
    a->release = NULL;
}

/* Child slots a fixed-size list restore may create for null lists stored
 * empty (see OP_FSL in apply_node). */
#define FSL_MAX_FILLER ((int64_t)1 << 27)

static bool bit(const void* b, int64_t i) { return (((const uint8_t*)b)[i >> 3] >> (i & 7)) & 1; }

/* Copy of the validity of rows [off, off + n) re-based to bit 0 (NULL if none). */
static carquet_status_t copy_validity(const struct ArrowArray* a, uint8_t** out) {
    *out = NULL;
    const void* v = a->n_buffers > 0 ? a->buffers[0] : NULL;
    if (!v || a->null_count == 0) return CARQUET_OK;
    uint8_t* o = (uint8_t*)calloc((size_t)((a->length + 7) / 8 + 1), 1);
    if (!o) return CARQUET_ERROR_OUT_OF_MEMORY;
    for (int64_t i = 0; i < a->length; i++) {
        if (bit(v, a->offset + i)) o[i >> 3] |= (uint8_t)(1u << (i & 7));
    }
    *out = o;
    return CARQUET_OK;
}

/* A new owned node of `len` rows with `nbuf` buffers (all NULL). */
static struct ArrowArray* new_node(int64_t len, int64_t null_count, int64_t nbuf) {
    struct ArrowArray* a = (struct ArrowArray*)calloc(1, sizeof(*a));
    if (!a) return NULL;
    a->length = len;
    a->null_count = null_count;
    a->n_buffers = nbuf;
    a->release = owned_release;
    if (nbuf > 0) {
        a->buffers = (const void**)calloc((size_t)nbuf, sizeof(void*));
        if (!a->buffers) { free(a); return NULL; }
    }
    return a;
}

static int64_t offset_at(const struct ArrowArray* a, int64_t i) {
    int32_t v;
    memcpy(&v, (const uint8_t*)a->buffers[1] + 4 * (size_t)i, 4);
    return v;
}

static void release_slot(struct ArrowArray** slot) {
    struct ArrowArray* old = *slot;
    if (old) {
        if (old->release) old->release(old);
        free(old);
    }
}

static carquet_status_t apply_node(const carquet_restore_plan_t* p, struct ArrowArray** slot);

static carquet_status_t apply_kids(const carquet_restore_plan_t* p, struct ArrowArray* a) {
    for (int32_t i = 0; i < p->nkids && i < a->n_children; i++) {
        carquet_status_t st = apply_node(&p->kids[i], &a->children[i]);
        if (st != CARQUET_OK) return st;
    }
    return CARQUET_OK;
}

/* Dictionary encoding. */
typedef struct { int64_t row; uint64_t hash; } dslot_t;

static uint64_t hash_bytes(const uint8_t* p, size_t n) {
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h ^ (n * 0x9E3779B97F4A7C15ULL);
}

static void value_at(const struct ArrowArray* a, const carquet_restore_plan_t* p, int64_t j,
                     const uint8_t** v, size_t* n) {
    if (p->value_kind == 'f') {
        *v = (const uint8_t*)a->buffers[1] + (size_t)j * (size_t)p->value_width;
        *n = (size_t)p->value_width;
    } else {
        int64_t s = offset_at(a, j), e = offset_at(a, j + 1);
        *v = (const uint8_t*)a->buffers[2] + s;
        *n = (size_t)(e - s);
    }
}

static carquet_status_t apply_dict(const carquet_restore_plan_t* p, struct ArrowArray** slot) {
    struct ArrowArray* a = *slot;
    int64_t n = a->length;
    const void* valid = a->n_buffers > 0 ? a->buffers[0] : NULL;
    bool has_nulls = valid && a->null_count != 0;
    size_t cap = 16;
    while (cap < (size_t)n * 2) cap <<= 1;
    dslot_t* table = (dslot_t*)malloc(cap * sizeof(dslot_t));
    int64_t* uniq = (int64_t*)malloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));   /* physical rows */
    uint8_t* idx = (uint8_t*)calloc((size_t)(n > 0 ? n : 1), (size_t)p->arg);
    if (!table || !uniq || !idx) { free(table); free(uniq); free(idx); return CARQUET_ERROR_OUT_OF_MEMORY; }
    for (size_t i = 0; i < cap; i++) table[i].row = -1;
    int64_t nu = 0;
    uint64_t limit = p->arg == 8 ? UINT64_MAX : (1ULL << (8 * p->arg - (p->arg2 ? 1 : 0)));
    for (int64_t i = 0; i < n; i++) {
        int64_t j = a->offset + i;
        if (has_nulls && !bit(valid, j)) continue;
        const uint8_t* v;
        size_t vn;
        value_at(a, p, j, &v, &vn);
        uint64_t h = hash_bytes(v, vn);
        size_t k = (size_t)h & (cap - 1);
        int64_t code = -1;
        while (table[k].row >= 0) {
            if (table[k].hash == h) {
                const uint8_t* w;
                size_t wn;
                value_at(a, p, uniq[table[k].row], &w, &wn);
                if (wn == vn && (vn == 0 || memcmp(v, w, vn) == 0)) { code = table[k].row; break; }
            }
            k = (k + 1) & (cap - 1);
        }
        if (code < 0) {
            if ((uint64_t)nu >= limit) {
                /* More distinct values than the stored index type counts: the
                 * stream's schema is fixed (index type + dictionary), so the
                 * batch cannot fall back to plain values under it. */
                free(table); free(uniq); free(idx);
                return CARQUET_ERROR_DECODE;
            }
            code = nu;
            uniq[nu++] = j;
            table[k].row = code;
            table[k].hash = h;
        }
        uint64_t c = (uint64_t)code;
        memcpy(idx + (size_t)i * (size_t)p->arg, &c, (size_t)p->arg);   /* little-endian truncation */
    }
    free(table);

    /* Dictionary values. */
    struct ArrowArray* dict;
    if (p->value_kind == 'f') {
        dict = new_node(nu, 0, 2);
        uint8_t* d = dict ? (uint8_t*)malloc((size_t)(nu > 0 ? nu : 1) * (size_t)p->value_width) : NULL;
        if (!dict || !d) { free(d); if (dict) owned_release(dict); free(dict); free(uniq); free(idx); return CARQUET_ERROR_OUT_OF_MEMORY; }
        for (int64_t u = 0; u < nu; u++) {
            const uint8_t* v; size_t vn;
            value_at(a, p, uniq[u], &v, &vn);
            memcpy(d + (size_t)u * vn, v, vn);
        }
        dict->buffers[1] = d;
    } else {
        bool large = p->value_kind == 'V';
        int64_t total = 0;
        for (int64_t u = 0; u < nu; u++) {
            const uint8_t* v; size_t vn;
            value_at(a, p, uniq[u], &v, &vn);
            total += (int64_t)vn;
        }
        dict = new_node(nu, 0, 3);
        size_t ow = large ? 8 : 4;
        uint8_t* offs = dict ? (uint8_t*)malloc((size_t)(nu + 1) * ow) : NULL;
        uint8_t* data = dict ? (uint8_t*)malloc((size_t)(total > 0 ? total : 1)) : NULL;
        if (!dict || !offs || !data || (!large && total > INT32_MAX)) {
            free(offs); free(data); if (dict) owned_release(dict); free(dict); free(uniq); free(idx);
            return dict && offs && data ? CARQUET_ERROR_DECODE : CARQUET_ERROR_OUT_OF_MEMORY;
        }
        int64_t at = 0;
        for (int64_t u = 0; u <= nu; u++) {
            if (large) memcpy(offs + 8 * (size_t)u, &at, 8);
            else { int32_t t = (int32_t)at; memcpy(offs + 4 * (size_t)u, &t, 4); }
            if (u == nu) break;
            const uint8_t* v; size_t vn;
            value_at(a, p, uniq[u], &v, &vn);
            if (vn) memcpy(data + at, v, vn);
            at += (int64_t)vn;
        }
        dict->buffers[1] = offs;
        dict->buffers[2] = data;
    }
    free(uniq);

    struct ArrowArray* o = new_node(n, a->null_count, 2);
    uint8_t* vb = NULL;
    if (!o || copy_validity(a, &vb) != CARQUET_OK) {
        if (o) owned_release(o);
        free(o); owned_release(dict); free(dict); free(idx);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    o->buffers[0] = vb;
    o->buffers[1] = idx;
    o->dictionary = dict;
    release_slot(slot);
    *slot = o;
    return CARQUET_OK;
}

static carquet_status_t apply_node(const carquet_restore_plan_t* p, struct ArrowArray** slot) {
    struct ArrowArray* a = *slot;
    if (!a) return CARQUET_ERROR_INTERNAL;
    int64_t n = a->length;
    switch (p->op) {
    case OP_KEEP:
        return apply_kids(p, a);

    case OP_DICT:
        return apply_dict(p, slot);

    case OP_NULL: {
        /* The stored type says the column holds no value: a value in it means
         * the file disagrees with its own ARROW:schema, and must not vanish. */
        if (a->n_buffers > 0 && n > 0) {
            const void* valid = a->buffers[0];
            if (!valid) return CARQUET_ERROR_DECODE;
            for (int64_t i = 0; i < n; i++) {
                if (bit(valid, a->offset + i)) return CARQUET_ERROR_DECODE;
            }
        }
        struct ArrowArray* o = new_node(n, n, 0);
        if (!o) return CARQUET_ERROR_OUT_OF_MEMORY;
        release_slot(slot);
        *slot = o;
        return CARQUET_OK;
    }

    case OP_DATE64: case OP_DIV1000_32: case OP_DIV1000_64: case OP_DECIMAL: {
        size_t in_w = p->op == OP_DATE64 || p->op == OP_DIV1000_32 ? 4
                    : p->op == OP_DIV1000_64 ? 8 : (size_t)p->arg2;
        size_t out_w = p->op == OP_DATE64 || p->op == OP_DIV1000_64 ? 8
                     : p->op == OP_DIV1000_32 ? 4 : (size_t)p->arg;
        struct ArrowArray* o = new_node(n, a->null_count, 2);
        uint8_t* vb = NULL;
        uint8_t* d = o ? (uint8_t*)malloc((size_t)(n > 0 ? n : 1) * out_w) : NULL;
        if (!o || !d || copy_validity(a, &vb) != CARQUET_OK) {
            free(d); if (o) owned_release(o); free(o);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        const uint8_t* src = (const uint8_t*)a->buffers[1] + (size_t)a->offset * in_w;
        for (int64_t i = 0; i < n; i++) {
            const uint8_t* s = src + (size_t)i * in_w;
            uint8_t* t = d + (size_t)i * out_w;
            if (p->op == OP_DATE64) {
                int32_t v; memcpy(&v, s, 4);
                int64_t ms = (int64_t)v * 86400000LL;
                memcpy(t, &ms, 8);
            } else if (p->op == OP_DIV1000_32) {
                int32_t v; memcpy(&v, s, 4);
                v = v >= 0 ? v / 1000 : -((-v + 999) / 1000);
                memcpy(t, &v, 4);
            } else if (p->op == OP_DIV1000_64) {
                int64_t v; memcpy(&v, s, 8);
                v = v >= 0 ? v / 1000 : -((-(v + 1)) / 1000) - 1;
                memcpy(t, &v, 8);
            } else {
                /* Two's complement, little-endian: sign-extend or truncate. */
                uint8_t fill = (s[in_w - 1] & 0x80) ? 0xFF : 0x00;
                for (size_t k = 0; k < out_w; k++) t[k] = k < in_w ? s[k] : fill;
            }
        }
        o->buffers[0] = vb;
        o->buffers[1] = d;
        release_slot(slot);
        *slot = o;
        return CARQUET_OK;
    }

    case OP_OFFSETS64:
    case OP_VIEW: {
        bool is_list = a->n_children == 1;
        if (p->op == OP_VIEW && is_list) return CARQUET_ERROR_INTERNAL;
        uint8_t* vb = NULL;
        if (copy_validity(a, &vb) != CARQUET_OK) return CARQUET_ERROR_OUT_OF_MEMORY;
        int64_t first = n > 0 ? offset_at(a, a->offset) : 0;
        int64_t last = n > 0 ? offset_at(a, a->offset + n) : 0;
        if (p->op == OP_OFFSETS64) {
            struct ArrowArray* o = new_node(n, a->null_count, is_list ? 2 : 3);
            int64_t* offs = o ? (int64_t*)malloc((size_t)(n + 1) * 8) : NULL;
            uint8_t* data = (!is_list && o) ? (uint8_t*)malloc((size_t)(last - first > 0 ? last - first : 1)) : NULL;
            if (!o || !offs || (!is_list && !data)) {
                free(vb); free(offs); free(data); if (o) owned_release(o); free(o);
                return CARQUET_ERROR_OUT_OF_MEMORY;
            }
            for (int64_t i = 0; i <= n; i++) offs[i] = offset_at(a, a->offset + i) - first;
            o->buffers[0] = vb;
            o->buffers[1] = offs;
            if (!is_list) {
                if (last > first) memcpy(data, (const uint8_t*)a->buffers[2] + first, (size_t)(last - first));
                o->buffers[2] = data;
            } else {
                o->children = (struct ArrowArray**)calloc(1, sizeof(void*));
                if (!o->children) { owned_release(o); free(o); return CARQUET_ERROR_OUT_OF_MEMORY; }
                o->n_children = 1;
                o->children[0] = a->children[0];       /* move the child, re-based */
                a->children[0] = NULL;
                o->children[0]->offset += first;
                o->children[0]->length = last - first;
                if (o->children[0]->null_count != 0) o->children[0]->null_count = -1;
                carquet_status_t st = p->nkids ? apply_node(&p->kids[0], &o->children[0]) : CARQUET_OK;
                if (st != CARQUET_OK) { owned_release(o); free(o); return st; }
            }
            release_slot(slot);
            *slot = o;
            return CARQUET_OK;
        }
        /* Views: 16 bytes each; values over 12 bytes point into one data buffer. */
        struct ArrowArray* o = new_node(n, a->null_count, 4);
        uint8_t* views = o ? (uint8_t*)calloc((size_t)(n > 0 ? n : 1), 16) : NULL;
        uint8_t* data = o ? (uint8_t*)malloc((size_t)(last - first > 0 ? last - first : 1)) : NULL;
        int64_t* sizes = o ? (int64_t*)malloc(sizeof(int64_t)) : NULL;
        if (!o || !views || !data || !sizes || last - first > INT32_MAX) {
            free(vb); free(views); free(data); free(sizes); if (o) owned_release(o); free(o);
            return CARQUET_ERROR_OUT_OF_MEMORY;
        }
        if (last > first) memcpy(data, (const uint8_t*)a->buffers[2] + first, (size_t)(last - first));
        for (int64_t i = 0; i < n; i++) {
            int64_t s = offset_at(a, a->offset + i), e = offset_at(a, a->offset + i + 1);
            int32_t len = (int32_t)(e - s);
            uint8_t* v = views + 16 * (size_t)i;
            memcpy(v, &len, 4);
            if (len <= 12) {
                if (len > 0) memcpy(v + 4, data + (s - first), (size_t)len);
            } else {
                int32_t buf = 0, off = (int32_t)(s - first);
                memcpy(v + 4, data + off, 4);
                memcpy(v + 8, &buf, 4);
                memcpy(v + 12, &off, 4);
            }
        }
        sizes[0] = last - first;
        o->buffers[0] = vb;
        o->buffers[1] = views;
        o->buffers[2] = data;
        o->buffers[3] = sizes;
        release_slot(slot);
        *slot = o;
        return CARQUET_OK;
    }

    case OP_FSL: {
        int64_t size = p->arg;
        const void* valid = a->n_buffers > 0 ? a->buffers[0] : NULL;
        bool has_nulls = valid && a->null_count != 0;
        struct ArrowArray* child = a->children[0];
        /* Every list must hold `size` elements; a null list may also be empty
         * (how Parquet stores it), and then needs `size` child slots made up. */
        int64_t filler = -1, nshort = 0;
        for (int64_t i = 0; i < n; i++) {
            int64_t s = offset_at(a, a->offset + i), e = offset_at(a, a->offset + i + 1);
            bool is_null = has_nulls && !bit(valid, a->offset + i);
            if (e - s == size) { if (filler < 0) filler = s; continue; }
            if (!(is_null && e == s)) return CARQUET_ERROR_DECODE;   /* not a fixed-size list */
            nshort++;
        }
        /* The size comes from the file's ARROW:schema and each made-up slot
         * costs memory no decoded value backs: bound them (a 500-byte file
         * with 3 null lists of size 2^30 asked for 25 GB). */
        if (nshort > 0 && nshort > FSL_MAX_FILLER / size) return CARQUET_ERROR_DECODE;
        struct ArrowArray* o = new_node(n, a->null_count, 1);
        uint8_t* vb = NULL;
        if (!o || copy_validity(a, &vb) != CARQUET_OK) { if (o) owned_release(o); free(o); return CARQUET_ERROR_OUT_OF_MEMORY; }
        o->buffers[0] = vb;
        o->children = (struct ArrowArray**)calloc(1, sizeof(void*));
        if (!o->children) { owned_release(o); free(o); return CARQUET_ERROR_OUT_OF_MEMORY; }
        o->n_children = 1;
        carquet_status_t st = CARQUET_OK;
        if (nshort == 0) {
            /* The lists are back to back: the child is the list's, re-based. */
            int64_t first = n > 0 ? offset_at(a, a->offset) : 0;
            o->children[0] = child;
            a->children[0] = NULL;
            child->offset += first;
            child->length = n * size;
            if (child->null_count != 0) child->null_count = -1;
        } else {
            int64_t total = n * size;
            struct ArrowArray* nc = (struct ArrowArray*)calloc(1, sizeof(*nc));
            if (!nc) { owned_release(o); free(o); return CARQUET_ERROR_OUT_OF_MEMORY; }
            if (filler < 0) {
                /* Only null lists: no child row to repeat, the child is all nulls. */
                st = carquet_arrow_null_array(&p->raw_child, total, nc);
            } else {
                /* Repeat a child row under the empty null lists (any value does there). */
                int64_t* rows = (int64_t*)malloc((size_t)total * sizeof(int64_t));
                if (!rows) { free(nc); owned_release(o); free(o); return CARQUET_ERROR_OUT_OF_MEMORY; }
                for (int64_t i = 0; i < n; i++) {
                    int64_t s = offset_at(a, a->offset + i), e = offset_at(a, a->offset + i + 1);
                    for (int64_t k = 0; k < size; k++) rows[i * size + k] = (e - s == size) ? s + k : filler;
                }
                st = carquet_arrow_take(child, &p->raw_child, rows, total, nc);
                free(rows);
            }
            if (st != CARQUET_OK) { free(nc); owned_release(o); free(o); return st; }
            o->children[0] = nc;
        }
        st = p->nkids ? apply_node(&p->kids[0], &o->children[0]) : CARQUET_OK;
        if (st != CARQUET_OK) { owned_release(o); free(o); return st; }
        release_slot(slot);
        *slot = o;
        return CARQUET_OK;
    }
    }
    return CARQUET_ERROR_INTERNAL;
}

carquet_status_t carquet_restore_apply(const carquet_restore_plan_t* p, struct ArrowArray* batch) {
    if (!p) return CARQUET_OK;
    return apply_kids(p, batch);
}
