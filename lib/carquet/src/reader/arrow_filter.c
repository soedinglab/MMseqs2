/**
 * @file arrow_filter.c
 * @brief Row predicates (carquet_predicate_t) for the Arrow stream: row-group
 *        pruning, page-filter clauses, exact evaluation on exported batches,
 *        and row compaction of Arrow arrays (carquet_arrow_take).
 *
 * A clause names a leaf column; it is evaluated on the value that column
 * exports to Arrow (arrow_leaf.c), so the literal, given in the column's
 * physical encoding, is converted once to the same representation and both
 * are compared in a common domain: int64, uint64, double, wide decimal,
 * bytes, fixed bytes or bool.
 */

#include <carquet/carquet.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "reader_internal.h"
#include "arrow_leaf.h"

/* ============================================================================
 * Compiled predicate
 * ============================================================================
 */

typedef enum { CMP_I64, CMP_U64, CMP_F64, CMP_DEC, CMP_BYTES, CMP_FIXED, CMP_BOOL } cmp_kind_t;

typedef struct {
    int64_t i;
    uint64_t u;
    double f;
    uint8_t* bytes;      /* owned: decimal (LE, width), fixed or variable bytes */
    int32_t len;
} lit_t;

typedef struct {
    carquet_filter_op_t op;
    int32_t leaf;
    cmp_kind_t kind;
    size_t width;          /* Arrow value width (fixed kinds) */
    carquet_arrow_leaf_t plan;
    int32_t depth;         /* path length: top-level field, then struct children */
    int32_t* path;         /* child index at each level; path[0] = batch child */
    lit_t* lits;           /* 1 (scalar), 2 (RANGE: lo, hi) or value_count (IN / NOT_IN) */
    int32_t nlits;
    bool has_lo, has_hi;
    bool float_like;       /* float / double / float16 column */
    bool lit_nan;          /* a literal is NaN */
} cclause_t;

struct carquet_compiled_predicate {
    cclause_t* clauses;
    int32_t nclauses;
    int32_t* conj_sizes;
    int32_t nconj;
    /* Page filter: the first conjunction's prunable clauses, values owned. */
    carquet_filter_clause_t* page;
    int32_t npage;
    void** owned;          /* buffers the page clauses point to */
    int32_t nowned, capowned;
};

static void* own(carquet_compiled_predicate_t* p, size_t n) {
    if (p->nowned == p->capowned) {
        int32_t cap = p->capowned ? 2 * p->capowned : 16;
        void** o = (void**)realloc(p->owned, (size_t)cap * sizeof(void*));
        if (!o) return NULL;
        p->owned = o;
        p->capowned = cap;
    }
    void* b = malloc(n ? n : 1);
    if (b) p->owned[p->nowned++] = b;
    return b;
}

void carquet_predicate_free(carquet_compiled_predicate_t* p) {
    if (!p) return;
    for (int32_t i = 0; i < p->nclauses; i++) {
        cclause_t* c = &p->clauses[i];
        for (int32_t k = 0; k < c->nlits; k++) free(c->lits[k].bytes);
        free(c->lits);
        free(c->path);
    }
    for (int32_t i = 0; i < p->nowned; i++) free(p->owned[i]);
    free(p->owned);
    free(p->clauses);
    free(p->conj_sizes);
    free(p->page);
    free(p);
}

static size_t physical_value_size(carquet_physical_type_t pt, int32_t type_length) {
    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN: return 1;
    case CARQUET_PHYSICAL_INT32: case CARQUET_PHYSICAL_FLOAT: return 4;
    case CARQUET_PHYSICAL_INT64: case CARQUET_PHYSICAL_DOUBLE: return 8;
    case CARQUET_PHYSICAL_INT96: return 12;
    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return type_length > 0 ? (size_t)type_length : 0;
    default: return 0;   /* BYTE_ARRAY: value_size */
    }
}

static double half_to_double(uint16_t h) {
    int sign = (h >> 15) & 1, exp = (h >> 10) & 0x1F, mant = h & 0x3FF;
    double v;
    if (exp == 0) v = ldexp((double)mant, -24);
    else if (exp == 31) v = mant ? NAN : INFINITY;
    else v = ldexp((double)(mant | 0x400), exp - 25);
    return sign ? -v : v;
}

static cmp_kind_t kind_of(const carquet_arrow_leaf_t* plan, size_t* width) {
    const char* f = plan->format;
    *width = plan->dst_width;
    if (plan->kind == CARQUET_ARROW_LEAF_BOOL) return CMP_BOOL;
    if (plan->kind == CARQUET_ARROW_LEAF_BINARY) return CMP_BYTES;
    switch (f[0]) {
    case 'd': return CMP_DEC;
    case 'f': case 'g': case 'e': return CMP_F64;
    case 'w': return CMP_FIXED;
    case 'C': case 'S': case 'I': case 'L': return CMP_U64;
    default: return CMP_I64;
    }
}

/* Literal from a physical-encoded value (`v`, `size` for BYTE_ARRAY). */
static carquet_status_t make_lit(const cclause_t* c, carquet_physical_type_t pt,
                                 const carquet_logical_type_t* lt, const void* v, int32_t size,
                                 lit_t* out) {
    memset(out, 0, sizeof(*out));
    if (!v && !(pt == CARQUET_PHYSICAL_BYTE_ARRAY && size == 0)) return CARQUET_ERROR_INVALID_ARGUMENT;
    bool unsig = lt && lt->id == CARQUET_LOGICAL_INTEGER && !lt->params.integer.is_signed;
    switch (c->kind) {
    case CMP_BOOL:
        out->i = *(const uint8_t*)v != 0;
        return CARQUET_OK;
    case CMP_I64:
    case CMP_U64:
        if (pt == CARQUET_PHYSICAL_INT96) {
            uint8_t ns[8];
            carquet_arrow_leaf_store(&c->plan, ns, v);
            memcpy(&out->i, ns, 8);
        } else if (pt == CARQUET_PHYSICAL_INT32) {
            int32_t x; memcpy(&x, v, 4);
            out->i = unsig ? (int64_t)(uint32_t)x : x;
        } else {
            int64_t x; memcpy(&x, v, 8);
            out->i = x;
        }
        out->u = (uint64_t)out->i;
        return CARQUET_OK;
    case CMP_F64:
        if (pt == CARQUET_PHYSICAL_FLOAT) { float x; memcpy(&x, v, 4); out->f = x; }
        else if (pt == CARQUET_PHYSICAL_DOUBLE) { memcpy(&out->f, v, 8); }
        else { uint16_t h; memcpy(&h, v, 2); out->f = half_to_double(h); }
        return CARQUET_OK;
    case CMP_DEC: {
        out->bytes = (uint8_t*)malloc(c->width);
        if (!out->bytes) return CARQUET_ERROR_OUT_OF_MEMORY;
        out->len = (int32_t)c->width;
        bool ok;
        if (pt == CARQUET_PHYSICAL_BYTE_ARRAY) {
            carquet_byte_array_t ba = { (uint8_t*)v, size };
            ok = carquet_arrow_leaf_store(&c->plan, out->bytes, &ba);
        } else {
            ok = carquet_arrow_leaf_store(&c->plan, out->bytes, v);
        }
        return ok ? CARQUET_OK : CARQUET_ERROR_INVALID_ARGUMENT;
    }
    case CMP_FIXED:
    case CMP_BYTES: {
        int32_t n = c->kind == CMP_FIXED ? (int32_t)c->width : size;
        if (n < 0) return CARQUET_ERROR_INVALID_ARGUMENT;
        out->bytes = (uint8_t*)malloc(n ? (size_t)n : 1);
        if (!out->bytes) return CARQUET_ERROR_OUT_OF_MEMORY;
        if (n) memcpy(out->bytes, v, (size_t)n);
        out->len = n;
        return CARQUET_OK;
    }
    }
    return CARQUET_ERROR_INVALID_ARGUMENT;
}

/* Copy `size` bytes of a page-filter value into memory the predicate owns. */
static const void* own_copy(carquet_compiled_predicate_t* p, const void* v, size_t size) {
    void* b = own(p, size);
    if (b && size) memcpy(b, v, size);
    return b;
}

carquet_status_t carquet_predicate_compile(
    const carquet_reader_t* reader, const carquet_predicate_t* pred,
    const int32_t* top_elems, int32_t ntop,
    carquet_compiled_predicate_t** out, carquet_error_t* error) {

    *out = NULL;
    const carquet_schema_t* cs = carquet_reader_schema(reader);
    if (!pred->clauses || !pred->conjunction_sizes || pred->num_conjunctions <= 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "predicate: no conjunction");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int64_t total = 0;
    for (int32_t k = 0; k < pred->num_conjunctions; k++) {
        if (pred->conjunction_sizes[k] <= 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "predicate: empty conjunction %d", k);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        total += pred->conjunction_sizes[k];
    }
    if (total > 65536) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "predicate: too many clauses");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_compiled_predicate_t* p = (carquet_compiled_predicate_t*)calloc(1, sizeof(*p));
    if (!p) return CARQUET_ERROR_OUT_OF_MEMORY;
    p->clauses = (cclause_t*)calloc((size_t)total, sizeof(cclause_t));
    p->conj_sizes = (int32_t*)malloc((size_t)pred->num_conjunctions * sizeof(int32_t));
    p->page = (carquet_filter_clause_t*)calloc((size_t)pred->conjunction_sizes[0], sizeof(carquet_filter_clause_t));
    if (!p->clauses || !p->conj_sizes || !p->page) {
        carquet_predicate_free(p);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    memcpy(p->conj_sizes, pred->conjunction_sizes, (size_t)pred->num_conjunctions * sizeof(int32_t));
    p->nconj = pred->num_conjunctions;

    carquet_status_t st = CARQUET_OK;
    for (int32_t i = 0; i < (int32_t)total && st == CARQUET_OK; i++) {
        const carquet_filter_clause_t* in = &pred->clauses[i];
        cclause_t* c = &p->clauses[i];
        p->nclauses = i + 1;
        c->op = in->op;
        c->leaf = in->column_index;
        if (in->column_index < 0 || in->column_index >= cs->num_leaves) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "predicate: column %d out of range", in->column_index);
            st = CARQUET_ERROR_INVALID_ARGUMENT;
            break;
        }
        if (cs->max_rep_levels[c->leaf] != 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "predicate: column %d is inside a LIST or MAP", c->leaf);
            st = CARQUET_ERROR_INVALID_ARGUMENT;
            break;
        }
        int32_t e = cs->leaf_indices[c->leaf];
        const parquet_schema_element_t* el = &cs->elements[e];
        const carquet_logical_type_t* lt = el->has_logical_type ? &el->logical_type : NULL;
        st = carquet_arrow_leaf_plan(el->type, el->type_length, lt, &c->plan);
        if (st != CARQUET_OK) break;
        c->kind = kind_of(&c->plan, &c->width);
        c->float_like = c->kind == CMP_F64;

        /* Path: the element's ancestors up to a top-level field. */
        int32_t chain[128], depth = 0;
        for (int32_t x = e; x > 0 && depth < 128; x = cs->parent_indices[x]) chain[depth++] = x;
        int32_t top = depth > 0 ? chain[depth - 1] : -1;
        int32_t pos = -1;
        for (int32_t t = 0; t < ntop; t++) if (top_elems[t] == top) pos = t;
        if (pos < 0 || depth >= 128) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "predicate: column %d is not among the fields read", c->leaf);
            st = CARQUET_ERROR_INVALID_ARGUMENT;
            break;
        }
        c->path = (int32_t*)malloc((size_t)depth * sizeof(int32_t));
        if (!c->path) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
        c->depth = depth;
        c->path[0] = pos;
        for (int32_t d = 1; d < depth; d++) {
            int32_t parent = chain[depth - d], child = chain[depth - d - 1], k = 0;
            for (int32_t y = 1; y < child; y++) if (cs->parent_indices[y] == parent) k++;
            c->path[d] = k;
        }

        /* Literals. */
        size_t vsize = physical_value_size(el->type, el->type_length);
        bool is_ba = el->type == CARQUET_PHYSICAL_BYTE_ARRAY;
        /* A fixed-width literal is read at the column's width: a size given
         * for it must say the same (a short buffer would be over-read). */
        bool sized_wrong = !is_ba && vsize > 0 &&
            ((in->value && in->value_size > 0 && (size_t)in->value_size != vsize) ||
             (in->has_lo && in->lo_size > 0 && (size_t)in->lo_size != vsize) ||
             (in->has_hi && in->hi_size > 0 && (size_t)in->hi_size != vsize));
        if (sized_wrong) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "predicate: clause %d gives a %d-byte value for a %d-byte column", i,
                (int)in->value_size, (int)vsize);
            st = CARQUET_ERROR_INVALID_ARGUMENT;
            break;
        }
        switch (c->op) {
        case CARQUET_FILTER_EQ: case CARQUET_FILTER_NE: case CARQUET_FILTER_LT:
        case CARQUET_FILTER_LE: case CARQUET_FILTER_GT: case CARQUET_FILTER_GE:
            if (el->type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && in->value_size != el->type_length) {
                st = CARQUET_ERROR_INVALID_ARGUMENT;
                break;
            }
            c->lits = (lit_t*)calloc(1, sizeof(lit_t));
            if (!c->lits) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            c->nlits = 1;
            st = make_lit(c, el->type, lt, in->value, in->value_size, &c->lits[0]);
            break;
        case CARQUET_FILTER_RANGE:
            c->has_lo = in->has_lo;
            c->has_hi = in->has_hi;
            if (!c->has_lo && !c->has_hi) { st = CARQUET_ERROR_INVALID_ARGUMENT; break; }
            /* As for the unary ops (and the page filter, which would refuse
             * the clause later, on files with a page index only). */
            if (el->type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY &&
                ((c->has_lo && in->lo_size != el->type_length) ||
                 (c->has_hi && in->hi_size != el->type_length))) {
                st = CARQUET_ERROR_INVALID_ARGUMENT;
                break;
            }
            c->lits = (lit_t*)calloc(2, sizeof(lit_t));
            if (!c->lits) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            c->nlits = 2;
            if (c->has_lo) st = make_lit(c, el->type, lt, in->lo, in->lo_size, &c->lits[0]);
            if (st == CARQUET_OK && c->has_hi) st = make_lit(c, el->type, lt, in->hi, in->hi_size, &c->lits[1]);
            break;
        case CARQUET_FILTER_IN:
        case CARQUET_FILTER_NOT_IN:
            if (in->value_count <= 0 || !in->values) { st = CARQUET_ERROR_INVALID_ARGUMENT; break; }
            c->lits = (lit_t*)calloc((size_t)in->value_count, sizeof(lit_t));
            if (!c->lits) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            c->nlits = in->value_count;
            for (int32_t k = 0; k < in->value_count && st == CARQUET_OK; k++) {
                if (is_ba) {
                    const carquet_byte_array_t* b = &((const carquet_byte_array_t*)in->values)[k];
                    st = make_lit(c, el->type, lt, b->data, b->length, &c->lits[k]);
                } else {
                    st = make_lit(c, el->type, lt, (const uint8_t*)in->values + (size_t)k * vsize,
                                  (int32_t)vsize, &c->lits[k]);
                }
            }
            break;
        case CARQUET_FILTER_IS_NULL:
        case CARQUET_FILTER_IS_NOT_NULL:
            break;
        default:
            st = CARQUET_ERROR_INVALID_ARGUMENT;
        }
        if (st != CARQUET_OK) {
            CARQUET_SET_ERROR(error, st, "predicate: clause %d on column %d is malformed", i, c->leaf);
            break;
        }
        for (int32_t k = 0; k < c->nlits; k++) {
            if (c->kind == CMP_F64 && isnan(c->lits[k].f)) c->lit_nan = true;
        }

        /* Page filter copy (first conjunction only, prunable ops). */
        if (i < pred->conjunction_sizes[0] && el->type != CARQUET_PHYSICAL_INT96 && !c->lit_nan &&
            c->op != CARQUET_FILTER_NOT_IN && !(c->op == CARQUET_FILTER_NE && c->float_like) &&
            !((c->op == CARQUET_FILTER_IN) && in->value_count > 256)) {
            carquet_filter_clause_t* pc = &p->page[p->npage];
            /* Only the members the op reads are copied: the others may hold
             * anything in a clause the caller did not zero. */
            memset(pc, 0, sizeof(*pc));
            pc->column_index = in->column_index;
            pc->op = in->op;
            size_t sz = is_ba ? (size_t)(in->value_size > 0 ? in->value_size : 0) : vsize;
            bool ok = true;
            bool is_range = in->op == CARQUET_FILTER_RANGE;
            bool is_set = in->op == CARQUET_FILTER_IN;
            bool is_unary = !is_range && !is_set && in->op != CARQUET_FILTER_IS_NULL &&
                            in->op != CARQUET_FILTER_IS_NOT_NULL;
            if (is_unary) {
                pc->value_size = in->value_size;
                ok = ok && (pc->value = own_copy(p, in->value, sz)) != NULL;
            }
            if (is_range && in->has_lo) {
                pc->has_lo = true;
                pc->lo_size = in->lo_size;
                ok = ok && (pc->lo = own_copy(p, in->lo, is_ba ? (size_t)(in->lo_size > 0 ? in->lo_size : 0) : vsize)) != NULL;
            }
            if (is_range && in->has_hi) {
                pc->has_hi = true;
                pc->hi_size = in->hi_size;
                ok = ok && (pc->hi = own_copy(p, in->hi, is_ba ? (size_t)(in->hi_size > 0 ? in->hi_size : 0) : vsize)) != NULL;
            }
            if (is_set) {
                pc->value_count = in->value_count;
                if (is_ba) {
                    carquet_byte_array_t* arr = (carquet_byte_array_t*)own(p, (size_t)in->value_count * sizeof(*arr));
                    ok = ok && arr;
                    for (int32_t k = 0; ok && k < in->value_count; k++) {
                        const carquet_byte_array_t* b = &((const carquet_byte_array_t*)in->values)[k];
                        arr[k].length = b->length;
                        arr[k].data = (uint8_t*)own_copy(p, b->data, (size_t)(b->length > 0 ? b->length : 0));
                        ok = arr[k].data != NULL;
                    }
                    pc->values = arr;
                } else {
                    ok = ok && (pc->values = own_copy(p, in->values, (size_t)in->value_count * vsize)) != NULL;
                }
            }
            if (!ok) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            p->npage++;
        }
    }
    if (st != CARQUET_OK) {
        carquet_predicate_free(p);
        return st;
    }
    *out = p;
    return CARQUET_OK;
}

/* ============================================================================
 * Row-group pruning (statistics)
 * ============================================================================
 */

static bool clause_rg_might_match(const carquet_reader_t* r, int32_t rg, const cclause_t* c,
                                  const carquet_filter_clause_t* in) {
    const carquet_schema_t* cs = carquet_reader_schema(r);
    const parquet_schema_element_t* el = &cs->elements[cs->leaf_indices[c->leaf]];
    if (c->op == CARQUET_FILTER_IS_NULL || c->op == CARQUET_FILTER_IS_NOT_NULL) {
        carquet_column_statistics_t s;
        if (carquet_reader_column_statistics(r, rg, c->leaf, &s) != CARQUET_OK || !s.has_null_count) return true;
        return c->op == CARQUET_FILTER_IS_NULL ? s.null_count > 0
                                                : (s.num_values <= 0 || s.null_count < s.num_values);
    }
    if (el->type == CARQUET_PHYSICAL_INT96 || c->lit_nan) return true;
    bool m = true;
    switch (c->op) {
    case CARQUET_FILTER_EQ: case CARQUET_FILTER_LT: case CARQUET_FILTER_LE:
    case CARQUET_FILTER_GT: case CARQUET_FILTER_GE:
        if (carquet_reader_row_group_matches(r, rg, c->leaf, (carquet_compare_op_t)c->op,
                                             in->value, in->value_size, &m) != CARQUET_OK) return true;
        return m;
    case CARQUET_FILTER_NE:
        if (c->float_like) return true;   /* statistics skip NaN, which NE keeps */
        if (carquet_reader_row_group_matches(r, rg, c->leaf, CARQUET_COMPARE_NE,
                                             in->value, in->value_size, &m) != CARQUET_OK) return true;
        return m;
    case CARQUET_FILTER_RANGE:
        if (in->has_lo && carquet_reader_row_group_matches(r, rg, c->leaf, CARQUET_COMPARE_GE,
                in->lo, in->lo_size, &m) == CARQUET_OK && !m) return false;
        if (in->has_hi && carquet_reader_row_group_matches(r, rg, c->leaf, CARQUET_COMPARE_LE,
                in->hi, in->hi_size, &m) == CARQUET_OK && !m) return false;
        return true;
    case CARQUET_FILTER_IN: {
        size_t vsize = physical_value_size(el->type, el->type_length);
        for (int32_t k = 0; k < in->value_count; k++) {
            const void* v;
            int32_t sz;
            if (el->type == CARQUET_PHYSICAL_BYTE_ARRAY) {
                const carquet_byte_array_t* b = &((const carquet_byte_array_t*)in->values)[k];
                v = b->data; sz = b->length;
                if (!v) return true;
            } else {
                v = (const uint8_t*)in->values + (size_t)k * vsize; sz = (int32_t)vsize;
            }
            if (carquet_reader_row_group_matches(r, rg, c->leaf, CARQUET_COMPARE_EQ, v, sz, &m) != CARQUET_OK || m)
                return true;
        }
        return false;
    }
    default:
        return true;
    }
}

bool carquet_predicate_row_group_might_match(const carquet_compiled_predicate_t* p,
                                             const carquet_predicate_t* pred,
                                             const carquet_reader_t* r, int32_t rg) {
    int32_t base = 0;
    for (int32_t k = 0; k < p->nconj; k++) {
        bool all = true;
        for (int32_t j = 0; j < p->conj_sizes[k] && all; j++) {
            if (!clause_rg_might_match(r, rg, &p->clauses[base + j], &pred->clauses[base + j])) all = false;
        }
        if (all) return true;
        base += p->conj_sizes[k];
    }
    return false;
}

const carquet_filter_clause_t* carquet_predicate_page_clauses(const carquet_compiled_predicate_t* p,
                                                              int32_t* count) {
    *count = p->nconj == 1 ? p->npage : 0;
    return p->page;
}

void carquet_predicate_columns(const carquet_compiled_predicate_t* p, int32_t* leaves, int32_t* n) {
    *n = 0;
    for (int32_t i = 0; i < p->nclauses; i++) {
        bool seen = false;
        for (int32_t j = 0; j < *n; j++) if (leaves[j] == p->clauses[i].leaf) seen = true;
        if (!seen) leaves[(*n)++] = p->clauses[i].leaf;
    }
}

/* ============================================================================
 * Exact evaluation
 * ============================================================================
 */

static bool bit_at(const void* bits, int64_t i) {
    return (((const uint8_t*)bits)[i >> 3] >> (i & 7)) & 1;
}

static int cmp_dec(const uint8_t* a, const uint8_t* b, size_t w) {
    /* Little-endian two's complement: most significant byte signed. */
    int8_t ha = (int8_t)a[w - 1], hb = (int8_t)b[w - 1];
    if (ha != hb) return ha < hb ? -1 : 1;
    for (size_t i = w - 1; i-- > 0;) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    return 0;
}

static int cmp_bytes(const uint8_t* a, int32_t na, const uint8_t* b, int32_t nb) {
    int32_t n = na < nb ? na : nb;
    int c = n ? memcmp(a, b, (size_t)n) : 0;
    if (c) return c < 0 ? -1 : 1;
    return na == nb ? 0 : (na < nb ? -1 : 1);
}

/* Three-way comparison of the value at physical index `j` of `a` with `l`;
 * 2 when unordered (a NaN is involved). */
static int compare_at(const cclause_t* c, const struct ArrowArray* a, int64_t j, const lit_t* l) {
    const uint8_t* d = (const uint8_t*)a->buffers[1];
    switch (c->kind) {
    case CMP_BOOL: {
        int64_t v = bit_at(d, j);
        return v == l->i ? 0 : (v < l->i ? -1 : 1);
    }
    case CMP_I64: {
        int64_t v;
        switch (c->width) {
        case 1: v = ((const int8_t*)d)[j]; break;
        case 2: { int16_t x; memcpy(&x, d + 2 * j, 2); v = x; break; }
        case 4: { int32_t x; memcpy(&x, d + 4 * j, 4); v = x; break; }
        default: memcpy(&v, d + 8 * j, 8); break;
        }
        return v == l->i ? 0 : (v < l->i ? -1 : 1);
    }
    case CMP_U64: {
        uint64_t v;
        switch (c->width) {
        case 1: v = d[j]; break;
        case 2: { uint16_t x; memcpy(&x, d + 2 * j, 2); v = x; break; }
        case 4: { uint32_t x; memcpy(&x, d + 4 * j, 4); v = x; break; }
        default: memcpy(&v, d + 8 * j, 8); break;
        }
        return v == l->u ? 0 : (v < l->u ? -1 : 1);
    }
    case CMP_F64: {
        double v;
        if (c->width == 4) { float x; memcpy(&x, d + 4 * j, 4); v = x; }
        else if (c->width == 8) memcpy(&v, d + 8 * j, 8);
        else { uint16_t h; memcpy(&h, d + 2 * j, 2); v = half_to_double(h); }
        if (isnan(v) || isnan(l->f)) return (isnan(v) && isnan(l->f)) ? 3 : 2;   /* 3: both NaN */
        return v == l->f ? 0 : (v < l->f ? -1 : 1);
    }
    case CMP_DEC:
        return cmp_dec(d + c->width * (size_t)j, l->bytes, c->width);
    case CMP_FIXED:
        return cmp_bytes(d + c->width * (size_t)j, (int32_t)c->width, l->bytes, l->len);
    case CMP_BYTES: {
        const int32_t* off = (const int32_t*)a->buffers[1];
        const uint8_t* data = (const uint8_t*)a->buffers[2];
        return cmp_bytes(data + off[j], off[j + 1] - off[j], l->bytes, l->len);
    }
    }
    return 2;
}

static bool test_value(const cclause_t* c, const struct ArrowArray* a, int64_t j) {
    int r;
    switch (c->op) {
    case CARQUET_FILTER_EQ: return compare_at(c, a, j, &c->lits[0]) == 0;
    case CARQUET_FILTER_NE: r = compare_at(c, a, j, &c->lits[0]); return r != 0;
    case CARQUET_FILTER_LT: r = compare_at(c, a, j, &c->lits[0]); return r == -1;
    case CARQUET_FILTER_LE: r = compare_at(c, a, j, &c->lits[0]); return r == -1 || r == 0;
    case CARQUET_FILTER_GT: r = compare_at(c, a, j, &c->lits[0]); return r == 1;
    case CARQUET_FILTER_GE: r = compare_at(c, a, j, &c->lits[0]); return r == 1 || r == 0;
    case CARQUET_FILTER_RANGE:
        if (c->has_lo) { r = compare_at(c, a, j, &c->lits[0]); if (r != 0 && r != 1) return false; }
        if (c->has_hi) { r = compare_at(c, a, j, &c->lits[1]); if (r != 0 && r != -1) return false; }
        return true;
    case CARQUET_FILTER_IN:
    case CARQUET_FILTER_NOT_IN: {
        bool found = false;
        for (int32_t k = 0; k < c->nlits && !found; k++) {
            r = compare_at(c, a, j, &c->lits[k]);
            found = r == 0 || r == 3;   /* NaN is a member of a set holding NaN */
        }
        return c->op == CARQUET_FILTER_IN ? found : !found;
    }
    default:
        return false;
    }
}

/* Mark in `match[i]` (rows 0..n-1 of the batch) whether clause `c` holds.
 * Row i of a struct with offset o is row i + o of each child, so the index
 * of row i in a node's buffers is i + shift + node->offset, where `shift`
 * sums the offsets of the node's ancestors (the batch's included). */
static carquet_status_t eval_clause(const cclause_t* c, const struct ArrowArray* batch, uint8_t* match) {
    int64_t n = batch->length;
    const struct ArrowArray* node = batch;
    int64_t shift = batch->offset;
    memset(match, 1, (size_t)n);      /* 1 = not null so far */
    for (int32_t d = 0; d < c->depth; d++) {
        if (c->path[d] >= node->n_children || !node->children[c->path[d]]) return CARQUET_ERROR_INTERNAL;
        node = node->children[c->path[d]];
        if (node->offset < 0 || node->length < shift + n) return CARQUET_ERROR_INTERNAL;
        const void* valid = node->n_buffers > 0 ? node->buffers[0] : NULL;
        if (node->n_buffers == 0) {
            memset(match, 0, (size_t)n);   /* the null type: every row is null */
        } else if (valid && node->null_count != 0) {
            for (int64_t i = 0; i < n; i++) {
                if (match[i] && !bit_at(valid, i + shift + node->offset)) match[i] = 0;
            }
        }
        shift += node->offset;
    }
    if (c->op == CARQUET_FILTER_IS_NULL || c->op == CARQUET_FILTER_IS_NOT_NULL) {
        bool want_null = c->op == CARQUET_FILTER_IS_NULL;
        for (int64_t i = 0; i < n; i++) match[i] = (match[i] == 0) == want_null;
        return CARQUET_OK;
    }
    if (n == 0 || node->n_buffers == 0) {   /* value ops on the null type: false */
        memset(match, 0, (size_t)n);
        return CARQUET_OK;
    }
    if (node->n_buffers < 2 || !node->buffers[1] ||
        (c->kind == CMP_BYTES && (node->n_buffers < 3 || !node->buffers[2]))) {
        return CARQUET_ERROR_INTERNAL;
    }
    for (int64_t i = 0; i < n; i++) {
        if (match[i]) match[i] = test_value(c, node, i + shift);
    }
    return CARQUET_OK;
}

carquet_status_t carquet_predicate_select(const carquet_compiled_predicate_t* p,
                                          const struct ArrowArray* batch,
                                          uint8_t* keep, int64_t* kept) {
    int64_t n = batch->length;
    uint8_t* conj = (uint8_t*)malloc((size_t)(n > 0 ? n : 1));
    uint8_t* m = (uint8_t*)malloc((size_t)(n > 0 ? n : 1));
    if (!conj || !m) { free(conj); free(m); return CARQUET_ERROR_OUT_OF_MEMORY; }
    memset(keep, 0, (size_t)n);
    int32_t base = 0;
    carquet_status_t st = CARQUET_OK;
    for (int32_t k = 0; k < p->nconj && st == CARQUET_OK; k++) {
        memset(conj, 1, (size_t)n);
        for (int32_t j = 0; j < p->conj_sizes[k] && st == CARQUET_OK; j++) {
            st = eval_clause(&p->clauses[base + j], batch, m);
            for (int64_t i = 0; i < n; i++) conj[i] &= m[i];
        }
        for (int64_t i = 0; i < n; i++) keep[i] |= conj[i];
        base += p->conj_sizes[k];
    }
    int64_t c = 0;
    for (int64_t i = 0; i < n; i++) c += keep[i];
    *kept = c;
    free(conj);
    free(m);
    return st;
}

/* ============================================================================
 * Compaction: keep rows of an Arrow array
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
    a->release = NULL;
}

/* Bytes per value of a fixed-width Arrow format, 0 for anything else. */
static size_t fixed_width(const char* f) {
    switch (f[0]) {
    case 'c': case 'C': return 1;
    case 's': case 'S': case 'e': return 2;
    case 'i': case 'I': case 'f': return 4;
    case 'l': case 'L': case 'g': return 8;
    case 'w': return f[1] == ':' ? (size_t)atoi(f + 2) : 0;
    case 'd': {
        const char* last = strrchr(f, ',');
        int commas = 0;
        for (const char* q = f; *q; q++) commas += *q == ',';
        if (commas >= 2 && last) {
            int bits = atoi(last + 1);
            return bits == 256 ? 32 : bits == 64 ? 8 : bits == 32 ? 4 : 16;
        }
        return 16;
    }
    case 't':
        if (f[1] == 'd') return f[2] == 'D' ? 4 : 8;
        if (f[1] == 't') return (f[2] == 's' || f[2] == 'm') ? 4 : 8;
        return 8;   /* timestamps, durations */
    default: return 0;
    }
}

static carquet_status_t take_node(const struct ArrowArray* in, const struct ArrowSchema* s,
                                  const int64_t* pos, int64_t n, struct ArrowArray* out, int depth);

static uint8_t* take_validity(const struct ArrowArray* in, const int64_t* pos, int64_t n,
                              int64_t* nulls, bool* oom) {
    *nulls = 0;
    *oom = false;
    const void* v = in->n_buffers > 0 ? in->buffers[0] : NULL;
    if (!v || in->null_count == 0) return NULL;
    uint8_t* o = (uint8_t*)calloc((size_t)((n + 7) / 8 + 1), 1);
    if (!o) { *oom = true; return NULL; }
    for (int64_t k = 0; k < n; k++) {
        if (bit_at(v, pos[k] + in->offset)) o[k >> 3] |= (uint8_t)(1u << (k & 7));
        else (*nulls)++;
    }
    if (*nulls == 0) { free(o); return NULL; }
    return o;
}

static carquet_status_t take_node(const struct ArrowArray* in, const struct ArrowSchema* s,
                                  const int64_t* pos, int64_t n, struct ArrowArray* out, int depth) {
    memset(out, 0, sizeof(*out));
    if (depth > 64 || !s->format) return CARQUET_ERROR_INVALID_ARGUMENT;
    const char* f = s->format;
    out->release = owned_release;
    out->length = n;
    bool oom;
    int64_t nulls;
    uint8_t* valid = take_validity(in, pos, n, &nulls, &oom);
    if (oom) return CARQUET_ERROR_OUT_OF_MEMORY;
    out->null_count = nulls;

    if (f[0] == '+' && (f[1] == 's')) {
        out->n_buffers = 1;
        out->buffers = (const void**)calloc(1, sizeof(void*));
        if (!out->buffers) { free(valid); out->release = NULL; return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->buffers[0] = valid;
        if (in->n_children != s->n_children) return CARQUET_ERROR_INVALID_ARGUMENT;
        int64_t* cpos = (int64_t*)malloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));
        out->children = (struct ArrowArray**)calloc((size_t)(in->n_children > 0 ? in->n_children : 1), sizeof(void*));
        if (!cpos || !out->children) { free(cpos); return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->n_children = in->n_children;
        for (int64_t k = 0; k < n; k++) cpos[k] = pos[k] + in->offset;
        for (int64_t c = 0; c < in->n_children; c++) {
            out->children[c] = (struct ArrowArray*)malloc(sizeof(struct ArrowArray));
            if (!out->children[c]) { free(cpos); return CARQUET_ERROR_OUT_OF_MEMORY; }
            carquet_status_t st = take_node(in->children[c], s->children[c], cpos, n, out->children[c], depth + 1);
            if (st != CARQUET_OK) { free(cpos); return st; }
        }
        free(cpos);
        return CARQUET_OK;
    }

    if (f[0] == '+' && (f[1] == 'l' || f[1] == 'L' || f[1] == 'm')) {
        bool large = f[1] == 'L';
        out->n_buffers = 2;
        out->buffers = (const void**)calloc(2, sizeof(void*));
        if (!out->buffers) { free(valid); out->release = NULL; return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->buffers[0] = valid;
        size_t ow = large ? 8 : 4;
        uint8_t* offs = (uint8_t*)malloc((size_t)(n + 1) * ow);
        out->buffers[1] = offs;
        if (!offs) return CARQUET_ERROR_OUT_OF_MEMORY;
        const uint8_t* io = (const uint8_t*)in->buffers[1];
        int64_t total = 0;
        for (int64_t k = 0; k < n; k++) {
            int64_t j = pos[k] + in->offset, a, b;
            if (large) { memcpy(&a, io + 8 * j, 8); memcpy(&b, io + 8 * (j + 1), 8); }
            else { int32_t x, y; memcpy(&x, io + 4 * j, 4); memcpy(&y, io + 4 * (j + 1), 4); a = x; b = y; }
            if (large) memcpy(offs + 8 * k, &total, 8);
            else { int32_t t = (int32_t)total; memcpy(offs + 4 * k, &t, 4); }
            total += b - a;
        }
        if (large) memcpy(offs + 8 * n, &total, 8);
        else { int32_t t = (int32_t)total; memcpy(offs + 4 * n, &t, 4); }
        int64_t* cpos = (int64_t*)malloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
        out->children = (struct ArrowArray**)calloc(1, sizeof(void*));
        if (!cpos || !out->children) { free(cpos); return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->n_children = 1;
        int64_t q = 0;
        for (int64_t k = 0; k < n; k++) {
            int64_t j = pos[k] + in->offset, a, b;
            if (large) { memcpy(&a, io + 8 * j, 8); memcpy(&b, io + 8 * (j + 1), 8); }
            else { int32_t x, y; memcpy(&x, io + 4 * j, 4); memcpy(&y, io + 4 * (j + 1), 4); a = x; b = y; }
            for (int64_t x = a; x < b; x++) cpos[q++] = x;
        }
        out->children[0] = (struct ArrowArray*)malloc(sizeof(struct ArrowArray));
        carquet_status_t st = out->children[0]
            ? take_node(in->children[0], s->children[0], cpos, total, out->children[0], depth + 1)
            : CARQUET_ERROR_OUT_OF_MEMORY;
        free(cpos);
        return st;
    }

    if (f[0] == 'n' && f[1] == '\0') {
        out->null_count = n;
        out->n_buffers = 0;
        free(valid);
        return CARQUET_OK;
    }

    if ((f[0] == 'u' || f[0] == 'z' || f[0] == 'U' || f[0] == 'Z') && f[1] == '\0') {
        bool large = f[0] == 'U' || f[0] == 'Z';
        size_t ow = large ? 8 : 4;
        out->n_buffers = 3;
        out->buffers = (const void**)calloc(3, sizeof(void*));
        if (!out->buffers) { free(valid); out->release = NULL; return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->buffers[0] = valid;
        const uint8_t* io = (const uint8_t*)in->buffers[1];
        const uint8_t* id = (const uint8_t*)in->buffers[2];
        int64_t total = 0;
        for (int64_t k = 0; k < n; k++) {
            int64_t j = pos[k] + in->offset, a, b;
            if (large) { memcpy(&a, io + 8 * j, 8); memcpy(&b, io + 8 * (j + 1), 8); }
            else { int32_t x, y; memcpy(&x, io + 4 * j, 4); memcpy(&y, io + 4 * (j + 1), 4); a = x; b = y; }
            total += b - a;
        }
        uint8_t* offs = (uint8_t*)malloc((size_t)(n + 1) * ow);
        uint8_t* data = (uint8_t*)malloc((size_t)(total > 0 ? total : 1));
        out->buffers[1] = offs;
        out->buffers[2] = data;
        if (!offs || !data) return CARQUET_ERROR_OUT_OF_MEMORY;
        int64_t at = 0;
        for (int64_t k = 0; k < n; k++) {
            int64_t j = pos[k] + in->offset, a, b;
            if (large) { memcpy(&a, io + 8 * j, 8); memcpy(&b, io + 8 * (j + 1), 8); }
            else { int32_t x, y; memcpy(&x, io + 4 * j, 4); memcpy(&y, io + 4 * (j + 1), 4); a = x; b = y; }
            if (large) memcpy(offs + 8 * k, &at, 8);
            else { int32_t t = (int32_t)at; memcpy(offs + 4 * k, &t, 4); }
            if (b > a) memcpy(data + at, id + a, (size_t)(b - a));
            at += b - a;
        }
        if (large) memcpy(offs + 8 * n, &at, 8);
        else { int32_t t = (int32_t)at; memcpy(offs + 4 * n, &t, 4); }
        return CARQUET_OK;
    }

    out->n_buffers = 2;
    out->buffers = (const void**)calloc(2, sizeof(void*));
    if (!out->buffers) { free(valid); out->release = NULL; return CARQUET_ERROR_OUT_OF_MEMORY; }
    out->buffers[0] = valid;
    if (f[0] == 'b' && f[1] == '\0') {
        uint8_t* bits = (uint8_t*)calloc((size_t)((n + 7) / 8 + 1), 1);
        out->buffers[1] = bits;
        if (!bits) return CARQUET_ERROR_OUT_OF_MEMORY;
        for (int64_t k = 0; k < n; k++) {
            if (bit_at(in->buffers[1], pos[k] + in->offset)) bits[k >> 3] |= (uint8_t)(1u << (k & 7));
        }
        return CARQUET_OK;
    }
    size_t w = fixed_width(f);
    if (w == 0) return CARQUET_ERROR_NOT_IMPLEMENTED;
    uint8_t* data = (uint8_t*)malloc((size_t)(n > 0 ? n : 1) * w);
    out->buffers[1] = data;
    if (!data) return CARQUET_ERROR_OUT_OF_MEMORY;
    const uint8_t* id = (const uint8_t*)in->buffers[1];
    for (int64_t k = 0; k < n; k++) memcpy(data + (size_t)k * w, id + (size_t)(pos[k] + in->offset) * w, w);
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_take(const struct ArrowArray* in, const struct ArrowSchema* schema,
                                    const int64_t* rows, int64_t n, struct ArrowArray* out) {
    carquet_status_t st = take_node(in, schema, rows, n, out, 0);
    if (st != CARQUET_OK && out->release) out->release(out);
    return st;
}

/* ============================================================================
 * All-null arrays
 * ============================================================================
 */

static carquet_status_t null_node(const struct ArrowSchema* s, int64_t len, struct ArrowArray* out, int depth) {
    memset(out, 0, sizeof(*out));
    if (depth > 64 || !s->format) return CARQUET_ERROR_INVALID_ARGUMENT;
    const char* f = s->format;
    out->release = owned_release;
    out->length = len;
    out->null_count = len;
    size_t vbytes = (size_t)((len + 7) / 8 + 1);
    if (f[0] == 'n' && f[1] == '\0') return CARQUET_OK;   /* no buffers */
    int64_t nbuf;
    size_t data_bytes = 0, off_bytes = 0;
    int nkids = 0;
    int64_t kid_len = 0;
    if (f[0] == '+' && f[1] == 's') {
        nbuf = 1; nkids = (int)s->n_children; kid_len = len;
    } else if (f[0] == '+' && (f[1] == 'l' || f[1] == 'm')) {
        nbuf = 2; off_bytes = (size_t)(len + 1) * 4; nkids = 1;
    } else if (f[0] == '+' && f[1] == 'L') {
        nbuf = 2; off_bytes = (size_t)(len + 1) * 8; nkids = 1;
    } else if (f[0] == '+' && f[1] == 'w') {
        nbuf = 1; nkids = 1; kid_len = len * (int64_t)atoi(f + 3);
    } else if ((f[0] == 'u' || f[0] == 'z') && f[1] == '\0') {
        nbuf = 3; off_bytes = (size_t)(len + 1) * 4; data_bytes = 1;
    } else if ((f[0] == 'U' || f[0] == 'Z') && f[1] == '\0') {
        nbuf = 3; off_bytes = (size_t)(len + 1) * 8; data_bytes = 1;
    } else if (f[0] == 'b' && f[1] == '\0') {
        nbuf = 2; data_bytes = vbytes;
    } else {
        size_t w = fixed_width(f);
        if (w == 0) return CARQUET_ERROR_NOT_IMPLEMENTED;
        nbuf = 2; data_bytes = (size_t)(len > 0 ? len : 1) * w;
    }
    out->n_buffers = nbuf;
    out->buffers = (const void**)calloc((size_t)nbuf, sizeof(void*));
    if (!out->buffers) { out->release = NULL; return CARQUET_ERROR_OUT_OF_MEMORY; }
    out->buffers[0] = calloc(vbytes, 1);
    if (!out->buffers[0]) return CARQUET_ERROR_OUT_OF_MEMORY;
    if (off_bytes) {
        out->buffers[1] = calloc(off_bytes, 1);
        if (!out->buffers[1]) return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (data_bytes) {
        int64_t slot = off_bytes ? 2 : 1;
        out->buffers[slot] = calloc(data_bytes, 1);
        if (!out->buffers[slot]) return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (nkids > 0) {
        if (s->n_children < nkids) return CARQUET_ERROR_INVALID_ARGUMENT;
        out->children = (struct ArrowArray**)calloc((size_t)nkids, sizeof(void*));
        if (!out->children) return CARQUET_ERROR_OUT_OF_MEMORY;
        out->n_children = nkids;
        for (int i = 0; i < nkids; i++) {
            out->children[i] = (struct ArrowArray*)malloc(sizeof(struct ArrowArray));
            if (!out->children[i]) return CARQUET_ERROR_OUT_OF_MEMORY;
            carquet_status_t st = null_node(s->children[i], kid_len, out->children[i], depth + 1);
            if (st != CARQUET_OK) return st;
        }
    }
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_null_array(const struct ArrowSchema* schema, int64_t len, struct ArrowArray* out) {
    carquet_status_t st = null_node(schema, len, out, 0);
    if (st != CARQUET_OK && out->release) out->release(out);
    return st;
}
