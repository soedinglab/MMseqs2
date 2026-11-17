/**
 * @file arrow_c_stream.c
 * @brief Arrow C Stream Interface over a reader (carquet_reader_export_arrow_stream).
 *
 * Two sources, chosen once from the fields read:
 *   - the batch reader + carquet_arrow_export_batch(), when every field is a
 *     primitive or a single-level LIST / MAP of primitives (what a row batch
 *     represents): parallel decode, batches of up to batch_size rows;
 *   - carquet_reader_read_arrow_fields(), one row group per batch, otherwise.
 * Both type the arrays the same way (arrow_leaf.c), and get_schema() is the
 * tree schema restricted to the selected fields, as each batch's.
 *
 * Row selection, in order: the row groups asked for (or all), minus those a
 * predicate's statistics rule out; then the page filter (single conjunction,
 * batch path) or the row window; then the exact predicate on each batch
 * (arrow_filter.c), after which the fields read only for it are dropped.
 */

#include <carquet/carquet.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader_internal.h"

typedef struct {
    carquet_reader_t* reader;
    bool close_reader;
    int32_t* fields;              /* top-level field ordinals: output ones, then predicate-only */
    int32_t num_fields;           /* output fields */
    int32_t num_all;              /* fields read */
    struct ArrowSchema read_schema;   /* schema of the batches as read (num_all fields) */
    struct ArrowSchema out_schema;    /* schema of the batches returned */
    carquet_restore_plan_t* restore;
    carquet_compiled_predicate_t* pred;

    int32_t* rgs;                 /* row groups to visit, in order */
    int32_t nrgs;
    bool window;
    int64_t row_offset, row_end;  /* [offset, end) file rows when `window` */

    /* Batch path */
    bool batch_path;
    int32_t* leaves;
    int32_t nleaves;
    int32_t batch_size, num_threads;
    carquet_batch_reader_t* batch_reader;
    int32_t run_end;              /* rgs[next_rg .. run_end) is the current run */
    bool* in_run;                 /* [num_row_groups]: row-group filter of the run */
    const carquet_filter_clause_t* page;
    int32_t npage;

    int32_t next_rg;              /* next index into rgs */
    char message[CARQUET_ERROR_MESSAGE_MAX + 64];
} stream_t;

/* Element indices of the top-level fields (children of the root), in order.
 * Returns the count, or -1 when out of memory. */
static int32_t top_fields(const carquet_schema_t* cs, int32_t** out) {
    int32_t n = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == 0) n++;
    }
    *out = (int32_t*)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
    if (!*out) return -1;
    n = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == 0) (*out)[n++] = i;
    }
    return n;
}

static int32_t only_child(const carquet_schema_t* cs, int32_t g, int32_t* count) {
    int32_t c = -1;
    *count = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] == g) {
            if (*count == 0) c = i;
            (*count)++;
        }
    }
    return c;
}

/* Leaf ordinal of element `e`, or -1. */
static int32_t leaf_of(const carquet_schema_t* cs, int32_t e) {
    for (int32_t l = 0; l < cs->num_leaves; l++) {
        if (cs->leaf_indices[l] == e) return l;
    }
    return -1;
}

/* True when top-level element `e` is what a row batch can carry: a primitive
 * (not repeated), or a single-level LIST / MAP of primitives. Its leaf
 * ordinals are appended to leaves[*n]; *repeated is set for a LIST / MAP. */
static bool batchable(const carquet_schema_t* cs, int32_t e, int32_t* leaves, int32_t* n,
                      bool* repeated) {
    const parquet_schema_element_t* el = &cs->elements[e];
    if (el->has_type) {
        int32_t l = leaf_of(cs, e);
        if (l < 0 || cs->max_rep_levels[l] != 0) return false;
        leaves[(*n)++] = l;
        return true;
    }
    if (!el->has_logical_type) return false;
    bool is_list = el->logical_type.id == CARQUET_LOGICAL_LIST;
    bool is_map = el->logical_type.id == CARQUET_LOGICAL_MAP;
    if (!is_list && !is_map) return false;
    int32_t nrep = 0;
    int32_t rep = only_child(cs, e, &nrep);
    if (nrep != 1 || cs->elements[rep].has_type) return false;
    if (is_list && cs->elements[rep].name) {
        /* Legacy shape: a repeated group named "array" or "<list>_tuple" is the
         * element itself (a struct), which read_arrow_fields assembles. */
        const char* rn = cs->elements[rep].name;
        const char* ln = el->name;
        size_t n = ln ? strlen(ln) : 0;
        if (strcmp(rn, "array") == 0 ||
            (ln && strncmp(rn, ln, n) == 0 && strcmp(rn + n, "_tuple") == 0)) return false;
    }
    int32_t want = is_list ? 1 : 2, got = 0;
    for (int32_t i = 1; i < cs->num_elements; i++) {
        if (cs->parent_indices[i] != rep) continue;
        int32_t l = cs->elements[i].has_type ? leaf_of(cs, i) : -1;
        if (l < 0 || cs->max_rep_levels[l] != 1 || got == want) return false;
        leaves[*n + got++] = l;
    }
    if (got != want) return false;
    *n += got;
    *repeated = true;
    return true;
}

/* Top-level field ordinal holding leaf `l`. */
static int32_t field_of_leaf(const carquet_schema_t* cs, const int32_t* top, int32_t ntop, int32_t l) {
    int32_t e = cs->leaf_indices[l];
    while (e > 0 && cs->parent_indices[e] != 0) e = cs->parent_indices[e];
    for (int32_t f = 0; f < ntop; f++) if (top[f] == e) return f;
    return -1;
}

static int stream_errno(carquet_status_t st) {
    switch (st) {
    case CARQUET_ERROR_OUT_OF_MEMORY: return ENOMEM;
    case CARQUET_ERROR_INVALID_ARGUMENT:
    case CARQUET_ERROR_INVALID_STATE: return EINVAL;
    default: return EIO;
    }
}

static int stream_fail(stream_t* s, carquet_status_t st, const char* detail) {
    snprintf(s->message, sizeof(s->message), "%s: %s",
             carquet_status_string(st), detail && detail[0] ? detail : "read failed");
    return stream_errno(st);
}

static carquet_status_t stream_schema(const stream_t* s, int32_t nfields, struct ArrowSchema* out,
                                     carquet_error_t* err) {
    carquet_status_t st = carquet_arrow_export_schema(carquet_reader_schema(s->reader), out, err);
    if (st == CARQUET_OK) {
        st = carquet_arrow_schema_select_children(out, s->fields, nfields);
        if (st != CARQUET_OK) out->release(out);
    }
    return st;
}

static int stream_get_schema(struct ArrowArrayStream* self, struct ArrowSchema* out) {
    stream_t* s = (stream_t*)self->private_data;
    carquet_status_t st = carquet_arrow_schema_copy(&s->out_schema, out);
    return st == CARQUET_OK ? 0 : stream_fail(s, st, "schema copy failed");
}

static bool in_run_filter(const carquet_reader_t* r, int32_t rg, void* ctx) {
    (void)r;
    return ((const bool*)ctx)[rg];
}

/* Open the batch reader for the run starting at rgs[next_rg] (row groups in
 * ascending order), or for the whole window. */
static carquet_status_t open_run(stream_t* s, carquet_error_t* err) {
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = s->batch_size;
    cfg.num_threads = s->num_threads;
    cfg.use_mmap = carquet_reader_is_mmap(s->reader);
    cfg.column_indices = s->leaves;
    cfg.num_columns = s->nleaves;
    if (!s->window) {
        int32_t nrg = carquet_reader_num_row_groups(s->reader);
        memset(s->in_run, 0, (size_t)(nrg > 0 ? nrg : 1) * sizeof(bool));
        int32_t e = s->next_rg;
        s->in_run[s->rgs[e++]] = true;
        while (e < s->nrgs && s->rgs[e] > s->rgs[e - 1]) s->in_run[s->rgs[e++]] = true;
        s->run_end = e;
        cfg.row_group_filter = in_run_filter;
        cfg.row_group_filter_ctx = s->in_run;
    } else {
        s->run_end = s->nrgs;
    }
    s->batch_reader = carquet_batch_reader_create(s->reader, &cfg, err);
    if (!s->batch_reader) return (err && err->code != CARQUET_OK) ? err->code : CARQUET_ERROR_INTERNAL;
    carquet_status_t st = CARQUET_OK;
    if (s->window) {
        st = carquet_batch_reader_set_row_range(s->batch_reader, s->row_offset,
                                                s->row_end - s->row_offset);
    } else if (s->npage > 0) {
        st = carquet_batch_reader_set_page_filter(s->batch_reader, s->page, s->npage);
    }
    if (st != CARQUET_OK) {
        carquet_batch_reader_free(s->batch_reader);
        s->batch_reader = NULL;
        CARQUET_SET_ERROR(err, st, "Arrow stream: configuring the batch reader failed");
    }
    return st;
}

/* Next raw batch (all fields read), release == NULL at the end. */
static carquet_status_t next_raw(stream_t* s, struct ArrowArray* out, carquet_error_t* err) {
    memset(out, 0, sizeof(*out));
    if (s->batch_path) {
        for (;;) {
            if (!s->batch_reader) {
                if (s->next_rg >= s->nrgs) return CARQUET_OK;
                carquet_status_t st = open_run(s, err);
                if (st != CARQUET_OK) return st;
            }
            carquet_row_batch_t* batch = NULL;
            carquet_status_t st = carquet_batch_reader_next(s->batch_reader, &batch);
            if (st == CARQUET_ERROR_END_OF_DATA || (st == CARQUET_OK && !batch)) {
                carquet_batch_reader_free(s->batch_reader);
                s->batch_reader = NULL;
                s->next_rg = s->run_end;
                continue;
            }
            if (st != CARQUET_OK) {
                const carquet_error_t* last = carquet_batch_reader_last_error(s->batch_reader);
                if (err && last) *err = *last;
                return st;
            }
            st = carquet_arrow_export_batch(batch, carquet_reader_schema(s->reader), NULL, out, err);
            carquet_row_batch_free(batch);
            return st;
        }
    }

    while (s->next_rg < s->nrgs) {
        int32_t rg = s->rgs[s->next_rg++];
        int64_t start = 0, len = -1;
        if (s->window) {
            int64_t first = 0;
            for (int32_t g = 0; g < rg; g++) {
                carquet_row_group_metadata_t md;
                if (carquet_reader_row_group_metadata(s->reader, g, &md) == CARQUET_OK) first += md.num_rows;
            }
            carquet_row_group_metadata_t md;
            if (carquet_reader_row_group_metadata(s->reader, rg, &md) != CARQUET_OK) continue;
            int64_t lo = s->row_offset > first ? s->row_offset : first;
            int64_t hi = s->row_end < first + md.num_rows ? s->row_end : first + md.num_rows;
            if (hi <= lo) continue;
            start = lo - first;
            len = hi - lo;
        }
        carquet_status_t st = carquet_reader_read_arrow_fields(
            s->reader, rg, s->fields, s->num_all, NULL, out, err);
        if (st != CARQUET_OK) return st;
        if (len >= 0) {
            /* The window goes into the fields, not the batch: a record batch
             * with an offset of its own is refused by Arrow C++ / PyArrow
             * ("struct has non-zero offset, cannot be imported as RecordBatch"). */
            for (int64_t i = 0; i < out->n_children; i++) {
                struct ArrowArray* c = out->children[i];
                c->offset += start;
                c->length = len;
                /* Unknown for the slice; the null type has no bitmap to count from. */
                if (c->null_count != 0) c->null_count = c->n_buffers == 0 ? len : -1;
            }
            out->length = len;
        }
        return CARQUET_OK;
    }
    return CARQUET_OK;
}

/* Drop the fields read for the predicate only. */
static void drop_extra_fields(stream_t* s, struct ArrowArray* a) {
    for (int64_t i = s->num_fields; i < a->n_children; i++) {
        struct ArrowArray* c = a->children[i];
        if (c) {
            if (c->release) c->release(c);
            free(c);
            a->children[i] = NULL;
        }
    }
    if (a->n_children > s->num_fields) a->n_children = s->num_fields;
}

static int stream_get_next(struct ArrowArrayStream* self, struct ArrowArray* out) {
    stream_t* s = (stream_t*)self->private_data;
    carquet_error_t err = CARQUET_ERROR_INIT;
    for (;;) {
        struct ArrowArray raw;
        carquet_status_t st = next_raw(s, &raw, &err);
        if (st != CARQUET_OK) {
            memset(out, 0, sizeof(*out));
            return stream_fail(s, st, err.message);
        }
        if (!raw.release) {
            memset(out, 0, sizeof(*out));
            return 0;   /* end of stream */
        }
        if (s->pred && raw.length > 0) {
            uint8_t* keep = (uint8_t*)malloc((size_t)raw.length);
            int64_t kept = 0;
            st = keep ? carquet_predicate_select(s->pred, &raw, keep, &kept) : CARQUET_ERROR_OUT_OF_MEMORY;
            if (st == CARQUET_OK && kept == 0) {
                free(keep);
                raw.release(&raw);
                continue;
            }
            if (st == CARQUET_OK && kept < raw.length) {
                int64_t* rows = (int64_t*)malloc((size_t)kept * sizeof(int64_t));
                if (!rows) {
                    st = CARQUET_ERROR_OUT_OF_MEMORY;
                } else {
                    int64_t k = 0;
                    for (int64_t i = 0; i < raw.length; i++) if (keep[i]) rows[k++] = i;
                    struct ArrowArray taken;
                    st = carquet_arrow_take(&raw, &s->read_schema, rows, kept, &taken);
                    free(rows);
                    if (st == CARQUET_OK) {
                        raw.release(&raw);
                        raw = taken;
                    }
                }
            }
            free(keep);
            if (st != CARQUET_OK) {
                raw.release(&raw);
                memset(out, 0, sizeof(*out));
                return stream_fail(s, st, "filtering a batch failed");
            }
        }
        drop_extra_fields(s, &raw);
        st = carquet_restore_apply(s->restore, &raw);
        if (st != CARQUET_OK) {
            raw.release(&raw);
            memset(out, 0, sizeof(*out));
            return stream_fail(s, st, "restoring the stored Arrow types failed");
        }
        *out = raw;
        return 0;
    }
}

static const char* stream_get_last_error(struct ArrowArrayStream* self) {
    stream_t* s = (stream_t*)self->private_data;
    return s->message[0] ? s->message : NULL;
}

static void stream_free(stream_t* s) {
    if (!s) return;
    if (s->batch_reader) carquet_batch_reader_free(s->batch_reader);
    if (s->read_schema.release) s->read_schema.release(&s->read_schema);
    if (s->out_schema.release) s->out_schema.release(&s->out_schema);
    carquet_restore_free(s->restore);
    carquet_predicate_free(s->pred);
    if (s->close_reader && s->reader) carquet_reader_close(s->reader);
    free(s->fields);
    free(s->rgs);
    free(s->leaves);
    free(s->in_run);
    free(s);
}

static void stream_release(struct ArrowArrayStream* self) {
    if (!self || !self->release) return;
    stream_free((stream_t*)self->private_data);
    self->release = NULL;
    self->private_data = NULL;
}

void carquet_arrow_stream_options_init(carquet_arrow_stream_options_t* options) {
    memset(options, 0, sizeof(*options));
    options->batch_size = 65536;
    options->row_limit = -1;
}

#define FAIL(code, ...) do { CARQUET_SET_ERROR(error, code, __VA_ARGS__); st = code; goto fail; } while (0)

carquet_status_t carquet_reader_export_arrow_stream(
    carquet_reader_t* reader,
    const carquet_arrow_stream_options_t* options,
    struct ArrowArrayStream* out,
    carquet_error_t* error) {

    if (!reader || !out) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader or out");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    memset(out, 0, sizeof(*out));
    carquet_arrow_stream_options_t defaults;
    if (!options) {
        carquet_arrow_stream_options_init(&defaults);
        options = &defaults;
    }

    const carquet_schema_t* cs = carquet_reader_schema(reader);
    int32_t nrg = carquet_reader_num_row_groups(reader);
    carquet_status_t st = CARQUET_OK;
    int32_t* top = NULL;
    int32_t* top_read = NULL;
    stream_t* s = (stream_t*)calloc(1, sizeof(*s));
    int32_t ntop = top_fields(cs, &top);
    bool explicit_fields = options->field_indices && options->num_fields > 0;
    int32_t nf = explicit_fields ? options->num_fields : ntop;
    int32_t cap = (nf > 0 ? nf : 0) + (cs->num_leaves > 0 ? cs->num_leaves : 1);
    if (ntop < 0 || !s ||
        !(s->fields = (int32_t*)malloc((size_t)cap * sizeof(int32_t))) ||
        !(s->leaves = (int32_t*)malloc((size_t)(cs->num_leaves > 0 ? cs->num_leaves : 1) * sizeof(int32_t))) ||
        !(s->rgs = (int32_t*)malloc((size_t)(nrg > 0 ? nrg : 1) * sizeof(int32_t))) ||
        !(s->in_run = (bool*)calloc((size_t)(nrg > 0 ? nrg : 1), sizeof(bool)))) {
        free(top);
        stream_free(s);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Arrow stream: alloc failed");
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    s->reader = reader;
    s->num_fields = nf;
    s->batch_size = options->batch_size > 0 ? options->batch_size : 65536;
    s->num_threads = options->num_threads;

    /* Fields. */
    for (int32_t i = 0; i < nf; i++) {
        int32_t f = explicit_fields ? options->field_indices[i] : i;
        if (f < 0 || f >= ntop)
            FAIL(CARQUET_ERROR_INVALID_ARGUMENT,
                 "Arrow stream: field index %d out of range (schema has %d top-level fields)", f, ntop);
        for (int32_t j = 0; j < i; j++)
            if (s->fields[j] == f) FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: field index %d selected twice", f);
        s->fields[i] = f;
    }
    s->num_all = nf;

    /* Predicate: its fields are read too. */
    if (options->predicate) {
        const carquet_predicate_t* p = options->predicate;
        int64_t nclauses = 0;
        for (int32_t k = 0; p->conjunction_sizes && k < p->num_conjunctions; k++) nclauses += p->conjunction_sizes[k];
        for (int64_t i = 0; p->clauses && i < nclauses; i++) {
            int32_t l = p->clauses[i].column_index;
            if (l < 0 || l >= cs->num_leaves)
                FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "predicate: column %d out of range", l);
            int32_t f = field_of_leaf(cs, top, ntop, l);
            bool have = false;
            for (int32_t j = 0; j < s->num_all; j++) if (s->fields[j] == f) have = true;
            if (!have && f >= 0 && s->num_all < cap) s->fields[s->num_all++] = f;
        }
        top_read = (int32_t*)malloc((size_t)(s->num_all > 0 ? s->num_all : 1) * sizeof(int32_t));
        if (!top_read) FAIL(CARQUET_ERROR_OUT_OF_MEMORY, "Arrow stream: alloc failed");
        for (int32_t i = 0; i < s->num_all; i++) top_read[i] = top[s->fields[i]];
        st = carquet_predicate_compile(reader, p, top_read, s->num_all, &s->pred, error);
        if (st != CARQUET_OK) goto fail;
    }

    /* Source: batch reader when every field read is batchable. */
    s->batch_path = s->num_all > 0;
    bool repeated = false;
    for (int32_t i = 0; i < s->num_all && s->batch_path; i++) {
        if (!batchable(cs, top[s->fields[i]], s->leaves, &s->nleaves, &repeated)) s->batch_path = false;
    }
    if (!s->batch_path) s->nleaves = 0;

    /* Row groups. */
    if (options->row_groups && options->num_row_groups > 0) {
        if (options->row_offset != 0 || options->row_limit >= 0)
            FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: row_groups and a row window are exclusive");
        if (options->num_row_groups > nrg)
            FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: %d row groups requested, the file has %d",
                 options->num_row_groups, nrg);
        for (int32_t i = 0; i < options->num_row_groups; i++) {
            int32_t rg = options->row_groups[i];
            if (rg < 0 || rg >= nrg) FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: row group %d out of range", rg);
            if (s->in_run[rg]) FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: row group %d selected twice", rg);
            s->in_run[rg] = true;
            s->rgs[s->nrgs++] = rg;
        }
    } else {
        for (int32_t rg = 0; rg < nrg; rg++) s->rgs[s->nrgs++] = rg;
    }
    if (options->row_offset < 0) FAIL(CARQUET_ERROR_INVALID_ARGUMENT, "Arrow stream: negative row_offset");
    if (options->row_offset != 0 || options->row_limit >= 0) {
        s->window = true;
        s->row_offset = options->row_offset;
        int64_t total = carquet_reader_num_rows(reader);
        s->row_end = options->row_limit < 0 || options->row_limit > total - s->row_offset
                   ? total : s->row_offset + options->row_limit;
        if (s->row_end < s->row_offset) s->row_end = s->row_offset;
        /* The batch reader's row range does not take repeated columns. */
        if (repeated) s->batch_path = false;
    }

    /* Statistics pruning. */
    if (s->pred) {
        int32_t k = 0;
        for (int32_t i = 0; i < s->nrgs; i++) {
            if (carquet_predicate_row_group_might_match(s->pred, options->predicate, reader, s->rgs[i]))
                s->rgs[k++] = s->rgs[i];
        }
        s->nrgs = k;
        /* With a window the batch reader reads every row group: keep it exact
         * and simple by leaving pruning to the row-group path then. */
        if (s->window && s->batch_path) {
            s->nrgs = 0;
            for (int32_t rg = 0; rg < nrg; rg++) s->rgs[s->nrgs++] = rg;
        }
    }

    /* Page filter: single conjunction, batch path, offset index everywhere. */
    if (s->pred && s->batch_path && !s->window) {
        int32_t npage = 0;
        const carquet_filter_clause_t* page = carquet_predicate_page_clauses(s->pred, &npage);
        bool ok = npage > 0;
        for (int32_t i = 0; ok && i < s->nrgs; i++) {
            for (int32_t c = 0; ok && c < npage; c++) {
                carquet_column_chunk_metadata_t md;
                if (carquet_reader_column_chunk_metadata(reader, s->rgs[i], page[c].column_index, &md) != CARQUET_OK ||
                    !md.has_offset_index) ok = false;
            }
        }
        if (ok) {
            s->page = page;
            s->npage = npage;
        }
    }

    /* Schema of the batches as read, for compaction. */
    if (s->pred) {
        carquet_error_t e2 = CARQUET_ERROR_INIT;
        st = stream_schema(s, s->num_all, &s->read_schema, &e2);
        if (st != CARQUET_OK) FAIL(st, "Arrow stream: %s", e2.message);
    }

    /* Schema of the batches returned, with the stored Arrow types. */
    {
        carquet_error_t e2 = CARQUET_ERROR_INIT;
        st = stream_schema(s, s->num_fields, &s->out_schema, &e2);
        if (st != CARQUET_OK) FAIL(st, "Arrow stream: %s", e2.message);
        st = carquet_restore_plan(reader, options->restore_arrow_schema, &s->out_schema, &s->restore);
        if (st != CARQUET_OK) FAIL(st, "Arrow stream: restoring the stored Arrow schema failed");
    }

    /* Validate the batch-reader configuration now rather than at the first read. */
    if (s->batch_path && s->nrgs > 0) {
        st = open_run(s, error);
        if (st != CARQUET_OK) goto fail;
    }

    free(top);
    free(top_read);
    s->close_reader = options->close_reader;
    out->get_schema = stream_get_schema;
    out->get_next = stream_get_next;
    out->get_last_error = stream_get_last_error;
    out->release = stream_release;
    out->private_data = s;
    return CARQUET_OK;

fail:
    free(top);
    free(top_read);
    stream_free(s);   /* close_reader not set yet: the reader stays the caller's */
    return st;
}

int32_t carquet_reader_predicate_row_groups(const carquet_reader_t* reader, const carquet_predicate_t* predicate,
                                            int32_t* matching, int32_t max_matching, carquet_error_t* error) {
    if (!reader || !predicate || (max_matching > 0 && !matching)) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL reader, predicate or output");
        return -CARQUET_ERROR_INVALID_ARGUMENT;
    }
    const carquet_schema_t* cs = carquet_reader_schema(reader);
    int32_t* top = NULL;
    int32_t ntop = top_fields(cs, &top);
    if (ntop < 0) return -CARQUET_ERROR_OUT_OF_MEMORY;
    carquet_compiled_predicate_t* p = NULL;
    carquet_status_t st = carquet_predicate_compile(reader, predicate, top, ntop, &p, error);
    free(top);
    if (st != CARQUET_OK) return -(int32_t)st;
    int32_t n = 0;
    for (int32_t rg = 0; rg < carquet_reader_num_row_groups(reader); rg++) {
        if (carquet_predicate_row_group_might_match(p, predicate, reader, rg)) {
            if (n < max_matching) matching[n] = rg;
            n++;
        }
    }
    carquet_predicate_free(p);
    return n;
}
