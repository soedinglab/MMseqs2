/**
 * @file metricstore.c
 * @brief MetricStore implementation — a realistic carquet consumer.
 *
 * Column layout (leaf indices are stable and referenced throughout):
 *
 *   0  event_id    INT64                    REQUIRED   DELTA_BINARY_PACKED, bloom
 *   1  ts          INT64  TIMESTAMP(us,utc) REQUIRED
 *   2  host        BYTE_ARRAY STRING        REQUIRED   RLE_DICTIONARY, bloom
 *   3  region      BYTE_ARRAY STRING        REQUIRED   RLE_DICTIONARY
 *   4  metric      BYTE_ARRAY STRING        REQUIRED   RLE_DICTIONARY, bloom
 *   5  value       DOUBLE                   REQUIRED   BYTE_STREAM_SPLIT
 *   6  error_code  INT32  INTEGER(32,signed) OPTIONAL
 *   7  session_id  FIXED_LEN_BYTE_ARRAY(16) UUID  REQUIRED
 */
#include "ms_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Small helpers (shared, declared in ms_internal.h)                          */
/* ------------------------------------------------------------------------- */

void ms_set_err(char* buf, size_t len, const char* msg) {
    if (buf && len) {
        strncpy(buf, msg, len - 1);
        buf[len - 1] = '\0';
    }
}

carquet_compression_t ms_to_carquet_codec(ms_codec_t c) {
    switch (c) {
        case MS_CODEC_NONE:   return CARQUET_COMPRESSION_UNCOMPRESSED;
        case MS_CODEC_SNAPPY: return CARQUET_COMPRESSION_SNAPPY;
        case MS_CODEC_ZSTD:   return CARQUET_COMPRESSION_ZSTD;
        case MS_CODEC_GZIP:   return CARQUET_COMPRESSION_GZIP;
        case MS_CODEC_LZ4:    return CARQUET_COMPRESSION_LZ4;
        default:              return CARQUET_COMPRESSION_ZSTD;
    }
}

/* Materialize @p n events as a single row group. */
carquet_status_t ms_write_group(carquet_writer_t* w, const ms_event_t* ev,
                                size_t n, bool first, const char** errmsg) {
    if (n == 0) return CARQUET_OK;
    carquet_status_t st = CARQUET_ERROR_OUT_OF_MEMORY;
    *errmsg = "out of memory building row group";

    int64_t* ids  = malloc(n * sizeof(int64_t));
    int64_t* tss  = malloc(n * sizeof(int64_t));
    double*  vals = malloc(n * sizeof(double));
    carquet_byte_array_t* hosts   = malloc(n * sizeof(carquet_byte_array_t));
    carquet_byte_array_t* regions = malloc(n * sizeof(carquet_byte_array_t));
    carquet_byte_array_t* metrics = malloc(n * sizeof(carquet_byte_array_t));
    int32_t* errs   = malloc(n * sizeof(int32_t));   /* packed non-null */
    int16_t* errdef = malloc(n * sizeof(int16_t));   /* one per row */
    uint8_t* sessions = malloc(n * MS_SESSION_ID_LEN);

    if (!ids || !tss || !vals || !hosts || !regions || !metrics ||
        !errs || !errdef || !sessions) goto done;

    size_t err_count = 0;
    for (size_t i = 0; i < n; i++) {
        const ms_event_t* e = &ev[i];
        const char* host   = e->host   ? e->host   : "";
        const char* region = e->region ? e->region : "";
        const char* metric = e->metric ? e->metric : "";
        ids[i]  = e->event_id;
        tss[i]  = e->ts_micros;
        vals[i] = e->value;
        hosts[i].data   = (uint8_t*)host;   hosts[i].length   = (int32_t)strlen(host);
        regions[i].data = (uint8_t*)region; regions[i].length = (int32_t)strlen(region);
        metrics[i].data = (uint8_t*)metric; metrics[i].length = (int32_t)strlen(metric);
        if (e->has_error) { errs[err_count++] = e->error_code; errdef[i] = 1; }
        else              { errdef[i] = 0; }
        memcpy(sessions + i * MS_SESSION_ID_LEN, e->session_id, MS_SESSION_ID_LEN);
    }

    if (!first) {
        st = carquet_writer_new_row_group(w);
        if (st != CARQUET_OK) { *errmsg = "new_row_group failed"; goto done; }
    }

    if ((st = carquet_writer_write_batch(w, COL_EVENT_ID, ids, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_TS, tss, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_HOST, hosts, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_REGION, regions, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_METRIC, metrics, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_VALUE, vals, n, NULL, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_ERROR, errs, n, errdef, NULL)) ||
        (st = carquet_writer_write_batch(w, COL_SESSION, sessions, n, NULL, NULL))) {
        *errmsg = carquet_status_string(st);
        goto done;
    }
    st = CARQUET_OK;

done:
    free(ids); free(tss); free(vals);
    free(hosts); free(regions); free(metrics);
    free(errs); free(errdef); free(sessions);
    return st;
}

/* ------------------------------------------------------------------------- */
/* Writer                                                                     */
/* ------------------------------------------------------------------------- */

/* Deep-copied event: strings are owned so the caller's buffers can go away. */
typedef struct buffered_event {
    int64_t event_id;
    int64_t ts_micros;
    char*   host;
    char*   region;
    char*   metric;
    double  value;
    bool    has_error;
    int32_t error_code;
    uint8_t session_id[MS_SESSION_ID_LEN];
} buffered_event_t;

struct metricstore_writer {
    carquet_writer_t* writer;
    carquet_schema_t* schema;
    ms_writer_config_t cfg;

    buffered_event_t* rows;   /* pending, not yet flushed to a row group */
    size_t            count;
    size_t            capacity;
    bool              first_group; /* have we written the initial group yet? */
};

void ms_writer_config_init(ms_writer_config_t* cfg) {
    if (!cfg) return;
    cfg->codec = MS_CODEC_ZSTD;
    cfg->compression_level = 0;
    cfg->rows_per_group = 50000;
    cfg->bloom_filters = true;
    cfg->page_index = true;
    cfg->statistics = true;
    cfg->arrow_schema = true;
}

carquet_schema_t* ms_build_schema(char* errbuf, size_t errlen) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    if (!s) { ms_set_err(errbuf, errlen, err.message); return NULL; }

    carquet_logical_type_t str  = { .id = CARQUET_LOGICAL_STRING };
    carquet_logical_type_t ts   = { .id = CARQUET_LOGICAL_TIMESTAMP };
    ts.params.timestamp.unit = CARQUET_TIME_UNIT_MICROS;
    ts.params.timestamp.is_adjusted_to_utc = true;
    carquet_logical_type_t i32  = { .id = CARQUET_LOGICAL_INTEGER };
    i32.params.integer.bit_width = 32;
    i32.params.integer.is_signed = true;
    carquet_logical_type_t uuid = { .id = CARQUET_LOGICAL_UUID };

    carquet_schema_add_column(s, "event_id", CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "ts",       CARQUET_PHYSICAL_INT64,  &ts,  CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "host",     CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "region",   CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "metric",   CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "value",    CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "error_code", CARQUET_PHYSICAL_INT32, &i32, CARQUET_REPETITION_OPTIONAL, 0, 0);
    carquet_schema_add_column(s, "session_id", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, &uuid, CARQUET_REPETITION_REQUIRED, MS_SESSION_ID_LEN, 0);

    /* Arrow-style field labels — surfaced by any Arrow-aware reader. */
    carquet_schema_set_field_metadata(s, COL_TS + 1,     "Label", "Event timestamp (UTC microseconds)");
    carquet_schema_set_field_metadata(s, COL_VALUE + 1,  "Label", "Metric sample value");
    carquet_schema_set_field_metadata(s, COL_ERROR + 1,  "Label", "Error code (null when healthy)");

    return s;
}

void ms_tune_writer(carquet_writer_t* w, const ms_writer_config_t* cfg) {
    /* Per-column encoding a real consumer would pick for this data shape. */
    carquet_writer_set_column_encoding(w, COL_EVENT_ID, CARQUET_ENCODING_DELTA_BINARY_PACKED);
    carquet_writer_set_column_encoding(w, COL_HOST,     CARQUET_ENCODING_RLE_DICTIONARY);
    carquet_writer_set_column_encoding(w, COL_REGION,   CARQUET_ENCODING_RLE_DICTIONARY);
    carquet_writer_set_column_encoding(w, COL_METRIC,   CARQUET_ENCODING_RLE_DICTIONARY);
    carquet_writer_set_column_encoding(w, COL_VALUE,    CARQUET_ENCODING_BYTE_STREAM_SPLIT);

    /* Per-column compression + page size overrides. */
    carquet_writer_set_column_compression(w, COL_SESSION, CARQUET_COMPRESSION_UNCOMPRESSED, 0);
    carquet_writer_set_column_page_size(w, COL_TS, 256 * 1024);
    carquet_writer_set_column_statistics(w, COL_SESSION, false);

    /* We write ts in ascending order — advertise it. */
    carquet_sorting_column_t sort = { COL_TS, /*descending=*/false, /*nulls_first=*/false };
    carquet_writer_set_sorting_columns(w, &sort, 1);

    /* Cap truncation of variable-length min/max in statistics. */
    carquet_writer_set_max_statistics_size(w, 64);

    if (cfg->bloom_filters) {
        /* Deliberately mix both bloom APIs: event_id uses the ndv/fpp options
           variant (explicit sizing), host/metric use the plain per-column
           setter, region is explicitly disabled. The two APIs compose (carquet
           honors an explicit legacy enable even after the options API takes
           per-column control — see test_writer_extensions:test_bloom_api_mix). */
        carquet_writer_set_column_bloom_filter_options(w, COL_EVENT_ID, true,
                                                       /*ndv=*/100000, /*fpp=*/0.01);
        carquet_writer_set_column_bloom_filter(w, COL_HOST,   true);
        carquet_writer_set_column_bloom_filter(w, COL_METRIC, true);
        carquet_writer_set_column_bloom_filter(w, COL_REGION, false);
    }

    carquet_writer_add_metadata(w, "app", "metricstore");
    carquet_writer_add_metadata(w, "schema_version", "1");
    carquet_writer_add_metadata(w, "carquet_version", carquet_version());
}

/* Fill writer options from a store config. */
static void fill_options(carquet_writer_options_t* opts, const ms_writer_config_t* cfg) {
    carquet_writer_options_init(opts);
    opts->compression         = ms_to_carquet_codec(cfg->codec);
    opts->compression_level   = cfg->compression_level;
    opts->write_statistics    = cfg->statistics;
    opts->write_bloom_filters = cfg->bloom_filters;
    opts->write_page_index    = cfg->page_index;
    opts->write_arrow_schema  = cfg->arrow_schema;
    opts->write_crc           = true;
}

metricstore_writer_t* ms_writer_open(const char* path,
                                     const ms_writer_config_t* cfg,
                                     char* errbuf, size_t errlen) {
    if (!path) { ms_set_err(errbuf, errlen, "path is NULL"); return NULL; }

    metricstore_writer_t* w = calloc(1, sizeof(*w));
    if (!w) { ms_set_err(errbuf, errlen, "out of memory"); return NULL; }
    ms_writer_config_init(&w->cfg);
    if (cfg) w->cfg = *cfg;
    w->first_group = true;

    w->schema = ms_build_schema(errbuf, errlen);
    if (!w->schema) { free(w); return NULL; }

    carquet_writer_options_t opts;
    fill_options(&opts, &w->cfg);

    carquet_error_t err = CARQUET_ERROR_INIT;
    w->writer = carquet_writer_create(path, w->schema, &opts, &err);
    if (!w->writer) {
        ms_set_err(errbuf, errlen, err.message);
        carquet_schema_free(w->schema);
        free(w);
        return NULL;
    }

    ms_tune_writer(w->writer, &w->cfg);

    /* Sanity: writer must expose the same leaf-column count as our schema. */
    if (carquet_writer_num_columns(w->writer) != COL_COUNT) {
        ms_set_err(errbuf, errlen, "unexpected writer column count");
        carquet_writer_abort(w->writer);
        carquet_schema_free(w->schema);
        free(w);
        return NULL;
    }

    return w;
}

static char* dup_str(const char* s) {
    if (!s) s = "";
    size_t n = strlen(s) + 1;
    char* p = malloc(n);
    if (p) memcpy(p, s, n);
    return p;
}

static void free_row(buffered_event_t* r) {
    free(r->host);
    free(r->region);
    free(r->metric);
    r->host = r->region = r->metric = NULL;
}

/* Materialize the buffered rows as a row group and write it. Reuses the
   shared ms_write_group() builder by exposing the buffered rows as a
   temporary ms_event_t view (strings are borrowed, valid for the call). */
static int flush_group(metricstore_writer_t* w, const char** errmsg) {
    if (w->count == 0) return 0;
    const size_t n = w->count;
    int rc = -1;

    ms_event_t* view = malloc(n * sizeof(*view));
    if (!view) { *errmsg = "out of memory"; goto done; }
    for (size_t i = 0; i < n; i++) {
        buffered_event_t* r = &w->rows[i];
        ms_event_t* e = &view[i];
        e->event_id = r->event_id;
        e->ts_micros = r->ts_micros;
        e->host = r->host; e->region = r->region; e->metric = r->metric;
        e->value = r->value;
        e->has_error = r->has_error; e->error_code = r->error_code;
        memcpy(e->session_id, r->session_id, MS_SESSION_ID_LEN);
    }

    if (ms_write_group(w->writer, view, n, w->first_group, errmsg) == CARQUET_OK)
        rc = 0;
    w->first_group = false;

done:
    free(view);
    for (size_t i = 0; i < w->count; i++) free_row(&w->rows[i]);
    w->count = 0;
    return rc;
}

int ms_writer_append(metricstore_writer_t* w,
                     const ms_event_t* events, size_t count) {
    if (!w || (!events && count)) return -1;

    for (size_t i = 0; i < count; i++) {
        if (w->count == w->capacity) {
            size_t nc = w->capacity ? w->capacity * 2 : 4096;
            buffered_event_t* nr = realloc(w->rows, nc * sizeof(*nr));
            if (!nr) return -1;
            w->rows = nr;
            w->capacity = nc;
        }
        const ms_event_t* e = &events[i];
        buffered_event_t* r = &w->rows[w->count];
        r->event_id  = e->event_id;
        r->ts_micros = e->ts_micros;
        r->value     = e->value;
        r->has_error = e->has_error;
        r->error_code = e->error_code;
        r->host   = dup_str(e->host);
        r->region = dup_str(e->region);
        r->metric = dup_str(e->metric);
        if (!r->host || !r->region || !r->metric) { free_row(r); return -1; }
        memcpy(r->session_id, e->session_id, MS_SESSION_ID_LEN);
        w->count++;

        if ((int64_t)w->count >= w->cfg.rows_per_group) {
            const char* msg = NULL;
            if (flush_group(w, &msg) != 0) return -1;
        }
    }
    return 0;
}

int ms_writer_close(metricstore_writer_t* w) {
    if (!w) return -1;
    int rc = 0;
    const char* msg = NULL;
    if (flush_group(w, &msg) != 0) rc = -1;

    if (rc == 0) {
        if (carquet_writer_close(w->writer) != CARQUET_OK) rc = -1;
    } else {
        carquet_writer_abort(w->writer);
    }
    carquet_schema_free(w->schema);
    free(w->rows);
    free(w);
    return rc;
}

void ms_writer_abort(metricstore_writer_t* w) {
    if (!w) return;
    for (size_t i = 0; i < w->count; i++) free_row(&w->rows[i]);
    carquet_writer_abort(w->writer);
    carquet_schema_free(w->schema);
    free(w->rows);
    free(w);
}

/* ------------------------------------------------------------------------- */
/* Query — projection + predicate pushdown via row-group statistics          */
/* ------------------------------------------------------------------------- */

typedef struct {
    int64_t lo, hi;
} ts_filter_ctx_t;

/* Does row group @p rg overlap the [lo,hi] ts window? Groups with no usable
   statistics are considered overlapping (cannot be pruned). */
static bool ts_group_overlaps(const carquet_reader_t* reader, int32_t rg,
                              int64_t lo, int64_t hi) {
    carquet_column_statistics_t s;
    if (carquet_reader_column_statistics(reader, rg, COL_TS, &s) != CARQUET_OK ||
        !s.has_min_max || !s.min_value || !s.max_value) {
        return true;
    }
    int64_t rg_min = *(const int64_t*)s.min_value;
    int64_t rg_max = *(const int64_t*)s.max_value;
    return !(rg_max < lo || rg_min > hi);
}

/* Pushdown callback: carquet may invoke this several times per row group, so it
   must be side-effect free and purely a function of (reader, rg). */
static bool ts_row_group_filter(const carquet_reader_t* reader,
                                int32_t rg, void* ud) {
    ts_filter_ctx_t* c = (ts_filter_ctx_t*)ud;
    return ts_group_overlaps(reader, rg, c->lo, c->hi);
}

static bool str_eq(const char* want, const char* data, int32_t len) {
    if (!want) return true;
    size_t wl = strlen(want);
    return (size_t)len == wl && memcmp(want, data, wl) == 0;
}

/* Core scan: works on an already-open reader (file- or buffer-backed). */
int ms_scan_reader(carquet_reader_t* reader, const ms_query_t* q,
                   ms_sample_fn on_row, void* ctx,
                   ms_query_result_t* out, char* errbuf, size_t errlen) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    bool seen = false;

    out->row_groups_total = carquet_reader_num_row_groups(reader);

    ts_filter_ctx_t fctx = { q->ts_lo, q->ts_hi };
    bool bounded = (q->ts_lo != MS_TS_MIN || q->ts_hi != MS_TS_MAX);

    /* Deterministically account how many groups the pushdown filter will skip
       (independent of how many times carquet calls the filter callback). */
    if (bounded) {
        for (int32_t rg = 0; rg < out->row_groups_total; rg++) {
            if (!ts_group_overlaps(reader, rg, q->ts_lo, q->ts_hi))
                out->row_groups_pruned++;
        }
    }

    int32_t proj[] = { COL_TS, COL_HOST, COL_METRIC, COL_VALUE };
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 8192;
    cfg.column_indices = proj;
    cfg.num_columns = 4;
    cfg.use_mmap = q->use_mmap;
    /* Only install the pushdown filter when the range is actually bounded. */
    if (bounded) {
        cfg.row_group_filter = ts_row_group_filter;
        cfg.row_group_filter_ctx = &fctx;
    }

    /* Optionally decode columns on a shared worker pool. */
    carquet_thread_pool_t* pool = NULL;
    if (q->parallel) {
        pool = carquet_thread_pool_create(0 /* auto */);
        cfg.thread_pool = pool;
    }

    /* Warm the OS cache for the projected columns of the first row group. */
    (void)carquet_reader_prebuffer(reader, 0, proj, 4, &err);

    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        ms_set_err(errbuf, errlen, err.message);
        if (pool) carquet_thread_pool_destroy(pool);
        carquet_reader_release_prebuffer(reader);
        return -1;
    }

    /* Page-level pruning: when the ts range is bounded, push a RANGE clause so
       whole pages outside the window are skipped before materialization. */
    if (bounded) {
        carquet_filter_clause_t clause;
        memset(&clause, 0, sizeof(clause));
        clause.column_index = COL_TS;
        clause.op = CARQUET_FILTER_RANGE;
        if (q->ts_lo != MS_TS_MIN) { clause.lo = &q->ts_lo; clause.lo_size = sizeof(q->ts_lo); clause.has_lo = true; }
        if (q->ts_hi != MS_TS_MAX) { clause.hi = &q->ts_hi; clause.hi_size = sizeof(q->ts_hi); clause.has_hi = true; }
        (void)carquet_batch_reader_set_page_filter(br, &clause, 1);
    }

    int rc = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        /* Cross-check the batch's own row/column counts. */
        (void)carquet_row_batch_num_columns(batch);
        int64_t batch_rows = carquet_row_batch_num_rows(batch);

        const void *tsd, *hostd, *metricd, *vald;
        const uint8_t *nb;
        int64_t n = 0, tmp;
        carquet_row_batch_column(batch, 0, &tsd, &nb, &n);
        if (n != batch_rows) n = batch_rows; /* stay consistent */
        carquet_row_batch_column(batch, 1, &hostd, &nb, &tmp);
        carquet_row_batch_column(batch, 2, &metricd, &nb, &tmp);
        carquet_row_batch_column(batch, 3, &vald, &nb, &tmp);

        const int64_t* ts = (const int64_t*)tsd;
        const carquet_byte_array_t* host = (const carquet_byte_array_t*)hostd;
        const carquet_byte_array_t* metric = (const carquet_byte_array_t*)metricd;
        const double* value = (const double*)vald;

        out->rows_scanned += n;
        for (int64_t i = 0; i < n; i++) {
            if (ts[i] < q->ts_lo || ts[i] > q->ts_hi) continue;
            if (!str_eq(q->metric, (const char*)metric[i].data, metric[i].length)) continue;
            if (!str_eq(q->host, (const char*)host[i].data, host[i].length)) continue;

            out->rows_matched++;
            out->value_sum += value[i];
            if (!seen || value[i] < out->value_min) out->value_min = value[i];
            if (!seen || value[i] > out->value_max) out->value_max = value[i];
            seen = true;

            if (on_row) {
                ms_sample_t s = {
                    ts[i],
                    (const char*)host[i].data, host[i].length,
                    (const char*)metric[i].data, metric[i].length,
                    value[i]
                };
                on_row(&s, ctx);
            }
        }
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    /* Account for rows the page filter skipped without materializing. */
    out->rows_page_skipped = carquet_batch_reader_rows_skipped(br);

    carquet_batch_reader_free(br);
    carquet_reader_release_prebuffer(reader);
    if (pool) carquet_thread_pool_destroy(pool);
    return rc;
}

int ms_query(const char* path, const ms_query_t* q,
             ms_sample_fn on_row, void* ctx,
             ms_query_result_t* out, char* errbuf, size_t errlen) {
    if (!path || !q || !out) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }
    memset(out, 0, sizeof(*out));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ropts;
    carquet_reader_options_init(&ropts);
    ropts.use_mmap = q->use_mmap;

    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    int rc = ms_scan_reader(reader, q, on_row, ctx, out, errbuf, errlen);
    carquet_reader_close(reader);
    return rc;
}

/* ------------------------------------------------------------------------- */
/* Bloom-filter membership                                                    */
/* ------------------------------------------------------------------------- */

int ms_might_contain_host(const char* path, const char* host,
                          char* errbuf, size_t errlen) {
    if (!path || !host) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    int result = 0; /* definitely absent until a filter says "maybe" */
    bool any_filter = false;
    int32_t rgs = carquet_reader_num_row_groups(reader);
    size_t hlen = strlen(host);

    for (int32_t rg = 0; rg < rgs; rg++) {
        carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(reader, rg, COL_HOST, &err);
        if (!bf) continue; /* no filter on this chunk -> inconclusive */
        any_filter = true;
        if (carquet_bloom_filter_check_bytes(bf, (const uint8_t*)host, hlen)) {
            result = 1;
            carquet_bloom_filter_destroy(bf);
            break;
        }
        carquet_bloom_filter_destroy(bf);
    }

    /* If no bloom filter exists at all, we cannot claim absence. */
    if (!any_filter) result = 1;

    carquet_reader_close(reader);
    return result;
}

/* ------------------------------------------------------------------------- */
/* Introspection                                                              */
/* ------------------------------------------------------------------------- */

int ms_describe(const char* path, FILE* out, char* errbuf, size_t errlen) {
    if (!path) { ms_set_err(errbuf, errlen, "path is NULL"); return -1; }
    if (!out) out = stdout;

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_file_info_t info;
    if (carquet_get_file_info(path, &info, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        return -1;
    }

    fprintf(out, "MetricStore file: %s\n", path);
    fprintf(out, "  size=%lld bytes  rows=%lld  row_groups=%d  columns=%d\n",
            (long long)info.file_size, (long long)info.num_rows,
            info.num_row_groups, info.num_columns);
    fprintf(out, "  created_by=%s\n", info.created_by[0] ? info.created_by : "(unset)");

    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t ncols = carquet_schema_num_columns(schema);

    fprintf(out, "  reader: rows=%lld columns=%d row_groups=%d mmap=%s\n",
            (long long)carquet_reader_num_rows(reader),
            carquet_reader_num_columns(reader),
            carquet_reader_num_row_groups(reader),
            carquet_reader_is_mmap(reader) ? "yes" : "no");
    fprintf(out, "  schema elements=%d leaf columns=%d  (find 'value' -> col %d)\n",
            carquet_schema_num_elements(schema), ncols,
            carquet_schema_find_column(schema, "value"));

    /* Predicate-pushdown preview independent of the batch reader: which row
       groups could hold ts >= the midpoint of the file's range? */
    {
        carquet_column_statistics_t s0;
        if (carquet_reader_column_statistics(reader, 0, COL_TS, &s0) == CARQUET_OK &&
            s0.has_min_max && s0.min_value) {
            int64_t probe = *(const int64_t*)s0.min_value;
            int32_t idx[64];
            int32_t nmatch = carquet_reader_filter_row_groups(
                reader, COL_TS, CARQUET_COMPARE_GE, &probe, sizeof(probe), idx, 64);
            bool rg0 = false;
            carquet_reader_row_group_matches(
                reader, 0, COL_TS, CARQUET_COMPARE_GE, &probe, sizeof(probe), &rg0);
            fprintf(out, "  pushdown: %d/%d row groups have ts>=min (rg0 matches=%s)\n",
                    nmatch, carquet_reader_num_row_groups(reader), rg0 ? "yes" : "no");
        }
    }

    /* Row-group level metadata. */
    fprintf(out, "\n  Row groups:\n");
    for (int32_t g = 0; g < carquet_reader_num_row_groups(reader); g++) {
        carquet_row_group_metadata_t rg;
        if (carquet_reader_row_group_metadata(reader, g, &rg) == CARQUET_OK) {
            fprintf(out, "    [%d] rows=%lld  %lld -> %lld B (%.1fx)\n", g,
                    (long long)rg.num_rows, (long long)rg.total_byte_size,
                    (long long)rg.total_compressed_size,
                    rg.total_compressed_size
                        ? (double)rg.total_byte_size / (double)rg.total_compressed_size
                        : 0.0);
        }
    }

    fprintf(out, "\n  Key-value metadata:\n");
    int32_t nkv = carquet_reader_num_metadata(reader);
    for (int32_t i = 0; i < nkv; i++) {
        const char *k = NULL, *v = NULL;
        carquet_reader_get_metadata(reader, i, &k, &v);
        fprintf(out, "    %s = %s\n", k ? k : "?", v ? v : "(null)");
    }
    const char* app = carquet_reader_find_metadata(reader, "app");
    fprintf(out, "    (lookup 'app' -> %s)\n", app ? app : "not found");

    fprintf(out, "\n  Columns:\n");
    for (int32_t c = 0; c < ncols; c++) {
        const char* name = carquet_schema_column_name(schema, c);
        carquet_physical_type_t pt = carquet_schema_column_type(schema, c);

        /* Schema-node view: path, nesting, def/rep levels. */
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, c + 1);
        const char* pathbuf[8];
        int32_t depth = carquet_schema_column_path(schema, c, pathbuf, 8);
        char pathstr[128]; pathstr[0] = '\0';
        for (int32_t d = 0; d < depth; d++) {
            if (d) strncat(pathstr, ".", sizeof(pathstr) - strlen(pathstr) - 1);
            strncat(pathstr, pathbuf[d], sizeof(pathstr) - strlen(pathstr) - 1);
        }

        carquet_column_chunk_metadata_t m;
        carquet_reader_column_chunk_metadata(reader, 0, c, &m);

        fprintf(out, "    [%d] %-11s %-22s codec=%-6s values=%lld comp=%lld/%lld B",
                c, name, carquet_physical_type_name(pt),
                carquet_compression_name(m.codec),
                (long long)m.num_values,
                (long long)m.total_compressed_size,
                (long long)m.total_uncompressed_size);
        if (m.num_encodings > 0) {
            fprintf(out, " enc=");
            for (int32_t e = 0; e < m.num_encodings; e++)
                fprintf(out, "%s%s", e ? "," : "", carquet_encoding_name(m.encodings[e]));
        }
        fprintf(out, " %s%s\n",
                m.has_bloom_filter ? "[bloom]" : "",
                m.has_column_index ? "[index]" : "");

        if (node) {
            const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
            fprintf(out, "         path=%s node=%s ptype=%d ltype=%d flen=%d "
                    "leaf=%s rep=%d def_level=%d/%d rep_level=%d/%d\n",
                    pathstr[0] ? pathstr : name,
                    carquet_schema_node_name(node),
                    (int)carquet_schema_node_physical_type(node),
                    lt ? (int)lt->id : -1,
                    carquet_schema_node_type_length(node),
                    carquet_schema_node_is_leaf(node) ? "yes" : "no",
                    (int)carquet_schema_node_repetition(node),
                    (int)carquet_schema_node_max_def_level(node),
                    (int)carquet_schema_max_def_level(schema, c),
                    (int)carquet_schema_node_max_rep_level(node),
                    (int)carquet_schema_max_rep_level(schema, c));
        }

        /* Arrow per-field metadata (the labels we set at write time). */
        int32_t nfm = carquet_reader_column_num_metadata(reader, c);
        for (int32_t f = 0; f < nfm; f++) {
            const char *fk = NULL, *fv = NULL;
            if (carquet_reader_column_get_metadata(reader, c, f, &fk, &fv) == CARQUET_OK)
                fprintf(out, "         field-meta: %s = %s\n", fk ? fk : "?", fv ? fv : "");
        }
        const char* label = carquet_reader_column_find_metadata(reader, c, "Label");
        if (label) fprintf(out, "         label = \"%s\"\n", label);

        carquet_column_statistics_t s;
        if (carquet_reader_column_statistics(reader, 0, c, &s) == CARQUET_OK &&
            s.has_min_max) {
            if (pt == CARQUET_PHYSICAL_INT64 && s.min_value && s.max_value) {
                fprintf(out, "         stats: min=%lld max=%lld nulls=%lld\n",
                        (long long)*(const int64_t*)s.min_value,
                        (long long)*(const int64_t*)s.max_value,
                        (long long)s.null_count);
            } else if (pt == CARQUET_PHYSICAL_DOUBLE && s.min_value && s.max_value) {
                fprintf(out, "         stats: min=%g max=%g nulls=%lld\n",
                        *(const double*)s.min_value, *(const double*)s.max_value,
                        (long long)s.null_count);
            } else {
                fprintf(out, "         stats: nulls=%lld\n", (long long)s.null_count);
            }
        }

        /* Page index (column index + offset index) for the first row group. */
        carquet_column_index_t* ci = carquet_reader_get_column_index(reader, 0, c, &err);
        carquet_offset_index_t* oi = carquet_reader_get_offset_index(reader, 0, c, &err);
        if (ci) {
            int32_t np = carquet_column_index_num_pages(ci);
            fprintf(out, "         page index: %d page(s), boundary_order=%d",
                    np, carquet_column_index_boundary_order(ci));
            if (oi) fprintf(out, ", offset pages=%d", carquet_offset_index_num_pages(oi));
            fprintf(out, "\n");
            if (np > 0) {
                carquet_page_stats_t ps;
                carquet_page_location_t loc;
                if (carquet_column_index_get_page_stats(ci, 0, &ps) == CARQUET_OK)
                    fprintf(out, "           page0: nulls=%lld null_page=%s\n",
                            (long long)ps.null_count, ps.is_null_page ? "yes" : "no");
                if (oi && carquet_offset_index_get_page_location(oi, 0, &loc) == CARQUET_OK)
                    fprintf(out, "           page0: offset=%lld size=%d first_row=%lld\n",
                            (long long)loc.offset, loc.compressed_size,
                            (long long)loc.first_row_index);
            }
        }
        if (ci) carquet_column_index_free(ci);
        if (oi) carquet_offset_index_free(oi);
    }

    carquet_reader_close(reader);
    return 0;
}

int ms_validate(const char* path, char* errbuf, size_t errlen) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_validate_file(path, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        return -1;
    }
    return 0;
}

const char* ms_backend_version(void) {
    return carquet_version();
}
