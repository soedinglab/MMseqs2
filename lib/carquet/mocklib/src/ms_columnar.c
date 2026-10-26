/**
 * @file ms_columnar.c
 * @brief Buffer I/O, append, low-level column reads and int64 bloom probes.
 *
 * Exercises the carquet surface a consumer reaches for beyond the basic
 * open/write/query path: in-memory writer/reader, append mode, the streaming
 * column reader (skip/remaining/has_next), and typed bloom-filter checks.
 */
#include "ms_internal.h"

#include <stdlib.h>
#include <string.h>

/* -- append a new row group to an existing store ------------------------- */

int ms_append(const char* path, const ms_event_t* events, size_t count,
              char* errbuf, size_t errlen) {
    if (!path || (!events && count)) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }

    carquet_schema_t* schema = ms_build_schema(errbuf, errlen);
    if (!schema) return -1;

    ms_writer_config_t cfg;
    ms_writer_config_init(&cfg);
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = ms_to_carquet_codec(cfg.codec);
    opts.write_statistics = cfg.statistics;
    opts.write_bloom_filters = cfg.bloom_filters;
    opts.write_page_index = cfg.page_index;
    opts.write_arrow_schema = cfg.arrow_schema;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_writer_t* w = carquet_writer_open_append(path, schema, &opts, &err);
    if (!w) { ms_set_err(errbuf, errlen, err.message); carquet_schema_free(schema); return -1; }

    ms_tune_writer(w, &cfg);

    /* Append always starts a fresh row group (first=false). */
    const char* msg = NULL;
    int rc = 0;
    if (ms_write_group(w, events, count, /*first=*/false, &msg) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, msg ? msg : "append write failed");
        carquet_writer_abort(w);
        rc = -1;
    } else if (carquet_writer_close(w) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, "append close failed");
        rc = -1;
    }
    carquet_schema_free(schema);
    return rc;
}

/* -- in-memory (buffer) writer / reader round trip ----------------------- */

int ms_pack(const ms_event_t* events, size_t count, const ms_writer_config_t* cfg,
            void** buf, size_t* size, char* errbuf, size_t errlen) {
    if (!events && count) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }
    if (!buf || !size) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }
    *buf = NULL; *size = 0;

    ms_writer_config_t local;
    ms_writer_config_init(&local);
    if (cfg) local = *cfg;

    carquet_schema_t* schema = ms_build_schema(errbuf, errlen);
    if (!schema) return -1;

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = ms_to_carquet_codec(local.codec);
    opts.write_statistics = local.statistics;
    opts.write_bloom_filters = local.bloom_filters;
    opts.write_page_index = local.page_index;
    opts.write_arrow_schema = local.arrow_schema;

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_writer_t* w = carquet_writer_create_buffer(schema, &opts, &err);
    if (!w) { ms_set_err(errbuf, errlen, err.message); carquet_schema_free(schema); return -1; }

    ms_tune_writer(w, &local);

    const char* msg = NULL;
    if (ms_write_group(w, events, count, /*first=*/true, &msg) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, msg ? msg : "pack write failed");
        carquet_writer_abort(w);
        carquet_schema_free(schema);
        return -1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, "pack close failed");
        carquet_schema_free(schema);
        return -1;
    }

    /* Take an owning copy of the writer's internal buffer. */
    void* internal = NULL; size_t isz = 0;
    int rc = -1;
    if (carquet_writer_get_buffer(w, &internal, &isz) == CARQUET_OK && internal) {
        void* copy = malloc(isz);
        if (copy) { memcpy(copy, internal, isz); *buf = copy; *size = isz; rc = 0; }
        else ms_set_err(errbuf, errlen, "out of memory copying buffer");
        free(internal);
    } else {
        ms_set_err(errbuf, errlen, "get_buffer failed");
    }
    carquet_schema_free(schema);
    return rc;
}

void ms_free_buffer(void* buf) { free(buf); }

int ms_query_buffer(const void* buf, size_t size, const ms_query_t* q,
                    ms_query_result_t* out, char* errbuf, size_t errlen) {
    if (!buf || !q || !out) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }
    memset(out, 0, sizeof(*out));

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open_buffer(buf, size, NULL, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    int rc = ms_scan_reader(reader, q, NULL, NULL, out, errbuf, errlen);
    carquet_reader_close(reader);
    return rc;
}

/* -- low-level streaming column reader ----------------------------------- */

int64_t ms_read_event_ids(const char* path, int32_t row_group, int64_t skip,
                          int64_t* out, int64_t max, char* errbuf, size_t errlen) {
    if (!path || !out || max <= 0) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    carquet_column_reader_t* col =
        carquet_reader_get_column(reader, row_group, COL_EVENT_ID, &err);
    if (!col) { ms_set_err(errbuf, errlen, err.message); carquet_reader_close(reader); return -1; }

    if (skip > 0) carquet_column_skip(col, skip);

    int64_t got = 0;
    if (carquet_column_has_next(col)) {
        /* Use the detailed variant so its error channel is exercised. */
        carquet_error_t rerr = CARQUET_ERROR_INIT;
        got = carquet_column_read_batch_ex(col, out, max, NULL, NULL, &rerr);
        if (got < 0) {
            char b[256];
            carquet_error_format(&rerr, b, sizeof(b));
            ms_set_err(errbuf, errlen, b);
            carquet_column_reader_free(col);
            carquet_reader_close(reader);
            return -1;
        }
    }

    /* Touch column_remaining so the accessor is exercised too. */
    (void)carquet_column_remaining(col);

    carquet_column_reader_free(col);
    carquet_reader_close(reader);
    return got;
}

/* -- typed bloom membership on the event_id (INT64) column --------------- */

int ms_might_contain_event_id(const char* path, int64_t id,
                              char* errbuf, size_t errlen) {
    if (!path) { ms_set_err(errbuf, errlen, "path is NULL"); return -1; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);
    if (!reader) { ms_set_err(errbuf, errlen, err.message); return -1; }

    int result = 0;
    bool any = false;
    int32_t rgs = carquet_reader_num_row_groups(reader);
    for (int32_t rg = 0; rg < rgs; rg++) {
        carquet_bloom_filter_t* bf =
            carquet_reader_get_bloom_filter(reader, rg, COL_EVENT_ID, &err);
        if (!bf) continue;
        any = true;
        (void)carquet_bloom_filter_size(bf);
        if (carquet_bloom_filter_check_i64(bf, id)) {
            result = 1;
            carquet_bloom_filter_destroy(bf);
            break;
        }
        carquet_bloom_filter_destroy(bf);
    }
    if (!any) result = 1;
    carquet_reader_close(reader);
    return result;
}
