/**
 * @file ms_internal.h
 * @brief Shared internals across MetricStore translation units. Not installed.
 */
#ifndef MS_INTERNAL_H
#define MS_INTERNAL_H

#include "metricstore/metricstore.h"
#include <carquet/carquet.h>

/* Stable leaf-column indices for the MetricStore schema. */
enum {
    COL_EVENT_ID = 0,
    COL_TS,
    COL_HOST,
    COL_REGION,
    COL_METRIC,
    COL_VALUE,
    COL_ERROR,
    COL_SESSION,
    COL_COUNT
};

/* Copy a message into a caller buffer (both may be NULL). */
void ms_set_err(char* buf, size_t len, const char* msg);

/* Map the public codec enum to a carquet codec. */
carquet_compression_t ms_to_carquet_codec(ms_codec_t c);

/* Build the canonical MetricStore schema (caller frees). */
carquet_schema_t* ms_build_schema(char* errbuf, size_t errlen);

/* Apply the standard per-column tuning to a freshly created writer. */
void ms_tune_writer(carquet_writer_t* w, const ms_writer_config_t* cfg);

/* Write @p n events as a single row group into an open writer. When @p first
   is false a new row group is started first. On error returns non-OK and sets
   *errmsg to a static string. */
carquet_status_t ms_write_group(carquet_writer_t* w, const ms_event_t* ev,
                                size_t n, bool first, const char** errmsg);

/* Core scan over an already-open reader; fills @p out. Does not close reader. */
int ms_scan_reader(carquet_reader_t* reader, const ms_query_t* q,
                   ms_sample_fn on_row, void* ctx,
                   ms_query_result_t* out, char* errbuf, size_t errlen);

#endif /* MS_INTERNAL_H */
