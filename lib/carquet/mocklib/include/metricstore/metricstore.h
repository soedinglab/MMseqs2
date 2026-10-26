/**
 * @file metricstore.h
 * @brief MetricStore — a small telemetry/observability store built on carquet.
 *
 * MetricStore is a *mockup* downstream library: it does not re-implement any
 * Parquet logic itself. Instead it models a realistic domain (time-series
 * telemetry events) and drives a broad slice of the public carquet API the way
 * a real consumer would — schema construction with logical types, per-column
 * encoding/compression tuning, bloom filters, page indexes and statistics on
 * the write side; column projection, predicate pushdown, bloom-filter
 * membership tests and metadata introspection on the read side.
 *
 * Its purpose is real-world integration testing of carquet: exercise the API
 * surface end to end, round-trip data, and assert the results.
 *
 * The API here is deliberately carquet-agnostic (no carquet types leak through)
 * so it reads like an independent product built on top of the library.
 */
#ifndef METRICSTORE_H
#define METRICSTORE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Length in bytes of an opaque session identifier (stored as UUID). */
#define MS_SESSION_ID_LEN 16

/** Sentinels for an unbounded timestamp range in ::ms_query_t. */
#define MS_TS_MIN INT64_MIN
#define MS_TS_MAX INT64_MAX

/** Compression codec selector, mapped internally to a carquet codec. */
typedef enum ms_codec {
    MS_CODEC_NONE = 0,
    MS_CODEC_SNAPPY,
    MS_CODEC_ZSTD,
    MS_CODEC_GZIP,
    MS_CODEC_LZ4
} ms_codec_t;

/** A single telemetry sample (one logical row). */
typedef struct ms_event {
    int64_t     event_id;                     /**< Monotonic-ish event id.       */
    int64_t     ts_micros;                     /**< Unix timestamp, microseconds. */
    const char* host;                          /**< Emitting host (required).     */
    const char* region;                        /**< Datacenter region (required). */
    const char* metric;                        /**< Metric name (required).       */
    double      value;                         /**< Sample value.                 */
    bool        has_error;                     /**< false => error_code is NULL.  */
    int32_t     error_code;                    /**< Optional error code.          */
    uint8_t     session_id[MS_SESSION_ID_LEN]; /**< Opaque session id (UUID).     */
} ms_event_t;

/** Writer tuning. Initialise with ::ms_writer_config_init, then override. */
typedef struct ms_writer_config {
    ms_codec_t codec;             /**< Default codec for all columns.        */
    int32_t    compression_level; /**< 0 = codec default.                    */
    int64_t    rows_per_group;    /**< Flush a new row group past this many.  */
    bool       bloom_filters;     /**< Emit bloom filters on id/host/metric. */
    bool       page_index;        /**< Emit page (column+offset) index.      */
    bool       statistics;        /**< Emit column statistics.               */
    bool       arrow_schema;      /**< Embed ARROW:schema for field labels.  */
} ms_writer_config_t;

/** Opaque append-only writer handle. */
typedef struct metricstore_writer metricstore_writer_t;

/** Populate @p cfg with sensible defaults (zstd, bloom+index+stats on). */
void ms_writer_config_init(ms_writer_config_t* cfg);

/**
 * Open a store for writing at @p path. On failure returns NULL and, if
 * @p errbuf is non-NULL, writes a message.
 */
metricstore_writer_t* ms_writer_open(const char* path,
                                     const ms_writer_config_t* cfg,
                                     char* errbuf, size_t errlen);

/** Append @p count events. Returns 0 on success, -1 on error. */
int ms_writer_append(metricstore_writer_t* w,
                     const ms_event_t* events, size_t count);

/** Flush, finalize the file and free the writer. Returns 0 on success. */
int ms_writer_close(metricstore_writer_t* w);

/** Discard any buffered data and free the writer without finalizing. */
void ms_writer_abort(metricstore_writer_t* w);

/** A projected row delivered to a scan callback. Pointers are borrowed. */
typedef struct ms_sample {
    int64_t     ts_micros;
    const char* host;
    int32_t     host_len;
    const char* metric;
    int32_t     metric_len;
    double      value;
} ms_sample_t;

/** Predicate for a scan. NULL string / sentinel ts means "unbounded". */
typedef struct ms_query {
    const char* metric;  /**< Match metric exactly, or NULL for any.  */
    const char* host;    /**< Match host exactly, or NULL for any.    */
    int64_t     ts_lo;   /**< Inclusive lower bound (MS_TS_MIN = any).*/
    int64_t     ts_hi;   /**< Inclusive upper bound (MS_TS_MAX = any).*/
    bool        use_mmap;/**< Read via memory-mapped I/O.             */
    bool        parallel;/**< Use a worker pool for column decode.    */
} ms_query_t;

/** Aggregates and pushdown accounting returned by ::ms_query. */
typedef struct ms_query_result {
    int64_t rows_scanned;      /**< Rows materialized by the batch reader.  */
    int64_t rows_matched;      /**< Rows passing the predicate.             */
    int32_t row_groups_total;  /**< Row groups in the file.                 */
    int32_t row_groups_pruned; /**< Row groups skipped via statistics.      */
    int64_t rows_page_skipped; /**< Rows skipped by the page filter.        */
    double  value_sum;         /**< Sum of value over matched rows.         */
    double  value_min;         /**< Min value over matched rows.            */
    double  value_max;         /**< Max value over matched rows.            */
} ms_query_result_t;

/** Called once per matched row during ::ms_query (may be NULL). */
typedef void (*ms_sample_fn)(const ms_sample_t* sample, void* ctx);

/**
 * Scan @p path applying @p q, invoking @p on_row for each match and filling
 * @p out with aggregates. Returns 0 on success, -1 on error.
 */
int ms_query(const char* path, const ms_query_t* q,
             ms_sample_fn on_row, void* ctx,
             ms_query_result_t* out, char* errbuf, size_t errlen);

/**
 * Bloom-filter membership probe on the host column.
 * Returns 1 if @p host might be present, 0 if it is definitely absent,
 * and -1 on error.
 */
int ms_might_contain_host(const char* path, const char* host,
                          char* errbuf, size_t errlen);

/** Print a human-readable description (schema, chunks, stats). Returns 0/-1. */
int ms_describe(const char* path, FILE* out, char* errbuf, size_t errlen);

/** Structurally validate the file. Returns 0 if valid, -1 otherwise. */
int ms_validate(const char* path, char* errbuf, size_t errlen);

/** Version string of the carquet library backing this store. */
const char* ms_backend_version(void);

/* ========================================================================= */
/* Extended surface — broader carquet coverage for integration testing.      */
/* ========================================================================= */

/**
 * Append @p count events as a new row group to an existing store at @p path
 * (opens the file in append mode). Returns 0 on success, -1 on error.
 */
int ms_append(const char* path, const ms_event_t* events, size_t count,
              char* errbuf, size_t errlen);

/**
 * Serialize @p count events to an in-memory Parquet blob. On success the caller
 * owns @p *buf and must release it with ::ms_free_buffer. Returns 0 / -1.
 */
int ms_pack(const ms_event_t* events, size_t count, const ms_writer_config_t* cfg,
            void** buf, size_t* size, char* errbuf, size_t errlen);

/** Run a query directly against an in-memory blob from ::ms_pack. */
int ms_query_buffer(const void* buf, size_t size, const ms_query_t* q,
                    ms_query_result_t* out, char* errbuf, size_t errlen);

/** Release a buffer returned by ::ms_pack. */
void ms_free_buffer(void* buf);

/**
 * Read up to @p max event ids from @p row_group using the low-level column
 * reader (streaming, with skip). Returns the count read, or -1 on error.
 */
int64_t ms_read_event_ids(const char* path, int32_t row_group, int64_t skip,
                          int64_t* out, int64_t max, char* errbuf, size_t errlen);

/**
 * Bloom-filter membership probe on the event_id column (INT64 filter).
 * Returns 1 (maybe present), 0 (definitely absent), or -1 on error.
 */
int ms_might_contain_event_id(const char* path, int64_t id,
                              char* errbuf, size_t errlen);

/**
 * Self-contained exercise of nested data: writes a small file with an id
 * column and a list<string> "tags" column, reads it back reconstructing list
 * boundaries, and verifies the reconstruction. Returns 0 on success, -1 on
 * error (message in @p errbuf).
 */
int ms_nested_selfcheck(const char* path, char* errbuf, size_t errlen);

/**
 * Read one row batch from @p src, bridge it through the Arrow C Data Interface,
 * and import it straight back into a new Parquet file @p dst. Returns the number
 * of rows moved, or -1 on error.
 */
int64_t ms_arrow_roundtrip(const char* src, const char* dst,
                           char* errbuf, size_t errlen);

/** Print backend diagnostics (version components, detected CPU features). */
void ms_print_diagnostics(FILE* out);

/** Release backend thread-local resources (call once before program exit). */
void ms_shutdown(void);

#ifdef __cplusplus
}
#endif

#endif /* METRICSTORE_H */
