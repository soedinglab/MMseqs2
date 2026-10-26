/**
 * @file write_sample.c
 * @brief Emit a MetricStore Parquet file for external inspection.
 *
 * Writes a small, fully-deterministic dataset (default 500 rows across a few
 * hosts/regions/metrics, split into a couple of row groups) and leaves it on
 * disk so it can be opened with pyarrow, parquet-tools, DuckDB, etc. Then
 * prints a description of what was written.
 *
 * Usage: write_sample [output.parquet] [num_rows]
 */
#include "metricstore/metricstore.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* HOSTS[]   = { "web-01", "web-02", "db-01" };
static const char* REGIONS[] = { "us-east", "eu-central" };
static const char* METRICS[] = { "cpu.util", "mem.used", "req.latency" };
#define NELEM(a) ((int)(sizeof(a) / sizeof((a)[0])))

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "sample.parquet";
    long rows = (argc > 2) ? strtol(argv[2], NULL, 10) : 500;
    if (rows <= 0) rows = 500;

    char errbuf[512] = {0};

    ms_writer_config_t cfg;
    ms_writer_config_init(&cfg);
    cfg.codec = MS_CODEC_ZSTD;
    cfg.rows_per_group = 200; /* force multiple row groups even for small N */

    metricstore_writer_t* w = ms_writer_open(path, &cfg, errbuf, sizeof(errbuf));
    if (!w) { fprintf(stderr, "open: %s\n", errbuf); return 1; }

    ms_event_t* buf = malloc((size_t)rows * sizeof(*buf));
    int64_t base_ts = 1700000000000000LL;
    for (long i = 0; i < rows; i++) {
        ms_event_t* e = &buf[i];
        e->event_id  = i;
        e->ts_micros = base_ts + i * 1000000LL; /* 1 second apart */
        e->host   = HOSTS[i % NELEM(HOSTS)];
        e->region = REGIONS[i % NELEM(REGIONS)];
        e->metric = METRICS[i % NELEM(METRICS)];
        e->value  = (double)((i * 37) % 1000) + 0.5;
        e->has_error = (i % 10 == 0);
        e->error_code = e->has_error ? (int32_t)(500 + (i % 4)) : 0;
        for (int b = 0; b < MS_SESSION_ID_LEN; b++)
            e->session_id[b] = (uint8_t)(i + b);
    }
    if (ms_writer_append(w, buf, (size_t)rows) != 0) {
        fprintf(stderr, "append failed\n"); free(buf); ms_writer_abort(w); return 1;
    }
    free(buf);
    if (ms_writer_close(w) != 0) { fprintf(stderr, "close failed\n"); return 1; }

    printf("Wrote %ld rows to %s (backend carquet %s)\n\n",
           rows, path, ms_backend_version());
    ms_describe(path, stdout, errbuf, sizeof(errbuf));
    printf("\nInspect with e.g.:\n");
    printf("  python3 -c \"import pyarrow.parquet as pq; print(pq.read_table('%s').to_pandas().head())\"\n",
           path);
    return 0;
}
