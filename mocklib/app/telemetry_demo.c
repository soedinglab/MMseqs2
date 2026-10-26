/**
 * @file telemetry_demo.c
 * @brief End-to-end MetricStore scenario: ingest synthetic telemetry, then
 *        introspect and query it. Exercises the carquet-backed store the way a
 *        real observability pipeline would.
 *
 * Usage: telemetry_demo [output.parquet] [num_events]
 */
#include "metricstore/metricstore.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* HOSTS[]   = { "web-01", "web-02", "web-03", "db-01", "cache-01" };
static const char* REGIONS[] = { "us-east", "us-west", "eu-central", "ap-south" };
static const char* METRICS[] = { "cpu.util", "mem.used", "req.latency", "disk.io", "net.rx" };

#define NELEM(a) ((int)(sizeof(a) / sizeof((a)[0])))

/* Deterministic LCG so runs are reproducible without touching global rand(). */
static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static uint32_t rng_next(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (uint32_t)(rng_state >> 33);
}

static void row_printer(const ms_sample_t* s, void* ctx) {
    int* shown = (int*)ctx;
    if (*shown < 5) {
        printf("      ts=%lld host=%.*s metric=%.*s value=%.3f\n",
               (long long)s->ts_micros,
               s->host_len, s->host, s->metric_len, s->metric, s->value);
        (*shown)++;
    }
}

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "/tmp/metricstore_demo.parquet";
    long num_events  = (argc > 2) ? strtol(argv[2], NULL, 10) : 200000;
    if (num_events <= 0) num_events = 200000;

    char errbuf[512] = {0};

    printf("MetricStore demo\n");
    ms_print_diagnostics(stdout);
    printf("Ingesting %ld telemetry events -> %s\n\n", num_events, path);

    ms_writer_config_t cfg;
    ms_writer_config_init(&cfg);
    cfg.codec = MS_CODEC_ZSTD;
    cfg.rows_per_group = 40000; /* several row groups so pushdown is visible */

    metricstore_writer_t* w = ms_writer_open(path, &cfg, errbuf, sizeof(errbuf));
    if (!w) { fprintf(stderr, "open failed: %s\n", errbuf); return 1; }

    /* Ingest in chunks the way a streaming collector would. */
    enum { CHUNK = 4096 };
    ms_event_t* buf = malloc(CHUNK * sizeof(*buf));
    int64_t base_ts = 1700000000000000LL; /* ~2023-11-14 in micros */
    long produced = 0;
    while (produced < num_events) {
        int m = 0;
        for (; m < CHUNK && produced < num_events; m++, produced++) {
            ms_event_t* e = &buf[m];
            e->event_id  = produced;
            e->ts_micros = base_ts + produced * 1000LL; /* 1ms apart, monotonic */
            e->host   = HOSTS[rng_next() % NELEM(HOSTS)];
            e->region = REGIONS[rng_next() % NELEM(REGIONS)];
            e->metric = METRICS[rng_next() % NELEM(METRICS)];
            e->value  = (double)(rng_next() % 100000) / 100.0;
            e->has_error = (rng_next() % 20 == 0); /* ~5% carry an error code */
            e->error_code = e->has_error ? (int32_t)(500 + rng_next() % 4) : 0;
            for (int b = 0; b < MS_SESSION_ID_LEN; b++)
                e->session_id[b] = (uint8_t)rng_next();
        }
        if (ms_writer_append(w, buf, m) != 0) {
            fprintf(stderr, "append failed\n");
            free(buf); ms_writer_abort(w); return 1;
        }
    }
    free(buf);

    if (ms_writer_close(w) != 0) { fprintf(stderr, "close failed\n"); return 1; }
    printf("Ingest complete.\n\n");

    /* ---- introspection ---- */
    printf("=== describe ===\n");
    if (ms_describe(path, stdout, errbuf, sizeof(errbuf)) != 0)
        fprintf(stderr, "describe failed: %s\n", errbuf);

    printf("\n=== validate ===\n");
    printf("  %s\n", ms_validate(path, errbuf, sizeof(errbuf)) == 0 ? "OK" : errbuf);

    /* ---- queries ---- */
    printf("\n=== query: metric=cpu.util, all time ===\n");
    ms_query_t q1 = { "cpu.util", NULL, MS_TS_MIN, MS_TS_MAX, false };
    ms_query_result_t r1;
    if (ms_query(path, &q1, NULL, NULL, &r1, errbuf, sizeof(errbuf)) == 0) {
        printf("  matched=%lld scanned=%lld sum=%.1f min=%.2f max=%.2f\n",
               (long long)r1.rows_matched, (long long)r1.rows_scanned,
               r1.value_sum, r1.value_min, r1.value_max);
    } else fprintf(stderr, "query failed: %s\n", errbuf);

    printf("\n=== query: host=db-01, first 10%% of the time window (pushdown) ===\n");
    int64_t window_hi = base_ts + (num_events / 10) * 1000LL;
    ms_query_t q2 = { NULL, "db-01", MS_TS_MIN, window_hi, false };
    ms_query_result_t r2;
    int shown = 0;
    if (ms_query(path, &q2, row_printer, &shown, &r2, errbuf, sizeof(errbuf)) == 0) {
        printf("  matched=%lld scanned=%lld row_groups=%d pruned=%d\n",
               (long long)r2.rows_matched, (long long)r2.rows_scanned,
               r2.row_groups_total, r2.row_groups_pruned);
    } else fprintf(stderr, "query failed: %s\n", errbuf);

    /* ---- bloom membership ---- */
    printf("\n=== bloom membership (host column) ===\n");
    const char* probes[] = { "db-01", "web-02", "nonexistent-host", "ghost-99" };
    for (int i = 0; i < NELEM(probes); i++) {
        int m = ms_might_contain_host(path, probes[i], errbuf, sizeof(errbuf));
        printf("  %-18s -> %s\n", probes[i],
               m < 0 ? errbuf : (m ? "maybe present" : "definitely absent"));
    }

    /* ---- parallel scan ---- */
    printf("\n=== parallel scan (worker pool) ===\n");
    ms_query_t qp = { "mem.used", NULL, MS_TS_MIN, MS_TS_MAX, false, true };
    ms_query_result_t rp;
    if (ms_query(path, &qp, NULL, NULL, &rp, errbuf, sizeof(errbuf)) == 0)
        printf("  metric=mem.used matched=%lld sum=%.1f\n",
               (long long)rp.rows_matched, rp.value_sum);

    /* ---- nested list<string> round-trip ---- */
    printf("\n=== nested LIST<STRING> self-check ===\n");
    printf("  %s\n", ms_nested_selfcheck("/tmp/metricstore_tags.parquet",
                                         errbuf, sizeof(errbuf)) == 0 ? "OK" : errbuf);

    /* ---- Arrow C Data Interface round-trip ---- */
    printf("\n=== Arrow C Data round-trip ===\n");
    int64_t moved = ms_arrow_roundtrip(path, "/tmp/metricstore_arrow.parquet",
                                       errbuf, sizeof(errbuf));
    if (moved > 0) printf("  bridged %lld rows through Arrow -> /tmp/metricstore_arrow.parquet\n",
                          (long long)moved);
    else fprintf(stderr, "  arrow roundtrip failed: %s\n", errbuf);

    printf("\nDone. (%s left on disk)\n", path);
    ms_shutdown();
    return 0;
}
