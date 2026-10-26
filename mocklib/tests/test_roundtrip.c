/**
 * @file test_roundtrip.c
 * @brief Self-checking round-trip test for MetricStore.
 *
 * Writes a small, fully-known dataset and asserts that queries, aggregates,
 * predicate pushdown and bloom-filter membership all return exact,
 * hand-computed answers. Exits non-zero on the first failure.
 */
#include "metricstore/metricstore.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("[FAIL] %s\n", (msg)); failures++; } \
    else         { printf("[ok]   %s\n", (msg)); } \
} while (0)

#define CHECK_EQ_LL(a, b, msg) do { \
    long long _a = (long long)(a), _b = (long long)(b); \
    if (_a != _b) { printf("[FAIL] %s (got %lld, want %lld)\n", (msg), _a, _b); failures++; } \
    else          { printf("[ok]   %s\n", (msg)); } \
} while (0)

static ms_event_t make_event(int64_t id, int64_t ts, const char* host,
                             const char* metric, double v, bool err) {
    ms_event_t e;
    memset(&e, 0, sizeof(e));
    e.event_id = id;
    e.ts_micros = ts;
    e.host = host;
    e.region = "test-region";
    e.metric = metric;
    e.value = v;
    e.has_error = err;
    e.error_code = err ? 503 : 0;
    for (int b = 0; b < MS_SESSION_ID_LEN; b++) e.session_id[b] = (uint8_t)(id + b);
    return e;
}

int main(void) {
    const char* path = "/tmp/metricstore_roundtrip.parquet";
    char errbuf[512] = {0};

    printf("MetricStore round-trip test (backend %s)\n\n", ms_backend_version());

    /* Known dataset: 30 rows, ts = 1000..1029, 3 hosts round-robin,
       metric alternating cpu/mem, value = id. rows_per_group=10 => 3 groups. */
    ms_writer_config_t cfg;
    ms_writer_config_init(&cfg);
    cfg.rows_per_group = 10;

    metricstore_writer_t* w = ms_writer_open(path, &cfg, errbuf, sizeof(errbuf));
    CHECK(w != NULL, "writer opens");
    if (!w) return 1;

    const char* hosts[] = { "alpha", "beta", "gamma" };
    ms_event_t rows[30];
    for (int i = 0; i < 30; i++) {
        rows[i] = make_event(i, 1000 + i, hosts[i % 3],
                             (i % 2 == 0) ? "cpu" : "mem",
                             (double)i, (i % 5 == 0));
    }
    CHECK(ms_writer_append(w, rows, 30) == 0, "append 30 rows");
    CHECK(ms_writer_close(w) == 0, "writer closes");

    CHECK(ms_validate(path, errbuf, sizeof(errbuf)) == 0, "file validates");

    /* Query all cpu rows (ids 0,2,4,...,28 => 15 rows, values 0..28 even). */
    ms_query_t q_cpu = { "cpu", NULL, MS_TS_MIN, MS_TS_MAX, false };
    ms_query_result_t r;
    CHECK(ms_query(path, &q_cpu, NULL, NULL, &r, errbuf, sizeof(errbuf)) == 0, "query cpu ok");
    CHECK_EQ_LL(r.rows_matched, 15, "cpu matched == 15");
    /* sum of even numbers 0..28 = 2*(0+1+..+14) = 2*105 = 210 */
    CHECK(r.value_sum == 210.0, "cpu value_sum == 210");
    CHECK(r.value_min == 0.0 && r.value_max == 28.0, "cpu min/max == 0/28");
    CHECK_EQ_LL(r.row_groups_total, 3, "file has 3 row groups");

    /* Query host=alpha (ids 0,3,6,...,27 => 10 rows). */
    ms_query_t q_alpha = { NULL, "alpha", MS_TS_MIN, MS_TS_MAX, false };
    ms_query_result_t ra;
    CHECK(ms_query(path, &q_alpha, NULL, NULL, &ra, errbuf, sizeof(errbuf)) == 0, "query alpha ok");
    CHECK_EQ_LL(ra.rows_matched, 10, "alpha matched == 10");

    /* Bounded ts range hitting only the last group (ts 1020..1029 => ids 20..29).
       With 3 groups of 10 and monotonic ts, the first two groups must prune. */
    ms_query_t q_range = { NULL, NULL, 1020, 2000, false };
    ms_query_result_t rr;
    CHECK(ms_query(path, &q_range, NULL, NULL, &rr, errbuf, sizeof(errbuf)) == 0, "range query ok");
    CHECK_EQ_LL(rr.rows_matched, 10, "range matched == 10");
    CHECK_EQ_LL(rr.rows_scanned, 10, "range scanned == 10 (pushdown worked)");
    CHECK_EQ_LL(rr.row_groups_pruned, 2, "range pruned 2 row groups");

    /* Bloom membership: present hosts maybe-present, bogus host absent. */
    CHECK(ms_might_contain_host(path, "alpha", errbuf, sizeof(errbuf)) == 1, "bloom: alpha maybe present");
    CHECK(ms_might_contain_host(path, "gamma", errbuf, sizeof(errbuf)) == 1, "bloom: gamma maybe present");
    CHECK(ms_might_contain_host(path, "zzz-not-here", errbuf, sizeof(errbuf)) == 0, "bloom: bogus host absent");

    /* Typed (INT64) bloom on event_id. */
    CHECK(ms_might_contain_event_id(path, 7, errbuf, sizeof(errbuf)) == 1, "bloom: event_id 7 maybe present");
    CHECK(ms_might_contain_event_id(path, 999999, errbuf, sizeof(errbuf)) == 0, "bloom: event_id 999999 absent");

    /* Parallel query returns identical aggregates to the serial path. */
    ms_query_t q_par = { "cpu", NULL, MS_TS_MIN, MS_TS_MAX, false, true };
    ms_query_result_t rp;
    CHECK(ms_query(path, &q_par, NULL, NULL, &rp, errbuf, sizeof(errbuf)) == 0, "parallel query ok");
    CHECK_EQ_LL(rp.rows_matched, 15, "parallel cpu matched == 15");
    CHECK(rp.value_sum == 210.0, "parallel cpu value_sum == 210");

    /* Low-level column reader: read event ids with a skip. */
    int64_t ids_out[8];
    int64_t got = ms_read_event_ids(path, /*row_group=*/0, /*skip=*/3, ids_out, 8, errbuf, sizeof(errbuf));
    CHECK(got >= 1, "low-level read returned ids");
    CHECK(got >= 1 && ids_out[0] == 3, "low-level skip landed on id 3");

    /* In-memory pack + query on the blob. */
    void* blob = NULL; size_t blob_sz = 0;
    CHECK(ms_pack(rows, 30, NULL, &blob, &blob_sz, errbuf, sizeof(errbuf)) == 0, "pack to buffer ok");
    CHECK(blob != NULL && blob_sz > 0, "buffer non-empty");
    ms_query_t q_all = { "cpu", NULL, MS_TS_MIN, MS_TS_MAX, false, false };
    ms_query_result_t rb;
    CHECK(ms_query_buffer(blob, blob_sz, &q_all, &rb, errbuf, sizeof(errbuf)) == 0, "query buffer ok");
    CHECK_EQ_LL(rb.rows_matched, 15, "buffer cpu matched == 15");
    ms_free_buffer(blob);

    /* Append a fourth row group and re-query the whole store. */
    ms_event_t extra[5];
    for (int i = 0; i < 5; i++)
        extra[i] = make_event(100 + i, 5000 + i, "delta", "cpu", 1.0, false);
    CHECK(ms_append(path, extra, 5, errbuf, sizeof(errbuf)) == 0, "append row group ok");
    ms_query_t q_delta = { NULL, "delta", MS_TS_MIN, MS_TS_MAX, false, false };
    ms_query_result_t rd;
    CHECK(ms_query(path, &q_delta, NULL, NULL, &rd, errbuf, sizeof(errbuf)) == 0, "query appended ok");
    CHECK_EQ_LL(rd.rows_matched, 5, "appended delta matched == 5");
    CHECK_EQ_LL(rd.row_groups_total, 4, "store now has 4 row groups");

    /* Nested LIST<STRING> round-trip self-check. */
    CHECK(ms_nested_selfcheck("/tmp/metricstore_nested.parquet", errbuf, sizeof(errbuf)) == 0, "nested list round-trip");

    /* Arrow C Data Interface round-trip. */
    int64_t moved = ms_arrow_roundtrip(path, "/tmp/metricstore_arrow.parquet", errbuf, sizeof(errbuf));
    CHECK(moved > 0, "arrow C-data round-trip moved rows");

    printf("\n%s (%d failure%s)\n",
           failures == 0 ? "ALL PASSED" : "FAILED",
           failures, failures == 1 ? "" : "s");
    ms_shutdown();
    return failures == 0 ? 0 : 1;
}
