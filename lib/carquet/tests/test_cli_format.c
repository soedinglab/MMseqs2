/**
 * @file test_cli_format.c
 * @brief Tests for the CLI value/type formatting helpers (src/cli/commands.c)
 *
 * Covers the temporal logical types the CLI renders specially: DATE,
 * TIME(MILLIS/MICROS/NANOS) and TIMESTAMP(MILLIS/MICROS/NANOS), including
 * pre-epoch and out-of-range values, and the row printers on LIST columns
 * (one list per row, aligned with flat columns across row groups).
 */

#include "cli/cli.h"
#include "test_helpers.h"
#include <stdint.h>
#if defined(_WIN32)
#include <io.h>
#define cap_dup _dup
#define cap_dup2 _dup2
#define cap_fileno _fileno
#define cap_close _close
#else
#include <unistd.h>
#define cap_dup dup
#define cap_dup2 dup2
#define cap_fileno fileno
#define cap_close close
#endif

static int failures = 0;

static void check(const char* name, const char* got, const char* want) {
    if (strcmp(got, want) != 0) {
        printf("[FAIL] %s: got \"%s\", want \"%s\"\n", name, got, want);
        failures++;
    }
}

static carquet_logical_type_t time_lt(carquet_time_unit_t unit) {
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_TIME;
    lt.params.time.unit = unit;
    lt.params.time.is_adjusted_to_utc = true;
    return lt;
}

static carquet_logical_type_t timestamp_lt(carquet_time_unit_t unit) {
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_TIMESTAMP;
    lt.params.timestamp.unit = unit;
    lt.params.timestamp.is_adjusted_to_utc = true;
    return lt;
}

static const char* fmt_time32(carquet_time_unit_t unit, int32_t v, char* buf, size_t n) {
    carquet_logical_type_t lt = time_lt(unit);
    return cli_format_value(CARQUET_PHYSICAL_INT32, &v, 0, &lt, buf, n);
}

static const char* fmt_time64(carquet_time_unit_t unit, int64_t v, char* buf, size_t n) {
    carquet_logical_type_t lt = time_lt(unit);
    return cli_format_value(CARQUET_PHYSICAL_INT64, &v, 0, &lt, buf, n);
}

static int test_time_values(void) {
    char buf[256];

    /* TIME(MILLIS) is INT32-backed; fraction always at the unit's precision. */
    check("time_ms zero",     fmt_time32(CARQUET_TIME_UNIT_MILLIS, 0, buf, sizeof buf),
          "00:00:00.000");
    check("time_ms end_of_day", fmt_time32(CARQUET_TIME_UNIT_MILLIS, 86399999, buf, sizeof buf),
          "23:59:59.999");
    check("time_ms one_hour",  fmt_time32(CARQUET_TIME_UNIT_MILLIS, 3600001, buf, sizeof buf),
          "01:00:00.001");

    /* TIME(MICROS) / TIME(NANOS) are INT64-backed. */
    check("time_us zero",     fmt_time64(CARQUET_TIME_UNIT_MICROS, 0, buf, sizeof buf),
          "00:00:00.000000");
    check("time_us small",    fmt_time64(CARQUET_TIME_UNIT_MICROS, 86399999, buf, sizeof buf),
          "00:01:26.399999");
    check("time_us end_of_day",
          fmt_time64(CARQUET_TIME_UNIT_MICROS, 86399999999LL, buf, sizeof buf),
          "23:59:59.999999");
    check("time_ns zero",     fmt_time64(CARQUET_TIME_UNIT_NANOS, 0, buf, sizeof buf),
          "00:00:00.000000000");
    check("time_ns end_of_day",
          fmt_time64(CARQUET_TIME_UNIT_NANOS, 86399999999999LL, buf, sizeof buf),
          "23:59:59.999999999");

    /* Out-of-spec values (negative, past 24h) print readably rather than
     * wrapping or overflowing. INT64_MIN must not trip signed negation UB. */
    check("time_us negative",  fmt_time64(CARQUET_TIME_UNIT_MICROS, -1, buf, sizeof buf),
          "-00:00:00.000001");
    check("time_us past_24h",  fmt_time64(CARQUET_TIME_UNIT_MICROS, 90000000000LL, buf, sizeof buf),
          "25:00:00.000000");
    check("time_us int64_min",
          fmt_time64(CARQUET_TIME_UNIT_MICROS, INT64_MIN, buf, sizeof buf),
          "-2562047788:00:54.775808");
    check("time_ms int32_min", fmt_time32(CARQUET_TIME_UNIT_MILLIS, INT32_MIN, buf, sizeof buf),
          "-596:31:23.648");

    if (failures) return 1;
    TEST_PASS("cli TIME value formatting");
    return 0;
}

static int test_timestamp_and_date_values(void) {
    char buf[256];
    int before = failures;

    carquet_logical_type_t lt = timestamp_lt(CARQUET_TIME_UNIT_MICROS);
    int64_t us = -62135596799137728LL; /* 0001-01-01T00:00:00.862272 */
    check("ts_us year_one",
          cli_format_value(CARQUET_PHYSICAL_INT64, &us, 0, &lt, buf, sizeof buf),
          "0001-01-01T00:00:00.862272");

    lt = timestamp_lt(CARQUET_TIME_UNIT_MILLIS);
    int64_t ms = 253402300798991LL; /* 9999-12-31T23:59:58.991 */
    check("ts_ms year_9999",
          cli_format_value(CARQUET_PHYSICAL_INT64, &ms, 0, &lt, buf, sizeof buf),
          "9999-12-31T23:59:58.991");

    ms = 0;
    check("ts_ms epoch",
          cli_format_value(CARQUET_PHYSICAL_INT64, &ms, 0, &lt, buf, sizeof buf),
          "1970-01-01T00:00:00");

    /* Pre-epoch: the seconds/fraction split must floor, not truncate. */
    ms = -1;
    check("ts_ms pre_epoch",
          cli_format_value(CARQUET_PHYSICAL_INT64, &ms, 0, &lt, buf, sizeof buf),
          "1969-12-31T23:59:59.999");
    ms = -86400000LL;
    check("ts_ms one_day_before",
          cli_format_value(CARQUET_PHYSICAL_INT64, &ms, 0, &lt, buf, sizeof buf),
          "1969-12-31T00:00:00");

    lt = timestamp_lt(CARQUET_TIME_UNIT_NANOS);
    int64_t ns = 1234567890123456789LL;
    check("ts_ns far_future",
          cli_format_value(CARQUET_PHYSICAL_INT64, &ns, 0, &lt, buf, sizeof buf),
          "2009-02-13T23:31:30.123456789");

    carquet_logical_type_t dlt;
    memset(&dlt, 0, sizeof(dlt));
    dlt.id = CARQUET_LOGICAL_DATE;
    struct { int32_t days; const char* want; } dates[] = {
        { 0,       "1970-01-01" },
        { -1,      "1969-12-31" },
        { -719162, "0001-01-01" },  /* Parquet DATE minimum */
        { 2932896, "9999-12-31" },  /* Parquet DATE maximum */
        { 59,      "1970-03-01" },
        { 719468,  "3939-11-03" },  /* past the MSVC CRT's year-3000 ceiling */
        { 11016,   "2000-02-29" },  /* leap year divisible by 400 */
        { -25567,  "1900-01-01" },  /* year divisible by 100, not a leap year */
    };
    for (size_t i = 0; i < sizeof dates / sizeof dates[0]; i++) {
        check("date", cli_format_value(CARQUET_PHYSICAL_INT32, &dates[i].days, 0,
                                       &dlt, buf, sizeof buf),
              dates[i].want);
    }

    if (failures != before) return 1;
    TEST_PASS("cli DATE/TIMESTAMP value formatting");
    return 0;
}

static int test_type_names(void) {
    char buf[64];
    int before = failures;

    carquet_logical_type_t lt = time_lt(CARQUET_TIME_UNIT_MILLIS);
    cli_format_type(CARQUET_PHYSICAL_INT32, &lt, buf, sizeof buf);
    check("type TIME(ms)", buf, "TIME(ms,UTC)");

    lt = time_lt(CARQUET_TIME_UNIT_MICROS);
    lt.params.time.is_adjusted_to_utc = false;
    cli_format_type(CARQUET_PHYSICAL_INT64, &lt, buf, sizeof buf);
    check("type TIME(us) local", buf, "TIME(us)");

    if (failures != before) return 1;
    TEST_PASS("cli type naming");
    return 0;
}

/* `carquet stat` renders raw footer statistics: a malformed file can give an
 * INT64 / DOUBLE / INT96 / FLBA column a stat shorter than the type, which
 * used to be dereferenced as a full value (heap over-read). Exact-size heap
 * buffers make the old over-read visible under ASan. */
static int test_stat_value_width(void) {
    int before = failures;
    char buf[64];
    uint8_t* one = (uint8_t*)malloc(1);
    one[0] = 7;
    check("stat int64 1 byte",
          cli_format_stat_value(CARQUET_PHYSICAL_INT64, one, 1, 0, NULL, buf, sizeof(buf)),
          "(invalid)");
    check("stat double 1 byte",
          cli_format_stat_value(CARQUET_PHYSICAL_DOUBLE, one, 1, 0, NULL, buf, sizeof(buf)),
          "(invalid)");
    check("stat int96 1 byte",
          cli_format_stat_value(CARQUET_PHYSICAL_INT96, one, 1, 0, NULL, buf, sizeof(buf)),
          "(invalid)");
    check("stat flba(16) 1 byte",
          cli_format_stat_value(CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, one, 1, 16, NULL,
                                buf, sizeof(buf)),
          "(invalid)");
    check("stat negative size",
          cli_format_stat_value(CARQUET_PHYSICAL_INT32, one, -1, 0, NULL, buf, sizeof(buf)),
          "(invalid)");
    free(one);

    /* Correct widths still render, from deliberately misaligned storage. */
    uint8_t* raw = (uint8_t*)malloc(9);
    int64_t v = -42;
    memcpy(raw + 1, &v, 8);
    check("stat int64 misaligned",
          cli_format_stat_value(CARQUET_PHYSICAL_INT64, raw + 1, 8, 0, NULL, buf, sizeof(buf)),
          "-42");
    int32_t v32 = 123;
    memcpy(raw + 1, &v32, 4);
    check("stat int32",
          cli_format_stat_value(CARQUET_PHYSICAL_INT32, raw + 1, 4, 0, NULL, buf, sizeof(buf)),
          "123");
    raw[1] = 0xAB; raw[2] = 0xCD;
    check("stat flba(2)",
          cli_format_stat_value(CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, raw + 1, 2, 2, NULL,
                                buf, sizeof(buf)),
          "abcd");
    free(raw);

    if (failures == before) TEST_PASS("stat_value_width");
    return failures != before;
}

/* Regression (CLI sweep over parquet-testing, nonnullable.impala.parquet):
 * the row printers took "nullable" from the leaf's own repetition and the
 * null threshold from carquet_schema_node_max_def_level() (the node's 0/1
 * share), so a REQUIRED leaf under an OPTIONAL group read no def levels, the
 * dense value index ran past the decoded values and `sample` formatted a stale
 * BYTE_ARRAY pointer (bus error). Every printer must handle such a column. */
static int test_required_leaf_under_optional_group(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "cli_req_under_opt");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    int32_t g = carquet_schema_add_group(s, "g", CARQUET_REPETITION_OPTIONAL, 0);
    carquet_schema_add_column(s, "name", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, g);
    carquet_writer_t* w = carquet_writer_create(path, s, NULL, &err);
    enum { N = 64 };
    carquet_byte_array_t vals[N];
    int16_t def[N];
    int nv = 0;
    for (int i = 0; i < N; i++) {
        def[i] = (i % 2) ? 1 : 0;              /* every other group is null */
        if (def[i]) { vals[nv].data = (uint8_t*)"abc"; vals[nv].length = 3; nv++; }
    }
    int ok = w && carquet_writer_write_batch(w, 0, vals, N, def, NULL) == CARQUET_OK &&
             carquet_writer_close(w) == CARQUET_OK;
    carquet_schema_free(s);
    if (!ok) { printf("[FAIL] required_leaf_under_optional_group: write\n"); failures++; return 1; }

    /* Output goes to stdout; only the return codes / absence of a crash or
     * sanitizer report are checked. */
    fflush(stdout);
    int rc = cmd_head(path, 10, NULL) | cmd_tail(path, 10, NULL) |
             cmd_sample(path, N, NULL);
    carquet_test_cleanup(path);
    if (rc != 0) {
        printf("[FAIL] required_leaf_under_optional_group: printer returned %d\n", rc);
        failures++;
        return 1;
    }
    printf("[PASS] required_leaf_under_optional_group\n");
    return 0;
}

/* Run `fn` with stdout redirected to a temp file; return the output. */
static char* capture_stdout(int (*fn)(const char*, void*), const char* path, void* arg,
                            int* rc) {
    char out_path[512];
    carquet_test_temp_path_ext(out_path, sizeof(out_path), "cli_capture", "txt");
    fflush(stdout);
    int saved = cap_dup(cap_fileno(stdout));
    FILE* f = fopen(out_path, "w+b");
    if (saved < 0 || !f) { if (f) fclose(f); *rc = -1; return NULL; }
    cap_dup2(cap_fileno(f), cap_fileno(stdout));
    *rc = fn(path, arg);
    fflush(stdout);
    cap_dup2(saved, cap_fileno(stdout));
    cap_close(saved);
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    rewind(f);
    char* buf = malloc((size_t)n + 1);
    size_t got = buf ? fread(buf, 1, (size_t)n, f) : 0;
    if (buf) buf[got] = '\0';
    fclose(f);
    remove(out_path);
    return buf;
}

static int run_export(const char* path, void* arg) {
    return cmd_export(path, (const row_select_opts_t*)arg, CLI_EXPORT_CSV);
}
static int run_tail3(const char* path, void* arg) { (void)arg; return cmd_tail(path, 3, NULL); }
static int run_head_huge_filter(const char* path, void* arg) {
    (void)arg;
    return cmd_head(path, INT64_MAX / 2, "id >= 6");
}

/* Rows r = 0..9 in two row groups of 5: id = r, o = r * 10 (null when r % 3
 * == 0), l by r % 4: [r, r+100, r+200] / [] / null / [null, r]. */
static int write_list_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    int ok = s &&
        carquet_schema_add_column(s, "id", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_list(s, "l", CARQUET_PHYSICAL_INT32, NULL,
                                CARQUET_REPETITION_OPTIONAL, 0, 0) >= 0 &&
        carquet_schema_add_column(s, "o", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK;
    carquet_writer_options_t o;
    carquet_writer_options_init(&o);
    o.write_page_index = true;
    carquet_writer_t* w = ok ? carquet_writer_create(path, s, &o, &err) : NULL;
    ok = w != NULL;
    for (int g = 0; ok && g < 2; g++) {
        int32_t ids[5], lv[20], ov[5];
        int16_t ld[20], lr[20], od[5];
        int nl = 0, nv = 0, no = 0;
        for (int k = 0; k < 5; k++) {
            int r = g * 5 + k;
            ids[k] = r;
            if (r % 3 == 0) od[k] = 0; else { od[k] = 1; ov[no++] = r * 10; }
            switch (r % 4) {
                case 0:
                    for (int j = 0; j < 3; j++) {
                        lr[nl] = j ? 1 : 0; ld[nl] = 3; nl++; lv[nv++] = r + 100 * j;
                    }
                    break;
                case 1: lr[nl] = 0; ld[nl] = 1; nl++; break;            /* [] */
                case 2: lr[nl] = 0; ld[nl] = 0; nl++; break;            /* null */
                default:
                    lr[nl] = 0; ld[nl] = 2; nl++;                       /* null elem */
                    lr[nl] = 1; ld[nl] = 3; nl++; lv[nv++] = r;
                    break;
            }
        }
        ok = carquet_writer_write_batch(w, 0, ids, 5, NULL, NULL) == CARQUET_OK &&
             carquet_writer_write_batch(w, 1, lv, nl, ld, lr) == CARQUET_OK &&
             carquet_writer_write_batch(w, 2, ov, 5, od, NULL) == CARQUET_OK &&
             (g == 1 || carquet_writer_new_row_group(w) == CARQUET_OK);
    }
    if (w && carquet_writer_close(w) != CARQUET_OK) ok = 0;
    if (s) carquet_schema_free(s);
    return ok;
}

static int test_list_rows(void) {
    const char* name = "list_rows_aligned";
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "cli_list_rows");
    if (!write_list_file(path)) {
        printf("[FAIL] %s: write\n", name);
        failures++;
        return 1;
    }
    int bad = 0, rc = 0;

    row_select_opts_t all = {0, -1, NULL, NULL};
    char* out = capture_stdout(run_export, path, &all, &rc);
    check("export all", out ? out : "",
          "id,element,o\n"
          "0,\"[0, 100, 200]\",\n"
          "1,[],10\n"
          "2,,20\n"
          "3,\"[null, 3]\",\n"
          "4,\"[4, 104, 204]\",40\n"
          "5,[],50\n"
          "6,,\n"
          "7,\"[null, 7]\",70\n"
          "8,\"[8, 108, 208]\",80\n"
          "9,[],\n");
    bad |= rc;
    free(out);

    /* Slice crossing the row-group boundary (skip by rows, not slots). */
    row_select_opts_t slice = {3, 4, "id,element", NULL};
    out = capture_stdout(run_export, path, &slice, &rc);
    check("export slice", out ? out : "",
          "id,element\n3,\"[null, 3]\"\n4,\"[4, 104, 204]\"\n5,[]\n6,\n");
    bad |= rc;
    free(out);

    out = capture_stdout(run_tail3, path, NULL, &rc);
    bad |= rc;
    if (!out || !strstr(out, "7   [null, 7]      70") ||
        !strstr(out, "9   []             null") || strstr(out, "  6 ")) {
        printf("[FAIL] %s: tail output:\n%s\n", name, out ? out : "(none)");
        failures++;
        bad = 1;
    }
    free(out);

    /* head -n <huge> --filter: the row count used to size the output. */
    out = capture_stdout(run_head_huge_filter, path, NULL, &rc);
    bad |= rc;
    if (!out || !strstr(out, "8   [8, 108, 208]  80")) {
        printf("[FAIL] %s: head --filter output:\n%s\n", name, out ? out : "(none)");
        failures++;
        bad = 1;
    }
    free(out);

    carquet_test_cleanup(path);
    if (bad) {
        printf("[FAIL] %s\n", name);
        failures++;
        return 1;
    }
    printf("[PASS] %s\n", name);
    return 0;
}

int main(void) {
    printf("=== CLI formatting tests ===\n");
    int rc = 0;
    rc |= test_time_values();
    rc |= test_timestamp_and_date_values();
    rc |= test_type_names();
    rc |= test_stat_value_width();
    rc |= test_required_leaf_under_optional_group();
    rc |= test_list_rows();
    if (rc) {
        printf("=== %d failure(s) ===\n", failures);
        return 1;
    }
    printf("=== All CLI formatting tests passed ===\n");
    return 0;
}
