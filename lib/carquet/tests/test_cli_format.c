/**
 * @file test_cli_format.c
 * @brief Tests for the CLI value/type formatting helpers (src/cli/commands.c)
 *
 * Covers the temporal logical types the CLI renders specially: DATE,
 * TIME(MILLIS/MICROS/NANOS) and TIMESTAMP(MILLIS/MICROS/NANOS), including
 * pre-epoch and out-of-range values.
 */

#include "cli/cli.h"
#include "test_helpers.h"
#include <stdint.h>

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

int main(void) {
    printf("=== CLI formatting tests ===\n");
    int rc = 0;
    rc |= test_time_values();
    rc |= test_timestamp_and_date_values();
    rc |= test_type_names();
    if (rc) {
        printf("=== %d failure(s) ===\n", failures);
        return 1;
    }
    printf("=== All CLI formatting tests passed ===\n");
    return 0;
}
