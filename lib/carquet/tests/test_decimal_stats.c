/**
 * @file test_decimal_stats.c
 * @brief Regression tests for DECIMAL min/max statistics ordering.
 *
 * A DECIMAL column backed by FIXED_LEN_BYTE_ARRAY or BYTE_ARRAY is ordered as a
 * SIGNED big-endian two's-complement integer, not unsigned-lexicographic. The
 * writer previously compared these stats with an unsigned memcmp, so a negative
 * value (leading 0xFF...) recorded as the maximum and inverted the recorded
 * min/max; the reader's predicate pushdown made the same mistake and could
 * prune row groups that actually matched. These tests pin the signed ordering
 * on both the write (recorded stats) and read (row-group pruning) paths, for
 * both physical backings.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"
#include <string.h>

/* DECIMAL logical type helper (precision/scale chosen to fit the backing). */
static carquet_logical_type_t decimal_lt(int32_t precision, int32_t scale) {
    carquet_logical_type_t lt;
    memset(&lt, 0, sizeof(lt));
    lt.id = CARQUET_LOGICAL_DECIMAL;
    lt.params.decimal.precision = precision;
    lt.params.decimal.scale = scale;
    return lt;
}

/* ---- FLBA DECIMAL: mixed signs, recorded min/max must be signed ----------- */
static int test_flba_decimal_stats(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "dec_flba_stats");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t lt = decimal_lt(9, 0);
    if (carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
                                  &lt, CARQUET_REPETITION_REQUIRED, 4, 0) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL("flba_decimal_stats", "add_column");
    }

    /* 4-byte big-endian two's complement: 1, -1, 100, -100.
     * Signed order => min = -100 (FF FF FF 9C), max = 100 (00 00 00 64).
     * Unsigned lex (the bug) would give min = 1, max = -1. */
    uint8_t vals[4 * 4] = {
        0x00, 0x00, 0x00, 0x01,  /* 1    */
        0xFF, 0xFF, 0xFF, 0xFF,  /* -1   */
        0x00, 0x00, 0x00, 0x64,  /* 100  */
        0xFF, 0xFF, 0xFF, 0x9C,  /* -100 */
    };
    const uint8_t exp_min[4] = {0xFF, 0xFF, 0xFF, 0x9C};  /* -100 */
    const uint8_t exp_max[4] = {0x00, 0x00, 0x00, 0x64};  /*  100 */

    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.write_statistics = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("flba_decimal_stats", "create"); }
    /* Two batches so page stats are merged at the column-chunk level too. */
    if (carquet_writer_write_batch(w, 0, vals, 2, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 0, vals + 8, 2, NULL, NULL) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL("flba_decimal_stats", "write");
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); TEST_FAIL("flba_decimal_stats", "close"); }
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("flba_decimal_stats", "open"); }

    carquet_column_statistics_t st;
    if (carquet_reader_column_statistics(r, 0, 0, &st) != CARQUET_OK || !st.has_min_max) {
        carquet_reader_close(r); carquet_test_cleanup(path);
        TEST_FAIL("flba_decimal_stats", "no stats");
    }
    int ok = st.min_value_size == 4 && st.max_value_size == 4 &&
             memcmp(st.min_value, exp_min, 4) == 0 &&
             memcmp(st.max_value, exp_max, 4) == 0;
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("flba_decimal_stats", "min/max not signed-ordered");
    TEST_PASS("flba_decimal_stats");
    return 0;
}

/* ---- BYTE_ARRAY DECIMAL: variable length, signed min/max ------------------ */
static int test_byte_array_decimal_stats(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "dec_ba_stats");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t lt = decimal_lt(18, 0);
    if (carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_BYTE_ARRAY,
                                  &lt, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL("ba_decimal_stats", "add_column");
    }

    /* Minimal big-endian two's complement of 1, -1, 256, -256.
     * Signed => min = -256 (FF 00), max = 256 (01 00). */
    static const uint8_t b_one[]     = {0x01};
    static const uint8_t b_negone[]  = {0xFF};
    static const uint8_t b_256[]     = {0x01, 0x00};
    static const uint8_t b_neg256[]  = {0xFF, 0x00};
    carquet_byte_array_t vals[4] = {
        { (uint8_t*)b_one,    1 },
        { (uint8_t*)b_negone, 1 },
        { (uint8_t*)b_256,    2 },
        { (uint8_t*)b_neg256, 2 },
    };

    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.write_statistics = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); TEST_FAIL("ba_decimal_stats", "create"); }
    if (carquet_writer_write_batch(w, 0, vals, 2, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 0, vals + 2, 2, NULL, NULL) != CARQUET_OK) {
        carquet_schema_free(s); TEST_FAIL("ba_decimal_stats", "write");
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); TEST_FAIL("ba_decimal_stats", "close"); }
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("ba_decimal_stats", "open"); }

    carquet_column_statistics_t st;
    if (carquet_reader_column_statistics(r, 0, 0, &st) != CARQUET_OK || !st.has_min_max) {
        carquet_reader_close(r); carquet_test_cleanup(path);
        TEST_FAIL("ba_decimal_stats", "no stats");
    }
    int ok = st.min_value_size == 2 && st.max_value_size == 2 &&
             memcmp(st.min_value, b_neg256, 2) == 0 &&
             memcmp(st.max_value, b_256, 2) == 0;
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!ok) TEST_FAIL("ba_decimal_stats", "min/max not signed-ordered");
    TEST_PASS("ba_decimal_stats");
    return 0;
}

/* ---- Reader pruning: signed compare in row-group predicate pushdown -------- */
static int test_decimal_row_group_pruning(void) {
    char path[512];
    carquet_test_temp_path(path, sizeof(path), "dec_prune");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_logical_type_t lt = decimal_lt(9, 0);
    carquet_schema_add_column(s, "d", CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY,
                              &lt, CARQUET_REPETITION_REQUIRED, 4, 0);

    /* All-negative values: min = -100, max = -1. Both signed and unsigned
     * writers agree here (negatives are lex-ordered the same as signed), so the
     * recorded stats are identical — this isolates the READER comparison. */
    uint8_t vals[2 * 4] = {
        0xFF, 0xFF, 0xFF, 0xFF,  /* -1   (max) */
        0xFF, 0xFF, 0xFF, 0x9C,  /* -100 (min) */
    };
    carquet_writer_options_t wo; carquet_writer_options_init(&wo);
    wo.write_statistics = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    carquet_writer_write_batch(w, 0, vals, 2, NULL, NULL);
    carquet_writer_close(w);
    carquet_schema_free(s);

    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { carquet_test_cleanup(path); TEST_FAIL("decimal_pruning", "open"); }

    /* Predicate: d >= 0. Every value is negative (max = -1 < 0), so the row
     * group cannot match and must be pruned. A lexicographic compare would rank
     * the predicate 0x00000000 below the stored max 0xFFFFFFFF and (wrongly)
     * keep the group. */
    const uint8_t zero[4] = {0x00, 0x00, 0x00, 0x00};
    bool might_match = true;
    carquet_status_t cs = carquet_reader_row_group_matches(
        r, 0, 0, CARQUET_COMPARE_GE, zero, 4, &might_match);
    int prune_ok = (cs == CARQUET_OK) && (might_match == false);

    /* Sanity: d >= -50 straddles the range, so it must NOT be pruned. */
    const uint8_t neg50[4] = {0xFF, 0xFF, 0xFF, 0xCE};  /* -50 */
    bool mm2 = false;
    carquet_status_t cs2 = carquet_reader_row_group_matches(
        r, 0, 0, CARQUET_COMPARE_GE, neg50, 4, &mm2);
    int keep_ok = (cs2 == CARQUET_OK) && (mm2 == true);

    carquet_reader_close(r);
    carquet_test_cleanup(path);
    if (!prune_ok) TEST_FAIL("decimal_pruning", "GE 0 not pruned (lexicographic compare?)");
    if (!keep_ok)  TEST_FAIL("decimal_pruning", "GE -50 wrongly pruned");
    TEST_PASS("decimal_pruning");
    return 0;
}

int main(void) {
    int rc = 0;
    rc |= test_flba_decimal_stats();
    rc |= test_byte_array_decimal_stats();
    rc |= test_decimal_row_group_pruning();
    if (rc == 0) printf("\nAll decimal-stats tests passed.\n");
    return rc;
}
