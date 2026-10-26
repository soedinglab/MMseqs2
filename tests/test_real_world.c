/**
 * @file test_real_world.c
 * @brief End-to-end "real dataset" tests that mirror how the library is used.
 *
 * Rather than probing one function in isolation, each test here builds a
 * realistic multi-column table — mixed physical types, nullable columns,
 * several row groups, compression, statistics — writes it through the public
 * writer, and reads it back through the public reader in the ways a real
 * application would: column-by-column, via the streaming batch reader, and
 * through the statistics / predicate-pushdown / key-value-metadata surfaces.
 *
 * The intent is a regression net for feature work: if a new encoding, codec,
 * or reader optimization silently corrupts values, drops nulls, miscounts row
 * groups, or breaks statistics, one of these full-pipeline assertions fails
 * with a concrete, reproducible message instead of surfacing later as a
 * mysterious data bug in production.
 */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

/* ---- an "events" dataset spread over several row groups ---------------- */

enum { RG = 4, PER_RG = 2500, TOTAL = RG * PER_RG };

/* Column layout:
 *   0 id       INT64    REQUIRED   monotonic primary key
 *   1 user_id  INT32    REQUIRED   low-cardinality (dictionary-friendly)
 *   2 score    DOUBLE   OPTIONAL   ~1/4 nulls
 *   3 name     BYTE_ARR OPTIONAL   ~1/3 nulls, shared prefixes
 *   4 active   BOOLEAN  REQUIRED
 *   5 price    FLOAT    REQUIRED
 */

typedef struct {
    int64_t id[TOTAL];
    int32_t user_id[TOTAL];
    double  score_present[TOTAL];  int16_t score_def[TOTAL]; int score_n;
    char    name_buf[TOTAL][24];
    carquet_byte_array_t name_present[TOTAL]; int16_t name_def[TOTAL]; int name_n;
    unsigned char active[TOTAL];
    float   price[TOTAL];
    /* expected reference for null-aware comparison */
    int     score_is_null[TOTAL];
    double  score_full[TOTAL];
    int     name_is_null[TOTAL];
} dataset_t;

static void build_dataset(dataset_t* d) {
    d->score_n = 0;
    d->name_n = 0;
    for (int i = 0; i < TOTAL; i++) {
        d->id[i] = 1000000 + (int64_t)i;
        d->user_id[i] = 100 + (i % 37);
        d->active[i] = (unsigned char)((i % 5) != 0);
        d->price[i] = (float)((i % 1000) * 0.01 + 0.5);

        /* score: null when i%4==0 */
        if (i % 4 == 0) {
            d->score_is_null[i] = 1;
            d->score_def[i] = 0;
        } else {
            double v = (double)i * 1.5 - 3000.0;
            d->score_is_null[i] = 0;
            d->score_full[i] = v;
            d->score_def[i] = 1;
            d->score_present[d->score_n++] = v;
        }

        /* name: null when i%3==0 */
        if (i % 3 == 0) {
            d->name_is_null[i] = 1;
            d->name_def[i] = 0;
        } else {
            d->name_is_null[i] = 0;
            d->name_def[i] = 1;
            int len = snprintf(d->name_buf[i], sizeof(d->name_buf[i]), "user/evt/%d", i);
            d->name_present[d->name_n].data = (uint8_t*)d->name_buf[i];
            d->name_present[d->name_n].length = len;
            d->name_n++;
        }
    }
}

static carquet_schema_t* make_schema(carquet_error_t* err) {
    carquet_schema_t* s = carquet_schema_create(err);
    if (!s) return NULL;
    int ok =
        carquet_schema_add_column(s, "id",      CARQUET_PHYSICAL_INT64,   NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "user_id", CARQUET_PHYSICAL_INT32,   NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "score",   CARQUET_PHYSICAL_DOUBLE,  NULL, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "name",    CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "active",  CARQUET_PHYSICAL_BOOLEAN, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK &&
        carquet_schema_add_column(s, "price",   CARQUET_PHYSICAL_FLOAT,   NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK;
    if (!ok) { carquet_schema_free(s); return NULL; }
    return s;
}

/* Write the dataset, one row group at a time, with the given compression. */
static int write_dataset(const char* path, const dataset_t* d,
                         carquet_compression_t comp) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = make_schema(&err);
    if (!s) return -1;

    carquet_writer_options_t wo;
    carquet_writer_options_init(&wo);
    wo.compression = comp;
    wo.write_statistics = true;
    carquet_writer_t* w = carquet_writer_create(path, s, &wo, &err);
    if (!w) { carquet_schema_free(s); return -1; }

    (void)carquet_writer_add_metadata(w, "producer", "carquet-test");
    (void)carquet_writer_add_metadata(w, "rows", "10000");

    int rc = -1;
    /* Track the sparse offsets into the present-value arrays per row group. */
    int score_off = 0, name_off = 0;
    for (int g = 0; g < RG; g++) {
        int base = g * PER_RG;
        /* Count present values in this row group's slice. */
        int score_here = 0, name_here = 0;
        for (int i = 0; i < PER_RG; i++) {
            if (!d->score_is_null[base + i]) score_here++;
            if (!d->name_is_null[base + i]) name_here++;
        }
        if (carquet_writer_write_batch(w, 0, d->id + base, PER_RG, NULL, NULL) != CARQUET_OK) goto done;
        if (carquet_writer_write_batch(w, 1, d->user_id + base, PER_RG, NULL, NULL) != CARQUET_OK) goto done;
        if (carquet_writer_write_batch(w, 2, d->score_present + score_off, PER_RG,
                                       d->score_def + base, NULL) != CARQUET_OK) goto done;
        if (carquet_writer_write_batch(w, 3, d->name_present + name_off, PER_RG,
                                       d->name_def + base, NULL) != CARQUET_OK) goto done;
        if (carquet_writer_write_batch(w, 4, d->active + base, PER_RG, NULL, NULL) != CARQUET_OK) goto done;
        if (carquet_writer_write_batch(w, 5, d->price + base, PER_RG, NULL, NULL) != CARQUET_OK) goto done;
        score_off += score_here;
        name_off += name_here;
        if (g + 1 < RG) {
            if (carquet_writer_new_row_group(w) != CARQUET_OK) goto done;
        }
    }
    if (carquet_writer_close(w) != CARQUET_OK) { carquet_schema_free(s); return -1; }
    carquet_schema_free(s);
    return 0;
done:
    carquet_writer_abort(w);
    carquet_schema_free(s);
    return rc;
}

/* ---- verification passes ----------------------------------------------- */

/* carquet returns nullable columns with DENSE (present-only) value buffers and
 * a per-row definition-level array — the same contract test_encoding_roundtrip
 * relies on. So value[k] is the k-th present value in the group, while def[i]
 * is per logical row: def==0 means row i is null. */
static int test_roundtrip_all_columns(const dataset_t* d, const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("roundtrip_all_columns", "reopen failed");
    if (carquet_reader_num_rows(r) != TOTAL) {
        carquet_reader_close(r);
        TEST_FAIL("roundtrip_all_columns", "num_rows mismatch");
    }

    int fail = 0;
    static int64_t got_id[PER_RG];
    static int32_t got_uid[PER_RG];
    static double  got_score[PER_RG]; static int16_t got_score_def[PER_RG];
    static carquet_byte_array_t got_name[PER_RG]; static int16_t got_name_def[PER_RG];
    static unsigned char got_active[PER_RG];
    static float   got_price[PER_RG];

    int32_t ngroups = TOTAL / PER_RG;
    for (int32_t g = 0; g < ngroups && !fail; g++) {
        int base = g * PER_RG;
        carquet_column_reader_t* c;
        /* BYTE_ARRAY values are zero-copy views into the column reader's
         * retained pages, so the name column's reader must outlive the
         * verification loop that reads got_name[].data. */
        carquet_column_reader_t* c_name = NULL;

        c = carquet_reader_get_column(r, g, 0, &err);
        if (!c || carquet_column_read_batch(c, got_id, PER_RG, NULL, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); break; }
        carquet_column_reader_free(c);
        c = carquet_reader_get_column(r, g, 1, &err);
        if (!c || carquet_column_read_batch(c, got_uid, PER_RG, NULL, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); break; }
        carquet_column_reader_free(c);
        c = carquet_reader_get_column(r, g, 2, &err);
        if (!c || carquet_column_read_batch(c, got_score, PER_RG, got_score_def, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); break; }
        carquet_column_reader_free(c);
        c_name = carquet_reader_get_column(r, g, 3, &err);
        if (!c_name || carquet_column_read_batch(c_name, got_name, PER_RG, got_name_def, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c_name); break; }
        /* Keep c_name alive: got_name[].data points into its retained pages. */
        c = carquet_reader_get_column(r, g, 4, &err);
        if (!c || carquet_column_read_batch(c, got_active, PER_RG, NULL, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); carquet_column_reader_free(c_name); break; }
        carquet_column_reader_free(c);
        c = carquet_reader_get_column(r, g, 5, &err);
        if (!c || carquet_column_read_batch(c, got_price, PER_RG, NULL, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); carquet_column_reader_free(c_name); break; }
        carquet_column_reader_free(c);

        int sk = 0, nk = 0;  /* dense present indices for score / name */
        for (int i = 0; i < PER_RG && !fail; i++) {
            if (got_id[i] != d->id[base + i]) fail = 1;
            else if (got_uid[i] != d->user_id[base + i]) fail = 1;
            else if (got_active[i] != d->active[base + i]) fail = 1;
            else if (got_price[i] != d->price[base + i]) fail = 1;

            int score_null = (got_score_def[i] == 0);
            if (score_null != d->score_is_null[base + i]) fail = 1;
            else if (!score_null) { if (got_score[sk++] != d->score_full[base + i]) fail = 1; }

            int name_null = (got_name_def[i] == 0);
            if (name_null != d->name_is_null[base + i]) fail = 1;
            else if (!name_null) {
                const char* exp = d->name_buf[base + i];
                if (got_name[nk].length != (int32_t)strlen(exp) ||
                    memcmp(got_name[nk].data, exp, strlen(exp)) != 0) fail = 1;
                nk++;
            }
        }
        carquet_column_reader_free(c_name);
    }
    carquet_reader_close(r);
    if (fail) TEST_FAIL("roundtrip_all_columns", "value mismatch");
    TEST_PASS("roundtrip_all_columns");
    return 0;
}

static int test_present_names_sequence(const dataset_t* d, const char* path) {
    /* Concatenate the dense present-name buffers across all row groups and
     * match them against the source's global present-value stream — a distinct
     * check that catches cross-page/cross-row-group value-order desync. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("present_names_sequence", "reopen failed");
    int32_t ngroups = TOTAL / PER_RG;
    static carquet_byte_array_t out[PER_RG];
    static int16_t def[PER_RG];
    int present_seen = 0, fail = 0;
    for (int32_t g = 0; g < ngroups && !fail; g++) {
        carquet_column_reader_t* c = carquet_reader_get_column(r, g, 3, &err);
        if (!c || carquet_column_read_batch(c, out, PER_RG, def, NULL) != PER_RG) { fail = 1; carquet_column_reader_free(c); break; }
        int present_in_group = 0;
        for (int i = 0; i < PER_RG; i++) if (def[i] != 0) present_in_group++;
        for (int k = 0; k < present_in_group; k++) {
            const carquet_byte_array_t* exp = &d->name_present[present_seen];
            if (out[k].length != exp->length ||
                memcmp(out[k].data, exp->data, (size_t)exp->length) != 0) { fail = 1; break; }
            present_seen++;
        }
        carquet_column_reader_free(c);
    }
    carquet_reader_close(r);
    if (fail || present_seen != d->name_n) TEST_FAIL("present_names_sequence", "present name stream mismatch");
    TEST_PASS("present_names_sequence");
    return 0;
}

static int test_statistics(const dataset_t* d, const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("statistics", "reopen failed");

    int32_t ngroups = TOTAL / PER_RG;
    int64_t total_score_nulls = 0;
    int fail = 0;
    for (int32_t g = 0; g < ngroups; g++) {
        /* id column: min/max must bracket the row-group id range exactly. */
        carquet_column_statistics_t st;
        if (carquet_reader_column_statistics(r, g, 0, &st) != CARQUET_OK) { fail = 1; break; }
        if (st.has_min_max) {
            int64_t mn, mx;
            memcpy(&mn, st.min_value, sizeof(mn));
            memcpy(&mx, st.max_value, sizeof(mx));
            int64_t exp_mn = d->id[g * PER_RG];
            int64_t exp_mx = d->id[g * PER_RG + PER_RG - 1];
            if (mn != exp_mn || mx != exp_mx) { fail = 1; break; }
        }
        /* score column: accumulate reported null counts. */
        if (carquet_reader_column_statistics(r, g, 2, &st) == CARQUET_OK && st.has_null_count)
            total_score_nulls += st.null_count;
    }
    carquet_reader_close(r);
    if (fail) TEST_FAIL("statistics", "id min/max wrong");

    /* Expected total score nulls. */
    int64_t exp_nulls = 0;
    for (int i = 0; i < TOTAL; i++) exp_nulls += d->score_is_null[i];
    if (total_score_nulls != exp_nulls) TEST_FAIL("statistics", "score null count wrong");
    TEST_PASS("statistics");
    return 0;
}

static int test_predicate_pushdown(const dataset_t* d, const char* path) {
    (void)d;
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("predicate_pushdown", "reopen failed");
    int32_t ngroups = TOTAL / PER_RG;

    /* id >= a value inside the last row group: only groups whose max can
     * satisfy it should report might_match. Since ids are globally monotonic,
     * exactly the groups from the target onward must match. */
    int64_t threshold = 1000000 + (int64_t)(TOTAL - 1);  /* the very last id */
    int matches = 0, fail = 0;
    for (int32_t g = 0; g < ngroups; g++) {
        bool might = false;
        if (carquet_reader_row_group_matches(r, g, 0, CARQUET_COMPARE_GE,
                                             &threshold, sizeof(threshold), &might) != CARQUET_OK) { fail = 1; break; }
        if (might) matches++;
    }
    /* A value that is below every id must prune nothing... */
    int64_t low = 0;
    bool all_match = true;
    for (int32_t g = 0; g < ngroups && all_match; g++) {
        bool might = false;
        (void)carquet_reader_row_group_matches(r, g, 0, CARQUET_COMPARE_GE, &low, sizeof(low), &might);
        if (!might) all_match = false;
    }
    /* ...and a value above every id must prune the earliest groups. */
    carquet_reader_close(r);
    if (fail) TEST_FAIL("predicate_pushdown", "row_group_matches errored");
    if (matches < 1) TEST_FAIL("predicate_pushdown", "last-id GE matched no group");
    if (matches == ngroups && ngroups > 1) TEST_FAIL("predicate_pushdown", "GE last id failed to prune any group");
    if (!all_match) TEST_FAIL("predicate_pushdown", "GE below-min pruned a group it shouldn't");
    TEST_PASS("predicate_pushdown");
    return 0;
}

static int test_kv_metadata(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("kv_metadata", "reopen failed");
    const char* producer = carquet_reader_find_metadata(r, "producer");
    const char* rows = carquet_reader_find_metadata(r, "rows");
    const char* missing = carquet_reader_find_metadata(r, "nope");
    int ok = producer && strcmp(producer, "carquet-test") == 0 &&
             rows && strcmp(rows, "10000") == 0 &&
             missing == NULL;
    carquet_reader_close(r);
    if (!ok) TEST_FAIL("kv_metadata", "metadata mismatch");
    TEST_PASS("kv_metadata");
    return 0;
}

static int test_batch_reader_crosscheck(const dataset_t* d, const char* path) {
    /* Stream the whole file through the batch reader (projecting id + price)
     * and confirm it agrees with the source, exercising a different read path
     * than the per-column reader above. */
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) TEST_FAIL("batch_reader_crosscheck", "reopen failed");

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    const char* cols[] = {"id", "price"};
    cfg.column_names = cols;
    cfg.num_column_names = 2;

    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { carquet_reader_close(r); TEST_FAIL("batch_reader_crosscheck", "batch reader create failed"); }

    int64_t seen = 0;
    int fail = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        int64_t n = carquet_row_batch_num_rows(batch);
        const void* id_data = NULL; const uint8_t* nb = NULL; int64_t nv = 0;
        const void* price_data = NULL; int64_t pv = 0;
        if (carquet_row_batch_column(batch, 0, &id_data, &nb, &nv) != CARQUET_OK ||
            carquet_row_batch_column(batch, 1, &price_data, &nb, &pv) != CARQUET_OK) { fail = 1; carquet_row_batch_free(batch); batch = NULL; break; }
        const int64_t* ids = (const int64_t*)id_data;
        const float* prices = (const float*)price_data;
        for (int64_t i = 0; i < n; i++) {
            if (ids[i] != d->id[seen + i] || prices[i] != d->price[seen + i]) { fail = 1; break; }
        }
        seen += n;
        carquet_row_batch_free(batch);
        batch = NULL;
        if (fail) break;
    }
    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    if (fail) TEST_FAIL("batch_reader_crosscheck", "batch value mismatch");
    if (seen != TOTAL) TEST_FAIL("batch_reader_crosscheck", "batch row count mismatch");
    TEST_PASS("batch_reader_crosscheck");
    return 0;
}

/* ---- driver: run the full suite under several compression codecs ------- */

static dataset_t g_dataset;  /* large; keep off the stack */

static int run_suite(carquet_compression_t comp, const char* tag) {
    char path[512];
    char base[64];
    snprintf(base, sizeof(base), "realworld_%s", tag);
    carquet_test_temp_path(path, sizeof(path), base);

    if (write_dataset(path, &g_dataset, comp) != 0) {
        carquet_test_cleanup(path);
        printf("[FAIL] real_world[%s]: write failed\n", tag);
        return 1;
    }
    int failures = 0;
    failures += test_roundtrip_all_columns(&g_dataset, path);
    failures += test_present_names_sequence(&g_dataset, path);
    failures += test_statistics(&g_dataset, path);
    failures += test_predicate_pushdown(&g_dataset, path);
    failures += test_kv_metadata(path);
    failures += test_batch_reader_crosscheck(&g_dataset, path);
    carquet_test_cleanup(path);
    if (failures) printf("  (%d failures under %s)\n", failures, tag);
    return failures;
}

int main(void) {
    build_dataset(&g_dataset);
    int failures = 0;
    failures += run_suite(CARQUET_COMPRESSION_UNCOMPRESSED, "none");
    failures += run_suite(CARQUET_COMPRESSION_SNAPPY, "snappy");
    failures += run_suite(CARQUET_COMPRESSION_ZSTD, "zstd");
    if (failures) { printf("\n%d real-world test(s) FAILED\n", failures); return 1; }
    printf("\nAll real-world integration tests passed\n");
    return 0;
}
