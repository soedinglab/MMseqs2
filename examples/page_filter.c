/**
 * @file page_filter.c
 * @brief Page-level predicate pushdown with carquet_batch_reader_set_page_filter()
 *
 * Writes a file with a page index (required for page filtering), then attaches
 * a conjunctive filter so only pages whose column-index min/max range could
 * match are decompressed and decoded. carquet_batch_reader_rows_skipped()
 * reports how many rows were pruned without being materialized.
 *
 * Page filters are conservative: rows inside a matching page that do not
 * satisfy the predicate are still returned, so we re-check exactly here.
 */

#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOTAL_ROWS    10000
#define ROWS_PER_PAGE 500
#define CHECK(expr) do { if ((expr) != CARQUET_OK) { \
    fprintf(stderr, "FAIL at %s:%d\n", __FILE__, __LINE__); return 1; } } while(0)

static int write_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    carquet_schema_add_column(schema, "id",     CARQUET_PHYSICAL_INT64,      NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "region", CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression      = CARQUET_COMPRESSION_ZSTD;
    opts.write_statistics = true;
    opts.write_page_index = true;            /* required for page filtering */
    opts.max_rows_per_page = ROWS_PER_PAGE;  /* small pages => fine-grained pruning */
    opts.write_batch_size  = ROWS_PER_PAGE;  /* chunk writes so each flush honors max_rows_per_page */
    /* Page splitting by row count happens on the eager per-batch flush; the
     * dictionary path instead flushes a single page on close, so force PLAIN. */
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { fprintf(stderr, "create: %s\n", err.message); carquet_schema_free(schema); return 1; }
    CHECK(carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_PLAIN));
    CHECK(carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_PLAIN));

    /* id is monotonic so each page covers a tight [min,max] range. region
     * alternates "east"/"west" per page so a region filter prunes half. */
    int64_t* ids = malloc(TOTAL_ROWS * sizeof(int64_t));
    carquet_byte_array_t* regions = malloc(TOTAL_ROWS * sizeof(carquet_byte_array_t));
    for (int i = 0; i < TOTAL_ROWS; i++) {
        ids[i] = i;
        int page = i / ROWS_PER_PAGE;
        const char* reg = (page % 2 == 0) ? "east" : "west";
        regions[i].data   = (uint8_t*)reg;
        regions[i].length = 4;
    }

    CHECK(carquet_writer_write_batch(w, 0, ids,     TOTAL_ROWS, NULL, NULL));
    CHECK(carquet_writer_write_batch(w, 1, regions, TOTAL_ROWS, NULL, NULL));
    CHECK(carquet_writer_close(w));

    free(ids); free(regions);
    carquet_schema_free(schema);
    return 0;
}

/* Run a filter and report rows returned, rows matching exactly, and pruned. */
static int run_filter(const char* path, const char* label,
                      const carquet_filter_clause_t* clauses, int32_t count,
                      int (*matches)(int64_t id, const char* region, int rlen)) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);
    if (!r) { fprintf(stderr, "open: %s\n", err.message); return 1; }

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);

    carquet_status_t st = carquet_batch_reader_set_page_filter(br, clauses, count);
    if (st != CARQUET_OK) {
        fprintf(stderr, "  set_page_filter: %s\n", carquet_status_string(st));
        carquet_batch_reader_free(br); carquet_reader_close(r); return 1;
    }

    int64_t returned = 0, exact = 0;
    carquet_row_batch_t* batch = NULL;
    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* iddata; const uint8_t* idnull; int64_t n;
        carquet_row_batch_column(batch, 0, &iddata, &idnull, &n);
        const int64_t* ids = (const int64_t*)iddata;

        const void* regdata; const uint8_t* regnull; int64_t rn;
        carquet_row_batch_column(batch, 1, &regdata, &regnull, &rn);
        const carquet_byte_array_t* regions = (const carquet_byte_array_t*)regdata;

        for (int64_t i = 0; i < n; i++) {
            returned++;
            if (matches(ids[i], (const char*)regions[i].data, (int)regions[i].length))
                exact++;
        }
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    int64_t skipped = carquet_batch_reader_rows_skipped(br);
    printf("  %-28s returned=%5lld  exact_match=%5lld  pruned=%5lld\n",
           label, (long long)returned, (long long)exact, (long long)skipped);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return 0;
}

static int match_id_gt(int64_t id, const char* region, int rlen) {
    (void)region; (void)rlen; return id > 7500;
}
static int match_region_east(int64_t id, const char* region, int rlen) {
    (void)id; return rlen == 4 && memcmp(region, "east", 4) == 0;
}
static int match_east_and_gt(int64_t id, const char* region, int rlen) {
    return match_id_gt(id, region, rlen) && match_region_east(id, region, rlen);
}

int main(void) {
    carquet_init();
    const char* path = "/tmp/carquet_page_filter_example.parquet";

    printf("Writing %d rows (%d per page, %d pages)...\n",
           TOTAL_ROWS, ROWS_PER_PAGE, TOTAL_ROWS / ROWS_PER_PAGE);
    if (write_file(path)) return 1;

    printf("\nPage-filtered scans (only matching pages are decoded):\n");

    /* id > 7500 — prunes the lower pages entirely. */
    int64_t threshold = 7500;
    carquet_filter_clause_t c_id = {0};
    c_id.column_index = 0;
    c_id.op = CARQUET_FILTER_GT;
    c_id.value = &threshold;
    c_id.value_size = (int32_t)sizeof(threshold);
    if (run_filter(path, "id > 7500", &c_id, 1, match_id_gt)) return 1;

    /* region == "east" — BYTE_ARRAY: raw bytes + length. */
    carquet_filter_clause_t c_reg = {0};
    c_reg.column_index = 1;
    c_reg.op = CARQUET_FILTER_EQ;
    c_reg.value = "east";
    c_reg.value_size = 4;
    if (run_filter(path, "region == 'east'", &c_reg, 1, match_region_east)) return 1;

    /* Conjunction: id > 7500 AND region == "east". */
    carquet_filter_clause_t conj[2] = { c_id, c_reg };
    if (run_filter(path, "id > 7500 AND east", conj, 2, match_east_and_gt)) return 1;

    printf("\nNote: 'exact_match' <= 'returned' because page filtering is\n"
           "conservative at page granularity; apply the predicate per row to\n"
           "drop the surviving non-matching rows.\n");

    remove(path);
    printf("\nDone.\n");
    return 0;
}
