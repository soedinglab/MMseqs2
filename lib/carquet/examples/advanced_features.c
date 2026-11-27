/**
 * @file advanced_features.c
 * @brief Compact examples of carquet advanced API features
 *
 * Demonstrates: bloom filters, page indexes, key-value metadata,
 * column chunk metadata, per-column writer options, and buffer writer.
 */

#include <carquet/carquet.h>
#include "example_paths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_ROWS 10000
#define CHECK(expr) do { if ((expr) != CARQUET_OK) { \
    fprintf(stderr, "FAIL at %s:%d\n", __FILE__, __LINE__); return 1; } } while(0)

/* ---------- helpers --------------------------------------------------- */

static carquet_schema_t* make_schema(void) {
    carquet_schema_t* s = carquet_schema_create(NULL);
    carquet_schema_add_column(s, "id",    CARQUET_PHYSICAL_INT64,      NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "value", CARQUET_PHYSICAL_DOUBLE,     NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "name",  CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    return s;
}

static int write_test_file(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = make_schema();

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression        = CARQUET_COMPRESSION_ZSTD;
    opts.write_statistics   = true;
    opts.write_bloom_filters = true;
    opts.write_page_index   = true;

    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) { fprintf(stderr, "writer: %s\n", err.message); exit(1); }

    /* -- per-column overrides ------------------------------------------ */
    carquet_writer_set_column_compression(w, 0, CARQUET_COMPRESSION_UNCOMPRESSED, 0);
    carquet_writer_set_column_bloom_filter(w, 2, false);  /* no bloom on 'name' */

    /* -- key-value metadata -------------------------------------------- */
    carquet_writer_add_metadata(w, "app", "carquet-example");
    carquet_writer_add_metadata(w, "version", CARQUET_VERSION_STRING);

    /* -- write data ---------------------------------------------------- */
    int64_t  ids[NUM_ROWS];
    double   vals[NUM_ROWS];
    for (int i = 0; i < NUM_ROWS; i++) { ids[i] = i; vals[i] = i * 0.1; }

    CHECK(carquet_writer_write_batch(w, 0, ids,  NUM_ROWS, NULL, NULL));
    CHECK(carquet_writer_write_batch(w, 1, vals, NUM_ROWS, NULL, NULL));

    /* name column: all non-null for simplicity */
    const char* tag = "row";
    carquet_byte_array_t names[NUM_ROWS];
    int16_t def[NUM_ROWS];
    for (int i = 0; i < NUM_ROWS; i++) {
        names[i].data   = (uint8_t*)tag;
        names[i].length = 3;
        def[i] = 1;
    }
    CHECK(carquet_writer_write_batch(w, 2, names, NUM_ROWS, def, NULL));
    CHECK(carquet_writer_close(w));
    carquet_schema_free(schema);
    (void)0; /* suppress -Wreturn-type: all CHECK calls above return 1 on failure */
    return 0;
}

/* ---------- demo: bloom filter --------------------------------------- */

static void demo_bloom_filter(const char* path) {
    printf("\n=== Bloom Filter ===\n");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);

    carquet_bloom_filter_t* bf = carquet_reader_get_bloom_filter(r, 0, 0, &err);
    if (bf) {
        printf("  filter size: %llu bytes\n", (unsigned long long)carquet_bloom_filter_size(bf));
        int64_t v42 = 42, v_neg = -999;
        printf("  might contain 42:   %s\n", carquet_bloom_filter_check_i64(bf, v42)  ? "yes" : "no");
        printf("  might contain -999: %s\n", carquet_bloom_filter_check_i64(bf, v_neg) ? "yes" : "no");
        carquet_bloom_filter_destroy(bf);
    } else {
        printf("  (no bloom filter on column 0)\n");
    }
    carquet_reader_close(r);
}

/* ---------- demo: page index ----------------------------------------- */

static void demo_page_index(const char* path) {
    printf("\n=== Page Index ===\n");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(path, NULL, &err);

    carquet_column_index_t* ci = carquet_reader_get_column_index(r, 0, 0, &err);
    if (ci) {
        int32_t np = carquet_column_index_num_pages(ci);
        printf("  column 0: %d pages, boundary_order=%d\n",
               np, carquet_column_index_boundary_order(ci));
        for (int32_t p = 0; p < np && p < 3; p++) {
            carquet_page_stats_t st;
            carquet_column_index_get_page_stats(ci, p, &st);
            int64_t mn = st.min_value ? *(const int64_t*)st.min_value : -1;
            int64_t mx = st.max_value ? *(const int64_t*)st.max_value : -1;
            printf("    page %d: min=%lld max=%lld nulls=%lld\n",
                   p, (long long)mn, (long long)mx, (long long)st.null_count);
        }
        carquet_column_index_free(ci);
    } else {
        printf("  (no column index)\n");
    }

    carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, 0, 0, &err);
    if (oi) {
        int32_t np = carquet_offset_index_num_pages(oi);
        printf("  offset index: %d pages\n", np);
        for (int32_t p = 0; p < np && p < 3; p++) {
            carquet_page_location_t loc;
            carquet_offset_index_get_page_location(oi, p, &loc);
            printf("    page %d: offset=%lld size=%d first_row=%lld\n",
                   p, (long long)loc.offset, loc.compressed_size,
                   (long long)loc.first_row_index);
        }
        carquet_offset_index_free(oi);
    }
    carquet_reader_close(r);
}

/* ---------- demo: key-value metadata --------------------------------- */

static void demo_kv_metadata(const char* path) {
    printf("\n=== Key-Value Metadata ===\n");
    carquet_reader_t* r = carquet_reader_open(path, NULL, NULL);
    int32_t n = carquet_reader_num_metadata(r);
    printf("  %d entries:\n", n);
    for (int32_t i = 0; i < n; i++) {
        const char *k, *v;
        carquet_reader_get_metadata(r, i, &k, &v);
        printf("    %s = %s\n", k, v ? v : "(null)");
    }
    const char* app = carquet_reader_find_metadata(r, "app");
    printf("  find 'app': %s\n", app ? app : "(not found)");
    carquet_reader_close(r);
}

/* ---------- demo: column chunk metadata ------------------------------ */

static void demo_chunk_metadata(const char* path) {
    printf("\n=== Column Chunk Metadata ===\n");
    carquet_reader_t* r = carquet_reader_open(path, NULL, NULL);
    for (int32_t c = 0; c < carquet_reader_num_columns(r); c++) {
        carquet_column_chunk_metadata_t m;
        carquet_reader_column_chunk_metadata(r, 0, c, &m);
        printf("  col %d: codec=%s  values=%lld  compressed=%lld  "
               "bloom=%s  col_idx=%s\n",
               c, carquet_compression_name(m.codec),
               (long long)m.num_values, (long long)m.total_compressed_size,
               m.has_bloom_filter ? "yes" : "no",
               m.has_column_index ? "yes" : "no");
    }
    carquet_reader_close(r);
}

/* ---------- demo: buffer writer -------------------------------------- */

static void demo_buffer_writer(void) {
    printf("\n=== Buffer Writer ===\n");
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(NULL);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) { printf("  buffer writer not available: %s\n", err.message); goto done; }

    int32_t data[] = {10, 20, 30};
    carquet_writer_write_batch(w, 0, data, 3, NULL, NULL);
    carquet_writer_close(w);

    void* buf = NULL;
    size_t sz = 0;
    if (carquet_writer_get_buffer(w, &buf, &sz) == CARQUET_OK) {
        printf("  wrote %llu bytes to memory\n", (unsigned long long)sz);
        /* Round-trip: read from buffer */
        carquet_reader_t* r = carquet_reader_open_buffer(buf, sz, NULL, NULL);
        if (r) {
            printf("  round-trip: %lld rows, %d cols\n",
                   (long long)carquet_reader_num_rows(r),
                   carquet_reader_num_columns(r));
            carquet_reader_close(r);
        }
        carquet_free(buf);
    }
done:
    carquet_schema_free(schema);
}

/* ---------- main ----------------------------------------------------- */

int main(void) {
    carquet_init();
    char path[1024];
    example_temp_path(path, sizeof(path), "carquet_advanced_example.parquet");

    printf("Writing test file with advanced features...\n");
    if (write_test_file(path)) return 1;

    demo_bloom_filter(path);
    demo_page_index(path);
    demo_kv_metadata(path);
    demo_chunk_metadata(path);
    demo_buffer_writer();

    remove(path);
    printf("\nDone.\n");
    return 0;
}
