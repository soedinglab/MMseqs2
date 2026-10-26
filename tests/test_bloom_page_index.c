/**
 * @file test_bloom_page_index.c
 * @brief Tests for bloom filter writing and page index writing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <carquet/carquet.h>
#include "reader/reader_internal.h"
#include "core/compat.h"
#include "test_helpers.h"

/* Bloom filter API (from bloom_filter.c) */
extern carquet_bloom_filter_t* carquet_bloom_filter_create(size_t num_bytes);
extern carquet_bloom_filter_t* carquet_bloom_filter_from_data(const uint8_t* data, size_t size);
extern void carquet_bloom_filter_insert_hash(carquet_bloom_filter_t* filter, uint64_t hash);
extern const uint8_t* carquet_bloom_filter_data(const carquet_bloom_filter_t* filter);
extern size_t carquet_bloom_filter_size(const carquet_bloom_filter_t* filter);
extern void carquet_bloom_filter_destroy(carquet_bloom_filter_t* filter);
extern bool carquet_bloom_filter_check_i32(const carquet_bloom_filter_t* filter, int32_t value);
extern bool carquet_bloom_filter_check_i64(const carquet_bloom_filter_t* filter, int64_t value);
extern bool carquet_bloom_filter_check_double(const carquet_bloom_filter_t* filter, double value);

typedef struct carquet_column_index_builder carquet_column_index_builder_t;
extern carquet_column_index_builder_t* carquet_column_index_builder_create(
    carquet_physical_type_t type, const carquet_logical_type_t* logical_type,
    int32_t type_length);
extern void carquet_column_index_builder_destroy(carquet_column_index_builder_t* builder);
extern carquet_status_t carquet_column_index_add_page(
    carquet_column_index_builder_t* builder,
    int64_t null_count, const void* min_value, int32_t min_value_len,
    const void* max_value, int32_t max_value_len, bool is_null_page,
    const int64_t* rep_level_hist, int32_t rep_level_hist_len,
    const int64_t* def_level_hist, int32_t def_level_hist_len);
extern carquet_status_t carquet_column_index_page_might_match(
    const carquet_column_index_builder_t* builder,
    int32_t page_idx,
    const void* min_value,
    const void* max_value,
    int32_t value_len,
    bool* might_match);

/* ============================================================================
 * Test: Bloom filter write + verify via metadata offsets
 * ============================================================================ */

static int test_bloom_filter_write(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "bloom_write");
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create schema: two columns */
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Write with bloom filters enabled */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_bloom_filters = true;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    int32_t ids[100];
    double values[100];
    for (int i = 0; i < 100; i++) {
        ids[i] = i * 10;
        values[i] = i * 1.5;
    }

    assert(carquet_writer_write_batch(writer, 0, ids, 100, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(writer, 1, values, 100, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back and verify metadata has bloom filter offsets */
    memset(&err, 0, sizeof(err));
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        char msg[512];
        snprintf(msg, sizeof(msg), "reader open failed: code=%d msg=%s", err.code, err.message);
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("bloom_filter_write", msg);
    }

    assert(carquet_reader_num_rows(reader) == 100);

    /* Check metadata: both columns should have bloom filter offsets */
    const parquet_file_metadata_t* meta = &reader->metadata;
    assert(meta->num_row_groups == 1);

    parquet_column_chunk_t* col0 = &meta->row_groups[0].columns[0];
    parquet_column_chunk_t* col1 = &meta->row_groups[0].columns[1];

    assert(col0->metadata.has_bloom_filter_offset);
    assert(col0->metadata.bloom_filter_offset > 0);
    assert(col0->metadata.has_bloom_filter_length);
    assert(col0->metadata.bloom_filter_length > 0);

    assert(col1->metadata.has_bloom_filter_offset);
    assert(col1->metadata.bloom_filter_offset > col0->metadata.bloom_filter_offset);
    assert(col1->metadata.has_bloom_filter_length);
    assert(col1->metadata.bloom_filter_length > 0);

    /* Read the raw bloom filter data from file and verify it works */
    FILE* f = fopen(test_file, "rb");
    assert(f);

    /* Read bloom filter for column 0 (skip Thrift header, read raw bitset) */
    /* The bloom filter header contains: numBytes (i32), algorithm, hash, compression
     * We can approximate the header size or just read the whole blob and use the
     * bloom filter API. For simplicity, let's seek to the offset and read. */
    int64_t bf0_offset = col0->metadata.bloom_filter_offset;
    int32_t bf0_length = col0->metadata.bloom_filter_length;

    uint8_t* bf0_data = malloc(bf0_length);
    assert(bf0_data);
    carquet_fseek64(f, bf0_offset, SEEK_SET);
    assert(fread(bf0_data, 1, bf0_length, f) == (size_t)bf0_length);

    /* The bloom filter data starts after the Thrift header.
     * Parse header to find numBytes, then create filter from the bitset portion.
     * The header is a Thrift compact struct. For a simpler test, we just verify
     * the data was written and the file is valid. */
    fclose(f);
    free(bf0_data);

    /* Verify data is still readable */
    carquet_column_reader_t* col_reader = carquet_reader_get_column(reader, 0, 0, &err);
    assert(col_reader);
    int32_t read_ids[100];
    int64_t count = carquet_column_read_batch(col_reader, read_ids, 100, NULL, NULL);
    assert(count == 100);
    assert(read_ids[0] == 0);
    assert(read_ids[50] == 500);
    assert(read_ids[99] == 990);
    carquet_column_reader_free(col_reader);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("bloom_filter_write");
    return 0;
}

/* ============================================================================
 * Test: Page index write + verify via metadata offsets
 * ============================================================================ */

static int test_page_index_write(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "page_index_write");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Write with page index enabled, small page size to force multiple pages */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 256;  /* Very small to force multiple pages */

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    int32_t data[1000];
    for (int i = 0; i < 1000; i++) {
        data[i] = i;
    }

    assert(carquet_writer_write_batch(writer, 0, data, 1000, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back and check metadata */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("page_index_write", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 1000);

    const parquet_file_metadata_t* meta = &reader->metadata;
    parquet_column_chunk_t* col = &meta->row_groups[0].columns[0];

    /* Verify column index and offset index metadata present */
    assert(col->has_column_index_offset);
    assert(col->column_index_offset > 0);
    assert(col->has_column_index_length);
    assert(col->column_index_length > 0);

    assert(col->has_offset_index_offset);
    assert(col->offset_index_offset > 0);
    assert(col->has_offset_index_length);
    assert(col->offset_index_length > 0);

    /* Column index should come before offset index in the file */
    assert(col->column_index_offset < col->offset_index_offset);

    /* Verify data integrity */
    carquet_column_reader_t* col_reader = carquet_reader_get_column(reader, 0, 0, &err);
    assert(col_reader);
    int32_t read_data[1000];
    int64_t count = carquet_column_read_batch(col_reader, read_data, 1000, NULL, NULL);
    assert(count == 1000);
    assert(read_data[0] == 0);
    assert(read_data[999] == 999);
    carquet_column_reader_free(col_reader);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("page_index_write");
    return 0;
}

/* ============================================================================
 * Test: Both bloom filter and page index together
 * ============================================================================ */

static int test_bloom_and_page_index(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "bloom_page_both");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT64,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_FLOAT,
                              NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_bloom_filters = true;
    opts.write_page_index = true;
    opts.page_size = 512;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    int64_t a_vals[500];
    float b_vals[500];
    int16_t b_def[500];
    for (int i = 0; i < 500; i++) {
        a_vals[i] = (int64_t)i * 100;
        b_vals[i] = (float)i * 0.5f;
        b_def[i] = (i % 5 == 0) ? 0 : 1;  /* Every 5th value is NULL */
    }

    assert(carquet_writer_write_batch(writer, 0, a_vals, 500, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(writer, 1, b_vals, 500, b_def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("bloom_and_page_index", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 500);

    const parquet_file_metadata_t* meta = &reader->metadata;

    /* Both columns should have bloom filter and page index metadata */
    for (int i = 0; i < 2; i++) {
        parquet_column_chunk_t* col = &meta->row_groups[0].columns[i];

        assert(col->metadata.has_bloom_filter_offset);
        assert(col->metadata.bloom_filter_offset > 0);
        assert(col->metadata.has_bloom_filter_length);
        assert(col->metadata.bloom_filter_length > 0);

        assert(col->has_column_index_offset);
        assert(col->column_index_offset > 0);
        assert(col->has_column_index_length);
        assert(col->column_index_length > 0);

        assert(col->has_offset_index_offset);
        assert(col->offset_index_offset > 0);
        assert(col->has_offset_index_length);
        assert(col->offset_index_length > 0);
    }

    /* Verify data round-trip for column 0 */
    carquet_column_reader_t* col_reader = carquet_reader_get_column(reader, 0, 0, &err);
    assert(col_reader);
    int64_t read_a[500];
    int64_t count = carquet_column_read_batch(col_reader, read_a, 500, NULL, NULL);
    assert(count == 500);
    assert(read_a[0] == 0);
    assert(read_a[499] == 49900);
    carquet_column_reader_free(col_reader);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("bloom_and_page_index");
    return 0;
}

/* ============================================================================
 * Test: Bloom filter with compression
 * ============================================================================ */

static int test_bloom_filter_compressed(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "bloom_compressed");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "val", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_bloom_filters = true;
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    int32_t data[200];
    for (int i = 0; i < 200; i++) data[i] = i * 7;

    assert(carquet_writer_write_batch(writer, 0, data, 200, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("bloom_filter_compressed", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 200);

    parquet_column_chunk_t* col = &reader->metadata.row_groups[0].columns[0];
    assert(col->metadata.has_bloom_filter_offset);
    assert(col->metadata.bloom_filter_length > 0);

    /* Read data back */
    carquet_column_reader_t* col_reader = carquet_reader_get_column(reader, 0, 0, &err);
    assert(col_reader);
    int32_t read_data[200];
    int64_t count = carquet_column_read_batch(col_reader, read_data, 200, NULL, NULL);
    assert(count == 200);
    assert(read_data[100] == 700);
    carquet_column_reader_free(col_reader);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("bloom_filter_compressed");
    return 0;
}

/* ============================================================================
 * Test: Page index with multiple row groups
 * ============================================================================ */

static int test_page_index_multi_row_groups(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "page_idx_multi_rg");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.page_size = 256;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    /* Write first row group */
    int32_t data1[500];
    for (int i = 0; i < 500; i++) data1[i] = i;
    assert(carquet_writer_write_batch(writer, 0, data1, 500, NULL, NULL) == CARQUET_OK);

    /* Force new row group */
    assert(carquet_writer_new_row_group(writer) == CARQUET_OK);

    /* Write second row group */
    int32_t data2[300];
    for (int i = 0; i < 300; i++) data2[i] = 1000 + i;
    assert(carquet_writer_write_batch(writer, 0, data2, 300, NULL, NULL) == CARQUET_OK);

    assert(carquet_writer_close(writer) == CARQUET_OK);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("page_index_multi_row_groups", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 800);
    assert(reader->metadata.num_row_groups == 2);

    /* Both row groups should have page index metadata */
    for (int rg = 0; rg < 2; rg++) {
        parquet_column_chunk_t* col = &reader->metadata.row_groups[rg].columns[0];
        assert(col->has_column_index_offset);
        assert(col->column_index_offset > 0);
        assert(col->has_offset_index_offset);
        assert(col->offset_index_offset > 0);
    }

    /* Second row group offsets should be after first */
    parquet_column_chunk_t* rg0 = &reader->metadata.row_groups[0].columns[0];
    parquet_column_chunk_t* rg1 = &reader->metadata.row_groups[1].columns[0];
    assert(rg1->column_index_offset > rg0->offset_index_offset);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("page_index_multi_row_groups");
    return 0;
}

/* ============================================================================
 * Test: Bloom filter disabled by default
 * ============================================================================ */

static int test_bloom_filter_disabled_default(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "bloom_disabled");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "val", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* Default options - no bloom filters */
    carquet_writer_t* writer = carquet_writer_create(test_file, schema, NULL, &err);
    assert(writer);

    int32_t data[50];
    for (int i = 0; i < 50; i++) data[i] = i;
    assert(carquet_writer_write_batch(writer, 0, data, 50, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    assert(reader);

    parquet_column_chunk_t* col = &reader->metadata.row_groups[0].columns[0];
    assert(!col->metadata.has_bloom_filter_offset);
    assert(!col->has_column_index_offset);
    assert(!col->has_offset_index_offset);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("bloom_filter_disabled_default");
    return 0;
}

/* ============================================================================
 * Test: Bloom filter with nullable column
 * ============================================================================ */

static int test_bloom_filter_nullable(void) {
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "bloom_nullable");
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema);
    carquet_schema_add_column(schema, "val", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_bloom_filters = true;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer);

    /* Sparse: only non-null values in array */
    int32_t data[] = {10, 20, 30, 40};
    int16_t def[] = {1, 0, 1, 1, 0, 1};  /* 6 rows: values at indices 0,2,3,5 */
    assert(carquet_writer_write_batch(writer, 0, data, 6, def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        TEST_FAIL("bloom_filter_nullable", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 6);

    parquet_column_chunk_t* col = &reader->metadata.row_groups[0].columns[0];
    assert(col->metadata.has_bloom_filter_offset);
    assert(col->metadata.bloom_filter_length > 0);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("bloom_filter_nullable");
    return 0;
}

/* ============================================================================
 * Test: Bloom filter block selection follows Parquet spec
 * ============================================================================ */

static int test_bloom_filter_spec_block_index(void) {
    carquet_bloom_filter_t* filter = carquet_bloom_filter_create(96);
    if (!filter) {
        TEST_FAIL("bloom_filter_spec_block_index", "create failed");
    }

    /* With three blocks and top hash bits 0x80000000, the Parquet formula
     * (top_bits * num_blocks) >> 32 selects block 1. A modulo-based selector
     * would select block 2, so this catches that compatibility bug. */
    carquet_bloom_filter_insert_hash(filter, 0x8000000000000001ULL);

    const uint8_t* data = carquet_bloom_filter_data(filter);
    size_t size = carquet_bloom_filter_size(filter);
    if (!data || size != 96) {
        carquet_bloom_filter_destroy(filter);
        TEST_FAIL("bloom_filter_spec_block_index", "unexpected filter size");
    }

    bool block0_set = false;
    bool block1_set = false;
    bool block2_set = false;
    for (int i = 0; i < 32; i++) {
        block0_set = block0_set || data[i] != 0;
        block1_set = block1_set || data[32 + i] != 0;
        block2_set = block2_set || data[64 + i] != 0;
    }

    carquet_bloom_filter_destroy(filter);

    if (block0_set || !block1_set || block2_set) {
        TEST_FAIL("bloom_filter_spec_block_index", "hash inserted into wrong block");
    }

    TEST_PASS("bloom_filter_spec_block_index");
    return 0;
}

static int test_page_index_unsigned_logical_order(void) {
    const char* tname = "page_index_unsigned_logical_order";
    carquet_logical_type_t u32_type = {
        .id = CARQUET_LOGICAL_INTEGER,
        .params.integer = { .bit_width = 32, .is_signed = false }
    };

    carquet_column_index_builder_t* builder =
        carquet_column_index_builder_create(CARQUET_PHYSICAL_INT32, &u32_type, 0);
    if (!builder) TEST_FAIL(tname, "builder create failed");

    uint32_t min_v = 1u;
    uint32_t max_v = UINT32_MAX;
    if (carquet_column_index_add_page(builder, 0, &min_v, 4, &max_v, 4, false,
                                      NULL, 0, NULL, 0) != CARQUET_OK) {
        carquet_column_index_builder_destroy(builder);
        TEST_FAIL(tname, "add page failed");
    }

    bool might_match = false;
    uint32_t query_min = UINT32_MAX - 1u;
    uint32_t query_max = UINT32_MAX;
    if (carquet_column_index_page_might_match(
            builder, 0, &query_min, &query_max, 4, &might_match) != CARQUET_OK ||
        !might_match) {
        carquet_column_index_builder_destroy(builder);
        TEST_FAIL(tname, "unsigned high range was incorrectly rejected");
    }

    query_min = 0u;
    query_max = 0u;
    if (carquet_column_index_page_might_match(
            builder, 0, &query_min, &query_max, 4, &might_match) != CARQUET_OK ||
        might_match) {
        carquet_column_index_builder_destroy(builder);
        TEST_FAIL(tname, "unsigned low range was not rejected");
    }

    carquet_column_index_builder_destroy(builder);
    TEST_PASS(tname);
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("=== Bloom Filter & Page Index Tests ===\n\n");

    int failures = 0;
    failures += test_bloom_filter_write();
    failures += test_page_index_write();
    failures += test_bloom_and_page_index();
    failures += test_bloom_filter_compressed();
    failures += test_page_index_multi_row_groups();
    failures += test_bloom_filter_disabled_default();
    failures += test_bloom_filter_nullable();
    failures += test_bloom_filter_spec_block_index();
    failures += test_page_index_unsigned_logical_order();

    printf("\n");
    if (failures == 0) {
        printf("All bloom filter & page index tests passed!\n");
    } else {
        printf("%d test(s) FAILED\n", failures);
    }

    return failures;
}
