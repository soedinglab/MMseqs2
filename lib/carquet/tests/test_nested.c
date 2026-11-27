/**
 * @file test_nested.c
 * @brief Comprehensive tests for nested schema support
 *
 * Tests schema creation, def/rep level computation, write/read round-trip,
 * and interop for nested Parquet schemas including optional groups,
 * repeated groups, and deeply nested structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <carquet/carquet.h>
#include "test_helpers.h"

/* ============================================================================
 * Test 1: Schema-level accumulated def/rep level accessors
 * ============================================================================
 *
 * Schema:
 *   schema (root, required)
 *   ├── id (required, INT32)               -> def=0, rep=0
 *   ├── name (optional, BYTE_ARRAY)        -> def=1, rep=0
 *   ├── address (optional, group)
 *   │   ├── street (required, BYTE_ARRAY)  -> def=1, rep=0
 *   │   └── city (optional, BYTE_ARRAY)    -> def=2, rep=0
 *   └── phones (repeated, group)
 *       ├── number (required, BYTE_ARRAY)  -> def=1, rep=1
 *       └── type (optional, BYTE_ARRAY)    -> def=2, rep=1
 */
static int test_accumulated_levels(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* id: required at root */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    /* name: optional at root */
    assert(carquet_schema_add_column(schema, "name", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, 0) == CARQUET_OK);

    /* address: optional group */
    int32_t address_idx = carquet_schema_add_group(schema, "address",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(address_idx >= 0);

    /* street: required under optional group */
    assert(carquet_schema_add_column(schema, "street", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_REQUIRED, 0, address_idx) == CARQUET_OK);

    /* city: optional under optional group */
    assert(carquet_schema_add_column(schema, "city", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, address_idx) == CARQUET_OK);

    /* phones: repeated group */
    int32_t phones_idx = carquet_schema_add_group(schema, "phones",
        CARQUET_REPETITION_REPEATED, 0);
    assert(phones_idx >= 0);

    /* number: required under repeated group */
    assert(carquet_schema_add_column(schema, "number", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_REQUIRED, 0, phones_idx) == CARQUET_OK);

    /* type: optional under repeated group */
    assert(carquet_schema_add_column(schema, "type", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, phones_idx) == CARQUET_OK);

    /* Verify leaf count */
    assert(carquet_schema_num_columns(schema) == 6);

    /* Verify accumulated levels via new schema-level API */
    /* Column 0 (id): required at root -> def=0, rep=0 */
    assert(carquet_schema_max_def_level(schema, 0) == 0);
    assert(carquet_schema_max_rep_level(schema, 0) == 0);

    /* Column 1 (name): optional at root -> def=1, rep=0 */
    assert(carquet_schema_max_def_level(schema, 1) == 1);
    assert(carquet_schema_max_rep_level(schema, 1) == 0);

    /* Column 2 (street): required under optional group -> def=1, rep=0 */
    assert(carquet_schema_max_def_level(schema, 2) == 1);
    assert(carquet_schema_max_rep_level(schema, 2) == 0);

    /* Column 3 (city): optional under optional group -> def=2, rep=0 */
    assert(carquet_schema_max_def_level(schema, 3) == 2);
    assert(carquet_schema_max_rep_level(schema, 3) == 0);

    /* Column 4 (number): required under repeated group -> def=1, rep=1 */
    assert(carquet_schema_max_def_level(schema, 4) == 1);
    assert(carquet_schema_max_rep_level(schema, 4) == 1);

    /* Column 5 (type): optional under repeated group -> def=2, rep=1 */
    assert(carquet_schema_max_def_level(schema, 5) == 2);
    assert(carquet_schema_max_rep_level(schema, 5) == 1);

    /* Verify column names */
    assert(strcmp(carquet_schema_column_name(schema, 0), "id") == 0);
    assert(strcmp(carquet_schema_column_name(schema, 1), "name") == 0);
    assert(strcmp(carquet_schema_column_name(schema, 2), "street") == 0);
    assert(strcmp(carquet_schema_column_name(schema, 3), "city") == 0);
    assert(strcmp(carquet_schema_column_name(schema, 4), "number") == 0);
    assert(strcmp(carquet_schema_column_name(schema, 5), "type") == 0);

    /* Verify out-of-bounds returns */
    assert(carquet_schema_max_def_level(schema, -1) == -1);
    assert(carquet_schema_max_def_level(schema, 6) == -1);
    assert(carquet_schema_column_name(schema, -1) == NULL);

    carquet_schema_free(schema);
    TEST_PASS("accumulated_levels");
    return 0;
}

/* ============================================================================
 * Test 2: Column path computation
 * ============================================================================
 */
static int test_column_paths(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* Flat column at root */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    /* Nested: address.street */
    int32_t addr = carquet_schema_add_group(schema, "address",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(addr >= 0);
    assert(carquet_schema_add_column(schema, "street", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_REQUIRED, 0, addr) == CARQUET_OK);

    /* Deeply nested: contact.info.email */
    int32_t contact = carquet_schema_add_group(schema, "contact",
        CARQUET_REPETITION_OPTIONAL, 0);
    int32_t info = carquet_schema_add_group(schema, "info",
        CARQUET_REPETITION_REQUIRED, contact);
    assert(carquet_schema_add_column(schema, "email", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, info) == CARQUET_OK);

    const char* path[16];
    int32_t depth;

    /* Column 0 (id): path = ["id"] */
    depth = carquet_schema_column_path(schema, 0, path, 16);
    assert(depth == 1);
    assert(strcmp(path[0], "id") == 0);

    /* Column 1 (street): path = ["address", "street"] */
    depth = carquet_schema_column_path(schema, 1, path, 16);
    assert(depth == 2);
    assert(strcmp(path[0], "address") == 0);
    assert(strcmp(path[1], "street") == 0);

    /* Column 2 (email): path = ["contact", "info", "email"] */
    depth = carquet_schema_column_path(schema, 2, path, 16);
    assert(depth == 3);
    assert(strcmp(path[0], "contact") == 0);
    assert(strcmp(path[1], "info") == 0);
    assert(strcmp(path[2], "email") == 0);

    carquet_schema_free(schema);
    TEST_PASS("column_paths");
    return 0;
}

/* ============================================================================
 * Test 3: Deeply nested levels (3+ levels)
 * ============================================================================
 *
 * Schema:
 *   schema (root)
 *   └── a (optional group)
 *       └── b (repeated group)
 *           └── c (optional group)
 *               └── value (required, INT32)  -> def=3, rep=1
 */
static int test_deep_nesting_levels(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    int32_t a = carquet_schema_add_group(schema, "a", CARQUET_REPETITION_OPTIONAL, 0);
    int32_t b = carquet_schema_add_group(schema, "b", CARQUET_REPETITION_REPEATED, a);
    int32_t c = carquet_schema_add_group(schema, "c", CARQUET_REPETITION_OPTIONAL, b);
    assert(carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, c) == CARQUET_OK);

    assert(carquet_schema_num_columns(schema) == 1);

    /* def = optional(a)=1 + repeated(b)=1 + optional(c)=1 + required(value)=0 = 3 */
    assert(carquet_schema_max_def_level(schema, 0) == 3);
    /* rep = repeated(b)=1 = 1 */
    assert(carquet_schema_max_rep_level(schema, 0) == 1);

    carquet_schema_free(schema);
    TEST_PASS("deep_nesting_levels");
    return 0;
}

/* ============================================================================
 * Test 4: Multiple repeated ancestors
 * ============================================================================
 *
 * Schema:
 *   schema (root)
 *   └── list_a (repeated group)
 *       └── list_b (repeated group)
 *           └── val (optional, INT32)  -> def=3, rep=2
 */
static int test_multiple_repeated_ancestors(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    int32_t la = carquet_schema_add_group(schema, "list_a", CARQUET_REPETITION_REPEATED, 0);
    int32_t lb = carquet_schema_add_group(schema, "list_b", CARQUET_REPETITION_REPEATED, la);
    assert(carquet_schema_add_column(schema, "val", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, lb) == CARQUET_OK);

    /* def = repeated(la)=1 + repeated(lb)=1 + optional(val)=1 = 3 */
    assert(carquet_schema_max_def_level(schema, 0) == 3);
    /* rep = repeated(la)=1 + repeated(lb)=1 = 2 */
    assert(carquet_schema_max_rep_level(schema, 0) == 2);

    carquet_schema_free(schema);
    TEST_PASS("multiple_repeated_ancestors");
    return 0;
}

/* ============================================================================
 * Test 5: Write/read round-trip with optional group (nullable nested data)
 * ============================================================================
 *
 * Schema:
 *   schema (root)
 *   ├── id (required, INT32)
 *   └── address (optional group)
 *       └── zip (required, INT32)
 *
 * Data (5 rows):
 *   Row 0: id=1, address.zip=10001
 *   Row 1: id=2, address=NULL
 *   Row 2: id=3, address.zip=10003
 *   Row 3: id=4, address=NULL
 *   Row 4: id=5, address.zip=10005
 */
static int test_roundtrip_optional_group(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_optional");

    /* Create schema */
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t addr = carquet_schema_add_group(schema, "address",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(addr >= 0);

    assert(carquet_schema_add_column(schema, "zip", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, addr) == CARQUET_OK);

    /* Verify levels: zip has def=1 (from optional address), rep=0 */
    assert(carquet_schema_max_def_level(schema, 1) == 1);
    assert(carquet_schema_max_rep_level(schema, 1) == 0);

    /* Write file */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("  Writer creation failed: %s\n", err.message);
        TEST_FAIL("roundtrip_optional_group", "writer creation failed");
    }

    /* Column 0: id (required, 5 values) */
    int32_t ids[] = {1, 2, 3, 4, 5};
    assert(carquet_writer_write_batch(writer, 0, ids, 5, NULL, NULL) == CARQUET_OK);

    /* Column 1: zip (under optional address group, def_level max=1)
     * def_level=1 means address is present (zip has a value)
     * def_level=0 means address is NULL (zip has no value)
     * Non-null values only: 10001, 10003, 10005 (3 values, but 5 rows with levels) */
    int32_t zips[] = {10001, 10003, 10005};
    int16_t zip_def[] = {1, 0, 1, 0, 1};  /* present, null, present, null, present */
    assert(carquet_writer_write_batch(writer, 1, zips, 5, zip_def, NULL) == CARQUET_OK);

    carquet_status_t status = carquet_writer_close(writer);
    assert(status == CARQUET_OK);

    /* Read file back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("roundtrip_optional_group", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 5);
    assert(carquet_reader_num_columns(reader) == 2);

    /* Verify schema levels from reader */
    const carquet_schema_t* read_schema = carquet_reader_schema(reader);
    assert(carquet_schema_max_def_level(read_schema, 0) == 0);  /* id: required */
    assert(carquet_schema_max_def_level(read_schema, 1) == 1);  /* zip: under optional group */
    assert(carquet_schema_max_rep_level(read_schema, 1) == 0);

    /* Read id column */
    carquet_column_reader_t* id_col = carquet_reader_get_column(reader, 0, 0, &err);
    assert(id_col != NULL);
    int32_t read_ids[5];
    int64_t count = carquet_column_read_batch(id_col, read_ids, 5, NULL, NULL);
    assert(count == 5);
    for (int i = 0; i < 5; i++) {
        assert(read_ids[i] == ids[i]);
    }
    carquet_column_reader_free(id_col);

    /* Read zip column with def levels */
    carquet_column_reader_t* zip_col = carquet_reader_get_column(reader, 0, 1, &err);
    assert(zip_col != NULL);
    int32_t read_zips[5];
    int16_t read_def[5];
    count = carquet_column_read_batch(zip_col, read_zips, 5, read_def, NULL);
    assert(count == 5);

    /* Verify definition levels */
    assert(read_def[0] == 1);  /* present */
    assert(read_def[1] == 0);  /* null */
    assert(read_def[2] == 1);  /* present */
    assert(read_def[3] == 0);  /* null */
    assert(read_def[4] == 1);  /* present */

    /* Verify non-null values */
    int val_idx = 0;
    int32_t expected_zips[] = {10001, 10003, 10005};
    for (int i = 0; i < 5; i++) {
        if (read_def[i] == 1) {
            assert(read_zips[val_idx] == expected_zips[val_idx]);
            val_idx++;
        }
    }

    carquet_column_reader_free(zip_col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("roundtrip_optional_group");
    return 0;
}

/* ============================================================================
 * Test 6: Write/read round-trip with flat required columns + nested
 * ============================================================================
 *
 * Schema:
 *   schema (root)
 *   ├── x (required, INT32)
 *   ├── y (required, DOUBLE)
 *   └── meta (optional group)
 *       └── tag (optional, INT32)
 *
 * "tag" has def=2: optional(meta)=1 + optional(tag)=1
 */
static int test_roundtrip_double_optional(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_double_opt");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    assert(carquet_schema_add_column(schema, "x", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);
    assert(carquet_schema_add_column(schema, "y", CARQUET_PHYSICAL_DOUBLE,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t meta = carquet_schema_add_group(schema, "meta",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(meta >= 0);
    assert(carquet_schema_add_column(schema, "tag", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, meta) == CARQUET_OK);

    /* tag: def=2 (optional group + optional field), rep=0 */
    assert(carquet_schema_max_def_level(schema, 2) == 2);
    assert(carquet_schema_max_rep_level(schema, 2) == 0);

    /* Write */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        TEST_FAIL("roundtrip_double_optional", "writer creation failed");
    }

    int32_t xs[] = {10, 20, 30, 40};
    double ys[] = {1.1, 2.2, 3.3, 4.4};
    assert(carquet_writer_write_batch(writer, 0, xs, 4, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(writer, 1, ys, 4, NULL, NULL) == CARQUET_OK);

    /* tag with def=2 means both meta group and tag value are present
     * def=1 means meta group is present but tag is null
     * def=0 means meta group is null (so tag is also null)
     *
     * Row 0: meta present, tag=100    (def=2)
     * Row 1: meta null                (def=0)
     * Row 2: meta present, tag=null   (def=1)
     * Row 3: meta present, tag=400    (def=2) */
    int32_t tags[] = {100, 400};  /* only non-null values */
    int16_t tag_def[] = {2, 0, 1, 2};
    assert(carquet_writer_write_batch(writer, 2, tags, 4, tag_def, NULL) == CARQUET_OK);

    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("roundtrip_double_optional", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 4);
    assert(carquet_reader_num_columns(reader) == 3);

    /* Verify levels from reader schema */
    const carquet_schema_t* rs = carquet_reader_schema(reader);
    assert(carquet_schema_max_def_level(rs, 2) == 2);

    /* Read tag column */
    carquet_column_reader_t* tag_col = carquet_reader_get_column(reader, 0, 2, &err);
    assert(tag_col != NULL);
    int32_t read_tags[4];
    int16_t read_def[4];
    int64_t count = carquet_column_read_batch(tag_col, read_tags, 4, read_def, NULL);
    assert(count == 4);

    /* Verify def levels preserved */
    assert(read_def[0] == 2);
    assert(read_def[1] == 0);
    assert(read_def[2] == 1);
    assert(read_def[3] == 2);

    carquet_column_reader_free(tag_col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("roundtrip_double_optional");
    return 0;
}

/* ============================================================================
 * Test 7: Schema with mixed flat and nested columns, verify element count
 * ============================================================================
 */
static int test_mixed_schema_structure(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* root (1) + id (1) + address group (1) + street (1) + city (1) = 5 elements */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t addr = carquet_schema_add_group(schema, "address",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(carquet_schema_add_column(schema, "street", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_REQUIRED, 0, addr) == CARQUET_OK);
    assert(carquet_schema_add_column(schema, "city", CARQUET_PHYSICAL_BYTE_ARRAY,
        NULL, CARQUET_REPETITION_OPTIONAL, 0, addr) == CARQUET_OK);

    assert(carquet_schema_num_elements(schema) == 5);  /* root + id + address + street + city */
    assert(carquet_schema_num_columns(schema) == 3);    /* id, street, city */

    /* Verify element types */
    const carquet_schema_node_t* root = carquet_schema_get_element(schema, 0);
    assert(!carquet_schema_node_is_leaf(root));

    const carquet_schema_node_t* id_node = carquet_schema_get_element(schema, 1);
    assert(carquet_schema_node_is_leaf(id_node));
    assert(carquet_schema_node_physical_type(id_node) == CARQUET_PHYSICAL_INT32);

    const carquet_schema_node_t* addr_node = carquet_schema_get_element(schema, 2);
    assert(!carquet_schema_node_is_leaf(addr_node));
    assert(carquet_schema_node_repetition(addr_node) == CARQUET_REPETITION_OPTIONAL);

    carquet_schema_free(schema);
    TEST_PASS("mixed_schema_structure");
    return 0;
}

/* ============================================================================
 * Test 8: Write/read nested schema and verify schema is preserved in metadata
 * ============================================================================
 */
static int test_schema_metadata_roundtrip(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_metadata");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t group = carquet_schema_add_group(schema, "data",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(carquet_schema_add_column(schema, "val", CARQUET_PHYSICAL_DOUBLE,
        NULL, CARQUET_REPETITION_REQUIRED, 0, group) == CARQUET_OK);

    /* Write minimal data */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer != NULL);

    int32_t id = 42;
    double val = 3.14;
    int16_t val_def[] = {1};
    assert(carquet_writer_write_batch(writer, 0, &id, 1, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(writer, 1, &val, 1, val_def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back and verify schema structure */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("schema_metadata_roundtrip", "reader open failed");
    }

    const carquet_schema_t* rs = carquet_reader_schema(reader);
    assert(carquet_schema_num_columns(rs) == 2);
    assert(carquet_schema_num_elements(rs) == 4);  /* root + id + data group + val */

    /* Verify element names */
    const carquet_schema_node_t* root = carquet_schema_get_element(rs, 0);
    assert(strcmp(carquet_schema_node_name(root), "schema") == 0);

    /* Verify group is preserved */
    const carquet_schema_node_t* data_group = carquet_schema_get_element(rs, 2);
    assert(!carquet_schema_node_is_leaf(data_group));
    assert(strcmp(carquet_schema_node_name(data_group), "data") == 0);
    assert(carquet_schema_node_repetition(data_group) == CARQUET_REPETITION_OPTIONAL);

    /* Verify levels preserved through serialization */
    assert(carquet_schema_max_def_level(rs, 0) == 0);  /* id */
    assert(carquet_schema_max_def_level(rs, 1) == 1);  /* val under optional data */

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("schema_metadata_roundtrip");
    return 0;
}

/* ============================================================================
 * Test 9: Validate invalid schema operations
 * ============================================================================
 */
static int test_schema_validation(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* Add a leaf column */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    /* Try to add a child under a leaf (should fail) */
    carquet_status_t status = carquet_schema_add_column(schema, "bad", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 1);
    assert(status == CARQUET_ERROR_INVALID_ARGUMENT);

    /* Try to add a group under a leaf (should fail) */
    int32_t bad_group = carquet_schema_add_group(schema, "bad_group",
        CARQUET_REPETITION_REQUIRED, 1);
    assert(bad_group == -1);

    /* Try invalid parent index */
    status = carquet_schema_add_column(schema, "bad2", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 999);
    assert(status == CARQUET_ERROR_INVALID_ARGUMENT);

    carquet_schema_free(schema);
    TEST_PASS("schema_validation");
    return 0;
}

/* ============================================================================
 * Test 10: Write/read with compression and nested schema
 * ============================================================================
 */
static int test_roundtrip_nested_compressed(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_compressed");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    assert(carquet_schema_add_column(schema, "key", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t grp = carquet_schema_add_group(schema, "payload",
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(grp >= 0);
    assert(carquet_schema_add_column(schema, "amount", CARQUET_PHYSICAL_DOUBLE,
        NULL, CARQUET_REPETITION_REQUIRED, 0, grp) == CARQUET_OK);

    /* Write with Snappy compression */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_SNAPPY;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("  Writer failed: %s\n", err.message);
        TEST_FAIL("roundtrip_nested_compressed", "writer creation failed");
    }

    const int N = 1000;
    int32_t keys[1000];
    double amounts[1000];
    int16_t amount_def[1000];

    for (int i = 0; i < N; i++) {
        keys[i] = i;
        amounts[i] = (double)i * 0.5;
        amount_def[i] = (i % 3 == 0) ? 0 : 1;  /* every 3rd row: payload is null */
    }

    assert(carquet_writer_write_batch(writer, 0, keys, N, NULL, NULL) == CARQUET_OK);
    assert(carquet_writer_write_batch(writer, 1, amounts, N, amount_def, NULL) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("roundtrip_nested_compressed", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == N);

    /* Read keys */
    carquet_column_reader_t* key_col = carquet_reader_get_column(reader, 0, 0, &err);
    assert(key_col != NULL);
    int32_t* read_keys = malloc(N * sizeof(int32_t));
    int64_t count = carquet_column_read_batch(key_col, read_keys, N, NULL, NULL);
    assert(count == N);
    for (int i = 0; i < N; i++) {
        assert(read_keys[i] == i);
    }
    free(read_keys);
    carquet_column_reader_free(key_col);

    /* Read amounts with def levels */
    carquet_column_reader_t* amt_col = carquet_reader_get_column(reader, 0, 1, &err);
    assert(amt_col != NULL);
    double* read_amounts = malloc(N * sizeof(double));
    int16_t* read_def = malloc(N * sizeof(int16_t));
    count = carquet_column_read_batch(amt_col, read_amounts, N, read_def, NULL);
    assert(count == N);

    /* Verify def levels */
    for (int i = 0; i < N; i++) {
        int16_t expected_def = (i % 3 == 0) ? 0 : 1;
        assert(read_def[i] == expected_def);
    }

    free(read_amounts);
    free(read_def);
    carquet_column_reader_free(amt_col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("roundtrip_nested_compressed");
    return 0;
}

/* ============================================================================
 * Test 11: carquet_count_rows helper
 * ============================================================================
 */
static int test_count_rows(void) {
    /* NULL rep_levels = flat column, num_values == num_rows */
    assert(carquet_count_rows(NULL, 5) == 5);
    assert(carquet_count_rows(NULL, 0) == 0);

    /* Repeated column: 3 rows with variable-length lists
     * Row 0: [10, 20, 30]  rep=[0, 1, 1]
     * Row 1: [40]          rep=[0]
     * Row 2: [50, 60]      rep=[0, 1]
     */
    int16_t rep[] = {0, 1, 1, 0, 0, 1};
    assert(carquet_count_rows(rep, 6) == 3);

    /* Single row */
    int16_t rep_single[] = {0, 1, 1, 1, 1};
    assert(carquet_count_rows(rep_single, 5) == 1);

    /* All separate rows (no repetition) */
    int16_t rep_flat[] = {0, 0, 0, 0};
    assert(carquet_count_rows(rep_flat, 4) == 4);

    TEST_PASS("count_rows");
    return 0;
}

/* ============================================================================
 * Test 12: carquet_list_offsets helper
 * ============================================================================
 */
static int test_list_offsets(void) {
    /* 3 lists:
     * List 0: values[0..2] = [10, 20, 30]
     * List 1: values[3]    = [40]
     * List 2: values[4..5] = [50, 60]
     */
    int16_t rep[] = {0, 1, 1, 0, 0, 1};
    int64_t offsets[10];

    int64_t num_lists = carquet_list_offsets(rep, 6, 1, offsets, 10);
    assert(num_lists == 3);
    assert(offsets[0] == 0);  /* list 0 starts at index 0 */
    assert(offsets[1] == 3);  /* list 1 starts at index 3 */
    assert(offsets[2] == 4);  /* list 2 starts at index 4 */
    assert(offsets[3] == 6);  /* sentinel: total values */

    /* Verify list lengths */
    assert(offsets[1] - offsets[0] == 3);  /* list 0 has 3 elements */
    assert(offsets[2] - offsets[1] == 1);  /* list 1 has 1 element */
    assert(offsets[3] - offsets[2] == 2);  /* list 2 has 2 elements */

    TEST_PASS("list_offsets");
    return 0;
}

/* ============================================================================
 * Test 13: LIST schema helper
 * ============================================================================
 *
 * carquet_schema_add_list creates:
 *   <name> (OPTIONAL/REQUIRED, LIST) {
 *     list (REPEATED) {
 *       element (OPTIONAL, <type>)
 *     }
 *   }
 *
 * Leaf "element" levels:
 *   For OPTIONAL list: def = optional(list_container)=1 + repeated(list)=1 + optional(element)=1 = 3
 *   For REQUIRED list: def = repeated(list)=1 + optional(element)=1 = 2
 *   rep = repeated(list) = 1
 */
static int test_list_schema_helper(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t list_idx = carquet_schema_add_list(schema, "scores",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    assert(list_idx >= 0);

    /* Should have: root + id + scores(group) + list(group) + element(leaf) = 5 elements */
    assert(carquet_schema_num_elements(schema) == 5);
    /* 2 leaf columns: id, element */
    assert(carquet_schema_num_columns(schema) == 2);

    /* element leaf: def = optional(scores)=1 + repeated(list)=1 + optional(element)=1 = 3, rep = 1 */
    assert(carquet_schema_max_def_level(schema, 1) == 3);
    assert(carquet_schema_max_rep_level(schema, 1) == 1);

    /* Verify path: ["scores", "list", "element"] */
    const char* path[8];
    int32_t depth = carquet_schema_column_path(schema, 1, path, 8);
    assert(depth == 3);
    assert(strcmp(path[0], "scores") == 0);
    assert(strcmp(path[1], "list") == 0);
    assert(strcmp(path[2], "element") == 0);

    /* Verify LIST logical type on outer group */
    const carquet_schema_node_t* scores_node = carquet_schema_get_element(schema, list_idx);
    assert(!carquet_schema_node_is_leaf(scores_node));
    const carquet_logical_type_t* lt = carquet_schema_node_logical_type(scores_node);
    assert(lt != NULL);
    assert(lt->id == CARQUET_LOGICAL_LIST);

    carquet_schema_free(schema);
    TEST_PASS("list_schema_helper");
    return 0;
}

/* ============================================================================
 * Test 14: MAP schema helper
 * ============================================================================
 *
 * carquet_schema_add_map creates:
 *   <name> (MAP) {
 *     key_value (REPEATED) {
 *       key (REQUIRED, <key_type>)
 *       value (OPTIONAL, <value_type>)
 *     }
 *   }
 */
static int test_map_schema_helper(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    int32_t map_idx = carquet_schema_add_map(schema, "attributes",
        CARQUET_PHYSICAL_BYTE_ARRAY, NULL, 0,   /* key: string */
        CARQUET_PHYSICAL_INT32, NULL, 0,          /* value: int32 */
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(map_idx >= 0);

    /* root + attributes(group) + key_value(group) + key(leaf) + value(leaf) = 5 */
    assert(carquet_schema_num_elements(schema) == 5);
    assert(carquet_schema_num_columns(schema) == 2);  /* key, value */

    /* key: required under repeated key_value under optional attributes
     * def = optional(attributes)=1 + repeated(key_value)=1 = 2, rep = 1 */
    assert(carquet_schema_max_def_level(schema, 0) == 2);
    assert(carquet_schema_max_rep_level(schema, 0) == 1);

    /* value: optional under repeated key_value under optional attributes
     * def = optional(attributes)=1 + repeated(key_value)=1 + optional(value)=1 = 3, rep = 1 */
    assert(carquet_schema_max_def_level(schema, 1) == 3);
    assert(carquet_schema_max_rep_level(schema, 1) == 1);

    /* Verify MAP logical type */
    const carquet_schema_node_t* map_node = carquet_schema_get_element(schema, map_idx);
    const carquet_logical_type_t* lt = carquet_schema_node_logical_type(map_node);
    assert(lt != NULL);
    assert(lt->id == CARQUET_LOGICAL_MAP);

    /* Verify key path: ["attributes", "key_value", "key"] */
    const char* path[8];
    int32_t depth = carquet_schema_column_path(schema, 0, path, 8);
    assert(depth == 3);
    assert(strcmp(path[0], "attributes") == 0);
    assert(strcmp(path[1], "key_value") == 0);
    assert(strcmp(path[2], "key") == 0);

    carquet_schema_free(schema);
    TEST_PASS("map_schema_helper");
    return 0;
}

/* ============================================================================
 * Test 15: Write/read repeated column with correct row count
 * ============================================================================
 *
 * This tests the row count fix: when column 0 is a repeated field,
 * the writer must count rows by rep_level==0, not by num_values.
 *
 * Schema:
 *   schema (root)
 *   └── items (repeated group)
 *       └── value (required, INT32)  -> def=1, rep=1
 *
 * Data (3 logical rows):
 *   Row 0: items = [10, 20, 30]
 *   Row 1: items = [40]
 *   Row 2: items = [50, 60]
 */
static int test_repeated_column_row_count(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_repeated");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    int32_t items = carquet_schema_add_group(schema, "items",
        CARQUET_REPETITION_REPEATED, 0);
    assert(items >= 0);
    assert(carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, items) == CARQUET_OK);

    /* value: def=1 (repeated items), rep=1 */
    assert(carquet_schema_max_def_level(schema, 0) == 1);
    assert(carquet_schema_max_rep_level(schema, 0) == 1);

    /* Write */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("  Writer failed: %s\n", err.message);
        TEST_FAIL("repeated_column_row_count", "writer creation failed");
    }

    /* 6 values, 3 logical rows */
    int32_t values[] = {10, 20, 30, 40, 50, 60};
    int16_t def_levels[] = {1, 1, 1, 1, 1, 1};  /* all present */
    int16_t rep_levels[] = {0, 1, 1, 0, 0, 1};   /* row boundaries */

    assert(carquet_writer_write_batch(writer, 0, values, 6, def_levels, rep_levels) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back and verify row count */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("repeated_column_row_count", "reader open failed");
    }

    /* KEY ASSERTION: metadata must say 3 rows, not 6 */
    int64_t num_rows = carquet_reader_num_rows(reader);
    printf("  Repeated column: 6 values, %lld rows in metadata\n", (long long)num_rows);
    assert(num_rows == 3);

    /* Read the values back with rep levels */
    carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, 0, &err);
    assert(col != NULL);

    int32_t read_vals[6];
    int16_t read_def[6];
    int16_t read_rep[6];
    int64_t count = carquet_column_read_batch(col, read_vals, 6, read_def, read_rep);
    assert(count == 6);

    /* Verify values */
    for (int i = 0; i < 6; i++) {
        assert(read_vals[i] == values[i]);
    }

    /* Verify rep levels preserved */
    assert(read_rep[0] == 0);  /* new row */
    assert(read_rep[1] == 1);  /* continuation */
    assert(read_rep[2] == 1);
    assert(read_rep[3] == 0);  /* new row */
    assert(read_rep[4] == 0);  /* new row */
    assert(read_rep[5] == 1);  /* continuation */

    /* Use carquet_count_rows to verify */
    assert(carquet_count_rows(read_rep, count) == 3);

    /* Use carquet_list_offsets to reconstruct list boundaries */
    int64_t offsets[10];
    int64_t num_lists = carquet_list_offsets(read_rep, count, 1, offsets, 10);
    assert(num_lists == 3);
    assert(offsets[0] == 0);
    assert(offsets[1] == 3);
    assert(offsets[2] == 4);
    assert(offsets[3] == 6);

    carquet_column_reader_free(col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("repeated_column_row_count");
    return 0;
}

/* ============================================================================
 * Test 16: LIST schema helper write/read round-trip
 * ============================================================================
 *
 * Write a file with list<int32> using the LIST helper, read it back,
 * and reconstruct list boundaries.
 */
static int test_list_roundtrip(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_list_rt");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* id column + list<int32> column */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t list_idx = carquet_schema_add_list(schema, "tags",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
    assert(list_idx >= 0);

    /* element: def=3 (optional tags + repeated list + optional element), rep=1 */
    assert(carquet_schema_max_def_level(schema, 1) == 3);
    assert(carquet_schema_max_rep_level(schema, 1) == 1);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("  Writer failed: %s\n", err.message);
        TEST_FAIL("list_roundtrip", "writer creation failed");
    }

    /* 4 rows */
    int32_t ids[] = {1, 2, 3, 4};
    assert(carquet_writer_write_batch(writer, 0, ids, 4, NULL, NULL) == CARQUET_OK);

    /* List column (leaf "element"):
     * max_def=3, max_rep=1
     *
     * Row 0: tags = [100, 200]         -> def=[3,3], rep=[0,1]
     * Row 1: tags = NULL               -> def=[0],   rep=[0]
     * Row 2: tags = [300]              -> def=[3],   rep=[0]
     * Row 3: tags = [400, 500, 600]    -> def=[3,3,3], rep=[0,1,1]
     *
     * Total: 8 entries in def/rep arrays, 6 non-null values
     */
    /* For a list<int32> with max_def=3:
     * def=0: tags group is null (entire list absent)
     * def=1: tags group is present but list is empty (repeated group absent)
     * def=2: list element exists but value is null
     * def=3: list element with non-null value
     *
     * Row 0: [100, 200]       -> entries: (def=3,rep=0), (def=3,rep=1)
     * Row 1: NULL              -> entries: (def=0,rep=0)
     * Row 2: [300]             -> entries: (def=3,rep=0)
     * Row 3: [400, 500, 600]  -> entries: (def=3,rep=0), (def=3,rep=1), (def=3,rep=1)
     *
     * Total: 7 entries, 6 non-null values
     */
    int32_t tag_vals2[] = {100, 200, 300, 400, 500, 600};
    int16_t tag_def2[] = {3, 3, 0, 3, 3, 3, 3};
    int16_t tag_rep2[] = {0, 1, 0, 0, 0, 1, 1};

    assert(carquet_writer_write_batch(writer, 1, tag_vals2, 7, tag_def2, tag_rep2) == CARQUET_OK);
    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("list_roundtrip", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 4);

    /* Read tag column with levels */
    carquet_column_reader_t* tag_col = carquet_reader_get_column(reader, 0, 1, &err);
    assert(tag_col != NULL);

    int32_t read_vals[10];
    int16_t read_def2[10];
    int16_t read_rep2[10];
    int64_t count = carquet_column_read_batch(tag_col, read_vals, 10, read_def2, read_rep2);
    assert(count == 7);

    /* Verify def levels */
    assert(read_def2[0] == 3);
    assert(read_def2[1] == 3);
    assert(read_def2[2] == 0);  /* NULL list */
    assert(read_def2[3] == 3);
    assert(read_def2[4] == 3);
    assert(read_def2[5] == 3);
    assert(read_def2[6] == 3);

    /* Verify rep levels */
    assert(read_rep2[0] == 0);
    assert(read_rep2[1] == 1);
    assert(read_rep2[2] == 0);
    assert(read_rep2[3] == 0);
    assert(read_rep2[4] == 0);
    assert(read_rep2[5] == 1);
    assert(read_rep2[6] == 1);

    /* Use list_offsets to reconstruct */
    int64_t offsets[10];
    int64_t num_lists = carquet_list_offsets(read_rep2, count, 1, offsets, 10);
    assert(num_lists == 4);  /* 4 rows */
    assert(offsets[0] == 0);
    assert(offsets[1] == 2);  /* row 1 starts at index 2 (the NULL entry) */
    assert(offsets[2] == 3);  /* row 2 starts at index 3 */
    assert(offsets[3] == 4);  /* row 3 starts at index 4 */
    assert(offsets[4] == 7);  /* sentinel */

    carquet_column_reader_free(tag_col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("list_roundtrip");
    return 0;
}

/* ============================================================================
 * Test 17: MAP write/read round-trip
 * ============================================================================
 */
static int test_map_roundtrip(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_map_rt");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* id + map<int32, int32> */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    int32_t map_idx = carquet_schema_add_map(schema, "props",
        CARQUET_PHYSICAL_INT32, NULL, 0,
        CARQUET_PHYSICAL_INT32, NULL, 0,
        CARQUET_REPETITION_OPTIONAL, 0);
    assert(map_idx >= 0);

    /* 3 leaf columns: id, key, value */
    assert(carquet_schema_num_columns(schema) == 3);

    /* key: def=2 (optional props + repeated key_value), rep=1 */
    assert(carquet_schema_max_def_level(schema, 1) == 2);
    assert(carquet_schema_max_rep_level(schema, 1) == 1);

    /* value: def=3 (optional props + repeated key_value + optional value), rep=1 */
    assert(carquet_schema_max_def_level(schema, 2) == 3);
    assert(carquet_schema_max_rep_level(schema, 2) == 1);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    if (!writer) {
        carquet_schema_free(schema);
        printf("  Writer failed: %s\n", err.message);
        TEST_FAIL("map_roundtrip", "writer creation failed");
    }

    /* 3 rows */
    int32_t ids[] = {1, 2, 3};
    assert(carquet_writer_write_batch(writer, 0, ids, 3, NULL, NULL) == CARQUET_OK);

    /* Map data:
     * Row 0: {1: 10, 2: 20}     -> key entries: (def=2,rep=0), (def=2,rep=1)
     * Row 1: NULL map            -> key entries: (def=0,rep=0)
     * Row 2: {3: 30}            -> key entries: (def=2,rep=0)
     * Total key entries: 4, non-null keys: 3
     */
    int32_t keys[] = {1, 2, 3};
    int16_t key_def[] = {2, 2, 0, 2};
    int16_t key_rep[] = {0, 1, 0, 0};
    assert(carquet_writer_write_batch(writer, 1, keys, 4, key_def, key_rep) == CARQUET_OK);

    /* Values:
     * Row 0: {1:10, 2:20}  -> (def=3,rep=0), (def=3,rep=1)
     * Row 1: NULL map       -> (def=0,rep=0)
     * Row 2: {3:30}         -> (def=3,rep=0)
     */
    int32_t vals[] = {10, 20, 30};
    int16_t val_def[] = {3, 3, 0, 3};
    int16_t val_rep[] = {0, 1, 0, 0};
    assert(carquet_writer_write_batch(writer, 2, vals, 4, val_def, val_rep) == CARQUET_OK);

    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Read back */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    if (!reader) {
        carquet_schema_free(schema);
        remove(test_file);
        printf("  Reader failed: %s\n", err.message);
        TEST_FAIL("map_roundtrip", "reader open failed");
    }

    assert(carquet_reader_num_rows(reader) == 3);

    /* Read key column and verify */
    carquet_column_reader_t* key_col = carquet_reader_get_column(reader, 0, 1, &err);
    assert(key_col != NULL);
    int32_t read_keys[10];
    int16_t read_key_def[10];
    int16_t read_key_rep[10];
    int64_t count = carquet_column_read_batch(key_col, read_keys, 10, read_key_def, read_key_rep);
    assert(count == 4);

    /* Verify def and rep levels */
    assert(read_key_def[0] == 2 && read_key_rep[0] == 0);
    assert(read_key_def[1] == 2 && read_key_rep[1] == 1);
    assert(read_key_def[2] == 0 && read_key_rep[2] == 0);
    assert(read_key_def[3] == 2 && read_key_rep[3] == 0);

    carquet_column_reader_free(key_col);
    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("map_roundtrip");
    return 0;
}

/* ============================================================================
 * Test 18: Repeated field with multiple write batches
 * ============================================================================
 *
 * Verifies that row counting works correctly across multiple write_batch calls.
 */
static int test_repeated_multi_batch(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_multi_batch");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    int32_t items = carquet_schema_add_group(schema, "items",
        CARQUET_REPETITION_REPEATED, 0);
    assert(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, items) == CARQUET_OK);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer != NULL);

    /* Batch 1: 2 rows with 4 values */
    int32_t vals1[] = {1, 2, 3, 4};
    int16_t def1[] = {1, 1, 1, 1};
    int16_t rep1[] = {0, 1, 0, 1};  /* row 0: [1,2], row 1: [3,4] */
    assert(carquet_writer_write_batch(writer, 0, vals1, 4, def1, rep1) == CARQUET_OK);

    /* Batch 2: 1 row with 3 values */
    int32_t vals2[] = {5, 6, 7};
    int16_t def2[] = {1, 1, 1};
    int16_t rep2[] = {0, 1, 1};  /* row 2: [5,6,7] */
    assert(carquet_writer_write_batch(writer, 0, vals2, 3, def2, rep2) == CARQUET_OK);

    assert(carquet_writer_close(writer) == CARQUET_OK);

    /* Verify 3 rows total, not 7 */
    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    assert(reader != NULL);
    assert(carquet_reader_num_rows(reader) == 3);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("repeated_multi_batch");
    return 0;
}

/* ============================================================================
 * Test 19: Mixed flat + repeated columns row count
 * ============================================================================
 *
 * When column 0 is flat and column 1 is repeated, row count
 * should be based on column 0 (flat), which is correct by default.
 */
static int test_mixed_flat_repeated(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    char test_file[512];
    carquet_test_temp_path(test_file, sizeof(test_file), "nested_mixed_rep");

    carquet_schema_t* schema = carquet_schema_create(&err);
    assert(schema != NULL);

    /* Flat column first */
    assert(carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, 0) == CARQUET_OK);

    /* Repeated column second */
    int32_t items = carquet_schema_add_group(schema, "items",
        CARQUET_REPETITION_REPEATED, 0);
    assert(carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT32,
        NULL, CARQUET_REPETITION_REQUIRED, 0, items) == CARQUET_OK);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;

    carquet_writer_t* writer = carquet_writer_create(test_file, schema, &opts, &err);
    assert(writer != NULL);

    /* 3 flat rows */
    int32_t ids[] = {1, 2, 3};
    assert(carquet_writer_write_batch(writer, 0, ids, 3, NULL, NULL) == CARQUET_OK);

    /* 3 rows with varying list lengths (7 total values) */
    int32_t vals[] = {10, 20, 30, 40, 50, 60, 70};
    int16_t def[] = {1, 1, 1, 1, 1, 1, 1};
    int16_t rep[] = {0, 1, 1, 0, 0, 1, 1};  /* [10,20,30], [40], [50,60,70] */
    assert(carquet_writer_write_batch(writer, 1, vals, 7, def, rep) == CARQUET_OK);

    assert(carquet_writer_close(writer) == CARQUET_OK);

    carquet_reader_t* reader = carquet_reader_open(test_file, NULL, &err);
    assert(reader != NULL);
    assert(carquet_reader_num_rows(reader) == 3);

    /* Read both columns */
    carquet_column_reader_t* id_col = carquet_reader_get_column(reader, 0, 0, &err);
    assert(id_col != NULL);
    int32_t read_ids[3];
    assert(carquet_column_read_batch(id_col, read_ids, 3, NULL, NULL) == 3);
    assert(read_ids[0] == 1 && read_ids[1] == 2 && read_ids[2] == 3);
    carquet_column_reader_free(id_col);

    carquet_column_reader_t* v_col = carquet_reader_get_column(reader, 0, 1, &err);
    assert(v_col != NULL);
    int32_t read_vals[10];
    int16_t read_rep[10];
    int64_t count = carquet_column_read_batch(v_col, read_vals, 10, NULL, read_rep);
    assert(count == 7);
    assert(carquet_count_rows(read_rep, count) == 3);
    carquet_column_reader_free(v_col);

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(test_file);

    TEST_PASS("mixed_flat_repeated");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */
int main(void) {
    int failures = 0;

    printf("=== Nested Schema Tests ===\n\n");

    /* Schema level tests */
    failures += test_accumulated_levels();
    failures += test_column_paths();
    failures += test_deep_nesting_levels();
    failures += test_multiple_repeated_ancestors();

    /* Round-trip tests (optional groups) */
    failures += test_roundtrip_optional_group();
    failures += test_roundtrip_double_optional();
    failures += test_mixed_schema_structure();
    failures += test_schema_metadata_roundtrip();
    failures += test_schema_validation();
    failures += test_roundtrip_nested_compressed();

    /* Helper function tests */
    failures += test_count_rows();
    failures += test_list_offsets();

    /* LIST/MAP schema helpers */
    failures += test_list_schema_helper();
    failures += test_map_schema_helper();

    /* Repeated field round-trip tests */
    failures += test_repeated_column_row_count();
    failures += test_list_roundtrip();
    failures += test_map_roundtrip();
    failures += test_repeated_multi_batch();
    failures += test_mixed_flat_repeated();

    printf("\n");
    if (failures == 0) {
        printf("All nested tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
