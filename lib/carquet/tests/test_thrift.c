/**
 * @file test_thrift.c
 * @brief Tests for Thrift encoding/decoding
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <carquet/carquet.h>
#include "thrift/thrift_decode.h"
#include "thrift/thrift_encode.h"
#include "thrift/parquet_types.h"
#include "core/buffer.h"
#include "core/arena.h"

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

static int test_thrift_varint_roundtrip(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    /* Write various integers */
    thrift_write_i32(&enc, 0);
    thrift_write_i32(&enc, 1);
    thrift_write_i32(&enc, -1);
    thrift_write_i32(&enc, 127);
    thrift_write_i32(&enc, 128);
    thrift_write_i32(&enc, 12345);
    thrift_write_i32(&enc, -12345);
    thrift_write_i32(&enc, 0x7FFFFFFF);
    thrift_write_i32(&enc, (int32_t)0x80000000);

    /* Read them back */
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, carquet_buffer_data_const(&buf), carquet_buffer_size(&buf));

    assert(thrift_read_i32(&dec) == 0);
    assert(thrift_read_i32(&dec) == 1);
    assert(thrift_read_i32(&dec) == -1);
    assert(thrift_read_i32(&dec) == 127);
    assert(thrift_read_i32(&dec) == 128);
    assert(thrift_read_i32(&dec) == 12345);
    assert(thrift_read_i32(&dec) == -12345);
    assert(thrift_read_i32(&dec) == 0x7FFFFFFF);
    assert(thrift_read_i32(&dec) == (int32_t)0x80000000);

    assert(!thrift_decoder_has_error(&dec));

    carquet_buffer_destroy(&buf);
    TEST_PASS("thrift_varint_roundtrip");
    return 0;
}

static int test_thrift_string_roundtrip(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    thrift_write_string(&enc, "Hello, World!");
    thrift_write_string(&enc, "");
    thrift_write_string(&enc, "A longer string with more characters");

    thrift_decoder_t dec;
    thrift_decoder_init(&dec, carquet_buffer_data_const(&buf), carquet_buffer_size(&buf));

    int32_t len;
    const uint8_t* data = NULL;
    (void)data;

    data = thrift_read_binary(&dec, &len);
    assert(len == 13);
    assert(memcmp(data, "Hello, World!", 13) == 0);

    data = thrift_read_binary(&dec, &len);
    assert(len == 0);

    data = thrift_read_binary(&dec, &len);
    assert(len == 36);

    assert(!thrift_decoder_has_error(&dec));

    carquet_buffer_destroy(&buf);
    TEST_PASS("thrift_string_roundtrip");
    return 0;
}

static int test_thrift_struct(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    /* Write a simple struct */
    thrift_write_struct_begin(&enc);
    THRIFT_WRITE_FIELD_I32(&enc, 1, 42);
    THRIFT_WRITE_FIELD_STRING(&enc, 2, "test");
    THRIFT_WRITE_FIELD_BOOL(&enc, 3, true);
    thrift_write_struct_end(&enc);

    /* Read it back */
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, carquet_buffer_data_const(&buf), carquet_buffer_size(&buf));

    thrift_read_struct_begin(&dec);

    thrift_type_t type;
    (void)type;
    int16_t field_id;
    (void)field_id;

    assert(thrift_read_field_begin(&dec, &type, &field_id));
    assert(field_id == 1);
    assert(type == THRIFT_TYPE_I32);
    assert(thrift_read_i32(&dec) == 42);

    assert(thrift_read_field_begin(&dec, &type, &field_id));
    assert(field_id == 2);
    assert(type == THRIFT_TYPE_BINARY);
    int32_t len;
    const uint8_t* data = thrift_read_binary(&dec, &len);
    (void)data;
    assert(len == 4 && memcmp(data, "test", 4) == 0);

    assert(thrift_read_field_begin(&dec, &type, &field_id));
    assert(field_id == 3);
    assert(type == THRIFT_TYPE_TRUE);
    assert(thrift_read_bool(&dec) == true);

    assert(!thrift_read_field_begin(&dec, &type, &field_id));  /* STOP */

    thrift_read_struct_end(&dec);

    assert(!thrift_decoder_has_error(&dec));

    carquet_buffer_destroy(&buf);
    TEST_PASS("thrift_struct");
    return 0;
}

static int test_thrift_list(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    /* Write a list of integers */
    thrift_write_list_begin(&enc, THRIFT_TYPE_I32, 5);
    thrift_write_i32(&enc, 1);
    thrift_write_i32(&enc, 2);
    thrift_write_i32(&enc, 3);
    thrift_write_i32(&enc, 4);
    thrift_write_i32(&enc, 5);

    /* Read it back */
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, carquet_buffer_data_const(&buf), carquet_buffer_size(&buf));

    thrift_type_t elem_type;
    int32_t count;
    thrift_read_list_begin(&dec, &elem_type, &count);

    assert(elem_type == THRIFT_TYPE_I32);
    assert(count == 5);

    for (int i = 1; i <= 5; i++) {
        assert(thrift_read_i32(&dec) == i);
    }

    assert(!thrift_decoder_has_error(&dec));

    carquet_buffer_destroy(&buf);
    TEST_PASS("thrift_list");
    return 0;
}

static int test_logical_type_converted_type_compat(void) {
    parquet_schema_element_t schema[20];
    memset(schema, 0, sizeof(schema));

    schema[0].name = "schema";
    schema[0].num_children = 19;

    schema[1].name = "str";
    schema[1].has_type = true;
    schema[1].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[1].has_logical_type = true;
    schema[1].logical_type.id = CARQUET_LOGICAL_STRING;

    schema[2].name = "date";
    schema[2].has_type = true;
    schema[2].type = CARQUET_PHYSICAL_INT32;
    schema[2].has_logical_type = true;
    schema[2].logical_type.id = CARQUET_LOGICAL_DATE;

    schema[3].name = "decimal";
    schema[3].has_type = true;
    schema[3].type = CARQUET_PHYSICAL_INT64;
    schema[3].has_logical_type = true;
    schema[3].logical_type.id = CARQUET_LOGICAL_DECIMAL;
    schema[3].logical_type.params.decimal.scale = 2;
    schema[3].logical_type.params.decimal.precision = 9;

    schema[4].name = "time_millis";
    schema[4].has_type = true;
    schema[4].type = CARQUET_PHYSICAL_INT32;
    schema[4].has_logical_type = true;
    schema[4].logical_type.id = CARQUET_LOGICAL_TIME;
    schema[4].logical_type.params.time.is_adjusted_to_utc = false;
    schema[4].logical_type.params.time.unit = CARQUET_TIME_UNIT_MILLIS;

    schema[5].name = "time_nanos";
    schema[5].has_type = true;
    schema[5].type = CARQUET_PHYSICAL_INT64;
    schema[5].has_logical_type = true;
    schema[5].logical_type.id = CARQUET_LOGICAL_TIME;
    schema[5].logical_type.params.time.unit = CARQUET_TIME_UNIT_NANOS;

    schema[6].name = "ts_micros";
    schema[6].has_type = true;
    schema[6].type = CARQUET_PHYSICAL_INT64;
    schema[6].has_logical_type = true;
    schema[6].logical_type.id = CARQUET_LOGICAL_TIMESTAMP;
    schema[6].logical_type.params.timestamp.is_adjusted_to_utc = false;
    schema[6].logical_type.params.timestamp.unit = CARQUET_TIME_UNIT_MICROS;

    schema[7].name = "ts_nanos";
    schema[7].has_type = true;
    schema[7].type = CARQUET_PHYSICAL_INT64;
    schema[7].has_logical_type = true;
    schema[7].logical_type.id = CARQUET_LOGICAL_TIMESTAMP;
    schema[7].logical_type.params.timestamp.unit = CARQUET_TIME_UNIT_NANOS;

    schema[8].name = "int16";
    schema[8].has_type = true;
    schema[8].type = CARQUET_PHYSICAL_INT32;
    schema[8].has_logical_type = true;
    schema[8].logical_type.id = CARQUET_LOGICAL_INTEGER;
    schema[8].logical_type.params.integer.bit_width = 16;
    schema[8].logical_type.params.integer.is_signed = true;

    schema[9].name = "uint64";
    schema[9].has_type = true;
    schema[9].type = CARQUET_PHYSICAL_INT64;
    schema[9].has_logical_type = true;
    schema[9].logical_type.id = CARQUET_LOGICAL_INTEGER;
    schema[9].logical_type.params.integer.bit_width = 64;
    schema[9].logical_type.params.integer.is_signed = false;

    schema[10].name = "json";
    schema[10].has_type = true;
    schema[10].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[10].has_logical_type = true;
    schema[10].logical_type.id = CARQUET_LOGICAL_JSON;

    schema[11].name = "bson";
    schema[11].has_type = true;
    schema[11].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[11].has_logical_type = true;
    schema[11].logical_type.id = CARQUET_LOGICAL_BSON;

    schema[12].name = "enum";
    schema[12].has_type = true;
    schema[12].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[12].has_logical_type = true;
    schema[12].logical_type.id = CARQUET_LOGICAL_ENUM;

    schema[13].name = "list";
    schema[13].has_logical_type = true;
    schema[13].logical_type.id = CARQUET_LOGICAL_LIST;

    schema[14].name = "map";
    schema[14].has_logical_type = true;
    schema[14].logical_type.id = CARQUET_LOGICAL_MAP;

    schema[15].name = "legacy_int8";
    schema[15].has_type = true;
    schema[15].type = CARQUET_PHYSICAL_INT32;
    schema[15].has_converted_type = true;
    schema[15].converted_type = CARQUET_CONVERTED_INT_8;

    schema[16].name = "legacy_map_key_value";
    schema[16].has_converted_type = true;
    schema[16].converted_type = CARQUET_CONVERTED_MAP_KEY_VALUE;

    schema[17].name = "variant";
    schema[17].has_logical_type = true;
    schema[17].logical_type.id = CARQUET_LOGICAL_VARIANT;
    schema[17].logical_type.params.variant.specification_version = 1;

    schema[18].name = "geometry";
    schema[18].has_type = true;
    schema[18].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[18].has_logical_type = true;
    schema[18].logical_type.id = CARQUET_LOGICAL_GEOMETRY;

    schema[19].name = "geography";
    schema[19].has_type = true;
    schema[19].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    schema[19].has_logical_type = true;
    schema[19].logical_type.id = CARQUET_LOGICAL_GEOGRAPHY;
    schema[19].logical_type.params.geography.has_algorithm = true;
    schema[19].logical_type.params.geography.algorithm = CARQUET_GEOSPATIAL_EDGE_VINCENTY;

    parquet_file_metadata_t metadata;
    memset(&metadata, 0, sizeof(metadata));
    metadata.version = 1;
    metadata.schema = schema;
    metadata.num_schema_elements = 20;
    metadata.num_rows = 0;
    metadata.created_by = "test";
    metadata.num_column_orders = 18;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    assert(parquet_write_file_metadata(&metadata, &buf, NULL) == CARQUET_OK);

    carquet_arena_t arena;
    assert(carquet_arena_init_size(&arena, 4096) == CARQUET_OK);

    parquet_file_metadata_t parsed;
    memset(&parsed, 0, sizeof(parsed));
    assert(parquet_parse_file_metadata(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        &arena, &parsed, NULL) == CARQUET_OK);

    assert(parsed.schema[1].has_converted_type);
    assert(parsed.schema[1].converted_type == CARQUET_CONVERTED_UTF8);
    assert(parsed.schema[2].converted_type == CARQUET_CONVERTED_DATE);
    assert(parsed.schema[3].converted_type == CARQUET_CONVERTED_DECIMAL);
    assert(parsed.schema[3].scale == 2);
    assert(parsed.schema[3].precision == 9);
    assert(parsed.schema[4].converted_type == CARQUET_CONVERTED_TIME_MILLIS);
    assert(!parsed.schema[5].has_converted_type);
    assert(parsed.schema[6].converted_type == CARQUET_CONVERTED_TIMESTAMP_MICROS);
    assert(!parsed.schema[7].has_converted_type);
    assert(parsed.schema[8].converted_type == CARQUET_CONVERTED_INT_16);
    assert(parsed.schema[9].converted_type == CARQUET_CONVERTED_UINT_64);
    assert(parsed.schema[10].converted_type == CARQUET_CONVERTED_JSON);
    assert(parsed.schema[11].converted_type == CARQUET_CONVERTED_BSON);
    assert(parsed.schema[12].converted_type == CARQUET_CONVERTED_ENUM);
    assert(parsed.schema[13].converted_type == CARQUET_CONVERTED_LIST);
    assert(parsed.schema[14].converted_type == CARQUET_CONVERTED_MAP);

    assert(parsed.schema[15].has_logical_type);
    assert(parsed.schema[15].logical_type.id == CARQUET_LOGICAL_INTEGER);
    assert(parsed.schema[15].logical_type.params.integer.is_signed);
    assert(parsed.schema[15].logical_type.params.integer.bit_width == 8);
    assert(parsed.schema[16].has_logical_type);
    assert(parsed.schema[16].logical_type.id == CARQUET_LOGICAL_MAP);
    assert(parsed.schema[17].has_logical_type);
    assert(parsed.schema[17].logical_type.id == CARQUET_LOGICAL_VARIANT);
    assert(parsed.schema[17].logical_type.params.variant.specification_version == 1);
    assert(parsed.schema[18].has_logical_type);
    assert(parsed.schema[18].logical_type.id == CARQUET_LOGICAL_GEOMETRY);
    assert(!parsed.schema[18].has_converted_type);
    assert(parsed.schema[19].has_logical_type);
    assert(parsed.schema[19].logical_type.id == CARQUET_LOGICAL_GEOGRAPHY);
    assert(parsed.schema[19].logical_type.params.geography.has_algorithm);
    assert(parsed.schema[19].logical_type.params.geography.algorithm == CARQUET_GEOSPATIAL_EDGE_VINCENTY);
    assert(!parsed.schema[19].has_converted_type);
    assert(parsed.num_column_orders == 18);

    carquet_arena_destroy(&arena);
    carquet_buffer_destroy(&buf);

    TEST_PASS("logical_type_converted_type_compat");
    return 0;
}

/*
 * A foreign file may legitimately carry an empty (zero-length) BYTE_ARRAY
 * min_value/max_value when the column's extreme value is the empty string.
 * The parser must record presence from the field itself, not from its length,
 * so that predicate pushdown stays enabled. carquet's own writer never emits an
 * empty min/max, so this is hand-encoded to reproduce the foreign-file case.
 */
static int test_statistics_empty_min_value_presence(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    /* FileMetaData */
    thrift_write_struct_begin(&enc);
    thrift_write_field_header(&enc, THRIFT_TYPE_I32, 1);   /* version */
    thrift_write_i32(&enc, 1);
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 2);  /* schema (empty) */
    thrift_write_list_begin(&enc, THRIFT_TYPE_STRUCT, 0);
    thrift_write_field_header(&enc, THRIFT_TYPE_I64, 3);   /* num_rows */
    thrift_write_i64(&enc, 3);
    thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 4);  /* row_groups */
    thrift_write_list_begin(&enc, THRIFT_TYPE_STRUCT, 1);
    {
        /* RowGroup */
        thrift_write_struct_begin(&enc);
        thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 1);  /* columns */
        thrift_write_list_begin(&enc, THRIFT_TYPE_STRUCT, 1);
        {
            /* ColumnChunk */
            thrift_write_struct_begin(&enc);
            thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 3);  /* meta_data */
            {
                /* ColumnMetaData */
                thrift_write_struct_begin(&enc);
                thrift_write_field_header(&enc, THRIFT_TYPE_I32, 1);   /* type */
                thrift_write_i32(&enc, (int32_t)CARQUET_PHYSICAL_BYTE_ARRAY);
                thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 2);  /* encodings */
                thrift_write_list_begin(&enc, THRIFT_TYPE_I32, 0);
                thrift_write_field_header(&enc, THRIFT_TYPE_LIST, 3);  /* path_in_schema */
                thrift_write_list_begin(&enc, THRIFT_TYPE_BINARY, 1);
                thrift_write_string(&enc, "s");
                thrift_write_field_header(&enc, THRIFT_TYPE_I32, 4);   /* codec */
                thrift_write_i32(&enc, 0);
                thrift_write_field_header(&enc, THRIFT_TYPE_I64, 5);   /* num_values */
                thrift_write_i64(&enc, 3);
                thrift_write_field_header(&enc, THRIFT_TYPE_I64, 6);   /* total_uncompressed */
                thrift_write_i64(&enc, 0);
                thrift_write_field_header(&enc, THRIFT_TYPE_I64, 7);   /* total_compressed */
                thrift_write_i64(&enc, 0);
                thrift_write_field_header(&enc, THRIFT_TYPE_I64, 9);   /* data_page_offset */
                thrift_write_i64(&enc, 4);
                thrift_write_field_header(&enc, THRIFT_TYPE_STRUCT, 12);  /* statistics */
                {
                    /* Statistics with EMPTY min_value and max_value */
                    thrift_write_struct_begin(&enc);
                    thrift_write_field_header(&enc, THRIFT_TYPE_I64, 3);  /* null_count */
                    thrift_write_i64(&enc, 0);
                    thrift_write_field_header(&enc, THRIFT_TYPE_BINARY, 5);  /* max_value = "" */
                    thrift_write_binary(&enc, (const uint8_t*)"", 0);
                    thrift_write_field_header(&enc, THRIFT_TYPE_BINARY, 6);  /* min_value = "" */
                    thrift_write_binary(&enc, (const uint8_t*)"", 0);
                    thrift_write_field_stop(&enc);
                }
                thrift_write_field_stop(&enc);  /* end ColumnMetaData */
            }
            thrift_write_field_stop(&enc);  /* end ColumnChunk */
        }
        thrift_write_field_header(&enc, THRIFT_TYPE_I64, 3);  /* RowGroup num_rows */
        thrift_write_i64(&enc, 3);
        thrift_write_field_stop(&enc);  /* end RowGroup */
    }
    thrift_write_field_stop(&enc);  /* end FileMetaData */

    carquet_arena_t arena;
    assert(carquet_arena_init_size(&arena, 4096) == CARQUET_OK);

    parquet_file_metadata_t parsed;
    memset(&parsed, 0, sizeof(parsed));
    assert(parquet_parse_file_metadata(
        carquet_buffer_data_const(&buf), carquet_buffer_size(&buf),
        &arena, &parsed, NULL) == CARQUET_OK);

    assert(parsed.num_row_groups == 1);
    assert(parsed.row_groups[0].num_columns == 1);
    const parquet_column_chunk_t* chunk = &parsed.row_groups[0].columns[0];
    assert(chunk->has_metadata);
    assert(chunk->metadata.has_statistics);
    const parquet_statistics_t* st = &chunk->metadata.statistics;

    /* The core assertion: presence is detected from the field, not its length. */
    assert(st->has_min_value);
    assert(st->has_max_value);
    assert(st->min_value_len == 0);
    assert(st->max_value_len == 0);

    carquet_arena_destroy(&arena);
    carquet_buffer_destroy(&buf);

    TEST_PASS("statistics_empty_min_value_presence");
    return 0;
}

int main(void) {
    int failures = 0;

    printf("=== Thrift Tests ===\n\n");

    failures += test_thrift_varint_roundtrip();
    failures += test_thrift_string_roundtrip();
    failures += test_thrift_struct();
    failures += test_thrift_list();
    failures += test_logical_type_converted_type_compat();
    failures += test_statistics_empty_min_value_presence();

    printf("\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
