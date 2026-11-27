/**
 * @file fuzz_thrift.c
 * @brief Fuzz target for carquet Thrift compact protocol decoder
 *
 * Tests the Thrift decoder which parses Parquet file metadata —
 * a critical attack surface for malicious files.
 *
 * Modes: primitives, struct parsing, containers (list/map),
 *        file metadata, page headers, and encoder roundtrips.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <carquet/carquet.h>

/* Internal Thrift headers */
#include "thrift/parquet_types.h"
#include "thrift/thrift_decode.h"
#include "thrift/thrift_encode.h"
#include "core/arena.h"
#include "core/buffer.h"

/**
 * Mode 0: Low-level Thrift primitives
 */
static void fuzz_thrift_primitives(const uint8_t* data, size_t size) {
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, size);

    while (!thrift_decoder_has_error(&dec) && thrift_decoder_remaining(&dec) > 0) {
        size_t remaining = thrift_decoder_remaining(&dec);

        if (remaining >= 1) { (void)thrift_read_byte(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 1) { (void)thrift_read_varint(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 1) { (void)thrift_read_zigzag(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 2) { (void)thrift_read_i16(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 4) { (void)thrift_read_i32(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 8) { (void)thrift_read_i64(&dec); if (thrift_decoder_has_error(&dec)) break; }
        if (remaining >= 8) { (void)thrift_read_double(&dec); if (thrift_decoder_has_error(&dec)) break; }

        { int32_t len = 0; (void)thrift_read_binary(&dec, &len); }
        break;
    }
}

/**
 * Mode 1: Thrift struct parsing with field iteration
 */
static void fuzz_thrift_struct(const uint8_t* data, size_t size) {
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, size);

    thrift_read_struct_begin(&dec);
    if (thrift_decoder_has_error(&dec)) return;

    thrift_type_t type;
    int16_t field_id;
    int count = 0;

    while (!thrift_decoder_has_error(&dec) &&
           thrift_read_field_begin(&dec, &type, &field_id) &&
           count < 100) {
        thrift_skip_field(&dec, type);
        count++;
    }
    thrift_read_struct_end(&dec);
}

/**
 * Mode 2: Thrift containers (list/set/map)
 */
static void fuzz_thrift_containers(const uint8_t* data, size_t size) {
    thrift_decoder_t dec;
    thrift_decoder_init(&dec, data, size);

    /* List */
    thrift_type_t elem_type;
    int32_t count;
    thrift_read_list_begin(&dec, &elem_type, &count);
    if (!thrift_decoder_has_error(&dec)) {
        int max_iter = count > 100 ? 100 : count;
        for (int i = 0; i < max_iter && !thrift_decoder_has_error(&dec); i++)
            thrift_skip(&dec, elem_type);
    }

    /* Map — fresh decoder */
    thrift_decoder_init(&dec, data, size);
    thrift_type_t key_type, value_type;
    thrift_read_map_begin(&dec, &key_type, &value_type, &count);
    if (!thrift_decoder_has_error(&dec)) {
        int max_iter = count > 50 ? 50 : count;
        for (int i = 0; i < max_iter && !thrift_decoder_has_error(&dec); i++) {
            thrift_skip(&dec, key_type);
            thrift_skip(&dec, value_type);
        }
    }
}

/**
 * Mode 3: Full Parquet file metadata parsing
 */
static void fuzz_parquet_metadata(const uint8_t* data, size_t size) {
    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) return;

    parquet_file_metadata_t metadata;
    memset(&metadata, 0, sizeof(metadata));
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_status_t status = parquet_parse_file_metadata(
        data, size, &arena, &metadata, &err);

    if (status == CARQUET_OK) {
        (void)metadata.version;
        (void)metadata.num_rows;
        (void)metadata.num_row_groups;

        for (int32_t i = 0; i < metadata.num_schema_elements && i < 100; i++) {
            if (metadata.schema) {
                (void)metadata.schema[i].name;
                (void)metadata.schema[i].type;
                (void)metadata.schema[i].repetition_type;
                (void)metadata.schema[i].num_children;
            }
        }

        for (int32_t i = 0; i < metadata.num_row_groups && i < 10; i++) {
            if (metadata.row_groups) {
                (void)metadata.row_groups[i].num_rows;
                (void)metadata.row_groups[i].num_columns;
                /* Walk column chunks */
                for (int32_t j = 0; j < metadata.row_groups[i].num_columns && j < 50; j++) {
                    if (metadata.row_groups[i].columns) {
                        (void)metadata.row_groups[i].columns[j].file_offset;
                        (void)metadata.row_groups[i].columns[j].has_metadata;
                    }
                }
            }
        }

        /* Access key-value metadata */
        for (int32_t i = 0; i < metadata.num_key_value && i < 50; i++) {
            if (metadata.key_value_metadata) {
                (void)metadata.key_value_metadata[i].key;
                (void)metadata.key_value_metadata[i].value;
            }
        }
    }

    parquet_file_metadata_free(&metadata);
    carquet_arena_destroy(&arena);
}

/**
 * Mode 4: Parquet page header parsing
 */
static void fuzz_parquet_page_header(const uint8_t* data, size_t size) {
    parquet_page_header_t header;
    memset(&header, 0, sizeof(header));
    size_t bytes_read = 0;
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_status_t status = parquet_parse_page_header(
        data, size, &header, &bytes_read, &err);

    if (status == CARQUET_OK) {
        (void)header.type;
        (void)header.compressed_page_size;
        (void)header.uncompressed_page_size;
    }
}

/**
 * Mode 5: Thrift encoder-decoder roundtrip
 */
static void fuzz_thrift_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4) return;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    thrift_encoder_t enc;
    thrift_encoder_init(&enc, &buf);

    /* Use fuzz data to generate values to encode */
    int32_t val_i32 = (int32_t)((uint32_t)data[0] | ((uint32_t)data[1] << 8) |
                               ((uint32_t)data[2] << 16) | ((uint32_t)data[3] << 24));

    thrift_write_struct_begin(&enc);
    thrift_write_field_header(&enc, THRIFT_TYPE_I32, 1);
    thrift_write_i32(&enc, val_i32);

    if (size >= 8) {
        int32_t len = (int32_t)(size - 4);
        if (len > 256) len = 256;
        thrift_write_field_header(&enc, THRIFT_TYPE_BINARY, 2);
        thrift_write_binary(&enc, data + 4, len);
    }

    thrift_write_struct_end(&enc);

    /* Decode what we just encoded */
    if (buf.size > 0) {
        thrift_decoder_t dec;
        thrift_decoder_init(&dec, buf.data, buf.size);
        thrift_read_struct_begin(&dec);

        thrift_type_t type;
        int16_t field_id;
        while (!thrift_decoder_has_error(&dec) &&
               thrift_read_field_begin(&dec, &type, &field_id)) {
            if (field_id == 1 && type == THRIFT_TYPE_I32) {
                int32_t decoded = thrift_read_i32(&dec);
                if (decoded != val_i32) __builtin_trap();
            } else {
                thrift_skip_field(&dec, type);
            }
        }
        thrift_read_struct_end(&dec);
    }

    carquet_buffer_destroy(&buf);
}

/**
 * Mode 6: ColumnIndex Thrift deserialization (page-level statistics)
 */
extern carquet_column_index_t* carquet_column_index_parse(const uint8_t* data, size_t size);
extern void carquet_column_index_free(carquet_column_index_t* index);

static void fuzz_column_index(const uint8_t* data, size_t size) {
    carquet_column_index_t* ci = carquet_column_index_parse(data, size);
    if (ci) {
        int32_t np = carquet_column_index_num_pages(ci);
        (void)carquet_column_index_boundary_order(ci);
        for (int32_t i = 0; i < np && i < 20; i++) {
            carquet_page_stats_t ps;
            (void)carquet_column_index_get_page_stats(ci, i, &ps);
            (void)ps.null_count;
            (void)ps.min_value;
            (void)ps.max_value;
            (void)ps.is_null_page;
        }
        carquet_column_index_free(ci);
    }
}

/**
 * Mode 7: OffsetIndex Thrift deserialization (page locations)
 */
extern carquet_offset_index_t* carquet_offset_index_parse(const uint8_t* data, size_t size);
extern void carquet_offset_index_free(carquet_offset_index_t* index);

static void fuzz_offset_index(const uint8_t* data, size_t size) {
    carquet_offset_index_t* oi = carquet_offset_index_parse(data, size);
    if (oi) {
        int32_t np = carquet_offset_index_num_pages(oi);
        for (int32_t i = 0; i < np && i < 20; i++) {
            carquet_page_location_t loc;
            (void)carquet_offset_index_get_page_location(oi, i, &loc);
            (void)loc.offset;
            (void)loc.compressed_size;
            (void)loc.first_row_index;
        }
        carquet_offset_index_free(oi);
    }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 2) return 0;
    (void)carquet_init();

    uint8_t mode = data[0] % 8;
    const uint8_t* payload = data + 1;
    size_t payload_size = size - 1;

    switch (mode) {
        case 0: fuzz_thrift_primitives(payload, payload_size); break;
        case 1: fuzz_thrift_struct(payload, payload_size); break;
        case 2: fuzz_thrift_containers(payload, payload_size); break;
        case 3: fuzz_parquet_metadata(payload, payload_size); break;
        case 4: fuzz_parquet_page_header(payload, payload_size); break;
        case 5: fuzz_thrift_roundtrip(payload, payload_size); break;
        case 6: fuzz_column_index(payload, payload_size); break;
        case 7: fuzz_offset_index(payload, payload_size); break;
    }
    return 0;
}

#ifdef AFL_MAIN
#include <stdio.h>
#include <sys/stat.h>
int main(int argc, char** argv) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <input_file>\n", argv[0]); return 1; }
    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 1; }
    struct stat st; fstat(fileno(f), &st);
    uint8_t* d = malloc((size_t)st.st_size);
    if (!d) { fclose(f); return 1; }
    fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, (size_t)st.st_size);
    free(d); return r;
}
#endif
