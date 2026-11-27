/**
 * @file test_dictionary_fallback.c
 * @brief preserve_dictionaries over a chunk falling back to PLAIN pages
 *
 * A dictionary-encoded column chunk may fall back to PLAIN pages after its
 * dictionary pages (as writers do once the dictionary is full). Read with
 * preserve_dictionaries, the batch reader returns the dictionary part's
 * indices, then fails at the PLAIN page: never the PLAIN values' bytes as
 * indices. carquet's writer falls back for whole chunks only, so the mixed
 * chunk is spliced from two files' chunks: a dictionary-encoded one, then a
 * PLAIN one's data pages, under the first's footer, grown to both.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"
#include "core/arena.h"
#include "core/buffer.h"
#include "thrift/parquet_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DICT_ROWS 3000
#define PLAIN_ROWS 2000
#define DICT_VALUES 7

static int64_t dict_value(int64_t i) { return 100 + i % DICT_VALUES; }
static int64_t plain_value(int64_t i) { return 1000000 + i; }

/* one INT64 REQUIRED column v, uncompressed, without statistics: rows of
 * value(i), dictionary-encoded or PLAIN */
static int write_file(const char* path, int64_t rows, int64_t (*value)(int64_t), bool dict) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema || carquet_schema_add_column(schema, "v", CARQUET_PHYSICAL_INT64, NULL,
                                             CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        return 1;
    }
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    opts.write_statistics = false;
    carquet_writer_t* w = carquet_writer_create(path, schema, &opts, &err);
    if (!w) {
        carquet_schema_free(schema);
        return 1;
    }
    int64_t* vs = malloc(sizeof(int64_t) * (size_t)rows);
    for (int64_t i = 0; i < rows; i++) vs[i] = value(i);
    int rc = (dict && carquet_writer_set_column_encoding(w, 0, CARQUET_ENCODING_RLE_DICTIONARY) != CARQUET_OK) ||
             carquet_writer_write_batch(w, 0, vs, rows, NULL, NULL) != CARQUET_OK;
    rc |= carquet_writer_close(w) != CARQUET_OK;
    free(vs);
    carquet_schema_free(schema);
    return rc;
}

static uint8_t* read_file(const char* path, size_t* size) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t* b = malloc((size_t)n);
    if (b && fread(b, 1, (size_t)n, f) != (size_t)n) {
        free(b);
        b = NULL;
    }
    fclose(f);
    *size = (size_t)n;
    return b;
}

/* a file's footer parsed into arena */
static int parse_footer(const uint8_t* file, size_t size, carquet_arena_t* arena,
                        parquet_file_metadata_t* md) {
    uint32_t len = (uint32_t)file[size - 8] | ((uint32_t)file[size - 7] << 8) |
                   ((uint32_t)file[size - 6] << 16) | ((uint32_t)file[size - 5] << 24);
    carquet_error_t err = CARQUET_ERROR_INIT;
    return parquet_parse_file_metadata(file + size - 8 - len, len, arena, md, &err) != CARQUET_OK;
}

/* the mixed file: the dictionary file's chunk, then the PLAIN file's data
 * pages, under the dictionary file's footer counting both */
static int splice(const char* dict_path, const char* plain_path, const char* out_path) {
    size_t ds = 0, ps = 0;
    uint8_t* d = read_file(dict_path, &ds);
    uint8_t* p = read_file(plain_path, &ps);
    carquet_arena_t arena;
    parquet_file_metadata_t dm, pm;
    int rc = 1;
    if (!d || !p || carquet_arena_init(&arena) != CARQUET_OK) goto out;
    if (parse_footer(d, ds, &arena, &dm) || parse_footer(p, ps, &arena, &pm)) goto out_arena;
    {
        parquet_column_metadata_t* dc = &dm.row_groups[0].columns[0].metadata;
        const parquet_column_metadata_t* pc = &pm.row_groups[0].columns[0].metadata;
        if (!dc->has_dictionary_page_offset || pc->has_dictionary_page_offset) goto out_arena;
        size_t dict_end = (size_t)(dc->dictionary_page_offset + dc->total_compressed_size);
        size_t plain_begin = (size_t)pc->data_page_offset;
        size_t plain_bytes = (size_t)pc->total_compressed_size;

        dc->num_values += pc->num_values;
        dc->total_compressed_size += pc->total_compressed_size;
        dc->total_uncompressed_size += pc->total_uncompressed_size;
        dc->num_encoding_stats = 0;
        dm.row_groups[0].num_rows += pm.row_groups[0].num_rows;
        dm.row_groups[0].total_byte_size += pm.row_groups[0].total_byte_size;
        if (dm.row_groups[0].has_total_compressed_size) {
            dm.row_groups[0].total_compressed_size += pc->total_compressed_size;
        }
        dm.num_rows += pm.num_rows;

        carquet_buffer_t footer;
        carquet_buffer_init(&footer);
        carquet_error_t err = CARQUET_ERROR_INIT;
        FILE* f = fopen(out_path, "wb");
        if (f && parquet_write_file_metadata(&dm, &footer, &err) == CARQUET_OK) {
            size_t n = carquet_buffer_size(&footer);
            uint8_t len[4] = {(uint8_t)n, (uint8_t)(n >> 8), (uint8_t)(n >> 16), (uint8_t)(n >> 24)};
            rc = fwrite(d, 1, dict_end, f) != dict_end ||
                 fwrite(p + plain_begin, 1, plain_bytes, f) != plain_bytes ||
                 fwrite(carquet_buffer_data(&footer), 1, n, f) != n ||
                 fwrite(len, 1, 4, f) != 4 || fwrite("PAR1", 1, 4, f) != 4;
        }
        if (f && fclose(f) != 0) rc = 1;
        carquet_buffer_destroy(&footer);
    }
out_arena:
    carquet_arena_destroy(&arena);
out:
    free(d);
    free(p);
    return rc;
}

int main(void) {
    carquet_init();
    char dict_path[512], plain_path[512], mixed_path[512];
    carquet_test_temp_path(dict_path, sizeof(dict_path), "test_dictionary_fallback_dict");
    carquet_test_temp_path(plain_path, sizeof(plain_path), "test_dictionary_fallback_plain");
    carquet_test_temp_path(mixed_path, sizeof(mixed_path), "test_dictionary_fallback_mixed");
    if (write_file(dict_path, DICT_ROWS, dict_value, true) ||
        write_file(plain_path, PLAIN_ROWS, plain_value, false) ||
        splice(dict_path, plain_path, mixed_path))
        TEST_FAIL("dictionary_fallback", "cannot write the files");

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(mixed_path, NULL, &err);
    if (!r) TEST_FAIL("dictionary_fallback", err.message);

    /* values: the spliced chunk reads whole */
    static int64_t values[DICT_ROWS + PLAIN_ROWS];
    carquet_column_reader_t* c = carquet_reader_get_column(r, 0, 0, &err);
    if (!c || carquet_column_read_batch(c, values, DICT_ROWS + PLAIN_ROWS, NULL, NULL) != DICT_ROWS + PLAIN_ROWS)
        TEST_FAIL("dictionary_fallback", "values of the mixed chunk");
    for (int64_t i = 0; i < DICT_ROWS + PLAIN_ROWS; i++) {
        if (values[i] != (i < DICT_ROWS ? dict_value(i) : plain_value(i - DICT_ROWS)))
            TEST_FAIL("dictionary_fallback", "a value of the mixed chunk");
    }
    carquet_column_reader_free(c);

    /* preserved dictionaries: the dictionary part's indices, then an error */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.preserve_dictionaries = true;
    config.batch_size = 1000;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &config, &err);
    if (!br) TEST_FAIL("dictionary_fallback", err.message);
    int64_t rows = 0;
    carquet_status_t status;
    for (;;) {
        carquet_row_batch_t* batch = NULL;
        status = carquet_batch_reader_next(br, &batch);
        if (status != CARQUET_OK || !batch) break;
        const uint32_t* indices;
        const uint8_t* nulls;
        const uint8_t* dict;
        const uint32_t* offsets;
        int64_t n = 0;
        int32_t count = 0;
        if (carquet_row_batch_column_dictionary(batch, 0, &indices, &nulls, &n, &dict, &count, &offsets) != CARQUET_OK)
            TEST_FAIL("dictionary_fallback", "a batch not dictionary-preserved");
        if (count != DICT_VALUES) TEST_FAIL("dictionary_fallback", "the dictionary");
        for (int64_t i = 0; i < n; i++) {
            if (rows + i >= DICT_ROWS) TEST_FAIL("dictionary_fallback", "PLAIN values read as indices");
            if (indices[i] >= (uint32_t)count)
                TEST_FAIL("dictionary_fallback", "an index past the dictionary");
            int64_t v;
            memcpy(&v, dict + (size_t)indices[i] * sizeof(int64_t), sizeof(v));
            if (v != dict_value(rows + i)) TEST_FAIL("dictionary_fallback", "an index");
        }
        rows += n;
        carquet_row_batch_free(batch);
    }
    if (rows != DICT_ROWS) TEST_FAIL("dictionary_fallback", "not every dictionary-encoded row read");
    if (status == CARQUET_OK || status == CARQUET_ERROR_END_OF_DATA)
        TEST_FAIL("dictionary_fallback", "no error at the PLAIN page");

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    remove(dict_path);
    remove(plain_path);
    remove(mixed_path);
    TEST_PASS("dictionary_fallback");
    printf("\nAll dictionary fallback tests passed.\n");
    return 0;
}
