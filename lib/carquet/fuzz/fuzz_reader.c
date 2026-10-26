/**
 * @file fuzz_reader.c
 * @brief Fuzz target for carquet Parquet reader
 *
 * Tests the full Parquet file reader with arbitrary input via buffer API.
 * Exercises: schema inspection, batch reader, low-level column reader,
 * statistics access, and type name utilities.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <carquet/carquet.h>

static bool fuzz_row_group_filter(
    const carquet_reader_t* reader,
    int32_t row_group_index,
    void* user_data) {

    (void)reader;
    uint8_t mask = user_data ? *(const uint8_t*)user_data : 0;
    return ((uint8_t)row_group_index & mask) == 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 12) return 0;
    (void)carquet_init();

    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_reader_t* reader = carquet_reader_open_buffer(data, size, NULL, &err);
    if (!reader) return 0;

    int64_t num_rows = carquet_reader_num_rows(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_row_groups = carquet_reader_num_row_groups(reader);
    (void)num_rows;

    /* Schema inspection */
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    if (schema) {
        int32_t num_elements = carquet_schema_num_elements(schema);
        for (int32_t i = 0; i < num_elements && i < 100; i++) {
            const carquet_schema_node_t* node = carquet_schema_get_element(schema, i);
            if (node) {
                const char* name = carquet_schema_node_name(node);
                carquet_physical_type_t ptype = carquet_schema_node_physical_type(node);
                /* Exercise type name utilities */
                (void)carquet_physical_type_name(ptype);
                (void)name;
            }
        }
    }

    /* Batch reader API */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = (int32_t)((data[0] % 64) + 1) * 32;
    config.num_threads = data[1] % 4;
    config.use_mmap = (data[2] & 1) != 0;
    config.preserve_dictionaries = (data[2] & 2) != 0;
    if (data[2] & 4) {
        config.row_group_filter = fuzz_row_group_filter;
        config.row_group_filter_ctx = (void*)&data[3];
    }

    carquet_batch_reader_t* batch_reader = carquet_batch_reader_create(reader, &config, &err);
    if (batch_reader) {
        carquet_row_batch_t* batch = NULL;
        int batch_count = 0;
        while (batch_count < 10 &&
               carquet_batch_reader_next(batch_reader, &batch) == CARQUET_OK && batch) {
            for (int32_t col = 0; col < num_cols && col < 100; col++) {
                const void* data_ptr;
                const uint8_t* nulls;
                int64_t count;
                (void)carquet_row_batch_column(batch, col, &data_ptr, &nulls, &count);
                if (config.preserve_dictionaries) {
                    const uint32_t* indices;
                    const uint8_t* dict_data;
                    int32_t dict_count;
                    const uint32_t* dict_offsets;
                    (void)carquet_row_batch_column_dictionary(
                        batch, col, &indices, &nulls, &count,
                        &dict_data, &dict_count, &dict_offsets);
                }
                /* Repeated (LIST) columns: exercise the reconstruction accessor,
                 * walking every reconstructed offset span. */
                const int32_t* offsets;
                int64_t num_lists;
                const void* values;
                const uint8_t* vvalid;
                const uint8_t* lvalid;
                if (carquet_row_batch_column_list(
                        batch, col, &offsets, &num_lists, &values,
                        &vvalid, &count, &lvalid) == CARQUET_OK) {
                    volatile int64_t acc = 0;
                    for (int64_t li = 0; li < num_lists; li++) {
                        acc += (int64_t)offsets[li + 1] - offsets[li];
                    }
                    (void)acc;
                }
            }
            /* Drive the Arrow C Data export over the whole batch. This walks the
             * flattened leaf/list/map children — including BYTE_ARRAY value
             * buffers — so it dereferences reconstructed value pointers that the
             * column accessors above only count. Release consumes what it built. */
            if (schema) {
                struct ArrowSchema as;
                struct ArrowArray aa;
                memset(&as, 0, sizeof(as));
                memset(&aa, 0, sizeof(aa));
                if (carquet_arrow_export_batch(batch, schema, &as, &aa, NULL)
                        == CARQUET_OK) {
                    if (aa.release) aa.release(&aa);
                    if (as.release) as.release(&as);
                }
            }
            carquet_row_batch_free(batch);
            batch = NULL;
            batch_count++;
        }
        carquet_batch_reader_free(batch_reader);
    }

    /* Row-range read: drive the [offset, limit) driver with attacker-derived
     * bounds over the (possibly malformed) file metadata. */
    carquet_batch_reader_t* rr = carquet_batch_reader_create(reader, &config, &err);
    if (rr) {
        int64_t off = (int64_t)((uint32_t)data[3] | ((uint32_t)data[4] << 8));
        int64_t lim = (data[5] & 1) ? -1
                    : (int64_t)((uint32_t)data[6] | ((uint32_t)data[7] << 8));
        if (carquet_batch_reader_set_row_range(rr, off, lim) == CARQUET_OK) {
            carquet_row_batch_t* batch = NULL;
            int batch_count = 0;
            while (batch_count < 10 &&
                   carquet_batch_reader_next(rr, &batch) == CARQUET_OK && batch) {
                for (int32_t col = 0; col < num_cols && col < 100; col++) {
                    const void* dptr; const uint8_t* nulls; int64_t count;
                    (void)carquet_row_batch_column(batch, col, &dptr, &nulls, &count);
                }
                carquet_row_batch_free(batch);
                batch = NULL;
                batch_count++;
            }
        }
        carquet_batch_reader_free(rr);
    }

    /* Low-level column reader API */
    for (int32_t rg = 0; rg < num_row_groups && rg < 5; rg++) {
        carquet_row_group_metadata_t rg_meta;
        if (carquet_reader_row_group_metadata(reader, rg, &rg_meta) != CARQUET_OK)
            continue;

        int64_t rg_rows = rg_meta.num_rows;
        if (rg_rows <= 0 || rg_rows > 10000) rg_rows = 10000;

        for (int32_t col = 0; col < num_cols && col < 50; col++) {
            carquet_column_reader_t* col_reader =
                carquet_reader_get_column(reader, rg, col, &err);
            if (!col_reader) continue;

            void* values = malloc((size_t)rg_rows * 16);
            int16_t* def_levels = malloc((size_t)rg_rows * sizeof(int16_t));
            int16_t* rep_levels = malloc((size_t)rg_rows * sizeof(int16_t));

            if (values && def_levels && rep_levels) {
                /* Drive the _ex variant with a live error out-param so the
                 * error-propagation branch gets malformed-input coverage, and
                 * assert its return/error contract: a hard failure (-1) must set
                 * the error, and a non-negative count must never report a hard
                 * status code (a partial-read truncation may set one). */
                carquet_error_t read_err = CARQUET_ERROR_INIT;
                int64_t values_read = carquet_column_read_batch_ex(
                    col_reader, values, rg_rows, def_levels, rep_levels, &read_err);
                if (values_read < 0) {
                    assert(read_err.code != CARQUET_OK);
                }
                (void)values_read;
            }

            free(values);
            free(def_levels);
            free(rep_levels);
            carquet_column_reader_free(col_reader);
        }
    }

    /* Exercise encoding/compression name utilities */
    (void)carquet_encoding_name(CARQUET_ENCODING_PLAIN);
    (void)carquet_encoding_name(CARQUET_ENCODING_RLE_DICTIONARY);
    (void)carquet_compression_name(CARQUET_COMPRESSION_SNAPPY);
    (void)carquet_compression_name(CARQUET_COMPRESSION_ZSTD);
    (void)carquet_status_string(CARQUET_OK);
    (void)carquet_status_string(CARQUET_ERROR_DECODE);

    /* Exercise v0.4.0 APIs: bloom filter, page index, kv metadata,
     * column chunk metadata.  These parse Thrift from untrusted file
     * offsets, so they are important fuzzing targets. */
    for (int32_t rg = 0; rg < num_row_groups && rg < 3; rg++) {
        for (int32_t col = 0; col < num_cols && col < 10; col++) {
            /* Bloom filter — parses Thrift BloomFilterHeader */
            carquet_bloom_filter_t* bf =
                carquet_reader_get_bloom_filter(reader, rg, col, NULL);
            if (bf) {
                (void)carquet_bloom_filter_size(bf);
                (void)carquet_bloom_filter_check_i64(bf, 42);
                (void)carquet_bloom_filter_check_i32(bf, 0);
                (void)carquet_bloom_filter_check_double(bf, 3.14);
                (void)carquet_bloom_filter_check_bytes(bf,
                    (const uint8_t*)"test", 4);
                carquet_bloom_filter_destroy(bf);
            }

            /* Column index — parses Thrift ColumnIndex */
            carquet_column_index_t* ci =
                carquet_reader_get_column_index(reader, rg, col, NULL);
            if (ci) {
                int32_t np = carquet_column_index_num_pages(ci);
                (void)carquet_column_index_boundary_order(ci);
                for (int32_t p = 0; p < np && p < 5; p++) {
                    carquet_page_stats_t ps;
                    (void)carquet_column_index_get_page_stats(ci, p, &ps);
                    (void)ps.min_value;
                    (void)ps.null_count;
                }
                carquet_column_index_free(ci);
            }

            /* Offset index — parses Thrift OffsetIndex */
            carquet_offset_index_t* oi =
                carquet_reader_get_offset_index(reader, rg, col, NULL);
            if (oi) {
                int32_t np = carquet_offset_index_num_pages(oi);
                for (int32_t p = 0; p < np && p < 5; p++) {
                    carquet_page_location_t loc;
                    (void)carquet_offset_index_get_page_location(oi, p, &loc);
                    (void)loc.offset;
                }
                carquet_offset_index_free(oi);
            }

            /* Column chunk metadata */
            carquet_column_chunk_metadata_t ccm;
            (void)carquet_reader_column_chunk_metadata(reader, rg, col, &ccm);
        }
    }

    /* Key-value metadata accessors */
    int32_t nkv = carquet_reader_num_metadata(reader);
    for (int32_t i = 0; i < nkv && i < 50; i++) {
        const char *k, *v;
        (void)carquet_reader_get_metadata(reader, i, &k, &v);
    }
    (void)carquet_reader_find_metadata(reader, "pandas");
    (void)carquet_reader_find_metadata(reader, "ARROW:schema");

    carquet_reader_close(reader);
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
