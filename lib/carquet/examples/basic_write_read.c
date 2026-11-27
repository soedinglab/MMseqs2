/**
 * @file basic_write_read.c
 * @brief Minimal example: write a Parquet file with Carquet and read it back
 */

#include <stdio.h>
#include <stdlib.h>
#include <carquet/carquet.h>
#include "example_paths.h"

#define NUM_ROWS 1000

int main(void) {
    char filename[1024];
    example_temp_path(filename, sizeof(filename), "example_basic.parquet");
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* ── Schema ── */
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",       CARQUET_PHYSICAL_INT32,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "price",    CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "quantity", CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    /* ── Generate data ── */
    int32_t ids[NUM_ROWS];
    double  prices[NUM_ROWS];
    int64_t quantities[NUM_ROWS];

    for (int i = 0; i < NUM_ROWS; i++) {
        ids[i]        = i + 1;
        prices[i]     = 9.99 + (i % 100) * 0.50;
        quantities[i] = (i % 50) + 1;
    }

    /* ── Write ── */
    carquet_writer_options_t wopts;
    carquet_writer_options_init(&wopts);
    wopts.compression = CARQUET_COMPRESSION_SNAPPY;

    carquet_writer_t* writer = carquet_writer_create(filename, schema, &wopts, &err);
    if (!writer) { fprintf(stderr, "write error: %s\n", err.message); return 1; }

    carquet_writer_write_batch(writer, 0, ids,        NUM_ROWS, NULL, NULL);
    carquet_writer_write_batch(writer, 1, prices,     NUM_ROWS, NULL, NULL);
    carquet_writer_write_batch(writer, 2, quantities, NUM_ROWS, NULL, NULL);
    carquet_writer_close(writer);

    printf("Wrote %d rows to %s\n", NUM_ROWS, filename);

    /* ── Read ── */
    carquet_reader_t* reader = carquet_reader_open(filename, NULL, &err);
    if (!reader) { fprintf(stderr, "read error: %s\n", err.message); return 1; }

    printf("Read back: %lld rows, %d columns, %d row groups\n",
           (long long)carquet_reader_num_rows(reader),
           carquet_reader_num_columns(reader),
           carquet_reader_num_row_groups(reader));

    /* Read first 5 ids from row group 0 */
    carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, 0, &err);
    if (col) {
        int32_t buf[5];
        int64_t n = carquet_column_read_batch(col, buf, 5, NULL, NULL);
        printf("First ids: ");
        for (int i = 0; i < n; i++) printf("%d ", buf[i]);
        printf("\n");
        carquet_column_reader_free(col);
    }

    carquet_reader_close(reader);
    carquet_schema_free(schema);
    remove(filename);
    return 0;
}
