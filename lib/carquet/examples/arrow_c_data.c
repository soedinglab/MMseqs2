/**
 * @file arrow_c_data.c
 * @brief Example: bridge Carquet to the Arrow C Data Interface.
 *
 * Reads a Parquet file into a row batch, exports it as a standard
 * ArrowSchema/ArrowArray (the zero-dependency bridge any Arrow consumer —
 * PyArrow, DuckDB, nanoarrow — understands), then imports it straight back into
 * a second Parquet file via the writer. No Arrow library is required to build
 * or run this; the ABI structs are declared inline in carquet.h.
 */

#include <stdio.h>
#include <stdlib.h>
#include <carquet/carquet.h>

#define NUM_ROWS 100

int main(void) {
    const char* in_path = "/tmp/example_arrow_in.parquet";
    const char* out_path = "/tmp/example_arrow_out.parquet";
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* ── Write a small source file (id REQUIRED, value OPTIONAL) ── */
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL,
                              CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL,
                              CARQUET_REPETITION_OPTIONAL, 0, 0);

    int64_t ids[NUM_ROWS];
    double values[NUM_ROWS];       /* packed non-null values */
    int16_t value_def[NUM_ROWS];   /* 1 = present, 0 = null */
    int packed = 0;
    for (int i = 0; i < NUM_ROWS; i++) {
        ids[i] = i;
        if (i % 5 == 0) {          /* every 5th row is null */
            value_def[i] = 0;
        } else {
            value_def[i] = 1;
            values[packed++] = i * 1.5;
        }
    }

    carquet_writer_t* w = carquet_writer_create(in_path, schema, NULL, &err);
    if (!w) { fprintf(stderr, "write: %s\n", err.message); return 1; }
    carquet_writer_write_batch(w, 0, ids, NUM_ROWS, NULL, NULL);
    carquet_writer_write_batch(w, 1, values, NUM_ROWS, value_def, NULL);
    carquet_writer_close(w);

    /* ── Read a batch and export it to the Arrow C Data Interface ── */
    carquet_reader_t* r = carquet_reader_open(in_path, NULL, &err);
    if (!r) { fprintf(stderr, "read: %s\n", err.message); return 1; }

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = NUM_ROWS;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    carquet_row_batch_t* batch = NULL;
    carquet_batch_reader_next(br, &batch);

    struct ArrowSchema aschema;
    struct ArrowArray aarray;
    const carquet_schema_t* rschema = carquet_reader_schema(r);
    if (carquet_arrow_export_batch(batch, rschema, &aschema, &aarray, &err) != CARQUET_OK) {
        fprintf(stderr, "export: %s\n", err.message);
        return 1;
    }
    printf("Exported %lld rows as an Arrow struct array (%lld columns)\n",
           (long long)aarray.length, (long long)aarray.n_children);
    printf("  column 0 format = %s (%s)\n", aschema.children[0]->format,
           (aschema.children[0]->flags & ARROW_FLAG_NULLABLE) ? "nullable" : "required");
    printf("  column 1 format = %s, null_count = %lld\n",
           aschema.children[1]->format, (long long)aarray.children[1]->null_count);

    /* The export owns independent copies — safe to drop the source now. */
    carquet_batch_reader_free(br);
    carquet_reader_close(r);

    /* ── Import the Arrow array straight into a new Parquet file ── */
    carquet_writer_t* w2 = carquet_writer_create(out_path, schema, NULL, &err);
    if (!w2) { fprintf(stderr, "write2: %s\n", err.message); return 1; }
    /* Consumes (releases) aarray and aschema. */
    if (carquet_writer_write_arrow(w2, &aarray, &aschema, &err) != CARQUET_OK) {
        fprintf(stderr, "write_arrow: %s\n", err.message);
        return 1;
    }
    carquet_writer_close(w2);

    carquet_reader_t* r2 = carquet_reader_open(out_path, NULL, &err);
    printf("Round-tripped through Arrow: %lld rows written to %s\n",
           (long long)carquet_reader_num_rows(r2), out_path);
    carquet_reader_close(r2);

    carquet_schema_free(schema);
    remove(in_path);
    remove(out_path);
    return 0;
}
