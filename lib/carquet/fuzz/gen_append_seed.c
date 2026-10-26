/* Generate a small append-able Parquet file whose schema matches the fixed
 * schema fuzz_append.c installs. Used to seed fuzz_append's corpus so that
 * mutated inputs can still pass append_validate_schema_matches() and exercise
 * the restore-row-groups / close-time footer-rewrite paths. Not part of the
 * normal build. */
#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    const char* out = argc > 1 ? argv[1] : "append_seed.parquet";
    carquet_init();
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT64,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BYTE_ARRAY,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_statistics = true;
    opts.write_page_index = true;  /* exercise the preserved-index region too */

    carquet_writer_t* w = carquet_writer_create(out, schema, &opts, &err);
    if (!w) { fprintf(stderr, "create: %s\n", err.message); return 1; }

    enum { N = 500 };
    int64_t* a = malloc(N * sizeof(int64_t));
    carquet_byte_array_t* b = malloc(N * sizeof(carquet_byte_array_t));
    static char buf[N][16];
    for (int i = 0; i < N; i++) {
        a[i] = (int64_t)i;
        int len = snprintf(buf[i], sizeof(buf[i]), "row%05d", i);
        b[i].data = (uint8_t*)buf[i];
        b[i].length = len;
    }
    if (carquet_writer_write_batch(w, 0, a, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, b, N, NULL, NULL) != CARQUET_OK) {
        fprintf(stderr, "write\n"); return 1;
    }
    carquet_writer_add_metadata(w, "seed_key", "seed_value");
    if (carquet_writer_close(w) != CARQUET_OK) { fprintf(stderr, "close\n"); return 1; }
    carquet_schema_free(schema);
    printf("wrote %s\n", out);
    return 0;
}
