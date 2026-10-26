/* Generate a small Parquet file with a page index, for seeding
 * fuzz_page_filter's corpus. Not part of the normal build. */
#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    const char* out = argc > 1 ? argv[1] : "seed.parquet";
    carquet_init();
    carquet_error_t err = CARQUET_ERROR_INIT;

    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "i32", CARQUET_PHYSICAL_INT32,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "i64", CARQUET_PHYSICAL_INT64,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "dbl", CARQUET_PHYSICAL_DOUBLE,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "str", CARQUET_PHYSICAL_BYTE_ARRAY,
                              NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.write_page_index = true;
    opts.write_statistics = true;
    opts.page_size = 256;  /* force several pages so the page index has entries */

    carquet_writer_t* w = carquet_writer_create(out, schema, &opts, &err);
    if (!w) { fprintf(stderr, "create: %s\n", err.message); return 1; }

    enum { N = 4000 };
    int32_t* i32 = malloc(N * sizeof(int32_t));
    int64_t* i64 = malloc(N * sizeof(int64_t));
    double* dbl = malloc(N * sizeof(double));
    carquet_byte_array_t* str = malloc(N * sizeof(carquet_byte_array_t));
    static char strbuf[N][16];
    for (int i = 0; i < N; i++) {
        i32[i] = i;
        i64[i] = (int64_t)i * 1000;
        dbl[i] = (double)i * 0.5;
        int len = snprintf(strbuf[i], sizeof(strbuf[i]), "v%06d", i);
        str[i].data = (uint8_t*)strbuf[i];
        str[i].length = len;
    }

    if (carquet_writer_write_batch(w, 0, i32, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, i64, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 2, dbl, N, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 3, str, N, NULL, NULL) != CARQUET_OK) {
        fprintf(stderr, "write\n"); return 1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        fprintf(stderr, "close\n"); return 1;
    }
    carquet_schema_free(schema);
    printf("wrote %s\n", out);
    return 0;
}
