/**
 * @file fuzz_append.c
 * @brief Fuzz target for the v0.6.0 append path (carquet_writer_open_append).
 *
 * open_append() parses the trailing footer of an *existing on-disk file*,
 * validates the caller schema against it, restores the prior row groups and
 * key-value metadata into the writer, and (on close) rewrites a fresh footer
 * listing the restored row groups followed by any new ones. None of that is
 * exercised by fuzz_reader (which reads from a buffer) or fuzz_writer (which
 * builds files from scratch).
 *
 * Because open_append takes a path, each input is written to a temp file and
 * opened. The harness installs a fixed schema that matches gen_append_seed.c;
 * mutations of that seed can still pass append_validate_schema_matches() and
 * reach the restore / close-time rewrite paths. On a successful open we also
 * append a row group and close, exercising the footer rewrite over restored
 * (attacker-influenced) row-group metadata.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <carquet/carquet.h>

static carquet_schema_t* make_schema(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) return NULL;
    /* Must match gen_append_seed.c. */
    if (carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT64,
            NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_BYTE_ARRAY,
            NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }
    return schema;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 12) return 0;
    (void)carquet_init();

    /* Write the candidate file to a fixed temp path (libFuzzer is
     * single-threaded by default, so a fixed name is safe). */
    char path[] = "/tmp/carquet_fuzz_append_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) return 0;
    if (write(fd, data, size) != (ssize_t)size) { close(fd); unlink(path); return 0; }
    close(fd);

    carquet_schema_t* schema = make_schema();
    if (!schema) { unlink(path); return 0; }

    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_writer_t* w = carquet_writer_open_append(path, schema, &opts, &err);
    if (w) {
        /* Append a small row group, then close — exercises the footer
         * rewrite over the restored row-group metadata. */
        int64_t a[4] = {1, 2, 3, 4};
        carquet_byte_array_t b[4];
        const char* s[4] = {"x", "yy", "zzz", "wwww"};
        for (int i = 0; i < 4; i++) {
            b[i].data = (uint8_t*)s[i];
            b[i].length = (int32_t)strlen(s[i]);
        }
        if (carquet_writer_write_batch(w, 0, a, 4, NULL, NULL) == CARQUET_OK &&
            carquet_writer_write_batch(w, 1, b, 4, NULL, NULL) == CARQUET_OK) {
            (void)carquet_writer_close(w);
        } else {
            (void)carquet_writer_close(w);
        }
    }

    carquet_schema_free(schema);
    unlink(path);
    return 0;
}

#ifdef AFL_MAIN
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
