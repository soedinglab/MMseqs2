/**
 * @file fuzz_arrow_schema.c
 * @brief Fuzz target for the "ARROW:schema" FlatBuffer field-metadata reader.
 *
 * carquet_apply_arrow_field_metadata() parses an untrusted, base64-encoded
 * Arrow IPC Schema message (the footer "ARROW:schema" blob) with a hand-rolled
 * FlatBuffer reader to recover per-field custom_metadata. This target drives
 * that reader with arbitrary bytes.
 *
 * The fuzz input is base64-encoded here and handed to the parser, so the
 * decoded FlatBuffer buffer is exactly the fuzz bytes — libFuzzer gets direct
 * control over every offset, vtable and vector span the reader walks.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "reader/arrow_schema_read.h"
#include "thrift/parquet_types.h"
#include "core/arena.h"

static char* b64_encode(const uint8_t* in, size_t n) {
    static const char T[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    size_t out_len = ((n + 2) / 3) * 4;
    char* out = (char*)malloc(out_len + 1);
    if (!out) return NULL;
    size_t o = 0;
    for (size_t i = 0; i < n; i += 3) {
        uint32_t v = (uint32_t)in[i] << 16;
        if (i + 1 < n) v |= (uint32_t)in[i + 1] << 8;
        if (i + 2 < n) v |= in[i + 2];
        out[o++] = T[(v >> 18) & 0x3F];
        out[o++] = T[(v >> 12) & 0x3F];
        out[o++] = (i + 1 < n) ? T[(v >> 6) & 0x3F] : '=';
        out[o++] = (i + 2 < n) ? T[v & 0x3F] : '=';
    }
    out[o] = '\0';
    return out;
}

/* A small fixed schema (root + 3 named leaves) so the reader has real targets
 * to match Arrow fields against. */
static void make_schema(parquet_schema_element_t* elems, int32_t* parents) {
    memset(elems, 0, 4 * sizeof(elems[0]));
    static char root[] = "schema";
    static char c0[] = "a", c1[] = "b", c2[] = "c";
    elems[0].name = root; elems[0].num_children = 3;
    elems[1].name = c0; elems[1].has_type = 1; elems[1].type = CARQUET_PHYSICAL_INT32;
    elems[2].name = c1; elems[2].has_type = 1; elems[2].type = CARQUET_PHYSICAL_INT32;
    elems[3].name = c2; elems[3].has_type = 1; elems[3].type = CARQUET_PHYSICAL_BYTE_ARRAY;
    parents[0] = -1; parents[1] = 0; parents[2] = 0; parents[3] = 0;
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size > 1u << 20) return 0;  /* keep inputs bounded */

    char* b64 = b64_encode(data, size);
    if (!b64) return 0;

    parquet_schema_element_t elems[4];
    int32_t parents[4];
    make_schema(elems, parents);

    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) { free(b64); return 0; }

    carquet_apply_arrow_field_metadata(b64, elems, 4, parents, &arena);

    carquet_arena_destroy(&arena);
    free(b64);
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
    uint8_t* d = malloc((size_t)st.st_size ? (size_t)st.st_size : 1);
    if (!d) { fclose(f); return 1; }
    size_t got = fread(d, 1, (size_t)st.st_size, f); fclose(f);
    int r = LLVMFuzzerTestOneInput(d, got);
    free(d); return r;
}
#endif
