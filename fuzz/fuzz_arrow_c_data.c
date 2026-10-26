/**
 * @file fuzz_arrow_c_data.c
 * @brief Fuzz target for the Arrow C Data Interface import parser.
 *
 * carquet_arrow_import_schema() maps arbitrary Arrow C Data Interface `format`
 * strings (from an external, possibly other-language producer) back to Carquet
 * physical/logical types. The format string is the genuine untrusted-input
 * parse surface of the bridge (the array *buffers* are a trusted ABI contract,
 * not a parse target). This harness builds a struct ArrowSchema whose child
 * `format` strings are carved from the fuzz input and drives the parser.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <carquet/carquet.h>

#define MAX_CHILDREN 64

/* Hand-built ArrowSchemas use stack/fuzzer-owned memory; releasing must free
 * nothing, only mark the struct released. */
static void noop_release(struct ArrowSchema* s) { s->release = NULL; }

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size > (1u << 20)) return 0;  /* keep inputs bounded */

    /* Mutable, NUL-terminated copy so tokens are valid C strings. */
    char* buf = (char*)malloc(size + 1);
    if (!buf) return 0;
    if (size) memcpy(buf, data, size);
    buf[size] = '\0';

    struct ArrowSchema children[MAX_CHILDREN];
    struct ArrowSchema* child_ptrs[MAX_CHILDREN];
    int nchild = 0;

    /* Split the buffer into NUL-separated tokens; each becomes a child format. */
    size_t i = 0;
    while (i <= size && nchild < MAX_CHILDREN) {
        char* tok = buf + i;
        size_t start = i;
        while (i < size && buf[i] != '\0') i++;
        buf[i] = '\0';  /* terminate token */
        i++;            /* skip separator */

        memset(&children[nchild], 0, sizeof(children[nchild]));
        children[nchild].format = tok;
        children[nchild].name = "c";
        /* Use the first byte of the token to toggle nullability for coverage. */
        children[nchild].flags = (tok[0] & 1) ? ARROW_FLAG_NULLABLE : 0;
        children[nchild].release = noop_release;
        child_ptrs[nchild] = &children[nchild];
        nchild++;

        if (start == size) break;  /* consumed everything */
    }

    struct ArrowSchema root;
    memset(&root, 0, sizeof(root));
    /* Alternate a valid struct top-level with a fuzzer-controlled one so both
     * the "must be +s" rejection and the child loop get exercised. */
    root.format = (size > 0 && (data[0] & 1)) ? buf : "+s";
    root.n_children = nchild;
    root.children = child_ptrs;
    root.release = noop_release;

    carquet_schema_t* out = NULL;
    carquet_error_t err = CARQUET_ERROR_INIT;
    if (carquet_arrow_import_schema(&root, &out, &err) == CARQUET_OK && out) {
        carquet_schema_free(out);
    }

    free(buf);
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
