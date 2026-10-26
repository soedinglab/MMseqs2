/**
 * @file fuzz_compression.c
 * @brief Fuzz target for carquet compression codecs
 *
 * Tests all 4 compression codecs (Snappy, LZ4, GZIP, ZSTD) with:
 *   - Decompression of arbitrary (likely invalid) data
 *   - Compress-then-decompress roundtrip verification
 *   - Decompression with undersized output buffer
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <carquet/carquet.h>

/* Internal compression functions — correct signatures from implementation */
carquet_status_t carquet_snappy_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
carquet_status_t carquet_snappy_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_snappy_compress_bound(size_t src_size);
carquet_status_t carquet_snappy_get_uncompressed_length(
    const uint8_t* src, size_t src_size, size_t* length);

/* Detailed-diagnostic entry point (src/compression/snappy.h). Declared here
 * rather than included so this target keeps its self-contained style. */
typedef struct carquet_snappy_diag {
    int reason;
    size_t input_pos;
    size_t output_pos;
    uint64_t declared_length;
    int32_t tag;
    uint64_t operand;
} carquet_snappy_diag_t;
carquet_status_t carquet_snappy_decompress_diag(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size,
    carquet_snappy_diag_t* diag);
const char* carquet_snappy_reason_string(int reason);
const char* carquet_snappy_diag_format(const carquet_snappy_diag_t* diag,
                                       char* buf, size_t buf_size);

carquet_status_t carquet_lz4_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
carquet_status_t carquet_lz4_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_lz4_compress_bound(size_t src_size);

int carquet_zstd_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
int carquet_zstd_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);
size_t carquet_zstd_compress_bound(size_t src_size);

int carquet_gzip_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
int carquet_gzip_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);
size_t carquet_gzip_compress_bound(size_t src_size);

/**
 * Decompress arbitrary data — tests robustness against malformed streams.
 */
static void fuzz_decompress(int codec, const uint8_t* data, size_t size) {
    size_t dst_capacity = size * 10 + 256;
    if (dst_capacity > 1024 * 1024) dst_capacity = 1024 * 1024;
    uint8_t* dst = malloc(dst_capacity);
    if (!dst) return;
    size_t dst_size = 0;

    switch (codec) {
        case 0: {
            (void)carquet_snappy_decompress(data, size, dst, dst_capacity, &dst_size);
            { size_t len = 0; (void)carquet_snappy_get_uncompressed_length(data, size, &len); }
            /* Same decode through the diagnostic entry point: exercises every
             * diag-filling branch and the formatter over malformed input. */
            carquet_snappy_diag_t diag;
            size_t diag_dst_size = 0;
            carquet_status_t st = carquet_snappy_decompress_diag(
                data, size, dst, dst_capacity, &diag_dst_size, &diag);
            if (st != CARQUET_OK) {
                char line[256];
                (void)carquet_snappy_diag_format(&diag, line, sizeof(line));
                (void)carquet_snappy_reason_string(diag.reason);
            }
            break;
        }
        case 1: (void)carquet_lz4_decompress(data, size, dst, dst_capacity, &dst_size); break;
        case 2: (void)carquet_gzip_decompress(data, size, dst, dst_capacity, &dst_size); break;
        case 3: (void)carquet_zstd_decompress(data, size, dst, dst_capacity, &dst_size); break;
    }
    free(dst);
}

/**
 * Compress-then-decompress roundtrip verification.
 */
static void fuzz_compress_roundtrip(int codec, const uint8_t* data, size_t size) {
    if (size < 1 || size > 100000) return;

    size_t comp_cap;
    switch (codec) {
        case 0: comp_cap = carquet_snappy_compress_bound(size); break;
        case 1: comp_cap = carquet_lz4_compress_bound(size); break;
        case 2: comp_cap = carquet_gzip_compress_bound(size); break;
        case 3: comp_cap = carquet_zstd_compress_bound(size); break;
        default: return;
    }

    uint8_t* compressed = malloc(comp_cap);
    if (!compressed) return;
    size_t comp_size = 0;
    int status;

    switch (codec) {
        case 0: status = carquet_snappy_compress(data, size, compressed, comp_cap, &comp_size); break;
        case 1: status = carquet_lz4_compress(data, size, compressed, comp_cap, &comp_size); break;
        case 2: status = carquet_gzip_compress(data, size, compressed, comp_cap, &comp_size, 6); break;
        case 3: status = carquet_zstd_compress(data, size, compressed, comp_cap, &comp_size, 1); break;
        default: free(compressed); return;
    }

    if (status == 0 && comp_size > 0) {
        uint8_t* decompressed = malloc(size);
        if (decompressed) {
            size_t dec_size = 0;
            int dec_status;
            switch (codec) {
                case 0: dec_status = carquet_snappy_decompress(compressed, comp_size, decompressed, size, &dec_size); break;
                case 1: dec_status = carquet_lz4_decompress(compressed, comp_size, decompressed, size, &dec_size); break;
                case 2: dec_status = carquet_gzip_decompress(compressed, comp_size, decompressed, size, &dec_size); break;
                case 3: dec_status = carquet_zstd_decompress(compressed, comp_size, decompressed, size, &dec_size); break;
                default: dec_status = -1; break;
            }
            if (dec_status == 0 && (dec_size != size || memcmp(data, decompressed, size) != 0)) {
                __builtin_trap();
            }
            free(decompressed);
        }
    }
    free(compressed);
}

/**
 * Decompress with undersized buffer — tests error handling paths.
 */
static void fuzz_decompress_small_buffer(int codec, const uint8_t* data, size_t size) {
    if (size < 2) return;
    uint8_t dst[16];
    size_t dst_size = 0;
    switch (codec) {
        case 0: (void)carquet_snappy_decompress(data, size, dst, sizeof(dst), &dst_size); break;
        case 1: (void)carquet_lz4_decompress(data, size, dst, sizeof(dst), &dst_size); break;
        case 2: (void)carquet_gzip_decompress(data, size, dst, sizeof(dst), &dst_size); break;
        case 3: (void)carquet_zstd_decompress(data, size, dst, sizeof(dst), &dst_size); break;
    }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 2) return 0;
    (void)carquet_init();

    /* High nibble = mode (0-2), low nibble = codec (0-3) */
    uint8_t mode = (data[0] >> 4) & 0x03;
    int codec = data[0] & 0x03;
    const uint8_t* payload = data + 1;
    size_t payload_size = size - 1;

    switch (mode) {
        case 0: fuzz_decompress(codec, payload, payload_size); break;
        case 1: fuzz_compress_roundtrip(codec, payload, payload_size); break;
        case 2: fuzz_decompress_small_buffer(codec, payload, payload_size); break;
        default: fuzz_decompress(codec, payload, payload_size); break;
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
