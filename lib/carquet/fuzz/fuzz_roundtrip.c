/**
 * @file fuzz_roundtrip.c
 * @brief Fuzz target for encode-decode roundtrip consistency
 *
 * Verifies that encode(data) -> decode == data for:
 *   - Delta INT32/INT64
 *   - LZ4, Snappy, GZIP, ZSTD compression
 *   - Byte stream split FLOAT/DOUBLE
 *   - Generic byte stream split INT32/INT64/FIXED_LEN_BYTE_ARRAY payloads
 *   - RLE encode/decode
 *   - Dictionary encode/decode (INT32)
 *   - Delta Length Byte Array / Delta Byte Array
 *   - Plain encode/decode (INT32)
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <carquet/carquet.h>
#include "core/buffer.h"
#include "encoding/rle.h"
#include "encoding/plain.h"

/* Delta encoding */
carquet_status_t carquet_delta_encode_int32(
    const int32_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);
carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data, size_t data_size,
    int32_t* values, int32_t num_values, size_t* bytes_consumed);
carquet_status_t carquet_delta_encode_int64(
    const int64_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);
carquet_status_t carquet_delta_decode_int64(
    const uint8_t* data, size_t data_size,
    int64_t* values, int32_t num_values, size_t* bytes_consumed);

/* Compression */
carquet_status_t carquet_lz4_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
carquet_status_t carquet_lz4_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_lz4_compress_bound(size_t src_size);

carquet_status_t carquet_snappy_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
carquet_status_t carquet_snappy_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_snappy_compress_bound(size_t src_size);

int carquet_gzip_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);
int carquet_gzip_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_gzip_compress_bound(size_t src_size);

int carquet_zstd_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);
int carquet_zstd_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);
size_t carquet_zstd_compress_bound(size_t src_size);

/* Byte stream split */
carquet_status_t carquet_byte_stream_split_encode_float(
    const float* values, int64_t count,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);
carquet_status_t carquet_byte_stream_split_decode_float(
    const uint8_t* data, size_t data_size, float* values, int64_t count);
carquet_status_t carquet_byte_stream_split_encode_double(
    const double* values, int64_t count,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);
carquet_status_t carquet_byte_stream_split_decode_double(
    const uint8_t* data, size_t data_size, double* values, int64_t count);
carquet_status_t carquet_byte_stream_split_encode(
    const uint8_t* values, int64_t count, int32_t type_length,
    uint8_t* output, size_t output_capacity, size_t* bytes_written);
carquet_status_t carquet_byte_stream_split_decode(
    const uint8_t* data, size_t data_size, int32_t type_length,
    uint8_t* values, int64_t count);

/* Dictionary encoding */
carquet_status_t carquet_dictionary_encode_int32(
    const int32_t* values, int64_t count,
    carquet_buffer_t* dict_output, carquet_buffer_t* indices_output);
carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    int32_t* output, int64_t output_count);

/* Byte array encodings */
carquet_status_t carquet_delta_length_encode(
    const carquet_byte_array_t* values, int32_t num_values,
    carquet_buffer_t* output);
carquet_status_t carquet_delta_length_decode(
    const uint8_t* data, size_t data_size,
    carquet_byte_array_t* values, int32_t num_values, size_t* bytes_consumed);
carquet_status_t carquet_delta_strings_encode(
    const carquet_byte_array_t* values, int32_t num_values,
    carquet_buffer_t* output);
carquet_status_t carquet_delta_strings_decode(
    const uint8_t* data, size_t data_size,
    carquet_byte_array_t* values, int32_t num_values,
    uint8_t* work_buffer, size_t work_buffer_size, size_t* bytes_consumed);

/* ── Roundtrip helpers ────────────────────────────────────────────────── */

static void fuzz_delta_int32_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4 || size > 4000) return;
    int32_t count = (int32_t)(size / 4);
    if (count < 1 || count > 1000) return;

    int32_t* input = malloc((size_t)count * sizeof(int32_t));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(int32_t));

    size_t enc_cap = (size_t)count * 10 + 100;
    uint8_t* encoded = malloc(enc_cap);
    if (!encoded) { free(input); return; }

    size_t enc_size = 0;
    if (carquet_delta_encode_int32(input, count, encoded, enc_cap, &enc_size) == CARQUET_OK && enc_size > 0) {
        int32_t* decoded = malloc((size_t)count * sizeof(int32_t));
        if (decoded) {
            size_t consumed = 0;
            if (carquet_delta_decode_int32(encoded, enc_size, decoded, count, &consumed) == CARQUET_OK) {
                for (int32_t i = 0; i < count; i++)
                    if (input[i] != decoded[i]) __builtin_trap();
            }
            free(decoded);
        }
    }
    free(encoded);
    free(input);
}

static void fuzz_delta_int64_roundtrip(const uint8_t* data, size_t size) {
    if (size < 8 || size > 8000) return;
    int32_t count = (int32_t)(size / 8);
    if (count < 1 || count > 1000) return;

    int64_t* input = malloc((size_t)count * sizeof(int64_t));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(int64_t));

    size_t enc_cap = (size_t)count * 20 + 100;
    uint8_t* encoded = malloc(enc_cap);
    if (!encoded) { free(input); return; }

    size_t enc_size = 0;
    if (carquet_delta_encode_int64(input, count, encoded, enc_cap, &enc_size) == CARQUET_OK && enc_size > 0) {
        int64_t* decoded = malloc((size_t)count * sizeof(int64_t));
        if (decoded) {
            size_t consumed = 0;
            if (carquet_delta_decode_int64(encoded, enc_size, decoded, count, &consumed) == CARQUET_OK) {
                for (int32_t i = 0; i < count; i++)
                    if (input[i] != decoded[i]) __builtin_trap();
            }
            free(decoded);
        }
    }
    free(encoded);
    free(input);
}

static void fuzz_lz4_roundtrip(const uint8_t* data, size_t size) {
    if (size < 1 || size > 100000) return;
    size_t comp_cap = carquet_lz4_compress_bound(size);
    uint8_t* compressed = malloc(comp_cap);
    if (!compressed) return;

    size_t comp_size = 0;
    if (carquet_lz4_compress(data, size, compressed, comp_cap, &comp_size) == CARQUET_OK && comp_size > 0) {
        uint8_t* dec = malloc(size);
        if (dec) {
            size_t dec_size = 0;
            if (carquet_lz4_decompress(compressed, comp_size, dec, size, &dec_size) == CARQUET_OK) {
                if (dec_size != size || memcmp(data, dec, size) != 0) __builtin_trap();
            }
            free(dec);
        }
    }
    free(compressed);
}

static void fuzz_bss_float_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4 || size > 40000) return;
    int64_t count = (int64_t)(size / 4);
    if (count < 1 || count > 10000) return;

    float* input = malloc((size_t)count * sizeof(float));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(float));

    size_t out_cap = (size_t)count * sizeof(float);
    uint8_t* encoded = malloc(out_cap);
    if (!encoded) { free(input); return; }

    size_t written = 0;
    if (carquet_byte_stream_split_encode_float(input, count, encoded, out_cap, &written) == CARQUET_OK && written > 0) {
        float* decoded = malloc((size_t)count * sizeof(float));
        if (decoded) {
            if (carquet_byte_stream_split_decode_float(encoded, (size_t)count * sizeof(float), decoded, count) == CARQUET_OK) {
                if (memcmp(input, decoded, (size_t)count * sizeof(float)) != 0) __builtin_trap();
            }
            free(decoded);
        }
    }
    free(encoded);
    free(input);
}

static void fuzz_bss_double_roundtrip(const uint8_t* data, size_t size) {
    if (size < 8 || size > 80000) return;
    int64_t count = (int64_t)(size / 8);
    if (count < 1 || count > 10000) return;

    double* input = malloc((size_t)count * sizeof(double));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(double));

    size_t out_cap = (size_t)count * sizeof(double);
    uint8_t* encoded = malloc(out_cap);
    if (!encoded) { free(input); return; }

    size_t written = 0;
    if (carquet_byte_stream_split_encode_double(input, count, encoded, out_cap, &written) == CARQUET_OK && written > 0) {
        double* decoded = malloc((size_t)count * sizeof(double));
        if (decoded) {
            if (carquet_byte_stream_split_decode_double(encoded, (size_t)count * sizeof(double), decoded, count) == CARQUET_OK) {
                if (memcmp(input, decoded, (size_t)count * sizeof(double)) != 0) __builtin_trap();
            }
            free(decoded);
        }
    }
    free(encoded);
    free(input);
}

static void fuzz_bss_generic_roundtrip(const uint8_t* data, size_t size, int32_t type_len) {
    if (type_len <= 0 || size < (size_t)type_len || size > 80000) return;
    int64_t count = (int64_t)(size / (size_t)type_len);
    if (count < 1 || count > 10000) return;
    size_t bytes = (size_t)count * (size_t)type_len;

    uint8_t* encoded = malloc(bytes);
    uint8_t* decoded = malloc(bytes);
    if (!encoded || !decoded) {
        free(encoded);
        free(decoded);
        return;
    }

    size_t written = 0;
    if (carquet_byte_stream_split_encode(
            data, count, type_len, encoded, bytes, &written) == CARQUET_OK &&
        written == bytes &&
        carquet_byte_stream_split_decode(
            encoded, bytes, type_len, decoded, count) == CARQUET_OK &&
        memcmp(data, decoded, bytes) != 0) {
        __builtin_trap();
    }

    free(encoded);
    free(decoded);
}

static void fuzz_rle_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4 || size > 4000) return;
    int64_t count = (int64_t)(size / 4);
    if (count < 1 || count > 1000) return;

    /* Interpret fuzz data as uint32 values, mask to valid bit width range */
    uint32_t* input = malloc((size_t)count * sizeof(uint32_t));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(uint32_t));

    int bit_width = (data[0] % 16) + 1; /* 1-16 bits */
    uint32_t mask = (1u << bit_width) - 1;
    for (int64_t i = 0; i < count; i++)
        input[i] &= mask;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    if (carquet_rle_encode_all(input, count, bit_width, &buf) == CARQUET_OK && buf.size > 0) {
        uint32_t* decoded = malloc((size_t)count * sizeof(uint32_t));
        if (decoded) {
            int64_t decoded_count = carquet_rle_decode_all(
                buf.data, buf.size, bit_width, decoded, count);
            if (decoded_count == count) {
                for (int64_t i = 0; i < count; i++)
                    if (input[i] != decoded[i]) __builtin_trap();
            }
            free(decoded);
        }
    }
    carquet_buffer_destroy(&buf);
    free(input);
}

static void fuzz_snappy_roundtrip(const uint8_t* data, size_t size) {
    if (size < 1 || size > 100000) return;
    size_t comp_cap = carquet_snappy_compress_bound(size);
    uint8_t* compressed = malloc(comp_cap);
    if (!compressed) return;

    size_t comp_size = 0;
    if (carquet_snappy_compress(data, size, compressed, comp_cap, &comp_size) == CARQUET_OK && comp_size > 0) {
        uint8_t* dec = malloc(size);
        if (dec) {
            size_t dec_size = 0;
            if (carquet_snappy_decompress(compressed, comp_size, dec, size, &dec_size) == CARQUET_OK) {
                if (dec_size != size || memcmp(data, dec, size) != 0) __builtin_trap();
            }
            free(dec);
        }
    }
    free(compressed);
}

static void fuzz_dict_int32_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4 || size > 4000) return;
    int64_t count = (int64_t)(size / 4);
    if (count < 1 || count > 1000) return;

    int32_t* input = malloc((size_t)count * sizeof(int32_t));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(int32_t));

    /* Reduce cardinality so dictionary encoding works well */
    for (int64_t i = 0; i < count; i++)
        input[i] = input[i] % 64;

    carquet_buffer_t dict_buf, idx_buf;
    carquet_buffer_init(&dict_buf);
    carquet_buffer_init(&idx_buf);

    if (carquet_dictionary_encode_int32(input, count, &dict_buf, &idx_buf) == CARQUET_OK) {
        int32_t dict_count = (int32_t)(dict_buf.size / sizeof(int32_t));
        int32_t* decoded = malloc((size_t)count * sizeof(int32_t));
        if (decoded && dict_count > 0) {
            if (carquet_dictionary_decode_int32(
                    dict_buf.data, dict_buf.size, dict_count,
                    idx_buf.data, idx_buf.size,
                    decoded, count) == CARQUET_OK) {
                for (int64_t i = 0; i < count; i++)
                    if (input[i] != decoded[i]) __builtin_trap();
            }
        }
        free(decoded);
    }
    carquet_buffer_destroy(&dict_buf);
    carquet_buffer_destroy(&idx_buf);
    free(input);
}

static void fuzz_plain_int32_roundtrip(const uint8_t* data, size_t size) {
    if (size < 4 || size > 4000) return;
    int64_t count = (int64_t)(size / 4);
    if (count < 1 || count > 1000) return;

    int32_t* input = malloc((size_t)count * sizeof(int32_t));
    if (!input) return;
    memcpy(input, data, (size_t)count * sizeof(int32_t));

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    if (carquet_encode_plain_int32(input, count, &buf) == CARQUET_OK && buf.size > 0) {
        int32_t* decoded = malloc((size_t)count * sizeof(int32_t));
        if (decoded) {
            int64_t consumed = carquet_decode_plain_int32(buf.data, buf.size, decoded, count);
            if (consumed >= 0) {
                for (int64_t i = 0; i < count; i++)
                    if (input[i] != decoded[i]) __builtin_trap();
            }
            free(decoded);
        }
    }
    carquet_buffer_destroy(&buf);
    free(input);
}

static int32_t build_byte_arrays(
    const uint8_t* data,
    size_t size,
    carquet_byte_array_t** arrays_out,
    uint8_t** storage_out,
    bool shared_prefixes) {

    if (size < 2) return 0;
    int32_t count = (int32_t)((data[0] % 32) + 1);
    uint8_t* storage = malloc((size_t)count * 40);
    carquet_byte_array_t* arrays = calloc((size_t)count, sizeof(carquet_byte_array_t));
    if (!storage || !arrays) {
        free(storage);
        free(arrays);
        return 0;
    }

    size_t pos = 1;
    for (int32_t i = 0; i < count; i++) {
        uint8_t* dst = storage + (size_t)i * 40;
        int32_t len = (int32_t)((pos < size ? data[pos++] : (uint8_t)i) % 24);
        if (shared_prefixes) {
            static const char prefix[] = "shared/prefix/";
            size_t prefix_len = sizeof(prefix) - 1;
            memcpy(dst, prefix, prefix_len);
            int32_t suffix_len = len > 8 ? 8 : len;
            for (int32_t j = 0; j < suffix_len; j++) {
                dst[prefix_len + j] = (uint8_t)((pos < size ? data[pos++] : (uint8_t)(i + j)) % 26 + 'a');
            }
            arrays[i].length = (int32_t)prefix_len + suffix_len;
        } else {
            for (int32_t j = 0; j < len; j++) {
                dst[j] = (uint8_t)(pos < size ? data[pos++] : (uint8_t)(i + j));
            }
            arrays[i].length = len;
        }
        arrays[i].data = dst;
    }

    *arrays_out = arrays;
    *storage_out = storage;
    return count;
}

static void fuzz_delta_length_ba_roundtrip(const uint8_t* data, size_t size) {
    carquet_byte_array_t* input = NULL;
    uint8_t* storage = NULL;
    int32_t count = build_byte_arrays(data, size, &input, &storage, false);
    if (count <= 0) return;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    if (carquet_delta_length_encode(input, count, &buf) == CARQUET_OK && buf.size > 0) {
        carquet_byte_array_t* decoded = calloc((size_t)count, sizeof(carquet_byte_array_t));
        if (decoded) {
            size_t consumed = 0;
            if (carquet_delta_length_decode(
                    buf.data, buf.size, decoded, count, &consumed) == CARQUET_OK) {
                for (int32_t i = 0; i < count; i++) {
                    if (decoded[i].length != input[i].length ||
                        memcmp(decoded[i].data, input[i].data, (size_t)input[i].length) != 0) {
                        __builtin_trap();
                    }
                }
            }
            free(decoded);
        }
    }
    carquet_buffer_destroy(&buf);
    free(input);
    free(storage);
}

static void fuzz_delta_byte_array_roundtrip(const uint8_t* data, size_t size) {
    carquet_byte_array_t* input = NULL;
    uint8_t* storage = NULL;
    int32_t count = build_byte_arrays(data, size, &input, &storage, true);
    if (count <= 0) return;

    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    if (carquet_delta_strings_encode(input, count, &buf) == CARQUET_OK && buf.size > 0) {
        carquet_byte_array_t* decoded = calloc((size_t)count, sizeof(carquet_byte_array_t));
        size_t work_size = 0;
        for (int32_t i = 0; i < count; i++) work_size += (size_t)input[i].length;
        uint8_t* work = malloc(work_size ? work_size : 1);
        if (decoded && work) {
            size_t consumed = 0;
            if (carquet_delta_strings_decode(
                    buf.data, buf.size, decoded, count,
                    work, work_size ? work_size : 1, &consumed) == CARQUET_OK) {
                for (int32_t i = 0; i < count; i++) {
                    if (decoded[i].length != input[i].length ||
                        memcmp(decoded[i].data, input[i].data, (size_t)input[i].length) != 0) {
                        __builtin_trap();
                    }
                }
            }
        }
        free(decoded);
        free(work);
    }
    carquet_buffer_destroy(&buf);
    free(input);
    free(storage);
}

static void fuzz_gzip_roundtrip(const uint8_t* data, size_t size) {
    if (size < 1 || size > 100000) return;
    size_t comp_cap = carquet_gzip_compress_bound(size);
    uint8_t* compressed = malloc(comp_cap);
    if (!compressed) return;

    size_t comp_size = 0;
    if (carquet_gzip_compress(data, size, compressed, comp_cap, &comp_size, 6) == 0 && comp_size > 0) {
        uint8_t* dec = malloc(size);
        if (dec) {
            size_t dec_size = 0;
            if (carquet_gzip_decompress(compressed, comp_size, dec, size, &dec_size) == 0) {
                if (dec_size != size || memcmp(data, dec, size) != 0) __builtin_trap();
            }
            free(dec);
        }
    }
    free(compressed);
}

static void fuzz_zstd_roundtrip(const uint8_t* data, size_t size) {
    if (size < 1 || size > 100000) return;
    size_t comp_cap = carquet_zstd_compress_bound(size);
    uint8_t* compressed = malloc(comp_cap);
    if (!compressed) return;

    size_t comp_size = 0;
    if (carquet_zstd_compress(data, size, compressed, comp_cap, &comp_size, 1) == 0 && comp_size > 0) {
        uint8_t* dec = malloc(size);
        if (dec) {
            size_t dec_size = 0;
            if (carquet_zstd_decompress(compressed, comp_size, dec, size, &dec_size) == 0) {
                if (dec_size != size || memcmp(data, dec, size) != 0) __builtin_trap();
            }
            free(dec);
        }
    }
    free(compressed);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size < 2) return 0;
    (void)carquet_init();

    uint8_t mode = data[0] % 16;
    const uint8_t* payload = data + 1;
    size_t payload_size = size - 1;

    switch (mode) {
        case 0:  fuzz_delta_int32_roundtrip(payload, payload_size); break;
        case 1:  fuzz_delta_int64_roundtrip(payload, payload_size); break;
        case 2:  fuzz_lz4_roundtrip(payload, payload_size); break;
        case 3:  fuzz_bss_float_roundtrip(payload, payload_size); break;
        case 4:  fuzz_bss_double_roundtrip(payload, payload_size); break;
        case 5:  fuzz_rle_roundtrip(payload, payload_size); break;
        case 6:  fuzz_snappy_roundtrip(payload, payload_size); break;
        case 7:  fuzz_dict_int32_roundtrip(payload, payload_size); break;
        case 8:  fuzz_plain_int32_roundtrip(payload, payload_size); break;
        case 9:  fuzz_gzip_roundtrip(payload, payload_size); break;
        case 10: fuzz_zstd_roundtrip(payload, payload_size); break;
        case 11: fuzz_bss_generic_roundtrip(payload, payload_size, 4); break;
        case 12: fuzz_bss_generic_roundtrip(payload, payload_size, 8); break;
        case 13: fuzz_bss_generic_roundtrip(payload, payload_size, 6); break;
        case 14: fuzz_delta_length_ba_roundtrip(payload, payload_size); break;
        case 15: fuzz_delta_byte_array_roundtrip(payload, payload_size); break;
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
