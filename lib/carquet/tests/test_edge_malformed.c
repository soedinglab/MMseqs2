/**
 * @file test_edge_malformed.c
 * @brief Malformed input tests for carquet
 *
 * Tests that the library handles invalid/malformed input gracefully
 * without crashes, memory corruption, or undefined behavior.
 * These tests simulate fuzzer-found inputs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <carquet/error.h>
#include <carquet/types.h>
#include "core/buffer.h"

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * Function Declarations
 * ============================================================================
 */

/* Compression decompression */
carquet_status_t carquet_lz4_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

carquet_status_t carquet_snappy_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

int carquet_gzip_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

int carquet_zstd_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

/* Delta decoding */
carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data, size_t data_size,
    int32_t* values, int32_t num_values, size_t* bytes_consumed);

carquet_status_t carquet_delta_decode_int64(
    const uint8_t* data, size_t data_size,
    int64_t* values, int32_t num_values, size_t* bytes_consumed);

/* RLE decoding */
int64_t carquet_rle_decode_all(
    const uint8_t* input, size_t input_size, int bit_width,
    uint32_t* output, int64_t max_values);

/* Dictionary decoding */
carquet_status_t carquet_dictionary_decode_int32(
    const uint8_t* dict_data, size_t dict_size, int32_t dict_count,
    const uint8_t* indices_data, size_t indices_size,
    int32_t* output, int64_t output_count);

/* ============================================================================
 * LZ4 Malformed Input Tests
 * ============================================================================
 */

static int test_lz4_truncated(void) {
    /* Truncated LZ4 frame */
    uint8_t truncated[] = {0x04, 0x22, 0x4D, 0x18};  /* Just magic */
    uint8_t output[1024];
    size_t out_size;

    carquet_status_t status = carquet_lz4_decompress(truncated, sizeof(truncated),
                                                      output, sizeof(output), &out_size);
    printf("  [DEBUG] LZ4 truncated: status=%d\n", status);
    /* Should fail gracefully, not crash */

    TEST_PASS("lz4_truncated");
    return 0;
}

static int test_lz4_garbage(void) {
    /* Complete garbage */
    uint8_t garbage[] = {0xDE, 0xAD, 0xBE, 0xEF, 0xCA, 0xFE, 0xBA, 0xBE};
    uint8_t output[1024];
    size_t out_size;

    carquet_status_t status = carquet_lz4_decompress(garbage, sizeof(garbage),
                                                      output, sizeof(output), &out_size);
    printf("  [DEBUG] LZ4 garbage: status=%d\n", status);

    TEST_PASS("lz4_garbage");
    return 0;
}

static int test_lz4_oversized_literal(void) {
    /* LZ4 with claimed oversized literal */
    uint8_t bad[] = {0xFF, 0xFF, 0xFF, 0xFF, 0x00};  /* Claims huge literal */
    uint8_t output[64];
    size_t out_size;

    carquet_status_t status = carquet_lz4_decompress(bad, sizeof(bad),
                                                      output, sizeof(output), &out_size);
    printf("  [DEBUG] LZ4 oversized literal: status=%d\n", status);

    TEST_PASS("lz4_oversized_literal");
    return 0;
}

/* ============================================================================
 * Snappy Malformed Input Tests
 * ============================================================================
 */

static int test_snappy_bad_varint(void) {
    /* Invalid varint (never terminates) */
    uint8_t bad_varint[] = {0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
    uint8_t output[1024];
    size_t out_size;

    carquet_status_t status = carquet_snappy_decompress(bad_varint, sizeof(bad_varint),
                                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] Snappy bad varint: status=%d\n", status);

    TEST_PASS("snappy_bad_varint");
    return 0;
}

static int test_snappy_oversized_length(void) {
    /* Claims to decompress to huge size */
    uint8_t oversized[] = {0xFF, 0xFF, 0xFF, 0xFF, 0x0F};  /* Max uint32 */
    uint8_t output[64];
    size_t out_size;

    carquet_status_t status = carquet_snappy_decompress(oversized, sizeof(oversized),
                                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] Snappy oversized length: status=%d\n", status);

    TEST_PASS("snappy_oversized_length");
    return 0;
}

static int test_snappy_bad_copy_offset(void) {
    /* Snappy with bad copy offset (references before start) */
    uint8_t bad_offset[] = {
        0x05,              /* Uncompressed length = 5 */
        0x00, 'H',         /* Literal 'H' */
        0x01, 0xFF, 0xFF   /* Copy with huge offset */
    };
    uint8_t output[64];
    size_t out_size;

    carquet_status_t status = carquet_snappy_decompress(bad_offset, sizeof(bad_offset),
                                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] Snappy bad copy offset: status=%d\n", status);

    TEST_PASS("snappy_bad_copy_offset");
    return 0;
}

/* ============================================================================
 * GZIP/DEFLATE Malformed Input Tests
 * ============================================================================
 */

static int test_gzip_bad_block_type(void) {
    /* DEFLATE with invalid block type (3 is reserved) */
    uint8_t bad_block[] = {0x07};  /* BFINAL=1, BTYPE=3 (reserved) */
    uint8_t output[64];
    size_t out_size;

    int status = carquet_gzip_decompress(bad_block, sizeof(bad_block),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] GZIP bad block type: status=%d\n", status);

    TEST_PASS("gzip_bad_block_type");
    return 0;
}

static int test_gzip_truncated_dynamic(void) {
    /* Start of dynamic huffman block, truncated */
    uint8_t truncated[] = {0x05, 0x00};  /* BFINAL=1, BTYPE=2, then truncated */
    uint8_t output[64];
    size_t out_size;

    int status = carquet_gzip_decompress(truncated, sizeof(truncated),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] GZIP truncated dynamic: status=%d\n", status);

    TEST_PASS("gzip_truncated_dynamic");
    return 0;
}

static int test_gzip_bad_lengths(void) {
    /* Stored block with mismatched LEN/NLEN */
    uint8_t bad_len[] = {
        0x01,                    /* BFINAL=1, BTYPE=0 (stored) */
        0x05, 0x00,              /* LEN = 5 */
        0x00, 0x00,              /* NLEN should be ~5 = 0xFFFA, but is 0 */
    };
    uint8_t output[64];
    size_t out_size;

    int status = carquet_gzip_decompress(bad_len, sizeof(bad_len),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] GZIP bad lengths: status=%d\n", status);

    TEST_PASS("gzip_bad_lengths");
    return 0;
}

/* ============================================================================
 * ZSTD Malformed Input Tests
 * ============================================================================
 */

static int test_zstd_bad_magic(void) {
    /* Wrong magic number */
    uint8_t bad_magic[] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t output[64];
    size_t out_size;

    int status = carquet_zstd_decompress(bad_magic, sizeof(bad_magic),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] ZSTD bad magic: status=%d\n", status);

    TEST_PASS("zstd_bad_magic");
    return 0;
}

static int test_zstd_truncated_frame(void) {
    /* Valid magic but truncated */
    uint8_t truncated[] = {0x28, 0xB5, 0x2F, 0xFD};  /* Just magic */
    uint8_t output[64];
    size_t out_size;

    int status = carquet_zstd_decompress(truncated, sizeof(truncated),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] ZSTD truncated frame: status=%d\n", status);

    TEST_PASS("zstd_truncated_frame");
    return 0;
}

static int test_zstd_fuzzer_crash_1(void) {
    /* Crash case from Reddit fuzzing feedback */
    static uint8_t crash_input[] = {
        0x28, 0xb5, 0x2f, 0xfd, 0x30, 0x30, 0xfd, 0x00,
        0x00, 0xfd, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30,
        0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30,
        0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30,
        0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30
    };
    uint8_t output[64];
    size_t out_size;

    int status = carquet_zstd_decompress(crash_input, sizeof(crash_input),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] ZSTD fuzzer crash 1: status=%d\n", status);

    TEST_PASS("zstd_fuzzer_crash_1");
    return 0;
}

static int test_zstd_oversized_window(void) {
    /* ZSTD frame claiming huge window size */
    uint8_t oversized[] = {
        0x28, 0xB5, 0x2F, 0xFD,  /* Magic */
        0xFF,                     /* Frame header with max window */
    };
    uint8_t output[64];
    size_t out_size;

    int status = carquet_zstd_decompress(oversized, sizeof(oversized),
                                         output, sizeof(output), &out_size);
    printf("  [DEBUG] ZSTD oversized window: status=%d\n", status);

    TEST_PASS("zstd_oversized_window");
    return 0;
}

/* ============================================================================
 * Delta Encoding Malformed Input Tests
 * ============================================================================
 */

static int test_delta_empty_input(void) {
    int32_t output[10];
    size_t consumed;

    carquet_status_t status = carquet_delta_decode_int32(NULL, 0, output, 10, &consumed);
    printf("  [DEBUG] Delta empty input: status=%d\n", status);

    TEST_PASS("delta_empty_input");
    return 0;
}

static int test_delta_truncated_header(void) {
    /* Delta header requires at least 2 bytes */
    uint8_t truncated[] = {0x80};
    int32_t output[10];
    size_t consumed;

    carquet_status_t status = carquet_delta_decode_int32(truncated, sizeof(truncated),
                                                          output, 10, &consumed);
    printf("  [DEBUG] Delta truncated header: status=%d\n", status);

    TEST_PASS("delta_truncated_header");
    return 0;
}

static int test_delta_bad_block_size(void) {
    /* Block size claiming more values than buffer */
    uint8_t bad[] = {
        0x80, 0x01,  /* Block size = 128 in ULEB128 */
        0x01,        /* Mini-block count = 1 */
        0xFF, 0xFF, 0xFF, 0xFF, 0x0F,  /* Huge total count */
        0x00         /* First value = 0 */
    };
    int32_t output[10];
    size_t consumed;

    carquet_status_t status = carquet_delta_decode_int32(bad, sizeof(bad),
                                                          output, 10, &consumed);
    printf("  [DEBUG] Delta bad block size: status=%d\n", status);

    TEST_PASS("delta_bad_block_size");
    return 0;
}

/* ============================================================================
 * RLE Encoding Malformed Input Tests
 * ============================================================================
 */

static int test_rle_zero_bit_width(void) {
    uint8_t data[] = {0x02, 0x00};  /* Run of 1, value 0 */
    uint32_t output[10];

    /* bit_width=0 is edge case */
    int64_t decoded = carquet_rle_decode_all(data, sizeof(data), 0, output, 10);
    printf("  [DEBUG] RLE zero bit width: decoded=%lld\n", (long long)decoded);

    TEST_PASS("rle_zero_bit_width");
    return 0;
}

static int test_rle_oversized_bit_width(void) {
    uint8_t data[] = {0x02, 0xFF, 0xFF, 0xFF, 0xFF};
    uint32_t output[10];

    /* bit_width > 32 is invalid for int32 */
    int64_t decoded = carquet_rle_decode_all(data, sizeof(data), 64, output, 10);
    printf("  [DEBUG] RLE oversized bit width: decoded=%lld\n", (long long)decoded);

    TEST_PASS("rle_oversized_bit_width");
    return 0;
}

static int test_rle_truncated_run(void) {
    /* Starts a run but data ends */
    uint8_t truncated[] = {0xFE};  /* Large run count, but no value */
    uint32_t output[1000];

    int64_t decoded = carquet_rle_decode_all(truncated, sizeof(truncated), 8, output, 1000);
    printf("  [DEBUG] RLE truncated run: decoded=%lld\n", (long long)decoded);

    TEST_PASS("rle_truncated_run");
    return 0;
}

/* ============================================================================
 * Dictionary Encoding Malformed Input Tests
 * ============================================================================
 */

static int test_dictionary_empty_dict(void) {
    /* Empty dictionary with indices pointing to it */
    uint8_t indices[] = {0x02, 0x00};  /* RLE: run of 1, index 0 */
    int32_t output[10];

    carquet_status_t status = carquet_dictionary_decode_int32(
        NULL, 0, 0,              /* Empty dictionary */
        indices, sizeof(indices),
        output, 10);
    printf("  [DEBUG] Dictionary empty dict: status=%d\n", status);

    TEST_PASS("dictionary_empty_dict");
    return 0;
}

static int test_dictionary_index_out_of_bounds(void) {
    /* Dictionary with 2 values, indices referencing index 100 */
    int32_t dict[] = {42, 84};
    uint8_t indices[] = {0x02, 0x64};  /* RLE: run of 1, index 100 */
    int32_t output[10];

    carquet_status_t status = carquet_dictionary_decode_int32(
        (uint8_t*)dict, sizeof(dict), 2,
        indices, sizeof(indices),
        output, 10);
    printf("  [DEBUG] Dictionary OOB index: status=%d\n", status);

    TEST_PASS("dictionary_index_out_of_bounds");
    return 0;
}

/* ============================================================================
 * Random Garbage Tests
 * ============================================================================
 */

static int test_random_garbage_all_codecs(void) {
    /* Feed random garbage to all decoders */
    uint8_t garbage[256];
    srand(12345);
    for (int i = 0; i < 256; i++) {
        garbage[i] = (uint8_t)(rand() % 256);
    }

    uint8_t output[1024];
    size_t out_size;
    int32_t int_output[100];
    uint32_t uint_output[100];
    size_t consumed;

    printf("  Testing random garbage against all decoders...\n");

    /* LZ4 */
    carquet_lz4_decompress(garbage, sizeof(garbage), output, sizeof(output), &out_size);

    /* Snappy */
    carquet_snappy_decompress(garbage, sizeof(garbage), output, sizeof(output), &out_size);

    /* GZIP */
    carquet_gzip_decompress(garbage, sizeof(garbage), output, sizeof(output), &out_size);

    /* ZSTD */
    carquet_zstd_decompress(garbage, sizeof(garbage), output, sizeof(output), &out_size);

    /* Delta */
    carquet_delta_decode_int32(garbage, sizeof(garbage), int_output, 100, &consumed);

    /* RLE */
    carquet_rle_decode_all(garbage, sizeof(garbage), 8, uint_output, 100);

    printf("  [DEBUG] All decoders survived random garbage\n");
    TEST_PASS("random_garbage_all_codecs");
    return 0;
}

static int test_all_zeros_input(void) {
    /* All zeros - edge case for many decoders */
    uint8_t zeros[256] = {0};
    uint8_t output[1024];
    size_t out_size;

    printf("  Testing all-zeros input...\n");

    carquet_lz4_decompress(zeros, sizeof(zeros), output, sizeof(output), &out_size);
    carquet_snappy_decompress(zeros, sizeof(zeros), output, sizeof(output), &out_size);
    carquet_gzip_decompress(zeros, sizeof(zeros), output, sizeof(output), &out_size);
    carquet_zstd_decompress(zeros, sizeof(zeros), output, sizeof(output), &out_size);

    printf("  [DEBUG] All decoders survived all-zeros input\n");
    TEST_PASS("all_zeros_input");
    return 0;
}

static int test_all_0xff_input(void) {
    /* All 0xFF - another edge case */
    uint8_t ones[256];
    memset(ones, 0xFF, sizeof(ones));
    uint8_t output[1024];
    size_t out_size;

    printf("  Testing all-0xFF input...\n");

    carquet_lz4_decompress(ones, sizeof(ones), output, sizeof(output), &out_size);
    carquet_snappy_decompress(ones, sizeof(ones), output, sizeof(output), &out_size);
    carquet_gzip_decompress(ones, sizeof(ones), output, sizeof(output), &out_size);
    carquet_zstd_decompress(ones, sizeof(ones), output, sizeof(output), &out_size);

    printf("  [DEBUG] All decoders survived all-0xFF input\n");
    TEST_PASS("all_0xff_input");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    printf("=== Malformed Input Tests ===\n\n");

    int failures = 0;

    printf("--- LZ4 Malformed Input ---\n");
    failures += test_lz4_truncated();
    failures += test_lz4_garbage();
    failures += test_lz4_oversized_literal();

    printf("\n--- Snappy Malformed Input ---\n");
    failures += test_snappy_bad_varint();
    failures += test_snappy_oversized_length();
    failures += test_snappy_bad_copy_offset();

    printf("\n--- GZIP Malformed Input ---\n");
    failures += test_gzip_bad_block_type();
    failures += test_gzip_truncated_dynamic();
    failures += test_gzip_bad_lengths();

    printf("\n--- ZSTD Malformed Input ---\n");
    failures += test_zstd_bad_magic();
    failures += test_zstd_truncated_frame();
    failures += test_zstd_fuzzer_crash_1();
    failures += test_zstd_oversized_window();

    printf("\n--- Delta Encoding Malformed Input ---\n");
    failures += test_delta_empty_input();
    failures += test_delta_truncated_header();
    failures += test_delta_bad_block_size();

    printf("\n--- RLE Encoding Malformed Input ---\n");
    failures += test_rle_zero_bit_width();
    failures += test_rle_oversized_bit_width();
    failures += test_rle_truncated_run();

    printf("\n--- Dictionary Encoding Malformed Input ---\n");
    failures += test_dictionary_empty_dict();
    failures += test_dictionary_index_out_of_bounds();

    printf("\n--- Random Garbage Tests ---\n");
    failures += test_random_garbage_all_codecs();
    failures += test_all_zeros_input();
    failures += test_all_0xff_input();

    printf("\n");
    if (failures == 0) {
        printf("All malformed input tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed.\n", failures);
        return 1;
    }
}
