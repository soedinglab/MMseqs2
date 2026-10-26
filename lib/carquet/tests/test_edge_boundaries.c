/**
 * @file test_edge_boundaries.c
 * @brief Boundary value and stress tests for carquet
 *
 * Tests size limits, memory boundaries, large allocations,
 * and numeric boundary conditions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include <carquet/carquet.h>
#include <carquet/error.h>
#include <carquet/types.h>
#include "core/buffer.h"
#include "core/arena.h"

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * External Function Declarations
 * ============================================================================
 */

carquet_status_t carquet_lz4_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

carquet_status_t carquet_lz4_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

size_t carquet_lz4_compress_bound(size_t src_size);

int carquet_zstd_compress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size, int level);

int carquet_zstd_decompress(
    const uint8_t* src, size_t src_size,
    uint8_t* dst, size_t dst_capacity, size_t* dst_size);

carquet_status_t carquet_delta_encode_int32(
    const int32_t* values, int32_t num_values,
    uint8_t* data, size_t data_capacity, size_t* bytes_written);

carquet_status_t carquet_delta_decode_int32(
    const uint8_t* data, size_t data_size,
    int32_t* values, int32_t num_values, size_t* bytes_consumed);

/* ============================================================================
 * Buffer Tests
 * ============================================================================
 */

static int test_buffer_grow_stress(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Grow buffer in small increments to stress reallocation */
    for (int i = 0; i < 10000; i++) {
        uint8_t byte = (uint8_t)(i & 0xFF);
        carquet_status_t status = carquet_buffer_append(&buf, &byte, 1);
        if (status != CARQUET_OK) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("buffer_grow_stress", "Buffer append failed");
        }
    }

    if (carquet_buffer_size(&buf) != 10000) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_grow_stress", "Wrong final size");
    }

    /* Verify contents */
    const uint8_t* data = carquet_buffer_data_const(&buf);
    for (int i = 0; i < 10000; i++) {
        if (data[i] != (uint8_t)(i & 0xFF)) {
            carquet_buffer_destroy(&buf);
            TEST_FAIL("buffer_grow_stress", "Data corruption");
        }
    }

    printf("  [DEBUG] Buffer grew to %zu bytes successfully\n", carquet_buffer_size(&buf));
    carquet_buffer_destroy(&buf);
    TEST_PASS("buffer_grow_stress");
    return 0;
}

static int test_buffer_large_append(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Append large chunk */
    size_t chunk_size = 1024 * 1024;  /* 1MB */
    uint8_t* chunk = malloc(chunk_size);
    if (!chunk) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_large_append", "Failed to allocate test data");
    }
    memset(chunk, 0xAB, chunk_size);

    carquet_status_t status = carquet_buffer_append(&buf, chunk, chunk_size);
    if (status != CARQUET_OK) {
        free(chunk);
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_large_append", "Large append failed");
    }

    printf("  [DEBUG] Buffer accepted %zu bytes\n", carquet_buffer_size(&buf));

    free(chunk);
    carquet_buffer_destroy(&buf);
    TEST_PASS("buffer_large_append");
    return 0;
}

static int test_buffer_reserve_exact(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Reserve exact capacity */
    carquet_status_t status = carquet_buffer_reserve(&buf, 4096);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_reserve_exact", "Reserve failed");
    }

    /* Size should still be 0 */
    if (carquet_buffer_size(&buf) != 0) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_reserve_exact", "Size should be 0 after reserve");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("buffer_reserve_exact");
    return 0;
}

/* ============================================================================
 * Arena Allocator Tests
 * ============================================================================
 */

static int test_arena_many_small_allocs(void) {
    carquet_arena_t arena;
    carquet_arena_init_size(&arena, 4096);

    /* Many small allocations */
    int count = 0;
    for (int i = 0; i < 10000; i++) {
        void* ptr = carquet_arena_alloc(&arena, 16);
        if (ptr) {
            memset(ptr, 0xCC, 16);  /* Touch memory */
            count++;
        }
    }

    printf("  [DEBUG] Arena: %d small allocations succeeded\n", count);

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_many_small_allocs");
    return 0;
}

static int test_arena_alignment(void) {
    carquet_arena_t arena;
    carquet_arena_init_size(&arena, 4096);

    /* Test alignment of allocations */
    for (int i = 0; i < 100; i++) {
        void* ptr = carquet_arena_alloc(&arena, 1 + (i % 32));
        if (ptr) {
            /* Check 8-byte alignment (arena uses 16-byte alignment) */
            if (((uintptr_t)ptr & 7) != 0) {
                carquet_arena_destroy(&arena);
                TEST_FAIL("arena_alignment", "Allocation not 8-byte aligned");
            }
        }
    }

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_alignment");
    return 0;
}

static int test_arena_save_restore(void) {
    carquet_arena_t arena;
    carquet_arena_init_size(&arena, 4096);

    /* Allocate some memory */
    void* ptr1 = carquet_arena_alloc(&arena, 100);
    void* ptr2 = carquet_arena_alloc(&arena, 200);
    (void)ptr1;
    (void)ptr2;

    /* Save state */
    carquet_arena_mark_t mark = carquet_arena_save(&arena);

    /* Allocate more */
    void* ptr3 = carquet_arena_alloc(&arena, 300);
    void* ptr4 = carquet_arena_alloc(&arena, 400);
    (void)ptr3;
    (void)ptr4;

    /* Restore - ptr3 and ptr4 should be "freed" */
    carquet_arena_restore(&arena, mark);

    /* Should be able to allocate from restored position */
    void* ptr5 = carquet_arena_alloc(&arena, 500);
    if (!ptr5) {
        carquet_arena_destroy(&arena);
        TEST_FAIL("arena_save_restore", "Allocation after restore failed");
    }

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_save_restore");
    return 0;
}

/* ============================================================================
 * Large Data Tests
 * ============================================================================
 */

static int test_large_compression_roundtrip(void) {
    /* Test with 10MB of data */
    size_t size = 10 * 1024 * 1024;
    uint8_t* input = malloc(size);
    if (!input) {
        printf("  [SKIP] Could not allocate 10MB for test\n");
        TEST_PASS("large_compression_roundtrip");
        return 0;
    }

    /* Fill with pattern */
    for (size_t i = 0; i < size; i++) {
        input[i] = (uint8_t)((i * 7 + i / 256) & 0xFF);
    }

    size_t bound = carquet_lz4_compress_bound(size);
    uint8_t* compressed = malloc(bound);
    if (!compressed) {
        free(input);
        printf("  [SKIP] Could not allocate compression buffer\n");
        TEST_PASS("large_compression_roundtrip");
        return 0;
    }

    size_t compressed_size;
    carquet_status_t status = carquet_lz4_compress(input, size, compressed, bound, &compressed_size);
    if (status != CARQUET_OK) {
        free(input);
        free(compressed);
        TEST_FAIL("large_compression_roundtrip", "Compression failed");
    }

    printf("  [DEBUG] Compressed 10MB: %zu -> %zu bytes (%.1f%%)\n",
           size, compressed_size, 100.0 * compressed_size / size);

    uint8_t* output = malloc(size);
    if (!output) {
        free(input);
        free(compressed);
        TEST_FAIL("large_compression_roundtrip", "Could not allocate output buffer");
    }

    size_t output_size;
    status = carquet_lz4_decompress(compressed, compressed_size, output, size, &output_size);
    if (status != CARQUET_OK || output_size != size) {
        free(input);
        free(compressed);
        free(output);
        TEST_FAIL("large_compression_roundtrip", "Decompression failed");
    }

    if (memcmp(input, output, size) != 0) {
        free(input);
        free(compressed);
        free(output);
        TEST_FAIL("large_compression_roundtrip", "Data mismatch");
    }

    free(input);
    free(compressed);
    free(output);
    TEST_PASS("large_compression_roundtrip");
    return 0;
}

static int test_large_delta_encoding(void) {
    /* Test delta encoding with many values */
    int32_t count = 100000;
    int32_t* input = malloc(count * sizeof(int32_t));
    if (!input) {
        printf("  [SKIP] Could not allocate input array\n");
        TEST_PASS("large_delta_encoding");
        return 0;
    }

    /* Sequential values - ideal for delta */
    for (int32_t i = 0; i < count; i++) {
        input[i] = i * 3 + 100;
    }

    size_t buffer_size = count * 8;  /* Generous estimate */
    uint8_t* buffer = malloc(buffer_size);
    if (!buffer) {
        free(input);
        TEST_FAIL("large_delta_encoding", "Could not allocate buffer");
    }

    size_t written;
    carquet_status_t status = carquet_delta_encode_int32(input, count, buffer, buffer_size, &written);
    if (status != CARQUET_OK) {
        free(input);
        free(buffer);
        TEST_FAIL("large_delta_encoding", "Encoding failed");
    }

    printf("  [DEBUG] Delta encoded %d int32s: %zu bytes (%.1f bytes/value)\n",
           count, written, (double)written / count);

    int32_t* output = malloc(count * sizeof(int32_t));
    if (!output) {
        free(input);
        free(buffer);
        TEST_FAIL("large_delta_encoding", "Could not allocate output");
    }

    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, count, &consumed);
    if (status != CARQUET_OK) {
        free(input);
        free(buffer);
        free(output);
        TEST_FAIL("large_delta_encoding", "Decoding failed");
    }

    for (int32_t i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            free(input);
            free(buffer);
            free(output);
            TEST_FAIL("large_delta_encoding", "Value mismatch");
        }
    }

    free(input);
    free(buffer);
    free(output);
    TEST_PASS("large_delta_encoding");
    return 0;
}

/* ============================================================================
 * Numeric Boundary Tests
 * ============================================================================
 */

static int test_int32_boundary_delta(void) {
    /* Test INT32_MIN to INT32_MAX transition */
    int32_t input[] = {
        INT32_MIN,
        INT32_MIN + 1,
        -1,
        0,
        1,
        INT32_MAX - 1,
        INT32_MAX
    };
    int count = sizeof(input) / sizeof(input[0]);

    uint8_t buffer[256];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(input, count, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("int32_boundary_delta", "Encoding failed");
    }

    int32_t output[7];
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, count, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("int32_boundary_delta", "Decoding failed");
    }

    for (int i = 0; i < count; i++) {
        if (output[i] != input[i]) {
            printf("  [DEBUG] Mismatch at %d: expected %d, got %d\n", i, input[i], output[i]);
            TEST_FAIL("int32_boundary_delta", "Value mismatch at boundary");
        }
    }

    TEST_PASS("int32_boundary_delta");
    return 0;
}

static int test_max_delta_jump(void) {
    /* Maximum possible delta: from INT32_MIN to INT32_MAX */
    int32_t input[] = {INT32_MIN, INT32_MAX};

    uint8_t buffer[256];
    size_t written;

    carquet_status_t status = carquet_delta_encode_int32(input, 2, buffer, sizeof(buffer), &written);
    if (status != CARQUET_OK) {
        TEST_FAIL("max_delta_jump", "Encoding failed");
    }

    int32_t output[2];
    size_t consumed;
    status = carquet_delta_decode_int32(buffer, written, output, 2, &consumed);
    if (status != CARQUET_OK) {
        TEST_FAIL("max_delta_jump", "Decoding failed");
    }

    if (output[0] != INT32_MIN || output[1] != INT32_MAX) {
        TEST_FAIL("max_delta_jump", "Value mismatch");
    }

    TEST_PASS("max_delta_jump");
    return 0;
}

/* ============================================================================
 * Size Calculation Tests
 * ============================================================================
 */

static int test_compress_bound_accuracy(void) {
    /* Verify compress_bound is always sufficient */
    for (size_t size = 1; size <= 65536; size *= 2) {
        uint8_t* input = malloc(size);
        if (!input) continue;

        /* Fill with random data (worst case) */
        for (size_t i = 0; i < size; i++) {
            input[i] = (uint8_t)(rand() % 256);
        }

        size_t bound = carquet_lz4_compress_bound(size);
        uint8_t* output = malloc(bound);
        if (!output) {
            free(input);
            continue;
        }

        size_t compressed_size;
        carquet_status_t status = carquet_lz4_compress(input, size, output, bound, &compressed_size);

        if (status != CARQUET_OK) {
            printf("  [DEBUG] LZ4 failed for size %zu with bound %zu\n", size, bound);
            free(input);
            free(output);
            TEST_FAIL("compress_bound_accuracy", "Compression failed with calculated bound");
        }

        if (compressed_size > bound) {
            printf("  [DEBUG] Compressed size %zu exceeds bound %zu\n", compressed_size, bound);
            free(input);
            free(output);
            TEST_FAIL("compress_bound_accuracy", "Compressed exceeded bound");
        }

        free(input);
        free(output);
    }

    printf("  [DEBUG] Compress bounds verified for sizes 1 to 65536\n");
    TEST_PASS("compress_bound_accuracy");
    return 0;
}

/* ============================================================================
 * Zero-size Tests
 * ============================================================================
 */

static int test_zero_size_operations(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Append zero bytes */
    carquet_status_t status = carquet_buffer_append(&buf, NULL, 0);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("zero_size_operations", "Zero append failed");
    }

    /* Reserve zero bytes */
    status = carquet_buffer_reserve(&buf, 0);
    if (status != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("zero_size_operations", "Zero reserve failed");
    }

    carquet_buffer_destroy(&buf);
    TEST_PASS("zero_size_operations");
    return 0;
}

/* ============================================================================
 * Caller-supplied values that index or size library memory
 * ============================================================================
 */

extern carquet_status_t carquet_encode_plain_byte_array(
    const carquet_byte_array_t* input, int64_t count, carquet_buffer_t* output);

/* carquet_byte_array_t.length is a signed int32. A negative one must be
 * rejected before anything is reserved: summed as unsigned it contributes 3
 * bytes to the reservation, then copies ~4GB. */
static int test_negative_byte_array_length(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);
    uint8_t payload[4] = {1, 2, 3, 4};
    carquet_byte_array_t in[2] = { { payload, 4 }, { payload, -1 } };
    carquet_status_t status = carquet_encode_plain_byte_array(in, 2, &buf);
    size_t size = buf.size;
    carquet_buffer_destroy(&buf);
    if (status != CARQUET_ERROR_INVALID_ARGUMENT)
        TEST_FAIL("negative_byte_array_length", "negative length was accepted");
    if (size != 0)
        TEST_FAIL("negative_byte_array_length", "output modified before validation");
    TEST_PASS("negative_byte_array_length");
    return 0;
}

/* Projection indices from the batch reader config index the schema's leaf
 * table directly, so an out-of-range one must fail creation. */
static int test_projection_index_bounds(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (!schema) TEST_FAIL("projection_index_bounds", "schema create failed");
    if (carquet_schema_add_column(schema, "a", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK ||
        carquet_schema_add_column(schema, "b", CARQUET_PHYSICAL_INT32, NULL,
                                  CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK)
        TEST_FAIL("projection_index_bounds", "add column failed");
    carquet_writer_t* w = carquet_writer_create_buffer(schema, NULL, &err);
    if (!w) TEST_FAIL("projection_index_bounds", "writer create failed");
    int32_t v[4] = {1, 2, 3, 4};
    if (carquet_writer_write_batch(w, 0, v, 4, NULL, NULL) != CARQUET_OK ||
        carquet_writer_write_batch(w, 1, v, 4, NULL, NULL) != CARQUET_OK ||
        carquet_writer_close(w) != CARQUET_OK)
        TEST_FAIL("projection_index_bounds", "write failed");
    void* data = NULL; size_t size = 0;
    if (carquet_writer_get_buffer(w, &data, &size) != CARQUET_OK)
        TEST_FAIL("projection_index_bounds", "get buffer failed");
    carquet_schema_free(schema);

    carquet_reader_t* r = carquet_reader_open_buffer(data, size, NULL, &err);
    if (!r) { free(data); TEST_FAIL("projection_index_bounds", "open failed"); }
    static const int32_t bad[3][2] = { {0, 2}, {-1, 1}, {INT32_MAX, 0} };
    int rc = 0;
    for (int k = 0; k < 3 && rc == 0; k++) {
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.column_indices = bad[k];
        cfg.num_columns = 2;
        carquet_error_t berr = CARQUET_ERROR_INIT;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &berr);
        if (br || berr.code != CARQUET_ERROR_INVALID_ARGUMENT) rc = 1;
        if (br) carquet_batch_reader_free(br);
    }
    /* A valid projection still works. */
    if (rc == 0) {
        static const int32_t good[2] = {1, 0};
        carquet_batch_reader_config_t cfg;
        carquet_batch_reader_config_init(&cfg);
        cfg.column_indices = good;
        cfg.num_columns = 2;
        carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
        if (!br) rc = 1; else carquet_batch_reader_free(br);
    }
    carquet_reader_close(r);
    free(data);
    if (rc) TEST_FAIL("projection_index_bounds", "out-of-range index not rejected");
    TEST_PASS("projection_index_bounds");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    printf("=== Boundary and Stress Tests ===\n\n");

    /* Initialize carquet for any internal state */
    carquet_init();

    int failures = 0;

    printf("--- Buffer Tests ---\n");
    failures += test_buffer_grow_stress();
    failures += test_buffer_large_append();
    failures += test_buffer_reserve_exact();

    printf("\n--- Arena Allocator Tests ---\n");
    failures += test_arena_many_small_allocs();
    failures += test_arena_alignment();
    failures += test_arena_save_restore();

    printf("\n--- Large Data Tests ---\n");
    failures += test_large_compression_roundtrip();
    failures += test_large_delta_encoding();

    printf("\n--- Numeric Boundary Tests ---\n");
    failures += test_int32_boundary_delta();
    failures += test_max_delta_jump();

    printf("\n--- Size Calculation Tests ---\n");
    failures += test_compress_bound_accuracy();

    printf("\n--- Zero Size Tests ---\n");
    failures += test_zero_size_operations();
    failures += test_negative_byte_array_length();
    failures += test_projection_index_bounds();

    printf("\n");
    if (failures == 0) {
        printf("All boundary and stress tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed.\n", failures);
        return 1;
    }
}
