/**
 * @file test_core.c
 * @brief Tests for core utilities
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include "core/arena.h"
#include "core/buffer.h"
#include "core/endian.h"
#include "core/bitpack.h"

extern int64_t carquet_dispatch_find_run_length_i32(const int32_t* values, int64_t count);
extern void carquet_dispatch_minmax_i32(const int32_t* values, int64_t count,
                                         int32_t* min_value, int32_t* max_value);
extern void carquet_dispatch_minmax_i64(const int64_t* values, int64_t count,
                                         int64_t* min_value, int64_t* max_value);
extern void carquet_dispatch_minmax_float(const float* values, int64_t count,
                                           float* min_value, float* max_value);
extern void carquet_dispatch_minmax_double(const double* values, int64_t count,
                                            double* min_value, double* max_value);
extern void carquet_dispatch_copy_minmax_i32(const int32_t* values, int64_t count, int32_t* output,
                                             int32_t* min_value, int32_t* max_value);
extern void carquet_dispatch_copy_minmax_float(const float* values, int64_t count, float* output,
                                               float* min_value, float* max_value);
extern void carquet_dispatch_copy_minmax_double(const double* values, int64_t count, double* output,
                                                double* min_value, double* max_value);
extern int64_t carquet_dispatch_count_non_nulls(const int16_t* def_levels, int64_t count,
                                                int16_t max_def_level);
extern void carquet_dispatch_build_null_bitmap(const int16_t* def_levels, int64_t count,
                                               int16_t max_def_level, uint8_t* null_bitmap);
extern void carquet_dispatch_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * Arena Tests
 * ============================================================================
 */

static int test_arena_basic(void) {
    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) {
        TEST_FAIL("arena_basic", "failed to initialize arena");
    }

    /* Test basic allocation */
    void* p1 = carquet_arena_alloc(&arena, 100);
    if (p1 == NULL) {
        carquet_arena_destroy(&arena);
        TEST_FAIL("arena_basic", "first allocation returned NULL");
    }

    void* p2 = carquet_arena_alloc(&arena, 200);
    (void)p2;
    assert(p2 != NULL);
    assert(p2 != p1);

    /* Test calloc */
    int* arr = carquet_arena_calloc(&arena, 10, sizeof(int));
    (void)arr;
    assert(arr != NULL);
    for (int i = 0; i < 10; i++) {
        assert(arr[i] == 0);
    }

    /* Test strdup */
    char* s = carquet_arena_strdup(&arena, "Hello, World!");
    (void)s;
    assert(s != NULL);
    assert(strcmp(s, "Hello, World!") == 0);

    /* Test reset */
    size_t allocated_before = carquet_arena_allocated(&arena);
    (void)allocated_before;
    assert(allocated_before > 0);

    carquet_arena_reset(&arena);
    assert(carquet_arena_allocated(&arena) == 0);

    /* Can still allocate after reset */
    void* p3 = carquet_arena_alloc(&arena, 50);
    (void)p3;
    assert(p3 != NULL);

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_basic");
    return 0;
}

static int test_arena_large_allocation(void) {
    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) {
        TEST_FAIL("arena_large_allocation", "failed to initialize arena");
    }

    /* Allocate something larger than default block size */
    void* p = carquet_arena_alloc(&arena, 256 * 1024);
    if (p == NULL) {
        carquet_arena_destroy(&arena);
        TEST_FAIL("arena_large_allocation", "large allocation returned NULL");
    }

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_large_allocation");
    return 0;
}

static int test_arena_save_restore(void) {
    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) {
        TEST_FAIL("arena_save_restore", "failed to initialize arena");
    }

    void* p1 = carquet_arena_alloc(&arena, 100);
    if (p1 == NULL) {
        carquet_arena_destroy(&arena);
        TEST_FAIL("arena_save_restore", "first allocation returned NULL");
    }

    carquet_arena_mark_t mark = carquet_arena_save(&arena);
    (void)mark;
    size_t allocated_at_mark = carquet_arena_allocated(&arena);
    (void)allocated_at_mark;

    void* p2 = carquet_arena_alloc(&arena, 200);
    (void)p2;
    assert(p2 != NULL);

    carquet_arena_restore(&arena, mark);
    assert(carquet_arena_allocated(&arena) == allocated_at_mark);

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_save_restore");
    return 0;
}

static int test_arena_overflow_rejected(void) {
    carquet_arena_t arena;
    if (carquet_arena_init(&arena) != CARQUET_OK) {
        TEST_FAIL("arena_overflow_rejected", "init failed");
    }

    void* p = carquet_arena_alloc(&arena, SIZE_MAX);
    if (p != NULL) {
        carquet_arena_destroy(&arena);
        TEST_FAIL("arena_overflow_rejected", "oversized allocation succeeded");
    }

    carquet_arena_destroy(&arena);
    TEST_PASS("arena_overflow_rejected");
    return 0;
}

/* ============================================================================
 * Buffer Tests
 * ============================================================================
 */

static int test_buffer_basic(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    /* Test append */
    uint8_t data[] = {1, 2, 3, 4, 5};
    if (carquet_buffer_append(&buf, data, 5) != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_basic", "failed to append data");
    }
    assert(carquet_buffer_size(&buf) == 5);

    /* Test append byte */
    assert(carquet_buffer_append_byte(&buf, 6) == CARQUET_OK);
    assert(carquet_buffer_size(&buf) == 6);

    /* Verify contents */
    const uint8_t* ptr = carquet_buffer_data_const(&buf);
    (void)ptr;
    assert(ptr[0] == 1 && ptr[5] == 6);

    carquet_buffer_destroy(&buf);
    TEST_PASS("buffer_basic");
    return 0;
}

static int test_buffer_integers(void) {
    carquet_buffer_t buf;
    carquet_buffer_init(&buf);

    if (carquet_buffer_append_u16_le(&buf, 0x1234) != CARQUET_OK) {
        carquet_buffer_destroy(&buf);
        TEST_FAIL("buffer_integers", "failed to append u16");
    }
    assert(carquet_buffer_append_u32_le(&buf, 0x12345678) == CARQUET_OK);
    assert(carquet_buffer_append_u64_le(&buf, 0x123456789ABCDEF0ULL) == CARQUET_OK);

    const uint8_t* ptr = carquet_buffer_data_const(&buf);
    (void)ptr;

    /* Verify u16 */
    assert(ptr[0] == 0x34 && ptr[1] == 0x12);

    /* Verify u32 */
    assert(carquet_read_u32_le(ptr + 2) == 0x12345678);

    /* Verify u64 */
    assert(carquet_read_u64_le(ptr + 6) == 0x123456789ABCDEF0ULL);

    carquet_buffer_destroy(&buf);
    TEST_PASS("buffer_integers");
    return 0;
}

static int test_buffer_reader(void) {
    uint8_t data[] = {0x34, 0x12, 0x78, 0x56, 0x34, 0x12};
    carquet_buffer_reader_t reader;
    carquet_buffer_reader_init_data(&reader, data, sizeof(data));

    uint16_t u16;
    (void)u16;
    assert(carquet_buffer_reader_read_u16_le(&reader, &u16) == CARQUET_OK);
    assert(u16 == 0x1234);

    uint32_t u32;
    (void)u32;
    assert(carquet_buffer_reader_read_u32_le(&reader, &u32) == CARQUET_OK);
    assert(u32 == 0x12345678);

    assert(carquet_buffer_reader_remaining(&reader) == 0);

    TEST_PASS("buffer_reader");
    return 0;
}

static int test_buffer_append_overflow_rejected(void) {
    carquet_buffer_t buf;
    uint8_t storage[8] = {0};
    uint8_t byte = 1;

    carquet_buffer_init_wrap(&buf, storage, sizeof(storage));
    if (carquet_buffer_append(&buf, &byte, SIZE_MAX) != CARQUET_ERROR_OUT_OF_MEMORY) {
        TEST_FAIL("buffer_append_overflow_rejected", "overflow append was not rejected");
    }
    if (carquet_buffer_size(&buf) != sizeof(storage)) {
        TEST_FAIL("buffer_append_overflow_rejected", "buffer size changed after failed append");
    }

    TEST_PASS("buffer_append_overflow_rejected");
    return 0;
}

/* ============================================================================
 * Endian Tests
 * ============================================================================
 */

static int test_endian_read_write(void) {
    uint8_t buf[8];

    /* Test 16-bit */
    carquet_write_u16_le(buf, 0x1234);
    assert(carquet_read_u16_le(buf) == 0x1234);

    /* Test 32-bit */
    carquet_write_u32_le(buf, 0x12345678);
    assert(carquet_read_u32_le(buf) == 0x12345678);

    /* Test 64-bit */
    carquet_write_u64_le(buf, 0x123456789ABCDEF0ULL);
    assert(carquet_read_u64_le(buf) == 0x123456789ABCDEF0ULL);

    TEST_PASS("endian_read_write");
    return 0;
}

static int test_varint(void) {
    uint8_t buf[10];
    uint32_t val32;
    uint64_t val64;

    /* Test small value */
    int len = carquet_encode_varint32(buf, 127);
    (void)len;
    assert(len == 1);
    int consumed = carquet_decode_varint32(buf, 10, &val32);
    (void)consumed;
    assert(consumed == 1 && val32 == 127);

    /* Test larger value */
    len = carquet_encode_varint32(buf, 300);
    assert(len == 2);
    consumed = carquet_decode_varint32(buf, 10, &val32);
    assert(consumed == 2 && val32 == 300);

    /* Test 64-bit */
    len = carquet_encode_varint64(buf, 0x123456789ULL);
    consumed = carquet_decode_varint64(buf, 10, &val64);
    assert(val64 == 0x123456789ULL);

    TEST_PASS("varint");
    return 0;
}

static int test_zigzag(void) {
    /* Positive numbers */
    assert(carquet_zigzag_encode32(0) == 0);
    assert(carquet_zigzag_encode32(1) == 2);
    assert(carquet_zigzag_encode32(2) == 4);

    /* Negative numbers */
    assert(carquet_zigzag_encode32(-1) == 1);
    assert(carquet_zigzag_encode32(-2) == 3);

    /* Roundtrip */
    for (int32_t i = -1000; i <= 1000; i++) {
        uint32_t encoded = carquet_zigzag_encode32(i);
        (void)encoded;
        int32_t decoded = carquet_zigzag_decode32(encoded);
        (void)decoded;
        assert(decoded == i);
    }

    TEST_PASS("zigzag");
    return 0;
}

/* ============================================================================
 * Bitpack Tests
 * ============================================================================
 */

static int test_bitpack_1bit(void) {
    uint8_t input = 0xAA;  /* 0b10101010 */
    uint32_t output[8];

    carquet_bitunpack8_1bit(&input, output);

    assert(output[0] == 0);
    assert(output[1] == 1);
    assert(output[2] == 0);
    assert(output[3] == 1);
    assert(output[4] == 0);
    assert(output[5] == 1);
    assert(output[6] == 0);
    assert(output[7] == 1);

    TEST_PASS("bitpack_1bit");
    return 0;
}

static int test_bitpack_4bit(void) {
    uint8_t input[4] = {0x21, 0x43, 0x65, 0x87};
    uint32_t output[8];

    carquet_bitunpack8_4bit(input, output);

    assert(output[0] == 1);
    assert(output[1] == 2);
    assert(output[2] == 3);
    assert(output[3] == 4);
    assert(output[4] == 5);
    assert(output[5] == 6);
    assert(output[6] == 7);
    assert(output[7] == 8);

    TEST_PASS("bitpack_4bit");
    return 0;
}

static int test_bitpack_roundtrip(void) {
    uint32_t original[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    uint8_t packed[32];
    uint32_t unpacked[8];

    for (int bit_width = 1; bit_width <= 8; bit_width++) {
        /* Clear buffers */
        memset(packed, 0, sizeof(packed));
        memset(unpacked, 0, sizeof(unpacked));

        /* Pack */
        carquet_bitpack8_32(original, bit_width, packed);

        /* Unpack */
        carquet_bitunpack8_32(packed, bit_width, unpacked);

        /* Verify */
        uint32_t mask = (1U << bit_width) - 1;
        (void)mask;
        for (int i = 0; i < 8; i++) {
            assert((original[i] & mask) == unpacked[i]);
        }
    }

    TEST_PASS("bitpack_roundtrip");
    return 0;
}

static int test_bitpack_16bit(void) {
    uint32_t original[8] = {
        0x0001, 0x0102, 0x1234, 0x4567,
        0x89AB, 0xBEEF, 0xCAFE, 0xFFFF
    };
    uint8_t packed[16];
    uint32_t unpacked[8];

    memset(packed, 0, sizeof(packed));
    memset(unpacked, 0, sizeof(unpacked));

    carquet_bitpack8_32(original, 16, packed);
    carquet_bitunpack8_32(packed, 16, unpacked);

    for (int i = 0; i < 8; i++) {
        assert(unpacked[i] == (original[i] & 0xFFFFu));
    }

    TEST_PASS("bitpack_16bit");
    return 0;
}

static int test_run_detection(void) {
    int32_t all_same[32];
    int32_t mismatch_early[16];
    int32_t mismatch_late[33];

    for (int i = 0; i < 32; i++) {
        all_same[i] = 7;
    }
    for (int i = 0; i < 16; i++) {
        mismatch_early[i] = 11;
    }
    for (int i = 0; i < 33; i++) {
        mismatch_late[i] = -3;
    }

    mismatch_early[5] = 12;
    mismatch_late[31] = 99;

    assert(carquet_dispatch_find_run_length_i32(NULL, 0) == 0);
    assert(carquet_dispatch_find_run_length_i32(all_same, 32) == 32);
    assert(carquet_dispatch_find_run_length_i32(mismatch_early, 16) == 5);
    assert(carquet_dispatch_find_run_length_i32(mismatch_late, 33) == 31);

    TEST_PASS("run_detection");
    return 0;
}

static int test_dispatch_minmax(void) {
    int32_t i32v[] = {7, -2, 9, 4, -5, 6, 0, 3, 12};
    int64_t i64v[] = {42, -100, 999, 8, -7, 77, 3};
    float fv[] = {3.5f, -1.0f, 9.25f, 4.0f, -7.5f, 6.25f, 0.0f};
    double dv[] = {8.0, -3.0, 15.5, 1.0, -9.0, 2.5, 14.0};
    int32_t i32min, i32max;
    int64_t i64min, i64max;
    float fmin, fmax;
    double dmin, dmax;

    carquet_dispatch_minmax_i32(i32v, 9, &i32min, &i32max);
    carquet_dispatch_minmax_i64(i64v, 7, &i64min, &i64max);
    carquet_dispatch_minmax_float(fv, 7, &fmin, &fmax);
    carquet_dispatch_minmax_double(dv, 7, &dmin, &dmax);

    assert(i32min == -5 && i32max == 12);
    assert(i64min == -100 && i64max == 999);
    assert(fmin == -7.5f && fmax == 9.25f);
    assert(dmin == -9.0 && dmax == 15.5);

    TEST_PASS("dispatch_minmax");
    return 0;
}

static int test_dispatch_copy_minmax(void) {
    int32_t src[] = {4, -2, 9, 1, -7, 3, 11, 0};
    int32_t dst[8] = {0};
    int32_t min_v = 0;
    int32_t max_v = 0;

    carquet_dispatch_copy_minmax_i32(src, 8, dst, &min_v, &max_v);

    assert(memcmp(src, dst, sizeof(src)) == 0);
    assert(min_v == -7);
    assert(max_v == 11);

    TEST_PASS("dispatch_copy_minmax");
    return 0;
}

/* The page writer relies on this contract for float/double statistics: the
 * dispatched min/max either returns NaN (it then rescans NaN-skipping) or the
 * exact min/max of the non-NaN values. A backend that lets a NaN lane swallow
 * the real extreme and then drops the NaN in its horizontal reduction writes
 * wrong statistics (#35: {1,9,2,NaN,4,3,2,1} gave max 4 on NEON). */
#define MINMAX_NAN_OK(mn, mx, rmn, rmx, found) \
    (isnan(mn) || isnan(mx) || ((found) && (mn) == (rmn) && (mx) == (rmx)))

static int test_dispatch_minmax_nan(void) {
    const double pr_example[8] = {1, 9, 2, NAN, 4, 3, 2, 1};
    double d[64], dcopy[64];
    float f[64], fcopy[64];
    double dmin, dmax;
    float fmin, fmax;
    uint32_t seed = 12345;

    carquet_dispatch_minmax_double(pr_example, 8, &dmin, &dmax);
    if (!MINMAX_NAN_OK(dmin, dmax, 1.0, 9.0, 1)) {
        TEST_FAIL("dispatch_minmax_nan", "minmax_double dropped a value next to a NaN");
    }
    carquet_dispatch_copy_minmax_double(pr_example, 8, dcopy, &dmin, &dmax);
    if (!MINMAX_NAN_OK(dmin, dmax, 1.0, 9.0, 1)) {
        TEST_FAIL("dispatch_minmax_nan", "copy_minmax_double dropped a value next to a NaN");
    }

    /* Every length up to 64 (covers each backend's vector body and tail),
     * NaN in roughly one value in six, including values[0]. */
    for (int iter = 0; iter < 20000; iter++) {
        int64_t n = 1 + (int64_t)(iter % 64);
        double rmin = 0, rmax = 0;
        float rfmin = 0, rfmax = 0;
        int found = 0;
        for (int64_t i = 0; i < n; i++) {
            seed = seed * 1103515245u + 12345u;
            double v = (double)((int32_t)((seed >> 8) % 2001u) - 1000) / 7.0;
            if (((seed >> 4) % 6u) == 0) v = NAN;
            d[i] = v;
            f[i] = (float)v;
            if (isnan(v)) continue;
            if (!found) {
                rmin = rmax = v;
                rfmin = rfmax = (float)v;
                found = 1;
            } else {
                if (v < rmin) rmin = v;
                if (v > rmax) rmax = v;
                if ((float)v < rfmin) rfmin = (float)v;
                if ((float)v > rfmax) rfmax = (float)v;
            }
        }

        carquet_dispatch_minmax_double(d, n, &dmin, &dmax);
        if (!MINMAX_NAN_OK(dmin, dmax, rmin, rmax, found)) {
            TEST_FAIL("dispatch_minmax_nan", "minmax_double violates NaN contract");
        }
        carquet_dispatch_copy_minmax_double(d, n, dcopy, &dmin, &dmax);
        if (!MINMAX_NAN_OK(dmin, dmax, rmin, rmax, found) ||
            memcmp(d, dcopy, (size_t)n * sizeof(double)) != 0) {
            TEST_FAIL("dispatch_minmax_nan", "copy_minmax_double violates NaN contract");
        }
        carquet_dispatch_minmax_float(f, n, &fmin, &fmax);
        if (!MINMAX_NAN_OK(fmin, fmax, rfmin, rfmax, found)) {
            TEST_FAIL("dispatch_minmax_nan", "minmax_float violates NaN contract");
        }
        carquet_dispatch_copy_minmax_float(f, n, fcopy, &fmin, &fmax);
        if (!MINMAX_NAN_OK(fmin, fmax, rfmin, rfmax, found) ||
            memcmp(f, fcopy, (size_t)n * sizeof(float)) != 0) {
            TEST_FAIL("dispatch_minmax_nan", "copy_minmax_float violates NaN contract");
        }
    }

    TEST_PASS("dispatch_minmax_nan");
    return 0;
}

static int test_dispatch_null_helpers(void) {
    int16_t levels[] = {1, 0, 1, 1, 0, 0, 1, 0, 1, 1};
    int16_t filled[10];
    uint8_t bitmap[2] = {0, 0};

    assert(carquet_dispatch_count_non_nulls(levels, 10, 1) == 6);

    carquet_dispatch_build_null_bitmap(levels, 10, 1, bitmap);
    assert(bitmap[0] == 0x4D);
    assert(bitmap[1] == 0x03);

    carquet_dispatch_fill_def_levels(filled, 10, 3);
    for (int i = 0; i < 10; i++) {
        assert(filled[i] == 3);
    }

    TEST_PASS("dispatch_null_helpers");
    return 0;
}

static int test_bit_reader(void) {
    uint8_t data[] = {0xD2, 0xB4};  /* 0b11010010, 0b10110100, LSB first */
    carquet_bit_reader_t reader;
    carquet_bit_reader_init(&reader, data, 2);

    /* Read individual bits (LSB first) */
    assert(carquet_bit_reader_read_bit(&reader) == 0);
    assert(carquet_bit_reader_read_bit(&reader) == 1);
    assert(carquet_bit_reader_read_bit(&reader) == 0);
    assert(carquet_bit_reader_read_bit(&reader) == 0);

    /* Read 4 bits */
    uint32_t nibble = carquet_bit_reader_read_bits(&reader, 4);
    (void)nibble;
    assert(nibble == 0xD);  /* bits 4-7 of first byte: 0b1101 = 13 */

    TEST_PASS("bit_reader");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    int failures = 0;

    printf("=== Core Utility Tests ===\n\n");

    /* Arena tests */
    failures += test_arena_basic();
    failures += test_arena_large_allocation();
    failures += test_arena_save_restore();
    failures += test_arena_overflow_rejected();

    /* Buffer tests */
    failures += test_buffer_basic();
    failures += test_buffer_integers();
    failures += test_buffer_reader();
    failures += test_buffer_append_overflow_rejected();

    /* Endian tests */
    failures += test_endian_read_write();
    failures += test_varint();
    failures += test_zigzag();

    /* Bitpack tests */
    failures += test_bitpack_1bit();
    failures += test_bitpack_4bit();
    failures += test_bitpack_roundtrip();
    failures += test_bitpack_16bit();
    failures += test_run_detection();
    failures += test_dispatch_minmax();
    failures += test_dispatch_copy_minmax();
    failures += test_dispatch_minmax_nan();
    failures += test_dispatch_null_helpers();
    failures += test_bit_reader();

    printf("\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
