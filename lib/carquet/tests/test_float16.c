/**
 * @file test_float16.c
 * @brief Unit tests for IEEE 754 binary16 (half) -> binary32 conversion.
 *
 * carquet_half_to_float() (src/core/float16.h) underpins FLOAT16 column
 * statistics, which the Parquet spec orders by the represented floating-point
 * value (NaNs excluded), not lexicographically. These tests pin down the exact
 * bit-pattern -> value mapping across every class of half-precision number:
 * signed zero, subnormals, normals, the largest finite value, and the inf/NaN
 * special cases. If the conversion ever regresses, min/max stats for FLOAT16
 * columns would silently corrupt, so this is a pure-function regression guard.
 */

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "core/float16.h"
#include "test_helpers.h"

/* Bit-exact equality via the integer representation, so that -0.0 vs +0.0 and
 * distinct NaN payloads compare correctly (plain == treats -0.0 == +0.0 and
 * NaN != NaN). */
static int bits_equal(float a, float b) {
    uint32_t ba, bb;
    memcpy(&ba, &a, 4);
    memcpy(&bb, &b, 4);
    return ba == bb;
}

static int close_enough(float a, float b) {
    if (a == b) return 1;
    float scale = fabsf(a) > fabsf(b) ? fabsf(a) : fabsf(b);
    if (scale == 0.0f) return 1;
    return fabsf(a - b) <= scale * 1e-6f;
}

static int test_signed_zero(void) {
    float pz = carquet_half_to_float(0x0000);
    float nz = carquet_half_to_float(0x8000);
    /* Both are numerically zero but carry distinct sign bits. */
    if (pz != 0.0f || nz != 0.0f) TEST_FAIL("signed_zero", "not zero-valued");
    if (!bits_equal(pz, 0.0f)) TEST_FAIL("signed_zero", "+0 wrong sign bit");
    if (!bits_equal(nz, -0.0f)) TEST_FAIL("signed_zero", "-0 wrong sign bit");
    if (bits_equal(pz, nz)) TEST_FAIL("signed_zero", "+0 and -0 collapsed");
    TEST_PASS("signed_zero");
    return 0;
}

static int test_exact_normals(void) {
    struct { uint16_t h; float f; } cases[] = {
        {0x3C00,  1.0f},
        {0xBC00, -1.0f},
        {0x4000,  2.0f},
        {0xC000, -2.0f},
        {0x4200,  3.0f},
        {0x3800,  0.5f},
        {0x3555,  0.333251953125f},   /* nearest half to 1/3 */
        {0x4900, 10.0f},
        {0x6400, 1024.0f},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        float got = carquet_half_to_float(cases[i].h);
        if (!bits_equal(got, cases[i].f)) {
            printf("  half=0x%04X got=%.9g want=%.9g\n", cases[i].h, got, cases[i].f);
            TEST_FAIL("exact_normals", "value mismatch");
        }
    }
    TEST_PASS("exact_normals");
    return 0;
}

static int test_largest_finite(void) {
    /* 0x7BFF is the largest finite half: 65504. */
    float got = carquet_half_to_float(0x7BFF);
    if (!bits_equal(got, 65504.0f)) TEST_FAIL("largest_finite", "max half != 65504");
    float ngot = carquet_half_to_float(0xFBFF);
    if (!bits_equal(ngot, -65504.0f)) TEST_FAIL("largest_finite", "min half != -65504");
    TEST_PASS("largest_finite");
    return 0;
}

static int test_subnormals(void) {
    /* Smallest positive subnormal: 2^-24. */
    float smallest = carquet_half_to_float(0x0001);
    if (!close_enough(smallest, ldexpf(1.0f, -24))) TEST_FAIL("subnormals", "smallest subnormal wrong");
    if (smallest <= 0.0f) TEST_FAIL("subnormals", "smallest subnormal not positive");

    /* Largest subnormal: (1023/1024) * 2^-14, just below smallest normal. */
    float largest_sub = carquet_half_to_float(0x03FF);
    float smallest_norm = carquet_half_to_float(0x0400);  /* 2^-14 */
    if (!close_enough(smallest_norm, ldexpf(1.0f, -14))) TEST_FAIL("subnormals", "smallest normal wrong");
    if (!(largest_sub < smallest_norm)) TEST_FAIL("subnormals", "subnormal >= normal boundary");
    if (!(largest_sub > 0.0f)) TEST_FAIL("subnormals", "largest subnormal not positive");

    /* Negative subnormal carries the sign. */
    float neg_sub = carquet_half_to_float(0x8001);
    if (!(neg_sub < 0.0f)) TEST_FAIL("subnormals", "negative subnormal not negative");
    if (!close_enough(neg_sub, -smallest)) TEST_FAIL("subnormals", "negative subnormal magnitude wrong");
    TEST_PASS("subnormals");
    return 0;
}

static int test_infinities(void) {
    float pinf = carquet_half_to_float(0x7C00);
    float ninf = carquet_half_to_float(0xFC00);
    if (!isinf(pinf) || pinf <= 0.0f) TEST_FAIL("infinities", "not +inf");
    if (!isinf(ninf) || ninf >= 0.0f) TEST_FAIL("infinities", "not -inf");
    TEST_PASS("infinities");
    return 0;
}

static int test_nan(void) {
    /* Any half with exp==0x1F and nonzero mantissa is NaN. */
    float qnan = carquet_half_to_float(0x7E00);
    float snan = carquet_half_to_float(0x7C01);
    float nnan = carquet_half_to_float(0xFE00);
    if (!isnan(qnan)) TEST_FAIL("nan", "quiet NaN not NaN");
    if (!isnan(snan)) TEST_FAIL("nan", "signaling NaN not NaN");
    if (!isnan(nnan)) TEST_FAIL("nan", "negative NaN not NaN");
    TEST_PASS("nan");
    return 0;
}

static int test_monotonic_ordering(void) {
    /* For non-negative, non-NaN halves, increasing bit pattern is increasing
     * value. This is the property FLOAT16 statistics ordering relies on when
     * comparing by represented value. Sweep the finite positive range. */
    float prev = carquet_half_to_float(0x0000);
    for (uint16_t h = 0x0001; h < 0x7C00; h++) {  /* stop before +inf */
        float cur = carquet_half_to_float(h);
        if (!(cur > prev)) {
            printf("  non-monotonic at 0x%04X: prev=%.9g cur=%.9g\n", h, prev, cur);
            TEST_FAIL("monotonic_ordering", "not strictly increasing");
        }
        prev = cur;
    }
    TEST_PASS("monotonic_ordering");
    return 0;
}

static int test_roundtrip_from_float(void) {
    /* Encode a handful of floats to half by hand, decode, and confirm they
     * survive when exactly representable. */
    struct { uint16_t h; float f; } exact[] = {
        {0x0000, 0.0f}, {0x3C00, 1.0f}, {0xC500, -5.0f}, {0x4980, 11.0f},
    };
    for (size_t i = 0; i < sizeof(exact) / sizeof(exact[0]); i++) {
        if (!bits_equal(carquet_half_to_float(exact[i].h), exact[i].f))
            TEST_FAIL("roundtrip_from_float", "exact value not preserved");
    }
    TEST_PASS("roundtrip_from_float");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_signed_zero();
    failures += test_exact_normals();
    failures += test_largest_finite();
    failures += test_subnormals();
    failures += test_infinities();
    failures += test_nan();
    failures += test_monotonic_ordering();
    failures += test_roundtrip_from_float();
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll float16 tests passed\n");
    return 0;
}
