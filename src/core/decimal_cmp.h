/**
 * @file decimal_cmp.h
 * @brief Signed big-endian two's-complement comparison for Parquet DECIMAL.
 *
 * Kept separate from endian.h on purpose: this header intentionally does NOT
 * define CARQUET_LITTLE_ENDIAN. page_writer.c gates a SIMD encode fast path on
 * `#if CARQUET_LITTLE_ENDIAN`, so pulling in endian.h there would silently flip
 * which branch compiles. Callers that only need the DECIMAL comparator include
 * this instead.
 */

#ifndef CARQUET_CORE_DECIMAL_CMP_H
#define CARQUET_CORE_DECIMAL_CMP_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compare two big-endian two's-complement integers of possibly different byte
 * lengths, as used by Parquet DECIMAL columns backed by FIXED_LEN_BYTE_ARRAY or
 * BYTE_ARRAY (whose sort order is SIGNED, not unsigned-lexicographic). Returns
 * a value < 0, 0, or > 0 for a < b, a == b, a > b.
 *
 * A negative value (sign bit of the most-significant byte set) always orders
 * below a non-negative one; two same-sign values compare by sign-extending the
 * shorter operand to the longer width (0x00 for non-negative, 0xFF for
 * negative). A zero-length operand is treated as the value 0.
 */
static inline int carquet_compare_decimal_be(const uint8_t* a, size_t alen,
                                             const uint8_t* b, size_t blen) {
    int a_neg = (alen > 0) && (a[0] & 0x80);
    int b_neg = (blen > 0) && (b[0] & 0x80);
    if (a_neg != b_neg) return a_neg ? -1 : 1;  /* negative < non-negative */

    /* Same sign: the shorter operand is virtually sign-extended on the left. */
    uint8_t ext = a_neg ? 0xFFu : 0x00u;
    size_t n = alen > blen ? alen : blen;
    size_t pa = n - alen;  /* count of leading extension bytes for a */
    size_t pb = n - blen;
    for (size_t i = 0; i < n; i++) {
        uint8_t av = (i < pa) ? ext : a[i - pa];
        uint8_t bv = (i < pb) ? ext : b[i - pb];
        if (av != bv) return av < bv ? -1 : 1;
    }
    return 0;
}

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_CORE_DECIMAL_CMP_H */
