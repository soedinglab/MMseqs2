/**
 * @file float16.h
 * @brief IEEE 754 binary16 (half) -> binary32 conversion.
 *
 * Used for FLOAT16 column statistics, which the Parquet spec orders by the
 * represented floating-point value (NaNs excluded), not lexicographically.
 */
#ifndef CARQUET_FLOAT16_H
#define CARQUET_FLOAT16_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline float carquet_half_to_float(uint16_t h) {
    uint32_t sign = (uint32_t)(h >> 15) & 1u;
    uint32_t exp  = (h >> 10) & 0x1Fu;
    uint32_t mant = h & 0x3FFu;
    uint32_t f;
    if (exp == 0) {
        if (mant == 0) {
            f = sign << 31;
        } else {
            exp = 1;
            while ((mant & 0x400u) == 0) { mant <<= 1; exp--; }
            mant &= 0x3FFu;
            f = (sign << 31) | ((exp + (127 - 15)) << 23) | (mant << 13);
        }
    } else if (exp == 0x1Fu) {
        f = (sign << 31) | (0xFFu << 23) | (mant << 13);
    } else {
        f = (sign << 31) | ((exp + (127 - 15)) << 23) | (mant << 13);
    }
    float out;
    memcpy(&out, &f, sizeof(out));
    return out;
}

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_FLOAT16_H */
