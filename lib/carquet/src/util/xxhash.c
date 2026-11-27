/**
 * @file xxhash.c
 * @brief xxHash implementation for bloom filters
 *
 * xxHash is a fast non-cryptographic hash algorithm.
 * Parquet uses xxHash64 for bloom filter hashing.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* xxHash64 constants */
#define XXH_PRIME64_1 0x9E3779B185EBCA87ULL
#define XXH_PRIME64_2 0xC2B2AE3D27D4EB4FULL
#define XXH_PRIME64_3 0x165667B19E3779F9ULL
#define XXH_PRIME64_4 0x85EBCA77C2B2AE63ULL
#define XXH_PRIME64_5 0x27D4EB2F165667C5ULL

static inline uint64_t xxh64_rotl(uint64_t x, int r) {
    return (x << r) | (x >> (64 - r));
}

static inline uint64_t xxh64_round(uint64_t acc, uint64_t input) {
    acc += input * XXH_PRIME64_2;
    acc = xxh64_rotl(acc, 31);
    acc *= XXH_PRIME64_1;
    return acc;
}

static inline uint64_t xxh64_merge_round(uint64_t acc, uint64_t val) {
    val = xxh64_round(0, val);
    acc ^= val;
    acc = acc * XXH_PRIME64_1 + XXH_PRIME64_4;
    return acc;
}

static inline uint64_t read64_le(const uint8_t* p) {
    return (uint64_t)p[0] |
           ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) |
           ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) |
           ((uint64_t)p[7] << 56);
}

static inline uint32_t read32_le(const uint8_t* p) {
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

uint64_t carquet_xxhash64(const void* data, size_t length, uint64_t seed) {
    const uint8_t* p = (const uint8_t*)data;
    const uint8_t* end = p + length;
    uint64_t h64;

    if (length >= 32) {
        const uint8_t* limit = end - 32;
        uint64_t v1 = seed + XXH_PRIME64_1 + XXH_PRIME64_2;
        uint64_t v2 = seed + XXH_PRIME64_2;
        uint64_t v3 = seed + 0;
        uint64_t v4 = seed - XXH_PRIME64_1;

        do {
            v1 = xxh64_round(v1, read64_le(p)); p += 8;
            v2 = xxh64_round(v2, read64_le(p)); p += 8;
            v3 = xxh64_round(v3, read64_le(p)); p += 8;
            v4 = xxh64_round(v4, read64_le(p)); p += 8;
        } while (p <= limit);

        h64 = xxh64_rotl(v1, 1) + xxh64_rotl(v2, 7) +
              xxh64_rotl(v3, 12) + xxh64_rotl(v4, 18);

        h64 = xxh64_merge_round(h64, v1);
        h64 = xxh64_merge_round(h64, v2);
        h64 = xxh64_merge_round(h64, v3);
        h64 = xxh64_merge_round(h64, v4);
    } else {
        h64 = seed + XXH_PRIME64_5;
    }

    h64 += (uint64_t)length;

    /* Process remaining 8-byte chunks */
    while (p + 8 <= end) {
        uint64_t k1 = xxh64_round(0, read64_le(p));
        h64 ^= k1;
        h64 = xxh64_rotl(h64, 27) * XXH_PRIME64_1 + XXH_PRIME64_4;
        p += 8;
    }

    /* Process remaining 4 bytes */
    if (p + 4 <= end) {
        h64 ^= (uint64_t)read32_le(p) * XXH_PRIME64_1;
        h64 = xxh64_rotl(h64, 23) * XXH_PRIME64_2 + XXH_PRIME64_3;
        p += 4;
    }

    /* Process remaining bytes */
    while (p < end) {
        h64 ^= (uint64_t)(*p) * XXH_PRIME64_5;
        h64 = xxh64_rotl(h64, 11) * XXH_PRIME64_1;
        p++;
    }

    /* Final mix */
    h64 ^= h64 >> 33;
    h64 *= XXH_PRIME64_2;
    h64 ^= h64 >> 29;
    h64 *= XXH_PRIME64_3;
    h64 ^= h64 >> 32;

    return h64;
}
