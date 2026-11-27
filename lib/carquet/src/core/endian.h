/**
 * @file endian.h
 * @brief Endianness handling utilities
 *
 * Parquet uses little-endian byte order for all multi-byte values.
 * These utilities handle reading and writing values in the correct byte order.
 */

#ifndef CARQUET_CORE_ENDIAN_H
#define CARQUET_CORE_ENDIAN_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Endianness Detection
 * ============================================================================
 */

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__)
    /* GCC/Clang - most reliable detection */
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        #define CARQUET_LITTLE_ENDIAN 1
    #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        #define CARQUET_LITTLE_ENDIAN 0
    #else
        #error "Unknown byte order"
    #endif
#elif defined(__LITTLE_ENDIAN__) || defined(__ARMEL__) || defined(__THUMBEL__) || \
      defined(__AARCH64EL__) || defined(_MIPSEL) || defined(__MIPSEL) || \
      defined(__MIPSEL__) || defined(__XTENSA_EL__) || defined(__RISCV__) || \
      defined(_WIN32) || defined(__x86_64__) || defined(__i386__) || \
      defined(__amd64__)
    /* Explicitly little-endian platforms */
    #define CARQUET_LITTLE_ENDIAN 1
#elif defined(__BIG_ENDIAN__) || defined(__ARMEB__) || defined(__THUMBEB__) || \
      defined(__AARCH64EB__) || defined(_MIPSEB) || defined(__MIPSEB) || \
      defined(__MIPSEB__) || defined(__XTENSA_EB__) || defined(__sparc__) || \
      defined(__s390__) || defined(__s390x__) || defined(__hppa__) || \
      defined(__HPPA__) || defined(__powerpc__) || defined(__ppc__) || \
      defined(__PPC__) || defined(_POWER)
    /* Explicitly big-endian platforms */
    #define CARQUET_LITTLE_ENDIAN 0
#else
    /* Fallback: use runtime check or compilation will fail if incorrect */
    #warning "Unknown endianness - assuming little-endian. Define CARQUET_LITTLE_ENDIAN=0 for big-endian."
    #define CARQUET_LITTLE_ENDIAN 1
#endif

/* ============================================================================
 * Byte Swap Intrinsics
 * ============================================================================
 */

#if defined(__GNUC__) || defined(__clang__)
    #define carquet_bswap16(x) __builtin_bswap16(x)
    #define carquet_bswap32(x) __builtin_bswap32(x)
    #define carquet_bswap64(x) __builtin_bswap64(x)
#elif defined(_MSC_VER)
    #include <stdlib.h>
    #define carquet_bswap16(x) _byteswap_ushort(x)
    #define carquet_bswap32(x) _byteswap_ulong(x)
    #define carquet_bswap64(x) _byteswap_uint64(x)
#else
    static inline uint16_t carquet_bswap16(uint16_t x) {
        return (x >> 8) | (x << 8);
    }
    static inline uint32_t carquet_bswap32(uint32_t x) {
        return ((x >> 24) & 0x000000FF) |
               ((x >> 8)  & 0x0000FF00) |
               ((x << 8)  & 0x00FF0000) |
               ((x << 24) & 0xFF000000);
    }
    static inline uint64_t carquet_bswap64(uint64_t x) {
        return ((x >> 56) & 0x00000000000000FFULL) |
               ((x >> 40) & 0x000000000000FF00ULL) |
               ((x >> 24) & 0x0000000000FF0000ULL) |
               ((x >> 8)  & 0x00000000FF000000ULL) |
               ((x << 8)  & 0x000000FF00000000ULL) |
               ((x << 24) & 0x0000FF0000000000ULL) |
               ((x << 40) & 0x00FF000000000000ULL) |
               ((x << 56) & 0xFF00000000000000ULL);
    }
#endif

/* ============================================================================
 * Little-Endian Read Functions
 * ============================================================================
 */

/**
 * Read a 16-bit unsigned integer from little-endian bytes.
 */
static inline uint16_t carquet_read_u16_le(const uint8_t* p) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    uint16_t v;
    memcpy(&v, p, sizeof(v));
    return v;
#else
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
#endif
}

/**
 * Read a 32-bit unsigned integer from little-endian bytes.
 */
static inline uint32_t carquet_read_u32_le(const uint8_t* p) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    uint32_t v;
    memcpy(&v, p, sizeof(v));
    return v;
#else
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
#endif
}

/**
 * Read a 64-bit unsigned integer from little-endian bytes.
 */
static inline uint64_t carquet_read_u64_le(const uint8_t* p) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    uint64_t v;
    memcpy(&v, p, sizeof(v));
    return v;
#else
    return (uint64_t)p[0] |
           ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) |
           ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) |
           ((uint64_t)p[7] << 56);
#endif
}

/**
 * Read a 16-bit signed integer from little-endian bytes.
 */
static inline int16_t carquet_read_i16_le(const uint8_t* p) {
    return (int16_t)carquet_read_u16_le(p);
}

/**
 * Read a 32-bit signed integer from little-endian bytes.
 */
static inline int32_t carquet_read_i32_le(const uint8_t* p) {
    return (int32_t)carquet_read_u32_le(p);
}

/**
 * Read a 64-bit signed integer from little-endian bytes.
 */
static inline int64_t carquet_read_i64_le(const uint8_t* p) {
    return (int64_t)carquet_read_u64_le(p);
}

/**
 * Read a 32-bit float from little-endian bytes.
 */
static inline float carquet_read_f32_le(const uint8_t* p) {
    uint32_t bits = carquet_read_u32_le(p);
    float f;
    memcpy(&f, &bits, sizeof(f));
    return f;
}

/**
 * Read a 64-bit double from little-endian bytes.
 */
static inline double carquet_read_f64_le(const uint8_t* p) {
    uint64_t bits = carquet_read_u64_le(p);
    double d;
    memcpy(&d, &bits, sizeof(d));
    return d;
}

/* ============================================================================
 * Little-Endian Write Functions
 * ============================================================================
 */

/**
 * Write a 16-bit unsigned integer as little-endian bytes.
 */
static inline void carquet_write_u16_le(uint8_t* p, uint16_t v) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(p, &v, sizeof(v));
#else
    p[0] = (uint8_t)(v);
    p[1] = (uint8_t)(v >> 8);
#endif
}

/**
 * Write a 32-bit unsigned integer as little-endian bytes.
 */
static inline void carquet_write_u32_le(uint8_t* p, uint32_t v) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(p, &v, sizeof(v));
#else
    p[0] = (uint8_t)(v);
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
#endif
}

/**
 * Write a 64-bit unsigned integer as little-endian bytes.
 */
static inline void carquet_write_u64_le(uint8_t* p, uint64_t v) {
#if CARQUET_LITTLE_ENDIAN && !defined(CARQUET_STRICT_ALIGN)
    memcpy(p, &v, sizeof(v));
#else
    p[0] = (uint8_t)(v);
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
    p[4] = (uint8_t)(v >> 32);
    p[5] = (uint8_t)(v >> 40);
    p[6] = (uint8_t)(v >> 48);
    p[7] = (uint8_t)(v >> 56);
#endif
}

/**
 * Write a 16-bit signed integer as little-endian bytes.
 */
static inline void carquet_write_i16_le(uint8_t* p, int16_t v) {
    carquet_write_u16_le(p, (uint16_t)v);
}

/**
 * Write a 32-bit signed integer as little-endian bytes.
 */
static inline void carquet_write_i32_le(uint8_t* p, int32_t v) {
    carquet_write_u32_le(p, (uint32_t)v);
}

/**
 * Write a 64-bit signed integer as little-endian bytes.
 */
static inline void carquet_write_i64_le(uint8_t* p, int64_t v) {
    carquet_write_u64_le(p, (uint64_t)v);
}

/**
 * Write a 32-bit float as little-endian bytes.
 */
static inline void carquet_write_f32_le(uint8_t* p, float f) {
    uint32_t bits;
    memcpy(&bits, &f, sizeof(bits));
    carquet_write_u32_le(p, bits);
}

/**
 * Write a 64-bit double as little-endian bytes.
 */
static inline void carquet_write_f64_le(uint8_t* p, double d) {
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));
    carquet_write_u64_le(p, bits);
}

/* ============================================================================
 * Varint Encoding (for Thrift)
 * ============================================================================
 */

/**
 * Encode a 32-bit unsigned integer as a varint.
 * Returns number of bytes written (1-5).
 */
static inline int carquet_encode_varint32(uint8_t* p, uint32_t v) {
    int i = 0;
    while (v >= 0x80) {
        p[i++] = (uint8_t)((v & 0x7F) | 0x80);
        v >>= 7;
    }
    p[i++] = (uint8_t)v;
    return i;
}

/**
 * Encode a 64-bit unsigned integer as a varint.
 * Returns number of bytes written (1-10).
 */
static inline int carquet_encode_varint64(uint8_t* p, uint64_t v) {
    int i = 0;
    while (v >= 0x80) {
        p[i++] = (uint8_t)((v & 0x7F) | 0x80);
        v >>= 7;
    }
    p[i++] = (uint8_t)v;
    return i;
}

/**
 * Decode a varint32 from bytes.
 * Returns number of bytes consumed, or -1 on error.
 */
static inline int carquet_decode_varint32(const uint8_t* p, size_t len, uint32_t* out) {
    uint32_t result = 0;
    int shift = 0;
    size_t i = 0;

    while (i < len && i < 5) {
        uint8_t byte = p[i];
        result |= (uint32_t)(byte & 0x7F) << shift;

        if ((byte & 0x80) == 0) {
            *out = result;
            return (int)(i + 1);
        }

        shift += 7;
        i++;
    }

    return -1;  /* Truncated or overflow */
}

/**
 * Decode a varint64 from bytes.
 * Returns number of bytes consumed, or -1 on error.
 */
static inline int carquet_decode_varint64(const uint8_t* p, size_t len, uint64_t* out) {
    uint64_t result = 0;
    int shift = 0;
    size_t i = 0;

    while (i < len && i < 10) {
        uint8_t byte = p[i];
        result |= (uint64_t)(byte & 0x7F) << shift;

        if ((byte & 0x80) == 0) {
            *out = result;
            return (int)(i + 1);
        }

        shift += 7;
        i++;
    }

    return -1;  /* Truncated or overflow */
}

/**
 * Zigzag encode a signed 32-bit integer for varint encoding.
 */
static inline uint32_t carquet_zigzag_encode32(int32_t v) {
    return ((uint32_t)v << 1) ^ ((uint32_t)((int32_t)v >> 31));
}

/**
 * Zigzag encode a signed 64-bit integer for varint encoding.
 */
static inline uint64_t carquet_zigzag_encode64(int64_t v) {
    return ((uint64_t)v << 1) ^ ((uint64_t)((int64_t)v >> 63));
}

/**
 * Zigzag decode a 32-bit varint to signed integer.
 */
static inline int32_t carquet_zigzag_decode32(uint32_t v) {
    return (int32_t)((v >> 1) ^ (-(int32_t)(v & 1)));
}

/**
 * Zigzag decode a 64-bit varint to signed integer.
 */
static inline int64_t carquet_zigzag_decode64(uint64_t v) {
    return (int64_t)((v >> 1) ^ (-(int64_t)(v & 1)));
}

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_CORE_ENDIAN_H */
