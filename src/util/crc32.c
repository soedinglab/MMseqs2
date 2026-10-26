/**
 * @file crc32.c
 * @brief CRC32 checksum for Parquet page checksums
 *
 * Parquet page CRCs use the standard IEEE CRC32 polynomial (0x04C11DB7),
 * computed over the serialized page body and excluding the page header.
 *
 * Delegates to zlib's crc32, which is already linked for GZIP/DEFLATE
 * support and ships with hardware-accelerated implementations
 * (PCLMULQDQ folding on x86, FEAT_CRC32 on ARMv8).
 */

#include <stdint.h>
#include <stddef.h>

#include <zlib.h>

uint32_t carquet_crc32_update(uint32_t crc, const uint8_t* data, size_t length) {
    /* zlib's crc32 takes uInt (32-bit) lengths; chunk if needed. */
    uLong c = (uLong)crc;
    while (length > 0) {
        uInt chunk = length > (size_t)0xFFFFFFF0u ? (uInt)0xFFFFFFF0u : (uInt)length;
        c = crc32(c, (const Bytef*)data, chunk);
        data += chunk;
        length -= chunk;
    }
    return (uint32_t)c;
}

uint32_t carquet_crc32(const uint8_t* data, size_t length) {
    return carquet_crc32_update(0, data, length);
}
