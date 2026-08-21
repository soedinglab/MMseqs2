#ifndef MMSEQS_VARINTCODEC_H
#define MMSEQS_VARINTCODEC_H

#include <cstddef>
#include <cstdint>

// LEB128 varints, zigzag mapping and narrow fixed-width integers, for the k-mer
// and candidate-edge bucket formats.
//
// Why those intermediates are encoded at all. They are the two largest things
// this pipeline writes, and at the target scales they are what decides whether a
// run fits its scratch filesystem. With fixed-width records, 1e12 sequences give
// ~504 TB of k-mer shuffle (21 x 24 B) and several hundred TB of candidate edges
// (~22 x 17 B measured at 1e9) against a budget of roughly 1 PB. Waves divide the
// shuffle but not the edges, because alignment runs once after every wave.
//
// Both formats are written once and read back sequentially exactly once, which is
// precisely the case where a variable-width encoding costs nothing but a little
// CPU: nothing seeks into them, nothing updates them in place.
//
// LEB128 rather than group- or stream-vbyte: the fields here are heterogeneous (a
// packed k-mer, an ascending id delta, a position, a zigzagged diagonal), so they
// do not form the uniform four-lane groups those layouts need, and the branchy
// decode is far from the bottleneck next to the sort and the alignment it feeds.
namespace VarintCodec {

// A 64-bit value is at most ten 7-bit groups.
const size_t MAX_BYTES = 10;

// Bytes write() will emit for this value. Used to size buffers up front, so the
// encoders never have to test for capacity per field.
inline size_t size(uint64_t value) {
    size_t n = 1;
    while (value >= 0x80) {
        value >>= 7;
        n++;
    }
    return n;
}

// Appends value and returns the new write position. The caller guarantees
// MAX_BYTES of room -- every producer here writes into a buffer sized with size()
// or with a MAX_BYTES-per-field bound.
inline uint8_t *write(uint8_t *out, uint64_t value) {
    while (value >= 0x80) {
        *out++ = static_cast<uint8_t>((value & 0x7F) | 0x80);
        value >>= 7;
    }
    *out++ = static_cast<uint8_t>(value);
    return out;
}

// Bounds-checked decode. Returns false on a truncated or over-long encoding and
// leaves `in` untouched, so a torn block tail is a clean failure rather than a
// read past the buffer.
//
// This is not defensive programming for its own sake: a worker killed mid-flush
// leaves exactly such a tail, and the readers are required to round down to whole
// records and carry on rather than treat the partition as poisoned.
inline bool read(const uint8_t *&in, const uint8_t *end, uint64_t &value) {
    uint64_t result = 0;
    unsigned int shift = 0;
    const uint8_t *p = in;
    while (p < end) {
        const uint8_t byte = *p++;
        // The tenth group carries only the single remaining bit; anything wider
        // is a corrupt encoding, not a large number.
        if (shift == 63 && (byte & 0xFE) != 0) {
            return false;
        }
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) {
            value = result;
            in = p;
            return true;
        }
        shift += 7;
        if (shift > 63) {
            return false;
        }
    }
    return false;
}

// Maps a signed value onto an unsigned one that is small when the input is small
// in magnitude, so a varint stays short for both signs. Diagonals are signed and
// cluster near zero; without this every negative one would encode as ten bytes.
inline uint64_t zigzag(int64_t value) {
    return (static_cast<uint64_t>(value) << 1) ^ static_cast<uint64_t>(-(value < 0 ? 1 : 0));
}

inline int64_t unzigzag(uint64_t value) {
    return static_cast<int64_t>(value >> 1) ^ -static_cast<int64_t>(value & 1);
}

// Narrowest byte count that can hold every value up to maxValue. The k-mer field
// uses this: at --min-seq-id 0.9 the index is base-13 over k = 14, so 13^14 is
// 2^51.8 and seven bytes are enough, while --min-seq-id >= 0.99 switches to the
// 21-letter alphabet (2^61.5) and needs eight.
inline unsigned int fixedWidthFor(uint64_t maxValue) {
    unsigned int width = 1;
    while (width < 8 && (maxValue >> (width * 8)) != 0) {
        width++;
    }
    return width;
}

// Little-endian fixed width. Used where the value has a known ceiling but no
// useful skew towards small numbers, so a varint would only add a continuation
// bit per byte.
inline uint8_t *writeFixed(uint8_t *out, uint64_t value, unsigned int width) {
    for (unsigned int i = 0; i < width; i++) {
        *out++ = static_cast<uint8_t>(value & 0xFF);
        value >>= 8;
    }
    return out;
}

inline bool readFixed(const uint8_t *&in, const uint8_t *end, unsigned int width, uint64_t &value) {
    if (static_cast<size_t>(end - in) < width) {
        return false;
    }
    uint64_t result = 0;
    for (unsigned int i = 0; i < width; i++) {
        result |= static_cast<uint64_t>(in[i]) << (i * 8);
    }
    in += width;
    value = result;
    return true;
}

}  // namespace VarintCodec

#endif
