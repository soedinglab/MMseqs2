// Tests for the varint / zigzag / fixed-width codec the k-mer and candidate-edge
// bucket formats are built on.
//
// Two properties matter, and they fail in opposite directions:
//
//   - round-trip exactness. These encodings sit under the two largest
//     intermediates in the pipeline, and a value that decodes to something else
//     does not crash -- it produces a k-mer that groups with the wrong partners,
//     or an edge on the wrong diagonal. That is a silently different clustering,
//     which is the failure mode this whole project is least able to detect.
//   - clean rejection of a truncated or corrupt encoding. A worker killed
//     mid-flush leaves a torn tail, and readers are required to round down to
//     whole records and carry on. That is only safe if the decoder reliably says
//     "no" instead of reading past the buffer.

#include "VarintCodec.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

const char* binary_name = "test_varintcodec";

static int failures = 0;

static void check(bool condition, const std::string &what) {
    if (condition == false) {
        fprintf(stderr, "FAIL: %s\n", what.c_str());
        failures++;
    } else {
        fprintf(stdout, "ok:   %s\n", what.c_str());
    }
}

// A cheap deterministic generator, so a failure is reproducible without carrying
// a seed corpus around.
static uint64_t nextRandom(uint64_t &state) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static void testUnsignedRoundTrip() {
    std::vector<uint64_t> values;
    // Every power-of-two boundary, and both sides of it: the byte count changes
    // at 2^7, 2^14, ... and off-by-one there is the classic varint defect.
    for (unsigned int bit = 0; bit < 64; bit++) {
        const uint64_t v = static_cast<uint64_t>(1) << bit;
        values.push_back(v - 1);
        values.push_back(v);
        values.push_back(v + 1);
    }
    values.push_back(0);
    values.push_back(UINT64_MAX);
    values.push_back(UINT64_MAX - 1);
    uint64_t state = 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 20000; i++) {
        values.push_back(nextRandom(state));
        // Also small values, which the random stream almost never produces and
        // which are the common case in every real field.
        values.push_back(nextRandom(state) % 256);
    }

    bool allOk = true;
    bool sizeOk = true;
    for (size_t i = 0; i < values.size(); i++) {
        uint8_t buffer[VarintCodec::MAX_BYTES];
        uint8_t *end = VarintCodec::write(buffer, values[i]);
        const size_t written = static_cast<size_t>(end - buffer);
        if (written != VarintCodec::size(values[i]) || written > VarintCodec::MAX_BYTES) {
            sizeOk = false;
        }
        const uint8_t *in = buffer;
        uint64_t decoded = 0;
        if (VarintCodec::read(in, end, decoded) == false || decoded != values[i] ||
            in != end) {
            allOk = false;
        }
    }
    check(allOk, "varint round-trips every boundary and 40000 random values");
    check(sizeOk, "size() agrees with write() and never exceeds MAX_BYTES");
}

static void testZigzagRoundTrip() {
    std::vector<int64_t> values;
    values.push_back(0);
    values.push_back(1);
    values.push_back(-1);
    values.push_back(INT64_MAX);
    values.push_back(INT64_MIN);
    // The diagonal field is an int16_t, so its whole range is worth covering
    // exhaustively rather than sampling.
    for (int v = -32768; v <= 32767; v++) {
        values.push_back(v);
    }

    bool allOk = true;
    bool smallOk = true;
    for (size_t i = 0; i < values.size(); i++) {
        const uint64_t mapped = VarintCodec::zigzag(values[i]);
        if (VarintCodec::unzigzag(mapped) != values[i]) {
            allOk = false;
        }
        // The point of zigzag: a small negative must stay short. Without it every
        // negative diagonal would encode as ten bytes.
        if (values[i] >= -63 && values[i] <= 63 && VarintCodec::size(mapped) != 1) {
            smallOk = false;
        }
    }
    check(allOk, "zigzag round-trips the full int16 range and the int64 extremes");
    check(smallOk, "zigzag keeps small-magnitude negatives to one varint byte");
}

static void testTruncationRejected() {
    // A ten-byte value, then every proper prefix of it. Each must be refused
    // rather than decoded as some shorter number.
    uint8_t buffer[VarintCodec::MAX_BYTES];
    uint8_t *end = VarintCodec::write(buffer, UINT64_MAX);
    const size_t full = static_cast<size_t>(end - buffer);

    bool allRejected = true;
    for (size_t cut = 0; cut < full; cut++) {
        const uint8_t *in = buffer;
        uint64_t decoded = 0;
        if (VarintCodec::read(in, buffer + cut, decoded)) {
            allRejected = false;
        }
        // A rejected read must not consume anything, or a caller that retries
        // with more data would resume from the wrong offset.
        if (in != buffer) {
            allRejected = false;
        }
    }
    check(allRejected, "every truncated prefix of a 10-byte varint is rejected, cursor unmoved");

    // An encoding that continues past 64 bits is corrupt, not a big number.
    uint8_t overlong[12];
    for (size_t i = 0; i < 11; i++) {
        overlong[i] = 0xFF;
    }
    overlong[11] = 0x01;
    const uint8_t *in = overlong;
    uint64_t decoded = 0;
    check(VarintCodec::read(in, overlong + 12, decoded) == false,
          "an over-long varint is rejected rather than silently truncated");

    // The tenth group may carry only bit 63.
    uint8_t tenth[VarintCodec::MAX_BYTES];
    for (size_t i = 0; i < 9; i++) {
        tenth[i] = 0x80;
    }
    tenth[9] = 0x02;  // bit 1 of the tenth group is past bit 64
    const uint8_t *in2 = tenth;
    check(VarintCodec::read(in2, tenth + VarintCodec::MAX_BYTES, decoded) == false,
          "a tenth group carrying more than one bit is rejected");
}

static void testFixedWidth() {
    bool widthOk = true;
    widthOk = widthOk && VarintCodec::fixedWidthFor(0) == 1;
    widthOk = widthOk && VarintCodec::fixedWidthFor(255) == 1;
    widthOk = widthOk && VarintCodec::fixedWidthFor(256) == 2;
    // The two cases the k-mer field actually takes: 13^14 at --min-seq-id 0.9 and
    // 21^14 at >= 0.99. Seven bytes and eight respectively -- the whole reason the
    // width is derived rather than fixed at 8.
    uint64_t base13 = 1;
    uint64_t base21 = 1;
    for (int i = 0; i < 14; i++) {
        base13 *= 13;
        base21 *= 21;
    }
    widthOk = widthOk && VarintCodec::fixedWidthFor(base13 - 1) == 7;
    widthOk = widthOk && VarintCodec::fixedWidthFor(base21 - 1) == 8;
    check(widthOk, "fixedWidthFor picks 7 bytes for a base-13 k=14 k-mer and 8 for base-21");

    bool roundTripOk = true;
    bool truncationOk = true;
    uint64_t state = 0xDEADBEEFCAFEBABEULL;
    for (unsigned int width = 1; width <= 8; width++) {
        const uint64_t mask =
            (width == 8) ? UINT64_MAX : ((static_cast<uint64_t>(1) << (width * 8)) - 1);
        for (int i = 0; i < 2000; i++) {
            const uint64_t value = nextRandom(state) & mask;
            uint8_t buffer[8];
            uint8_t *end = VarintCodec::writeFixed(buffer, value, width);
            if (static_cast<size_t>(end - buffer) != width) {
                roundTripOk = false;
            }
            const uint8_t *in = buffer;
            uint64_t decoded = 0;
            if (VarintCodec::readFixed(in, end, width, decoded) == false || decoded != value ||
                in != end) {
                roundTripOk = false;
            }
            // One byte short must be refused.
            const uint8_t *shortIn = buffer;
            uint64_t ignored = 0;
            if (VarintCodec::readFixed(shortIn, buffer + width - 1, width, ignored)) {
                truncationOk = false;
            }
        }
    }
    check(roundTripOk, "fixed-width round-trips every width from 1 to 8 bytes");
    check(truncationOk, "a short fixed-width field is rejected");
}

int main(int, const char**) {
    testUnsignedRoundTrip();
    testZigzagRoundTrip();
    testTruncationRejected();
    testFixedWidth();

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
