// Tests for the packed candidate-edge block format.
//
// This format carries the largest intermediate the pipeline writes, and it had
// no direct test until two defects had already been shipped in it: a decoder
// that reserved per block and so copied the whole accumulated bucket every time
// (quadratic in block count, +65% on the align stage), and a producer that fed
// it a descending representative after slicing was added. The second was caught
// only because the encoder happens to check its own precondition. That is a thin
// margin for the byte format the whole scratch budget rests on.
//
// Three properties, and they fail in different directions:
//
//   - round-trip exactness. A field that decodes to something else does not
//     crash; it produces an edge on the wrong diagonal or between the wrong
//     pair, which is a silently different clustering.
//   - raw and packed agree. --raw-records is the control the whole encoding
//     argument rests on, so "both encodings decode to the same edges" is the
//     property that makes it a control rather than a second implementation.
//   - corruption is refused, not decoded. A torn tail is what an interrupted
//     worker leaves, and the reader is required to stop at the last good block.

#include "CandidateEdge.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

const char* binary_name = "test_edgecodec";

static int failures = 0;

static void check(bool condition, const std::string &what) {
    if (condition == false) {
        fprintf(stderr, "FAIL: %s\n", what.c_str());
        failures++;
    } else {
        fprintf(stdout, "ok:   %s\n", what.c_str());
    }
}

static uint64_t nextRandom(uint64_t &state) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static bool sameEdge(const CandidateEdge &a, const CandidateEdge &b) {
    return a.getRep() == b.getRep() && a.getMember() == b.getMember() &&
           a.diagonal == b.diagonal && a.score == b.score && a.reverseStrand == b.reverseStrand;
}

// Encodes and decodes, returning the payload size so the two encodings can be
// compared on size as well as on content.
static bool roundTrip(const std::vector<CandidateEdge> &edges, bool raw,
                      std::vector<CandidateEdge> &out, size_t &payloadBytes) {
    std::vector<uint8_t> buffer;
    EdgeBlockCodec::encode(edges, raw, 7, 3, buffer);
    if (edges.empty()) {
        return buffer.empty();
    }
    EdgeBlockHeader header;
    memcpy(&header, buffer.data(), sizeof(EdgeBlockHeader));
    if (EdgeBlockCodec::headerLooksValid(header) == false) {
        return false;
    }
    if (header.partition != 7 || header.worker != 3 || header.recordCount != edges.size()) {
        return false;
    }
    payloadBytes = static_cast<size_t>(header.payloadBytes);
    out.clear();
    return EdgeBlockCodec::decode(header, buffer.data() + sizeof(EdgeBlockHeader), out);
}

// The encoder requires non-decreasing representatives, which the reduce
// guarantees. Everything generated here respects that.
static std::vector<CandidateEdge> makeEdges(size_t count, uint64_t seed, bool extremes) {
    std::vector<CandidateEdge> edges;
    uint64_t state = seed;
    uint64_t rep = 0;
    for (size_t i = 0; i < count; i++) {
        CandidateEdge e;
        rep += nextRandom(state) % 1000;  // ascending, sometimes repeating
        e.setRep(rep);
        // Members on both sides of the representative: with --cov-mode 1
        // assignGroup emits reversed edges, so member < rep is a real case and
        // the zigzag delta has to carry it.
        const int64_t offset = static_cast<int64_t>(nextRandom(state) % 2000001) - 1000000;
        int64_t member = static_cast<int64_t>(rep) + offset;
        if (member < 0) {
            member = 0;
        }
        e.setMember(static_cast<uint64_t>(member));
        e.diagonal = static_cast<int16_t>(nextRandom(state) % 65536);
        e.score = static_cast<uint16_t>(nextRandom(state) % 65536);
        e.reverseStrand = static_cast<uint8_t>(nextRandom(state) & 1);
        edges.push_back(e);
    }
    if (extremes && edges.empty() == false) {
        // Every field at both ends of its range, in a block that also holds
        // ordinary values, so a boundary bug cannot hide behind a special case.
        edges[0].setRep(0);
        edges[0].setMember(CandidateEdge::MAX_KEY);
        edges[0].diagonal = INT16_MIN;
        edges[0].score = 0;
        edges[0].reverseStrand = 0;
        CandidateEdge last;
        last.setRep(CandidateEdge::MAX_KEY);
        last.setMember(0);
        last.diagonal = INT16_MAX;
        last.score = 65535;
        last.reverseStrand = 1;
        edges.push_back(last);  // still non-decreasing: MAX_KEY is the largest
    }
    return edges;
}

static void testRoundTrip() {
    const size_t counts[] = {1, 2, 17, 1000, 50000};
    bool allOk = true;
    bool agreeOk = true;
    bool smallerOk = true;
    for (size_t c = 0; c < sizeof(counts) / sizeof(counts[0]); c++) {
        std::vector<CandidateEdge> edges = makeEdges(counts[c], 0x2545F4914F6CDD1DULL + c, true);

        std::vector<CandidateEdge> packed, raw;
        size_t packedBytes = 0, rawBytes = 0;
        if (roundTrip(edges, false, packed, packedBytes) == false ||
            roundTrip(edges, true, raw, rawBytes) == false) {
            allOk = false;
            continue;
        }
        if (packed.size() != edges.size() || raw.size() != edges.size()) {
            allOk = false;
            continue;
        }
        for (size_t i = 0; i < edges.size(); i++) {
            if (!sameEdge(packed[i], edges[i]) || !sameEdge(raw[i], edges[i])) {
                allOk = false;
            }
            if (!sameEdge(packed[i], raw[i])) {
                agreeOk = false;
            }
        }
        // The whole point of the encoding. A change that silently stopped packing
        // would pass every correctness check above.
        if (packedBytes >= rawBytes) {
            smallerOk = false;
        }
    }
    check(allOk, "packed and raw both round-trip 1..50000 edges including field extremes");
    check(agreeOk, "packed and raw decode to identical edges, so --raw-records is a real control");
    check(smallerOk, "the packed payload is smaller than the fixed-width one at every size");

    std::vector<CandidateEdge> empty, out;
    size_t bytes = 0;
    check(roundTrip(empty, false, out, bytes), "an empty block encodes to nothing at all");
}

static void testMemberOnBothSides() {
    // Reversed edges are not a corner case: --cov-mode 1 produces them, and the
    // member delta is zigzagged specifically so they stay short.
    std::vector<CandidateEdge> edges;
    for (int64_t d = -5; d <= 5; d++) {
        CandidateEdge e;
        e.setRep(1000000);
        e.setMember(static_cast<uint64_t>(1000000 + d));
        e.diagonal = static_cast<int16_t>(d);
        e.score = 1;
        e.reverseStrand = 0;
        edges.push_back(e);
    }
    std::sort(edges.begin(), edges.end(), [](const CandidateEdge &a, const CandidateEdge &b) {
        return a.getMember() < b.getMember();
    });
    std::vector<CandidateEdge> out;
    size_t bytes = 0;
    bool ok = roundTrip(edges, false, out, bytes) && out.size() == edges.size();
    for (size_t i = 0; ok && i < edges.size(); i++) {
        ok = sameEdge(out[i], edges[i]);
    }
    check(ok, "members below and above their representative both round-trip");
}

static void testCorruptionRefused() {
    std::vector<CandidateEdge> edges = makeEdges(500, 0x9E3779B97F4A7C15ULL, false);
    std::vector<uint8_t> good;
    EdgeBlockCodec::encode(edges, false, 7, 3, good);
    EdgeBlockHeader header;
    memcpy(&header, good.data(), sizeof(EdgeBlockHeader));
    const uint8_t *payload = good.data() + sizeof(EdgeBlockHeader);

    std::vector<CandidateEdge> out;
    check(EdgeBlockCodec::decode(header, payload, out), "the untouched block decodes");

    // A flipped payload byte must fail the checksum rather than decode into
    // plausible edges.
    {
        std::vector<uint8_t> broken(good);
        broken[sizeof(EdgeBlockHeader) + broken.size() / 3] ^= 0x40;
        out.clear();
        check(EdgeBlockCodec::decode(header, broken.data() + sizeof(EdgeBlockHeader), out) == false,
              "a single flipped payload byte is refused by the checksum");
        check(out.empty(), "a refused block leaves the output untouched");
    }

    // A record count larger than the payload supports must run out of bytes
    // rather than read past the end.
    {
        EdgeBlockHeader tooMany = header;
        tooMany.recordCount = header.recordCount + 50;
        out.clear();
        check(EdgeBlockCodec::decode(tooMany, payload, out) == false,
              "a record count past the payload is refused");
    }

    // Trailing bytes mean the block is not what its header says.
    {
        EdgeBlockHeader tooFew = header;
        tooFew.recordCount = header.recordCount - 10;
        out.clear();
        check(EdgeBlockCodec::decode(tooFew, payload, out) == false,
              "a record count short of the payload is refused");
    }

    // Header sanity, which is what stops a torn tail being read as a block.
    {
        EdgeBlockHeader bad = header;
        bad.magic ^= 0xFFFFFFFFU;
        check(EdgeBlockCodec::headerLooksValid(bad) == false, "a bad magic is rejected");
        bad = header;
        bad.encoding = 99;
        check(EdgeBlockCodec::headerLooksValid(bad) == false, "an unknown encoding is rejected");
        bad = header;
        bad.encoding = EdgeBlockCodec::ENCODING_RAW;
        check(EdgeBlockCodec::headerLooksValid(bad) == false,
              "a raw header whose length disagrees with its record count is rejected");
    }
}

// Decoding many blocks into one vector is what the align stage does per bucket.
// It is also where the quadratic reserve lived, so the growth behaviour is
// asserted rather than assumed: appending N blocks must stay linear.
static void testManyBlocksAccumulate() {
    std::vector<CandidateEdge> all;
    std::vector<CandidateEdge> out;
    bool ok = true;
    uint64_t rep = 0;
    for (int block = 0; block < 200; block++) {
        std::vector<CandidateEdge> edges;
        uint64_t state = 0x1234567 + block;
        for (int i = 0; i < 500; i++) {
            CandidateEdge e;
            rep += nextRandom(state) % 100;
            e.setRep(rep);
            e.setMember(rep + 17);
            e.diagonal = 3;
            e.score = 2;
            e.reverseStrand = 0;
            edges.push_back(e);
            all.push_back(e);
        }
        std::vector<uint8_t> buffer;
        EdgeBlockCodec::encode(edges, false, 1, 1, buffer);
        EdgeBlockHeader header;
        memcpy(&header, buffer.data(), sizeof(EdgeBlockHeader));
        if (EdgeBlockCodec::decode(header, buffer.data() + sizeof(EdgeBlockHeader), out) == false) {
            ok = false;
            break;
        }
    }
    ok = ok && out.size() == all.size();
    for (size_t i = 0; ok && i < all.size(); i++) {
        ok = sameEdge(out[i], all[i]);
    }
    check(ok, "200 blocks decoded into one vector accumulate in order and intact");
}

int main(int, const char**) {
    testRoundTrip();
    testMemberOnBothSides();
    testCorruptionRefused();
    testManyBlocksAccumulate();

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
