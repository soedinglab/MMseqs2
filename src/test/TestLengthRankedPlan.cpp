// Tests for the placement math behind the distributed, length-ranked createdb.
//
// Pass 2 writes sequences straight to their final byte offsets with no merge and
// no communication between workers, so the plan is the only thing keeping those
// writes from colliding. The properties under test are therefore the ones that
// make the output a valid database at all:
//
//   - keys tile [0, N) exactly: no key unassigned, none assigned twice;
//   - keys are length-ranked: reading in key order yields non-increasing lengths;
//   - data and header byte ranges tile their files exactly, so no worker can
//     overwrite another's bytes and no gap is left uninitialised;
//   - the plan depends only on the input, not on the order histograms arrive in.

#include "LengthRankedPlan.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <unistd.h>

const char* binary_name = "test_lengthrankedplan";

static int failures = 0;

static void check(bool condition, const std::string &what) {
    if (condition == false) {
        fprintf(stderr, "FAIL: %s\n", what.c_str());
        failures++;
    } else {
        fprintf(stdout, "ok:   %s\n", what.c_str());
    }
}

static std::string makeTempDir() {
    char tmpl[] = "/tmp/mmseqs_lrplan_testXXXXXX";
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) {
        perror("mkdtemp");
        exit(EXIT_FAILURE);
    }
    return std::string(dir);
}

static void removeTempDir(const std::string &dir) {
    std::string cmd = "rm -rf '" + dir + "'";
    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "could not remove %s\n", dir.c_str());
    }
}

static ChunkHistogram makeHistogram(uint64_t chunkIdx, uint64_t fileIdx,
                                    const std::vector<ChunkHistogram::Bucket> &buckets) {
    ChunkHistogram histogram;
    histogram.chunkIdx = chunkIdx;
    histogram.fileIdx = fileIdx;
    histogram.buckets = buckets;
    for (size_t i = 0; i < buckets.size(); i++) {
        histogram.seqCount += buckets[i].count;
    }
    return histogram;
}

static ChunkHistogram::Bucket bucket(uint64_t length, uint64_t count, uint64_t headerBytes) {
    ChunkHistogram::Bucket b;
    b.length = length;
    b.count = count;
    b.headerBytes = headerBytes;
    return b;
}

// A worked two-chunk example, checked entry by entry against hand-computed
// offsets. The randomised test below proves the invariants hold in general; this
// one pins down that the arithmetic itself is the intended arithmetic.
static void testWorkedExample() {
    std::vector<ChunkHistogram> histograms;
    histograms.push_back(makeHistogram(0, 0, {bucket(5, 2, 20), bucket(10, 1, 12)}));
    histograms.push_back(makeHistogram(1, 0, {bucket(5, 3, 33), bucket(7, 1, 9)}));

    std::vector<ChunkPlan> plans;
    LengthRankedTotals totals = buildLengthRankedPlan(histograms, plans);

    check(totals.seqCount == 7, "worked example totals 7 sequences");
    // 1*(10+2) + 1*(7+2) + 5*(5+2) = 12 + 9 + 35
    check(totals.dataBytes == 56, "worked example totals 56 data bytes");
    check(totals.headerBytes == 74, "worked example totals 74 header bytes");
    check(totals.maxSeqLen == 10, "worked example reports the longest sequence");

    check(plans.size() == 2, "worked example produces one plan per chunk");
    // Entries are stored ascending by length.
    check(plans[0].entries.size() == 2 && plans[1].entries.size() == 2,
          "worked example plans keep every bucket");

    const ChunkPlan::Entry &c0len10 = plans[0].entries[1];
    check(c0len10.length == 10 && c0len10.keyStart == 0 && c0len10.dataOffset == 0 && c0len10.hdrOffset == 0,
          "longest sequence takes key 0 at offset 0");

    const ChunkPlan::Entry &c1len7 = plans[1].entries[1];
    check(c1len7.length == 7 && c1len7.keyStart == 1 && c1len7.dataOffset == 12 && c1len7.hdrOffset == 12,
          "second-longest sequence follows the longest");

    // Both chunks hold length 5; chunk 0 must come first because ties break by
    // chunk index, which is input order.
    const ChunkPlan::Entry &c0len5 = plans[0].entries[0];
    const ChunkPlan::Entry &c1len5 = plans[1].entries[0];
    check(c0len5.keyStart == 2 && c0len5.dataOffset == 21 && c0len5.hdrOffset == 21,
          "lower chunk index wins a length tie");
    check(c1len5.keyStart == 4 && c1len5.dataOffset == 35 && c1len5.hdrOffset == 41,
          "higher chunk index follows within the same length");
}

// Materialises every (key, length, dataOffset, headerBytes) run the plan implies
// and checks the runs tile all three address spaces exactly.
static void testTilingInvariants() {
    srand(20260728);

    for (int trial = 0; trial < 20; trial++) {
        const size_t chunkCount = 1 + (rand() % 12);
        std::vector<ChunkHistogram> histograms;
        for (size_t c = 0; c < chunkCount; c++) {
            // Lengths are drawn from a small shared pool so chunks collide on
            // lengths often, which is what exercises the tie-breaking.
            std::vector<ChunkHistogram::Bucket> buckets;
            for (uint64_t length = 1; length <= 20; length++) {
                if (rand() % 3 == 0) {
                    continue;
                }
                const uint64_t count = 1 + (rand() % 5);
                buckets.push_back(bucket(length, count, count * (3 + (rand() % 7))));
            }
            histograms.push_back(makeHistogram(c, c % 3, buckets));
        }

        std::vector<ChunkHistogram> reference = histograms;
        std::vector<ChunkPlan> plans;
        LengthRankedTotals totals = buildLengthRankedPlan(histograms, plans);

        // Flatten the plan into runs ordered by key.
        struct Run {
            uint64_t keyStart, count, length, dataOffset, hdrOffset, hdrBytes;
        };
        std::vector<Run> runs;
        for (size_t i = 0; i < plans.size(); i++) {
            for (size_t e = 0; e < plans[i].entries.size(); e++) {
                const ChunkPlan::Entry &entry = plans[i].entries[e];
                uint64_t hdrBytes = 0;
                for (size_t b = 0; b < histograms[i].buckets.size(); b++) {
                    if (histograms[i].buckets[b].length == entry.length) {
                        hdrBytes = histograms[i].buckets[b].headerBytes;
                    }
                }
                Run run = {entry.keyStart, entry.count, entry.length,
                           entry.dataOffset, entry.hdrOffset, hdrBytes};
                runs.push_back(run);
            }
        }
        std::sort(runs.begin(), runs.end(),
                  [](const Run &a, const Run &b) { return a.keyStart < b.keyStart; });

        uint64_t expectedKey = 0, expectedData = 0, expectedHdr = 0;
        uint64_t previousLength = UINT64_MAX;
        bool keysTile = true, dataTiles = true, hdrTiles = true, lengthRanked = true;
        for (size_t r = 0; r < runs.size(); r++) {
            keysTile = keysTile && runs[r].keyStart == expectedKey;
            dataTiles = dataTiles && runs[r].dataOffset == expectedData;
            hdrTiles = hdrTiles && runs[r].hdrOffset == expectedHdr;
            lengthRanked = lengthRanked && runs[r].length <= previousLength;

            previousLength = runs[r].length;
            expectedKey += runs[r].count;
            expectedData += runs[r].count * (runs[r].length + 2);
            expectedHdr += runs[r].hdrBytes;
        }

        if (trial == 0) {
            check(keysTile, "key ranges tile [0, N) with no gap or overlap");
            check(dataTiles, "data byte ranges tile the data file exactly");
            check(hdrTiles, "header byte ranges tile the header file exactly");
            check(lengthRanked, "sequences are ordered longest first in key order");
            check(expectedKey == totals.seqCount, "totals agree with the tiled key count");
            check(expectedData == totals.dataBytes, "totals agree with the tiled data size");
            check(expectedHdr == totals.headerBytes, "totals agree with the tiled header size");
        } else if (!keysTile || !dataTiles || !hdrTiles || !lengthRanked ||
                   expectedKey != totals.seqCount || expectedData != totals.dataBytes ||
                   expectedHdr != totals.headerBytes) {
            check(false, "tiling invariants hold on randomised trial " + std::to_string(trial));
            return;
        }

        // Feeding the histograms in a different order must not change the plan:
        // workers finish in arbitrary order, so the planner must not depend on it.
        std::vector<ChunkHistogram> shuffled = reference;
        std::reverse(shuffled.begin(), shuffled.end());
        std::vector<ChunkPlan> shuffledPlans;
        buildLengthRankedPlan(shuffled, shuffledPlans);
        bool identical = shuffledPlans.size() == plans.size();
        for (size_t i = 0; identical && i < plans.size(); i++) {
            identical = shuffledPlans[i].chunkIdx == plans[i].chunkIdx &&
                        shuffledPlans[i].entries.size() == plans[i].entries.size();
            for (size_t e = 0; identical && e < plans[i].entries.size(); e++) {
                identical = shuffledPlans[i].entries[e].keyStart == plans[i].entries[e].keyStart &&
                            shuffledPlans[i].entries[e].dataOffset == plans[i].entries[e].dataOffset &&
                            shuffledPlans[i].entries[e].hdrOffset == plans[i].entries[e].hdrOffset;
            }
        }
        if (trial == 0) {
            check(identical, "plan is independent of the order histograms are collected in");
        } else if (!identical) {
            check(false, "plan stayed order-independent on randomised trial " + std::to_string(trial));
            return;
        }
    }
    check(true, "tiling invariants hold across 20 randomised chunk layouts");
}

static void testRoundTrip(const std::string &dir) {
    ChunkHistogram histogram = makeHistogram(7, 3, {bucket(4, 9, 40), bucket(11, 2, 18)});
    histogram.nuclVotes = 5;
    histogram.sampleCount = 11;
    const std::string histPath = dir + "/chunk7.hist";
    histogram.write(histPath);
    ChunkHistogram loaded = ChunkHistogram::read(histPath);

    bool same = loaded.chunkIdx == 7 && loaded.fileIdx == 3 && loaded.seqCount == 11 &&
                loaded.nuclVotes == 5 && loaded.sampleCount == 11 &&
                loaded.buckets.size() == 2;
    for (size_t i = 0; same && i < loaded.buckets.size(); i++) {
        same = loaded.buckets[i].length == histogram.buckets[i].length &&
               loaded.buckets[i].count == histogram.buckets[i].count &&
               loaded.buckets[i].headerBytes == histogram.buckets[i].headerBytes;
    }
    check(same, "chunk histogram survives a write/read round trip");

    std::vector<ChunkHistogram> histograms(1, histogram);
    std::vector<ChunkPlan> plans;
    buildLengthRankedPlan(histograms, plans);
    const std::string planPath = dir + "/chunk7.plan";
    plans[0].write(planPath);
    ChunkPlan loadedPlan = ChunkPlan::read(planPath);

    bool planSame = loadedPlan.chunkIdx == plans[0].chunkIdx &&
                    loadedPlan.fileIdx == plans[0].fileIdx &&
                    loadedPlan.entries.size() == plans[0].entries.size();
    for (size_t i = 0; planSame && i < loadedPlan.entries.size(); i++) {
        planSame = loadedPlan.entries[i].length == plans[0].entries[i].length &&
                   loadedPlan.entries[i].count == plans[0].entries[i].count &&
                   loadedPlan.entries[i].keyStart == plans[0].entries[i].keyStart &&
                   loadedPlan.entries[i].dataOffset == plans[0].entries[i].dataOffset &&
                   loadedPlan.entries[i].hdrOffset == plans[0].entries[i].hdrOffset;
    }
    check(planSame, "chunk plan survives a write/read round trip");
}

static void testEmptyInput() {
    std::vector<ChunkHistogram> histograms;
    histograms.push_back(makeHistogram(0, 0, {}));
    std::vector<ChunkPlan> plans;
    LengthRankedTotals totals = buildLengthRankedPlan(histograms, plans);
    check(totals.seqCount == 0 && totals.dataBytes == 0 && totals.headerBytes == 0,
          "an empty chunk plans to an empty database");
}

int main(int, char **) {
    const std::string dir = makeTempDir();

    testWorkedExample();
    testTilingInvariants();
    testRoundTrip(dir);
    testEmptyInput();

    removeTempDir(dir);

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
