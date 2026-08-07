// Tests for k-mer space partitioning, the basis of the distributed kmermatcher.
//
// The property everything rests on is that partitioning k-mer space is
// *lossless*: every occurrence of a given k-mer lands in exactly one partition,
// so a partition can be grouped in isolation and the union of the partitions is
// exactly the input. If that ever fails, candidate pairs go missing silently and
// the clustering is quietly wrong rather than visibly broken -- so it is checked
// here directly, on a corpus with heavy k-mer repetition.

#include "CandidateEdge.h"
#include "FileUtil.h"
#include "KmerPartition.h"
#include "kmermatcher.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <unistd.h>

const char* binary_name = "test_kmerpartition";

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
    char tmpl[] = "/tmp/mmseqs_kmerpart_testXXXXXX";
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

// Stands in for kmermatcher's `(unsigned short) hashUInt64(kmerIdx, hashShift)`.
// The only property the partitioner relies on is that it is a pure function of
// the k-mer, which this reproduces.
static unsigned short scoreOf(uint64_t kmer) {
    uint64_t h = kmer * 0x9E3779B97F4A7C15ULL;
    h ^= h >> 29;
    h *= 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 32;
    return static_cast<unsigned short>(h & 0xFFFF);
}

static void testRecordLayout() {
    check(sizeof(KmerRecord) == 24, "k-mer record is 24 packed bytes");

    KmerRecord record;
    record.kmer = 0xDEADBEEFCAFEF00DULL;
    record.pos = 4242;
    record.seqLen = 65535;
    for (int i = 0; i < 6; i++) {
        record.adjacent[i] = static_cast<uint8_t>(i + 1);
    }

    // Boundary ids: the 48-bit field must survive its extremes intact, since a
    // silent truncation here would corrupt sequence identity, not just size.
    const uint64_t ids[] = {0, 1, 255, 256, 1000000000000ULL, KmerRecord::MAX_ID};
    bool allRoundTrip = true;
    for (size_t i = 0; i < sizeof(ids) / sizeof(ids[0]); i++) {
        record.setId(ids[i]);
        allRoundTrip = allRoundTrip && record.getId() == ids[i];
    }
    check(allRoundTrip, "48-bit id round-trips at 0, 1, 255, 256, 1e12 and 2^48-1");
    check(KmerRecord::MAX_ID > 1000000000000ULL, "id field spans past the 1e12 target");

    record.setId(1000000000000ULL);
    check(record.kmer == 0xDEADBEEFCAFEF00DULL && record.pos == 4242 && record.seqLen == 65535 &&
              record.adjacent[5] == 6,
          "setting the id leaves the neighbouring packed fields untouched");
}

static void testPartitioner() {
    KmerPartitioner partitioner(8192);
    check(partitioner.getPartitionCount() == 8192, "partitioner reports its partition count");

    // The plan sizes P = 8192 over the 16-bit hash space, i.e. 8 hash values per
    // partition. Verify that is exactly what happens, and that the whole space
    // maps inside range. The 8 values a partition owns are strided rather than
    // contiguous, which is the point -- see the prefix test below.
    std::set<unsigned int> seen;
    bool inRange = true;
    std::map<unsigned int, int> perPartition;
    for (unsigned int score = 0; score <= 65535; score++) {
        const unsigned int p = partitioner.partitionOf(static_cast<unsigned short>(score));
        inRange = inRange && p < 8192;
        seen.insert(p);
        perPartition[p]++;
    }
    check(inRange, "every 16-bit score maps inside the partition range");
    check(seen.size() == 8192, "every partition is reachable");
    bool evenlySized = true;
    for (std::map<unsigned int, int>::const_iterator it = perPartition.begin();
         it != perPartition.end(); ++it) {
        evenlySized = evenlySized && it->second == 8;
    }
    check(evenlySized, "each of 8192 partitions covers exactly 8 hash values");

    KmerPartitioner single(1);
    bool allZero = true;
    for (unsigned int score = 0; score <= 65535; score += 97) {
        allZero = allZero && single.partitionOf(static_cast<unsigned short>(score)) == 0;
    }
    check(allZero, "a single partition collects the whole hash space");

    // The property the balance depends on, and the one an earlier high-bits
    // partitioner failed. Linclust keeps the bottom-scoring k-mers of a sequence,
    // so the scores that reach a bucket are a *prefix* of the space, not a sample
    // of all of it. A partitioner that spreads the whole space evenly can still
    // send an entire prefix to one partition; this checks the prefix itself
    // spreads. Measured on real data the high-bits version gave partition 0 68% of
    // 19.9M k-mers at P = 16.
    KmerPartitioner small(16);
    std::map<unsigned int, int> prefixCounts;
    const unsigned int prefixEnd = 4096;  // a plausible per-sequence threshold
    for (unsigned int score = 0; score < prefixEnd; score++) {
        prefixCounts[small.partitionOf(static_cast<unsigned short>(score))]++;
    }
    check(prefixCounts.size() == 16, "a prefix of the score space reaches every partition");
    bool prefixBalanced = true;
    for (std::map<unsigned int, int>::const_iterator it = prefixCounts.begin();
         it != prefixCounts.end(); ++it) {
        prefixBalanced = prefixBalanced && it->second == static_cast<int>(prefixEnd) / 16;
    }
    check(prefixBalanced, "a prefix of the score space spreads evenly over the partitions");
}

// The load-bearing test: run a repetitive corpus through the writer, read it
// back partition by partition, and prove nothing was lost, duplicated or split.
//
// Run twice, once per encoding. The packed form drops seqLen from the record and
// recovers it from the length-rank table, so running the same corpus both ways
// tests the table and the codec together against a single expectation -- which is
// exactly what --raw-records exists to make possible in a real run.
static void testLosslessRoundTrip(const std::string &dir, bool raw) {
    const unsigned int partitionCount = 64;
    const std::string suffix = raw ? "_raw" : "_packed";
    const std::string bucketDir = dir + "/buckets" + suffix;
    KmerBucketWriter::createLayout(bucketDir, partitionCount);
    KmerPartitioner partitioner(partitionCount);

    const uint64_t sequences = 20000;

    // Lengths must be *non-increasing in the key*, because that is what a
    // length-ranked database guarantees and what the table encodes. Fifty
    // sequences per length, so there are runs to binary-search rather than one
    // length per key.
    std::vector<LengthRankTable::Run> runs;
    for (uint64_t seqId = 0; seqId < sequences; seqId += 50) {
        LengthRankTable::Run run;
        run.length = 500 - (seqId / 50);
        run.firstKey = seqId;
        run.count = 50;
        runs.push_back(run);
    }
    const std::string tableDb = bucketDir + "/db";
    LengthRankTable::write(tableDb, runs, sequences);
    LengthRankTable lengths;
    lengths.open(tableDb);
    const LengthRankTable *lengthsForRead = raw ? NULL : &lengths;

    // Few distinct k-mers over many sequences, so most k-mers recur -- which is
    // the case where a partitioning bug would actually lose candidate pairs.
    const uint64_t distinctKmers = 500;
    std::map<uint64_t, std::vector<uint64_t> > expected;  // kmer -> sequence ids
    srand(20260728);

    // Three shards, standing in for three threads or three worker nodes.
    const char *shardNames[] = {"w0", "w1", "w2"};
    std::vector<KmerBucketWriter *> writers;
    for (int s = 0; s < 3; s++) {
        writers.push_back(
            new KmerBucketWriter(bucketDir, partitionCount, shardNames[s], 4096, 0, 0, raw));
    }

    uint64_t written = 0;
    for (uint64_t seqId = 0; seqId < sequences; seqId++) {
        const int shard = static_cast<int>(seqId % 3);
        const uint16_t seqLen = static_cast<uint16_t>(500 - (seqId / 50));
        for (int k = 0; k < 21; k++) {
            const uint64_t kmer = static_cast<uint64_t>(rand()) % distinctKmers;
            KmerRecord record;
            record.kmer = kmer;
            record.setId(seqId);
            record.pos = static_cast<uint16_t>(k);
            record.seqLen = seqLen;
            for (int i = 0; i < 6; i++) {
                // Residue indices inside a 21-letter alphabet, plus the UCHAR_MAX
                // "no adjacency" sentinel on every seventh sequence: the sentinel
                // is the one value the 5-bit packing has to special-case, and
                // assignGroup tests it directly.
                record.adjacent[i] = (seqId % 7 == 0)
                                         ? static_cast<uint8_t>(UCHAR_MAX)
                                         : static_cast<uint8_t>((kmer + i) % 21);
            }
            writers[shard]->append(partitioner.partitionOf(scoreOf(kmer)), record);
            expected[kmer].push_back(seqId);
            written++;
        }
    }
    for (int s = 0; s < 3; s++) {
        writers[s]->close();
        delete writers[s];
    }

    uint64_t counted = 0;
    for (unsigned int p = 0; p < partitionCount; p++) {
        counted += KmerBucketReader::countRecords(bucketDir, p);
    }
    check(counted == written,
          "every written record is counted back across partitions" + suffix);

    // Read each partition and check the k-mers it holds belong to it alone.
    std::map<uint64_t, unsigned int> kmerPartition;
    std::map<uint64_t, std::vector<uint64_t> > recovered;
    bool fieldsIntact = true;
    bool partitionPure = true;
    uint64_t readBack = 0;
    for (unsigned int p = 0; p < partitionCount; p++) {
        std::vector<KmerRecord> records;
        KmerBucketReader::readPartition(bucketDir, p, records, lengthsForRead);
        for (size_t i = 0; i < records.size(); i++) {
            const KmerRecord &record = records[i];
            // Every occurrence of this k-mer must be in this partition and no other.
            if (kmerPartition.count(record.kmer) == 0) {
                kmerPartition[record.kmer] = p;
            } else if (kmerPartition[record.kmer] != p) {
                partitionPure = false;
            }
            if (partitioner.partitionOf(scoreOf(record.kmer)) != p) {
                partitionPure = false;
            }
            const uint64_t seqId = record.getId();
            const uint8_t wantAdjacent =
                (seqId % 7 == 0) ? static_cast<uint8_t>(UCHAR_MAX)
                                 : static_cast<uint8_t>((record.kmer + 0) % 21);
            fieldsIntact = fieldsIntact && record.seqLen == 500 - (seqId / 50) &&
                           record.adjacent[0] == wantAdjacent && record.pos < 21;
            recovered[record.kmer].push_back(seqId);
            readBack++;
        }
    }

    check(readBack == written, "every record survives the write/read round trip" + suffix);
    check(partitionPure, "each k-mer appears in exactly one partition" + suffix);
    check(fieldsIntact, "id, seqLen, position and adjacency survive the round trip" + suffix);

    bool sameMultiset = recovered.size() == expected.size();
    for (std::map<uint64_t, std::vector<uint64_t> >::iterator it = recovered.begin();
         sameMultiset && it != recovered.end(); ++it) {
        std::vector<uint64_t> &got = it->second;
        std::vector<uint64_t> want = expected[it->first];
        std::sort(got.begin(), got.end());
        std::sort(want.begin(), want.end());
        sameMultiset = got == want;
    }
    check(sameMultiset,
          "grouping partition by partition sees exactly the same k-mer/sequence pairs as the "
          "whole input" + suffix);

    // A legitimately empty partition is one createLayout made and no worker wrote
    // to -- an *existing*, readable directory with no shards in it. A directory
    // that cannot be opened is a failure, not an empty partition, because the map
    // always creates the layout before the reduce runs; treating the two the same
    // let a worker record an item done having silently read nothing.
    const std::string emptyDir = bucketDir + "_empty";
    KmerBucketWriter::createLayout(emptyDir, 1);
    std::vector<KmerRecord> empty;
    KmerBucketReader::readPartition(emptyDir, 0, empty, lengthsForRead);
    check(empty.empty() && KmerBucketReader::countRecords(emptyDir, 0) == 0,
          "a partition nothing was written to reads back empty rather than failing" + suffix);
}

// The packed encoding must be smaller than the fixed-width one it replaces --
// that is the entire reason it exists, and a change that quietly stopped packing
// would otherwise pass every correctness check above.
static void testPackedIsSmaller(const std::string &dir) {
    uint64_t bytes[2] = {0, 0};
    for (int raw = 0; raw < 2; raw++) {
        const std::string bucketDir = dir + "/buckets" + (raw ? "_raw" : "_packed");
        for (unsigned int p = 0; p < 64; p++) {
            const std::vector<std::string> shards = KmerBucketReader::shardFiles(bucketDir, p);
            for (size_t i = 0; i < shards.size(); i++) {
                bytes[raw] += FileUtil::getFileSize(shards[i]);
            }
        }
    }
    const double perRecordRaw = static_cast<double>(bytes[1]) / (20000.0 * 21.0);
    const double perRecordPacked = static_cast<double>(bytes[0]) / (20000.0 * 21.0);
    fprintf(stdout, "      raw %.2f B/record, packed %.2f B/record (%.2fx)\n", perRecordRaw,
            perRecordPacked, perRecordRaw / perRecordPacked);
    // The fixed-width record is 24 B plus framing; the packed one should land well
    // under it even on this synthetic corpus, whose tiny k-mer values flatter the
    // k-mer field but whose 4096-byte buffers make the per-block header costly.
    check(bytes[0] * 10 < bytes[1] * 8, "the packed encoding is at least 20% smaller than raw");
}

// The reduce stage feeds records into assignGroup through KmerPosition, so the
// record has to carry everything that struct exposes. In particular seqLen must
// come back out of getSeqLen() without any global key-indexed table, which is the
// whole reason the field is in the record.
static void testKmerPositionConversion() {
    typedef KmerPosition<short, true, true> Position;

    // The widest id that survives *both* the record's 48-bit field and DBKeyType.
    // A hardcoded 48-bit constant silently truncates in the default build, where
    // MMSEQS_INT64_IDS is 0 and DBKeyType is uint32_t, so this check has to be
    // written against the key width the build actually has.
    const uint64_t wideId =
        std::min<uint64_t>(KmerRecord::MAX_ID, static_cast<uint64_t>(DB_KEY_INVALID) - 1);

    KmerRecord record;
    record.kmer = 0x0123456789ABCDEFULL;
    record.setId(wideId);
    record.pos = 517;
    record.seqLen = 30000;
    for (int i = 0; i < 6; i++) {
        record.adjacent[i] = static_cast<uint8_t>(i * 3 + 1);
    }

    Position position;
    kmerRecordToPosition(record, position);

    bool adjacencyKept = true;
    for (int i = 0; i < 6; i++) {
        adjacencyKept = adjacencyKept && position.getAdjacentSeq(i) == static_cast<unsigned char>(i * 3 + 1);
    }
    check(position.kmer == record.kmer && static_cast<uint64_t>(position.id) == wideId &&
              position.pos == 517,
          "record converts into the KmerPosition assignGroup consumes");
    check(position.getSeqLen() == 30000,
          "sequence length comes from the record, needing no key-indexed length table");
    check(adjacencyKept, "adjacency survives conversion into KmerPosition");

    KmerRecord back;
    kmerPositionToRecord(position, back);
    bool sameBytes = back.kmer == record.kmer && back.getId() == record.getId() &&
                     back.pos == record.pos && back.seqLen == record.seqLen;
    for (int i = 0; i < 6; i++) {
        sameBytes = sameBytes && back.adjacent[i] == record.adjacent[i];
    }
    check(sameBytes, "KmerPosition converts back into an identical record");
}

// Pins the derivation to the two scales the project actually targets, because
// the whole argument for computing P rather than exposing it is that the right
// answer differs between them.
static void testShuffleSizing() {
    const uint64_t TB = 1024ULL * 1024 * 1024 * 1024;
    const uint64_t GB = 1024ULL * 1024 * 1024;

    // 100B: 1e11 sequences must fit a hard 100 TB budget. Persistent load is the
    // 25 TB database plus ~7.6 TB of surviving edges.
    KmerShuffleSizing small = deriveKmerShuffleSizing(100000000000ULL, 21, 100 * TB,
                                                      33 * TB, 64 * GB);
    check(small.totalKmerBytes == 100000000000ULL * 21 * 24,
          "100B k-mer volume follows from 21 records of 24 bytes per sequence");
    check(small.waveCount == 1, "100B fits the 100 TB budget in a single wave");
    check(small.bytesPerWave <= 100 * TB - 33 * TB,
          "100B peak k-mer bytes stay inside the budget after the database and edges");
    check(small.partitionCount == 4096, "100B derives P = 4096");
    // P bounds the reduce's RESIDENT set, which is ~1.9x the bucket's on-disk
    // bytes before candidate edges; the factor of 3 is what keeps a partition
    // inside a worker rather than 1.9x outside it.
    check(small.bytesPerPartition * 3 <= 64 * GB, "100B partitions fit worker memory when expanded");
    check(small.bytesPerPartition <= 64 * GB, "100B buckets fit the per-worker memory");

    // 1T: 1e12 sequences, no hard ceiling, so only per-worker memory sets P.
    KmerShuffleSizing large = deriveKmerShuffleSizing(1000000000000ULL, 21, 0, 0, 64 * GB);
    check(large.waveCount == 1, "an unlimited budget means a single wave");
    check(large.partitionCount == 32768, "1T derives P = 32768");
    check(large.bytesPerPartition * 3 <= 64 * GB, "1T partitions fit worker memory when expanded");
    check(large.bytesPerPartition <= 64 * GB, "1T buckets fit the per-worker memory");
    check(large.partitionCount > small.partitionCount,
          "the two target scales genuinely want different P, which is why it is derived");

    // A budget that cannot hold the whole shuffle at once must split into waves,
    // and more waves must mean smaller buckets rather than a larger P.
    KmerShuffleSizing waved = deriveKmerShuffleSizing(1000000000000ULL, 21, 400 * TB,
                                                      100 * TB, 64 * GB);
    check(waved.waveCount > 1, "a budget below the k-mer volume forces multiple waves");
    check(waved.bytesPerWave <= 300 * TB, "each wave stays inside the budget");
    check(waved.bytesPerPartition <= 64 * GB, "waved buckets still fit the per-worker memory");
    check(waved.partitionCount <= large.partitionCount,
          "waves reduce the peak, so no more partitions are needed than without them");

    // A wave owns a contiguous slice of partitions, so the slices are only equal
    // -- and peak scratch only really 1/W of the whole -- if W divides P. Both
    // being powers of two is how that is guaranteed.
    bool wavesDivide = true;
    // From 8 TB up: below that, 504 TB of k-mers needs more than the 64 waves
    // deriveKmerShuffleSizing now refuses as a misconfigured budget.
    for (uint64_t budget = 16; budget <= 512; budget *= 2) {
        // Sweep budgets from far below the k-mer volume to above it. Per-worker
        // memory is large so nothing but the wave count can raise P, which is the
        // case where the two used to disagree.
        const KmerShuffleSizing s =
            deriveKmerShuffleSizing(1000000000000ULL, 21, budget * TB, 0, 1024 * GB);
        wavesDivide = wavesDivide && s.waveCount > 0 &&
                      (s.waveCount & (s.waveCount - 1)) == 0 &&
                      s.partitionCount >= s.waveCount &&
                      s.partitionCount % s.waveCount == 0 &&
                      // the largest slice really is within budget
                      (s.totalKmerBytes / s.partitionCount) *
                              (s.partitionCount / s.waveCount) <= budget * TB;
    }
    check(wavesDivide, "the wave count is a power of two that divides P at every budget");

    // P must always be a usable power of two.
    bool powerOfTwo = true;
    for (uint64_t seqs = 1000000; seqs <= 100000000000ULL; seqs *= 10) {
        const KmerShuffleSizing s = deriveKmerShuffleSizing(seqs, 21, 0, 0, 8 * GB);
        powerOfTwo = powerOfTwo && s.partitionCount > 0 &&
                     (s.partitionCount & (s.partitionCount - 1)) == 0 &&
                     s.partitionCount <= 65536;
    }
    check(powerOfTwo, "derived P is always a power of two inside the 16-bit hash space");

    KmerShuffleSizing tiny = deriveKmerShuffleSizing(1000, 21, 0, 0, 64 * GB);
    check(tiny.partitionCount == 1 && tiny.waveCount == 1,
          "a small input collapses to one partition and one wave");
}

// P is the reduce's unit of work as well as its unit of memory, so a per-node
// memory limit large enough to hold everything used to derive P = 1 and leave the
// whole allocation but one worker with nothing to claim. --workers raises it.
void testShuffleSizingWorkerFloor() {
    const uint64_t TB = 1024ULL * 1024 * 1024 * 1024;
    const uint64_t GB = 1024ULL * 1024 * 1024;
    // 1e9 sequences against a 800 GB node: the memory sizing alone gives a P far
    // below the worker count, which is exactly the measured failure.
    const KmerShuffleSizing untold = deriveKmerShuffleSizing(1000000000ULL, 21, 0, 0, 800 * GB);
    const KmerShuffleSizing told = deriveKmerShuffleSizing(1000000000ULL, 21, 0, 0, 800 * GB, 8);
    check(told.partitionCount >= 8,
          "--workers gives every reduce worker at least one partition to claim");
    check(told.partitionCount >= untold.partitionCount,
          "--workers never lowers P below what the memory limit required");

    // The hint must not break the invariants the rest of the shuffle relies on.
    bool powerOfTwo = true;
    bool wavesDivide = true;
    for (unsigned int workers = 1; workers <= 1024; workers *= 2) {
        const KmerShuffleSizing s =
            deriveKmerShuffleSizing(1000000000ULL, 21, 100 * TB, 0, 64 * GB, workers);
        powerOfTwo = powerOfTwo && s.partitionCount > 0 &&
                     (s.partitionCount & (s.partitionCount - 1)) == 0 &&
                     s.partitionCount <= 65536;
        wavesDivide = wavesDivide && s.partitionCount % s.waveCount == 0;
    }
    check(powerOfTwo, "P stays a power of two inside the hash space at every worker count");
    check(wavesDivide, "the wave count still divides P at every worker count");

    // A large --workers against a small database must not shatter it into
    // partitions holding almost nothing.
    const KmerShuffleSizing small =
        deriveKmerShuffleSizing(10000000ULL, 21, 0, 0, 64 * GB, 4096);
    check(small.totalKmerBytes / small.partitionCount >= 64ULL * 1024 * 1024,
          "--workers never derives partitions smaller than the useful floor");
    // And an absurd worker count is clamped rather than killing the run.
    const KmerShuffleSizing absurd =
        deriveKmerShuffleSizing(1000000000ULL, 21, 0, 0, 64 * GB, 4000000000U);
    check(absurd.partitionCount <= 65536,
          "an absurd --workers is clamped into the hash space rather than fatal");

    // Not told is the previous behaviour, exactly.
    const KmerShuffleSizing zero = deriveKmerShuffleSizing(1000000000ULL, 21, 0, 0, 800 * GB, 0);
    check(zero.partitionCount == untold.partitionCount,
          "--workers 0 leaves the memory-only sizing untouched");
}

// The align stage's parallelism is capped by the edge bucket count, so the same
// failure mode applies: a generous memory limit must not collapse it to one.
void testAlignBucketCount() {
    const uint64_t GiB = 1024ULL * 1024 * 1024;
    const uint64_t MiB = 1024ULL * 1024;

    // The measured 1e9 case: 125 GB of edges on a 800 GB node produced one bucket
    // and left alignparallel running on a single node for half the run.
    const unsigned int big = deriveAlignBucketCount(125 * GiB, 1000000000ULL, 800 * GiB, 0);
    check(big >= 64, "a large edge set is split into many buckets despite a large memory limit");

    // More memory must never mean less parallelism. Swept from 256 MiB, because
    // min(mem/4, 1 GiB) pins the target at 1 GiB for anything at or above 4 GiB --
    // a sweep that starts above that never leaves the flat region and would pass
    // whatever the code did.
    bool monotone = true;
    bool sawMemoryMatter = false;
    unsigned int previous = deriveAlignBucketCount(125 * GiB, 1000000000ULL, 256 * MiB, 0);
    for (uint64_t limit = 512ULL * 1024 * 1024; limit <= 1024 * GiB; limit *= 2) {
        const unsigned int now = deriveAlignBucketCount(125 * GiB, 1000000000ULL, limit, 0);
        monotone = monotone && now <= previous;
        sawMemoryMatter = sawMemoryMatter || now != previous;
        previous = now;
    }
    check(monotone, "raising the memory limit never raises the bucket count");
    check(sawMemoryMatter, "the sweep covers the range where the memory limit changes the count");
    check(deriveAlignBucketCount(125 * GiB, 1000000000ULL, 1024 * GiB, 0) >= 64,
          "even an unlimited-looking memory limit keeps the byte target's buckets");

    // The property that actually protects correctness: every representative key
    // lands in a bucket that exists. bucketSpan is what the reduce computes from
    // the returned count, so the two have to agree for all of them.
    bool everyKeyMaps = true;
    const uint64_t sizes[] = {1, 2, 3, 1000, 1048576, 100000000ULL, 1000000000ULL};
    const uint64_t limits[] = {0, 1, 1024, 256ULL * MiB, 8 * GiB, 800 * GiB};
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
        for (size_t j = 0; j < sizeof(limits) / sizeof(limits[0]); j++) {
            for (unsigned int w = 0; w <= 64; w = w ? w * 2 : 1) {
                const uint64_t entries = sizes[i];
                const unsigned int b =
                    deriveAlignBucketCount(entries * 189, entries, limits[j], w);
                const uint64_t span = (entries + b - 1) / b;
                everyKeyMaps = everyKeyMaps && b >= 1 && b <= MAX_ALIGN_BUCKETS && span >= 1 &&
                               (entries - 1) / span < b;
            }
        }
    }
    check(everyKeyMaps,
          "every representative key maps into an existing bucket at every count");

    // The worker hint only ever raises it.
    const unsigned int small = deriveAlignBucketCount(2 * GiB, 10000000ULL, 800 * GiB, 0);
    const unsigned int hinted = deriveAlignBucketCount(2 * GiB, 10000000ULL, 800 * GiB, 8);
    check(hinted >= 8, "--workers gives every align worker a bucket on a small input");
    check(hinted >= small, "--workers never lowers the bucket count");

    // Bounds.
    check(deriveAlignBucketCount(1024, 4, 800 * GiB, 64) <= 4,
          "there are never more buckets than representative keys");
    check(deriveAlignBucketCount(0, 0, 800 * GiB, 0) == 1,
          "an empty edge set is a single bucket");
    check(deriveAlignBucketCount(1ULL << 50, 1ULL << 40, 0, 1000000) <= MAX_ALIGN_BUCKETS,
          "the bucket count stays inside the 16-bit id space");
}

int main(int, char **) {
    const std::string dir = makeTempDir();

    testRecordLayout();
    testPartitioner();
    testShuffleSizing();
    testShuffleSizingWorkerFloor();
    testAlignBucketCount();
    testKmerPositionConversion();
    testLosslessRoundTrip(dir, true);
    testLosslessRoundTrip(dir, false);
    testPackedIsSmaller(dir);

    removeTempDir(dir);

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
