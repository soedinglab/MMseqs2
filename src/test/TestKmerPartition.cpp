// Tests for k-mer space partitioning, the basis of the distributed kmermatcher.
//
// The property everything rests on is that partitioning k-mer space is
// *lossless*: every occurrence of a given k-mer lands in exactly one partition,
// so a partition can be grouped in isolation and the union of the partitions is
// exactly the input. If that ever fails, candidate pairs go missing silently and
// the clustering is quietly wrong rather than visibly broken -- so it is checked
// here directly, on a corpus with heavy k-mer repetition.

#include "KmerPartition.h"
#include "kmermatcher.h"

#include <algorithm>
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
    // maps inside range.
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
}

// The load-bearing test: run a repetitive corpus through the writer, read it
// back partition by partition, and prove nothing was lost, duplicated or split.
static void testLosslessRoundTrip(const std::string &dir) {
    const unsigned int partitionCount = 64;
    const std::string bucketDir = dir + "/buckets";
    KmerBucketWriter::createLayout(bucketDir, partitionCount);
    KmerPartitioner partitioner(partitionCount);

    // Few distinct k-mers over many sequences, so most k-mers recur -- which is
    // the case where a partitioning bug would actually lose candidate pairs.
    const uint64_t distinctKmers = 500;
    std::map<uint64_t, std::vector<uint64_t> > expected;  // kmer -> sequence ids
    srand(20260728);

    // Three shards, standing in for three threads or three worker nodes.
    const char *shardNames[] = {"w0", "w1", "w2"};
    std::vector<KmerBucketWriter *> writers;
    for (int s = 0; s < 3; s++) {
        writers.push_back(new KmerBucketWriter(bucketDir, partitionCount, shardNames[s], 4096));
    }

    uint64_t written = 0;
    for (uint64_t seqId = 0; seqId < 20000; seqId++) {
        const int shard = static_cast<int>(seqId % 3);
        for (int k = 0; k < 21; k++) {
            const uint64_t kmer = static_cast<uint64_t>(rand()) % distinctKmers;
            KmerRecord record;
            record.kmer = kmer;
            record.setId(seqId);
            record.pos = static_cast<uint16_t>(k);
            record.seqLen = static_cast<uint16_t>(100 + (seqId % 400));
            for (int i = 0; i < 6; i++) {
                record.adjacent[i] = static_cast<uint8_t>((kmer + i) & 0xFF);
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
    check(counted == written, "every written record is counted back across partitions");

    // Read each partition and check the k-mers it holds belong to it alone.
    std::map<uint64_t, unsigned int> kmerPartition;
    std::map<uint64_t, std::vector<uint64_t> > recovered;
    bool fieldsIntact = true;
    bool partitionPure = true;
    uint64_t readBack = 0;
    for (unsigned int p = 0; p < partitionCount; p++) {
        std::vector<KmerRecord> records;
        KmerBucketReader::readPartition(bucketDir, p, records);
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
            fieldsIntact = fieldsIntact && record.seqLen == 100 + (seqId % 400) &&
                           record.adjacent[0] == static_cast<uint8_t>(record.kmer & 0xFF);
            recovered[record.kmer].push_back(seqId);
            readBack++;
        }
    }

    check(readBack == written, "every record survives the write/read round trip");
    check(partitionPure, "each k-mer appears in exactly one partition");
    check(fieldsIntact, "id, seqLen and adjacency survive the round trip");

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
          "grouping partition by partition sees exactly the same k-mer/sequence pairs as the whole input");

    std::vector<KmerRecord> empty;
    KmerBucketReader::readPartition(bucketDir + "_missing", 0, empty);
    check(empty.empty() && KmerBucketReader::countRecords(bucketDir + "_missing", 0) == 0,
          "a partition nothing was written to reads back empty rather than failing");
}

// The reduce stage feeds records into assignGroup through KmerPosition, so the
// record has to carry everything that struct exposes. In particular seqLen must
// come back out of getSeqLen() without any global key-indexed table, which is the
// whole reason the field is in the record.
static void testKmerPositionConversion() {
    typedef KmerPosition<short, true, true> Position;

    KmerRecord record;
    record.kmer = 0x0123456789ABCDEFULL;
    record.setId(999999999999ULL);
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
    check(position.kmer == record.kmer && static_cast<uint64_t>(position.id) == 999999999999ULL &&
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

int main(int, char **) {
    const std::string dir = makeTempDir();

    testRecordLayout();
    testPartitioner();
    testKmerPositionConversion();
    testLosslessRoundTrip(dir);

    removeTempDir(dir);

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
