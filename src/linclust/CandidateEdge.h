#ifndef MMSEQS_CANDIDATEEDGE_H
#define MMSEQS_CANDIDATEEDGE_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "VarintCodec.h"

// One (representative, member) pair proposed by the distributed reduce.
//
// Packed binary rather than a prefilter DB. A `pref` DB stores each hit as ASCII
// with a per-representative index entry, and at 1e11 sequences the index alone is
// per-key state no single node can hold. 17 bytes per edge is also ~4x smaller
// than the text form, which matters directly: this is the largest intermediate
// the pipeline writes.
struct __attribute__((__packed__)) CandidateEdge {
    // 48-bit keys, as in KmerRecord: 2.8e14 keys, past the 1e12 target, and two
    // bytes cheaper per side than a full 64-bit key.
    uint8_t repBytes[6];
    uint8_t memberBytes[6];
    // 16 bits deliberately, and *not* wide enough to hold every diagonal a 65535
    // residue sequence can produce. That is stock's convention, not an oversight:
    //
    //   - stock stores it in `short KmerEntry::diagonal` (kmermatcher.h:187,199)
    //     and truncates an int into it at kmermatcher.cpp:2050, exactly as
    //     collectRoundEdges does here;
    //   - the prefilter hands it on as `unsigned short hit_t::diagonal`
    //     (QueryMatcher.h:36);
    //   - the aligner *undoes* the truncation. ungappedAlign takes an
    //     `unsigned short` and DistanceCalculator::computeUngappedAlignment
    //     (DistanceCalculator.h:93-112) tries every real diagonal congruent to it
    //     mod 65536 that the two lengths allow, keeping the best-scoring one.
    //
    // So a diagonal of 45000 stored as -20536 is recovered by the aligner, and two
    // real diagonals that alias mod 65536 alias in stock too. Widening this would
    // diverge from stock rather than converge on it, and cost 6% of the largest
    // intermediate the pipeline writes.
    int16_t diagonal;
    // Nucleotide strand. Stock carries it in bit 63 of the representative key,
    // which a 48-bit key has no room for.
    uint8_t reverseStrand;
    // How many k-mers put this pair on this diagonal. Stock's prefilter score,
    // and the value the align stage ranks diagonals by.
    //
    // 16 bits, not 8. At 8 bits this saturated on 199 of 4.9M pass-2 records
    // (pass 2 uses --kmer-per-seq-scale aa:0.100, so a long sequence contributes
    // thousands), which turned a diagonal comparison into a tie and moved 31 of
    // 1,000,000 sequences; widening it took that residual to 8. The bound is
    // kmersPerSeq(65535) x rounds ~= 26k, which fits.
    //
    // Note this is not bit-for-bit stock behaviour, in either width. Stock's
    // KmerEntry::score is an unsigned char assigned from an int
    // (kmermatcher.h:188, kmermatcher.cpp:2048), so a single entry above 255
    // *wraps* rather than saturating, and the wrapped value is then summed into a
    // wider accumulator. Saturating is the more defensible reading of "how many
    // k-mers agree on this diagonal" and measured closer to stock, but the two
    // can still disagree on a pair whose per-entry count exceeds 255.
    uint16_t score;

    // 48-bit keys, matching KmerRecord::MAX_ID. Checked on decode: a packed edge
    // whose delta arithmetic produced a key past this is corrupt, not large.
    static const uint64_t MAX_KEY = (static_cast<uint64_t>(1) << 48) - 1;

    uint64_t getRep() const { return get(repBytes); }
    uint64_t getMember() const { return get(memberBytes); }
    void setRep(uint64_t key) { set(repBytes, key); }
    void setMember(uint64_t key) { set(memberBytes, key); }

private:
    static uint64_t get(const uint8_t *bytes) {
        uint64_t value = 0;
        for (int i = 5; i >= 0; i--) {
            value = (value << 8) | bytes[i];
        }
        return value;
    }
    static void set(uint8_t *bytes, uint64_t value) {
        for (int i = 0; i < 6; i++) {
            bytes[i] = static_cast<uint8_t>(value & 0xFF);
            value >>= 8;
        }
    }
};

// Buffered writer for one bucket's edge file, used for the *alignment* output:
// alignparallel writes one of these per bucket and greedycluster reads them back.
//
// One file per bucket, written under a per-process temporary name and renamed on
// close. The rename is atomic, so a bucket redone after a crash -- or run twice
// because a lease lapsed -- replaces the earlier file with an equally complete
// one rather than appending to it. That is why this stage needs no block framing,
// unlike EdgeBucketWriter below, whose output several producers share.
class EdgeWriter {
public:
    // workerId names the temporary file this writes before renaming. It has to be
    // the globally unique worker id rather than a pid, because the path is on the
    // shared filesystem and pids collide across nodes.
    EdgeWriter(const std::string &path, int64_t workerId, size_t bufferRecords = 1024 * 1024);
    ~EdgeWriter();

    void append(const CandidateEdge &edge);
    void close();

    uint64_t getEdgeCount() const { return edgeCount; }

    static std::string partitionPath(const std::string &dir, unsigned int partition);

private:
    EdgeWriter(const EdgeWriter &);
    EdgeWriter &operator=(const EdgeWriter &);

    void flush();

    std::string path;
    int64_t workerId;
    std::string tmpPath;
    FILE *file;
    std::vector<CandidateEdge> buffer;
    size_t bufferRecords;
    uint64_t edgeCount;
    bool closed;
};


// Framing for one flushed run of edges, so a record can be traced to the k-mer
// partition and the worker that produced it.
//
// This is what makes the reduce idempotent across a crash. A worker that flushes
// a partition's edges and dies before the queue records the item leaves that data
// on disk; the item is then re-claimed and redone by another worker, whose edges
// land in *its* shard. Both copies are present, and the align stage sums the
// support of matching (pair, diagonal) records, so the redone pair would count
// twice and can win a diagonal it should not.
//
// The duplicates cannot be recognised from the records alone: unlike a k-mer
// record, where (kmer, id, pos) names one occurrence, two different partitions
// may legitimately emit byte-identical edges that must both be summed. Naming the
// producer in the frame is what separates "the same partition written twice" from
// "two partitions that agree".
struct __attribute__((__packed__)) EdgeBlockHeader {
    // Distinguishes a real header from the tail of a torn write.
    uint32_t magic;
    uint32_t partition;
    uint32_t worker;
    uint32_t recordCount;
    // ENCODING_RAW or ENCODING_PACKED. Named per block, so --raw-records is a
    // writer-side switch that costs the reader nothing and a directory holding
    // both kinds still reads correctly.
    uint8_t encoding;
    uint8_t reserved[3];
    uint64_t payloadBytes;
    // FNV-1a over the payload. Length alone catches a truncated block; this also
    // catches a block boundary that has drifted, which with variable-width
    // records would otherwise decode into plausible nonsense.
    uint32_t checksum;
    uint32_t reserved2;

    static const uint32_t MAGIC = 0x45444745;  // "EDGE"
};

// Packs a sorted run of edges.
//
// This is the intermediate that decides whether the pipeline fits its scratch
// budget. Waves divide the k-mer shuffle but not the edges -- alignment runs once,
// after every wave -- so at 1e12 the edges are the term nothing else reduces:
// ~22 records per sequence measured at 1e9, at 17 fixed-width bytes each.
//
// Four fields, each encoded against what the surrounding design already
// guarantees:
//
//   - `rep` as a delta. The reduce sorts its edges by (rep, member, diagonal,
//     strand) and appends them in that order, and bucket == rep / bucketSpan, so
//     a bucket's appends are a subsequence of a sorted sequence and every block is
//     sorted by construction. No sort is needed here.
//   - `member` as a zigzag delta *from its own rep*, not from the previous
//     member. Keys are length-ranked and homologous sequences have similar
//     lengths, so |rep - member| is small for most edges -- 82.5% below 1e7 on
//     measured UniRef100 data. This is a direct dividend of the length-ranked key
//     design; with input-order keys it would be a full 6 bytes.
//   - `diagonal` zigzagged. It stays an int16_t: stock truncates an int into a
//     short (kmermatcher.cpp:2050) and DistanceCalculator undoes it by trying
//     every diagonal congruent mod 65536, so widening it would diverge from stock
//     rather than converge on it.
//   - `score` and `reverseStrand` share a varint, strand in the low bit. Strand
//     belongs in the sort key (two opposite-strand edges on one diagonal are
//     different alignments), but it is one bit and a byte of its own would be
//     6% of the record.
namespace EdgeBlockCodec {

const uint8_t ENCODING_RAW = 0;
const uint8_t ENCODING_PACKED = 1;

// Ceiling on the bytes one edge can encode to, used to size the writer's budget
// against both the edge buffer and the encode buffer it packs into.
const size_t MAX_ENCODED_BYTES_PER_EDGE = VarintCodec::MAX_BYTES * 2 + 3 + 3;

// Appends one framed block. Requires edges sorted by non-decreasing rep, which
// the reduce guarantees; the encoder checks rather than assumes, because a
// negative delta would wrap into a huge one and decode into a different edge.
void encode(const std::vector<CandidateEdge> &edges, bool raw, uint32_t partition, uint32_t worker,
            std::vector<uint8_t> &out);

// Decodes one block payload, appending to out. Returns false on corruption rather
// than exiting, so a caller can stop at the last good block.
bool decode(const EdgeBlockHeader &header, const uint8_t *payload,
            std::vector<CandidateEdge> &out);

bool headerLooksValid(const EdgeBlockHeader &header);

uint32_t checksum(const uint8_t *data, size_t size);

}  // namespace EdgeBlockCodec

// Writes edges into buckets by *representative key range*.
//
// This is the layout the alignment stage needs, and the reason is not obvious:
// aligning inside the k-mer partition fails because a partition's pairs are
// scattered over the whole key space, so every partition ends up re-reading the
// entire sequence database (measured: 52x amplification, see DESIGN_DECISIONS.md
// §9). Bucketing by representative key gives each align worker a contiguous slice
// of sequences instead.
//
// Two things fall out for free:
//   - every copy of a pair produced by different k-mer partitions lands in the
//     same bucket, so the cross-partition duplicates are removed here rather than
//     paid for in duplicate alignments;
//   - having all copies together means the per-(pair, diagonal) score
//     accumulation stock does in its global merge can be reproduced exactly,
//     instead of approximated per partition.
//
// The bucket id is a 16-bit value on disk, so this is a hard ceiling rather than
// a policy choice.
const unsigned int MAX_ALIGN_BUCKETS = 65536;

// Chooses how many edge buckets the reduce writes.
//
// A bucket is two things at once, and sizing it for only the first is what made
// this pipeline fail to scale:
//
//   1. a slice that has to fit an align worker, which loads its sequences whole;
//   2. the unit of work the align stage's queue hands out, which makes the count
//      a hard ceiling on how many workers that stage can use at all.
//
// Deriving it from (1) alone means a *generous* --split-memory-limit produces
// *fewer* buckets and less parallelism -- backwards, and silent, because the
// starved workers find nothing claimable and exit 0. Measured at 1e9 sequences
// with 800G per node: 125 GB of edges against a 200 GB target left bucketCount at
// 1, and alignparallel ran 3h16m on one node while three others slept, half the
// wall clock of the whole run.
//
// So memoryLimitBytes becomes a ceiling and the target is taken well below it,
// and workerCount (0 when not told) raises the count further so an allocation
// always has something to claim. The result never exceeds entryCount, since a
// bucket is a range of representative key.
unsigned int deriveAlignBucketCount(uint64_t dataSize, uint64_t entryCount,
                                    uint64_t memoryLimitBytes, unsigned int workerCount);

// One file per (worker, bucket), like the k-mer shards, so no locking is needed.
class EdgeBucketWriter {
public:
    EdgeBucketWriter(const std::string &dir, unsigned int bucketCount, const std::string &shardId,
                     size_t bufferBudgetBytes = 256 * 1024 * 1024, bool rawRecords = false);
    ~EdgeBucketWriter();

    // Names the partition every subsequent append belongs to, and the worker
    // doing the work. Flushes anything still buffered for the previous partition
    // first, so a block never spans two of them.
    void beginPartition(unsigned int partition, int64_t worker);

    void append(unsigned int bucket, const CandidateEdge &edge);
    // Pushes buffered edges to the OS. Call before marking a work item done, for
    // the same reason the k-mer writer does.
    void flushAll();
    void close();

    uint64_t getEdgeCount() const { return edgeCount; }

    static void createLayout(const std::string &dir, unsigned int bucketCount);
    static std::string bucketDir(const std::string &dir, unsigned int bucket);
    static std::vector<std::string> shardFiles(const std::string &dir, unsigned int bucket);

private:
    EdgeBucketWriter(const EdgeBucketWriter &);
    EdgeBucketWriter &operator=(const EdgeBucketWriter &);

    void flush(unsigned int bucket);
    void closeFile(unsigned int bucket);
    std::string shardPath(unsigned int bucket) const;

    std::string dir;
    std::string shardId;
    unsigned int bucketCount;
    size_t edgesPerBuffer;
    std::vector<std::vector<CandidateEdge> > buffers;
    // One reusable encode target per bucket, so a flush allocates nothing.
    std::vector<std::vector<uint8_t> > encodeBuffers;
    std::vector<FILE *> files;
    // Descriptors kept open across flushes, up to descriptorBudget. Plain counters:
    // one writer belongs to one process and its appends come from the single
    // thread that drains the work queue.
    size_t descriptorBudget;
    size_t openFiles;
    uint64_t edgeCount;
    bool closed;
    unsigned int currentPartition;
    int64_t currentWorker;
    // Writes fixed-width edges instead of packed blocks; see --raw-records.
    bool rawRecords;
};

// Reads the blocks EdgeBucketWriter wrote, keeping only those whose producer the
// reduce's work queue recorded as the one that completed the partition.
//
// `authority[partition]` is that worker (see WorkQueue::readCompletedWorkers). A
// block naming any other worker is a dead worker's copy of an item that was
// redone, and is skipped. An empty authority disables filtering, for callers that
// have no queue to consult.
class EdgeBucketReader {
public:
    static size_t readShard(const std::string &path, const std::vector<int64_t> &authority,
                            std::vector<CandidateEdge> &out);
};

#endif
