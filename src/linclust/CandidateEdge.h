#ifndef MMSEQS_CANDIDATEEDGE_H
#define MMSEQS_CANDIDATEEDGE_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

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

    static const uint32_t MAGIC = 0x45444745;  // "EDGE"
};

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
// One file per (worker, bucket), like the k-mer shards, so no locking is needed.
class EdgeBucketWriter {
public:
    EdgeBucketWriter(const std::string &dir, unsigned int bucketCount, const std::string &shardId,
                     size_t bufferBudgetBytes = 256 * 1024 * 1024);
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
