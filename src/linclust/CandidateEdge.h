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
    int16_t diagonal;
    // Nucleotide strand. Stock carries it in bit 63 of the representative key,
    // which a 48-bit key has no room for.
    uint8_t reverseStrand;
    // How many k-mers put this pair on this diagonal. Stock's prefilter score,
    // and the value the align stage ranks diagonals by.
    //
    // 16 bits, not 8. The align stage picks a pair's diagonal by comparing these
    // counts, so a ceiling low enough to be reached turns a comparison into a tie
    // and hands the decision to a tie-break stock does not have. Pass 1 extracts
    // 21 k-mers per sequence and never came close; pass 2 uses
    // --kmer-per-seq-scale aa:0.100, so a long sequence contributes thousands, and
    // at 255 this saturated on 199 of 4.9M records, which gave 23 pairs a
    // different diagonal than stock and moved 31 of 1,000,000 sequences into a
    // different cluster; widening it here took that residual to 8. The bound is kmersPerSeq(65535) x rounds ~= 26k, which
    // fits 16 bits; stock accumulates in an int and has no ceiling at all.
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

// Buffered writer for one partition's edge file.
//
// One file per partition, written whole, so a partition redone after a crash
// simply overwrites its file. That makes the reduce idempotent, unlike the map,
// whose per-worker shards can retain a dead worker's partial output.
//
// Used for the reference layout (`--align`, small scale). Production uses
// EdgeBucketWriter below.
class EdgeWriter {
public:
    EdgeWriter(const std::string &path, size_t bufferRecords = 1024 * 1024);
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
    std::string tmpPath;
    FILE *file;
    std::vector<CandidateEdge> buffer;
    size_t bufferRecords;
    uint64_t edgeCount;
    bool closed;
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

    std::string dir;
    std::string shardId;
    unsigned int bucketCount;
    size_t edgesPerBuffer;
    std::vector<std::vector<CandidateEdge> > buffers;
    std::vector<FILE *> files;
    uint64_t edgeCount;
    bool closed;
};

#endif
