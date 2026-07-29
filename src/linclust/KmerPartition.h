#ifndef MMSEQS_KMERPARTITION_H
#define MMSEQS_KMERPARTITION_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// K-mer space partitioning for the distributed linclust map/reduce.
//
// Stock kmermatcher handles a database too large for memory by *splitting*: it
// re-reads the whole database once per split and re-extracts every k-mer, then
// throws away the k-mers whose hash falls outside the split's range
// (kmermatcher.cpp:322, loop at :139-144). At 1e12 sequences on a 2 TB node that
// is ~68 splits, so the extraction work -- the single most expensive part of the
// pipeline -- is done ~68 times over.
//
// Here the database is scanned *once*. Every extracted k-mer is appended to the
// bucket file of its partition, so the extraction is paid for exactly once and
// the reduce stage later reads back one self-contained partition at a time.
//
// Why partitioning by the existing 16-bit hash is correct: kmermatcher already
// computes `score = (unsigned short) hashUInt64(kmerIdx, hashShift)` for every
// k-mer (kmermatcher.cpp:202,213,223). That is a *pure function of the k-mer
// index*, so two occurrences of the same k-mer always produce the same score and
// therefore always land in the same partition. Grouping only ever compares equal
// k-mers, so a partition can be grouped in isolation and the result is identical
// to grouping the whole database at once. This is what makes the partitioning
// lossless rather than an approximation -- unlike sequence-space sharding, which
// co-locates near-duplicates only about a third of the time at 90% identity.
class KmerPartitioner {
public:
    // partitionCount must be a power of two and at most 65536, since the score
    // it partitions is 16 bits wide.
    explicit KmerPartitioner(unsigned int partitionCount);

    unsigned int getPartitionCount() const { return partitionCount; }

    // Partition of a k-mer, from the 16-bit hash kmermatcher already computed.
    // Taking the high bits keeps whole contiguous runs of the score space
    // together, which is the same space the stock split ranges carve up, so a
    // partition is directly comparable to a stock split.
    unsigned int partitionOf(unsigned short score) const {
        return static_cast<unsigned int>(score) >> shift;
    }

private:
    unsigned int partitionCount;
    unsigned int shift;
};

// One k-mer occurrence as stored in a bucket file.
//
// 24 bytes, packed. Two fields are here specifically to delete arrays that stock
// kmermatcher sizes by *key space* rather than by entry count, and which
// therefore cannot exist at 1e12:
//   - seqLen removes `seqkey_to_len[dbKeySize]` (kmermatcher.cpp:1236),
//   - id is carried explicitly so nothing needs a dense key-indexed side table.
// `countTable` (:1256) and `repSequence` (:1357) fall the same way.
struct __attribute__((__packed__)) KmerRecord {
    uint64_t kmer;
    // 48 bits is 2.8e14 keys, comfortably past the 1e12 target, and saves 2 bytes
    // per record against a full 64-bit key -- 40 TB at 1e12 with 21 k-mers each.
    uint8_t idBytes[6];
    uint16_t pos;
    uint16_t seqLen;
    // Residues flanking the k-mer, as kmermatcher's --include-adjacency 1 keeps.
    uint8_t adjacent[6];

    uint64_t getId() const {
        uint64_t id = 0;
        for (int i = 5; i >= 0; i--) {
            id = (id << 8) | idBytes[i];
        }
        return id;
    }

    void setId(uint64_t id) {
        for (int i = 0; i < 6; i++) {
            idBytes[i] = static_cast<uint8_t>(id & 0xFF);
            id >>= 8;
        }
    }

    static const uint64_t MAX_ID = (static_cast<uint64_t>(1) << 48) - 1;
};

// Converts a bucket record into the KmerPosition that assignGroup consumes, and
// back. Templated rather than including kmermatcher.h so this header stays free
// of the whole DBReader/Parameters chain; the caller instantiates it with the
// KmerPosition variant it uses.
//
// The variant must be one with IncludeSeqLen = true. That is the point of the
// record carrying seqLen: with IncludeSeqLen = false, KmerPosition::getSeqLen()
// reads the *static* `SeqLenData<T,false>::seqkey_to_len` array, which is sized
// by key space (kmermatcher.cpp:1236) and so cannot exist at 1e12. Reading the
// length out of the record is what deletes it.
template <typename KmerPositionT>
void kmerRecordToPosition(const KmerRecord &record, KmerPositionT &out) {
    out.kmer = record.kmer;
    out.id = static_cast<decltype(out.id)>(record.getId());
    out.pos = static_cast<decltype(out.pos)>(record.pos);
    out.sl.setSeqLen(static_cast<decltype(out.pos)>(record.seqLen));
    for (int i = 0; i < 6; i++) {
        out.setAdjacentSeq(i, record.adjacent[i]);
    }
}

template <typename KmerPositionT>
void kmerPositionToRecord(KmerPositionT &in, KmerRecord &out) {
    out.kmer = in.kmer;
    out.setId(static_cast<uint64_t>(in.id));
    out.pos = static_cast<uint16_t>(in.pos);
    out.seqLen = static_cast<uint16_t>(in.getSeqLen());
    for (int i = 0; i < 6; i++) {
        out.adjacent[i] = in.getAdjacentSeq(i);
    }
}

// Appends k-mer records into per-partition bucket files.
//
// Each writer owns one shard id and writes <dir>/p<partition>/<shard>.kmers, so
// concurrent writers -- threads within a worker, and workers across nodes --
// never touch the same file and need no locking at all. Bucket files are laid
// out one directory per partition because the reduce stage reads exactly one
// partition, and because a flat directory would hold partitions x shards entries.
class KmerBucketWriter {
public:
    // bufferBudgetBytes is split evenly across partitions, so memory is bounded
    // regardless of partition count. Records accumulate per partition and are
    // flushed in large contiguous appends rather than one write per k-mer.
    KmerBucketWriter(const std::string &dir, unsigned int partitionCount,
                     const std::string &shardId, size_t bufferBudgetBytes = 32 * 1024 * 1024);
    ~KmerBucketWriter();

    void append(unsigned int partition, const KmerRecord &record);
    // Flushes and closes every open bucket. Called by the destructor, but call it
    // explicitly to see write errors before the object goes away.
    void close();

    uint64_t getRecordCount() const { return recordCount; }

    // Creates the per-partition directories. Safe to call from many workers.
    static void createLayout(const std::string &dir, unsigned int partitionCount);
    static std::string partitionDir(const std::string &dir, unsigned int partition);

private:
    KmerBucketWriter(const KmerBucketWriter &);
    KmerBucketWriter &operator=(const KmerBucketWriter &);

    void flush(unsigned int partition);

    std::string dir;
    std::string shardId;
    unsigned int partitionCount;
    size_t recordsPerBuffer;
    std::vector<std::vector<KmerRecord> > buffers;
    std::vector<FILE *> files;
    uint64_t recordCount;
};

// Reads every shard of one partition back.
class KmerBucketReader {
public:
    // Total records across all shards of the partition, so the reduce stage can
    // size its array in one allocation before reading.
    static uint64_t countRecords(const std::string &dir, unsigned int partition);

    // Appends every record of the partition to out.
    static void readPartition(const std::string &dir, unsigned int partition,
                              std::vector<KmerRecord> &out);

    static std::vector<std::string> shardFiles(const std::string &dir, unsigned int partition);
};

#endif
