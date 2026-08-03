#ifndef MMSEQS_KMERPARTITION_H
#define MMSEQS_KMERPARTITION_H

#include <cstdint>
#include <cstdio>
#include <mutex>
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
    //
    // The *low* bits, not the high ones. Linclust keeps the bottom-scoring k-mers
    // of each sequence (kmermatcher.cpp:318, `score < threshold`), so the selected
    // scores are not spread over the 16-bit space at all -- they pile up against
    // zero. Partitioning on the high bits therefore puts almost everything in
    // partition 0: measured on 1M sequences at P=16 it gave partition 0 68% of all
    // k-mers and partition 15 0.3%, which would make the reduce's per-partition
    // memory bound meaningless. The selected scores form a prefix [0, threshold)
    // of the space, and the low bits of a prefix are uniform, so masking balances
    // the partitions without needing the hash-distribution pass stock's
    // setupKmerSplits runs to carve its uneven split ranges.
    //
    // Any pure function of the score would be *correct*; only balance picks
    // between them, since equal k-mers always share a score.
    unsigned int partitionOf(unsigned short score) const {
        return static_cast<unsigned int>(score) & mask;
    }

private:
    unsigned int partitionCount;
    unsigned int mask;
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

// The two numbers that shape the k-mer shuffle, and where they came from.
struct KmerShuffleSizing {
    uint64_t totalKmerBytes;     // every k-mer record, summed over all waves
    // Both powers of two, and waveCount divides partitionCount: a wave owns the
    // contiguous slice [w * P / W, (w + 1) * P / W) of partition space, and the
    // slices must be equal for peak scratch to actually be 1 / W of the whole.
    unsigned int waveCount;      // extraction passes, so peak scratch stays in budget
    unsigned int partitionCount; // P
    uint64_t bytesPerWave;       // peak k-mer bytes on disk at any one time
    uint64_t bytesPerPartition;  // what one worker loads in the reduce
};

// Derives the wave count and P from the scratch budget and per-worker memory,
// rather than making either a knob.
//
// Both are over-determined by things already known: the k-mer volume follows
// from the sequence count at a measured 21 k-mers x 24 B per sequence, the wave
// count is whatever keeps peak scratch inside the budget, and P is the smallest
// power of two whose buckets still fit a worker. Exposing them raw invites two
// failures that only show up deep into a run -- P too small and workers die of
// memory, P too large and the reduce silently pays P/W sequential re-scans of
// the whole database for no benefit.
//
// The two scales this is sized for pull in opposite directions, which is exactly
// why it should be computed: 50.4 TB of k-mers at 100B against 504 TB at 1T. Note
// that P follows from workerMemoryBytes as much as from the budget, so a target
// scale alone does not pin it.
//
// persistentBytes is everything sharing the scratch budget with the k-mer wave:
// what is already on disk (the caller passes a measured figure, so pass 2 can
// account for what pass 1 left behind) plus the candidate edges this stage's own
// reduce will write.
//
// scratchBudgetBytes == 0 means unlimited, giving a single wave.
KmerShuffleSizing deriveKmerShuffleSizing(uint64_t sequenceCount,
                                          unsigned int kmersPerSequence,
                                          uint64_t scratchBudgetBytes,
                                          uint64_t persistentBytes,
                                          uint64_t workerMemoryBytes);

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
// One writer per worker *process*, shared by all its threads, writing
// <dir>/p<partition>/<shard>.kmers. Two constraints pick that granularity:
//
//   - File count. A writer per thread would leave workers x threads x partitions
//     shards -- 32 million at 500 workers, 64 threads and P = 1024 -- which no
//     shared filesystem enjoys. Per process it is workers x partitions, so 0.5 M
//     at the 100B target and 4 M at 1T, a few hundred files per directory.
//   - Write size. The buffer budget is split across partitions, so per process a
//     partition gets budget/P; per thread it would get budget/(threads x P),
//     which at these partition counts is a few hundred bytes -- the write size
//     that makes a parallel filesystem collapse.
//
// Threads therefore share the per-partition buffers under a per-partition mutex.
// There are P independent locks and the k-mer hash spreads threads uniformly
// across them, so the expected contention is threads/P, and the k-mer extraction
// that produces each record costs far more than the uncontended acquire.
//
// Workers still never share a file, so nothing is locked across nodes. Bucket
// files are laid out one directory per partition because the reduce stage reads
// exactly one partition, and because a flat directory would hold partitions x
// shards entries.
class KmerBucketWriter {
public:
    // bufferBudgetBytes is split evenly across partitions, so memory is bounded
    // regardless of partition count. Records accumulate per partition and are
    // flushed in large contiguous appends rather than one write per k-mer.
    // partitionFrom/partitionTo restrict writing to one wave's slice of the
    // partition space. A wave re-extracts every k-mer but keeps only its own
    // slice, so peak scratch is the whole shuffle divided by the wave count --
    // which is what makes 1e12 fit a budget that cannot hold 504 TB at once.
    // Appends outside the window are dropped; the wave that owns them writes them.
    KmerBucketWriter(const std::string &dir, unsigned int partitionCount,
                     const std::string &shardId, size_t bufferBudgetBytes = 1024 * 1024 * 1024,
                     unsigned int partitionFrom = 0, unsigned int partitionTo = 0);
    ~KmerBucketWriter();

    // Thread-safe.
    void append(unsigned int partition, const KmerRecord &record);

    // Pushes every buffered record out to the operating system.
    //
    // Must be called before a work item is marked done. Buffers span items, so
    // without it an item can be recorded as complete while its k-mers are still
    // in memory; a worker that then dies loses them, and because the item is
    // already done nobody redoes it. This is the same durability the work queue
    // itself has -- both survive the process dying, neither survives the node
    // going down, which is the level at which the whole stage restarts anyway.
    void flushAll();
    // Flushes and closes every open bucket. Called by the destructor, but call it
    // explicitly to see write errors before the object goes away.
    void close();

    // Records appended by this writer. Call it when the threads are quiescent.
    uint64_t getRecordCount();

    // Creates the per-partition directories. Safe to call from many workers.
    static void createLayout(const std::string &dir, unsigned int partitionCount);
    static std::string partitionDir(const std::string &dir, unsigned int partition);

private:
    KmerBucketWriter(const KmerBucketWriter &);
    KmerBucketWriter &operator=(const KmerBucketWriter &);

    // Caller must hold mutexes[partition].
    void flush(unsigned int partition);

    std::string dir;
    std::string shardId;
    unsigned int partitionCount;
    size_t recordsPerBuffer;
    std::vector<std::vector<KmerRecord> > buffers;
    std::vector<FILE *> files;
    std::vector<std::mutex> mutexes;
    std::vector<uint64_t> recordCounts;
    unsigned int partitionFrom;
    unsigned int partitionTo;  // exclusive
};

// Reads every shard of one partition back.
class KmerBucketReader {
public:
    // Total records across all shards of the partition, so the reduce stage can
    // size its array in one allocation before reading.
    static uint64_t countRecords(const std::string &dir, unsigned int partition);

    // Appends every record of the partition to out.
    //
    // The reduce must drop exact duplicate records after sorting. A worker that
    // dies mid-item has already flushed part of that item into its own shard, and
    // the worker that redoes the item once the lease expires writes the same
    // records into a *different* shard, so the partition can legitimately contain
    // a record twice. Per-item shard files would make the map idempotent instead,
    // but that is items x partitions files -- 8e9 at 1e12 sequences -- so the
    // duplicates are removed here rather than prevented there.
    //
    // Dropping them is exact, not a heuristic: (kmer, id, pos) identifies one
    // k-mer occurrence in one sequence, so two byte-identical records can only
    // come from the same occurrence being written twice.
    static void readPartition(const std::string &dir, unsigned int partition,
                              std::vector<KmerRecord> &out);

    static std::vector<std::string> shardFiles(const std::string &dir, unsigned int partition);
};

#endif
