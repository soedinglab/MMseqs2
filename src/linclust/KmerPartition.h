#ifndef MMSEQS_KMERPARTITION_H
#define MMSEQS_KMERPARTITION_H

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

#include "LengthRankTable.h"
#include "VarintCodec.h"

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

// Sub-slice of a partition, for a reduce that cannot hold the whole thing.
//
// P is derived from the *average* partition, and a partition that exceeds a
// worker's memory has to be groupable anyway: nodes in an allocation are not
// always identical, a cgroup limit can be lower than the flag says, a resumed run
// can land on smaller machines, and k-mer skew can push one partition above its
// peers. Failing there wastes the whole run.
//
// Slicing is exact for the same reason partitioning is. The slice is a pure
// function of the k-mer, so every occurrence of a k-mer lands in one slice and a
// group is never split -- and the group is the atomic unit of the reduce, since
// assignGroup only ever reads and swaps within one group's index range and
// buildThreadOffsets refuses to cut a group across threads. Grouping the slices
// separately therefore produces exactly the edges grouping the partition whole
// would, and the align stage sums their support the same way it already sums
// across partitions.
//
// A different mix from the partitioner's: that one uses the low bits of the
// 16-bit hash kmermatcher already computed, so reusing it would put every record
// of a partition in one slice.
//
// The case this does NOT cover is a single k-mer *group* too large for the
// budget. No slicing can help there -- the group cannot be divided without
// changing the answer. Measured on 10M MGnify sequences the largest group is
// 3,718 records (0.11% of its partition, ~171 KB resident) against a mean of
// 1.59, so that case is orders of magnitude away; it is diagnosed rather than
// handled.
inline unsigned int kmerSliceOf(uint64_t kmer, unsigned int sliceCount) {
    if (sliceCount <= 1) {
        return 0;
    }
    uint64_t mixed = kmer * 0x9E3779B97F4A7C15ULL;
    mixed ^= mixed >> 29;
    mixed *= 0xBF58476D1CE4E5B9ULL;
    mixed ^= mixed >> 32;
    return static_cast<unsigned int>(mixed % sliceCount);
}

// One k-mer occurrence, as held in memory.
//
// 24 bytes, packed. Two fields are here specifically to delete arrays that stock
// kmermatcher sizes by *key space* rather than by entry count, and which
// therefore cannot exist at 1e12:
//   - seqLen removes `seqkey_to_len[dbKeySize]` (kmermatcher.cpp:1236),
//   - id is carried explicitly so nothing needs a dense key-indexed side table.
// `countTable` (:1256) and `repSequence` (:1357) fall the same way.
//
// **This is the in-memory form, no longer the on-disk one.** Bucket files hold
// the packed block encoding below, which is ~13.5 B per record against these 24.
// `seqLen` in particular is never written: keys are length-ranked, so
// LengthRankTable recovers it from the id for the price of a binary search over a
// table of distinct lengths. Keeping the field resident costs nothing (the decode
// buffer is bounded) and keeps kmerRecordToPosition and everything downstream
// unchanged.
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
//
// workerCount is how many workers the runner launched for the stage, and it only
// ever raises P. A partition is the reduce's unit of *work* as well as its unit of
// memory, so deriving P from workerMemoryBytes alone makes it a ceiling on how
// many workers the stage can use, and a generous per-node limit then starves the
// very allocation it was raised for. 0 means "not told", which keeps the
// memory-only sizing and the previous behaviour.
KmerShuffleSizing deriveKmerShuffleSizing(uint64_t sequenceCount,
                                          unsigned int kmersPerSequence,
                                          uint64_t scratchBudgetBytes,
                                          uint64_t persistentBytes,
                                          uint64_t workerMemoryBytes,
                                          unsigned int workerCount = 0);

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

// Framing for one flushed run of k-mer records.
//
// Three things this buys, beyond the packing itself:
//
//   - a torn tail from a worker killed mid-flush is *recognisable*. The reader
//     rounds down to the last whole block and carries on, which is the behaviour
//     the fixed-width format could only approximate by rounding the file size
//     down to a record boundary -- and could not do at all once records vary in
//     width.
//   - the record count is readable without decoding, so the reduce can size its
//     array in one allocation by scanning headers rather than payloads.
//   - the encoding is named per block, so `--raw-records` is a writer-side switch
//     that costs the reader nothing: a directory holding both kinds still reads
//     correctly, which is what makes the two comparable in one run.
struct __attribute__((__packed__)) KmerBlockHeader {
    // Distinguishes a real header from the tail of a torn write.
    uint32_t magic;
    // ENCODING_RAW or ENCODING_PACKED.
    uint8_t encoding;
    // Bytes per k-mer field, packed mode only. Derived from the largest k-mer in
    // *this* block rather than from (k, alphabetSize): it needs no parameter
    // plumbing, it cannot overflow for a large nucleotide k, and it adapts when a
    // block happens to hold only small k-mers.
    uint8_t kmerWidth;
    uint16_t reserved;
    uint64_t recordCount;
    uint64_t payloadBytes;
    // FNV-1a over the payload, computed while the bytes are already in hand. A
    // truncated block is caught by payloadBytes alone; this also catches a block
    // boundary that has drifted, which would otherwise decode into plausible
    // nonsense.
    uint32_t checksum;
    uint32_t reserved2;
};

namespace KmerBlockCodec {

const uint32_t MAGIC = 0x4B424C4BU;  // "KBLK"
const uint8_t ENCODING_RAW = 0;
const uint8_t ENCODING_PACKED = 1;

// Six adjacent residues at five bits each. Stock uses UCHAR_MAX as the "no
// adjacency" sentinel that assignGroup tests, so it maps onto the one value five
// bits leaves spare.
const uint8_t ADJACENT_SENTINEL = 31;
const unsigned int ADJACENT_COUNT = 6;
const unsigned int ADJACENT_BITS = 5;
const unsigned int ADJACENT_BYTES = 4;

// Ceiling on the bytes one record can encode to: a 10-byte id delta, an 8-byte
// k-mer, a 3-byte position and the packed adjacency. The writer sizes its buffer
// budget against this so the encode buffer and the record buffer together stay
// inside it.
const size_t MAX_ENCODED_BYTES_PER_RECORD = VarintCodec::MAX_BYTES + 8 + 3 + ADJACENT_BYTES;

// Largest residue index the packing can carry. Checked against the alphabet at
// startup rather than per record.
const unsigned int MAX_ALPHABET_SIZE = 31;

// Sorts records by (id, kmer, pos) and appends one framed block to out.
//
// The sort is what makes the id field cost about one byte: the map scans a
// contiguous key range, so a flush holds ids from that range and, once ordered,
// consecutive deltas are the range divided by the record count. Threads share the
// per-partition buffer (that is deliberate -- a per-thread buffer would divide the
// write size by the thread count), so the ordering has to be restored here rather
// than relied on. Measured against the whole pipeline the sort is a rounding
// error; the field it shrinks is 6 bytes on every one of 2.1e13 records at 1e12.
void encode(std::vector<KmerRecord> &records, bool raw, std::vector<uint8_t> &out);

// Decodes one block payload, appending to out.
//
// `lengths` may be NULL only for a raw block, which carries seqLen itself. That
// asymmetry is deliberate: it makes --raw-records a control that does not depend
// on the length table, so comparing the two paths tests the table as well as the
// codec.
//
// Returns false on a corrupt payload rather than exiting, because the caller's
// correct response is to stop at the last good block, not to fail the run.
bool decode(const KmerBlockHeader &header, const uint8_t *payload,
            const LengthRankTable *lengths, std::vector<KmerRecord> &out);

// Whether a header could plausibly start a block. Cheap pre-check before
// trusting payloadBytes.
bool headerLooksValid(const KmerBlockHeader &header);

uint32_t checksum(const uint8_t *data, size_t size);

}  // namespace KmerBlockCodec

// Reads framed blocks from one shard file.
class KmerShardReader {
public:
    explicit KmerShardReader(const std::string &path);
    ~KmerShardReader();

    // Appends the next block's records to out. Returns false at end of file or at
    // the first block that does not decode, which is how a torn tail ends the
    // shard without failing the run.
    bool next(std::vector<KmerRecord> &out, const LengthRankTable *lengths);

    // True if the shard ended on a partial or corrupt block rather than cleanly.
    bool endedTorn() const { return torn; }

    // Sums recordCount over the shard's block headers, seeking past each payload.
    // O(blocks), not O(records): a shard holds ~1e7 records in a few hundred
    // blocks at the target scale, so this is a few hundred forward seeks.
    static uint64_t countRecords(const std::string &path, bool &torn);

private:
    KmerShardReader(const KmerShardReader &);
    KmerShardReader &operator=(const KmerShardReader &);

    FILE *file;
    std::string path;
    std::vector<uint8_t> payload;
    bool torn;
};

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
                     unsigned int partitionFrom = 0, unsigned int partitionTo = 0,
                     bool rawRecords = false);
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
    void closeFile(unsigned int partition);

    std::string dir;
    std::string shardId;
    unsigned int partitionCount;
    size_t recordsPerBuffer;
    std::vector<std::vector<KmerRecord> > buffers;
    // One reusable encode target per partition, so a flush allocates nothing.
    // Per partition rather than shared because flush() runs under the
    // per-partition mutex and a shared buffer would serialise every partition
    // behind one lock.
    std::vector<std::vector<uint8_t> > encodeBuffers;
    std::vector<FILE *> files;
    std::vector<std::mutex> mutexes;
    std::vector<uint64_t> recordCounts;
    unsigned int partitionFrom;
    unsigned int partitionTo;  // exclusive
    // Writes fixed-width records instead of packed blocks. The escape hatch that
    // separates a codec defect from a semantic one: a run with it on must produce
    // the same candidate edges as a run with it off, and if it does not, the
    // difference is in the encoding rather than in the clustering.
    bool rawRecords;
    // Descriptors are kept open across flushes up to this many; past it a flush
    // reverts to open-append-close. Counted atomically because the count is shared
    // across the per-partition mutexes.
    size_t descriptorBudget;
    std::atomic<size_t> openFiles;
};

// Reads every shard of one partition back.
class KmerBucketReader {
public:
    // Total records across all shards of the partition, so the reduce stage can
    // size its array in one allocation before reading. Reads block headers only.
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
    // lengths is required unless every block was written with --raw-records.
    static void readPartition(const std::string &dir, unsigned int partition,
                              std::vector<KmerRecord> &out, const LengthRankTable *lengths);

    static std::vector<std::string> shardFiles(const std::string &dir, unsigned int partition);
};

#endif
