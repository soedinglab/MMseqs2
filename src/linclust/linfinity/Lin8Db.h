#ifndef MMSEQS_CREATELINDB_H
#define MMSEQS_CREATELINDB_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#ifdef OPENMP
#include <omp.h>
#endif

class Lin8DbIndex {
public:
    static const unsigned int RANK_BITS = 44;
    static const uint64_t MAX_RANK = (1ull << RANK_BITS) - 1;
    static const uint64_t MAX_BYTE = (1ull << 48) - 1;
    static const uint32_t MAX_FILE = (1u << 16) - 1;
    static const uint32_t MAX_SEQ_LEN = 32764;
    static const uint32_t MAX_ENTRY_LEN = 65535;

    struct LengthRange {
        uint64_t rankAndLength;
        uint64_t offsetAndFile;
        uint64_t headerByte;

        uint64_t firstRank() const { return rankAndLength & MAX_RANK; }
        uint32_t getSeqLen() const { return static_cast<uint32_t>((rankAndLength >> RANK_BITS) & 0xFFFFu); }
        uint64_t dataOffset() const { return offsetAndFile & MAX_BYTE; }
        uint32_t fileIndex() const { return static_cast<uint32_t>(offsetAndFile >> 48); }
        uint64_t headerOffset() const { return headerByte; }
    };

    Lin8DbIndex();

    void reserve(size_t rangeCount);
    void append(uint64_t rankBase, uint32_t seqLen, uint64_t byteBase, uint32_t fileIdx,
                uint64_t hdrBase);

    size_t rangeCount() const { return ranges.size(); }
    uint64_t getSize() const { return sequenceCount; }
    const LengthRange *data() const { return ranges.data(); }
    const LengthRange &operator[](size_t at) const { return ranges[at]; }

    size_t rangeIndexOf(uint64_t rank) const;
    size_t rangeIndexFrom(uint64_t rank, size_t cursor) const;

    uint32_t getSeqLen(uint64_t rank) const { return ranges[rangeIndexOf(rank)].getSeqLen(); }
    uint32_t getMaxSeqLen() const { return ranges.empty() ? 0 : ranges[0].getSeqLen(); }
    uint32_t getFileIndex(uint64_t rank) const { return ranges[rangeIndexOf(rank)].fileIndex(); }
    uint64_t getOffset(uint64_t rank) const { return offsetInRange(rangeIndexOf(rank), rank); }
    uint64_t offsetInRange(size_t range, uint64_t rank) const;
    uint64_t rankAfter(size_t range) const {
        return (range + 1 < ranges.size()) ? ranges[range + 1].firstRank() : sequenceCount;
    }

    uint64_t rankAtByte(uint64_t globalByte) const;
    uint64_t byteAtRank(uint64_t rank) const;
    uint64_t getDataSize() const { return dataBytes; }

    void write(const std::string &path) const;
    void read(const std::string &path);
    unsigned int nodeCount() const { return nodes; }
    unsigned int filesPerNode() const { return perNodeFiles; }
    unsigned int fileCount() const { return nodes * perNodeFiles; }
    void setLayout(unsigned int nodeCount, unsigned int filesPerNode) {
        nodes = nodeCount;
        perNodeFiles = filesPerNode;
    }
    void finish(uint64_t totalEntries);
    void checkLengthsDescend() const;

private:
    std::vector<LengthRange> ranges;
    std::vector<uint64_t> byteStarts;
    uint64_t sequenceCount;
    uint64_t dataBytes;
    unsigned int nodes;
    unsigned int perNodeFiles;

    void rebuildByteStarts();
};

struct __attribute__((packed)) KmerRecord {
    uint64_t low;
    uint64_t high;

    static const unsigned int BUCKET_BITS = 13;
    static const unsigned int KEY_BITS = 51 - BUCKET_BITS;
    static const unsigned int RANK_BITS = Lin8DbIndex::RANK_BITS;
    static const unsigned int POS_BITS = 15;

    static const unsigned int SUB_BUCKET_BITS = 8;
    static const size_t SUB_BUCKET_COUNT = size_t(1) << SUB_BUCKET_BITS;

    static const unsigned int ADJACENT_COUNT = 6;
    static const unsigned int ADJACENT_BITS = 5;

    static const uint64_t BUCKET_COUNT = uint64_t(1) << BUCKET_BITS;
    static const uint64_t KEY_MAX = (uint64_t(1) << KEY_BITS) - 1;
    static const uint64_t RANK_MAX = (uint64_t(1) << RANK_BITS) - 1;
    static const uint64_t POS_MAX = (uint64_t(1) << POS_BITS) - 1;
    static const uint64_t ADJACENT_MAX = (uint64_t(1) << ADJACENT_BITS) - 1;

    static const unsigned int RANK_HIGH_BITS = 64 - KEY_BITS;
    static const unsigned int RANK_LOW_BITS = RANK_BITS - RANK_HIGH_BITS;
    static const unsigned int POS_SHIFT = RANK_LOW_BITS + POS_BITS;
    static const unsigned int ADJACENT_SHIFT =
        64 - RANK_LOW_BITS - POS_BITS - ADJACENT_COUNT * ADJACENT_BITS;
    static const uint64_t ADJACENT_ALL = (uint64_t(1) << (ADJACENT_COUNT * ADJACENT_BITS)) - 1;

    static const size_t DISK_BYTES = 16;
    void pack(unsigned char *to) const { memcpy(to, this, DISK_BYTES); }
    void unpack(const unsigned char *from) { memcpy(this, from, DISK_BYTES); }
    static_assert(RANK_LOW_BITS + POS_BITS + ADJACENT_COUNT * ADJACENT_BITS <= 64,
                  "the k-mer record's second word is over budget");

    void set(uint64_t key, uint64_t rank, uint64_t pos, uint64_t adjacent) {
        low = (key << RANK_HIGH_BITS) | (rank >> RANK_LOW_BITS);
        high = ((rank & ((uint64_t(1) << RANK_LOW_BITS) - 1)) << (64 - RANK_LOW_BITS))
               | (pos << (64 - POS_SHIFT)) | (adjacent << ADJACENT_SHIFT);
    }

    uint64_t key() const { return low >> RANK_HIGH_BITS; }
    uint64_t rank() const {
        return ((low & ((uint64_t(1) << RANK_HIGH_BITS) - 1)) << RANK_LOW_BITS)
               | (high >> (64 - RANK_LOW_BITS));
    }
    uint64_t pos() const { return (high >> (64 - POS_SHIFT)) & POS_MAX; }
    uint64_t adjacent() const { return (high >> ADJACENT_SHIFT) & ADJACENT_ALL; }
    unsigned int adjacentAt(unsigned int slot) const {
        return (high >> (ADJACENT_SHIFT + slot * ADJACENT_BITS)) & ADJACENT_MAX;
    }

    unsigned int subBucket() const {
        return static_cast<unsigned int>(key() >> (KEY_BITS - SUB_BUCKET_BITS));
    }

    static bool byKeyAndRank(const KmerRecord &first, const KmerRecord &second) {
        if (first.low != second.low) {
            return first.low < second.low;
        }
        return first.high < second.high;
    }
};

struct __attribute__((packed)) PairRecord {
    uint64_t low;
    uint64_t high;

    static const unsigned int RANK_BITS = KmerRecord::RANK_BITS;
    static const unsigned int DIAGONAL_BITS = 16;
    static const int DIAGONAL_BIAS = 1 << (DIAGONAL_BITS - 1);

    static const size_t DEFAULT_REP_RANK_BLOCKS = 512;
    static const size_t MAX_REP_RANK_BLOCKS = 4096;
    static_assert(MAX_REP_RANK_BLOCKS == (size_t(1) << 12),
                  "the fine repRankBlock packing below reserves 12 bits for the repRankBlock count");

    static size_t repRankBlockOf(uint64_t rep, uint64_t ranks, size_t repRankBlocks) {
        return ranks == 0 ? 0 : std::min(rep * repRankBlocks / ranks, repRankBlocks - 1);
    }

    static uint64_t firstRankOf(size_t repRankBlock, uint64_t ranks, size_t repRankBlocks) {
        return (repRankBlock * ranks + repRankBlocks - 1) / repRankBlocks;
    }

    static const unsigned int REP_RANK_SUB_BLOCK_BITS = 8;
    static_assert(Lin8DbIndex::RANK_BITS + 12 + REP_RANK_SUB_BLOCK_BITS <= 64,
                  "a rank times the repRankBlock count times the sub repRankBlock count is over 64 bits");
    static const size_t REP_RANK_SUB_BLOCKS = size_t(1) << REP_RANK_SUB_BLOCK_BITS;
    static size_t fineOf(uint64_t rep, uint64_t ranks, size_t repRankBlocks) {
        return ranks == 0 ? 0 : std::min(rep * repRankBlocks * REP_RANK_SUB_BLOCKS / ranks,
                                         repRankBlocks * REP_RANK_SUB_BLOCKS - 1);
    }
    static size_t repRankSubBlockOf(uint64_t rep, uint64_t ranks, size_t repRankBlocks) {
        return fineOf(rep, ranks, repRankBlocks) % REP_RANK_SUB_BLOCKS;
    }

    static const unsigned int MEMBER_HIGH_BITS = 64 - RANK_BITS;
    static const unsigned int MEMBER_LOW_BITS = RANK_BITS - MEMBER_HIGH_BITS;
    static const unsigned int DIAGONAL_SHIFT = 64 - MEMBER_LOW_BITS - DIAGONAL_BITS;

    void set(uint64_t rep, uint64_t member, int diagonal) {
        low = (rep << MEMBER_HIGH_BITS) | (member >> MEMBER_LOW_BITS);
        high = ((member & ((uint64_t(1) << MEMBER_LOW_BITS) - 1)) << (64 - MEMBER_LOW_BITS))
               | (uint64_t(uint16_t(diagonal + DIAGONAL_BIAS)) << DIAGONAL_SHIFT);
    }

    uint64_t rep() const { return low >> MEMBER_HIGH_BITS; }
    uint64_t member() const {
        return ((low & ((uint64_t(1) << MEMBER_HIGH_BITS) - 1)) << MEMBER_LOW_BITS)
               | (high >> (64 - MEMBER_LOW_BITS));
    }
    int diagonal() const {
        return int((high >> DIAGONAL_SHIFT) & 0xFFFF) - DIAGONAL_BIAS;
    }

    static bool sameRepAndMember(const PairRecord &first, const PairRecord &second) {
        return first.rep() == second.rep() && first.member() == second.member();
    }

    static bool byRepAndMember(const PairRecord &first, const PairRecord &second) {
        if (first.low != second.low) {
            return first.low < second.low;
        }
        return first.high < second.high;
    }

    static const size_t DISK_BYTES = 13;
    void pack(unsigned char *to) const {
        const uint64_t top = high >> SPARE_BITS;
        memcpy(to, &low, sizeof(uint64_t));
        memcpy(to + sizeof(uint64_t), &top, DISK_BYTES - sizeof(uint64_t));
    }
    void unpack(const unsigned char *from) {
        uint64_t top = 0;
        memcpy(&low, from, sizeof(uint64_t));
        memcpy(&top, from + sizeof(uint64_t), DISK_BYTES - sizeof(uint64_t));
        high = top << SPARE_BITS;
    }

private:
    static const unsigned int SPARE_BITS = 8 * (16 - DISK_BYTES);
    static_assert(2 * RANK_BITS + DIAGONAL_BITS <= 8 * DISK_BYTES,
                  "a pair record's fields are wider than the bytes that reach the disk");
};

class KSeqWrapper;

struct InputSplit {
    size_t file;
    uint64_t from;
    uint64_t until;
    bool compressed;
    uint64_t bytes() const { return until - from; }
};

std::vector<InputSplit> planInputSplits(const std::vector<std::string> &filenames, size_t want);

class InputSplitReader {
public:
    static const size_t PACKED_READ_BYTES = 4u << 20;

    InputSplitReader(const std::string &filename, const InputSplit &chunk);
    ~InputSplitReader();

    bool next(const char *&header, size_t &headerLength, const char *&sequence, size_t &length);

private:
    InputSplitReader(const InputSplitReader &);
    InputSplitReader &operator=(const InputSplitReader &);

    bool fill();

    int fd;
    InputSplit chunk;
    std::string name;
    std::vector<char> buffer;
    size_t at;
    size_t filled;
    uint64_t readTo;
    bool started;
    std::string header;
    std::string sequence;
    KSeqWrapper *whole;
    void *stream;
    std::vector<char> packed;
    size_t packedAt;
    size_t packedFilled;
    uint64_t packedTo;
    uint64_t fileSize;
    uint64_t endsAt;
};

// an encoder turns one flush into the bytes of the bucket file and, if it needs one, an index entry
template <typename Record>
struct PackRecords {
    void operator()(std::vector<Record> &records, std::vector<unsigned char> &out, std::vector<unsigned char> &index) const {
        index.clear();
        out.resize(records.size() * Record::DISK_BYTES);
        for (size_t i = 0; i < records.size(); i++) {
            records[i].pack(&out[i * Record::DISK_BYTES]);
        }
    }
};

template <typename Record, typename Encoder = PackRecords<Record> >
class BucketWriter {
public:
    BucketWriter(const std::string &prefix, size_t buckets, unsigned int threads, size_t budget,
                 const Encoder &encoder = Encoder())
        : prefix(prefix), buckets(buckets), encoder(encoder), written(buckets, 0),
          bytes(buckets, 0), indexBytes(buckets, 0), offsets(buckets, 0), indexOffsets(buckets, 0),
          pendingIndex(buckets), staged(threads, std::vector<std::vector<Record> >(buckets)) {
        const size_t share = budget / STAGING_SHARE;
        depth = std::max<size_t>(256, (share / 4) / (threads * buckets * sizeof(Record)));
        depth = std::min<size_t>(depth, FLUSH_BYTES / sizeof(Record));
        poolDepth = std::max<size_t>(depth, (share - share / 4) / (buckets * sizeof(Record)));
        poolDepth = std::min<size_t>(poolDepth, POOL_BYTES / sizeof(Record));
        pooled.resize(buckets);
        leaving.resize(threads);
        gate.assign(buckets, 0);
        held = (threads * depth + poolDepth + depth) * buckets * sizeof(Record);
    }

    size_t bytesHeld() const { return held; }

    ~BucketWriter() { close(); }

    int openBucket(size_t bucket) const {
        const std::string path = name(bucket);
        const int fd = open(path.c_str(), O_WRONLY | O_CREAT, 0666);
        if (fd < 0) {
            Debug(Debug::ERROR) << "Cannot open " << path << " for writing, error " << errno << "\n";
            EXIT(EXIT_FAILURE);
        }
        return fd;
    }

    void closeBucket(int fd, size_t bucket) const {
        if (::close(fd) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << name(bucket) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    static const size_t DESCRIPTOR_SLACK = 512;

    size_t openableBuckets() const {
        struct rlimit limit;
        if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {
            return 0;
        }
        const rlim_t want = (rlim_t) buckets + DESCRIPTOR_SLACK;
        if (limit.rlim_cur < want) {
            limit.rlim_cur = std::min(want, limit.rlim_max);
            if (setrlimit(RLIMIT_NOFILE, &limit) != 0 || getrlimit(RLIMIT_NOFILE, &limit) != 0) {
                return 0;
            }
        }
        const size_t lent = limit.rlim_cur > DESCRIPTOR_SLACK
                                ? std::min<size_t>(buckets, (size_t) limit.rlim_cur - DESCRIPTOR_SLACK)
                                : 0;
        Debug(Debug::INFO) << "Keeping " << lent << " of " << buckets << " buckets open"
                           << (lent < buckets ? "; the rest are opened for the write" : "") << "\n";
        return lent;
    }

    int fdOf(size_t bucket) const { return kept[bucket] >= 0 ? kept[bucket] : openBucket(bucket); }

    void release(int fd, size_t bucket) const {
        if (kept[bucket] < 0) {
            closeBucket(fd, bucket);
        }
    }

    void openAt(const std::vector<uint64_t> &keepBytes, const std::vector<uint64_t> &keepIndexBytes) {
        const size_t lent = openableBuckets();
        kept.assign(buckets, -1);
        for (size_t i = 0; i < buckets; i++) {
            const int fd = openBucket(i);
            offsets[i] = keepBytes[i];
            indexOffsets[i] = keepIndexBytes[i];
            if (ftruncate(fd, static_cast<off_t>(offsets[i])) != 0) {
                Debug(Debug::ERROR) << "Cannot cut " << name(i) << " to " << offsets[i] << " byte\n";
                EXIT(EXIT_FAILURE);
            }
            if (truncate(indexName(i).c_str(), static_cast<off_t>(indexOffsets[i])) != 0
                && (errno != ENOENT || indexOffsets[i] > 0)) {
                Debug(Debug::ERROR) << "Cannot cut " << indexName(i) << " to " << indexOffsets[i] << " byte\n";
                EXIT(EXIT_FAILURE);
            }
            if (i < lent) {
                kept[i] = fd;
            } else {
                closeBucket(fd, i);
            }
            written[i] = 0;
        }
    }

    void add(unsigned int thread, const Record &record, size_t bucket) {
        std::vector<Record> &buffer = staged[thread][bucket];
        if (buffer.capacity() < depth) {
            buffer.reserve(depth);
        }
        buffer.push_back(record);
        if (buffer.size() >= depth) {
            drain(bucket, buffer);
        }
    }

    void flushAll(unsigned int threads) {
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
        for (size_t bucket = 0; bucket < buckets; bucket++) {
            for (size_t thread = 0; thread < staged.size(); thread++) {
                drain(bucket, staged[thread][bucket]);
            }
            flush(bucket, pooled[bucket]);
            writeIndex(bucket);
        }
    }

    void finishChunk(unsigned int threads) {
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
        for (size_t bucket = 0; bucket < buckets; bucket++) {
            for (size_t thread = 0; thread < staged.size(); thread++) {
                drain(bucket, staged[thread][bucket]);
            }
            flush(bucket, pooled[bucket]);
            writeIndex(bucket);
            if (written[bucket] > 0) {
                const int fd = fdOf(bucket);
                sync_file_range(fd, 0, 0, SYNC_FILE_RANGE_WRITE);
                release(fd, bucket);
            }
        }
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
        for (size_t bucket = 0; bucket < buckets; bucket++) {
            if (written[bucket] == 0) {
                continue;
            }
            const int fd = fdOf(bucket);
            if (fdatasync(fd) != 0) {
                Debug(Debug::ERROR) << "Cannot flush " << name(bucket) << " to storage\n";
                EXIT(EXIT_FAILURE);
            }
            release(fd, bucket);
        }
    }

    void close() {
        for (size_t i = 0; i < kept.size(); i++) {
            if (kept[i] >= 0) {
                closeBucket(kept[i], i);
                kept[i] = -1;
            }
        }
    }

    const std::vector<uint64_t> &chunkCounts() const { return written; }
    const std::vector<uint64_t> &chunkBytes() const { return bytes; }
    const std::vector<uint64_t> &chunkIndexBytes() const { return indexBytes; }
    void resetCounts() {
        std::fill(written.begin(), written.end(), 0);
        std::fill(bytes.begin(), bytes.end(), 0);
        std::fill(indexBytes.begin(), indexBytes.end(), 0);
    }

private:
    std::string name(size_t bucket) const { return prefix + "." + SSTR(bucket); }
    std::string indexName(size_t bucket) const { return name(bucket) + ".idx"; }

    void drain(size_t bucket, std::vector<Record> &buffer) {
        if (buffer.empty()) {
            return;
        }
        while (__sync_lock_test_and_set(&gate[bucket], 1) != 0) {
            while (__atomic_load_n(&gate[bucket], __ATOMIC_RELAXED) != 0) {
            }
        }
        std::vector<Record> &pool = pooled[bucket];
        if (pool.capacity() < poolDepth + depth) {
            pool.reserve(poolDepth + depth);
        }
        pool.insert(pool.end(), buffer.begin(), buffer.end());
        size_t here = 0;
#ifdef OPENMP
        here = (size_t) omp_get_thread_num();
#endif
        if (here >= leaving.size()) {
            if (pool.size() >= poolDepth) {
                flush(bucket, pool);
            }
            __sync_lock_release(&gate[bucket]);
            buffer.clear();
            return;
        }
        std::vector<Record> &outbound = leaving[here];
        if (pool.size() >= poolDepth) {
            outbound.swap(pool);
        }
        __sync_lock_release(&gate[bucket]);
        buffer.clear();
        flush(bucket, outbound);
    }

    void flush(size_t bucket, std::vector<Record> &buffer) {
        if (buffer.empty()) {
            return;
        }
        std::vector<unsigned char> packed;
        std::vector<unsigned char> index;
        encoder(buffer, packed, index);
        const size_t size = packed.size();
        const uint64_t at = __sync_fetch_and_add(&offsets[bucket], size);
        __sync_fetch_and_add(&written[bucket], buffer.size());
        __sync_fetch_and_add(&bytes[bucket], size);
        const int fd = fdOf(bucket);
        const ssize_t wrote = pwrite(fd, packed.data(), size, static_cast<off_t>(at));
        if (wrote < 0 || static_cast<size_t>(wrote) != size) {
            Debug(Debug::ERROR) << "Cannot write " << size << " byte to " << name(bucket) << "\n";
            EXIT(EXIT_FAILURE);
        }
        release(fd, bucket);
        buffer.clear();
        if (index.empty() == false) {
            const uint64_t where[2] = {at, size};
#pragma omp critical(bucket_index)
            {
                std::vector<unsigned char> &pending = pendingIndex[bucket];
                pending.insert(pending.end(), reinterpret_cast<const unsigned char *>(where),
                               reinterpret_cast<const unsigned char *>(where) + sizeof(where));
                pending.insert(pending.end(), index.begin(), index.end());
            }
        }
    }

    // index entries wait in memory and reach the .idx file once a chunk, so a flush opens one file
    void writeIndex(size_t bucket) {
        std::vector<unsigned char> &pending = pendingIndex[bucket];
        if (pending.empty()) {
            return;
        }
        const int fd = open(indexName(bucket).c_str(), O_WRONLY | O_CREAT, 0666);
        const ssize_t wrote = fd < 0 ? -1 : pwrite(fd, pending.data(), pending.size(), static_cast<off_t>(indexOffsets[bucket]));
        if (wrote < 0 || static_cast<size_t>(wrote) != pending.size() || fdatasync(fd) != 0 || ::close(fd) != 0) {
            Debug(Debug::ERROR) << "Cannot write " << pending.size() << " byte to " << indexName(bucket) << "\n";
            EXIT(EXIT_FAILURE);
        }
        indexOffsets[bucket] += pending.size();
        indexBytes[bucket] += pending.size();
        pending.clear();
    }

    static const size_t STAGING_SHARE = 8;
    static const size_t FLUSH_BYTES = 4 * 1024 * 1024;
    static const size_t POOL_BYTES = 16 * 1024 * 1024;
    std::string prefix;
    size_t buckets;
    Encoder encoder;
    size_t depth;
    size_t poolDepth;
    std::vector<std::vector<Record> > pooled;
    std::vector<std::vector<Record> > leaving;
    std::vector<int> gate;
    mutable std::vector<int> kept;
    size_t held;
    std::vector<uint64_t> written;
    std::vector<uint64_t> bytes;
    std::vector<uint64_t> indexBytes;
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> indexOffsets;
    std::vector<std::vector<unsigned char> > pendingIndex;
    std::vector<std::vector<std::vector<Record> > > staged;
};

template <class Record>
size_t writeRecords(const Record *from, size_t count, FILE *out) {
    if (Record::DISK_BYTES == sizeof(Record)) {
        return fwrite(from, sizeof(Record), count, out);
    }
    thread_local std::vector<unsigned char> packed;
    packed.resize(count * Record::DISK_BYTES);
    for (size_t i = 0; i < count; i++) {
        from[i].pack(&packed[i * Record::DISK_BYTES]);
    }
    return fwrite(packed.data(), Record::DISK_BYTES, count, out);
}

template <class Record>
size_t readRecords(Record *into, size_t count, FILE *in) {
    if (Record::DISK_BYTES == sizeof(Record)) {
        return fread(into, sizeof(Record), count, in);
    }
    thread_local std::vector<unsigned char> packed;
    packed.resize(count * Record::DISK_BYTES);
    const size_t read = fread(packed.data(), Record::DISK_BYTES, count, in);
    for (size_t i = 0; i < read; i++) {
        into[i].unpack(&packed[i * Record::DISK_BYTES]);
    }
    return read;
}

void publishAllAtomically(std::vector<std::pair<std::string, std::string> > &pending,
                          unsigned int threads);

void requireMemory(const std::string &what, size_t bytes, size_t budget,
                  const std::string &narrower);
static const uint64_t PUBLISH_BATCH_BYTES = 4ull * 1024 * 1024 * 1024;
static const size_t PUBLISH_BATCH_FILES = 64;

void writeBucketManifest(const std::string &path, const std::vector<uint64_t> &counts,
                         const std::vector<uint64_t> &bytes, const std::vector<uint64_t> &indexBytes,
                         const std::string &spanKey, uint64_t spanBegin, uint64_t spanEnd);

inline std::string uniqueTmpSuffix() {
    char host[HOST_NAME_MAX + 1];
    memset(host, 0, sizeof(host));
    gethostname(host, HOST_NAME_MAX);
    return std::string(host) + "." + SSTR(getpid());
}

inline std::string nodeDonePath(const std::string &path, unsigned int node) {
    return path + "." + SSTR(node) + ".done";
}

// the alignments beside a pair file, one text line per pair, accepted ones under a "#\trep" line
inline std::string alnTextPath(const std::string &prefix, unsigned int node, size_t block) {
    return prefix + "_text." + SSTR(node) + "." + SSTR(block);
}
static const char ALN_TEXT_REP_MARK = '#';

class AlnTextReader {
public:
    AlnTextReader(const std::string &prefix, unsigned int node, size_t block, bool required)
        : path(alnTextPath(prefix, node, block)), line(NULL), cap(0) {
        file = fopen(path.c_str(), "r");
        if (file == NULL && required) {
            Debug(Debug::ERROR) << "Cannot open " << path
                                << ", which the aligning pass should have written\n";
            EXIT(EXIT_FAILURE);
        }
        if (file != NULL) {
            setvbuf(file, NULL, _IOFBF, 1u << 20);
        }
    }

    ~AlnTextReader() {
        if (file != NULL) {
            fclose(file);
        }
        free(line);
    }

    bool next(char *&begin, size_t &length) {
        if (file == NULL) {
            return false;
        }
        const ssize_t got = getline(&line, &cap, file);
        if (got > 0) {
            begin = line;
            length = (size_t) got;
            return true;
        }
        if (ferror(file) != 0) {
            Debug(Debug::ERROR) << "Cannot read " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        return false;
    }

    const std::string &name() const { return path; }

private:
    AlnTextReader(const AlnTextReader &);
    AlnTextReader &operator=(const AlnTextReader &);

    std::string path;
    FILE *file;
    char *line;
    size_t cap;
};
void markNodeDone(const std::string &path, unsigned int node);

// a counter rewritten in place, because a new name on NFS is a negative dentry another node caches
void publishProgress(const std::string &path, uint64_t value);

void removeConsumedBuckets(const std::string &prefix, unsigned int nodes, size_t from, size_t until,
                  size_t step);
void requireEveryNodeDone(const std::string &path, unsigned int nodes);
void waitNodeDone(const std::string &path, unsigned int node, unsigned int limitSeconds);
void waitEveryNodeDone(const std::string &path, unsigned int nodes, unsigned int limitSeconds);

void writeSubBucketCounts(const std::string &path, const std::vector<uint64_t> &base,
                          const std::vector<std::vector<uint64_t> > &perThread, size_t chunks);
std::vector<uint64_t> readSubBucketCounts(const std::string &path, size_t entries, size_t &chunks);

class SubBucketCounts {
public:
    SubBucketCounts(const std::string &prefix, unsigned int nodes, size_t subBuckets,
                 size_t buckets);
    ~SubBucketCounts();
    std::vector<uint64_t> of(size_t bucket) const;
    std::vector<uint64_t> of(size_t bucket, size_t node) const;

private:
    SubBucketCounts(const SubBucketCounts &);
    SubBucketCounts &operator=(const SubBucketCounts &);
    std::vector<FILE *> files;
    std::vector<std::string> paths;
    size_t subBuckets;
};

size_t readBucketManifests(const std::string &prefix, size_t chunks, std::vector<uint64_t> &bytesInto,
                           std::vector<uint64_t> &indexBytesInto, uint64_t *resumeAt = NULL);

template <typename T>
class RawArray {
public:
    RawArray() : items(NULL), count(0) {}
    ~RawArray() { delete[] items; }

    void resize(size_t n) {
        delete[] items;
        items = new T[n];
        count = n;
    }

    T *begin() const { return items; }
    T &operator[](size_t i) const { return items[i]; }
    size_t size() const { return count; }

private:
    RawArray(const RawArray &);
    RawArray &operator=(const RawArray &);
    T *items;
    size_t count;
};

void preadFully(int fd, void *into, size_t bytes, uint64_t at, const std::string &what);

// a segment is one sorted flush, bit-packed in blocks that never cross a sub-bucket, sliced by subStart
struct SegmentHeader {
    static const size_t SUBS = 256;
    static const size_t BLOCK_RECORDS = 128;

    uint64_t magic;
    uint64_t records;
    uint32_t rankBits;
    uint32_t classCount;
    uint64_t subStart[SUBS + 1];
};
static_assert(KmerRecord::SUB_BUCKET_COUNT == SegmentHeader::SUBS && PairRecord::REP_RANK_SUB_BLOCKS == SegmentHeader::SUBS,
              "segments are sliced by the same 256 sub-buckets both record types are counted in");

// a block's rows share a key column: a flag per row after the first says whether the key changed
struct BlockHeader {
    uint8_t records;
    uint8_t newKeys;
    uint8_t keyBits;
    uint8_t widthBits;

    size_t payloadBytes(unsigned int firstKeyBits, size_t perRecordBits) const {
        const size_t bits = (records - 1) + firstKeyBits + (size_t) (newKeys - 1) * keyBits + records * perRecordBits;
        return (bits + 7) / 8;
    }
};

inline unsigned int bitsFor(uint64_t value) {
    unsigned int bits = 0;
    while (value != 0) {
        bits++;
        value >>= 1;
    }
    return bits;
}

inline uint64_t zigzag(int64_t value) { return value < 0 ? uint64_t(-value) * 2 - 1 : uint64_t(value) * 2; }
inline int64_t unzigzag(uint64_t value) { return (value & 1) ? -int64_t((value + 1) / 2) : int64_t(value / 2); }

class BitWriter {
public:
    BitWriter(std::vector<unsigned char> &out) : out(out), acc(0), filled(0) {}
    void put(uint64_t value, unsigned int bits) {
        acc |= value << filled;
        filled += bits;
        while (filled >= 8) {
            out.push_back(static_cast<unsigned char>(acc));
            acc >>= 8;
            filled -= 8;
        }
    }
    void finish() {
        if (filled > 0) {
            out.push_back(static_cast<unsigned char>(acc));
            acc = 0;
            filled = 0;
        }
    }
private:
    std::vector<unsigned char> &out;
    uint64_t acc;
    unsigned int filled;
};

class BitReader {
public:
    BitReader(const unsigned char *at) : at(at), acc(0), filled(0) {}
    uint64_t get(unsigned int bits) {
        while (filled < bits) {
            acc |= static_cast<uint64_t>(*at++) << filled;
            filled += 8;
        }
        const uint64_t value = acc & ((uint64_t(1) << bits) - 1);
        acc >>= bits;
        filled -= bits;
        return value;
    }
private:
    const unsigned char *at;
    uint64_t acc;
    unsigned int filled;
};

inline void countKeys(const uint64_t *keys, size_t n, BlockHeader &block) {
    block.records = static_cast<uint8_t>(n);
    block.newKeys = 1;
    uint64_t widest = 0;
    for (size_t i = 1; i < n; i++) {
        if (keys[i] != keys[i - 1]) {
            block.newKeys++;
            widest = std::max(widest, keys[i] - keys[i - 1]);
        }
    }
    block.keyBits = static_cast<uint8_t>(bitsFor(widest));
}

inline void putKeys(BitWriter &bits, const uint64_t *keys, size_t n, unsigned int firstBits, unsigned int deltaBits) {
    for (size_t i = 1; i < n; i++) {
        bits.put(keys[i] != keys[i - 1], 1);
    }
    bits.put(keys[0], firstBits);
    for (size_t i = 1; i < n; i++) {
        if (keys[i] != keys[i - 1]) {
            bits.put(keys[i] - keys[i - 1], deltaBits);
        }
    }
}

inline void getKeys(BitReader &bits, uint64_t *keys, size_t n, unsigned int firstBits, unsigned int deltaBits) {
    bool fresh[SegmentHeader::BLOCK_RECORDS];
    for (size_t i = 1; i < n; i++) {
        fresh[i] = bits.get(1) != 0;
    }
    keys[0] = bits.get(firstBits);
    for (size_t i = 1; i < n; i++) {
        keys[i] = keys[i - 1] + (fresh[i] ? bits.get(deltaBits) : 0);
    }
}

inline void putBlockHeader(std::vector<unsigned char> &out, const BlockHeader &block) {
    out.insert(out.end(), reinterpret_cast<const unsigned char *>(&block),
               reinterpret_cast<const unsigned char *>(&block) + sizeof(block));
}

// call at every block start and once at the end to fill header.subStart
inline void markSubStart(SegmentHeader &header, unsigned int &sub, unsigned int upTo, size_t at) {
    while (sub <= upTo) {
        header.subStart[sub++] = at;
    }
}

template <class Codec>
struct SegmentEncoder {
    Codec codec;
    SegmentEncoder(const Codec &codec) : codec(codec) {}
    void operator()(std::vector<typename Codec::Record> &records, std::vector<unsigned char> &out,
                    std::vector<unsigned char> &index) const {
        std::sort(records.begin(), records.end(), Codec::less);
        SegmentHeader header;
        codec.encode(records, out, header);
        index.assign(reinterpret_cast<const unsigned char *>(&header),
                     reinterpret_cast<const unsigned char *>(&header) + sizeof(header));
    }
};

struct Segment {
    int fd;
    uint64_t at;
    uint64_t bytes;
    SegmentHeader header;
};

// the segments one bucket holds across the writer nodes, read from the .idx files
template <class Codec>
class BucketSegments {
public:
    BucketSegments(const std::string &prefix, unsigned int nodes, size_t bucket) : fds(nodes, -1) {
        for (unsigned int node = 0; node < nodes; node++) {
            const std::string path = prefix + "." + SSTR(node) + "." + SSTR(bucket);
            fds[node] = open(path.c_str(), O_RDONLY);
            if (fds[node] < 0) {
                Debug(Debug::ERROR) << "Cannot open " << path << ", which " << Codec::PRODUCER << " wrote and this"
                                    << " pass has not read yet. With --remove-tmp-files it drops them as it"
                                    << " consumes them, so rerun " << Codec::PRODUCER << " to make them again\n";
                EXIT(EXIT_FAILURE);
            }
            const std::string indexPath = path + ".idx";
            const size_t entry = 2 * sizeof(uint64_t) + sizeof(SegmentHeader);
            const uint64_t indexSize = FileUtil::fileExists(indexPath.c_str()) ? FileUtil::getFileSize(indexPath) : 0;
            const uint64_t bodySize = FileUtil::getFileSize(path);
            if (indexSize % entry != 0 || (indexSize == 0 && bodySize != 0)) {
                Debug(Debug::ERROR) << indexPath << " does not index " << path << ". Was it written by an older "
                                    << Codec::PRODUCER << ", or one still running?\n";
                EXIT(EXIT_FAILURE);
            }
            std::vector<unsigned char> index(indexSize);
            if (indexSize > 0) {
                const int indexFd = open(indexPath.c_str(), O_RDONLY);
                preadFully(indexFd, index.data(), indexSize, 0, indexPath);
                close(indexFd);
            }
            for (size_t at = 0; at < indexSize; at += entry) {
                Segment segment;
                segment.fd = fds[node];
                memcpy(&segment.at, &index[at], sizeof(uint64_t));
                memcpy(&segment.bytes, &index[at + sizeof(uint64_t)], sizeof(uint64_t));
                memcpy(&segment.header, &index[at + 2 * sizeof(uint64_t)], sizeof(SegmentHeader));
                if (segment.header.magic != Codec::MAGIC || segment.at + segment.bytes > bodySize
                    || segment.header.subStart[SegmentHeader::SUBS] != segment.bytes) {
                    Debug(Debug::ERROR) << indexPath << " entry " << (at / entry) << " does not fit " << path
                                        << ". Was " << Codec::PRODUCER << " still running?\n";
                    EXIT(EXIT_FAILURE);
                }
                segments.push_back(segment);
            }
        }
    }
    ~BucketSegments() {
        for (size_t i = 0; i < fds.size(); i++) {
            close(fds[i]);
        }
    }

    uint64_t packedBytes(size_t sub) const {
        uint64_t bytes = 0;
        for (size_t r = 0; r < segments.size(); r++) {
            bytes += segments[r].header.subStart[sub + 1] - segments[r].header.subStart[sub];
        }
        return bytes;
    }

    std::vector<Segment> segments;

private:
    std::vector<int> fds;
};

// cuts the sub-buckets into ranges whose packed bytes, records and merge scratch fit the budget
template <class Codec>
std::vector<size_t> planWindows(const BucketSegments<Codec> &bucketSegments, const std::vector<uint64_t> &subCounts,
                               size_t budget, unsigned int threads, const std::string &what, const char *narrower) {
    typedef typename Codec::Record Record;
    std::vector<size_t> cuts(1, 0);
    size_t widest = 0;
    for (size_t sub = 0; sub < SegmentHeader::SUBS; sub++) {
        widest = std::max<size_t>(widest, subCounts[sub]);
    }
    size_t held = 0;
    size_t heldRecords = 0;
    for (size_t sub = 0; sub < SegmentHeader::SUBS; sub++) {
        const size_t need = subCounts[sub] * sizeof(Record) + bucketSegments.packedBytes(sub);
        requireMemory(what + " prefix " + SSTR(sub), need + subCounts[sub] * sizeof(Record), budget, narrower);
        const size_t scratch = std::min<size_t>((size_t) threads * widest, heldRecords + subCounts[sub]) * sizeof(Record);
        if (held + need + scratch > budget) {
            cuts.push_back(sub);
            held = 0;
            heldRecords = 0;
        }
        held += need;
        heldRecords += subCounts[sub];
    }
    cuts.push_back(size_t(SegmentHeader::SUBS));
    return cuts;
}

// reads one slice of every segment for [from, to) and merges each sub-bucket into the order the sort made
template <class Codec>
void loadWindow(const BucketSegments<Codec> &bucketSegments, const std::vector<uint64_t> &subCounts, size_t from, size_t to,
               unsigned int threads, const std::string &what, RawArray<typename Codec::Record> &into) {
    typedef typename Codec::Record Record;
    const std::vector<Segment> &segments = bucketSegments.segments;
    std::vector<size_t> starts(to - from + 1, 0);
    for (size_t sub = from; sub < to; sub++) {
        starts[sub - from + 1] = starts[sub - from] + subCounts[sub];
    }
    into.resize(starts.back());

    std::vector<std::vector<unsigned char> > packed(segments.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t r = 0; r < segments.size(); r++) {
        const SegmentHeader &header = segments[r].header;
        packed[r].resize(header.subStart[to] - header.subStart[from]);
        preadFully(segments[r].fd, packed[r].data(), packed[r].size(), segments[r].at + header.subStart[from], what);
    }

#pragma omp parallel num_threads(threads)
    {
        std::vector<Record> scratch;
        std::vector<size_t> begin(segments.size()), end(segments.size()), heap;
        // the heap keeps the segment whose next record sorts first on top
        const auto later = [&](size_t a, size_t b) { return Codec::less(scratch[begin[b]], scratch[begin[a]]); };
#pragma omp for schedule(dynamic, 1)
        for (size_t sub = from; sub < to; sub++) {
            const size_t want = subCounts[sub];
            if (want == 0) {
                continue;
            }
            size_t got = 0;
            heap.clear();
            scratch.resize(want);
            for (size_t r = 0; r < segments.size(); r++) {
                const SegmentHeader &header = segments[r].header;
                const unsigned char *slice = packed[r].data() + (header.subStart[sub] - header.subStart[from]);
                begin[r] = got;
                got += Codec::decode(slice, slice + (header.subStart[sub + 1] - header.subStart[sub]), header,
                                     scratch.data() + got, want - got);
                end[r] = got;
                if (end[r] > begin[r]) {
                    heap.push_back(r);
                }
            }
            if (got != want) {
                Debug(Debug::ERROR) << what << " prefix " << sub << " holds " << got << " records and the counts say "
                                    << want << ". Was " << Codec::PRODUCER << " still running?\n";
                EXIT(EXIT_FAILURE);
            }
            Record *out = into.begin() + starts[sub - from];
            std::make_heap(heap.begin(), heap.end(), later);
            for (size_t i = 0; i < want; i++) {
                std::pop_heap(heap.begin(), heap.end(), later);
                const size_t r = heap.back();
                out[i] = scratch[begin[r]++];
                if (begin[r] == end[r]) {
                    heap.pop_back();
                } else {
                    std::push_heap(heap.begin(), heap.end(), later);
                }
            }
        }
    }
}


// the six flanking classes as one base-alphabet integer
unsigned int adjacentBitsFor(unsigned int classCount);

// k-mer records by key and rank, the key as a column and the rest bit-packed per block
struct KmerSegmentCodec {
    typedef KmerRecord Record;
    static const uint64_t MAGIC = 0x4C494E384B52554Eull;
    static constexpr const char *PRODUCER = "lin8-extractkmers";

    unsigned int rankBits;
    unsigned int classCount;
    KmerSegmentCodec(unsigned int rankBits, unsigned int classCount) : rankBits(rankBits), classCount(classCount) {}

    static bool less(const KmerRecord &a, const KmerRecord &b) { return KmerRecord::byKeyAndRank(a, b); }
    void encode(const std::vector<KmerRecord> &records, std::vector<unsigned char> &out, SegmentHeader &header) const;
    static size_t decode(const unsigned char *from, const unsigned char *to, const SegmentHeader &header,
                         KmerRecord *into, size_t capacity);
};

// pairs by representative and member, the rep as a column and the rest bit-packed per block
struct PairSegmentCodec {
    typedef PairRecord Record;
    static const uint64_t MAGIC = 0x4C494E385052554Eull;
    static constexpr const char *PRODUCER = "lin8-assignedpairs";

    unsigned int rankBits;
    uint64_t ranks;
    size_t repRankBlocks;
    PairSegmentCodec(unsigned int rankBits, uint64_t ranks, size_t repRankBlocks)
        : rankBits(rankBits), ranks(ranks), repRankBlocks(repRankBlocks) {}

    static bool less(const PairRecord &a, const PairRecord &b) { return PairRecord::byRepAndMember(a, b); }
    void encode(const std::vector<PairRecord> &records, std::vector<unsigned char> &out, SegmentHeader &header) const;
    static size_t decode(const unsigned char *from, const unsigned char *to, const SegmentHeader &header,
                         PairRecord *into, size_t capacity);
};

#endif
