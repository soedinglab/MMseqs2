#include "KmerPartition.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "VarintCodec.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstring>

#include <dirent.h>
#include <sys/resource.h>
#include <sys/stat.h>

KmerPartitioner::KmerPartitioner(unsigned int partitionCount) : partitionCount(partitionCount) {
    if (partitionCount == 0 || (partitionCount & (partitionCount - 1)) != 0) {
        Debug(Debug::ERROR) << "Partition count " << partitionCount << " is not a power of two\n";
        EXIT(EXIT_FAILURE);
    }
    if (partitionCount > 65536) {
        // The score being partitioned is 16 bits, so more partitions than that
        // would leave some permanently empty and silently unbalance the run.
        Debug(Debug::ERROR) << "Partition count " << partitionCount
                            << " exceeds the 16-bit k-mer hash space (max 65536)\n";
        EXIT(EXIT_FAILURE);
    }
    mask = partitionCount - 1;
}

namespace {

// Resident bytes the reduce needs per byte of k-mer bucket: 26/24 for the input
// array plus 20/24 for assignGroup's output array is 1.92, and the round-by-round
// candidate edges take the rest of the headroom.
const uint64_t REDUCE_MEMORY_FACTOR = 3;

// Partitions per worker when the caller passes --workers. More than one so the
// reduce's makespan is not set by whichever worker draws the heaviest partition;
// k-mer hashes spread evenly enough that a small factor suffices.
const uint64_t REDUCE_PARTITIONS_PER_WORKER = 2;

// Floor on what one partition should hold before --workers is allowed to make
// another. Only bounds the worker-derived floor; the memory sizing can still
// derive a smaller partition when the limit genuinely requires it.
const uint64_t MIN_USEFUL_PARTITION_BYTES = 64ULL * 1024 * 1024;

// More waves than this means the budget is wrong, not that the plan is clever.
const unsigned int MAX_SENSIBLE_WAVES = 64;

uint64_t divideRoundingUp(uint64_t value, uint64_t divisor) {
    return (value + divisor - 1) / divisor;
}

// How many bucket descriptors one writer may hold open across flushes.
//
// The reason the writers used to open, write and close on *every* flush is real:
// P can reach 65536 against the 8192 soft limit FileUtil::fixRlimitNoFile raises
// to, so one descriptor per bucket runs the process out of them partway through.
// But that made the map issue ~1e9 open+close pairs at 1e12 against a single
// metadata server -- 1e6 work items x 1024 partitions, plus intra-item flushes.
//
// Keeping a bounded set open costs nothing and removes almost all of it: a wave
// only writes P/W partitions, so the common case fits inside the budget entirely
// and every flush becomes a plain fwrite. Buckets past the budget fall back to
// open-write-close, which is exactly today's cost for those alone.
size_t deriveDescriptorBudget(unsigned int activeCount) {
    // Half the soft limit, so the rest of the process -- the sequence database,
    // the coordination files, the work queue -- keeps its share.
    size_t allowed = 256;
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) == 0 && limit.rlim_cur != RLIM_INFINITY
            && limit.rlim_cur > 64) {
        allowed = static_cast<size_t>(limit.rlim_cur) / 2;
    }
    return std::min<size_t>(std::max<size_t>(allowed, 16), std::max<unsigned int>(activeCount, 1));
}

unsigned int roundUpToPowerOfTwo(uint64_t value) {
    // Above 2^31 the shift below wraps to 0 and the loop never terminates, which
    // an extreme --split-memory-limit or --scratch-budget can reach. Nothing here
    // may legitimately exceed the 16-bit hash space anyway, so refuse early with a
    // number the caller can act on.
    if (value > 65536) {
        Debug(Debug::ERROR) << "Derived a partition or wave count of " << value
                            << ", far past the 65536 the 16-bit k-mer hash can address. "
                            << "--split-memory-limit or --scratch-budget is implausibly small.\n";
        EXIT(EXIT_FAILURE);
    }
    unsigned int result = 1;
    while (result < value) {
        result <<= 1;
    }
    return result;
}

}  // namespace

KmerShuffleSizing deriveKmerShuffleSizing(uint64_t sequenceCount, unsigned int kmersPerSequence,
                                          uint64_t scratchBudgetBytes, uint64_t persistentBytes,
                                          uint64_t workerMemoryBytes, unsigned int workerCount) {
    KmerShuffleSizing sizing;
    sizing.totalKmerBytes = sequenceCount * kmersPerSequence * sizeof(KmerRecord);

    if (scratchBudgetBytes == 0) {
        // No budget given: one wave, and let the memory constraint alone pick P.
        sizing.waveCount = 1;
    } else {
        if (persistentBytes >= scratchBudgetBytes) {
            Debug(Debug::ERROR) << "Scratch budget of " << scratchBudgetBytes
                                << " bytes is already exhausted by the sequence database and "
                                << "surviving edges (" << persistentBytes
                                << " bytes). No wave count can fit the k-mer shuffle.\n";
            EXIT(EXIT_FAILURE);
        }
        const uint64_t available = scratchBudgetBytes - persistentBytes;
        // Rounded up to a power of two so the wave count divides the partition
        // count exactly. A wave owns a contiguous slice of partitions, so if the
        // count did not divide P the largest slice would exceed totalKmerBytes /
        // waveCount and peak scratch would quietly exceed the budget -- and with
        // P below the wave count, the last waves would own nothing at all while
        // the first still held everything. At most this doubles the number of
        // extraction passes; overrunning the scratch filesystem kills the run.
        sizing.waveCount = roundUpToPowerOfTwo(
            std::max<uint64_t>(divideRoundingUp(sizing.totalKmerBytes, available), 1));
    }
    // A budget only slightly above the persistent footprint leaves a sliver for
    // the shuffle and derives a wave count in the hundreds -- a ~256x slowdown
    // presented as a normal run. Refuse instead: at this point the budget is
    // wrong, not the plan.
    if (sizing.waveCount > MAX_SENSIBLE_WAVES) {
        Debug(Debug::ERROR) << "A scratch budget of " << scratchBudgetBytes << " bytes leaves only "
                            << (scratchBudgetBytes > persistentBytes
                                    ? scratchBudgetBytes - persistentBytes : 0)
                            << " bytes for " << sizing.totalKmerBytes << " bytes of k-mers, which "
                            << "needs " << sizing.waveCount << " extraction waves. Each wave "
                            << "re-scans every sequence, so this would run roughly "
                            << sizing.waveCount << "x slower than a single pass. Raise "
                            << "--scratch-budget above " << MAX_SENSIBLE_WAVES << " waves' worth.\n";
        EXIT(EXIT_FAILURE);
    }
    sizing.bytesPerWave = divideRoundingUp(sizing.totalKmerBytes, sizing.waveCount);

    // Sized against totalKmerBytes, NOT bytesPerWave.
    //
    // A wave re-extracts every k-mer and keeps only its own contiguous slice of
    // partition space, so a partition that a wave owns receives *all* of its
    // k-mers -- there is no such thing as a partial partition. A partition
    // therefore holds totalKmerBytes / P regardless of the wave count; what waves
    // divide is how many partitions are on disk at once, which is why
    // bytesPerWave is still the right figure for the scratch budget.
    //
    // Deriving P from bytesPerWave made P a factor of waveCount too small and
    // reported a per-partition size the same factor too low, so a reduce worker
    // held waveCount times the memory the limit allowed. Measured at 100M with
    // two waves: reported 1.575 GB per partition against an actual 3.15 GB.
    if (workerMemoryBytes == 0) {
        sizing.partitionCount = 1;
    } else {
        // The reduce does not hold a partition at its on-disk size. It builds a
        // KmerPosition<T,true,true> array (26 B per 24 B record) and, alongside
        // it, the KmerPosition<T,false,true> array assignGroup writes into (20 B),
        // then accumulates candidate edges across all rounds. Sizing P against the
        // raw bucket bytes therefore overshot resident memory by ~1.9x before the
        // edges were counted at all: with the workflow's default
        // --split-memory-limit 0 on a 2 TB node at 1e11 that derived P = 32, a
        // 1.58 TB partition and ~3 TB of arrays.
        sizing.partitionCount = roundUpToPowerOfTwo(
            divideRoundingUp(sizing.totalKmerBytes * REDUCE_MEMORY_FACTOR, workerMemoryBytes));
    }
    // Both are powers of two, so this makes the wave count a divisor of P.
    sizing.partitionCount = std::max(sizing.partitionCount, sizing.waveCount);
    // A partition is also the reduce's unit of work, so P is a ceiling on how
    // many workers that stage can use -- deriving it from memory alone means a
    // generous --split-memory-limit starves the very allocation it was raised
    // for. Measured at 1e9 with 800G per node: P = 2, so six of eight reduce
    // workers had nothing to claim.
    //
    // Only ever raises P, so the memory sizing above still bounds what one worker
    // holds, and the result stays a power of two, which keeps waveCount a divisor
    // of P. A hint, not a contract: the manifest the caller writes stays
    // authoritative, so a restart with a different --workers cannot re-shape an
    // existing shuffle.
    if (workerCount > 0) {
        // Clamped before roundUpToPowerOfTwo, which EXITs above 65536 with a
        // message blaming --split-memory-limit. An absurd --workers is not that,
        // and the align path clamps the same input rather than dying, so this
        // matches it.
        const uint64_t wanted = std::min<uint64_t>(
            static_cast<uint64_t>(workerCount) * REDUCE_PARTITIONS_PER_WORKER, 65536);
        unsigned int floorP = roundUpToPowerOfTwo(wanted);
        // Never so many partitions that each holds almost nothing: a partition is
        // a directory plus one shard per worker plus a queue item, so a large
        // --workers against a small database -- pass 2 runs over just the
        // representatives -- would otherwise trade real work for millions of tiny
        // files. Measured at 1e7 sequences with --workers 1024: P = 2048, i.e.
        // 2.5 MB per partition.
        while (floorP > 1 && sizing.totalKmerBytes / floorP < MIN_USEFUL_PARTITION_BYTES) {
            floorP /= 2;
        }
        sizing.partitionCount = std::max(sizing.partitionCount, floorP);
    }
    if (sizing.partitionCount > 65536) {
        // Above the 16-bit hash space the partitioner cannot tell partitions
        // apart, so this is a real dead end rather than something to clamp: the
        // run needs a bigger node, more waves, or a smaller input.
        Debug(Debug::ERROR) << "A wave holds " << sizing.bytesPerWave << " bytes of k-mers, which "
                            << "needs more than 65536 partitions to fit " << workerMemoryBytes
                            << " bytes per worker. The 16-bit k-mer hash cannot address that many. "
                            << "Increase the per-worker memory or lower the scratch budget so more "
                            << "waves are used.\n";
        EXIT(EXIT_FAILURE);
    }
    sizing.bytesPerPartition = divideRoundingUp(sizing.totalKmerBytes, sizing.partitionCount);
    return sizing;
}

std::string KmerBucketWriter::partitionDir(const std::string &dir, unsigned int partition) {
    return dir + "/p" + SSTR(partition);
}

namespace KmerBlockCodec {

uint32_t checksum(const uint8_t *data, size_t size) {
    // FNV-1a. Chosen for being a few instructions per byte on data that is
    // already in cache from having just been written or read; this is not a
    // cryptographic check, it is a guard against a block boundary that has drifted.
    uint32_t hash = 2166136261U;
    for (size_t i = 0; i < size; i++) {
        hash ^= data[i];
        hash *= 16777619U;
    }
    return hash;
}

bool headerLooksValid(const KmerBlockHeader &header) {
    if (header.magic != MAGIC) {
        return false;
    }
    if (header.encoding != ENCODING_RAW && header.encoding != ENCODING_PACKED) {
        return false;
    }
    if (header.encoding == ENCODING_PACKED && (header.kmerWidth < 1 || header.kmerWidth > 8)) {
        return false;
    }
    if (header.encoding == ENCODING_RAW &&
        header.payloadBytes != header.recordCount * sizeof(KmerRecord)) {
        return false;
    }
    return true;
}

// Six residues at five bits each, little-endian within a 32-bit word.
static uint32_t packAdjacent(const uint8_t *adjacent) {
    uint32_t bits = 0;
    for (unsigned int i = 0; i < ADJACENT_COUNT; i++) {
        uint32_t value = adjacent[i];
        if (value == UCHAR_MAX) {
            value = ADJACENT_SENTINEL;
        } else if (value >= ADJACENT_SENTINEL) {
            // Only reachable if the alphabet is wider than the packing allows,
            // which is checked once at startup. Reaching it here means that check
            // was skipped, and silently truncating would corrupt the adjacency
            // rounds rather than fail.
            Debug(Debug::ERROR) << "Adjacent residue index " << value
                                << " does not fit the 5-bit packing (max "
                                << (ADJACENT_SENTINEL - 1) << ", or " << UCHAR_MAX
                                << " for absent)\n";
            EXIT(EXIT_FAILURE);
        }
        bits |= value << (i * ADJACENT_BITS);
    }
    return bits;
}

static void unpackAdjacent(uint32_t bits, uint8_t *adjacent) {
    for (unsigned int i = 0; i < ADJACENT_COUNT; i++) {
        const uint32_t value = (bits >> (i * ADJACENT_BITS)) & 0x1FU;
        adjacent[i] = (value == ADJACENT_SENTINEL) ? UCHAR_MAX : static_cast<uint8_t>(value);
    }
}

void encode(std::vector<KmerRecord> &records, bool raw, std::vector<uint8_t> &out) {
    if (records.empty()) {
        return;
    }

    KmerBlockHeader header;
    header.magic = MAGIC;
    header.reserved = 0;
    header.reserved2 = 0;
    header.recordCount = records.size();

    const size_t headerAt = out.size();
    out.resize(headerAt + sizeof(KmerBlockHeader));
    const size_t payloadAt = out.size();

    if (raw) {
        header.encoding = ENCODING_RAW;
        header.kmerWidth = 0;
        header.payloadBytes = records.size() * sizeof(KmerRecord);
        out.resize(payloadAt + header.payloadBytes);
        memcpy(&out[payloadAt], records.data(), header.payloadBytes);
    } else {
        // Ascending ids are what make the id delta a byte; see the note on
        // encode() in the header. Ties broken by (kmer, pos) so the encoding is a
        // pure function of the record multiset, which is what lets a redone work
        // item produce byte-identical output.
        std::sort(records.begin(), records.end(),
                  [](const KmerRecord &a, const KmerRecord &b) {
                      const uint64_t ida = a.getId();
                      const uint64_t idb = b.getId();
                      if (ida != idb) {
                          return ida < idb;
                      }
                      if (a.kmer != b.kmer) {
                          return a.kmer < b.kmer;
                      }
                      return a.pos < b.pos;
                  });

        uint64_t maxKmer = 0;
        for (size_t i = 0; i < records.size(); i++) {
            maxKmer = std::max(maxKmer, records[i].kmer);
        }
        const unsigned int kmerWidth = VarintCodec::fixedWidthFor(maxKmer);
        header.encoding = ENCODING_PACKED;
        header.kmerWidth = static_cast<uint8_t>(kmerWidth);

        // Sized exactly rather than reserved at the worst case. The loose bound is
        // 25 B per record against a ~13.5 B result, and this buffer is live
        // alongside the 24 B record buffer it encodes -- so reserving the bound
        // would put the writer's real footprint at nearly twice what the budget
        // was told to expect. Costing one extra pass over data already in cache is
        // the cheaper half of that trade.
        size_t payloadSize = 0;
        {
            uint64_t prevId = 0;
            for (size_t i = 0; i < records.size(); i++) {
                const uint64_t id = records[i].getId();
                payloadSize += VarintCodec::size(id - prevId);
                payloadSize += VarintCodec::size(records[i].pos);
                prevId = id;
            }
            payloadSize += records.size() * (kmerWidth + ADJACENT_BYTES);
        }
        out.resize(payloadAt + payloadSize);
        uint8_t *cursor = &out[payloadAt];

        // prevId starts at zero, so the first record's "delta" is its absolute id.
        // That keeps the decode loop uniform and costs a few bytes once per block.
        uint64_t prevId = 0;
        for (size_t i = 0; i < records.size(); i++) {
            const uint64_t id = records[i].getId();
            cursor = VarintCodec::write(cursor, id - prevId);
            prevId = id;
            cursor = VarintCodec::writeFixed(cursor, records[i].kmer, kmerWidth);
            cursor = VarintCodec::write(cursor, records[i].pos);
            cursor = VarintCodec::writeFixed(cursor, packAdjacent(records[i].adjacent),
                                             ADJACENT_BYTES);
        }
        header.payloadBytes = static_cast<uint64_t>(cursor - &out[payloadAt]);
        if (header.payloadBytes != payloadSize) {
            Debug(Debug::ERROR) << "K-mer block encoder wrote " << header.payloadBytes
                                << " bytes where it sized " << payloadSize << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    header.checksum = checksum(&out[payloadAt], header.payloadBytes);
    memcpy(&out[headerAt], &header, sizeof(KmerBlockHeader));
}

bool decode(const KmerBlockHeader &header, const uint8_t *payload,
            const LengthRankTable *lengths, std::vector<KmerRecord> &out) {
    if (checksum(payload, header.payloadBytes) != header.checksum) {
        return false;
    }
    if (header.encoding == ENCODING_RAW) {
        if (header.payloadBytes != header.recordCount * sizeof(KmerRecord)) {
            return false;
        }
        const size_t at = out.size();
        out.resize(at + header.recordCount);
        memcpy(&out[at], payload, header.payloadBytes);
        return true;
    }

    // Packed blocks do not carry seqLen; it comes from the key.
    if (lengths == NULL || lengths->isOpen() == false) {
        Debug(Debug::ERROR) << "A packed k-mer block needs the length-rank table to recover "
                            << "sequence lengths, and none was opened\n";
        EXIT(EXIT_FAILURE);
    }

    const uint8_t *cursor = payload;
    const uint8_t *end = payload + header.payloadBytes;
    uint64_t prevId = 0;
    const size_t at = out.size();
    // Grown geometrically, never to exactly the size needed.
    //
    // std::vector::reserve allocates *exactly* what is asked for, so reserving
    // `already + thisBlock` once per block defeats the geometric growth push_back
    // would otherwise get: every block reallocates and copies everything decoded
    // so far, making a shard quadratic in its block count. Measured on 10M, that
    // cost the align stage +65% (111.9 s -> 184.9 s) against the fixed-width path,
    // which uses resize() and is geometric already.
    if (out.capacity() < at + header.recordCount) {
        out.reserve(std::max(at + static_cast<size_t>(header.recordCount), out.capacity() * 2));
    }
    for (uint64_t i = 0; i < header.recordCount; i++) {
        uint64_t idDelta = 0;
        uint64_t kmer = 0;
        uint64_t pos = 0;
        uint64_t adjacentBits = 0;
        if (VarintCodec::read(cursor, end, idDelta) == false ||
            VarintCodec::readFixed(cursor, end, header.kmerWidth, kmer) == false ||
            VarintCodec::read(cursor, end, pos) == false ||
            VarintCodec::readFixed(cursor, end, ADJACENT_BYTES, adjacentBits) == false) {
            out.resize(at);
            return false;
        }
        prevId += idDelta;
        if (prevId > KmerRecord::MAX_ID || pos > UINT16_MAX) {
            out.resize(at);
            return false;
        }
        unsigned int length = 0;
        if (lengths->tryLengthOf(prevId, length) == false) {
            out.resize(at);
            return false;
        }

        KmerRecord record;
        record.kmer = kmer;
        record.setId(prevId);
        record.pos = static_cast<uint16_t>(pos);
        record.seqLen = static_cast<uint16_t>(std::min<unsigned int>(length, UINT16_MAX));
        unpackAdjacent(static_cast<uint32_t>(adjacentBits), record.adjacent);
        out.push_back(record);
    }
    // Trailing bytes mean the block is not what its header says it is.
    if (cursor != end) {
        out.resize(at);
        return false;
    }
    return true;
}

}  // namespace KmerBlockCodec

KmerShardReader::KmerShardReader(const std::string &path) : path(path), torn(false) {
    file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        // A shard that vanished under us is not an error: the reduce unlinks
        // consumed partitions, and a restart can race that sweep.
        torn = false;
    }
}

KmerShardReader::~KmerShardReader() {
    if (file != NULL) {
        fclose(file);
    }
}

bool KmerShardReader::next(std::vector<KmerRecord> &out, const LengthRankTable *lengths) {
    if (file == NULL) {
        return false;
    }
    KmerBlockHeader header;
    const size_t got = fread(&header, 1, sizeof(KmerBlockHeader), file);
    if (got == 0) {
        return false;  // clean end of shard
    }
    if (got != sizeof(KmerBlockHeader) || KmerBlockCodec::headerLooksValid(header) == false) {
        torn = true;
        return false;
    }
    payload.resize(static_cast<size_t>(header.payloadBytes));
    if (header.payloadBytes > 0 &&
        fread(payload.data(), 1, payload.size(), file) != payload.size()) {
        torn = true;
        return false;
    }
    if (KmerBlockCodec::decode(header, payload.data(), lengths, out) == false) {
        torn = true;
        return false;
    }
    return true;
}

uint64_t KmerShardReader::countRecords(const std::string &path, bool &torn) {
    torn = false;
    FILE *file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        return 0;
    }
    uint64_t total = 0;
    while (true) {
        KmerBlockHeader header;
        const size_t got = fread(&header, 1, sizeof(KmerBlockHeader), file);
        if (got == 0) {
            break;
        }
        if (got != sizeof(KmerBlockHeader) || KmerBlockCodec::headerLooksValid(header) == false) {
            torn = true;
            break;
        }
        // Seeking past the payload rather than reading it is what keeps this
        // O(blocks) instead of O(records).
        if (fseeko(file, static_cast<off_t>(header.payloadBytes), SEEK_CUR) != 0) {
            torn = true;
            break;
        }
        total += header.recordCount;
    }
    fclose(file);
    return total;
}

void KmerBucketWriter::createLayout(const std::string &dir, unsigned int partitionCount) {
    if (FileUtil::directoryExists(dir.c_str()) == false) {
        FileUtil::makeDir(dir.c_str());
    }
    for (unsigned int p = 0; p < partitionCount; p++) {
        const std::string path = partitionDir(dir, p);
        if (FileUtil::directoryExists(path.c_str()) == false) {
            // Racing workers may both try this; only a failure that also leaves
            // no directory behind is real.
            if (mkdir(path.c_str(), 0777) != 0 && FileUtil::directoryExists(path.c_str()) == false) {
                Debug(Debug::ERROR) << "Cannot create bucket directory " << path << ": "
                                    << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
}

KmerBucketWriter::KmerBucketWriter(const std::string &dir, unsigned int partitionCount,
                                   const std::string &shardId, size_t bufferBudgetBytes,
                                   unsigned int partitionFrom, unsigned int partitionTo,
                                   bool rawRecords)
    : dir(dir), shardId(shardId), partitionCount(partitionCount),
      mutexes(partitionCount), partitionFrom(partitionFrom),
      partitionTo(partitionTo == 0 ? partitionCount : partitionTo), rawRecords(rawRecords) {
    // At least a handful of records per partition even with a tiny budget, so a
    // large partition count degrades to more frequent flushes rather than to
    // one write syscall per k-mer.
    //
    // Divided by the partitions this wave actually writes, not by P. A wave owns
    // the slice [partitionFrom, partitionTo) and drops every append outside it, so
    // dividing by P gave each live partition a buffer waveCount times smaller than
    // the budget allowed -- and the write size is exactly what a parallel
    // filesystem cares about.
    const unsigned int activeCount =
        (this->partitionTo > partitionFrom) ? (this->partitionTo - partitionFrom) : 1;
    // Divided by both buffers a partition holds: the 24-byte records and the
    // encoded block they are packed into. Sizing against the records alone would
    // understate the writer's footprint by the encoded size, which is the same
    // class of quiet overshoot that reserving the encoder's worst case would be.
    const size_t perPartition =
        bufferBudgetBytes /
        (activeCount * (sizeof(KmerRecord) + KmerBlockCodec::MAX_ENCODED_BYTES_PER_RECORD));
    recordsPerBuffer = std::max<size_t>(perPartition, 16);
    buffers.resize(partitionCount);
    // Reserved up front: push_back until size() >= recordsPerBuffer lets a
    // std::vector's geometric growth land on twice the budgeted capacity, so a
    // 1 GB budget became 2 GB resident across the active partitions.
    for (unsigned int p = partitionFrom; p < this->partitionTo && p < partitionCount; p++) {
        buffers[p].reserve(recordsPerBuffer);
    }
    encodeBuffers.resize(partitionCount);
    files.assign(partitionCount, NULL);
    recordCounts.assign(partitionCount, 0);
    descriptorBudget = deriveDescriptorBudget(activeCount);
    openFiles.store(0);
}

KmerBucketWriter::~KmerBucketWriter() {
    close();
}

// Caller must hold mutexes[partition].
void KmerBucketWriter::closeFile(unsigned int partition) {
    if (files[partition] == NULL) {
        return;
    }
    // Checked: buffered I/O can report ENOSPC, a quota, or a remote filesystem
    // error only at close, and the queue marks the item done straight after
    // flushAll(). An unchecked close here let an item be recorded complete with
    // k-mers missing, which the reduce then reads as "this k-mer had no partner"
    // -- a wrong clustering with no diagnostic.
    FILE *file = files[partition];
    files[partition] = NULL;
    openFiles.fetch_sub(1);
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close k-mer bucket " << partition << " of " << dir << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void KmerBucketWriter::flush(unsigned int partition) {
    std::vector<KmerRecord> &buffer = buffers[partition];
    if (buffer.empty()) {
        return;
    }
    if (files[partition] == NULL) {
        // Opened lazily: with 8192 partitions and a sparse shard, most buckets
        // stay untouched and should not cost a file descriptor or an empty file.
        // Append is safe because this shard belongs to this worker alone.
        const std::string path = partitionDir(dir, partition) + "/" + shardId + ".kmers";
        files[partition] = fopen(path.c_str(), "ab");
        if (files[partition] == NULL) {
            Debug(Debug::ERROR) << "Cannot open bucket " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        openFiles.fetch_add(1);
    }
    // Framed and packed here rather than written as fixed-width structs. The
    // encode buffer is per partition and reused across flushes, so it allocates
    // once; recordsPerBuffer is derived against both buffers so the pair stays
    // inside the budget.
    std::vector<uint8_t> &encoded = encodeBuffers[partition];
    encoded.clear();
    KmerBlockCodec::encode(buffer, rawRecords, encoded);
    if (encoded.empty() == false &&
        fwrite(encoded.data(), 1, encoded.size(), files[partition]) != encoded.size()) {
        // Name the reason: a full scratch filesystem is by far the likeliest way
        // this stage fails, and "cannot write" alone sends you looking for a bug.
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " k-mer records to bucket "
                            << partition << " of " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Held open while the writer's descriptor budget allows, so a partition that
    // is flushed once per work item does not pay an open and a close each time
    // (see deriveDescriptorBudget). Past the budget this is the old
    // open-write-close, which is what makes a 65536-partition run safe.
    if (openFiles.load() > descriptorBudget) {
        closeFile(partition);
    }
    buffer.clear();
}

void KmerBucketWriter::append(unsigned int partition, const KmerRecord &record) {
    // Outside this wave's slice: another wave re-extracts and writes it.
    if (partition < partitionFrom || partition >= partitionTo) {
        return;
    }
    std::lock_guard<std::mutex> guard(mutexes[partition]);
    buffers[partition].push_back(record);
    if (buffers[partition].size() >= recordsPerBuffer) {
        flush(partition);
    }
    // Counted per partition rather than in one shared counter: this runs on every
    // k-mer, and a single counter would be the one point every thread contends on.
    recordCounts[partition]++;
}

uint64_t KmerBucketWriter::getRecordCount() {
    uint64_t total = 0;
    for (unsigned int p = 0; p < partitionCount; p++) {
        std::lock_guard<std::mutex> guard(mutexes[p]);
        total += recordCounts[p];
    }
    return total;
}

void KmerBucketWriter::flushAll() {
    // Descriptors that survive a flush now hold the records in a stdio buffer, so
    // getting them to the OS -- which is exactly what this call owes its caller
    // before the work item is marked done -- takes a checked fflush.
    for (unsigned int p = 0; p < partitionCount; p++) {
        std::lock_guard<std::mutex> guard(mutexes[p]);
        flush(p);
        if (files[p] != NULL && fflush(files[p]) != 0) {
            Debug(Debug::ERROR) << "Cannot flush k-mer bucket for partition " << p << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void KmerBucketWriter::close() {
    for (unsigned int p = 0; p < partitionCount; p++) {
        std::lock_guard<std::mutex> guard(mutexes[p]);
        flush(p);
        closeFile(p);
    }
}

std::vector<std::string> KmerBucketReader::shardFiles(const std::string &dir, unsigned int partition) {
    const std::string path = KmerBucketWriter::partitionDir(dir, partition);
    std::vector<std::string> shards;
    // createLayout pre-creates every partition/bucket directory, so a directory
    // that cannot be opened is a real failure, not a legitimately empty one --
    // treating it as empty let a worker record the item done having read nothing,
    // silently dropping every record in it. An empty partition is an *existing*,
    // readable directory holding no shards.
    DIR *handle = opendir(path.c_str());
    if (handle == NULL) {
        // errno captured before the Debug chain, which does its own I/O and would
        // otherwise overwrite it.
        const int err = errno;
        Debug(Debug::ERROR) << "Cannot open k-mer partition directory " << path << ": " << strerror(err) << "\n";
        EXIT(EXIT_FAILURE);
    }
    struct dirent *entry;
    errno = 0;
    while ((entry = readdir(handle)) != NULL) {
        const std::string name = entry->d_name;
        if (name.size() > 6 && name.compare(name.size() - 6, 6, ".kmers") == 0) {
            shards.push_back(path + "/" + name);
        }
    }
    // errno was set to 0 before the loop precisely so it could be read here, and
    // then never was. A directory read that fails part-way returns NULL exactly as
    // the end of the directory does, so an unchecked loop turns an I/O error into a
    // *partial* shard list -- which downstream reads as "these k-mers had no
    // partner", a smaller clustering with no diagnostic. EdgeBucketWriter::shardFiles
    // already does this; the two are now the same shape.
    const int readErr = errno;
    if (readErr != 0) {
        Debug(Debug::ERROR) << "Cannot read " << path << ": " << strerror(readErr) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (closedir(handle) != 0) {
        const int closeErr = errno;
        Debug(Debug::ERROR) << "Cannot close " << path << ": " << strerror(closeErr) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Sorted so a partition reads back in the same order on every run, which
    // keeps the whole pipeline reproducible regardless of directory order.
    std::sort(shards.begin(), shards.end());
    return shards;
}

uint64_t KmerBucketReader::countRecords(const std::string &dir, unsigned int partition) {
    const std::vector<std::string> shards = shardFiles(dir, partition);
    uint64_t total = 0;
    for (size_t i = 0; i < shards.size(); i++) {
        // Block headers only, seeking past each payload: O(blocks) rather than
        // O(records), and with variable-width records the file size no longer
        // gives the count at all.
        bool torn = false;
        total += KmerShardReader::countRecords(shards[i], torn);
        if (torn) {
            // Counted down to the last whole block rather than refused. A torn
            // tail is exactly what an interrupted worker leaves, and the map
            // redoes that item into a *different* shard, so the records are not
            // lost -- but the torn shard stays on disk, and making it fatal meant
            // every later reduce of that partition died on it forever, with no
            // recovery but deleting the file by hand.
            Debug(Debug::WARNING) << "Bucket " << shards[i]
                                  << " ends on a partial block, as an interrupted worker leaves "
                                     "it; counting the whole blocks it holds.\n";
        }
    }
    return total;
}

void KmerBucketReader::readPartition(const std::string &dir, unsigned int partition,
                                     std::vector<KmerRecord> &out, const LengthRankTable *lengths) {
    const std::vector<std::string> shards = shardFiles(dir, partition);
    // Sized in one go rather than grown per shard: resize() on a vector holding
    // hundreds of gigabytes reallocates and copies the whole array once per shard,
    // and a partition has one shard per worker.
    uint64_t total = 0;
    for (size_t i = 0; i < shards.size(); i++) {
        bool torn = false;
        total += KmerShardReader::countRecords(shards[i], torn);
    }
    out.reserve(out.size() + static_cast<size_t>(total));
    for (size_t i = 0; i < shards.size(); i++) {
        KmerShardReader reader(shards[i]);
        while (reader.next(out, lengths)) {
        }
        if (reader.endedTorn()) {
            Debug(Debug::WARNING) << "Bucket " << shards[i]
                                  << " ends on a partial block, as an interrupted worker leaves "
                                     "it; reading the whole blocks it holds.\n";
        }
    }
}

