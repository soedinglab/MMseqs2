#include "KmerPartition.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>

#include <dirent.h>
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

// More waves than this means the budget is wrong, not that the plan is clever.
const unsigned int MAX_SENSIBLE_WAVES = 64;

uint64_t divideRoundingUp(uint64_t value, uint64_t divisor) {
    return (value + divisor - 1) / divisor;
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
                                          uint64_t workerMemoryBytes) {
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
                                   unsigned int partitionFrom, unsigned int partitionTo)
    : dir(dir), shardId(shardId), partitionCount(partitionCount),
      mutexes(partitionCount), partitionFrom(partitionFrom),
      partitionTo(partitionTo == 0 ? partitionCount : partitionTo) {
    // At least a handful of records per partition even with a tiny budget, so a
    // large partition count degrades to more frequent flushes rather than to
    // one write syscall per k-mer.
    const size_t perPartition = bufferBudgetBytes / (partitionCount * sizeof(KmerRecord));
    recordsPerBuffer = std::max<size_t>(perPartition, 16);
    buffers.resize(partitionCount);
    files.assign(partitionCount, NULL);
    recordCounts.assign(partitionCount, 0);
}

KmerBucketWriter::~KmerBucketWriter() {
    close();
}

void KmerBucketWriter::flush(unsigned int partition) {
    std::vector<KmerRecord> &buffer = buffers[partition];
    if (buffer.empty()) {
        return;
    }
    if (files[partition] == NULL) {
        // Opened lazily: with 8192 partitions and a sparse shard, most buckets
        // stay untouched and should not cost a file descriptor or an empty file.
        // Opened in append mode and closed again after the write (see below).
        // Holding one descriptor per partition for the life of the writer needs P
        // of them -- 8192 at the 1e12 sizing -- against a soft limit
        // FileUtil::fixRlimitNoFile only raises to 8192, so the last partitions
        // would fail to open. Append is safe because this shard belongs to this
        // worker alone.
        const std::string path = partitionDir(dir, partition) + "/" + shardId + ".kmers";
        files[partition] = fopen(path.c_str(), "ab");
        if (files[partition] == NULL) {
            Debug(Debug::ERROR) << "Cannot open bucket " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fwrite(buffer.data(), sizeof(KmerRecord), buffer.size(), files[partition]) != buffer.size()) {
        // Name the reason: a full scratch filesystem is by far the likeliest way
        // this stage fails, and "cannot write" alone sends you looking for a bug.
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " k-mer records to bucket "
                            << partition << " of " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Checked, and before the buffer is dropped: buffered I/O can report ENOSPC,
    // a quota, or a remote filesystem error only at close, and the queue marks the
    // item done straight after flushAll(). An unchecked close here let an item be
    // recorded complete with k-mers missing, which the reduce then reads as "this
    // k-mer had no partner" -- a wrong clustering with no diagnostic.
    if (fclose(files[partition]) != 0) {
        Debug(Debug::ERROR) << "Cannot close k-mer bucket " << partition << " of " << dir << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    files[partition] = NULL;
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
        if (files[p] != NULL) {
            if (fclose(files[p]) != 0) {
                Debug(Debug::ERROR) << "Cannot close k-mer bucket for partition " << p << "\n";
                EXIT(EXIT_FAILURE);
            }
            files[p] = NULL;
        }
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
    closedir(handle);
    // Sorted so a partition reads back in the same order on every run, which
    // keeps the whole pipeline reproducible regardless of directory order.
    std::sort(shards.begin(), shards.end());
    return shards;
}

uint64_t KmerBucketReader::countRecords(const std::string &dir, unsigned int partition) {
    const std::vector<std::string> shards = shardFiles(dir, partition);
    uint64_t total = 0;
    for (size_t i = 0; i < shards.size(); i++) {
        const size_t bytes = FileUtil::getFileSize(shards[i]);
        if (bytes % sizeof(KmerRecord) != 0) {
            // Counted down to the last whole record rather than refused. A torn
            // tail is exactly what an interrupted worker leaves, and the map
            // redoes that item into a *different* shard, so the records are not
            // lost -- but the torn shard stays on disk, and making it fatal meant
            // every later reduce of that partition died on it forever, with no
            // recovery but deleting the file by hand. readPartitionAsPositions
            // already stops at the last whole record for the same reason.
            Debug(Debug::WARNING) << "Bucket " << shards[i] << " ends mid-record at " << bytes
                                  << " bytes, as an interrupted worker leaves it; reading the "
                                  << (bytes / sizeof(KmerRecord)) << " whole records it holds.\n";
        }
        total += bytes / sizeof(KmerRecord);
    }
    return total;
}

void KmerBucketReader::readPartition(const std::string &dir, unsigned int partition,
                                     std::vector<KmerRecord> &out) {
    const std::vector<std::string> shards = shardFiles(dir, partition);
    for (size_t i = 0; i < shards.size(); i++) {
        const size_t bytes = FileUtil::getFileSize(shards[i]);
        if (bytes % sizeof(KmerRecord) != 0) {
            Debug(Debug::ERROR) << "Bucket " << shards[i] << " is truncated\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t count = bytes / sizeof(KmerRecord);
        if (count == 0) {
            continue;
        }
        FILE *file = FileUtil::openFileOrDie(shards[i].c_str(), "rb", true);
        const size_t offset = out.size();
        out.resize(offset + count);
        if (fread(out.data() + offset, sizeof(KmerRecord), count, file) != count) {
            Debug(Debug::ERROR) << "Cannot read bucket " << shards[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(file);
    }
}
