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

uint64_t divideRoundingUp(uint64_t value, uint64_t divisor) {
    return (value + divisor - 1) / divisor;
}

unsigned int roundUpToPowerOfTwo(uint64_t value) {
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
        sizing.waveCount = static_cast<unsigned int>(
            std::max<uint64_t>(divideRoundingUp(sizing.totalKmerBytes, available), 1));
    }
    sizing.bytesPerWave = divideRoundingUp(sizing.totalKmerBytes, sizing.waveCount);

    if (workerMemoryBytes == 0) {
        sizing.partitionCount = 1;
    } else {
        sizing.partitionCount =
            roundUpToPowerOfTwo(divideRoundingUp(sizing.bytesPerWave, workerMemoryBytes));
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
    sizing.bytesPerPartition = divideRoundingUp(sizing.bytesPerWave, sizing.partitionCount);
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
                                   const std::string &shardId, size_t bufferBudgetBytes)
    : dir(dir), shardId(shardId), partitionCount(partitionCount),
      mutexes(partitionCount) {
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
        const std::string path = partitionDir(dir, partition) + "/" + shardId + ".kmers";
        files[partition] = fopen(path.c_str(), "wb");
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
    buffer.clear();
}

void KmerBucketWriter::append(unsigned int partition, const KmerRecord &record) {
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
    DIR *handle = opendir(path.c_str());
    if (handle == NULL) {
        // A partition no worker wrote to is empty, not an error.
        return shards;
    }
    struct dirent *entry;
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
            Debug(Debug::ERROR) << "Bucket " << shards[i] << " is " << bytes
                                << " bytes, not a whole number of k-mer records. "
                                << "It was probably written by an interrupted worker.\n";
            EXIT(EXIT_FAILURE);
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
