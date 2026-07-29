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
    unsigned int bits = 0;
    while ((1u << bits) < partitionCount) {
        bits++;
    }
    shift = 16 - bits;
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
    : dir(dir), shardId(shardId), partitionCount(partitionCount), recordCount(0) {
    // At least a handful of records per partition even with a tiny budget, so a
    // large partition count degrades to more frequent flushes rather than to
    // one write syscall per k-mer.
    const size_t perPartition = bufferBudgetBytes / (partitionCount * sizeof(KmerRecord));
    recordsPerBuffer = std::max<size_t>(perPartition, 16);
    buffers.resize(partitionCount);
    files.assign(partitionCount, NULL);
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
        Debug(Debug::ERROR) << "Cannot write k-mer bucket for partition " << partition << "\n";
        EXIT(EXIT_FAILURE);
    }
    buffer.clear();
}

void KmerBucketWriter::append(unsigned int partition, const KmerRecord &record) {
    buffers[partition].push_back(record);
    if (buffers[partition].size() >= recordsPerBuffer) {
        flush(partition);
    }
    recordCount++;
}

void KmerBucketWriter::close() {
    for (unsigned int p = 0; p < partitionCount; p++) {
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
