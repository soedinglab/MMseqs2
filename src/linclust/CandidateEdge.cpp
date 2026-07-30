#include "CandidateEdge.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <cerrno>
#include <algorithm>
#include <cstring>

#include <dirent.h>
#include <sys/stat.h>

std::string EdgeWriter::partitionPath(const std::string &dir, unsigned int partition) {
    return dir + "/p" + SSTR(partition) + ".edges";
}

EdgeWriter::EdgeWriter(const std::string &path, size_t bufferRecords)
    : path(path), file(NULL), bufferRecords(bufferRecords), edgeCount(0), closed(false) {
    buffer.reserve(bufferRecords);
}

EdgeWriter::~EdgeWriter() {
    close();
}

void EdgeWriter::flush() {
    if (buffer.empty()) {
        return;
    }
    if (file == NULL) {
        file = fopen(path.c_str(), "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fwrite(buffer.data(), sizeof(CandidateEdge), buffer.size(), file) != buffer.size()) {
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " edges to " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    buffer.clear();
}

void EdgeWriter::append(const CandidateEdge &edge) {
    buffer.push_back(edge);
    if (buffer.size() >= bufferRecords) {
        flush();
    }
    edgeCount++;
}

void EdgeWriter::close() {
    // The destructor calls this too. Without the guard a second call would fall
    // into the empty-file branch below and truncate what the first call wrote.
    if (closed) {
        return;
    }
    closed = true;
    flush();
    if (file != NULL) {
        // A partition that produced no edge still gets an empty file, so a reader
        // can tell "this partition was reduced and had nothing" from "this
        // partition was never reduced".
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        file = NULL;
    } else {
        FILE *empty = fopen(path.c_str(), "wb");
        if (empty == NULL) {
            Debug(Debug::ERROR) << "Cannot create edge file " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(empty);
    }
}

std::string EdgeBucketWriter::bucketDir(const std::string &dir, unsigned int bucket) {
    return dir + "/r" + SSTR(bucket);
}

void EdgeBucketWriter::createLayout(const std::string &dir, unsigned int bucketCount) {
    if (FileUtil::directoryExists(dir.c_str()) == false) {
        FileUtil::makeDir(dir.c_str());
    }
    for (unsigned int b = 0; b < bucketCount; b++) {
        const std::string path = bucketDir(dir, b);
        if (FileUtil::directoryExists(path.c_str()) == false) {
            // Racing workers may both try; only a failure that also leaves no
            // directory behind is real.
            if (mkdir(path.c_str(), 0777) != 0 && FileUtil::directoryExists(path.c_str()) == false) {
                Debug(Debug::ERROR) << "Cannot create edge bucket " << path << ": "
                                    << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
}

EdgeBucketWriter::EdgeBucketWriter(const std::string &dir, unsigned int bucketCount,
                                   const std::string &shardId, size_t bufferBudgetBytes)
    : dir(dir), shardId(shardId), bucketCount(bucketCount), edgeCount(0), closed(false) {
    const size_t perBucket = bufferBudgetBytes / (bucketCount * sizeof(CandidateEdge));
    edgesPerBuffer = std::max<size_t>(perBucket, 64);
    buffers.resize(bucketCount);
    files.assign(bucketCount, NULL);
}

EdgeBucketWriter::~EdgeBucketWriter() {
    close();
}

void EdgeBucketWriter::flush(unsigned int bucket) {
    std::vector<CandidateEdge> &buffer = buffers[bucket];
    if (buffer.empty()) {
        return;
    }
    if (files[bucket] == NULL) {
        // Opened lazily: a worker whose partitions produced nothing for a bucket
        // should not cost a descriptor or an empty file.
        const std::string path = bucketDir(dir, bucket) + "/" + shardId + ".edges";
        files[bucket] = fopen(path.c_str(), "wb");
        if (files[bucket] == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge bucket " << path << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fwrite(buffer.data(), sizeof(CandidateEdge), buffer.size(), files[bucket]) != buffer.size()) {
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " edges to bucket " << bucket
                            << " of " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    buffer.clear();
}

void EdgeBucketWriter::append(unsigned int bucket, const CandidateEdge &edge) {
    buffers[bucket].push_back(edge);
    if (buffers[bucket].size() >= edgesPerBuffer) {
        flush(bucket);
    }
    edgeCount++;
}

void EdgeBucketWriter::flushAll() {
    for (unsigned int b = 0; b < bucketCount; b++) {
        flush(b);
        if (files[b] != NULL && fflush(files[b]) != 0) {
            Debug(Debug::ERROR) << "Cannot flush edge bucket " << b << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void EdgeBucketWriter::close() {
    if (closed) {
        return;
    }
    closed = true;
    for (unsigned int b = 0; b < bucketCount; b++) {
        flush(b);
        if (files[b] != NULL) {
            if (fclose(files[b]) != 0) {
                Debug(Debug::ERROR) << "Cannot close edge bucket " << b << "\n";
                EXIT(EXIT_FAILURE);
            }
            files[b] = NULL;
        }
    }
}

std::vector<std::string> EdgeBucketWriter::shardFiles(const std::string &dir, unsigned int bucket) {
    const std::string path = bucketDir(dir, bucket);
    std::vector<std::string> shards;
    DIR *handle = opendir(path.c_str());
    if (handle == NULL) {
        return shards;  // a bucket nothing was written to is empty, not an error
    }
    struct dirent *entry;
    while ((entry = readdir(handle)) != NULL) {
        const std::string name = entry->d_name;
        if (name.size() > 6 && name.compare(name.size() - 6, 6, ".edges") == 0) {
            shards.push_back(path + "/" + name);
        }
    }
    closedir(handle);
    std::sort(shards.begin(), shards.end());
    return shards;
}
