#include "CandidateEdge.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <unistd.h>

#include <cerrno>
#include <algorithm>
#include <cstring>

#include <dirent.h>
#include <sys/stat.h>

std::string EdgeWriter::partitionPath(const std::string &dir, unsigned int partition) {
    return dir + "/p" + SSTR(partition) + ".edges";
}

EdgeWriter::EdgeWriter(const std::string &path, int64_t workerId, size_t bufferRecords)
    : path(path), workerId(workerId), file(NULL), bufferRecords(bufferRecords), edgeCount(0),
      closed(false) {
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
        // Written under a per-worker temporary name and renamed on close. The
        // output of this stage is named by bucket alone, so two workers that end
        // up holding the same bucket -- which a lease expiry can still cause --
        // would otherwise interleave into one file, and the result would stay a
        // whole number of records and pass every downstream integrity check.
        // rename(2) is atomic, so the loser simply replaces the winner with an
        // equally complete file.
        //
        // The worker id, not getpid(): pids are unique per node, and this path is
        // on a shared filesystem. Two workers on two nodes drawing the same pid --
        // routine across a homogeneous allocation -- would open the same "wb" path
        // and truncate each other, publishing a file with a hole of zero-filled
        // edges that every downstream size check still accepts. The bucket writers
        // alongside this one already key their shards on the worker id.
        tmpPath = path + ".w" + SSTR(workerId);
        file = fopen(tmpPath.c_str(), "wb");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge file " << tmpPath << ": " << strerror(errno)
                                << "\n";
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
        //
        // Not fsynced. Losing the node mid-stage is out of scope: the pipeline
        // resumes between stages, not within one, so a stage interrupted that way
        // is re-run from its start rather than trusted from its work queue.
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close edge file " << tmpPath << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        file = NULL;
        if (rename(tmpPath.c_str(), path.c_str()) != 0) {
            Debug(Debug::ERROR) << "Cannot rename " << tmpPath << " to " << path << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        const std::string tmp = path + ".w" + SSTR(workerId);
        FILE *empty = fopen(tmp.c_str(), "wb");
        if (empty == NULL) {
            Debug(Debug::ERROR) << "Cannot create edge file " << tmp << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(empty) != 0 || rename(tmp.c_str(), path.c_str()) != 0) {
            Debug(Debug::ERROR) << "Cannot publish empty edge file " << path << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
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
    : dir(dir), shardId(shardId), bucketCount(bucketCount), edgeCount(0), closed(false),
      currentPartition(0), currentWorker(-1) {
    const size_t perBucket = bufferBudgetBytes / (bucketCount * sizeof(CandidateEdge));
    edgesPerBuffer = std::max<size_t>(perBucket, 64);
    buffers.resize(bucketCount);
    files.assign(bucketCount, NULL);
}

EdgeBucketWriter::~EdgeBucketWriter() {
    close();
}

std::string EdgeBucketWriter::shardPath(unsigned int bucket) const {
    return bucketDir(dir, bucket) + "/" + shardId + ".edges";
}

void EdgeBucketWriter::flush(unsigned int bucket) {
    std::vector<CandidateEdge> &buffer = buffers[bucket];
    if (buffer.empty()) {
        return;
    }
    if (files[bucket] == NULL) {
        // Opened lazily: a worker whose partitions produced nothing for a bucket
        // should not cost a descriptor or an empty file.
        const std::string path = shardPath(bucket);
        // Append-and-close per flush, for the same reason as KmerBucketWriter:
        // one descriptor per bucket would need up to 65536 of them.
        files[bucket] = fopen(path.c_str(), "ab");
        if (files[bucket] == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge bucket " << path << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    // Header then records, as one write each. The header names the producer so a
    // crashed worker's superseded copy can be told apart from a second partition
    // that legitimately produced the same edges (see EdgeBlockHeader).
    if (currentWorker < 0) {
        // Refused rather than stamped with a placeholder: worker 0 is a real id,
        // so an unattributed block would be kept or dropped by coincidence.
        Debug(Debug::ERROR) << "Edges were appended to " << dir << " before beginPartition() named "
                            << "the partition and worker producing them.\n";
        EXIT(EXIT_FAILURE);
    }
    EdgeBlockHeader header;
    header.magic = EdgeBlockHeader::MAGIC;
    header.partition = currentPartition;
    header.worker = static_cast<uint32_t>(currentWorker);
    header.recordCount = static_cast<uint32_t>(buffer.size());
    if (fwrite(&header, sizeof(EdgeBlockHeader), 1, files[bucket]) != 1) {
        Debug(Debug::ERROR) << "Cannot write the block header of bucket " << bucket << " of " << dir
                            << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fwrite(buffer.data(), sizeof(CandidateEdge), buffer.size(), files[bucket]) != buffer.size()) {
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " edges to bucket " << bucket
                            << " of " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Closed here, not held until close(): with up to 65536 buckets against the
    // 8192 descriptors fixRlimitNoFile raises to, keeping one open per bucket runs
    // the process out of descriptors partway through the reduce. This is what the
    // comment above the open has always claimed; only the fclose was missing.
    if (fclose(files[bucket]) != 0) {
        Debug(Debug::ERROR) << "Cannot close edge bucket " << bucket << " of " << dir << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    files[bucket] = NULL;
    buffer.clear();
}

void EdgeBucketWriter::beginPartition(unsigned int partition, int64_t worker) {
    if (currentWorker >= 0 && (partition != currentPartition || worker != currentWorker)) {
        // A block must belong to exactly one partition, or the filter cannot
        // decide it. Anything still buffered is the previous one's.
        for (unsigned int b = 0; b < bucketCount; b++) {
            flush(b);
        }
    }
    currentPartition = partition;
    currentWorker = worker;
}

void EdgeBucketWriter::append(unsigned int bucket, const CandidateEdge &edge) {
    buffers[bucket].push_back(edge);
    if (buffers[bucket].size() >= edgesPerBuffer) {
        flush(bucket);
    }
    edgeCount++;
}

// Pushes every buffered edge to the OS before the caller marks the work item done.
//
// Deliberately not fsynced. That would make the data durable against losing the
// node, which the queue's own completion record already is -- but resume here is
// between stages, not within one, so an interrupted stage is re-run from its start
// rather than trusted item by item. Syncing each touched shard once per item cost
// real traffic on a parallel filesystem for a guarantee nothing consumes.
void EdgeBucketWriter::flushAll() {
    for (unsigned int b = 0; b < bucketCount; b++) {
        flush(b);
    }
}

void EdgeBucketWriter::close() {
    if (closed) {
        return;
    }
    closed = true;
    for (unsigned int b = 0; b < bucketCount; b++) {
        flush(b);
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

size_t EdgeBucketReader::readShard(const std::string &path, const std::vector<int64_t> &authority,
                                   std::vector<CandidateEdge> &out) {
    const size_t bytes = FileUtil::getFileSize(path);
    if (bytes == 0) {
        return 0;
    }
    FILE *file = FileUtil::openFileOrDie(path.c_str(), "rb", true);
    size_t offset = 0;
    size_t kept = 0;
    while (offset + sizeof(EdgeBlockHeader) <= bytes) {
        EdgeBlockHeader header;
        if (fread(&header, sizeof(EdgeBlockHeader), 1, file) != 1) {
            // The loop only enters with a whole header's worth of bytes left, so
            // a short read here is an I/O error rather than the end. Treating it
            // as the end would silently drop the rest of the shard's edges, and
            // missing edges are indistinguishable from pairs that never matched.
            Debug(Debug::ERROR) << "Cannot read the block header at byte " << offset << " of edge "
                                << "shard " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        offset += sizeof(EdgeBlockHeader);
        const size_t blockBytes = static_cast<size_t>(header.recordCount) * sizeof(CandidateEdge);
        // A worker killed mid-write leaves a partial block, always at the end of
        // its own shard: it appends, and a restarted worker takes a new id and so
        // a new shard. Stopping here discards exactly that tail. The partition it
        // belonged to is redone by another worker, whose copy is complete.
        if (header.magic != EdgeBlockHeader::MAGIC || offset + blockBytes > bytes) {
            Debug(Debug::WARNING) << "Edge shard " << path << " ends in a partial block at byte "
                                  << offset << "; it was written by an interrupted worker and the "
                                  << "partition it held was redone.\n";
            break;
        }
        const bool wanted =
            authority.empty() || header.partition >= authority.size() ||
            authority[header.partition] == static_cast<int64_t>(header.worker);
        if (wanted == false) {
            // A superseded copy: this worker did not record the partition done, so
            // another redid it and its edges are the ones that count.
            if (fseek(file, static_cast<long>(blockBytes), SEEK_CUR) != 0) {
                Debug(Debug::ERROR) << "Cannot skip a superseded block in " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            offset += blockBytes;
            continue;
        }
        if (header.recordCount > 0) {
            const size_t at = out.size();
            out.resize(at + header.recordCount);
            if (fread(out.data() + at, sizeof(CandidateEdge), header.recordCount, file) !=
                header.recordCount) {
                Debug(Debug::ERROR) << "Cannot read a block of edge shard " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            kept += header.recordCount;
        }
        offset += blockBytes;
    }
    fclose(file);
    return kept;
}
