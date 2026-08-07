#include "CandidateEdge.h"

#include "VarintCodec.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <unistd.h>

#include <cerrno>
#include <algorithm>
#include <cstring>

#include <dirent.h>
#include <sys/resource.h>
#include <sys/stat.h>

namespace {
// See deriveDescriptorBudget in KmerPartition.cpp: bucket descriptors are kept
// open across flushes up to a bounded set, so the common case (a bucket count
// well inside the descriptor limit) costs one open per bucket per stage rather
// than one per flush, while a 65536-bucket run still cannot exhaust descriptors.
size_t deriveEdgeDescriptorBudget(unsigned int bucketCount) {
    size_t allowed = 256;
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) == 0 && limit.rlim_cur != RLIM_INFINITY
            && limit.rlim_cur > 64) {
        allowed = static_cast<size_t>(limit.rlim_cur) / 2;
    }
    return std::min<size_t>(std::max<size_t>(allowed, 16), std::max<unsigned int>(bucketCount, 1));
}

// What one edge bucket should hold when memory would allow a bigger one.
//
// Deliberately far below any plausible per-node budget, because this count is
// what bounds align-stage parallelism. 1 GiB puts 1e8 sequences at 32 buckets and
// 1e9 at 128, which feeds any allocation those scales are run on, while staying
// large enough that the per-bucket write buffers EdgeBucketWriter carves out of a
// fixed budget are still worth having.
const uint64_t TARGET_ALIGN_BUCKET_BYTES = 1ULL << 30;

// Buckets per worker when --workers is given. More than one so a worker that
// draws a heavy bucket is not the whole stage's critical path: bucket boundaries
// are uniform in key space, not in edge count, and cluster sizes are heavy-tailed.
const unsigned int ALIGN_BUCKETS_PER_WORKER = 4;
}  // namespace

unsigned int deriveAlignBucketCount(uint64_t dataSize, uint64_t entryCount,
                                    uint64_t memoryLimitBytes, unsigned int workerCount) {
    const uint64_t targetBucketBytes =
        std::min<uint64_t>(memoryLimitBytes / 4, TARGET_ALIGN_BUCKET_BYTES);
    unsigned int bucketCount = 1;
    // A target of zero would double all the way to the ceiling on any input at
    // all, so an unusable memory limit falls back to the byte target alone.
    const uint64_t effectiveTarget = targetBucketBytes > 0 ? targetBucketBytes
                                                           : TARGET_ALIGN_BUCKET_BYTES;
    while (bucketCount < MAX_ALIGN_BUCKETS && dataSize / bucketCount > effectiveTarget) {
        bucketCount *= 2;
    }
    // Only ever raises the count, so the memory ceiling above still holds. A hint
    // rather than a contract: workers may join late, die or be restarted, and the
    // recorded manifest stays authoritative. It matters where the data is too
    // small for the byte target to make enough buckets by itself -- at 1e7 the
    // edges are ~2 GB, which is two buckets however many workers are waiting.
    if (workerCount > 0) {
        const uint64_t wanted = static_cast<uint64_t>(workerCount) * ALIGN_BUCKETS_PER_WORKER;
        while (bucketCount < MAX_ALIGN_BUCKETS && bucketCount < wanted) {
            bucketCount *= 2;
        }
    }
    // More buckets than there are keys leaves empty buckets and a bucketSpan of
    // zero.
    while (bucketCount > 1 && static_cast<uint64_t>(bucketCount) > entryCount) {
        bucketCount /= 2;
    }
    return bucketCount;
}

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

// C++14 has no inline variables, so the in-class initialiser is not a definition.
const uint64_t CandidateEdge::MAX_KEY;

namespace EdgeBlockCodec {

uint32_t checksum(const uint8_t *data, size_t size) {
    uint32_t hash = 2166136261U;
    for (size_t i = 0; i < size; i++) {
        hash ^= data[i];
        hash *= 16777619U;
    }
    return hash;
}

bool headerLooksValid(const EdgeBlockHeader &header) {
    if (header.magic != EdgeBlockHeader::MAGIC) {
        return false;
    }
    if (header.encoding != ENCODING_RAW && header.encoding != ENCODING_PACKED) {
        return false;
    }
    if (header.encoding == ENCODING_RAW &&
        header.payloadBytes != static_cast<uint64_t>(header.recordCount) * sizeof(CandidateEdge)) {
        return false;
    }
    return true;
}

void encode(const std::vector<CandidateEdge> &edges, bool raw, uint32_t partition, uint32_t worker,
            std::vector<uint8_t> &out) {
    if (edges.empty()) {
        return;
    }

    EdgeBlockHeader header;
    header.magic = EdgeBlockHeader::MAGIC;
    header.partition = partition;
    header.worker = worker;
    header.recordCount = static_cast<uint32_t>(edges.size());
    header.reserved[0] = header.reserved[1] = header.reserved[2] = 0;
    header.reserved2 = 0;

    const size_t headerAt = out.size();
    out.resize(headerAt + sizeof(EdgeBlockHeader));
    const size_t payloadAt = out.size();

    if (raw) {
        header.encoding = ENCODING_RAW;
        header.payloadBytes = static_cast<uint64_t>(edges.size()) * sizeof(CandidateEdge);
        out.resize(payloadAt + header.payloadBytes);
        memcpy(&out[payloadAt], edges.data(), static_cast<size_t>(header.payloadBytes));
    } else {
        header.encoding = ENCODING_PACKED;

        // Checked rather than assumed. The reduce sorts by (rep, member,
        // diagonal, strand) before appending and a bucket is a rep range, so a
        // block is sorted by construction -- but a rep that went backwards would
        // wrap its delta into a huge unsigned value and decode as a different
        // edge, silently.
        for (size_t i = 1; i < edges.size(); i++) {
            if (edges[i].getRep() < edges[i - 1].getRep()) {
                Debug(Debug::ERROR) << "Edge block for partition " << partition
                                    << " is not sorted by representative at record " << i
                                    << ". The packed encoding requires the order the reduce "
                                       "produces.\n";
                EXIT(EXIT_FAILURE);
            }
        }

        // Sized exactly, then written. Reserving the worst case instead would put
        // this buffer at 26 B per edge alongside the 17 B edge buffer it encodes,
        // which is nearly twice what the writer's budget was told to expect.
        size_t payloadSize = 0;
        uint64_t prevRep = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            const uint64_t rep = edges[i].getRep();
            const uint64_t member = edges[i].getMember();
            payloadSize += VarintCodec::size(rep - prevRep);
            payloadSize += VarintCodec::size(VarintCodec::zigzag(
                static_cast<int64_t>(member) - static_cast<int64_t>(rep)));
            payloadSize += VarintCodec::size(VarintCodec::zigzag(edges[i].diagonal));
            payloadSize += VarintCodec::size((static_cast<uint64_t>(edges[i].score) << 1) |
                                             (edges[i].reverseStrand ? 1U : 0U));
            prevRep = rep;
        }
        out.resize(payloadAt + payloadSize);
        uint8_t *cursor = &out[payloadAt];

        prevRep = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            const uint64_t rep = edges[i].getRep();
            const uint64_t member = edges[i].getMember();
            cursor = VarintCodec::write(cursor, rep - prevRep);
            cursor = VarintCodec::write(
                cursor, VarintCodec::zigzag(static_cast<int64_t>(member) -
                                            static_cast<int64_t>(rep)));
            cursor = VarintCodec::write(cursor, VarintCodec::zigzag(edges[i].diagonal));
            cursor = VarintCodec::write(cursor, (static_cast<uint64_t>(edges[i].score) << 1) |
                                                    (edges[i].reverseStrand ? 1U : 0U));
            prevRep = rep;
        }
        header.payloadBytes = static_cast<uint64_t>(cursor - &out[payloadAt]);
        if (header.payloadBytes != payloadSize) {
            Debug(Debug::ERROR) << "Edge block encoder wrote " << header.payloadBytes
                                << " bytes where it sized " << payloadSize << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    header.checksum = checksum(&out[payloadAt], static_cast<size_t>(header.payloadBytes));
    memcpy(&out[headerAt], &header, sizeof(EdgeBlockHeader));
}

bool decode(const EdgeBlockHeader &header, const uint8_t *payload,
            std::vector<CandidateEdge> &out) {
    if (checksum(payload, static_cast<size_t>(header.payloadBytes)) != header.checksum) {
        return false;
    }
    const size_t at = out.size();
    if (header.encoding == ENCODING_RAW) {
        if (header.payloadBytes != static_cast<uint64_t>(header.recordCount) * sizeof(CandidateEdge)) {
            return false;
        }
        out.resize(at + header.recordCount);
        memcpy(&out[at], payload, static_cast<size_t>(header.payloadBytes));
        return true;
    }

    const uint8_t *cursor = payload;
    const uint8_t *end = payload + header.payloadBytes;
    uint64_t prevRep = 0;
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
    for (uint32_t i = 0; i < header.recordCount; i++) {
        uint64_t repDelta = 0;
        uint64_t memberZigzag = 0;
        uint64_t diagonalZigzag = 0;
        uint64_t scoreStrand = 0;
        if (VarintCodec::read(cursor, end, repDelta) == false ||
            VarintCodec::read(cursor, end, memberZigzag) == false ||
            VarintCodec::read(cursor, end, diagonalZigzag) == false ||
            VarintCodec::read(cursor, end, scoreStrand) == false) {
            out.resize(at);
            return false;
        }
        prevRep += repDelta;
        const int64_t memberOffset = VarintCodec::unzigzag(memberZigzag);
        const int64_t member = static_cast<int64_t>(prevRep) + memberOffset;
        const int64_t diagonal = VarintCodec::unzigzag(diagonalZigzag);
        if (prevRep > CandidateEdge::MAX_KEY || member < 0 ||
            static_cast<uint64_t>(member) > CandidateEdge::MAX_KEY ||
            diagonal < INT16_MIN || diagonal > INT16_MAX || (scoreStrand >> 1) > UINT16_MAX) {
            out.resize(at);
            return false;
        }

        CandidateEdge edge;
        edge.setRep(prevRep);
        edge.setMember(static_cast<uint64_t>(member));
        edge.diagonal = static_cast<int16_t>(diagonal);
        edge.reverseStrand = static_cast<uint8_t>(scoreStrand & 1U);
        edge.score = static_cast<uint16_t>(scoreStrand >> 1);
        out.push_back(edge);
    }
    if (cursor != end) {
        out.resize(at);
        return false;
    }
    return true;
}

}  // namespace EdgeBlockCodec

EdgeBucketWriter::EdgeBucketWriter(const std::string &dir, unsigned int bucketCount,
                                   const std::string &shardId, size_t bufferBudgetBytes,
                                   bool rawRecords)
    : dir(dir), shardId(shardId), bucketCount(bucketCount), edgeCount(0), closed(false),
      currentPartition(0), currentWorker(-1), rawRecords(rawRecords) {
    // Divided by both buffers a bucket holds: the fixed-width edges and the
    // encoded block they are packed into.
    const size_t perBucket =
        bufferBudgetBytes /
        (bucketCount * (sizeof(CandidateEdge) + EdgeBlockCodec::MAX_ENCODED_BYTES_PER_EDGE));
    edgesPerBuffer = std::max<size_t>(perBucket, 64);
    buffers.resize(bucketCount);
    // Reserved up front, so a std::vector's geometric growth cannot leave the
    // buffers holding twice the budget once they settle at edgesPerBuffer.
    for (unsigned int b = 0; b < bucketCount; b++) {
        buffers[b].reserve(edgesPerBuffer);
    }
    encodeBuffers.resize(bucketCount);
    files.assign(bucketCount, NULL);
    descriptorBudget = deriveEdgeDescriptorBudget(bucketCount);
    openFiles = 0;
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
        files[bucket] = fopen(path.c_str(), "ab");
        if (files[bucket] == NULL) {
            Debug(Debug::ERROR) << "Cannot open edge bucket " << path << ": " << strerror(errno)
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        openFiles++;
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
    std::vector<uint8_t> &encoded = encodeBuffers[bucket];
    encoded.clear();
    EdgeBlockCodec::encode(buffer, rawRecords, currentPartition,
                           static_cast<uint32_t>(currentWorker), encoded);
    if (encoded.empty() == false &&
        fwrite(encoded.data(), 1, encoded.size(), files[bucket]) != encoded.size()) {
        Debug(Debug::ERROR) << "Cannot write " << buffer.size() << " edges to bucket " << bucket
                            << " of " << dir << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Held open while the descriptor budget allows, and closed past it. With up to
    // 65536 buckets against the 8192 descriptors fixRlimitNoFile raises to, keeping
    // one open per bucket unconditionally would run the process out of them
    // partway through the reduce; closing on every flush instead made the stage
    // issue an open and a close per buffer, which is what a metadata server
    // notices at 1e12.
    if (openFiles > descriptorBudget) {
        closeFile(bucket);
    }
    buffer.clear();
}

void EdgeBucketWriter::closeFile(unsigned int bucket) {
    if (files[bucket] == NULL) {
        return;
    }
    FILE *file = files[bucket];
    files[bucket] = NULL;
    openFiles--;
    // Checked: a buffered or remote filesystem reports ENOSPC and quota failures
    // here, and the queue marks the partition done right after flushAll().
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close edge bucket " << bucket << " of " << dir << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
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
        // Descriptors that survive a flush hold the block in a stdio buffer, so
        // getting it to the OS -- what this call owes its caller before the item
        // is marked done -- now takes an explicit, checked fflush.
        if (files[b] != NULL && fflush(files[b]) != 0) {
            Debug(Debug::ERROR) << "Cannot flush edge bucket " << b << " of " << dir << ": "
                                << strerror(errno) << "\n";
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
        closeFile(b);
    }
}

std::vector<std::string> EdgeBucketWriter::shardFiles(const std::string &dir, unsigned int bucket) {
    const std::string path = bucketDir(dir, bucket);
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
        Debug(Debug::ERROR) << "Cannot open edge bucket directory " << path << ": " << strerror(err) << "\n";
        EXIT(EXIT_FAILURE);
    }
    struct dirent *entry;
    errno = 0;
    while ((entry = readdir(handle)) != NULL) {
        const std::string name = entry->d_name;
        if (name.size() > 6 && name.compare(name.size() - 6, 6, ".edges") == 0) {
            shards.push_back(path + "/" + name);
        }
    }
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
    std::vector<uint8_t> payload;
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
        // A worker killed mid-write leaves a partial block, always at the end of
        // its own shard: it appends, and a restarted worker takes a new id and so
        // a new shard. Stopping here discards exactly that tail. The partition it
        // belonged to is redone by another worker, whose copy is complete.
        if (EdgeBlockCodec::headerLooksValid(header) == false ||
            offset + header.payloadBytes > bytes) {
            Debug(Debug::WARNING) << "Edge shard " << path << " ends in a partial block at byte "
                                  << offset << "; it was written by an interrupted worker and the "
                                  << "partition it held was redone.\n";
            break;
        }
        // Out of range is an error, not "keep it". With a non-empty authority the
        // vector covers every partition of the run, so a block naming a partition
        // past its end cannot be attributed to anyone -- it is a block from a
        // differently-partitioned run, or from a run whose queues are not all
        // present. Treating it as wanted silently kept stale or superseded blocks
        // and double-counted their edge support.
        if (authority.empty() == false && header.partition >= authority.size()) {
            Debug(Debug::ERROR) << "Edge shard " << path << " holds a block from partition "
                                << header.partition << ", but the reduce covered only "
                                << authority.size() << " partitions. The edge directory holds "
                                << "output from a run with a different partitioning.\n";
            EXIT(EXIT_FAILURE);
        }
        const bool wanted =
            authority.empty() || authority[header.partition] == static_cast<int64_t>(header.worker);
        if (wanted == false) {
            // A superseded copy: this worker did not record the partition done, so
            // another redid it and its edges are the ones that count.
            if (fseeko(file, static_cast<off_t>(header.payloadBytes), SEEK_CUR) != 0) {
                Debug(Debug::ERROR) << "Cannot skip a superseded block in " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            offset += static_cast<size_t>(header.payloadBytes);
            continue;
        }
        if (header.recordCount > 0) {
            payload.resize(static_cast<size_t>(header.payloadBytes));
            if (payload.empty() == false &&
                fread(payload.data(), 1, payload.size(), file) != payload.size()) {
                Debug(Debug::ERROR) << "Cannot read a block of edge shard " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (EdgeBlockCodec::decode(header, payload.data(), out) == false) {
                // Distinguished from the partial-block case above: the bytes are
                // all present but do not decode, which means the shard is not what
                // its headers claim rather than merely truncated.
                Debug(Debug::ERROR) << "Edge shard " << path << " has a block at byte " << offset
                                    << " that does not decode. The shard is corrupt.\n";
                EXIT(EXIT_FAILURE);
            }
            kept += header.recordCount;
        }
        offset += static_cast<size_t>(header.payloadBytes);
    }
    fclose(file);
    return kept;
}
