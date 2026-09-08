#include "Lin8DbReader.h"

#include "Debug.h"
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <linux/io_uring.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include "FileUtil.h"
#include "Util.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <zstd.h>

#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)

namespace {
struct Ring {
    int fd;
    unsigned *sqTail, *sqMask, *sqArray;
    unsigned *cqHead, *cqTail, *cqMask;
    struct io_uring_sqe *sqes;
    struct io_uring_cqe *cqes;
    void *sqPtr; size_t sqSize;
    void *cqPtr; size_t cqSize;
    void *sqePtr; size_t sqeSize;
    unsigned entries;
};

bool ringOpen(Ring &r, unsigned entries) {
    struct io_uring_params p;
    memset(&p, 0, sizeof(p));
    const int fd = syscall(__NR_io_uring_setup, entries, &p);
    if (fd < 0) {
        return false;
    }
    const size_t probeSize = sizeof(struct io_uring_probe)
                             + (IORING_OP_READ + 1) * sizeof(struct io_uring_probe_op);
    struct io_uring_probe *probe = (struct io_uring_probe *) calloc(1, probeSize);
    if (probe == NULL) {
        ::close(fd);
        return false;
    }
    const bool readSupported =
        syscall(__NR_io_uring_register, fd, IORING_REGISTER_PROBE, probe, IORING_OP_READ + 1) >= 0
        && probe->ops_len > IORING_OP_READ
        && (probe->ops[IORING_OP_READ].flags & IO_URING_OP_SUPPORTED) != 0;
    free(probe);
    if (readSupported == false) {
        ::close(fd);
        return false;
    }
    size_t sqSize = p.sq_off.array + p.sq_entries * sizeof(unsigned);
    size_t cqSize = p.cq_off.cqes + p.cq_entries * sizeof(struct io_uring_cqe);
    if (p.features & IORING_FEAT_SINGLE_MMAP) {
        sqSize = (cqSize > sqSize) ? cqSize : sqSize;
        cqSize = sqSize;
    }
    void *sq = mmap(NULL, sqSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd,
                    IORING_OFF_SQ_RING);
    if (sq == MAP_FAILED) {
        ::close(fd);
        return false;
    }
    void *cq = sq;
    if ((p.features & IORING_FEAT_SINGLE_MMAP) == 0) {
        cq = mmap(NULL, cqSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd,
                  IORING_OFF_CQ_RING);
        if (cq == MAP_FAILED) {
            munmap(sq, sqSize);
            ::close(fd);
            return false;
        }
    }
    const size_t sqeSize = p.sq_entries * sizeof(struct io_uring_sqe);
    void *sqes = mmap(NULL, sqeSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd,
                      IORING_OFF_SQES);
    if (sqes == MAP_FAILED) {
        if (cq != sq) {
            munmap(cq, cqSize);
        }
        munmap(sq, sqSize);
        ::close(fd);
        return false;
    }
    r.fd = fd;
    r.sqPtr = sq;    r.sqSize = sqSize;
    r.cqPtr = cq;    r.cqSize = cqSize;
    r.sqePtr = sqes; r.sqeSize = sqeSize;
    r.entries = p.sq_entries;
    r.sqTail  = (unsigned *) ((char *) sq + p.sq_off.tail);
    r.sqMask  = (unsigned *) ((char *) sq + p.sq_off.ring_mask);
    r.sqArray = (unsigned *) ((char *) sq + p.sq_off.array);
    r.cqHead  = (unsigned *) ((char *) cq + p.cq_off.head);
    r.cqTail  = (unsigned *) ((char *) cq + p.cq_off.tail);
    r.cqMask  = (unsigned *) ((char *) cq + p.cq_off.ring_mask);
    r.cqes    = (struct io_uring_cqe *) ((char *) cq + p.cq_off.cqes);
    r.sqes    = (struct io_uring_sqe *) sqes;
    return true;
}

void ringClose(Ring &r) {
    munmap(r.sqePtr, r.sqeSize);
    if (r.cqPtr != r.sqPtr) {
        munmap(r.cqPtr, r.cqSize);
    }
    munmap(r.sqPtr, r.sqSize);
    ::close(r.fd);
}
}
#endif

IoRing::IoRing() : ready(false), state(NULL), queued(0), done(0), inflight(0) {}

IoRing::~IoRing() {
#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)
    if (state != NULL) {
        Ring *r = static_cast<Ring *>(state);
        if (ready) {
            ringClose(*r);
        }
        delete r;
    }
#endif
}

bool IoRing::open(unsigned depth) {
#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)
    Ring *r = new Ring();
    ready = ringOpen(*r, depth);
    if (ready == false) {
        delete r;
        r = NULL;
    }
    state = r;
#else
    (void) depth;
#endif
    return ready;
}

void IoRing::preadAll(const char *what) {
    for (size_t i = 0; i < reads.size(); i++) {
        size_t got = 0;
        while (got < reads[i].length) {
            const ssize_t n = pread(reads[i].fd, static_cast<char *>(reads[i].into) + got,
                                    reads[i].length - got, reads[i].offset + got);
            if (n == 0) {
                break;
            }
            if (n < 0) {
                Debug(Debug::ERROR) << "Cannot read " << reads[i].length << " byte from " << what
                                    << " at " << reads[i].offset << "\n";
                EXIT(EXIT_FAILURE);
            }
            got += static_cast<size_t>(n);
        }
        if (got < reads[i].required) {
            Debug(Debug::ERROR) << "Short read of " << got << " byte from " << what << " at "
                                << reads[i].offset << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void IoRing::pump(const char *what, bool untilDone) {
#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)
    Ring &r = *static_cast<Ring *>(state);
    const unsigned sqMask = *r.sqMask;
    const unsigned cqMask = *r.cqMask;
    do {
        unsigned toQueue = 0;
        while (queued < reads.size() && inflight + toQueue < r.entries) {
            const unsigned sqeIdx = (unsigned) ((inflight + toQueue) % r.entries);
            struct io_uring_sqe *sqe = &r.sqes[sqeIdx];
            memset(sqe, 0, sizeof(*sqe));
            sqe->opcode = IORING_OP_READ;
            sqe->fd = reads[queued].fd;
            sqe->off = reads[queued].offset;
            sqe->addr = (uint64_t) (uintptr_t) reads[queued].into;
            sqe->len = static_cast<unsigned>(reads[queued].length);
            sqe->user_data = queued;
            const unsigned tail = *r.sqTail;
            r.sqArray[tail & sqMask] = sqeIdx;
            __atomic_store_n(r.sqTail, tail + 1, __ATOMIC_RELEASE);
            toQueue++;
            queued++;
        }
        unsigned waitFor = 0;
        if (untilDone) {
            waitFor = (queued >= reads.size()) ? (inflight + toQueue) : 1u;
        }
        long ret;
        do {
            ret = syscall(__NR_io_uring_enter, r.fd, toQueue, waitFor,
                          waitFor > 0 ? IORING_ENTER_GETEVENTS : 0u, NULL, 0);
        } while (ret < 0 && errno == EINTR);
        if (ret < 0) {
            Debug(Debug::ERROR) << "Cannot submit reads for " << what << ". Error " << errno << "\n";
            EXIT(EXIT_FAILURE);
        }
        inflight += toQueue;
        unsigned head = *r.cqHead;
        const unsigned cqTail = __atomic_load_n(r.cqTail, __ATOMIC_ACQUIRE);
        while (head != cqTail) {
            struct io_uring_cqe *cqe = &r.cqes[head & cqMask];
            if (cqe->res < 0) {
                Debug(Debug::ERROR) << "Cannot read " << what << ". Error " << -cqe->res << "\n";
                EXIT(EXIT_FAILURE);
            }
            const size_t which = static_cast<size_t>(cqe->user_data);
            if (static_cast<size_t>(cqe->res) < reads[which].required) {
                Debug(Debug::ERROR) << "Short read of " << cqe->res << " byte from " << what
                                    << " at " << reads[which].offset << "\n";
                EXIT(EXIT_FAILURE);
            }
            head++;
            inflight--;
            done++;
        }
        __atomic_store_n(r.cqHead, head, __ATOMIC_RELEASE);
    } while (untilDone && done < reads.size());
#else
    (void) what;
    (void) untilDone;
#endif
}

void IoRing::submit(const char *what) {
    if (inflight != 0) {
        Debug(Debug::ERROR) << "A batch for " << what << " was submitted with " << inflight
                            << " reads of the last one still with the kernel\n";
        EXIT(EXIT_FAILURE);
    }
    queued = 0;
    done = 0;
    inflight = 0;
    if (reads.empty()) {
        return;
    }
#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)
    if (ready == false) {
        preadAll(what);
        done = reads.size();
        return;
    }
    pump(what, false);
#else
    preadAll(what);
    done = reads.size();
#endif
}

void IoRing::await(const char *what) {
    if (done >= reads.size()) {
        return;
    }
#if defined(__linux__) && defined(HAVE_LINUX_IO_URING)
    pump(what, true);
#else
    (void) what;
#endif
}

const uint64_t Lin8DbReader::KEPT_BITMAP_MAGIC = 0x4C494E4356414C44ull;
const uint64_t Lin8DbReader::HEADER_FRAMES_MAGIC = 0x4C494E3848445246ull;

namespace {
size_t fileSizeIfExists(const std::string &path, bool &exists) {
    struct stat sb;
    exists = (stat(path.c_str(), &sb) == 0);
    return exists ? static_cast<size_t>(sb.st_size) : 0;
}

void unmapAndDrop(char *&at, int &fd, size_t length) {
    if (at != NULL) {
        munmap(at, length);
        at = NULL;
    }
    if (fd >= 0) {
#if defined(__linux__)
        posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
#endif
        close(fd);
        fd = -1;
    }
}

char *mapFileReadOnly(const std::string &path, size_t length, int &keptFd) {
    keptFd = -1;
    if (length == 0) {
        return NULL;
    }
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    char *at = static_cast<char *>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0));
    keptFd = fd;
    if (at == MAP_FAILED) {
        close(fd);
        Debug(Debug::ERROR) << "Cannot map " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    Util::madviseLogged(at, length, POSIX_MADV_SEQUENTIAL, path.c_str());
    return at;
}
}

const char *Lin8DbReader::KEPT_BITMAP_SUFFIX = ".clusthash_kept";

Lin8DbReader::Lin8DbReader(const std::string &db, bool withHeaders)
    : db(db), withHeaders(withHeaders), kept(NULL),
      keptMap(NULL), keptSize(0), keptCount(0), keptLoaded(false),
      wantDirect(true), batchAccess(ACCESS_UNHINTED) {}

Lin8DbReader::~Lin8DbReader() {
    close();
}

void Lin8DbReader::open() {
    index.read(db + ".ranges");
    index.checkLengthsDescend();
    for (unsigned int i = 0; i < index.fileCount(); i++) {
        bool exists = false;
        const size_t length = fileSizeIfExists(db + "." + SSTR(i), exists);
        if (exists == false) {
            Debug(Debug::ERROR) << "Data file " << (db + "." + SSTR(i)) << " is missing, the range "
                                << "table declares " << index.fileCount() << " of them\n";
            EXIT(EXIT_FAILURE);
        }
        data.push_back(NULL);
        dataFd.push_back(-1);
        dataSize.push_back(length);
    }
    if (withHeaders) {
        for (unsigned int i = 0; i < data.size(); i++) {
            bool exists = false;
            headerSize.push_back(fileSizeIfExists(db + "_h." + SSTR(i), exists));
            headers.push_back(NULL);
            headerFd.push_back(-1);
            headerFrames.push_back(std::vector<HeaderFrame>());
            headerRawSize.push_back(headerSize.back());
            if (exists == false) {
                Debug(Debug::ERROR) << "Header file " << (db + "_h." + SSTR(i)) << " is missing\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
    const std::string validPath = db + KEPT_BITMAP_SUFFIX;
    {
        const int fd = ::open(validPath.c_str(), O_RDONLY);
        if (fd < 0) {
            return;
        }
        uint64_t header[3];
        const uint64_t words = index.getSize() / 64 + (index.getSize() % 64 != 0);
        struct stat sb;
        if (fstat(fd, &sb) != 0
            || (size_t) sb.st_size < sizeof(header) + words * sizeof(uint64_t)) {
            Debug(Debug::ERROR) << "Bitmap " << validPath << " is truncated\n";
            EXIT(EXIT_FAILURE);
        }
        if (pread(fd, header, sizeof(header), 0) != (ssize_t) sizeof(header)
            || header[0] != KEPT_BITMAP_MAGIC) {
            Debug(Debug::ERROR) << "File " << validPath << " is not a kept bitmap\n";
            EXIT(EXIT_FAILURE);
        }
        if (header[1] != index.getSize()) {
            Debug(Debug::ERROR) << "Bitmap " << validPath << " covers " << header[1]
                                << " sequences, " << db << " holds " << index.getSize() << "\n";
            EXIT(EXIT_FAILURE);
        }
        keptSize = sizeof(header) + words * sizeof(uint64_t);
        void *at = mmap(NULL, keptSize, PROT_READ, MAP_PRIVATE, fd, 0);
        if (at == MAP_FAILED) {
            Debug(Debug::ERROR) << "Cannot map " << validPath << ", " << keptSize << " byte\n";
            EXIT(EXIT_FAILURE);
        }
        ::close(fd);
        keptMap = at;
        kept = reinterpret_cast<const uint64_t *>(static_cast<const char *>(at) + sizeof(header));
        keptCount = words;
        keptLoaded = true;
    }
}

void Lin8DbReader::close() {
    if (keptMap != NULL) {
        munmap(keptMap, keptSize);
        keptMap = NULL;
        kept = NULL;
        keptLoaded = false;
    }
    for (size_t i = 0; i < batch.size(); i++) {
        delete batch[i];
    }
    batch.clear();
    for (size_t i = 0; i < directFd.size(); i++) {
        if (directFd[i] >= 0) {
            ::close(directFd[i]);
            directFd[i] = -1;
        }
    }
    for (size_t i = 0; i < data.size(); i++) {
        unmapAndDrop(data[i], dataFd[i], dataSize[i]);
    }
    for (size_t i = 0; i < headers.size(); i++) {
        unmapAndDrop(headers[i], headerFd[i], headerSize[i]);
    }
    dataFd.clear();
    headerFd.clear();
    data.clear();
    dataSize.clear();
    headers.clear();
    headerSize.clear();
    headerFrames.clear();
    headerRawSize.clear();
}

void Lin8DbReader::mapFile(uint32_t file) const {
    char *at = NULL;
#pragma omp atomic read
    at = data[file];
    if (at != NULL || dataSize[file] == 0) {
        return;
    }
#pragma omp critical(rundb_map)
    {
        if (data[file] == NULL) {
            int fd = -1;
            char *mapped = mapFileReadOnly(db + "." + SSTR(file), dataSize[file], fd);
            dataFd[file] = fd;
#pragma omp atomic write
            data[file] = mapped;
        }
    }
}

void Lin8DbReader::mapHeader(uint32_t file) const {
    if (__atomic_load_n(&headers[file], __ATOMIC_ACQUIRE) != NULL || headerSize[file] == 0) {
        return;
    }
#pragma omp critical(rundb_map)
    {
        if (headers[file] == NULL) {
            int fd = -1;
            char *mapped = mapFileReadOnly(db + "_h." + SSTR(file), headerSize[file], fd);
            headerFd[file] = fd;
            loadHeaderFrames(file, mapped);
            __atomic_store_n(&headers[file], mapped, __ATOMIC_RELEASE);
        }
    }
}

void Lin8DbReader::loadHeaderFrames(uint32_t file, const char *mapped) const {
    const size_t footer = 2 * sizeof(uint64_t);
    if (headerSize[file] < footer) {
        return;
    }
    uint64_t tail[2];
    memcpy(tail, mapped + headerSize[file] - footer, footer);
    if (tail[1] != HEADER_FRAMES_MAGIC) {
        return;
    }
    const uint64_t count = tail[0];
    if (count > (headerSize[file] - footer) / sizeof(HeaderFrame)) {
        Debug(Debug::ERROR) << "Header file " << file << " of " << db << " declares " << count
                            << " frames but is only " << headerSize[file] << " byte long\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t table = count * sizeof(HeaderFrame);
    headerFrames[file].resize(count);
    memcpy(headerFrames[file].data(), mapped + headerSize[file] - footer - table, table);
    for (size_t i = 0; i < count; i++) {
        const HeaderFrame &f = headerFrames[file][i];
        if (f.packedAt + f.packedSize > headerSize[file] - footer - table
            || (i > 0 && f.rawAt != headerFrames[file][i - 1].rawAt + headerFrames[file][i - 1].rawSize)) {
            Debug(Debug::ERROR) << "Header file " << file << " of " << db << " has a broken frame table at " << i << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    headerRawSize[file] = count == 0 ? 0 : headerFrames[file].back().rawAt + headerFrames[file].back().rawSize;
}

const char *Lin8DbReader::fileData(uint32_t file, uint64_t offset) const {
    if (file >= data.size() || offset > dataSize[file]) {
        Debug(Debug::ERROR) << "Range table points at file " << file << " offset " << offset
                            << ", outside the " << data.size() << " data files of " << db << "\n";
        EXIT(EXIT_FAILURE);
    }
    mapFile(file);
    return data[file] + offset;
}

const char *Lin8DbReader::getData(uint64_t rank) const {
    const size_t range = index.rangeIndexOf(rank);
    return fileData(index[range].fileIndex(), index.offsetInRange(range, rank));
}

uint32_t Lin8DbReader::getSeqLen(uint64_t rank, Cursor &cursor) const {
    cursor.at = index.rangeIndexFrom(rank, cursor.at);
    return index[cursor.at].getSeqLen();
}

const char *Lin8DbReader::getData(uint64_t rank, Cursor &cursor) const {
    cursor.at = index.rangeIndexFrom(rank, cursor.at);
    return fileData(index[cursor.at].fileIndex(), index.offsetInRange(cursor.at, rank));
}

bool Lin8DbReader::isKept(uint64_t rank) const {
    if (keptLoaded == false) {
        return true;
    }
    return (kept[rank >> 6] & (uint64_t(1) << (rank & 63))) != 0;
}

void Lin8DbReader::releaseFileSlot(size_t fileSlot) {
    for (size_t file = fileSlot; file < data.size(); file += index.filesPerNode()) {
        unmapAndDrop(data[file], dataFd[file], dataSize[file]);
        if (file < headers.size()) {
            unmapAndDrop(headers[file], headerFd[file], headerSize[file]);
        }
    }
}

uint64_t Lin8DbReader::countKept() const {
    if (keptLoaded == false) {
        return index.getSize();
    }
    uint64_t count = 0;
    for (size_t i = 0; i < keptCount; i++) {
        count += static_cast<uint64_t>(__builtin_popcountll(kept[i]));
    }
    return count;
}

Lin8DbReader::HeaderStream::HeaderStream(const Lin8DbReader &owner)
    : owner(owner), range(0), left(0), at(0), frameFile(0), frame(0) {
    if (owner.headers.empty()) {
        Debug(Debug::ERROR) << "Headers of " << owner.db << " were not opened\n";
        EXIT(EXIT_FAILURE);
    }
}

// the text holding uncompressed offset `at`, straight from the mapping or from the frame that covers it
const char *Lin8DbReader::HeaderStream::frameText(uint32_t file, size_t &avail) {
    const std::vector<HeaderFrame> &frames = owner.headerFrames[file];
    if (frames.empty()) {
        avail = owner.headerSize[file] - at;
        return owner.headers[file] + at;
    }
    if (raw.empty() || file != frameFile || at < frames[frame].rawAt
        || at >= frames[frame].rawAt + frames[frame].rawSize) {
        frame = std::upper_bound(frames.begin(), frames.end(), at,
                                 [](uint64_t a, const HeaderFrame &f) { return a < f.rawAt; }) - frames.begin() - 1;
        frameFile = file;
        const HeaderFrame &f = frames[frame];
        raw.resize(f.rawSize);
        const size_t got = ZSTD_decompress(raw.data(), raw.size(), owner.headers[file] + f.packedAt, f.packedSize);
        if (ZSTD_isError(got) || got != f.rawSize) {
            Debug(Debug::ERROR) << "Header file " << file << " of " << owner.db << " frame " << frame
                                << " does not decompress: " << (ZSTD_isError(got) ? ZSTD_getErrorName(got) : "short") << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    const size_t in = at - frames[frame].rawAt;
    avail = frames[frame].rawSize - in;
    return raw.data() + in;
}

bool Lin8DbReader::HeaderStream::next(const char *&begin, size_t &length) {
    while (left == 0) {
        if (range >= owner.index.rangeCount()) {
            return false;
        }
        left = owner.index.rankAfter(range) - owner.index[range].firstRank();
        at = owner.index[range].headerOffset();
        range++;
    }
    const uint32_t file = owner.index[range - 1].fileIndex();
    if (file >= owner.headers.size()) {
        Debug(Debug::ERROR) << "Length range " << (range - 1) << " of " << owner.db
                            << " points past header file " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    owner.mapHeader(file);
    if (at >= owner.headerRawSize[file]) {
        Debug(Debug::ERROR) << "Length range " << (range - 1) << " of " << owner.db
                            << " points past header file " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t avail = 0;
    const char *from = frameText(file, avail);
    const char *end = static_cast<const char *>(memchr(from, '\n', avail));
    if (end == NULL) {
        Debug(Debug::ERROR) << "Header file " << file << " of " << owner.db
                            << " does not end with a newline\n";
        EXIT(EXIT_FAILURE);
    }
    begin = from;
    length = static_cast<size_t>(end - from);
    at += length + 1;
    left--;
    return true;
}

static const size_t DIRECT_BLOCK = 512;
static const unsigned RING_DEPTH = 1024;

void Lin8DbReader::openBatch(unsigned int threads, size_t arenaBytes,
                            size_t memoryBudget, int revisit, int access) {
    directFd.assign(data.size(), -1);
    batchAccess = access;
    const size_t arenaTotal = (size_t) threads * (arenaBytes + LANES * DIRECT_BLOCK);
    if (arenaTotal >= memoryBudget) {
        Debug(Debug::ERROR) << "Read arenas for " << threads << " thread need "
                            << (arenaTotal >> 20) << " MB, which is the whole "
                            << (memoryBudget >> 20) << " MB limit\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t budget = memoryBudget - arenaTotal;
    const uint64_t sequenceBytes = index.getDataSize();
    wantDirect = revisit == READ_ONCE && sequenceBytes > budget / 2;
    Debug(Debug::INFO) << "Sequence data: " << (sequenceBytes >> 30) << " GB, budget "
                       << (budget >> 30) << " GB after " << (arenaTotal >> 20)
                       << " MB read buffer, reading "
                       << (wantDirect ? "past the OS cache" : "through the OS cache")
                       << (access == ACCESS_RANDOM ? ", advising random access" : access == ACCESS_SEQUENTIAL ? ", advising sequential access" : "") << "\n";
    const size_t laneBytes = arenaBytes / LANES;
    const size_t longest = 2 * ((size_t) index.getMaxSeqLen() + DIRECT_BLOCK);
    if (laneBytes < longest) {
        Debug(Debug::ERROR) << "A read lane of " << laneBytes << " byte cannot hold two sequences of "
                            << index.getMaxSeqLen() << " byte, which needs " << longest << "\n";
        EXIT(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < threads; i++) {
        BatchWorker *worker = new BatchWorker();
        for (unsigned int l = 0; l < LANES; l++) {
            BatchLane &lane = worker->lane[l];
            lane.arena.resize(laneBytes + DIRECT_BLOCK);
            char *at = lane.arena.data();
            const size_t off = reinterpret_cast<uintptr_t>(at) % DIRECT_BLOCK;
            lane.aligned = at + (off == 0 ? 0 : DIRECT_BLOCK - off);
            lane.ring.open(RING_DEPTH);
        }
        batch.push_back(worker);
    }
}

int Lin8DbReader::directOf(uint32_t file) const {
    int fd = -1;
#pragma omp atomic read
    fd = directFd[file];
    if (fd >= 0) {
        return fd;
    }
#pragma omp critical(rundb_direct)
    {
        if (directFd[file] < 0) {
            const std::string path = db + "." + SSTR(file);
            int opened = ::open(path.c_str(), wantDirect ? (O_RDONLY | O_DIRECT) : O_RDONLY);
            if (opened < 0) {
                opened = ::open(path.c_str(), O_RDONLY);
                if (opened < 0) {
                    Debug(Debug::ERROR) << "Cannot open " << path << " for reading\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            if (batchAccess == ACCESS_RANDOM) {
                posix_fadvise(opened, 0, 0, POSIX_FADV_RANDOM);
            } else if (batchAccess == ACCESS_SEQUENTIAL) {
                posix_fadvise(opened, 0, 0, POSIX_FADV_SEQUENTIAL);
            }
#pragma omp atomic write
            directFd[file] = opened;
        }
    }
#pragma omp atomic read
    fd = directFd[file];
    return fd;
}

const char *Lin8DbReader::batchQueryAt(unsigned int thread, unsigned int lane) const {
    return batch[thread]->lane[lane].queryAt;
}

const char *Lin8DbReader::batchAt(unsigned int thread, unsigned int lane, size_t member) const {
    return batch[thread]->lane[lane].memberAt[member];
}

size_t Lin8DbReader::batchRoomFor(uint32_t seqLen) const {
    if (batch.empty()) {
        return 0;
    }
    const size_t room = batch[0]->lane[0].arena.size() - DIRECT_BLOCK;
    return room / ((size_t) seqLen + DIRECT_BLOCK);
}

void Lin8DbReader::awaitBatch(unsigned int thread, unsigned int lane) const {
    batch[thread]->lane[lane].ring.await(db.c_str());
}

bool Lin8DbReader::appendBatchRead(BatchLane &lane, uint64_t rank, Cursor &cursor,
                                  const char *&at) const {
    cursor.at = index.rangeIndexFrom(rank, cursor.at);
    const uint64_t offset = index.offsetInRange(cursor.at, rank);
    const size_t length = index[cursor.at].getSeqLen();
    const int fd = directOf(index[cursor.at].fileIndex());
    const uint64_t blockFrom = offset - offset % DIRECT_BLOCK;
    const uint64_t blockUntil = ((offset + length + DIRECT_BLOCK - 1) / DIRECT_BLOCK) * DIRECT_BLOCK;

    const size_t room = lane.arena.size() - DIRECT_BLOCK;
    size_t used = 0;
    std::vector<IoRing::Read> &reads = lane.ring.list();
    if (reads.empty() == false) {
        IoRing::Read &last = reads.back();
        char *into = static_cast<char *>(last.into);
        used = (size_t) (into - lane.aligned) + last.length;
        const uint64_t end = last.offset + last.length;
        if (last.fd == fd && last.offset <= blockFrom && blockFrom <= end) {
            const size_t extra = (blockUntil > end) ? (size_t) (blockUntil - end) : 0;
            if (used + extra > room) {
                return false;
            }
            last.length += extra;
            last.required = (size_t) (offset + length - last.offset);
            at = into + (size_t) (offset - last.offset);
            return true;
        }
    }
    const size_t span = (size_t) (blockUntil - blockFrom);
    if (used + span > room) {
        return false;
    }
    IoRing::Read read;
    read.into = lane.aligned + used;
    read.fd = fd;
    read.offset = blockFrom;
    read.length = span;
    read.required = (size_t) (offset + length - blockFrom);
    reads.push_back(read);
    at = lane.aligned + used + (size_t) (offset - blockFrom);
    return true;
}

size_t Lin8DbReader::startBatch(uint64_t queryRank, const uint64_t *members, size_t n,
                              unsigned int thread, unsigned int lane) const {
    BatchLane &at = batch[thread]->lane[lane];
    at.memberAt.clear();
    at.ring.list().clear();
    Cursor cursor;
    appendBatchRead(at, queryRank, cursor, at.queryAt);
    size_t loaded = 0;
    for (; loaded < n; loaded++) {
        const char *where = NULL;
        if (appendBatchRead(at, members[loaded], cursor, where) == false) {
            break;
        }
        at.memberAt.push_back(where);
    }
    if (loaded == 0 && n > 0) {
        Debug(Debug::ERROR) << "A read lane of " << at.arena.size() << " byte took none of "
                            << n << " member, so the pass would not move\n";
        EXIT(EXIT_FAILURE);
    }
    at.ring.submit(db.c_str());
    return loaded;
}

