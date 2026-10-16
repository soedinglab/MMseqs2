#include "BatchReader.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#ifdef __linux__
#include <sys/stat.h>
#include <sys/sysmacros.h>
#include <sys/vfs.h>
#endif

#ifdef HAVE_IO_URING
#include <fcntl.h>
#include <linux/io_uring.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#endif

static const size_t MAX_READ_REQUEST_SIZE = 16 * 1024 * 1024;
#ifdef HAVE_IO_URING
static const unsigned int MAX_RING_DEPTH = 4096;
static const unsigned int MIN_RING_DEPTH = 64;
#endif

static int storageConcurrency(int fd, int maxIoThreads) {
#ifdef __linux__
    struct statfs fs;
    if (fstatfs(fd, &fs) != 0) {
        return 16;
    }
    // GPFS, Lustre | NFS | ext2/3/4, xfs, btrfs, f2fs, zfs, tmpfs
    switch (static_cast<unsigned int>(fs.f_type)) {
        case 0x47504653: case 0x0BD00BD0:
            return 16;
        case 0x6969:
            return 4;
        case 0xEF53: case 0x58465342: case 0x9123683E: case 0xF2F52010: case 0x2FC12FC1: case 0x01021994:
            break;
        default:
            return 16;
    }
    struct stat st;
    if (fstat(fd, &st) == 0) {
        const char *queues[] = { "queue", "../queue" };
        for (size_t i = 0; i < 2; i++) {
            char path[64];
            snprintf(path, sizeof(path), "/sys/dev/block/%u:%u/%s/rotational", major(st.st_dev), minor(st.st_dev), queues[i]);
            FILE *file = fopen(path, "r");
            if (file != NULL) {
                const int rotational = fgetc(file);
                fclose(file);
                return (rotational == '1') ? 1 : maxIoThreads;
            }
        }
    }
#endif
    return maxIoThreads;
}

BatchReader::BatchReader(int maxIoThreads, int fdToDetectStorage, bool useIoUring) : nextSplitReadRequest(0), readError(0), readErrorFileName(NULL) {
    concurrencyRequests = std::max(1, std::min(storageConcurrency(fdToDetectStorage, maxIoThreads), maxIoThreads));
#ifdef HAVE_IO_URING
    ringFd = -1;
    ringWorkersLimited = false;
    inflightReadRequests = 0;
    inflightReadLimit = MIN_RING_DEPTH;
    if (useIoUring == false) {
        return;
    }
    io_uring_params params;
    memset(&params, 0, sizeof(params));
    int fd = syscall(__NR_io_uring_setup, MAX_RING_DEPTH, &params);
    // kernels before 5.12 charge the ring to RLIMIT_MEMLOCK, which may only fit MIN_RING_DEPTH
    if (fd < 0) {
        memset(&params, 0, sizeof(params));
        fd = syscall(__NR_io_uring_setup, MIN_RING_DEPTH, &params);
    }
    if (fd < 0) {
        return;
    }
    const unsigned int probeOps = IORING_OP_READ + 1;
    io_uring_probe *probe = static_cast<io_uring_probe *>(calloc(1, sizeof(io_uring_probe) + probeOps * sizeof(io_uring_probe_op)));
    const bool canRead = probe != NULL && syscall(__NR_io_uring_register, fd, IORING_REGISTER_PROBE, probe, probeOps) >= 0
                         && probe->ops_len > IORING_OP_READ && (probe->ops[IORING_OP_READ].flags & IO_URING_OP_SUPPORTED);
    free(probe);
    if (canRead == false || (params.features & IORING_FEAT_SINGLE_MMAP) == 0) {
        close(fd);
        return;
    }
    ringSize = std::max(params.sq_off.array + params.sq_entries * sizeof(unsigned int), params.cq_off.cqes + params.cq_entries * sizeof(io_uring_cqe));
    ring = static_cast<char *>(mmap(NULL, ringSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd, IORING_OFF_SQ_RING));
    if (ring == MAP_FAILED) {
        close(fd);
        return;
    }
    sqesSize = params.sq_entries * sizeof(io_uring_sqe);
    sqes = static_cast<io_uring_sqe *>(mmap(NULL, sqesSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd, IORING_OFF_SQES));
    if (sqes == MAP_FAILED) {
        munmap(ring, ringSize);
        close(fd);
        return;
    }
    sqHead = reinterpret_cast<unsigned int *>(ring + params.sq_off.head);
    sqTail = reinterpret_cast<unsigned int *>(ring + params.sq_off.tail);
    sqMask = *reinterpret_cast<unsigned int *>(ring + params.sq_off.ring_mask);
    cqHead = reinterpret_cast<unsigned int *>(ring + params.cq_off.head);
    cqTail = reinterpret_cast<unsigned int *>(ring + params.cq_off.tail);
    cqMask = *reinterpret_cast<unsigned int *>(ring + params.cq_off.ring_mask);
    cqes = reinterpret_cast<io_uring_cqe *>(ring + params.cq_off.cqes);
    unsigned int *sqArray = reinterpret_cast<unsigned int *>(ring + params.sq_off.array);
    for (unsigned int i = 0; i < params.sq_entries; i++) {
        sqArray[i] = i;
    }
    ringFd = fd;
#else
    (void)useIoUring;
#endif
}

BatchReader::~BatchReader() {
    waitAsyncRead();
#ifdef HAVE_IO_URING
    if (ringFd >= 0) {
        munmap(sqes, sqesSize);
        munmap(ring, ringSize);
        close(ringFd);
    }
#endif
}

void BatchReader::readAsync(const std::vector<ReadRequest> &requests) {
    startRead(requests, false);
}

void BatchReader::readSync(const std::vector<ReadRequest> &requests) {
    startRead(requests, true);
    waitAsyncRead();
}

void BatchReader::waitAsyncRead() {
#ifdef HAVE_IO_URING
    if (ringFd >= 0) {
        readWithRing(true);
    }
#endif
    readWithPread();
    for (size_t i = 0; i < preadThreads.size(); i++) {
        preadThreads[i].join();
    }
    preadThreads.clear();
    if (readError != 0) {
        Debug(Debug::ERROR) << "Cannot read " << readErrorFileName << ": " << (readError > 0 ? strerror(readError) : "unexpected end of file") << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void BatchReader::startRead(const std::vector<ReadRequest> &requests, bool callingThreadReads) {
    waitAsyncRead();
    splitReadRequests.clear();
    nextSplitReadRequest = 0;
    for (size_t i = 0; i < requests.size(); i++) {
        const ReadRequest &request = requests[i];
        const size_t end = request.fileOffset + request.bytesToRead;
        const size_t minEnd = request.fileOffset + request.minBytesToRead;
        for (size_t at = request.fileOffset; at < minEnd;) {
            const size_t next = std::min((at / MAX_READ_REQUEST_SIZE + 1) * MAX_READ_REQUEST_SIZE, end);
            ReadRequest splitReadRequest = request;
            splitReadRequest.destination += at - request.fileOffset;
            splitReadRequest.fileOffset = at;
            splitReadRequest.bytesToRead = next - at;
            splitReadRequest.minBytesToRead = std::min(next, minEnd) - at;
            splitReadRequests.push_back(splitReadRequest);
            at = next;
        }
    }
#ifdef HAVE_IO_URING
    if (ringFd >= 0) {
        if (ringWorkersLimited == false) {
            unsigned int maxWorkers[2] = { static_cast<unsigned int>(concurrencyRequests), static_cast<unsigned int>(concurrencyRequests) };
            // without a worker cap (Linux < 5.15) io-wq can start 4 workers per CPU
            ringWorkersLimited = syscall(__NR_io_uring_register, ringFd, IORING_REGISTER_IOWQ_MAX_WORKERS, maxWorkers, 2) >= 0;
        }
        inflightReadLimit = ringWorkersLimited ? sqMask + 1 : MIN_RING_DEPTH;
        for (size_t i = 0; i < requests.size(); i++) {
            if (callingThreadReads && (fcntl(requests[i].fileToRead, F_GETFL) & O_DIRECT) == 0) {
                inflightReadLimit = MIN_RING_DEPTH;
                break;
            }
        }
        readWithRing(false);
        return;
    }
#endif
    const size_t preadThreadCount = std::min(static_cast<size_t>(concurrencyRequests), splitReadRequests.size());
    for (size_t i = callingThreadReads ? 1 : 0; i < preadThreadCount; i++) {
        preadThreads.emplace_back(&BatchReader::readWithPread, this);
    }
}

void BatchReader::readWithPread() {
    const size_t pageSize = Util::getPageSize();
    while (__atomic_load_n(&readError, __ATOMIC_RELAXED) == 0) {
        const size_t i = __sync_fetch_and_add(&nextSplitReadRequest, 1);
        if (i >= splitReadRequests.size()) {
            return;
        }
        for (size_t bytesRead = 0; bytesRead < splitReadRequests[i].minBytesToRead;) {
            const ssize_t count = pread(splitReadRequests[i].fileToRead, splitReadRequests[i].destination + bytesRead, splitReadRequests[i].bytesToRead - bytesRead, splitReadRequests[i].fileOffset + bytesRead);
            if (count < 0 && errno == EINTR) {
                continue;
            }
            if (count < static_cast<ssize_t>(std::min(pageSize, splitReadRequests[i].minBytesToRead - bytesRead))) {
                if (__sync_bool_compare_and_swap(&readError, 0, count < 0 ? errno : -1)) {
                    readErrorFileName = splitReadRequests[i].fileName;
                }
                return;
            }
            bytesRead += count;
            if (bytesRead < splitReadRequests[i].minBytesToRead) {
                bytesRead -= bytesRead % pageSize;
            }
        }
    }
}

#ifdef HAVE_IO_URING
void BatchReader::readWithRing(bool waitForAll) {
    const size_t pageSize = Util::getPageSize();
    do {
        unsigned int tail = *sqTail;
        for (; readError == 0 && nextSplitReadRequest < splitReadRequests.size() && inflightReadRequests < inflightReadLimit; nextSplitReadRequest++, inflightReadRequests++, tail++) {
            io_uring_sqe &sqe = sqes[tail & sqMask];
            memset(&sqe, 0, sizeof(sqe));
            sqe.opcode = IORING_OP_READ;
            sqe.fd = splitReadRequests[nextSplitReadRequest].fileToRead;
            sqe.off = splitReadRequests[nextSplitReadRequest].fileOffset;
            sqe.addr = reinterpret_cast<unsigned long>(splitReadRequests[nextSplitReadRequest].destination);
            sqe.len = splitReadRequests[nextSplitReadRequest].bytesToRead;
            sqe.flags = (splitReadRequests[nextSplitReadRequest].bytesToRead == MAX_READ_REQUEST_SIZE) ? IOSQE_ASYNC : 0;
            sqe.user_data = nextSplitReadRequest;
        }
        __atomic_store_n(sqTail, tail, __ATOMIC_RELEASE);
        const unsigned int toSubmit = tail - __atomic_load_n(sqHead, __ATOMIC_RELAXED);
        const unsigned int toWait = (waitForAll && inflightReadRequests > 0) ? 1 : 0;
        if (syscall(__NR_io_uring_enter, ringFd, toSubmit, toWait, IORING_ENTER_GETEVENTS, NULL, 0) < 0 && errno != EAGAIN && errno != EINTR) {
            Debug(Debug::ERROR) << "Cannot read with io_uring: " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int head = *cqHead;
        for (const unsigned int end = __atomic_load_n(cqTail, __ATOMIC_ACQUIRE); head != end; head++, inflightReadRequests--) {
            const io_uring_cqe &cqe = cqes[head & cqMask];
            ReadRequest rest = splitReadRequests[cqe.user_data];
            if (cqe.res == -EAGAIN || cqe.res == -EINTR) {
                splitReadRequests.push_back(rest);
                continue;
            }
            if (cqe.res < static_cast<int>(std::min(pageSize, rest.minBytesToRead))) {
                if (readError == 0) {
                    readError = (cqe.res < 0) ? -cqe.res : -1;
                    readErrorFileName = rest.fileName;
                }
                continue;
            }
            if (static_cast<size_t>(cqe.res) < rest.minBytesToRead) {
                const size_t alignedBytesRead = cqe.res - cqe.res % pageSize;
                rest.destination += alignedBytesRead;
                rest.fileOffset += alignedBytesRead;
                rest.bytesToRead -= alignedBytesRead;
                rest.minBytesToRead -= alignedBytesRead;
                splitReadRequests.push_back(rest);
            }
        }
        __atomic_store_n(cqHead, head, __ATOMIC_RELEASE);
    } while (waitForAll && (inflightReadRequests > 0 || (readError == 0 && nextSplitReadRequest < splitReadRequests.size())));
}
#endif
