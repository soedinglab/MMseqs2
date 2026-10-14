#include "BatchReader.h"
#include "Debug.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#ifdef HAVE_IO_URING
#include <sys/mman.h>
#include <sys/syscall.h>
#endif

static const size_t PIECE_SIZE = 16 * 1024 * 1024;
#ifdef HAVE_IO_URING
static const unsigned int RING_DEPTH = 64;
#endif

BatchReader::BatchReader(int threads) : threads(threads) {
#ifdef HAVE_IO_URING
    ringFd = -1;
    io_uring_params params;
    memset(&params, 0, sizeof(params));
    const int fd = syscall(__NR_io_uring_setup, RING_DEPTH, &params);
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
#endif
}

BatchReader::~BatchReader() {
#ifdef HAVE_IO_URING
    if (ringFd >= 0) {
        munmap(sqes, sqesSize);
        munmap(ring, ringSize);
        close(ringFd);
    }
#endif
}

void BatchReader::read(const std::vector<Request> &requests, const char *name) {
    std::vector<Request> pieces;
    for (size_t i = 0; i < requests.size(); i++) {
        for (size_t at = 0; at < requests[i].minLength; at += PIECE_SIZE) {
            Request piece = requests[i];
            piece.data += at;
            piece.offset += at;
            piece.length = std::min(PIECE_SIZE, requests[i].length - at);
            piece.minLength = std::min(piece.length, requests[i].minLength - at);
            pieces.push_back(piece);
        }
    }
#ifdef HAVE_IO_URING
    if (ringFd >= 0) {
        size_t next = 0;
        size_t done = 0;
        unsigned int inflight = 0;
        while (done < pieces.size()) {
            unsigned int tail = *sqTail;
            for (; next < pieces.size() && inflight < RING_DEPTH; next++, inflight++, tail++) {
                io_uring_sqe &sqe = sqes[tail & sqMask];
                memset(&sqe, 0, sizeof(sqe));
                sqe.opcode = IORING_OP_READ;
                sqe.fd = pieces[next].fd;
                sqe.off = pieces[next].offset;
                sqe.addr = reinterpret_cast<unsigned long>(pieces[next].data);
                sqe.len = pieces[next].length;
                sqe.flags = (pieces[next].length == PIECE_SIZE) ? IOSQE_ASYNC : 0;
                sqe.user_data = next;
            }
            __atomic_store_n(sqTail, tail, __ATOMIC_RELEASE);
            const unsigned int toSubmit = tail - __atomic_load_n(sqHead, __ATOMIC_RELAXED);
            if (syscall(__NR_io_uring_enter, ringFd, toSubmit, 1, IORING_ENTER_GETEVENTS, NULL, 0) < 0 && errno != EAGAIN && errno != EINTR) {
                Debug(Debug::ERROR) << "Cannot read " << name << ": " << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
            unsigned int head = *cqHead;
            for (const unsigned int end = __atomic_load_n(cqTail, __ATOMIC_ACQUIRE); head != end; head++, inflight--, done++) {
                const io_uring_cqe &cqe = cqes[head & cqMask];
                Request rest = pieces[cqe.user_data];
                if (cqe.res == -EAGAIN || cqe.res == -EINTR) {
                    pieces.push_back(rest);
                    continue;
                }
                if (cqe.res <= 0) {
                    Debug(Debug::ERROR) << "Cannot read " << name << ": " << (cqe.res < 0 ? strerror(-cqe.res) : "unexpected end of file") << "\n";
                    EXIT(EXIT_FAILURE);
                }
                if (static_cast<size_t>(cqe.res) < rest.minLength) {
                    rest.data += cqe.res;
                    rest.offset += cqe.res;
                    rest.length -= cqe.res;
                    rest.minLength -= cqe.res;
                    pieces.push_back(rest);
                }
            }
            __atomic_store_n(cqHead, head, __ATOMIC_RELEASE);
        }
        return;
    }
#endif
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(pieces.size() > 1)
    for (size_t i = 0; i < pieces.size(); i++) {
        for (size_t got = 0; got < pieces[i].minLength;) {
            const ssize_t count = pread(pieces[i].fd, pieces[i].data + got, pieces[i].length - got, pieces[i].offset + got);
            if (count < 0 && errno == EINTR) {
                continue;
            }
            if (count <= 0) {
                Debug(Debug::ERROR) << "Cannot read " << name << ": " << (count < 0 ? strerror(errno) : "unexpected end of file") << "\n";
                EXIT(EXIT_FAILURE);
            }
            got += count;
        }
    }
}
