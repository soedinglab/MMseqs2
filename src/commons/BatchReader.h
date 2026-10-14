#ifndef MMSEQS_BATCHREADER_H
#define MMSEQS_BATCHREADER_H

#include <cstddef>
#include <vector>

#ifdef HAVE_IO_URING
#include <linux/io_uring.h>
#endif

class BatchReader {
public:
    struct Request {
        int fd;
        char *data;
        size_t offset;
        size_t length;
        size_t minLength;
    };

    explicit BatchReader(int threads);
    ~BatchReader();

    void read(const std::vector<Request> &requests, const char *name);

private:
    int threads;
#ifdef HAVE_IO_URING
    int ringFd;
    char *ring;
    size_t ringSize;
    io_uring_sqe *sqes;
    size_t sqesSize;
    unsigned int *sqHead;
    unsigned int *sqTail;
    unsigned int sqMask;
    unsigned int *cqHead;
    unsigned int *cqTail;
    unsigned int cqMask;
    io_uring_cqe *cqes;
#endif
};

#endif
