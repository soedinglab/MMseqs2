#ifndef MMSEQS_BATCHREADER_H
#define MMSEQS_BATCHREADER_H

#include <cstddef>
#include <thread>
#include <vector>

#ifdef HAVE_IO_URING
struct io_uring_sqe;
struct io_uring_cqe;
#endif

class BatchReader {
public:
    struct ReadRequest {
        int fileToRead;
        const char *fileName;
        char *destination;
        size_t fileOffset;
        size_t bytesToRead;
        size_t minBytesToRead;
    };

    BatchReader(int maxIoThreads, int fdToDetectStorage, bool useIoUring);
    ~BatchReader();

    // with io_uring, the thread calling readAsync must call waitAsyncRead
    void readAsync(const std::vector<ReadRequest> &requests);
    void waitAsyncRead();
    void readSync(const std::vector<ReadRequest> &requests);

private:
    void startRead(const std::vector<ReadRequest> &requests, bool callingThreadReads);
    void readWithPread();

    int concurrencyRequests;
    std::vector<ReadRequest> splitReadRequests;
    size_t nextSplitReadRequest;
    // errno of the first failed read, -1 for an unexpected end of file
    int readError;
    const char *readErrorFileName;
    std::vector<std::thread> preadThreads;
#ifdef HAVE_IO_URING
    void readWithRing(bool waitForAll);

    int ringFd;
    bool ringWorkersLimited;
    unsigned int inflightReadRequests;
    unsigned int inflightReadLimit;
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
