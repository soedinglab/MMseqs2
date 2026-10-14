#ifndef MMSEQS_BATCHREADER_H
#define MMSEQS_BATCHREADER_H

#include <cstddef>
#include <thread>
#include <vector>

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

    BatchReader(int maxIoThreads, int fdToDetectStorage);
    ~BatchReader();

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
};

#endif
