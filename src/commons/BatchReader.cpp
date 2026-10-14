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

static const size_t MAX_READ_REQUEST_SIZE = 16 * 1024 * 1024;

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

BatchReader::BatchReader(int maxIoThreads, int fdToDetectStorage) : nextSplitReadRequest(0), readError(0), readErrorFileName(NULL) {
    concurrencyRequests = std::max(1, std::min(storageConcurrency(fdToDetectStorage, maxIoThreads), maxIoThreads));
}

BatchReader::~BatchReader() {
    waitAsyncRead();
}

void BatchReader::readAsync(const std::vector<ReadRequest> &requests) {
    startRead(requests, false);
}

void BatchReader::readSync(const std::vector<ReadRequest> &requests) {
    startRead(requests, true);
    waitAsyncRead();
}

void BatchReader::waitAsyncRead() {
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

