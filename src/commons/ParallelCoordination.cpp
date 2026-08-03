#include "ParallelCoordination.h"

#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <ctime>

#include <fcntl.h>
#include <unistd.h>

namespace {

// Retries around EINTR, which a long F_SETLKW wait or a large pread can hit when
// the process takes a signal (Slurm sends plenty).
ssize_t preadFully(int fd, void *buffer, size_t size, off_t offset) {
    char *out = static_cast<char *>(buffer);
    size_t done = 0;
    while (done < size) {
        ssize_t got = pread(fd, out + done, size - done, offset + static_cast<off_t>(done));
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        if (got == 0) {
            break;
        }
        done += static_cast<size_t>(got);
    }
    return static_cast<ssize_t>(done);
}

ssize_t pwriteFully(int fd, const void *buffer, size_t size, off_t offset) {
    const char *in = static_cast<const char *>(buffer);
    size_t done = 0;
    while (done < size) {
        ssize_t put = pwrite(fd, in + done, size - done, offset + static_cast<off_t>(done));
        if (put < 0) {
            if (errno == EINTR) {
                continue;
            }
            return -1;
        }
        done += static_cast<size_t>(put);
    }
    return static_cast<ssize_t>(done);
}

int64_t nowSeconds() {
    return static_cast<int64_t>(time(NULL));
}

}  // namespace

FileLock::FileLock(const std::string &path) : path(path), fd(-1) {
    fd = open(path.c_str(), O_RDWR | O_CREAT, 0666);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Could not open coordination file " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

FileLock::~FileLock() {
    if (fd >= 0) {
        close(fd);
    }
}

void FileLock::lock() {
    threadMutex.lock();

    struct flock request;
    memset(&request, 0, sizeof(request));
    request.l_type = F_WRLCK;
    request.l_whence = SEEK_SET;
    request.l_start = 0;
    // len 0 means "to end of file, including anything appended later".
    request.l_len = 0;

    while (fcntl(fd, F_SETLKW, &request) < 0) {
        if (errno == EINTR) {
            continue;
        }
        Debug(Debug::ERROR) << "Could not lock coordination file " << path << ": "
                            << strerror(errno) << "\n";
        // A failure here means the filesystem is not honouring fcntl locks (a
        // Lustre mount without -o flock is the usual cause). Continuing would
        // silently corrupt shared state, so stop instead.
        threadMutex.unlock();
        EXIT(EXIT_FAILURE);
    }
}

void FileLock::unlock() {
    struct flock request;
    memset(&request, 0, sizeof(request));
    request.l_type = F_UNLCK;
    request.l_whence = SEEK_SET;
    request.l_start = 0;
    request.l_len = 0;

    while (fcntl(fd, F_SETLK, &request) < 0) {
        if (errno == EINTR) {
            continue;
        }
        Debug(Debug::ERROR) << "Could not unlock coordination file " << path << ": "
                            << strerror(errno) << "\n";
        break;
    }

    threadMutex.unlock();
}

SharedCounter::SharedCounter(const std::string &path) : lock(path) {
}

int64_t SharedCounter::readLocked() {
    int64_t value = 0;
    ssize_t got = preadFully(lock.getFd(), &value, sizeof(value), 0);
    if (got < 0) {
        Debug(Debug::ERROR) << "Could not read shared counter: " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // A short read means the file was just created and is still empty, which is
    // the same as a counter of zero.
    if (got != static_cast<ssize_t>(sizeof(value))) {
        return 0;
    }
    return value;
}

void SharedCounter::writeLocked(int64_t value) {
    if (pwriteFully(lock.getFd(), &value, sizeof(value), 0) < 0) {
        Debug(Debug::ERROR) << "Could not write shared counter: " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Push the update out before releasing the lock, so a reader on another node
    // that acquires the lock next is guaranteed to observe it.
    fsync(lock.getFd());
}

int64_t SharedCounter::fetchAdd(int64_t n) {
    lock.lock();
    int64_t previous = readLocked();
    writeLocked(previous + n);
    lock.unlock();
    return previous;
}

int64_t SharedCounter::get() {
    lock.lock();
    int64_t value = readLocked();
    lock.unlock();
    return value;
}

void SharedCounter::await(int64_t target, unsigned int pollSeconds) {
    while (get() < target) {
        sleep(pollSeconds);
    }
}

WorkQueue::WorkQueue(const std::string &path, int64_t itemCount)
        : path(path), itemCount(itemCount), lock(path) {
    lock.lock();
    initialiseLocked();
    lock.unlock();
}

void WorkQueue::initialiseLocked() {
    Header header;
    ssize_t got = preadFully(lock.getFd(), &header, sizeof(header), 0);
    if (got < 0) {
        Debug(Debug::ERROR) << "Could not read work queue " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (got == static_cast<ssize_t>(sizeof(header)) && header.magic == MAGIC) {
        if (header.version != VERSION) {
            Debug(Debug::ERROR) << "Work queue " << path << " has version " << header.version
                                << ", this build writes version " << VERSION << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (header.itemCount != static_cast<uint64_t>(itemCount)) {
            // Resuming a run that partitioned its work differently would silently
            // skip or duplicate items, so refuse rather than guess.
            Debug(Debug::ERROR) << "Work queue " << path << " was created for "
                                << header.itemCount << " items but this run has " << itemCount
                                << ". Remove the coordination directory to start over.\n";
            EXIT(EXIT_FAILURE);
        }
        // The records, not the counter, are the truth. completeLocked marks an item
        // DONE and then increments doneCount as two separate writes, so a node lost
        // between them leaves the count one short -- and since nothing is then
        // claimable and allDone() never becomes true, drain() would report the stage
        // stalled on every restart, permanently. Recounting once per open costs a
        // single bulk read and makes that unreachable.
        std::vector<Record> records;
        readRecordsLocked(records);
        uint64_t done = 0;
        for (size_t i = 0; i < records.size(); i++) {
            if (records[i].state == DONE) {
                done++;
            }
        }
        if (done != header.doneCount) {
            Debug(Debug::WARNING) << "Work queue " << path << " recorded " << header.doneCount
                                  << " completed items but holds " << done
                                  << "; repairing the count from the records.\n";
            header.doneCount = done;
            writeHeaderLocked(header);
            fsync(lock.getFd());
        }
        return;
    }

    // Zero-fill the record array so every later access is a plain offset read
    // rather than a short read off the end of a sparse file. PENDING is state 0,
    // so zeroing is also the correct initial state.
    //
    // Records first, header last. The header is what makes the file look like a
    // queue, so writing it first meant a death during the zero-fill left a file
    // that passed the magic and itemCount checks above over a record array that
    // was still short -- every later claim() then read past EOF and exited, for
    // every worker, forever, with no recovery but deleting the directory by hand.
    // The window is one pwrite for the linclust queues but ~100 for
    // createdbparallel at 1e11.
    const size_t batchSize = 65536;
    Record *blank = new Record[batchSize];
    memset(blank, 0, batchSize * sizeof(Record));
    for (int64_t written = 0; written < itemCount; written += static_cast<int64_t>(batchSize)) {
        size_t count = static_cast<size_t>(std::min<int64_t>(batchSize, itemCount - written));
        if (pwriteFully(lock.getFd(), blank, count * sizeof(Record), recordOffset(written)) < 0) {
            Debug(Debug::ERROR) << "Could not initialise work queue " << path << ": "
                                << strerror(errno) << "\n";
            delete[] blank;
            EXIT(EXIT_FAILURE);
        }
    }
    delete[] blank;
    fsync(lock.getFd());

    Header fresh;
    memset(&fresh, 0, sizeof(fresh));
    fresh.magic = MAGIC;
    fresh.version = VERSION;
    fresh.itemCount = static_cast<uint64_t>(itemCount);
    fresh.doneCount = 0;
    fresh.nextHint = 0;
    writeHeaderLocked(fresh);
    fsync(lock.getFd());
}

WorkQueue::Header WorkQueue::readHeaderLocked() {
    Header header;
    if (preadFully(lock.getFd(), &header, sizeof(header), 0) != static_cast<ssize_t>(sizeof(header))) {
        Debug(Debug::ERROR) << "Could not read work queue header " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    return header;
}

void WorkQueue::writeHeaderLocked(const Header &header) {
    if (pwriteFully(lock.getFd(), &header, sizeof(header), 0) < 0) {
        Debug(Debug::ERROR) << "Could not write work queue header " << path << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

WorkQueue::Record WorkQueue::readRecordLocked(int64_t index) {
    Record record;
    if (preadFully(lock.getFd(), &record, sizeof(record), recordOffset(index))
            != static_cast<ssize_t>(sizeof(record))) {
        Debug(Debug::ERROR) << "Could not read work queue record " << index << " in " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    return record;
}

void WorkQueue::readRecordsLocked(std::vector<Record> &out) {
    out.resize(static_cast<size_t>(itemCount));
    if (itemCount == 0) {
        return;
    }
    const size_t bytes = static_cast<size_t>(itemCount) * sizeof(Record);
    if (preadFully(lock.getFd(), &out[0], bytes, recordOffset(0)) != static_cast<ssize_t>(bytes)) {
        Debug(Debug::ERROR) << "Could not read the " << itemCount << " records of work queue "
                            << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void WorkQueue::writeRecordLocked(int64_t index, const Record &record) {
    if (pwriteFully(lock.getFd(), &record, sizeof(record), recordOffset(index)) < 0) {
        Debug(Debug::ERROR) << "Could not write work queue record " << index << " in " << path
                            << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

int64_t WorkQueue::claim(int64_t workerId, int64_t leaseSeconds) {
    lock.lock();

    Header header = readHeaderLocked();
    const int64_t now = nowSeconds();

    std::vector<Record> records;
    readRecordsLocked(records);

    int64_t firstUnfinished = -1;
    for (int64_t index = static_cast<int64_t>(header.nextHint); index < itemCount; index++) {
        Record record = records[static_cast<size_t>(index)];
        if (record.state == DONE) {
            continue;
        }
        if (firstUnfinished < 0) {
            firstUnfinished = index;
        }
        // Untouched, or held by a worker whose lease ran out -- in both cases the
        // item is ours to take. Re-claiming after a lease expiry is what makes a
        // killed job recoverable without an operator deciding what to redo.
        if (record.state == PENDING
                || (record.state == CLAIMED && static_cast<int64_t>(record.leaseExpiry) <= now)) {
            record.state = CLAIMED;
            record.worker = static_cast<uint32_t>(workerId);
            record.leaseExpiry = static_cast<uint64_t>(now + leaseSeconds);
            writeRecordLocked(index, record);

            if (firstUnfinished > static_cast<int64_t>(header.nextHint)) {
                header.nextHint = static_cast<uint64_t>(firstUnfinished);
                writeHeaderLocked(header);
            }
            fsync(lock.getFd());
            lock.unlock();
            return index;
        }
    }

    if (firstUnfinished > static_cast<int64_t>(header.nextHint)) {
        header.nextHint = static_cast<uint64_t>(firstUnfinished);
        writeHeaderLocked(header);
        fsync(lock.getFd());
    }

    lock.unlock();
    return -1;
}

void WorkQueue::completeLocked(int64_t index, int64_t workerId) {
    Record record = readRecordLocked(index);
    if (record.state == DONE) {
        // Already recorded, either by us before a crash or by another worker that
        // re-claimed the item after our lease lapsed. Nothing to do; making this
        // idempotent is what lets a worker redo an item it may or may not have
        // finished before it died.
        return;
    }

    record.state = DONE;
    record.worker = static_cast<uint32_t>(workerId);
    record.leaseExpiry = 0;
    writeRecordLocked(index, record);

    Header header = readHeaderLocked();
    header.doneCount += 1;
    writeHeaderLocked(header);
    fsync(lock.getFd());
}

void WorkQueue::renew(int64_t index, int64_t workerId, int64_t leaseSeconds) {
    lock.lock();
    Record record = readRecordLocked(index);
    // Only extend a claim we still hold. If the lease already lapsed and someone
    // else took the item, stealing it back would run it twice.
    if (record.state == CLAIMED && record.worker == static_cast<uint32_t>(workerId)) {
        record.leaseExpiry = static_cast<uint64_t>(nowSeconds() + leaseSeconds);
        writeRecordLocked(index, record);
        fsync(lock.getFd());
    }
    lock.unlock();
}

void WorkQueue::complete(int64_t index, int64_t workerId) {
    lock.lock();
    completeLocked(index, workerId);
    lock.unlock();
}

void WorkQueue::release(int64_t index, int64_t workerId) {
    lock.lock();
    Record record = readRecordLocked(index);
    if (record.state == CLAIMED && record.worker == static_cast<uint32_t>(workerId)) {
        record.state = PENDING;
        record.worker = 0;
        record.leaseExpiry = 0;
        writeRecordLocked(index, record);
        fsync(lock.getFd());
    }
    lock.unlock();
}

int64_t WorkQueue::getDoneCount() {
    lock.lock();
    Header header = readHeaderLocked();
    lock.unlock();
    return static_cast<int64_t>(header.doneCount);
}

// True while some item is held under a lease that has not yet expired.
//
// drain() uses this to tell "nobody has finished anything for a long time
// because the run is stuck" from "because one item legitimately takes hours".
// Heartbeats keep a live holder's expiry in the future, so a lease that is still
// valid means a worker is still on it.
bool WorkQueue::hasLiveClaim() {
    const int64_t now = nowSeconds();
    lock.lock();
    std::vector<Record> records;
    readRecordsLocked(records);
    bool live = false;
    for (int64_t i = 0; i < itemCount && live == false; i++) {
        const Record &record = records[static_cast<size_t>(i)];
        if (record.state == CLAIMED && static_cast<int64_t>(record.leaseExpiry) > now) {
            live = true;
        }
    }
    lock.unlock();
    return live;
}

bool WorkQueue::allDone() {
    return getDoneCount() >= itemCount;
}

bool WorkQueue::readCompletedWorkers(const std::string &path, std::vector<int64_t> &workers) {
    workers.clear();
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        return false;
    }
    Header header;
    if (preadFully(fd, &header, sizeof(Header), 0) != static_cast<ssize_t>(sizeof(Header))) {
        close(fd);
        Debug(Debug::ERROR) << "Cannot read the header of work queue " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (header.magic != MAGIC || header.version != VERSION) {
        close(fd);
        Debug(Debug::ERROR) << "File " << path << " is not a work queue\n";
        EXIT(EXIT_FAILURE);
    }
    workers.assign(static_cast<size_t>(header.itemCount), -1);
    for (uint64_t i = 0; i < header.itemCount; i++) {
        Record record;
        if (preadFully(fd, &record, sizeof(Record), recordOffset(static_cast<int64_t>(i))) !=
            static_cast<ssize_t>(sizeof(Record))) {
            close(fd);
            Debug(Debug::ERROR) << "Cannot read record " << i << " of work queue " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (record.state == DONE) {
            workers[static_cast<size_t>(i)] = static_cast<int64_t>(record.worker);
        }
    }
    close(fd);
    return true;
}

bool WorkQueue::awaitAll(unsigned int pollSeconds, unsigned int stallSeconds) {
    int64_t lastDone = -1;
    int64_t lastProgress = nowSeconds();
    while (true) {
        int64_t done = getDoneCount();
        if (done >= itemCount) {
            return true;
        }
        if (done != lastDone) {
            lastDone = done;
            lastProgress = nowSeconds();
        } else if (stallSeconds > 0
                   && nowSeconds() - lastProgress > static_cast<int64_t>(stallSeconds)) {
            Debug(Debug::WARNING) << "Work queue " << path << " made no progress for "
                                  << stallSeconds << " s (" << done << "/" << itemCount
                                  << " done)\n";
            return false;
        }
        sleep(pollSeconds);
    }
}
