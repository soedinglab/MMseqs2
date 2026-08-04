#include "ParallelCoordination.h"

#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <mutex>

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

namespace {

// One descriptor and one mutex per path, shared by every FileLock naming it.
//
// Two fcntl properties make this necessary rather than tidy:
//
//   - record locks are owned by the *process*, so two FileLock objects for one
//     path in two threads both believe they hold it. Their own std::mutexes are
//     different objects, so nothing serialises them and the shared state they
//     guard is updated concurrently.
//   - closing *any* descriptor to a file drops that process's locks on it. So a
//     short-lived FileLock going out of scope silently unlocked a long-lived one
//     naming the same path, mid-critical-section.
//
// Neither is hypothetical in this codebase: the stages construct SharedCounter,
// WorkQueue and bare FileLock objects over the same coordination directory, and
// createdbparallel nests a plan lock inside a queue drain.
//
// Keyed on the canonicalised *directory* plus the basename, because the lock file
// may not exist yet -- realpath() on the file itself would resolve differently
// before and after creation and hand out two entries for one file.
//
// The file is opened only when an entry is created, never speculatively. That is
// load-bearing rather than an optimisation: closing a duplicate descriptor for a
// file releases *all* of this process's record locks on it, so opening and then
// closing a second descriptor would silently drop a lock another thread was
// holding at that moment.
class LockRegistry {
public:
    struct Entry {
        int fd;
        std::mutex mutex;
        size_t refs;
        Entry() : fd(-1), refs(0) {}
    };

    static LockRegistry &instance() {
        static LockRegistry registry;
        return registry;
    }

    static std::string keyOf(const std::string &path) {
        const size_t slash = path.find_last_of('/');
        const std::string dir = slash == std::string::npos ? std::string(".")
                                                          : path.substr(0, slash);
        const std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
        char *resolved = realpath(dir.c_str(), NULL);
        if (resolved == NULL) {
            return path;
        }
        const std::string key = std::string(resolved) + "/" + base;
        free(resolved);
        return key;
    }

    Entry *acquire(const std::string &path) {
        const std::string key = keyOf(path);
        std::lock_guard<std::mutex> guard(tableMutex);
        Entry *&slot = table[key];
        if (slot == NULL) {
            const int fd = open(path.c_str(), O_RDWR | O_CREAT, 0666);
            if (fd < 0) {
                Debug(Debug::ERROR) << "Could not open coordination file " << path << ": "
                                    << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
            slot = new Entry();
            slot->fd = fd;
        }
        slot->refs++;
        return slot;
    }

    void release(Entry *entry) {
        std::lock_guard<std::mutex> guard(tableMutex);
        if (entry == NULL || entry->refs == 0) {
            return;
        }
        entry->refs--;
        // The descriptor is deliberately *not* closed at zero. Closing it would
        // drop any lock a later FileLock over the same path re-acquires, and the
        // count of distinct coordination files in a run is a handful.
    }

private:
    std::mutex tableMutex;
    std::map<std::string, Entry *> table;
};

}  // namespace

FileLock::FileLock(const std::string &path) : path(path), fd(-1), entry(NULL) {
    LockRegistry::Entry *shared = LockRegistry::instance().acquire(path);
    entry = shared;
    fd = shared->fd;
}

FileLock::~FileLock() {
    LockRegistry::instance().release(static_cast<LockRegistry::Entry *>(entry));
    // fd belongs to the registry, not to this object; see LockRegistry::release.
}

void FileLock::lock() {
    static_cast<LockRegistry::Entry *>(entry)->mutex.lock();

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
        static_cast<LockRegistry::Entry *>(entry)->mutex.unlock();
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

    static_cast<LockRegistry::Entry *>(entry)->mutex.unlock();
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
        uint64_t done = 0;
        for (int64_t from = 0; from < itemCount; from += SCAN_WINDOW) {
            const int64_t count = std::min<int64_t>(SCAN_WINDOW, itemCount - from);
            readRecordsLocked(from, count, records);
            for (size_t i = 0; i < records.size(); i++) {
                if (records[i].state == DONE) {
                    done++;
                }
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

void WorkQueue::readRecordsLocked(int64_t from, int64_t count, std::vector<Record> &out) {
    if (count < 0 || from < 0 || from + count > itemCount) {
        Debug(Debug::ERROR) << "Work queue " << path << " asked for records [" << from << ", "
                            << (from + count) << ") of " << itemCount << "\n";
        EXIT(EXIT_FAILURE);
    }
    out.resize(static_cast<size_t>(count));
    if (count == 0) {
        return;
    }
    const size_t bytes = static_cast<size_t>(count) * sizeof(Record);
    if (preadFully(lock.getFd(), &out[0], bytes, recordOffset(from)) != static_cast<ssize_t>(bytes)) {
        Debug(Debug::ERROR) << "Could not read records [" << from << ", " << (from + count)
                            << ") of work queue " << path << ": " << strerror(errno) << "\n";
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

    // Scanned in bounded windows from nextHint rather than reading the whole
    // record array. Both parts matter: starting at nextHint skips the finished
    // prefix, and stopping after SCAN_WINDOW records bounds the read even when
    // the prefix has not advanced -- which is the case while the queue is being
    // worked, since nextHint only moves past items that are *DONE*.
    std::vector<Record> records;
    int64_t firstUnfinished = -1;
    int64_t claimed = -1;
    for (int64_t base = static_cast<int64_t>(header.nextHint); base < itemCount && claimed < 0;
             base += SCAN_WINDOW) {
        const int64_t count = std::min<int64_t>(SCAN_WINDOW, itemCount - base);
        readRecordsLocked(base, count, records);
        for (int64_t offset = 0; offset < count; offset++) {
            const int64_t index = base + offset;
            Record record = records[static_cast<size_t>(offset)];
            if (record.state == DONE) {
                continue;
            }
            if (firstUnfinished < 0) {
                firstUnfinished = index;
            }
            // Untouched, or held by a worker whose lease ran out -- in both cases
            // the item is ours to take. Re-claiming after a lease expiry is what
            // makes a killed job recoverable without an operator deciding what to
            // redo.
            if (record.state == PENDING
                    || (record.state == CLAIMED
                        && static_cast<int64_t>(record.leaseExpiry) <= now)) {
                record.state = CLAIMED;
                record.worker = static_cast<uint32_t>(workerId);
                record.leaseExpiry = static_cast<uint64_t>(now + leaseSeconds);
                writeRecordLocked(index, record);
                claimed = index;
                break;
            }
        }
    }

    if (firstUnfinished > static_cast<int64_t>(header.nextHint)) {
        header.nextHint = static_cast<uint64_t>(firstUnfinished);
        writeHeaderLocked(header);
    }
    // Kept, against the reviewers' suggestion to drop it. It is not durability
    // against a crash -- the stage is re-run for that -- it is *visibility*: the
    // record was pwritten into the page cache, and a worker on another node that
    // takes the lock next must see it or it will claim the same item. fcntl
    // lock/unlock does imply a cache flush on a conformant NFSv4 client, but the
    // cost of that assumption being wrong is two workers silently running one item.
    // The O(itemCount) read this used to do alongside it was the real cost, and
    // that is what the windowed scan above removes.
    fsync(lock.getFd());

    lock.unlock();
    return claimed;
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
    if (record.state == CLAIMED && record.worker != static_cast<uint32_t>(workerId)) {
        // Someone else holds it. Our lease lapsed -- a stalled node, a long GC of
        // the filesystem, a clock skew -- and another worker re-claimed the item
        // and is running it *now*. Recording DONE here would tell the stage an
        // in-progress item is finished, and for the reduce it would also make this
        // worker the authority for that partition's edge blocks while the worker
        // actually producing them is not.
        //
        // renew() and release() have always guarded on ownership; only complete()
        // did not. Returning silently is right: our output for the item is
        // superseded by the holder's, which is exactly what the EdgeBlockHeader
        // filter is built to express.
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
    // From nextHint, in windows, stopping at the first live claim. Every idle
    // worker calls this every pollSeconds at the tail of a stage, so reading the
    // whole array here was the second O(itemCount) read under the global lock --
    // and the one that most directly crowded out the heartbeat renewals live
    // workers depend on.
    const Header header = readHeaderLocked();
    std::vector<Record> records;
    bool live = false;
    for (int64_t base = static_cast<int64_t>(header.nextHint); base < itemCount && live == false;
             base += SCAN_WINDOW) {
        const int64_t count = std::min<int64_t>(SCAN_WINDOW, itemCount - base);
        readRecordsLocked(base, count, records);
        for (int64_t i = 0; i < count; i++) {
            const Record &record = records[static_cast<size_t>(i)];
            if (record.state == CLAIMED && static_cast<int64_t>(record.leaseExpiry) > now) {
                live = true;
                break;
            }
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
    // Read in blocks, not one 16-byte pread per item: the align stage calls this
    // once per reduce wave over a queue that has up to 1e6 items at the 1e12
    // sizing, and a pread each was 1e6 round trips before any edge was read.
    std::vector<Record> block;
    for (uint64_t from = 0; from < header.itemCount; from += static_cast<uint64_t>(SCAN_WINDOW)) {
        const uint64_t count =
            std::min<uint64_t>(static_cast<uint64_t>(SCAN_WINDOW), header.itemCount - from);
        block.resize(static_cast<size_t>(count));
        const size_t bytes = static_cast<size_t>(count) * sizeof(Record);
        if (preadFully(fd, &block[0], bytes, recordOffset(static_cast<int64_t>(from)))
                != static_cast<ssize_t>(bytes)) {
            close(fd);
            Debug(Debug::ERROR) << "Cannot read records [" << from << ", " << (from + count)
                                << ") of work queue " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (uint64_t i = 0; i < count; i++) {
            if (block[static_cast<size_t>(i)].state == DONE) {
                workers[static_cast<size_t>(from + i)] =
                    static_cast<int64_t>(block[static_cast<size_t>(i)].worker);
            }
        }
    }
    close(fd);
    return true;
}

bool WorkQueue::awaitAll(unsigned int pollSeconds, unsigned int stallSeconds) {
    int64_t lastDone = -1;
    int64_t lastProgress = nowSeconds();
    // Backs off like drain(), and for the same reason: every waiting worker takes
    // the global lock once per poll, and at the tail of a stage that is thousands
    // of lock acquisitions competing with the heartbeats of the workers still
    // running.
    unsigned int wait = pollSeconds;
    const unsigned int maxWait = 60;
    while (true) {
        int64_t done = getDoneCount();
        if (done >= itemCount) {
            return true;
        }
        if (done != lastDone) {
            lastDone = done;
            lastProgress = nowSeconds();
            wait = pollSeconds;
        } else if (stallSeconds > 0
                   && nowSeconds() - lastProgress > static_cast<int64_t>(stallSeconds)) {
            Debug(Debug::WARNING) << "Work queue " << path << " made no progress for "
                                  << stallSeconds << " s (" << done << "/" << itemCount
                                  << " done)\n";
            return false;
        }
        sleep(wait);
        if (wait < maxWait) {
            wait = std::min(maxWait, wait * 2);
        }
    }
}
