#ifndef MMSEQS_PARALLELCOORDINATION_H
#define MMSEQS_PARALLELCOORDINATION_H

#include <cstdint>
#include <mutex>
#include <string>

// Shared-filesystem coordination for multi-node MMseqs2 stages.
//
// Every worker process runs the *same* command line and coordinates only through
// files in a shared directory: there is no MPI, no rank argument, and no direct
// node-to-node communication. Workers may join late, die, and be restarted; the
// on-disk state is the single source of truth.
//
// All mutual exclusion uses POSIX fcntl() whole-file write locks. That choice is
// deliberate: fcntl locks are honoured cluster-wide on GPFS, Lustre (mounted with
// -o flock) and NFSv4, whereas flock() is node-local on several of them and would
// silently give no protection at all. Two consequences of fcntl semantics drive
// the implementation below:
//   - locks are owned by the *process*, not the thread, so a process-local mutex
//     is needed as well to serialise threads within one worker;
//   - locks are dropped when *any* file descriptor to the file is closed, so the
//     descriptor is opened once and owned for the object's lifetime.
// Locks are released automatically when a process dies, so a crashed worker can
// never deadlock the run.
class FileLock {
public:
    // Opens (creating if needed) the lock file and keeps the descriptor for the
    // object's lifetime. Does not acquire the lock.
    explicit FileLock(const std::string &path);
    ~FileLock();

    // Blocks until the exclusive lock is held. Acquires the process-local mutex
    // first so that concurrent threads in this process serialise against each
    // other rather than all believing they hold the (per-process) fcntl lock.
    void lock();
    void unlock();

    int getFd() const { return fd; }

private:
    FileLock(const FileLock &);
    FileLock &operator=(const FileLock &);

    std::string path;
    int fd;
    std::mutex threadMutex;
};

// A single 64-bit integer in a shared file, updated atomically across nodes.
//
// Used for worker-id assignment (every worker calls fetchAdd once at startup and
// takes the returned value as its identity, so identities stay unique without any
// central authority) and as a completion counter for barriers.
class SharedCounter {
public:
    explicit SharedCounter(const std::string &path);

    // Adds n and returns the value *before* the addition.
    int64_t fetchAdd(int64_t n = 1);
    int64_t get();

    // Polls until the counter reaches at least target. Intended for coarse
    // stage barriers; callers that need failure detection should use WorkQueue,
    // which tracks leases, rather than blocking here indefinitely.
    void await(int64_t target, unsigned int pollSeconds = 1);

private:
    int64_t readLocked();
    void writeLocked(int64_t value);

    FileLock lock;
};

// A crash-tolerant work queue over a fixed set of items numbered [0, itemCount).
//
// Items are claimed under a lease rather than simply popped. A worker that dies
// holding a claim has its lease expire, after which another worker re-claims the
// item. This is the property a 24 h walltime makes mandatory: a stage that
// outlives its Slurm job must be resumable by a fresh set of workers with no
// manual intervention and no lost or duplicated work.
//
// The queue lives in one file: a fixed-size header followed by one fixed-width
// record per item. Fixed widths mean claiming never rewrites the whole file, and
// state is read back by offset rather than by scanning lines.
class WorkQueue {
public:
    enum State {
        PENDING = 0,
        CLAIMED = 1,
        DONE = 2
    };

    // Opens the queue at path, creating it for itemCount items if it does not
    // exist. Safe to call concurrently from every worker: the first to get the
    // lock initialises the file and the rest observe the finished state.
    //
    // Re-opening an existing queue keeps its progress, which is what makes a
    // stage resumable across jobs. Re-opening with a different itemCount is a
    // programming error and exits, since it would mean the work was partitioned
    // differently than in the run being resumed.
    WorkQueue(const std::string &path, int64_t itemCount);

    // Claims the lowest-numbered item that is either untouched or whose holder's
    // lease has expired. Returns the item index, or -1 when no item is currently
    // claimable (either everything is done, or every remaining item is held by a
    // worker with a live lease).
    int64_t claim(int64_t workerId, int64_t leaseSeconds = DEFAULT_LEASE_SECONDS);

    // Extends the lease on a held item. Long-running items must call this more
    // often than leaseSeconds, otherwise another worker will treat the item as
    // abandoned and start it again.
    void renew(int64_t index, int64_t workerId, int64_t leaseSeconds = DEFAULT_LEASE_SECONDS);

    // Marks a held item finished. Idempotent, so a worker that completes an item
    // and dies before recording it can safely redo and re-complete it.
    void complete(int64_t index, int64_t workerId);

    // Releases a held item back to PENDING without completing it, so it is
    // re-claimable immediately rather than after the lease expires.
    void release(int64_t index, int64_t workerId);

    int64_t getItemCount() const { return itemCount; }
    int64_t getDoneCount();
    bool allDone();

    // Polls until every item is DONE. Returns false if no progress was made and
    // no item was claimable for stallSeconds, which means the remaining work is
    // held by workers that are gone but whose leases have not yet expired, or
    // that the run is genuinely stuck.
    bool awaitAll(unsigned int pollSeconds = 5, unsigned int stallSeconds = 0);

    static const int64_t DEFAULT_LEASE_SECONDS = 1800;

private:
    // On-disk layout. Both structs are written and read verbatim; sizes are
    // asserted at construction so a mismatched build cannot corrupt a queue
    // written by another binary.
    struct Header {
        uint64_t magic;
        uint64_t version;
        uint64_t itemCount;
        uint64_t doneCount;
        // Lowest index that might not be DONE yet. Purely an optimisation: it
        // lets claim() skip a finished prefix instead of rescanning from 0.
        uint64_t nextHint;
        uint64_t reserved[3];
    };

    struct Record {
        uint32_t state;
        uint32_t worker;
        // Unix seconds after which the claim is considered abandoned. Node
        // clocks only need to agree to within a fraction of the lease, which
        // NTP-synchronised HPC nodes comfortably do.
        uint64_t leaseExpiry;
    };

    static const uint64_t MAGIC = 0x4d4d51554555453fULL;  // "MMQUEUE?"
    static const uint64_t VERSION = 1;

    void initialiseLocked();
    Header readHeaderLocked();
    void writeHeaderLocked(const Header &header);
    Record readRecordLocked(int64_t index);
    void writeRecordLocked(int64_t index, const Record &record);
    void completeLocked(int64_t index, int64_t workerId);

    static size_t recordOffset(int64_t index) {
        return sizeof(Header) + static_cast<size_t>(index) * sizeof(Record);
    }

    std::string path;
    int64_t itemCount;
    FileLock lock;
};

#endif
