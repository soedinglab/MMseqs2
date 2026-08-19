#ifndef MMSEQS_PARALLELCOORDINATION_H
#define MMSEQS_PARALLELCOORDINATION_H

#include <algorithm>
#include <cstdint>
#include <ctime>
#include <unistd.h>
#include <mutex>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <string>
#include <thread>
#include <vector>

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
//     descriptor must be owned in exactly one place.
// Locks are released automatically when a process dies, so a crashed worker can
// never deadlock the run.
//
// Both consequences are per *file*, not per object, so the descriptor and the
// mutex live in a process-global registry keyed on the canonical path (see
// LockRegistry in the .cpp). Two FileLock instances naming one path therefore
// serialise against each other, and neither can close the descriptor out from
// under the other.
class FileLock {
public:
    // Joins (creating the file if needed) the process-wide lock state for this
    // path. Does not acquire the lock.
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
    // Both borrowed from the per-path registry entry, which outlives this object.
    int fd;
    void *entry;
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

    // Reads a finished queue's completion records without opening it as a queue,
    // which would create it and needs an itemCount the caller does not know.
    //
    // workers[i] is the worker that recorded item i DONE, or -1 if it never was.
    // Exactly one worker is named per completed item -- complete() is idempotent
    // and keeps the first recorded worker -- which makes this usable as the
    // authority on whose output for an item counts, when an item may have been
    // run more than once. Returns false if the queue file does not exist.
    static bool readCompletedWorkers(const std::string &path, std::vector<int64_t> &workers);

    // True while some item is held under an unexpired lease.
    bool hasLiveClaim();

    // Polls until every item is DONE. Returns false if no progress was made and
    // no item was claimable for stallSeconds, which means the remaining work is
    // held by workers that are gone but whose leases have not yet expired, or
    // that the run is genuinely stuck.
    //
    // Prefer drain(): waiting here does not re-claim, so a worker that reaches
    // this call while a crashed worker's lease is still live will wait for items
    // that nobody is going to finish.
    bool awaitAll(unsigned int pollSeconds = 5, unsigned int stallSeconds = 0);

    // Claims items and runs body on each until every item in the queue is DONE.
    //
    // The loop re-enters claiming after waiting, which is the whole point: when
    // claim() returns -1 the work is not finished, it merely means every unfinished
    // item is held under a live lease. If the workers holding them are alive, they
    // will complete them; if they died, the leases expire and the items become this
    // worker's to redo. Stopping at the first -1 instead -- claim, then wait for
    // the rest -- leaves a crashed worker's items unfinished forever, since nobody
    // is left to complete them and waiting alone never re-claims.
    //
    // Returns false if no item completed anywhere in the run for stallSeconds,
    // which by then means the run really is stuck rather than merely waiting on a
    // lease.
    // pollMillis is the *first* idle wait, not the only one: it doubles up to a
    // minute, so a short value costs a few quick polls at the tail of a stage and
    // nothing thereafter.
    //
    // It is milliseconds rather than seconds because a second is not a rounding
    // error for every stage. The map, reduce and align run items for minutes to
    // hours and never notice; the stages that were single-node until recently
    // have phases that finish in well under a second, and there a one-second floor
    // on the straggler wake-up was the dominant cost -- four phases of it turned a
    // 0.34 s merge into 4.7 s at 1e6. The default is unchanged in effect.
    template <typename Body>
    bool drain(int64_t workerId, Body body, unsigned int pollMillis = 5000,
               unsigned int stallSeconds = 2 * DEFAULT_LEASE_SECONDS) {
        int64_t lastDone = -1;
        int64_t lastProgress = static_cast<int64_t>(time(NULL));
        // Idle polling backs off. At the tail of a stage every worker that has run
        // out of claimable items sits in the branch below, and each pass costs a
        // claim() and a hasLiveClaim() through the one global lock -- thousands of
        // workers hammering the lock exactly while the last few holders need it for
        // their heartbeats. Doubling up to a minute keeps a late joiner responsive
        // (the first waits are still pollMillis) without that pile-up. Reset on
        // every observed completion, so a queue that starts moving again is picked
        // up promptly.
        unsigned int wait = pollMillis;
        const unsigned int maxWait = 60 * 1000;
        while (true) {
            const int64_t item = claim(workerId);
            if (item >= 0) {
                wait = pollMillis;
                // Heartbeat for the duration of the item. Without it any item
                // taking longer than the lease -- which at 1e11 is every reduce
                // partition and every align bucket, both hours of work -- looks
                // abandoned, gets re-claimed, and is then run by two workers at
                // once. The stage's output is not written per worker, so that is
                // silent corruption rather than merely wasted effort.
                // Woken on a condition variable, not by polling sleep(1).
                //
                // The sleep loop this replaces cost ~0.5 s of dead time on *every*
                // work item: join() cannot return until the heartbeat thread comes
                // back from its current one-second sleep and notices the flag, so
                // the wait was uniform on [0, 1] s however short the item was.
                // Measured on the map at 1M -- 20 items, 128 threads -- the queue
                // accounted for 14.6 s of a 20.2 s stage against 2.8 s of actual
                // extraction, and the cost was identical at 1 thread and at 128,
                // for every partition count and every --max-seq-len, which is what
                // eventually identified it. Every stage that drains a queue paid
                // it; the map paid most because it has the most items.
                std::atomic<bool> running(true);
                std::mutex heartbeatMutex;
                std::condition_variable heartbeatStop;
                std::thread heartbeat([this, item, workerId, &running, &heartbeatMutex,
                                       &heartbeatStop]() {
                    const int64_t every = std::max<int64_t>(DEFAULT_LEASE_SECONDS / 3, 1);
                    while (true) {
                        {
                            std::unique_lock<std::mutex> guard(heartbeatMutex);
                            // Returns true only if the predicate held, i.e. we were
                            // told to stop; a timeout means the lease needs renewing.
                            if (heartbeatStop.wait_for(guard, std::chrono::seconds(every),
                                                       [&running]() { return running.load() == false; })) {
                                return;
                            }
                        }
                        // Renewed outside the mutex: it takes the queue's file lock,
                        // and holding a local mutex across that would make the
                        // shutdown notify below wait on unrelated I/O.
                        renew(item, workerId);
                    }
                });
                body(static_cast<size_t>(item));
                {
                    std::lock_guard<std::mutex> guard(heartbeatMutex);
                    running.store(false);
                }
                heartbeatStop.notify_all();
                heartbeat.join();
                complete(item, workerId);
                continue;
            }
            const int64_t done = getDoneCount();
            if (done >= itemCount) {
                return true;
            }
            if (done != lastDone) {
                lastDone = done;
                lastProgress = static_cast<int64_t>(time(NULL));
                wait = pollMillis;
            } else if (hasLiveClaim()) {
                // Someone is still working, and heartbeating to say so. A stage
                // whose last item takes hours must not be declared stalled.
                lastProgress = static_cast<int64_t>(time(NULL));
            } else if (stallSeconds > 0 &&
                       static_cast<int64_t>(time(NULL)) - lastProgress >
                           static_cast<int64_t>(stallSeconds)) {
                return false;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(wait));
            if (wait < maxWait) {
                wait = std::min(maxWait, wait * 2);
            }
        }
    }

    // Bounded by how long a dead worker's item should stay unclaimable, not by
    // how long an item takes: drain() heartbeats every DEFAULT_LEASE_SECONDS/3
    // for as long as the item runs, so a live holder never lets its lease lapse.
    // At 1800 s a single killed worker idled an entire measured run for 30
    // minutes before anyone could redo its item.
    static const int64_t DEFAULT_LEASE_SECONDS = 300;

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

    // The comment above has always said the sizes are asserted; they were not.
    // A queue is written by one worker and read by others that may have been
    // built at a different time, so a layout change has to fail the build rather
    // than reinterpret an existing file.
    static_assert(sizeof(Header) == 64, "WorkQueue::Header is written to disk verbatim");
    static_assert(sizeof(Record) == 16, "WorkQueue::Record is written to disk verbatim");

    static const uint64_t MAGIC = 0x4d4d51554555453fULL;  // "MMQUEUE?"
    static const uint64_t VERSION = 1;

    void initialiseLocked();
    Header readHeaderLocked();
    void writeHeaderLocked(const Header &header);
    Record readRecordLocked(int64_t index);
    // A contiguous run of records in one read, for the scans that would otherwise
    // issue one 16-byte pread per item while holding the lock. On a shared
    // filesystem each of those is a cross-node round trip, and at P = 8192 with
    // every idle worker polling every 5 s they crowd out the heartbeat renewals
    // that keep live workers' leases from lapsing.
    //
    // out is indexed *relative to `from`*. Reading [from, from + count) rather
    // than the whole array is what keeps claim() O(1) in the item count: this used
    // to pread from record 0 every time, so a 1,048,576-item queue -- the map's
    // size at 1e12 -- moved 16 MB under the global lock per claim, measured at
    // 13.2 ms against 4.2 ms for a 4096-item queue.
    void readRecordsLocked(int64_t from, int64_t count, std::vector<Record> &out);
    void writeRecordLocked(int64_t index, const Record &record);
    void completeLocked(int64_t index, int64_t workerId);

    // Records read per pass when scanning for a claimable item. Bounds the read
    // under the lock; the scan extends by another window only when the whole
    // window is DONE or live-leased, which happens on a queue whose head is a long
    // finished prefix -- and nextHint then skips it on the next call.
    static const int64_t SCAN_WINDOW = 4096;

    static size_t recordOffset(int64_t index) {
        return sizeof(Header) + static_cast<size_t>(index) * sizeof(Record);
    }

    std::string path;
    int64_t itemCount;
    FileLock lock;
};

#endif
