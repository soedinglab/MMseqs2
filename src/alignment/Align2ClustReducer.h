#ifndef MMSEQS_ALIGN2CLUST_REDUCER_H
#define MMSEQS_ALIGN2CLUST_REDUCER_H

#include "IndexTypes.h"

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

// One generated cluster candidate for a single global order position.
// Produced by candidate generation (OpenMP workers in doAlign2clust today;
// a distributed generator in a future multi-node setting) and consumed
// strictly in `sequenceIdx` order by Align2ClustReducer.
struct ClusterResult {
    size_t sequenceIdx;
    DBLocalId representativeId;
    size_t prefSize;
    std::vector<DBLocalId> memberIds;
};

// Deterministic, order-preserving reducer for align2clust's SET_COVER and
// GREEDY/GREEDY_MEM clustering modes.
//
// Candidate generation may run in any order, on any number of parallel
// workers (threads today; MPI ranks in a future distributed mode).
// Correctness only requires that every global order position in
// [0, dbSize) is pushed to pushResult() exactly once. Align2ClustReducer
// reorders candidates back into sequenceIdx order internally and performs
// the actual set-cover/greedy assignment strictly in that order, so the
// final clustering never depends on generation order, thread/rank
// scheduling, worker count, or (in a future distributed setting) how many
// nodes produced the candidates.
//
// isAssigned() lets a candidate generator skip alignments it already knows
// would be discarded (a pure performance optimization backed by a live,
// process-local view of the current assignment). It is never required for
// correctness: the reducer independently re-validates every candidate's
// representative and every member against the assignment state at the time
// the candidate is actually consumed, so a generator may safely produce a
// superset of the "true" final members/candidates (for example, because it
// has no live assignment state to consult, as would be the case for a
// distributed worker on a separate node). This is the property that a
// future MPI-distributed candidate generator relies on: worker ranks can
// generate candidates independently, without sharing live assignment
// state, and still feed the same deterministic reducer used here.
//
// This class has no dependency on OpenMP, MPI, DBReader, or DBWriter so it
// can be unit tested with synthetic candidates and reused unchanged by a
// future distributed candidate consumer.
class Align2ClustReducer {
public:
    enum Mode {
        SET_COVER,
        GREEDY
    };

    // dbSize is the total number of sequences (and the number of distinct
    // sequenceIdx values that must eventually be pushed exactly once).
    // reorderCapacity is the maximum out-of-order window; values less than
    // 1 are clamped to 1.
    //
    // distributedCollectOnly selects a lightweight second mode used by MPI
    // worker ranks generating candidates for a single checkpoint chunk: no
    // consumer thread is started, isAssigned() always returns false (a
    // worker rank has no live view of cross-rank assignment state -- see
    // the class comment), and pushResult() simply appends to an internal
    // vector retrievable via getCollectedResults() instead of driving the
    // ordered set-cover/greedy reduction. finish()/finalizeSingletons() are
    // no-ops in this mode; getAssignment() must not be called.
    Align2ClustReducer(Mode mode, size_t dbSize, size_t reorderCapacity, bool distributedCollectOnly = false);
    ~Align2ClustReducer();

    Align2ClustReducer(const Align2ClustReducer &) = delete;
    Align2ClustReducer &operator=(const Align2ClustReducer &) = delete;

    // Starts the internal consumer thread. Must be called exactly once,
    // before any call to pushResult(). No-op in distributed-collect-only
    // mode.
    void start();

    // Thread-safe; may be called concurrently by any number of worker
    // threads. Blocks if the reorder window is full until the reducer has
    // consumed enough in-order results to make room. Must not be called
    // after finish(). In distributed-collect-only mode, simply appends
    // under a mutex and never blocks.
    void pushResult(ClusterResult &&result);

    // Signals that no further results will be pushed for this run, then
    // waits for the internal consumer thread to finish and joins it. Must
    // be called exactly once, after every pushResult() call has returned.
    // No-op in distributed-collect-only mode.
    void finish();

    // True if sequenceId already has a representative assigned. See the
    // class comment: this is a pure performance optimization for candidate
    // generation and is never required for correctness. Always false in
    // distributed-collect-only mode.
    bool isAssigned(size_t sequenceId) const;

    // Must be called after finish(). Assigns every still-unassigned
    // sequence to itself (singleton cluster). Must not be called in
    // distributed-collect-only mode.
    void finalizeSingletons();

    // Valid only after finalizeSingletons(); never DB_LOCAL_ID_INVALID.
    DBLocalId getAssignment(size_t sequenceId) const;

    // Valid only in distributed-collect-only mode, after finish().
    const std::vector<ClusterResult> &getCollectedResults() const { return collected_; }

    size_t getDbSize() const { return dbSize_; }

private:
    struct SetCoverCandidate {
        DBLocalId representativeId;
        size_t memberCount;
        size_t memberOffset;
    };

    struct SetCoverComparator {
        bool operator()(const SetCoverCandidate &a, const SetCoverCandidate &b) const {
            if (a.memberCount < b.memberCount) {
                return true;
            }
            if (b.memberCount < a.memberCount) {
                return false;
            }
            if (a.representativeId < b.representativeId) {
                return true;
            }
            if (b.representativeId < a.representativeId) {
                return false;
            }
            return false;
        }
    };

    void consumerLoopSetCover();
    void consumerLoopGreedy();
    void compactSetCoverMemberPool();

    static DBLocalId loadAssignedCluster(const std::atomic<DBLocalId> *assignedCluster, size_t sequenceId) {
        return assignedCluster[sequenceId].load(std::memory_order_relaxed);
    }
    static void storeAssignedCluster(std::atomic<DBLocalId> *assignedCluster, size_t sequenceId, DBLocalId representativeId) {
        assignedCluster[sequenceId].store(representativeId, std::memory_order_relaxed);
    }

    Mode mode_;
    size_t dbSize_;
    std::atomic<DBLocalId> *assignedCluster_;

    mutable std::mutex mutex_;
    std::condition_variable clusterCondition_;
    std::condition_variable reorderSpaceCondition_;

    std::vector<ClusterResult> reorderSlots_;
    std::vector<unsigned char> reorderFilled_;
    size_t reorderCapacity_;
    size_t reorderBufferedCount_;

    size_t currentProcessPosition_;
    size_t currentPrefSize_;
    bool allCalculationsDone_;

    std::vector<SetCoverCandidate> setCoverCandidates_;
    std::vector<DBLocalId> setCoverMemberPool_;
    size_t setCoverLiveMemberCount_;

    std::thread consumerThread_;
    bool started_;
    bool finished_;
    bool distributedCollectOnly_;
    std::vector<ClusterResult> collected_;
};

#endif
