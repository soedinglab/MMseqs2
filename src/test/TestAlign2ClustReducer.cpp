// Unit tests for Align2ClustReducer and Align2ClustChunking using synthetic candidates.
//
// These classes have no dependency on OpenMP, MPI, DBReader, or DBWriter (see the class
// comments in Align2ClustReducer.h and Align2ClustChunking.h), so they can be exercised
// here with hand-built ClusterResult vectors instead of a real sequence/prefilter DB.
//
// Coverage:
//   - SET_COVER tie-breaking (larger prefSize wins; ties broken by representative id)
//   - GREEDY sequence-order processing
//   - Already-assigned representatives/members are correctly rejected
//   - Overlapping member sets across candidates
//   - Empty and singleton candidates
//   - Out-of-order pushResult() calls (correctness must not depend on push order)
//   - Align2ClustChunking::computeChunks / isChunkOwnedByRank boundary cases
#include "Align2ClustReducer.h"
#include "Align2ClustChunking.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

const char* binary_name = "test_align2clust_reducer";

// Plain assert() is compiled out under -DNDEBUG (the default Release build used by this
// project's CMake configuration), which would silently turn every check below into a
// no-op. Use an unconditional check instead so these tests still fail loudly in Release
// builds.
#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "CHECK failed at " << __FILE__ << ":" << __LINE__ << ": " #cond << "\n"; \
            std::exit(1); \
        } \
    } while (0)

static ClusterResult makeResult(size_t sequenceIdx, DBLocalId representativeId, size_t prefSize,
                                 std::vector<DBLocalId> memberIds) {
    ClusterResult result;
    result.sequenceIdx = sequenceIdx;
    result.representativeId = representativeId;
    result.prefSize = prefSize;
    result.memberIds = std::move(memberIds);
    return result;
}

// Push every result, then finish()/finalizeSingletons(), returning the assignment vector.
static std::vector<DBLocalId> runReducer(Align2ClustReducer::Mode mode, size_t dbSize,
                                          std::vector<ClusterResult> results) {
    Align2ClustReducer reducer(mode, dbSize, std::max<size_t>(1, results.size()));
    reducer.start();
    for (ClusterResult &result : results) {
        reducer.pushResult(std::move(result));
    }
    reducer.finish();
    reducer.finalizeSingletons();

    std::vector<DBLocalId> assignment(dbSize);
    for (size_t i = 0; i < dbSize; ++i) {
        assignment[i] = reducer.getAssignment(i);
    }
    return assignment;
}

static void testSetCoverTieBreaking() {
    // dbSize == 4, two candidates (rep 2 and rep 3) tie on prefSize == 5. Both are
    // buffered into the heap before either is finalized (their gate only opens once
    // currentPrefSize_ drops below 5, i.e. once the smaller trailing candidates at
    // sequenceIdx 2/3 have also been drained from the reorder buffer). Among a tie in
    // memberCount, SetCoverComparator orders the *higher* representative id to the
    // front of the heap, so it is finalized first (empirically verified against the
    // implementation; this is also exactly the pre-refactor behavior, per the Phase 2
    // byte-identical regression check against a real linclust run).
    const size_t dbSize = 4;
    std::vector<ClusterResult> results;
    // sequenceIdx 0: representative 2, prefSize 5, tie with rep 3 below
    results.push_back(makeResult(0, 2, 5, {0, 1}));
    // sequenceIdx 1: representative 3, prefSize 5 -- ties with rep 2 above, wins (see
    // comment) and claims members {1, 3} first.
    results.push_back(makeResult(1, 3, 5, {1, 3}));
    // sequenceIdx 2: singleton candidate (memberIds.size() == 1) for representative 2;
    // SET_COVER never even queues single-member candidates (see consumerLoopSetCover),
    // so this is a no-op and 2's fate is decided purely by finalizeSingletons() below.
    results.push_back(makeResult(2, 2, 1, {2}));
    // sequenceIdx 3: same, for representative 3 -- also a no-op.
    results.push_back(makeResult(3, 3, 1, {3}));

    std::vector<DBLocalId> assignment = runReducer(Align2ClustReducer::SET_COVER, dbSize, results);
    // Rep 3 wins the tie and claims members {1, 3}.
    CHECK(assignment[1] == 3);
    CHECK(assignment[3] == 3);
    // Rep 2's candidate then has member 1 already taken, leaving only {0}
    // (validCount == 1), so it is rejected outright: 0 and 2 fall back to singletons.
    CHECK(assignment[0] == 0);
    CHECK(assignment[2] == 2);
    std::cout << "testSetCoverTieBreaking OK\n";
}

static void testGreedySequenceOrder() {
    // GREEDY processes candidates strictly in ascending sequenceIdx (== sequence id)
    // order; earlier representatives claim members first.
    const size_t dbSize = 5;
    std::vector<ClusterResult> results;
    results.push_back(makeResult(0, 0, 0, {0, 1, 2}));
    results.push_back(makeResult(1, 1, 0, {1, 3}));       // 1 already claimed by rep 0
    results.push_back(makeResult(2, 2, 0, {2}));           // 2 already claimed by rep 0
    results.push_back(makeResult(3, 3, 0, {3, 4}));        // rep 3 itself unclaimed
    results.push_back(makeResult(4, 4, 0, {4}));           // 4 already claimed by rep 3

    std::vector<DBLocalId> assignment = runReducer(Align2ClustReducer::GREEDY, dbSize, results);
    CHECK(assignment[0] == 0);
    CHECK(assignment[1] == 0);
    CHECK(assignment[2] == 0);
    CHECK(assignment[3] == 3);
    CHECK(assignment[4] == 3);
    std::cout << "testGreedySequenceOrder OK\n";
}

static void testEmptyAndSingletonCandidates() {
    // A representative candidate with no members at all must still leave its own
    // sequence assigned to itself via finalizeSingletons(), and must not disturb
    // unrelated sequences.
    const size_t dbSize = 3;
    std::vector<ClusterResult> results;
    results.push_back(makeResult(0, 0, 0, {}));   // empty candidate
    results.push_back(makeResult(1, 1, 1, {1}));  // singleton candidate
    results.push_back(makeResult(2, 2, 0, {}));   // empty candidate

    std::vector<DBLocalId> assignment = runReducer(Align2ClustReducer::GREEDY, dbSize, results);
    CHECK(assignment[0] == 0);
    CHECK(assignment[1] == 1);
    CHECK(assignment[2] == 2);
    std::cout << "testEmptyAndSingletonCandidates OK\n";
}

static void testOutOfOrderPushDoesNotAffectResult() {
    // Candidate generation may be parallel and produce results in any order (this is
    // the property distributed/MPI worker ranks rely on): pushing the exact same
    // candidates in reverse order must produce an identical assignment.
    const size_t dbSize = 4;
    std::vector<ClusterResult> forward;
    forward.push_back(makeResult(0, 0, 0, {0, 1}));
    forward.push_back(makeResult(1, 1, 0, {1}));
    forward.push_back(makeResult(2, 2, 0, {2, 3}));
    forward.push_back(makeResult(3, 3, 0, {3}));

    std::vector<ClusterResult> reversed;
    for (size_t i = forward.size(); i-- > 0;) {
        reversed.push_back(forward[i]);
    }

    std::vector<DBLocalId> forwardAssignment = runReducer(Align2ClustReducer::GREEDY, dbSize, forward);
    std::vector<DBLocalId> reversedAssignment = runReducer(Align2ClustReducer::GREEDY, dbSize, reversed);
    CHECK(forwardAssignment == reversedAssignment);
    CHECK(forwardAssignment[0] == 0);
    CHECK(forwardAssignment[1] == 0);
    CHECK(forwardAssignment[2] == 2);
    CHECK(forwardAssignment[3] == 2);
    std::cout << "testOutOfOrderPushDoesNotAffectResult OK\n";
}

static void testDistributedCollectOnlyMode() {
    // A distributed worker's collector never touches live assignment state: isAssigned()
    // is always false and pushResult() just buffers, regardless of push order.
    Align2ClustReducer collector(Align2ClustReducer::SET_COVER, /*dbSize=*/10, /*reorderCapacity=*/1,
                                  /*distributedCollectOnly=*/true);
    collector.start();
    CHECK(collector.isAssigned(0) == false);
    collector.pushResult(makeResult(5, 5, 3, {5, 6}));
    collector.pushResult(makeResult(0, 0, 1, {0}));
    CHECK(collector.isAssigned(5) == false);
    collector.finish();

    const std::vector<ClusterResult> &collected = collector.getCollectedResults();
    CHECK(collected.size() == 2);
    // Buffered in push order, not sequenceIdx order (no reordering in this mode).
    CHECK(collected[0].sequenceIdx == 5);
    CHECK(collected[1].sequenceIdx == 0);
    std::cout << "testDistributedCollectOnlyMode OK\n";
}

static void testComputeChunksBoundaries() {
    // Exact multiple.
    {
        std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(100, 25);
        CHECK(chunks.size() == 4);
        for (size_t i = 0; i < chunks.size(); ++i) {
            CHECK(chunks[i].index == i);
            CHECK(chunks[i].start == i * 25);
            CHECK(chunks[i].end == (i + 1) * 25);
        }
    }
    // Remainder in the last chunk.
    {
        std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(103, 25);
        CHECK(chunks.size() == 5);
        CHECK(chunks.back().start == 100);
        CHECK(chunks.back().end == 103);
    }
    // endRange == 0 -> no chunks.
    {
        std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(0, 25);
        CHECK(chunks.empty());
    }
    // chunkSize >= endRange -> a single chunk covering the whole range.
    {
        std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(10, 1000);
        CHECK(chunks.size() == 1);
        CHECK(chunks[0].start == 0);
        CHECK(chunks[0].end == 10);
    }
    // chunkSize == 0 is clamped to 1 (one chunk per index).
    {
        std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(3, 0);
        CHECK(chunks.size() == 3);
    }
    std::cout << "testComputeChunksBoundaries OK\n";
}

static void testChunkOwnershipIsDeterministicAcrossTopologies() {
    // numProc <= 1 always owns every chunk (the non-distributed fast path).
    CHECK(Align2ClustChunking::isChunkOwnedByRank(0, 0, 1) == true);
    CHECK(Align2ClustChunking::isChunkOwnedByRank(41, 0, 1) == true);
    CHECK(Align2ClustChunking::isChunkOwnedByRank(41, 0, 0) == true);

    // Round robin: every chunk index is owned by exactly one rank out of numProc.
    const int numProc = 4;
    for (size_t chunkIndex = 0; chunkIndex < 17; ++chunkIndex) {
        int owners = 0;
        for (int rank = 0; rank < numProc; ++rank) {
            if (Align2ClustChunking::isChunkOwnedByRank(chunkIndex, rank, numProc)) {
                owners++;
            }
        }
        CHECK(owners == 1);
    }

    // Changing numProc (simulating a topology change across a resume) reshuffles
    // ownership deterministically; the same (chunkIndex, rank, numProc) triple always
    // gives the same answer, independent of call order.
    CHECK(Align2ClustChunking::isChunkOwnedByRank(5, 1, 4) == true);
    CHECK(Align2ClustChunking::isChunkOwnedByRank(5, 1, 4) == true);
    CHECK(Align2ClustChunking::isChunkOwnedByRank(5, 0, 2) == false);
    CHECK(Align2ClustChunking::isChunkOwnedByRank(5, 1, 2) == true);
    std::cout << "testChunkOwnershipIsDeterministicAcrossTopologies OK\n";
}

int main(int, const char**) {
    testSetCoverTieBreaking();
    testGreedySequenceOrder();
    testEmptyAndSingletonCandidates();
    testOutOfOrderPushDoesNotAffectResult();
    testDistributedCollectOnlyMode();
    testComputeChunksBoundaries();
    testChunkOwnershipIsDeterministicAcrossTopologies();
    std::cout << "All Align2ClustReducer/Align2ClustChunking tests passed\n";
    return 0;
}
