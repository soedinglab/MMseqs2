#include "Align2ClustReducer.h"
#include "Util.h"

#include <algorithm>

// Compact setCoverMemberPool_ when over half of it is dead. Only offsets
// change, so heap order holds.
static const size_t ALIGN2CLUST_MIN_COMPACTION_DEAD_MEMBERS = 16 * 1024 * 1024 / sizeof(DBLocalId);

Align2ClustReducer::Align2ClustReducer(Mode mode, size_t dbSize, size_t reorderCapacity, bool distributedCollectOnly)
    : mode_(mode), dbSize_(dbSize), assignedCluster_(nullptr),
      reorderCapacity_(std::max<size_t>(1, reorderCapacity)), reorderBufferedCount_(0),
      currentProcessPosition_(0), currentPrefSize_(0), allCalculationsDone_(false),
      setCoverLiveMemberCount_(0), started_(false), finished_(false),
      distributedCollectOnly_(distributedCollectOnly) {
    if (distributedCollectOnly_) {
        // No live assignment state is tracked in this mode (see class comment): a
        // distributed worker rank has no cross-rank view of assignment progress, so
        // isAssigned() always reports "unassigned" and no atomic array is needed.
        return;
    }
    assignedCluster_ = new(std::nothrow) std::atomic<DBLocalId>[dbSize_];
    Util::checkAllocation(assignedCluster_, "Can not allocate assignedCluster memory in Align2ClustReducer");
    for (size_t i = 0; i < dbSize_; ++i) {
        storeAssignedCluster(assignedCluster_, i, DB_LOCAL_ID_INVALID);
    }
    reorderSlots_.resize(reorderCapacity_);
    reorderFilled_.assign(reorderCapacity_, 0);
}

Align2ClustReducer::~Align2ClustReducer() {
    delete[] assignedCluster_;
}

void Align2ClustReducer::start() {
    if (distributedCollectOnly_ || started_) {
        started_ = true;
        return;
    }
    started_ = true;
    if (mode_ == SET_COVER) {
        consumerThread_ = std::thread(&Align2ClustReducer::consumerLoopSetCover, this);
    } else {
        consumerThread_ = std::thread(&Align2ClustReducer::consumerLoopGreedy, this);
    }
}

void Align2ClustReducer::pushResult(ClusterResult &&result) {
    if (distributedCollectOnly_) {
        std::lock_guard<std::mutex> lock(mutex_);
        collected_.push_back(std::move(result));
        return;
    }

    const size_t idx = result.sequenceIdx;
    bool shouldNotifyClusterThread = false;
    {
        std::unique_lock<std::mutex> lock(mutex_);
        // Wait for this result's slot to free up. The next-in-order slot is always
        // free, so the producer the consumer is waiting on never blocks (deadlock-free).
        reorderSpaceCondition_.wait(lock, [&] {
            return allCalculationsDone_ || idx < currentProcessPosition_ + reorderCapacity_;
        });
        const size_t slot = idx % reorderCapacity_;
        reorderSlots_[slot] = std::move(result);   // O(1) vector move, no heap sift
        reorderFilled_[slot] = 1;
        reorderBufferedCount_++;
        shouldNotifyClusterThread = (idx == currentProcessPosition_);
    }
    if (shouldNotifyClusterThread) {
        clusterCondition_.notify_one();
    }
}

void Align2ClustReducer::finish() {
    if (finished_) {
        return;
    }
    finished_ = true;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        allCalculationsDone_ = true;
    }
    clusterCondition_.notify_one();
    reorderSpaceCondition_.notify_all();
    if (consumerThread_.joinable()) {
        consumerThread_.join();
    }
}

bool Align2ClustReducer::isAssigned(size_t sequenceId) const {
    if (distributedCollectOnly_) {
        // A distributed worker rank has no live view of cross-rank assignment
        // progress; see the class comment. Always report "unassigned" so callers
        // fall back to full (safe superset) candidate generation.
        return false;
    }
    return loadAssignedCluster(assignedCluster_, sequenceId) != DB_LOCAL_ID_INVALID;
}

void Align2ClustReducer::finalizeSingletons() {
    if (distributedCollectOnly_) {
        return;
    }
    for (size_t i = 0; i < dbSize_; ++i) {
        if (loadAssignedCluster(assignedCluster_, i) == DB_LOCAL_ID_INVALID) {
            storeAssignedCluster(assignedCluster_, i, i);
        }
    }
}

DBLocalId Align2ClustReducer::getAssignment(size_t sequenceId) const {
    return loadAssignedCluster(assignedCluster_, sequenceId);
}

// Compact setCoverMemberPool_ when over half is dead. Only offsets change, so heap order holds.
void Align2ClustReducer::compactSetCoverMemberPool() {
    const size_t deadMemberCount = setCoverMemberPool_.size() - setCoverLiveMemberCount_;
    if (setCoverMemberPool_.empty() ||
        deadMemberCount < ALIGN2CLUST_MIN_COMPACTION_DEAD_MEMBERS ||
        deadMemberCount * 2 <= setCoverMemberPool_.size()) {
        return;
    }
    std::vector<DBLocalId> compactedPool;
    compactedPool.reserve(setCoverLiveMemberCount_);
    for (SetCoverCandidate &candidate : setCoverCandidates_) {
        const size_t newOffset = compactedPool.size();
        compactedPool.insert(compactedPool.end(),
                             setCoverMemberPool_.begin() + candidate.memberOffset,
                             setCoverMemberPool_.begin() + candidate.memberOffset + candidate.memberCount);
        candidate.memberOffset = newOffset;
    }
    setCoverMemberPool_.swap(compactedPool);
}

void Align2ClustReducer::consumerLoopSetCover() {
    while (true) {
        std::unique_lock<std::mutex> lock(mutex_);

        clusterCondition_.wait(lock, [this] {
            return reorderFilled_[currentProcessPosition_ % reorderCapacity_] != 0 ||
                   allCalculationsDone_;
        });

        // 1) reorder buffer -> setCoverCandidates_
        while (reorderFilled_[currentProcessPosition_ % reorderCapacity_] != 0) {
            const size_t slot = currentProcessPosition_ % reorderCapacity_;
            ClusterResult result = std::move(reorderSlots_[slot]);
            reorderFilled_[slot] = 0;
            reorderBufferedCount_--;
            currentProcessPosition_++;
            reorderSpaceCondition_.notify_all();
            currentPrefSize_ = result.prefSize;

            if (result.memberIds.size() > 1) {
                SetCoverCandidate candidate;
                candidate.representativeId = result.representativeId;
                candidate.memberCount = result.memberIds.size();
                candidate.memberOffset = setCoverMemberPool_.size();
                setCoverMemberPool_.insert(setCoverMemberPool_.end(),
                                          result.memberIds.begin(), result.memberIds.end());
                setCoverLiveMemberCount_ += candidate.memberCount;
                setCoverCandidates_.push_back(candidate);
                std::push_heap(setCoverCandidates_.begin(), setCoverCandidates_.end(), SetCoverComparator());
            }
        }

        // 2) assign candidates guaranteed to be the currently largest set
        while (setCoverCandidates_.empty() == false &&
               (allCalculationsDone_ ||
                setCoverCandidates_.front().memberCount > currentPrefSize_)) {

            std::pop_heap(setCoverCandidates_.begin(), setCoverCandidates_.end(), SetCoverComparator());
            SetCoverCandidate candidate = setCoverCandidates_.back();
            setCoverCandidates_.pop_back();
            setCoverLiveMemberCount_ -= candidate.memberCount;

            if (loadAssignedCluster(assignedCluster_, candidate.representativeId) != DB_LOCAL_ID_INVALID) {
                continue;
            }

            // Drop already-assigned members, compacting the survivors in place
            // within the pool region [memberOffset, memberOffset + memberCount).
            DBLocalId *members = setCoverMemberPool_.data() + candidate.memberOffset;
            size_t validCount = 0;
            for (size_t i = 0; i < candidate.memberCount; i++) {
                if (loadAssignedCluster(assignedCluster_, members[i]) == DB_LOCAL_ID_INVALID) {
                    members[validCount++] = members[i];
                }
            }

            if (validCount <= 1) {
                continue;
            }

            if (validCount != candidate.memberCount) {
                candidate.memberCount = validCount;
                setCoverLiveMemberCount_ += validCount;
                setCoverCandidates_.push_back(candidate);
                std::push_heap(setCoverCandidates_.begin(), setCoverCandidates_.end(), SetCoverComparator());
                continue;
            }

            for (size_t i = 0; i < candidate.memberCount; i++) {
                storeAssignedCluster(assignedCluster_, members[i], candidate.representativeId);
            }
        }

        // Compaction only touches consumer-private structures (setCoverCandidates_,
        // setCoverMemberPool_, setCoverLiveMemberCount_), never shared state, so release
        // the mutex during the copy so worker threads can keep pushing results.
        lock.unlock();
        compactSetCoverMemberPool();
        lock.lock();

        if (allCalculationsDone_ &&
            reorderBufferedCount_ == 0 &&
            setCoverCandidates_.empty()) {
            break;
        }
    }
}

void Align2ClustReducer::consumerLoopGreedy() {
    while (true) {
        std::unique_lock<std::mutex> lock(mutex_);

        clusterCondition_.wait(lock, [this] {
            return reorderFilled_[currentProcessPosition_ % reorderCapacity_] != 0 ||
                   allCalculationsDone_;
        });

        if (allCalculationsDone_ && reorderBufferedCount_ == 0) {
            break;
        }

        while (reorderFilled_[currentProcessPosition_ % reorderCapacity_] != 0) {
            const size_t slot = currentProcessPosition_ % reorderCapacity_;
            ClusterResult result = std::move(reorderSlots_[slot]);
            reorderFilled_[slot] = 0;
            reorderBufferedCount_--;
            currentProcessPosition_++;
            reorderSpaceCondition_.notify_all();

            if (loadAssignedCluster(assignedCluster_, result.representativeId) != DB_LOCAL_ID_INVALID) {
                continue;
            }

            std::vector<DBLocalId> validMemberIds;
            validMemberIds.reserve(result.memberIds.size());
            for (DBLocalId memberId : result.memberIds) {
                if (loadAssignedCluster(assignedCluster_, memberId) == DB_LOCAL_ID_INVALID) {
                    validMemberIds.push_back(memberId);
                }
            }

            if (validMemberIds.size() <= 1) {
                continue;
            }

            for (DBLocalId memberId : validMemberIds) {
                storeAssignedCluster(assignedCluster_, memberId, result.representativeId);
            }
        }
    }
}
