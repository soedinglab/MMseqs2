#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "FastSort.h"
#include "BlockAligner.h"
#include "Alignment.h"
#include <atomic>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef __linux__
#include <sys/mman.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#define MIN_SIZE 32

struct ClusterResult {
    size_t sequenceIdx;
    DBLocalId representativeId;
    size_t prefSize;
    std::vector<DBLocalId> memberIds;
};

struct PrefInfo {
    DBLocalId id;
    uint32_t size;

    static bool compareBySizeAndId(const PrefInfo &first, const PrefInfo &second){
        if(first.size > second.size)
            return true;
        if(second.size > first.size)
            return false;
        if(first.id < second.id)
            return true;
        if(second.id < first.id)
            return false;
        return false;
    }
};

// Lightweight entry stored in the set-cover ready queue. Instead of carrying the
// member id vector (as ClusterResult does) it only references the members, which
// are kept in a single shared pool (setCoverMemberPool), by offset and count.
struct SetCoverCandidate {
    DBLocalId representativeId;
    size_t memberCount;
    size_t memberOffset;
};

// entry length per local id, so the heap breaks ties by length like a length sorted reader did
static const unsigned int *setCoverLengths = nullptr;

struct SetCoverComparator {
    bool operator()(const SetCoverCandidate& a, const SetCoverCandidate& b) const {
        if (a.memberCount < b.memberCount) {
            return true;
        }
        if (b.memberCount < a.memberCount) {
            return false;
        }
        const unsigned int lenA = setCoverLengths[a.representativeId];
        const unsigned int lenB = setCoverLengths[b.representativeId];
        if (lenA != lenB) {
            return lenA > lenB;
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

static std::mutex clusterMutex;
static std::condition_variable clusterCondition;
static std::condition_variable reorderSpaceCondition;
// Reorders out-of-order worker results back to sequenceIdx order for the consumer,
// indexed by sequenceIdx % reorderCapacity. The next-in-order slot is always free,
// so producers never deadlock.
static std::vector<ClusterResult> reorderSlots;    // slot storage, size == reorderCapacity
static std::vector<unsigned char> reorderFilled;   // 1 == slot holds an unconsumed result
static size_t reorderCapacity = 0;                 // max out-of-order window
static size_t reorderBufferedCount = 0;
// max-heap by member count, kept in a plain vector so the pool stays compactable
static std::vector<SetCoverCandidate> setCoverCandidates;
static SetCoverComparator setCoverComparator;
// Shared backing store for all member ids referenced by setCoverCandidates.
static std::vector<DBLocalId> setCoverMemberPool;
static size_t setCoverLiveMemberCount = 0;

static size_t currentProcessPosition = 0;
static size_t currentPrefSize = 0;
static bool allCalculationsDone = false;

typedef std::atomic<DBLocalId> ClusterAssignment;
static std::atomic<uint64_t> *assignedFlags = nullptr;

static size_t assignedFlagWordCount(size_t dbSize) {
    return dbSize / 64 + static_cast<size_t>(dbSize % 64 != 0);
}

static DBLocalId loadAssignedCluster(const ClusterAssignment *assignedCluster, size_t sequenceId) {
    return assignedCluster[sequenceId].load(std::memory_order_relaxed);
}

static bool isAssigned(const ClusterAssignment *assignedCluster, size_t sequenceId) {
    (void) assignedCluster;
    const uint64_t bit = uint64_t(1) << (sequenceId & 63);
    return (assignedFlags[sequenceId >> 6].load(std::memory_order_relaxed) & bit) != 0;
}

static void storeAssignedCluster(ClusterAssignment *assignedCluster, size_t sequenceId, DBLocalId representativeId) {
    assignedCluster[sequenceId].store(representativeId, std::memory_order_relaxed);
    if (representativeId != DB_LOCAL_ID_INVALID) {
        std::atomic<uint64_t> &word = assignedFlags[sequenceId >> 6];
        word.store(word.load(std::memory_order_relaxed) | (uint64_t(1) << (sequenceId & 63)),
                   std::memory_order_relaxed);
    }
}

// Serialize a single alignment record and append it to the per-representative buffer.
static void appendAlignmentResult(std::string &alnResultBuffer, char *lineBuffer, const Matcher::result_t &result, bool addBacktrace) {
    size_t len = Matcher::resultToBuffer(lineBuffer, result, addBacktrace);
    alnResultBuffer.append(lineBuffer, len);
}

// Compact setCoverMemberPool when over half is dead. Only offsets change, so heap order holds.
static const size_t ALIGN2CLUST_MIN_COMPACTION_DEAD_MEMBERS = 16 * 1024 * 1024 / sizeof(DBLocalId);

static void compactSetCoverMemberPool() {
    const size_t deadMemberCount = setCoverMemberPool.size() - setCoverLiveMemberCount;
    if (setCoverMemberPool.empty() ||
        deadMemberCount < ALIGN2CLUST_MIN_COMPACTION_DEAD_MEMBERS ||
        deadMemberCount * 2 <= setCoverMemberPool.size()) {
        return;
    }
    std::vector<DBLocalId> compactedPool;
    compactedPool.reserve(setCoverLiveMemberCount);
    for (SetCoverCandidate &candidate : setCoverCandidates) {
        const size_t newOffset = compactedPool.size();
        compactedPool.insert(compactedPool.end(),
                             setCoverMemberPool.begin() + candidate.memberOffset,
                             setCoverMemberPool.begin() + candidate.memberOffset + candidate.memberCount);
        candidate.memberOffset = newOffset;
    }
    setCoverMemberPool.swap(compactedPool);
}

static const size_t ALIGN2CLUST_DEFAULT_REORDER_LIMIT = 2 * 1000 * 1000;

// one drop per this many bytes: DONTNEED shoots down every thread's TLB and a sub-page range frees nothing
static const size_t ALIGN2CLUST_CACHE_DROP_BYTES = 256 * 1024 * 1024;
static const size_t ALIGN2CLUST_CACHE_DROP_STRIDE = 65536;
static std::atomic<size_t> seqCacheDroppedTo(0);
static std::atomic<size_t> alnCacheDroppedTo(0);

// first id a query of this length can still cover; everything before it is out of the band for good
static size_t coverageBandStart(DBReader<DBKeyType> *seqDbr, size_t queryId, float covThr, int covMode) {
    const float queryLength = static_cast<float>(seqDbr->getSeqLen(queryId));
    size_t low = 0;
    size_t high = queryId;
    while (low < high) {
        const size_t mid = low + (high - low) / 2;
        if (Util::canBeCovered(covThr, covMode, queryLength, static_cast<float>(seqDbr->getSeqLen(mid)))) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    return low;
}

static void dropBehind(DBReader<DBKeyType> &reader, std::atomic<size_t> &droppedTo, size_t end) {
    size_t seen = droppedTo.load(std::memory_order_relaxed);
    if (end < seen + ALIGN2CLUST_CACHE_DROP_BYTES) {
        return;
    }
    if (droppedTo.compare_exchange_strong(seen, end, std::memory_order_relaxed)) {
        reader.dropCacheRange(seen, end);
    }
}

// Both uses of the split memory limit here only size internal buffers, and the mapped databases can
// legitimately be larger than the limit. Util::computeMemory exits in that case, which killed
// linclust --split-memory-limit 10M; this degrades to a zero budget instead, which makes the reorder
// buffer take its minimum and the io policy report a starved cache -- both correct under a tight limit.
static size_t availableBufferBudget(size_t splitMemoryLimit) {
    const size_t limit = (splitMemoryLimit > 0)
        ? splitMemoryLimit
        : static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    const size_t committed = MemoryTracker::getSize();
    return (committed < limit) ? (limit - committed) : 0;
}

static size_t computeReorderCapacity(const Parameters &par, size_t dbSize, int mode, size_t resultCount,
                                     bool ownsLengthOrder) {
    const size_t memoryLimit = availableBufferBudget(par.splitMemoryLimit);

    size_t fixedMemory = dbSize * sizeof(ClusterAssignment);
    fixedMemory += assignedFlagWordCount(dbSize) * sizeof(std::atomic<uint64_t>);
    if (mode == Parameters::SET_COVER) {
        fixedMemory += dbSize * sizeof(PrefInfo);
    } else if (ownsLengthOrder) {
        // GREEDY keeps this visit order until every producer has completed.
        fixedMemory += dbSize * sizeof(DBLocalId);
    }

    const size_t budget = (memoryLimit > fixedMemory)
        ? static_cast<size_t>(static_cast<double>(memoryLimit - fixedMemory) * 0.9)
        : 0;
    const size_t bytesPerResult = sizeof(ClusterResult) + sizeof(unsigned char)
                                  + (par.maxResListLen + 1) * sizeof(DBLocalId)
                                  + ((mode == Parameters::SET_COVER) ? 0 : sizeof(ClusterResult));

    size_t capacity = std::min(resultCount, std::max<size_t>(1, budget / bytesPerResult));
    capacity = std::min(capacity, ALIGN2CLUST_DEFAULT_REORDER_LIMIT);
    capacity = std::max<size_t>(1, capacity);

    return capacity;
}

static void discardLengthOrderPrefix(DBLocalId *lengthOrder, size_t begin, size_t end) {
#ifdef MADV_DONTNEED
    const size_t pageSize = Util::getPageSize();
    const uintptr_t beginAddress = reinterpret_cast<uintptr_t>(lengthOrder + begin);
    const uintptr_t endAddress = reinterpret_cast<uintptr_t>(lengthOrder + end);
    const uintptr_t discardBegin = (beginAddress + pageSize - 1) & ~(static_cast<uintptr_t>(pageSize) - 1);
    const uintptr_t discardEnd = endAddress & ~(static_cast<uintptr_t>(pageSize) - 1);
    if (discardEnd > discardBegin) {
        madvise(reinterpret_cast<void *>(discardBegin), discardEnd - discardBegin, MADV_DONTNEED);
    }
#else
    (void) lengthOrder;
    (void) begin;
    (void) end;
#endif
}

static void pushClusterResult(ClusterResult &&clusterResult) {
    const size_t idx = clusterResult.sequenceIdx;
    bool shouldNotifyClusterThread = false;
    {
        std::unique_lock<std::mutex> lock(clusterMutex);
        // Wait for this result's slot to free up. The next-in-order slot is always
        // free, so the producer the consumer is waiting on never blocks (deadlock-free).
        reorderSpaceCondition.wait(lock, [&] {
            return idx < currentProcessPosition + reorderCapacity;
        });
        const size_t slot = idx % reorderCapacity;
        reorderSlots[slot] = std::move(clusterResult);   // O(1) vector move, no heap sift
        reorderFilled[slot] = 1;
        reorderBufferedCount++;
        shouldNotifyClusterThread = (idx == currentProcessPosition);
    }
    if (shouldNotifyClusterThread) {
        clusterCondition.notify_one();
    }
}

static float parsePrecisionLib(const std::string &scoreFile, double targetSeqid, double targetCov, double targetPrecision) {
    std::stringstream in(scoreFile);
    std::string line;
    int intTargetSeqid = static_cast<int>((targetSeqid + 0.0001) * 100);
    int seqIdRest = (intTargetSeqid % 5);
    targetSeqid = static_cast<float>(intTargetSeqid - seqIdRest) / 100;
    targetCov = static_cast<float>(static_cast<int>((targetCov + 0.0001) * 10)) / 10;
    
    while (std::getline(in, line)) {
        std::vector<std::string> values = Util::split(line, " ");
        float cov = strtod(values[0].c_str(), NULL);
        float seqid = strtod(values[1].c_str(), NULL);
        float scorePerCol = strtod(values[2].c_str(), NULL);
        float precision = strtod(values[3].c_str(), NULL);
        if (MathUtil::AreSame(cov, targetCov) && MathUtil::AreSame(seqid, targetSeqid) && precision >= targetPrecision) {
            return scorePerCol;
        }
    }
    
    Debug(Debug::WARNING) << "Can not find any score per column for coverage "
                          << targetCov << " and sequence identity " << targetSeqid 
                          << ". No hit will be filtered.\n";
    return 0;
}

// memberOrder keeps a cluster contiguous, so nudging a split off the cluster it lands in is enough
static size_t clusterRangeStart(const ClusterAssignment *assignedCluster, const DBLocalId *memberOrder,
                                size_t dbSize, size_t at) {
    if (at >= dbSize) {
        return dbSize;
    }
    while (at > 0 && loadAssignedCluster(assignedCluster, memberOrder[at])
                     == loadAssignedCluster(assignedCluster, memberOrder[at - 1])) {
        at++;
        if (at >= dbSize) {
            return dbSize;
        }
    }
    return at;
}

// mirrors Clustering::writeData, but reads members through a local-id permutation
static void writeClustering(DBWriter *dbWriter, DBReader<DBKeyType> *seqDbr,
                            const ClusterAssignment *assignedCluster, const DBLocalId *memberOrder,
                            size_t dbSize, int threads) {
#pragma omp parallel num_threads(threads)
    {
    unsigned int thrIdx = 0;
    int threadCnt = 1;
#ifdef OPENMP
    thrIdx = static_cast<unsigned int>(omp_get_thread_num());
    threadCnt = omp_get_num_threads();
#endif
    // disjoint ascending ranges in ascending thread files, so the concatenation is serial identical
    const size_t stride = (dbSize + threadCnt - 1) / static_cast<size_t>(threadCnt);
    const size_t rangeBegin = clusterRangeStart(assignedCluster, memberOrder, dbSize, stride * thrIdx);
    const size_t rangeEnd = clusterRangeStart(assignedCluster, memberOrder, dbSize, stride * (thrIdx + 1));
    std::string resultString;
    resultString.reserve(1024 * 1024);
    char buffer[32];
    DBKeyType previousRepresentativeKey = DB_KEY_INVALID;

    for (size_t i = rangeBegin; i < rangeEnd; i++) {
        const DBLocalId memberId = memberOrder[i];
        const DBKeyType currentRepresentativeKey = seqDbr->getDbKey(loadAssignedCluster(assignedCluster, memberId));

        if (previousRepresentativeKey != currentRepresentativeKey) {
            if (previousRepresentativeKey != DB_KEY_INVALID) {
                dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey, thrIdx);
            }
            resultString.clear();
            char *outPos = Itoa::u64toa_sse2(static_cast<uint64_t>(currentRepresentativeKey), buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }

        const DBKeyType memberKey = seqDbr->getDbKey(memberId);
        if (memberKey != currentRepresentativeKey) {
            char *outPos = Itoa::u64toa_sse2(static_cast<uint64_t>(memberKey), buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }

        previousRepresentativeKey = currentRepresentativeKey;
    }

    if (previousRepresentativeKey != DB_KEY_INVALID) {
        dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey, thrIdx);
    }
    }
}

static size_t requireId(size_t id, const char *dbName, DBKeyType key) {
    if (id == DB_ENTRY_NOT_FOUND) {
        Debug(Debug::ERROR) << dbName << " has no entry for key " << key << "\n";
        EXIT(EXIT_FAILURE);
    }
    return id;
}

static void (*clusterThreadFunc)(ClusterAssignment*) = nullptr;

void clusterThreadFuncSetcover(ClusterAssignment* assignedCluster) {
    while (true) {
        bool drained = false;
        std::unique_lock<std::mutex> lock(clusterMutex);
        
        clusterCondition.wait(lock, [] {
            return reorderFilled[currentProcessPosition % reorderCapacity] != 0 ||
                   allCalculationsDone;
        });

        // 1) reorder buffer → setCoverCandidates
        while (reorderFilled[currentProcessPosition % reorderCapacity] != 0) {
            const size_t slot = currentProcessPosition % reorderCapacity;
            ClusterResult result = std::move(reorderSlots[slot]);
            reorderFilled[slot] = 0;
            reorderBufferedCount--;
            currentProcessPosition++;
            drained = true;
            currentPrefSize = result.prefSize;

            if (result.memberIds.size() > 1) {
                SetCoverCandidate candidate;
                candidate.representativeId = result.representativeId;
                candidate.memberCount = result.memberIds.size();
                candidate.memberOffset = setCoverMemberPool.size();
                setCoverMemberPool.insert(setCoverMemberPool.end(),
                                          result.memberIds.begin(), result.memberIds.end());
                setCoverLiveMemberCount += candidate.memberCount;
                setCoverCandidates.push_back(candidate);
                std::push_heap(setCoverCandidates.begin(), setCoverCandidates.end(), setCoverComparator);
            }
        }

        // 2) assign candidates guaranteed to be the currently largest set
        while (setCoverCandidates.empty() == false &&
               (allCalculationsDone ||
                setCoverCandidates.front().memberCount > currentPrefSize)) {

            std::pop_heap(setCoverCandidates.begin(), setCoverCandidates.end(), setCoverComparator);
            SetCoverCandidate candidate = setCoverCandidates.back();
            setCoverCandidates.pop_back();
            setCoverLiveMemberCount -= candidate.memberCount;

            if (isAssigned(assignedCluster, candidate.representativeId)) {
                continue;
            }

            // Drop already-assigned members, compacting the survivors in place
            // within the pool region [memberOffset, memberOffset + memberCount).
            DBLocalId *members = setCoverMemberPool.data() + candidate.memberOffset;
            size_t validCount = 0;
            for (size_t i = 0; i < candidate.memberCount; i++) {
                if (!isAssigned(assignedCluster, members[i])) {
                    members[validCount++] = members[i];
                }
            }

            if (validCount <= 1) {
                continue;
            }

            if (validCount != candidate.memberCount) {
                candidate.memberCount = validCount;
                setCoverLiveMemberCount += validCount;
                setCoverCandidates.push_back(candidate);
                std::push_heap(setCoverCandidates.begin(), setCoverCandidates.end(), setCoverComparator);
                continue;
            }

            for (size_t i = 0; i < candidate.memberCount; i++) {
                storeAssignedCluster(assignedCluster, members[i], candidate.representativeId);
            }
        }

        // Compaction only touches consumer-private structures (setCoverCandidates,
        // setCoverMemberPool, setCoverLiveMemberCount), never shared state, so release
        // the mutex during the copy so worker threads can keep pushing results.
        lock.unlock();
        // one wake per drained round, issued while the lock is free
        if (drained) {
            reorderSpaceCondition.notify_all();
        }
        compactSetCoverMemberPool();
        lock.lock();

        if (allCalculationsDone &&
            reorderBufferedCount == 0 &&
            setCoverCandidates.empty()) {
            break;
        }
    }
}

void clusterThreadFuncGreedy(ClusterAssignment* assignedCluster) {
    // consumer-private scratch, so both are reused instead of allocated once per result
    std::vector<ClusterResult> drainedResults;
    std::vector<DBLocalId> validMemberIds;
    while (true) {
        bool lastRound = false;
        drainedResults.clear();
        {
            std::unique_lock<std::mutex> lock(clusterMutex);

            clusterCondition.wait(lock, [] {
                return reorderFilled[currentProcessPosition % reorderCapacity] != 0 ||
                       allCalculationsDone;
            });

            lastRound = (allCalculationsDone && reorderBufferedCount == 0);

            // the mutex guards the ring, so hold it only for the moves out of it
            while (reorderFilled[currentProcessPosition % reorderCapacity] != 0) {
                const size_t slot = currentProcessPosition % reorderCapacity;
                drainedResults.push_back(std::move(reorderSlots[slot]));
                reorderFilled[slot] = 0;
                reorderBufferedCount--;
                currentProcessPosition++;
            }
        }

        // one wake per drained round, issued while the lock is free
        if (drainedResults.empty() == false) {
            reorderSpaceCondition.notify_all();
        }

        // only this thread writes assignedCluster while the producers run, so no mutex is needed
        for (ClusterResult &result : drainedResults) {
            if (isAssigned(assignedCluster, result.representativeId)) {
                continue;
            }

            validMemberIds.clear();
            validMemberIds.reserve(result.memberIds.size());
            for (DBLocalId memberId : result.memberIds) {
                if (!isAssigned(assignedCluster, memberId)) {
                    validMemberIds.push_back(memberId);
                }
            }

            if (validMemberIds.size() <= 1) {
                continue;
            }

            for (DBLocalId memberId : validMemberIds) {
                storeAssignedCluster(assignedCluster, memberId, result.representativeId);
            }
        }

        if (lastRound) {
            break;
        }
    }
}

int doAlign2clust(Parameters &par, DBWriter &resultWriter, DBReader<DBKeyType> &alnDbr, DBWriter *alnWriter) {
    DBReader<DBKeyType> *seqDbr = new DBReader<DBKeyType>(
        par.db1.c_str(), par.db1Index.c_str(), par.threads, 
        DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
    );
    seqDbr->setIoCacheAdvice(true);
    seqDbr->open(DBReader<DBKeyType>::NOSORT);

    const size_t ioMemoryLimit = availableBufferBudget(par.splitMemoryLimit);
    const size_t seqDataSize = seqDbr->getDataSize();
    const size_t alnDataSize = alnDbr.getDataSize();
    const bool cacheStarved = seqDataSize > ioMemoryLimit
                              || alnDataSize > ioMemoryLimit - seqDataSize;
    DBReader<DBKeyType> *cluDbr = nullptr;
    DBReader<DBKeyType> *cluSeqDbr = nullptr;
    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false) {
        std::string cluIndex = par.filterCluDBFile + ".index";
        cluDbr = new DBReader<DBKeyType>(
            par.filterCluDBFile.c_str(), cluIndex.c_str(), par.threads, 
            DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
        );
        // NOSORT: only ever read by key lookup, so the LINEAR_ACCCESS id mapping is dead weight
        cluDbr->open(DBReader<DBKeyType>::NOSORT);

        std::string cluSeqIndex = par.filterSeqDBFile + ".index";
        cluSeqDbr = new DBReader<DBKeyType>(
            par.filterSeqDBFile.c_str(), cluSeqIndex.c_str(), par.threads, 
            DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
        );
            
        cluSeqDbr->open(DBReader<DBKeyType>::NOSORT);
    } else if (par.filterCluDBFile.empty() != par.filterSeqDBFile.empty()) {
        Debug(Debug::ERROR)<< "Error: Both filterCluDBFile and filterSeqDBFile should be provided together.\n";
        EXIT(EXIT_FAILURE);
    }


    const size_t dbSize = seqDbr->getSize();
    const bool localPrefilterIds = Parameters::isEqualDbtype(
        alnDbr.getDbtype(), Parameters::DBTYPE_PREFILTER_LOCAL_RES);
    if (localPrefilterIds && alnDbr.getSize() != dbSize) {
        Debug(Debug::ERROR) << "Local-ID prefilter has " << alnDbr.getSize()
                            << " entries but its sequence DB has " << dbSize << ".\n";
        EXIT(EXIT_FAILURE);
    }
    BaseMatrix *subMat = new SubstitutionMatrix(
        par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0
    );
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
                                    
    float scorePerColThreshold = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.98);
    Debug(Debug::INFO) << "Score per column threshold for filtering: " << scorePerColThreshold << "\n";
    
    EvalueComputation evaluer(seqDbr->getAminoAcidDBSize(), subMat);
    int32_t xDrop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());
    
    ClusterAssignment *assignedCluster = new(std::nothrow) ClusterAssignment[dbSize];
    Util::checkAllocation(assignedCluster, "Can not allocate assignedCluster memory in Align2Clust");
    // parallel first touch, so the pages are not all faulted onto one node by one thread
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dbSize; ++i) {
        storeAssignedCluster(assignedCluster, i, DB_LOCAL_ID_INVALID);
    }

    int mode = par.clusteringMode;
    
    if (mode == Parameters::SET_COVER) {
        clusterThreadFunc = clusterThreadFuncSetcover;
        Debug(Debug::INFO) << "Using SET_COVER clustering mode\n";
    } else if (mode == Parameters::GREEDY || mode == Parameters::GREEDY_MEM) {
        clusterThreadFunc = clusterThreadFuncGreedy;
        Debug(Debug::INFO) << "Using GREEDY clustering mode\n";
    } else {
        Debug(Debug::ERROR) << "MMseqs2 align2clust doesn't support clustering mode: " << mode << "\n";
        delete[] assignedCluster;
        delete[] fastMatrix.matrix;
        delete[] fastMatrix.matrixData;
        delete subMat;
        seqDbr->close();
        delete seqDbr;
        return EXIT_FAILURE;
    }

    const size_t flagWordCount = assignedFlagWordCount(dbSize);
    assignedFlags = new(std::nothrow) std::atomic<uint64_t>[flagWordCount];
    Util::checkAllocation(assignedFlags, "Can not allocate assignedFlags memory in Align2Clust");
    unsigned int *entryLengths = nullptr;
    if (mode == Parameters::SET_COVER) {
        entryLengths = new(std::nothrow) unsigned int[dbSize];
        Util::checkAllocation(entryLengths, "Can not allocate entryLengths memory in Align2Clust");
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            entryLengths[i] = static_cast<unsigned int>(seqDbr->getEntryLen(i));
        }
        setCoverLengths = entryLengths;
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < flagWordCount; ++i) {
        assignedFlags[i].store(0, std::memory_order_relaxed);
    }

    bool identityOrder = false;
    if (mode != Parameters::SET_COVER) {
        Timer orderTimer;
        identityOrder = seqDbr->isSortedByEntryLengthDescending(par.threads);
        Debug(Debug::INFO) << "Sequence DB length order: "
                           << (identityOrder ? "descending, visit order is the local id" : "unsorted")
                           << " (" << orderTimer.lap() << ")\n";
    }

    // the sequence db is scattered either way, but a length ordered greedy run reads the prefilter in id order, where readahead wins
    const bool seqOnMapping = cacheStarved
                              ? seqDbr->useDescriptorIo(true) == false : true;
    if (cacheStarved && !identityOrder) {
        alnDbr.useDescriptorIo(true);
    }
    if (identityOrder) {
        alnDbr.setSequentialAdvice();
    }

    Timer timer;
    timer.reset();
    PrefInfo *prefRepSizePair = nullptr;
    DBLocalId *lengthOrder = nullptr;

    if (mode == Parameters::SET_COVER) {
        prefRepSizePair = new(std::nothrow) PrefInfo[dbSize];
        Util::checkAllocation(prefRepSizePair, "Can not allocate prefRepSizePair memory in ClusteringAlgorithms::execute");

#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 1000)
            for (size_t i = 0; i < seqDbr->getSize(); i++) {
                size_t alnId = i;
                if (!localPrefilterIds) {
                    const DBKeyType clusterId = seqDbr->getDbKey(i);
                    alnId = requireId(alnDbr.getId(clusterId), "Alignment DB", clusterId);
                }
                const char *data = alnDbr.getData(alnId, thread_idx);
                const size_t dataSize = alnDbr.getEntryLen(alnId);
                prefRepSizePair[i].id = i;
                prefRepSizePair[i].size = (*data == '\0') ? 1 : Util::countLines(data, dataSize);
            }
        }
        SORT_PARALLEL(prefRepSizePair, prefRepSizePair + dbSize, PrefInfo::compareBySizeAndId);
    } else if (identityOrder == false) {
        // NOSORT drops the reader's length order, so materialise it from a fused (length, id) array
        std::pair<unsigned int, DBLocalId> *sortForLength =
            new(std::nothrow) std::pair<unsigned int, DBLocalId>[dbSize];
        Util::checkAllocation(sortForLength, "Can not allocate sortForLength memory in Align2Clust");
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            sortForLength[i] = std::make_pair(static_cast<unsigned int>(seqDbr->getEntryLen(i)),
                                              static_cast<DBLocalId>(i));
        }
        // byte-for-byte DBReader::comparePairBySeqLength: entry length descending, local id ascending
        SORT_PARALLEL(sortForLength, sortForLength + dbSize,
                      [](const std::pair<unsigned int, DBLocalId> &first,
                         const std::pair<unsigned int, DBLocalId> &second) {
                          if (first.first != second.first) {
                              return first.first > second.first;
                          }
                          return first.second < second.second;
                      });
        lengthOrder = new(std::nothrow) DBLocalId[dbSize];
        Util::checkAllocation(lengthOrder, "Can not allocate lengthOrder memory in Align2Clust");
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            lengthOrder[i] = sortForLength[i].second;
        }
        delete[] sortForLength;
    }

    timer.reset();

    size_t endRange = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();
    // the visit index addresses the sequence db directly, so replace the bounds check getDbKey(i) did
    if (mode != Parameters::SET_COVER && endRange > dbSize) {
        Debug(Debug::ERROR) << "Alignment DB has " << endRange << " entries but the sequence DB has " << dbSize << "\n";
        EXIT(EXIT_FAILURE);
    }

    // built after the ordering peak, and the order stays resident, so reserve it against the cache budget too
    const size_t align2clustResultCount = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();
    const size_t reorderCapacityChosen = computeReorderCapacity(par, dbSize, mode, align2clustResultCount,
                                                                lengthOrder != nullptr);
    {
        std::lock_guard<std::mutex> lock(clusterMutex);
        reorderCapacity = reorderCapacityChosen;
        reorderSlots.clear();
        reorderSlots.resize(reorderCapacity);
        reorderFilled.assign(reorderCapacity, 0);
        reorderBufferedCount = 0;
        setCoverCandidates.clear();
        setCoverCandidates.shrink_to_fit();
        setCoverMemberPool.clear();
        setCoverMemberPool.shrink_to_fit();
        setCoverLiveMemberCount = 0;
        currentProcessPosition = 0;
        currentPrefSize = 0;
        allCalculationsDone = false;
    }

    std::thread clusterThread(clusterThreadFunc, assignedCluster);

    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
    bool dropAlnCache = cacheStarved && !identityOrder;
    const size_t alnEpoch = dropAlnCache
        ? std::max<size_t>(ioMemoryLimit / 4 / Util::getPageSize(), 1)
        : std::max<size_t>(endRange, 1);
    Debug::Progress progress(endRange);
    size_t db_maxseqlen = (cluSeqDbr != nullptr)
        ? std::max(seqDbr->getMaxSeqLen(), cluSeqDbr->getMaxSeqLen())
        : seqDbr->getMaxSeqLen();
    // queries are drawn in increasing order so the two held results cannot deadlock
    std::atomic<size_t> nextQuery(0);
#pragma omp parallel
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = (unsigned int) omp_get_thread_num();
#endif
        Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, db_maxseqlen, subMat, &evaluer, 
                       par.compBiasCorrection, par.compBiasCorrectionScale, 
                       par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 
                       0.0, par.zdrop);
        Sequence query(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        Sequence target(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        Sequence element(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        BlockAligner blockAligner(Parameters::DBTYPE_AMINO_ACIDS, db_maxseqlen, subMat, &fastMatrix, 
                                 &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, 
                                 -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
        struct TargetHit {
            size_t id;
            uint32_t slot;
            unsigned short diagonal;
        };
        // one parsed query per lane, kept while its member reads are in flight
        struct QueryWork {
            ClusterResult result;
            size_t representativeId;
            DBKeyType queryKey;
            int queryLength;
            std::vector<TargetHit> targets;
            std::vector<size_t> batchIds;
            size_t loaded;
            bool pending;
            bool skipped;
        };
        QueryWork works[DBReader<DBKeyType>::BATCH_LANES];
        for (unsigned int lane = 0; lane < DBReader<DBKeyType>::BATCH_LANES; lane++) {
            works[lane].targets.reserve(1000);
            works[lane].pending = false;
        }

        const bool includeAlignFiles = (alnWriter != nullptr);
        const bool needTargetKey = includeAlignFiles || cluDbr != nullptr;
        std::string queryCopy;
        std::string alnResultBuffer;
        // Staged member alignments; flushed only if the allpass-gate fully passes.
        std::string pendingMemberAln;
        std::vector<char> alnLineBuffer;
        if (includeAlignFiles) {
            alnLineBuffer.resize(1024 + 32768 * 4);
        }

        // parses one query and puts its member reads in flight
        auto prepare = [&](size_t i, QueryWork &work, unsigned int lane) {
            progress.updateProgress();
            ClusterResult &clusterResult = work.result;
            clusterResult.sequenceIdx = i;
            clusterResult.memberIds.clear();
            std::vector<TargetHit> &targetsWithDiagonal = work.targets;
            std::vector<size_t> &batchIds = work.batchIds;
            targetsWithDiagonal.clear();
            work.loaded = 0;
            work.pending = true;

            size_t representativeId;
            const bool needQueryKey = !localPrefilterIds || includeAlignFiles;
            DBKeyType queryKey = DB_KEY_INVALID;

            if (mode == Parameters::SET_COVER) {
                representativeId = prefRepSizePair[i].id;
                if (needQueryKey) {
                    queryKey = seqDbr->getDbKey(representativeId);
                }
                clusterResult.prefSize = prefRepSizePair[i].size;   // precomputed in the prefix pass
            } else { // GREEDY || GREEDY_MEM
                representativeId = identityOrder ? i : lengthOrder[i];
                if (needQueryKey) {
                    queryKey = seqDbr->getDbKey(representativeId);
                }
                clusterResult.prefSize = 0;                         // greedy has no currentPrefSize gate
            }
            clusterResult.representativeId = representativeId;
            work.representativeId = representativeId;
            work.queryKey = queryKey;

            // Representative already assigned to another cluster: this cluster is discarded
            // by the cluster thread anyway, so skip parsing and aligning it entirely. prefSize
            // is already set (precomputed for set-cover), so the currentPrefSize gate stays
            // correct.
            work.skipped = isAssigned(assignedCluster, representativeId);
            if (work.skipped) {
                return;
            }

            const size_t alignmentId = localPrefilterIds ? representativeId
                : requireId(alnDbr.getId(queryKey), "Alignment DB", queryKey);
            char *alignmentData = alnDbr.getData(alignmentId, threadIdx);
            const size_t queryId = representativeId;
            // index-only, so it is known without faulting the body in; Sequence::mapSequence sets L to it
            const int queryLength = static_cast<int>(seqDbr->getSeqLen(queryId));
            work.queryLength = queryLength;

            batchIds.clear();
            batchIds.push_back(queryId);

            size_t prefSize = 0;
            while (*alignmentData != '\0') {
                hit_t hit = QueryMatcher::parsePrefilterHit(alignmentData);
                const size_t targetId = localPrefilterIds ? static_cast<size_t>(hit.seqId)
                    : requireId(seqDbr->getId(hit.seqId), "Sequence DB", hit.seqId);
                if (localPrefilterIds && targetId >= dbSize) {
                    Debug(Debug::ERROR) << "Local-ID prefilter target " << hit.seqId
                                        << " is outside sequence DB size " << dbSize << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                const bool targetUnassigned = !isAssigned(assignedCluster, targetId);

                if (mode == Parameters::SET_COVER || targetUnassigned) {
                    uint32_t slot = std::numeric_limits<uint32_t>::max();
                    if (targetId != queryId && targetUnassigned &&
                        Util::canBeCovered(par.covThr, par.covMode, queryLength,
                                           seqDbr->getSeqLen(targetId))) {
                        if (batchIds.size() >= std::numeric_limits<uint32_t>::max()) {
                            Debug(Debug::ERROR) << "Too many batched targets for one query.\n";
                            EXIT(EXIT_FAILURE);
                        }
                        slot = static_cast<uint32_t>(batchIds.size());
                        batchIds.push_back(targetId);
                    }
                    targetsWithDiagonal.push_back({targetId, slot, hit.diagonal});
                }
                alignmentData = Util::skipLine(alignmentData);
                prefSize++;
            }
            clusterResult.prefSize = prefSize;
            // the identity hit needs no bytes, so a query without a candidate reads nothing
            if (batchIds.size() > 1) {
                work.loaded = seqDbr->startBatch(batchIds.data(), batchIds.size(), threadIdx, lane);
            }
        };

        auto align = [&](QueryWork &work, unsigned int lane) {
            if (work.pending == false) {
                return;
            }
            ClusterResult &clusterResult = work.result;
            // only the aligner pushes, so a thread never waits on the consumer while it holds another result
            if (work.skipped) {
                pushClusterResult(std::move(clusterResult));
                work.pending = false;
                return;
            }
            const size_t i = clusterResult.sequenceIdx;
            const size_t representativeId = work.representativeId;
            const size_t queryId = representativeId;
            const DBKeyType queryKey = work.queryKey;
            const int queryLength = work.queryLength;
            std::vector<TargetHit> &targetsWithDiagonal = work.targets;
            std::vector<size_t> &batchIds = work.batchIds;
            if (includeAlignFiles) {
                alnResultBuffer.clear();
            }
            // the body is mapped below, on the first target that actually reaches alignment
            const char *querySequence = nullptr;

            // [windowStart, windowEnd) indexes batchIds and describes the current batch arena.
            size_t windowStart = 0;
            size_t windowEnd = 0;

            for (size_t targetIdx = 0; targetIdx < targetsWithDiagonal.size(); targetIdx++) {
                // Representative assigned meanwhile: the cluster thread discards clusters
                // whose representative is assigned, so this result (partial or empty) is
                // never used; stop aligning the rest. Safe in set-cover too: prefSize is
                // already fully counted, so this is just the k-th-target form of the
                // (k=0) rep-skip above.
                if (isAssigned(assignedCluster, representativeId)) {
                    break;
                }

                const size_t targetId = targetsWithDiagonal[targetIdx].id;
                const unsigned short diagonal = targetsWithDiagonal[targetIdx].diagonal;

                const bool isIdentity = (targetId == queryId);
                if (isIdentity) {
                    clusterResult.memberIds.push_back(queryId);
                    if (includeAlignFiles) {
                        // mmseqs forces coverage and seqId to 1.0 for an identity hit (Alignment.cpp)
                        std::string backtrace = par.addBacktrace ? std::string(queryLength, 'M') : std::string();
                        Matcher::result_t selfResult(queryKey, queryLength, 1.0f, 1.0f, 1.0f, 0.0,
                            queryLength, 0, queryLength - 1, queryLength, 0, queryLength - 1, queryLength, backtrace);
                        appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), selfResult, par.addBacktrace);
                    }
                    continue;
                }

                // Skip the (expensive) alignment if the target was assigned meanwhile.
                // Safe in set-cover too: an assigned target is monotonic, so it would be
                // dropped by the cluster thread's re-evaluation anyway.
                if (isAssigned(assignedCluster, targetId)) {
                    continue;
                }

                const uint32_t slot = targetsWithDiagonal[targetIdx].slot;
                if (slot == std::numeric_limits<uint32_t>::max()) {
                    continue;
                }
                const size_t targetLength = seqDbr->getSeqLen(targetId);

                // batchIds[0] is the query, so copy it out before the next batch reuses the arena
                if (querySequence == nullptr) {
                    seqDbr->awaitBatch(threadIdx, lane);
                    windowStart = 0;
                    windowEnd = work.loaded;
                    if (windowEnd == 0) {
                        Debug(Debug::ERROR) << "Failed to batch-load query " << queryKey << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    queryCopy.assign(seqDbr->batchAt(threadIdx, lane, 0), queryLength);
                    querySequence = queryCopy.c_str();
                    query.mapSequence(queryId, queryKey, querySequence, queryLength);
                    blockAligner.initQuery(&query);
                    matcher.initQuery(&query);
                }

                if (slot < windowStart || slot >= windowEnd) {
                    windowStart = slot;
                    windowEnd = slot + seqDbr->startBatch(&batchIds[slot], batchIds.size() - slot, threadIdx, lane);
                    seqDbr->awaitBatch(threadIdx, lane);
                }
                const char *targetSequence = seqDbr->batchAt(threadIdx, lane, slot - windowStart);

                const DBKeyType targetKey = needTargetKey ? seqDbr->getDbKey(targetId) : DB_KEY_INVALID;
                target.mapSequence(targetId, targetKey, targetSequence, targetLength);

                BlockAligner::UngappedAln_res ungappedAlignment = blockAligner.ungappedAlign(&target, diagonal); 
                
                bool hasEvalue = (ungappedAlignment.eval <= par.evalThr);
                bool hasAlnLen = (ungappedAlignment.alnLen >= par.alnLenThr);
                bool hasCoverage = Util::hasCoverage(par.covThr, par.covMode, ungappedAlignment.qcov, ungappedAlignment.tcov);
                float seqId = 0;
                
                if (hasEvalue) {    
                    int identicalCount = 0;
                    for (int q = ungappedAlignment.qStart; q <= ungappedAlignment.qEnd; q++) {
                        char queryLetter = querySequence[q] & static_cast<unsigned char>(~0x20);
                        char targetLetter = targetSequence[ungappedAlignment.tStart + (q - ungappedAlignment.qStart)] & static_cast<unsigned char>(~0x20);
                        identicalCount += (queryLetter == targetLetter) ? 1 : 0;
                    }
                    seqId = Util::computeSeqId(par.seqIdMode, identicalCount, queryLength, target.L, ungappedAlignment.alnLen);
                }
                
                bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                if (isAssigned(assignedCluster, targetId)) continue;

                if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
                    if (isAssigned(assignedCluster, targetId)) continue;
                    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                        // check all the member from filtering file
                        const size_t cluId = requireId(cluDbr->getId(targetKey), "Filter cluster DB", targetKey);
                        char *cluData = cluDbr->getData(cluId, threadIdx);
                        const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                        size_t numClu = Util::countLines(cluData, cluDataSize);
                        bool allpass = true;
                        char buffer[1024];
                        if (includeAlignFiles) {
                            pendingMemberAln.clear();
                        }
                        if (numClu > 1) {
                            while (*cluData != '\0') {
                                Util::parseKey(cluData, buffer);

                                const DBKeyType elementKey = Util::fast_atoi<DBKeyType>(buffer);
                                if (elementKey == targetKey) {
                                    cluData = Util::skipLine(cluData);
                                    continue;
                                }
                                const size_t elementId = requireId(cluSeqDbr->getId(elementKey), "Filter sequence DB", elementKey);
                                const size_t elementLength = cluSeqDbr->getSeqLen(elementId);
                                if (Util::canBeCovered(par.covThr, par.covMode, queryLength, elementLength) == false) {
                                    allpass = false;
                                    break;
                                }
                                char *elementSequence = cluSeqDbr->getData(elementId, threadIdx);
                                short elementDiagonal = diagonal;

                                // 1. ungapped alignment
                                element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                BlockAligner::UngappedAln_res elementUngappedAlignment = blockAligner.ungappedAlign(&element, elementDiagonal);
                                
                                // 2. check the criteria
                                bool elementHasEvalue = (elementUngappedAlignment.eval <= par.evalThr);
                                bool elementHasAlnLen = (elementUngappedAlignment.alnLen >= par.alnLenThr);
                                bool elementHasCoverage = Util::hasCoverage(par.covThr, par.covMode, elementUngappedAlignment.qcov, elementUngappedAlignment.tcov);
                                int elementIdenticalCount = 0;
                                for (int q = elementUngappedAlignment.qStart; q <= elementUngappedAlignment.qEnd; q++) {
                                    char queryLetter = querySequence[q] & static_cast<unsigned char>(~0x20);
                                    char elementLetter = elementSequence[elementUngappedAlignment.tStart + (q - elementUngappedAlignment.qStart)] & static_cast<unsigned char>(~0x20);
                                    elementIdenticalCount += (queryLetter == elementLetter) ? 1 : 0;
                                }
                                float elementSeqId = Util::computeSeqId(par.seqIdMode, elementIdenticalCount, queryLength, elementLength, elementUngappedAlignment.alnLen);
                                bool elementHasSeqId = elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                                
                                if (!(elementHasAlnLen && elementHasCoverage && elementHasSeqId && elementHasEvalue)) {
                                    // 3. gapped alignment
                                    Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(elementDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                                                        swMode, par.seqIdMode, false);
                                    if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                        allpass = false;
                                        break;
                                    }
                                    // stage member alignment (flushed only if allpass holds)
                                    if (includeAlignFiles) {
                                        appendAlignmentResult(pendingMemberAln, alnLineBuffer.data(), res_element, par.addBacktrace);
                                    }
                                } else if (includeAlignFiles) {
                                    // member passed ungapped: gap-free (all 'M') record
                                    std::string elementBacktrace = par.addBacktrace ? std::string(elementUngappedAlignment.alnLen, 'M') : std::string();
                                    Matcher::result_t elementResult(elementKey, elementUngappedAlignment.score, elementUngappedAlignment.qcov,
                                        elementUngappedAlignment.tcov, elementSeqId, elementUngappedAlignment.eval, elementUngappedAlignment.alnLen,
                                        elementUngappedAlignment.qStart, elementUngappedAlignment.qEnd, queryLength,
                                        elementUngappedAlignment.tStart, elementUngappedAlignment.tEnd, elementLength, elementBacktrace);
                                    appendAlignmentResult(pendingMemberAln, alnLineBuffer.data(), elementResult, par.addBacktrace);
                                }
                                cluData = Util::skipLine(cluData);
                            }
                        }
                        if (allpass == false) {
                            continue;
                        }
                    }
                    if (includeAlignFiles) {
                        std::string backtrace = par.addBacktrace ? std::string(ungappedAlignment.alnLen, 'M') : std::string();
                        Matcher::result_t ungappedResult(targetKey, ungappedAlignment.score, ungappedAlignment.qcov,
                            ungappedAlignment.tcov, seqId, ungappedAlignment.eval, ungappedAlignment.alnLen,
                            ungappedAlignment.qStart, ungappedAlignment.qEnd, queryLength,
                            ungappedAlignment.tStart, ungappedAlignment.tEnd, targetLength, backtrace);
                        appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), ungappedResult, par.addBacktrace);
                        // flush staged member alignments (empty unless filter gate ran)
                        alnResultBuffer += pendingMemberAln;
                    }
                    clusterResult.memberIds.push_back(targetId);
                    continue;
                }

                float currentScorePerCol = static_cast<float>(ungappedAlignment.score) / static_cast<float>(ungappedAlignment.diagonalLen);
                if (currentScorePerCol < scorePerColThreshold) {
                    continue;
                }
                
                int alignmentLength = ungappedAlignment.alnLen;
                int queryStartPos = ungappedAlignment.qStart;
                int targetStartPos = ungappedAlignment.tStart;
                int newQueryStartPos = queryStartPos;
                int newTargetStartPos = targetStartPos;
                
                if (queryStartPos == -1 || targetStartPos == -1 || alignmentLength < 3) {
                    continue;
                }

                if (isAssigned(assignedCluster, targetId)) continue;

                bool foundConsecutiveMatchSeed = false;
                for (int blockIdx = 0; blockIdx <= alignmentLength - 3; ++blockIdx) {
                    int queryPos = queryStartPos + blockIdx;
                    int targetPos = targetStartPos + blockIdx;
                    
                    if (querySequence[queryPos] == targetSequence[targetPos] &&
                        querySequence[queryPos + 1] == targetSequence[targetPos + 1] &&
                        querySequence[queryPos + 2] == targetSequence[targetPos + 2]) {
                        newQueryStartPos = queryPos + 1; 
                        newTargetStartPos = targetPos + 1;
                        foundConsecutiveMatchSeed = true;
                        break;
                    }
                }
                
                if (foundConsecutiveMatchSeed) {
                    std::string gappedBacktrace;

                    s_align gappedAlignment = blockAligner.bandedalign(&target, newQueryStartPos, newTargetStartPos,
                                                                       gappedBacktrace, xDrop, par.covThr, par.covMode);
                    // bandedalign signals failure with evalue < 0 and an empty backtrace, which computeSeqId would divide by
                    if (gappedAlignment.evalue < 0 || gappedBacktrace.empty()) {
                        continue;
                    }
                    unsigned int gappedAlnLength = gappedBacktrace.size();
                    double gappedSeqId = Util::computeSeqId(par.seqIdMode, gappedAlignment.identicalAACnt,
                                                           queryLength, targetLength, gappedAlnLength);
                    Matcher::result_t result = Matcher::result_t(
                        targetKey, gappedAlignment.score1, gappedAlignment.qCov, gappedAlignment.tCov, 
                        gappedSeqId, gappedAlignment.evalue, gappedAlnLength,
                        gappedAlignment.qStartPos1, gappedAlignment.qEndPos1, queryLength,
                        gappedAlignment.dbStartPos1, gappedAlignment.dbEndPos1, targetLength, gappedBacktrace
                    );
                    if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, 
                                                par.alnLenThr, par.covMode, par.covThr)) {
                        if (isAssigned(assignedCluster, targetId)) continue;
                        if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                            // check all the member from filtering file
                            const size_t cluId = requireId(cluDbr->getId(targetKey), "Filter cluster DB", targetKey);
                            char *cluData = cluDbr->getData(cluId, threadIdx);
                            const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                            size_t numClu = Util::countLines(cluData, cluDataSize);
                            bool allpass = true;
                            char buffer[1024];
                            if (includeAlignFiles) {
                                pendingMemberAln.clear();
                            }
                            if (numClu > 1) {
                                while (*cluData != '\0') {
                                    Util::parseKey(cluData, buffer);
                                    const DBKeyType elementKey = Util::fast_atoi<DBKeyType>(buffer);
                                    if (elementKey == targetKey) {
                                        cluData = Util::skipLine(cluData);
                                        continue;
                                    }
                                    const size_t elementId = requireId(cluSeqDbr->getId(elementKey), "Filter sequence DB", elementKey);
                                    const size_t elementLength = cluSeqDbr->getSeqLen(elementId);
                                    if (Util::canBeCovered(par.covThr, par.covMode, queryLength, elementLength) == false) {
                                        allpass = false;
                                        break;
                                    }
                                    char *elementSequence = cluSeqDbr->getData(elementId, threadIdx);
                                    short elementDiagonal = 0;

                                    // 1. ungapped alignment
                                    element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                    BlockAligner::UngappedAln_res elementUngappedAlignment = blockAligner.ungappedAlign(&element, elementDiagonal);
                                    
                                    // 2. check the criteria
                                    bool elementHasEvalue = (elementUngappedAlignment.eval <= par.evalThr);
                                    bool elementHasAlnLen = (elementUngappedAlignment.alnLen >= par.alnLenThr);
                                    bool elementHasCoverage = Util::hasCoverage(par.covThr, par.covMode, elementUngappedAlignment.qcov, elementUngappedAlignment.tcov);
                                    int elementIdenticalCount = 0;
                                    for (int q = elementUngappedAlignment.qStart; q <= elementUngappedAlignment.qEnd; q++) {
                                        char queryLetter = querySequence[q] & static_cast<unsigned char>(~0x20);
                                        char elementLetter = elementSequence[elementUngappedAlignment.tStart + (q - elementUngappedAlignment.qStart)] & static_cast<unsigned char>(~0x20);
                                        elementIdenticalCount += (queryLetter == elementLetter) ? 1 : 0;
                                    }
                                    float elementSeqId = Util::computeSeqId(par.seqIdMode, elementIdenticalCount, queryLength, elementLength, elementUngappedAlignment.alnLen);
                                    bool elementHasSeqId = elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                                    
                                    if (!(elementHasAlnLen && elementHasCoverage && elementHasSeqId && elementHasEvalue)) {
                                        // 3. gapped alignment
                                        Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(elementDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                                                            swMode, par.seqIdMode, false);
                                        if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                            allpass = false;
                                            break;
                                        }
                                        // stage member alignment (flushed only if allpass holds)
                                        if (includeAlignFiles) {
                                            appendAlignmentResult(pendingMemberAln, alnLineBuffer.data(), res_element, par.addBacktrace);
                                        }
                                    } else if (includeAlignFiles) {
                                        // member passed ungapped: gap-free (all 'M') record
                                        std::string elementBacktrace = par.addBacktrace ? std::string(elementUngappedAlignment.alnLen, 'M') : std::string();
                                        Matcher::result_t elementResult(elementKey, elementUngappedAlignment.score, elementUngappedAlignment.qcov,
                                            elementUngappedAlignment.tcov, elementSeqId, elementUngappedAlignment.eval, elementUngappedAlignment.alnLen,
                                            elementUngappedAlignment.qStart, elementUngappedAlignment.qEnd, queryLength,
                                            elementUngappedAlignment.tStart, elementUngappedAlignment.tEnd, elementLength, elementBacktrace);
                                        appendAlignmentResult(pendingMemberAln, alnLineBuffer.data(), elementResult, par.addBacktrace);
                                    }
                                    cluData = Util::skipLine(cluData);
                                }
                            }
                            if (allpass == false) {
                                continue;
                            }
                        }
                        if (includeAlignFiles) {
                            appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), result, par.addBacktrace);
                            // flush staged member alignments (empty unless filter gate ran)
                            alnResultBuffer += pendingMemberAln;
                        }
                        clusterResult.memberIds.push_back(targetId);
                    }
                }
            }

            if (includeAlignFiles) {
                alnWriter->writeData(alnResultBuffer.c_str(), alnResultBuffer.length(), queryKey, threadIdx);
            }
            pushClusterResult(std::move(clusterResult));
            work.pending = false;

            if (identityOrder && i > reorderCapacity && (i % ALIGN2CLUST_CACHE_DROP_STRIDE) == 0) {
                // a producer here is not blocked, so every entry this far back is already consumed
                const size_t safeEnd = i - reorderCapacity;
                if (localPrefilterIds) {
                    dropBehind(alnDbr, alnCacheDroppedTo, alnDbr.getOffset(safeEnd));
                }
                if (seqOnMapping) {
                    dropBehind(*seqDbr, seqCacheDroppedTo,
                               seqDbr->getOffset(coverageBandStart(seqDbr, safeEnd, par.covThr, par.covMode)));
                }
            }
        };

        for (size_t epochStart = 0; epochStart < endRange || epochStart == 0; epochStart += alnEpoch) {
            const size_t epochEnd = std::min(endRange, epochStart + alnEpoch);
            // a thread submits the reads of the query it just drew, then aligns the one it still holds
            unsigned int lane = 0;
            bool holding = false;
            for (size_t i = nextQuery.fetch_add(1); i < epochEnd; i = nextQuery.fetch_add(1)) {
                prepare(i, works[lane], lane);
                if (holding) {
                    align(works[lane ^ 1u], lane ^ 1u);
                }
                holding = true;
                lane ^= 1u;
            }
            if (holding) {
                align(works[lane ^ 1u], lane ^ 1u);
            }
#pragma omp barrier
#pragma omp single
            {
                if (lengthOrder != nullptr) {
                    discardLengthOrderPrefix(lengthOrder, epochStart, epochEnd);
                }
                if (dropAlnCache) {
                    alnDbr.dropCacheAll();
                }
                nextQuery.store(epochEnd);
            }
        }
    }

    // nothing past the producer loop reads the visit order, so drop it before memberOrder is sized
    delete[] lengthOrder;
    // no sequence body is read past here either; writeClustering only needs the index
    seqDbr->dropCacheAll();
    seqDbr->unmapData();

    {
        std::lock_guard<std::mutex> lock(clusterMutex);
        allCalculationsDone = true;
    }
    clusterCondition.notify_one();
    reorderSpaceCondition.notify_all();
    
    if (clusterThread.joinable()) {
        clusterThread.join();
    }

    // ring and set-cover pool are dead once the consumer joins; free them before memberOrder is sized
    std::vector<ClusterResult>().swap(reorderSlots);
    std::vector<unsigned char>().swap(reorderFilled);
    std::vector<SetCoverCandidate>().swap(setCoverCandidates);
    std::vector<DBLocalId>().swap(setCoverMemberPool);

    // neither is read past the producer loop, so the output phase does not compete with them for memory
    delete[] prefRepSizePair;
    setCoverLengths = nullptr;
    delete[] entryLengths;
    alnDbr.close();

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dbSize; ++i) {
        if (!isAssigned(assignedCluster, i)) {
            assignedCluster[i].store(i, std::memory_order_relaxed);
        }
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < assignedFlagWordCount(dbSize); ++i) {
        assignedFlags[i].store(UINT64_MAX, std::memory_order_relaxed);
    }

    // group members by representative through a local-id permutation: 8 byte per sequence, not 16
    DBLocalId *memberOrder = new(std::nothrow) DBLocalId[dbSize];
    Util::checkAllocation(memberOrder, "Can not allocate memberOrder memory in Align2Clust");
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dbSize; i++) {
        memberOrder[i] = static_cast<DBLocalId>(i);
    }

    // comparator touches only the two arrays, so no getDbKey call per comparison
    SORT_PARALLEL(memberOrder, memberOrder + dbSize, [assignedCluster](DBLocalId first, DBLocalId second) {
        const DBLocalId firstRep = loadAssignedCluster(assignedCluster, first);
        const DBLocalId secondRep = loadAssignedCluster(assignedCluster, second);
        if (firstRep != secondRep) {
            return firstRep < secondRep;
        }
        return first < second;
    });

    size_t clusterCount = (dbSize > 0) ? 1 : 0;
#pragma omp parallel for schedule(static) reduction(+:clusterCount)
    for (size_t i = 1; i < dbSize; i++) {
        clusterCount += (loadAssignedCluster(assignedCluster, memberOrder[i])
                         != loadAssignedCluster(assignedCluster, memberOrder[i - 1]));
    }

    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << clusterCount << "\n";

    Timer writeTimer;
    writeClustering(&resultWriter, seqDbr, assignedCluster, memberOrder, dbSize, par.threads);
    Debug(Debug::INFO) << "Time for writing clusters: " << writeTimer.lap() << "\n";

    delete[] memberOrder;
    delete[] assignedFlags;
    assignedFlags = nullptr;
    delete[] assignedCluster;
    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
    delete subMat;
    seqDbr->close();
    delete seqDbr;

    if (cluDbr != nullptr) {
        cluDbr->close();
        delete cluDbr;
    }
    if (cluSeqDbr != nullptr) {
        cluSeqDbr->close();
        delete cluSeqDbr;
    }
    return 0;
}

int align2clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    
    Timer timer;
    timer.reset();
    
    DBReader<DBKeyType> alnDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    alnDbr.setIoCacheAdvice(true);
    alnDbr.open(DBReader<DBKeyType>::NOSORT);
    int dbtype =  Parameters::DBTYPE_CLUSTER_RES;

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, dbtype);
    resultWriter.open();

    // Optional alignment-result output; path derived from the cluster DB (db3 + "_aln").
    // Alignment output needs a backtrace (-a) and score+cov+seqid so member records carry a CIGAR.
    if (par.includeAlignFiles) {
        const unsigned int effectiveSwMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
        if (par.addBacktrace == false || effectiveSwMode != Matcher::SCORE_COV_SEQID) {
            Debug(Debug::ERROR) << "Writing alignment files requires backtrace and score+cov+seqid alignment.\n"
                                << "Please re-run with '-a 1' and '--alignment-mode "
                                << Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID << "'.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    DBWriter *alnWriter = nullptr;
    if (par.includeAlignFiles) {
        std::string alnDb = par.db3 + "_aln";
        std::string alnDbIndex = alnDb + ".index";
        alnWriter = new DBWriter(alnDb.c_str(), alnDbIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
        alnWriter->open();
    }

    int status = doAlign2clust(par, resultWriter, alnDbr, alnWriter);

    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";

    // memberOrder is sorted by representative and NOSORT keeps getDbKey ascending, so the entries leave key sorted
    resultWriter.close(false, false);
    if (alnWriter != nullptr) {
        alnWriter->close();
        delete alnWriter;
    }

    return status;
}
