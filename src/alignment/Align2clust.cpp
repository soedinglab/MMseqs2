#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "IndexReader.h"
#include "FastSort.h"
#include "BlockAligner.h"
#include "Alignment.h"
#include "AlignmentSymmetry.h"
#include <atomic>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

#define MAX_SIZE 4096 //change
#define MIN_SIZE 32

struct ClusterResult {
    size_t sequenceIdx;
    DBLocalId representativeId;
    size_t prefSize;
    std::vector<DBLocalId> memberIds;
};

struct PrefInfo {
    DBLocalId id;
    size_t size;

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

struct SetCoverComparator {
    bool operator()(const SetCoverCandidate& a, const SetCoverCandidate& b) const {
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

static std::mutex clusterMutex;
static std::condition_variable clusterCondition;
static std::condition_variable reorderSpaceCondition;
// In-order reorder buffer. Worker results carry a dense, contiguous sequenceIdx in
// [0, endRange); instead of a heap we use a fixed-size ring indexed by
// sequenceIdx % reorderCapacity. Publishing a result is O(1) under the lock, and the
// slot for the next-in-order result is always free -> deadlock-free throttling.
static std::vector<ClusterResult> reorderSlots;    // slot storage, size == reorderCapacity
static std::vector<unsigned char> reorderFilled;   // 1 == slot holds an unconsumed result
static size_t reorderCapacity = 0;                 // ring size == max out-of-order window
static size_t reorderBufferedCount = 0;            // number of filled (unconsumed) slots
// Binary heap (std::push_heap / pop_heap) of lightweight candidates. A plain
// vector is used instead of std::priority_queue so the member pool can be
// compacted by iterating over the live candidates.
static std::vector<SetCoverCandidate> setCoverReadyQueue;
static SetCoverComparator setCoverComparator;
// Shared backing store for all member ids referenced by setCoverReadyQueue.
static std::vector<DBLocalId> setCoverMemberPool;
static size_t setCoverLiveMemberCount = 0;

static size_t currentProcessPosition = 0;
static size_t currentPrefSize = 0;
static bool allCalculationsDone = false;

typedef std::atomic<DBLocalId> ClusterAssignment;

static DBLocalId loadAssignedCluster(const ClusterAssignment *assignedCluster, size_t sequenceId) {
    return assignedCluster[sequenceId].load(std::memory_order_relaxed);
}

static void storeAssignedCluster(ClusterAssignment *assignedCluster, size_t sequenceId, DBLocalId representativeId) {
    assignedCluster[sequenceId].store(representativeId, std::memory_order_relaxed);
}

// Serialize a single alignment record and append it to the per-representative buffer.
static void appendAlignmentResult(std::string &alnResultBuffer, char *lineBuffer, const Matcher::result_t &result, bool addBacktrace) {
    size_t len = Matcher::resultToBuffer(lineBuffer, result, addBacktrace);
    alnResultBuffer.append(lineBuffer, len);
}

// Reclaim dead space in setCoverMemberPool once more than half of it is dead.
// Live members are copied into a fresh pool and the candidates' offsets updated;
// the heap order is preserved because only offsets (not the sort keys) change.
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
    for (SetCoverCandidate &candidate : setCoverReadyQueue) {
        const size_t newOffset = compactedPool.size();
        compactedPool.insert(compactedPool.end(),
                             setCoverMemberPool.begin() + candidate.memberOffset,
                             setCoverMemberPool.begin() + candidate.memberOffset + candidate.memberCount);
        candidate.memberOffset = newOffset;
    }
    setCoverMemberPool.swap(compactedPool);
}

// User override for the reorder-buffer capacity. Returns 0 when unset/invalid, meaning
// "size automatically from the memory budget" (see computeReorderCapacity).
static size_t parseAlign2clustMaxQueuedResults() {
    const char *envValue = getenv("MMSEQS_ALIGN2CLUST_MAX_QUEUE");
    if (envValue == nullptr || *envValue == '\0') {
        return 0;
    }

    char *end = nullptr;
    unsigned long long parsedValue = strtoull(envValue, &end, 10);
    if (end == envValue || *end != '\0' || parsedValue == 0) {
        Debug(Debug::WARNING) << "Ignoring invalid MMSEQS_ALIGN2CLUST_MAX_QUEUE=" << envValue
                              << "; sizing the reorder buffer automatically\n";
        return 0;
    }
    return static_cast<size_t>(parsedValue);
}

// Choose the reorder-buffer capacity from the available memory budget (mirrors the
// kmermatcher OOM-aware split sizing). Reserve the per-sequence arrays that coexist
// with the buffer (assignedCluster, plus prefRepSizePair in set-cover), keep 10%
// headroom, and bound each buffered result by its worst case (a ClusterResult header
// plus up to maxResListLen member ids). MMSEQS_ALIGN2CLUST_MAX_QUEUE, if set, caps it
// further; the capacity never exceeds the number of results produced.
static size_t computeReorderCapacity(const Parameters &par, size_t dbSize, int mode, size_t resultCount) {
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);

    size_t fixedMemory = dbSize * sizeof(ClusterAssignment);
    if (mode == Parameters::SET_COVER) {
        fixedMemory += dbSize * sizeof(PrefInfo);
    }

    const size_t budget = (memoryLimit > fixedMemory)
        ? static_cast<size_t>(static_cast<double>(memoryLimit - fixedMemory) * 0.9)
        : 0;
    const size_t bytesPerResult = sizeof(ClusterResult) + (par.maxResListLen + 1) * sizeof(DBLocalId);

    size_t capacity = std::min(resultCount, std::max<size_t>(1, budget / bytesPerResult));
    const size_t userCap = parseAlign2clustMaxQueuedResults();  // 0 == auto
    if (userCap != 0) {
        capacity = std::min(capacity, userCap);
    }
    capacity = std::max<size_t>(1, capacity);

    Debug(Debug::INFO) << "Reorder buffer sizing: memory limit " << memoryLimit << " byte, reserved (fixed) "
                       << fixedMemory << " byte, budget " << budget << " byte, " << bytesPerResult
                       << " byte/result\n";
    Debug(Debug::INFO) << "Reorder buffer capacity: " << capacity << " results ("
                       << capacity * sizeof(ClusterResult) << " byte pre-allocated slots)\n";
    return capacity;
}

static void pushClusterResult(ClusterResult &&clusterResult) {
    const size_t idx = clusterResult.sequenceIdx;
    bool shouldNotifyClusterThread = false;
    {
        std::unique_lock<std::mutex> lock(clusterMutex);
        // Wait until this result's ring slot is free. Because reorderCapacity >= 1,
        // idx == currentProcessPosition always satisfies the predicate immediately, so
        // the producer of the result the consumer is waiting for never blocks (deadlock-free).
        reorderSpaceCondition.wait(lock, [&] {
            return allCalculationsDone || idx < currentProcessPosition + reorderCapacity;
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

static void writeData(DBWriter *dbWriter, const std::pair<DBKeyType, DBKeyType> * results, size_t dbSize) {
    std::string resultString;
    resultString.reserve(1024 * 1024 * 1024);
    char buffer[32];
    DBKeyType previousRepresentativeKey = DB_KEY_INVALID;
    
    for (size_t i = 0; i < dbSize; i++) {
        DBKeyType currentRepresentativeKey = results[i].first;
        
        if (previousRepresentativeKey != currentRepresentativeKey) {
            if (previousRepresentativeKey != DB_KEY_INVALID) {
                dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey);
            }
            resultString.clear();
            char *outPos = Itoa::u64toa_sse2(static_cast<uint64_t>(currentRepresentativeKey), buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }
        
        DBKeyType memberKey = results[i].second;
        if (memberKey != currentRepresentativeKey) {
            char *outPos = Itoa::u64toa_sse2(static_cast<uint64_t>(memberKey), buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }
        
        previousRepresentativeKey = currentRepresentativeKey;
    }
    
    if (previousRepresentativeKey != DB_KEY_INVALID) {
        dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey);
    }
}

static void (*clusterThreadFunc)(ClusterAssignment*) = nullptr;

void clusterThreadFuncSetcover(ClusterAssignment* assignedCluster) {
    while (true) {
        std::unique_lock<std::mutex> lock(clusterMutex);
        
        clusterCondition.wait(lock, [] {
            return reorderFilled[currentProcessPosition % reorderCapacity] != 0 ||
                   allCalculationsDone;
        });

        // 1) reorder buffer → setCoverReadyQueue
        while (reorderFilled[currentProcessPosition % reorderCapacity] != 0) {
            const size_t slot = currentProcessPosition % reorderCapacity;
            ClusterResult result = std::move(reorderSlots[slot]);
            reorderFilled[slot] = 0;
            reorderBufferedCount--;
            currentProcessPosition++;
            reorderSpaceCondition.notify_all();
            currentPrefSize = result.prefSize;

            if (result.memberIds.size() > 1) {
                SetCoverCandidate candidate;
                candidate.representativeId = result.representativeId;
                candidate.memberCount = result.memberIds.size();
                candidate.memberOffset = setCoverMemberPool.size();
                setCoverMemberPool.insert(setCoverMemberPool.end(),
                                          result.memberIds.begin(), result.memberIds.end());
                setCoverLiveMemberCount += candidate.memberCount;
                setCoverReadyQueue.push_back(candidate);
                std::push_heap(setCoverReadyQueue.begin(), setCoverReadyQueue.end(), setCoverComparator);
            }
        }

        // 2) assign candidates guaranteed to be the currently largest set
        while (setCoverReadyQueue.empty() == false &&
               (allCalculationsDone ||
                setCoverReadyQueue.front().memberCount > currentPrefSize)) {

            std::pop_heap(setCoverReadyQueue.begin(), setCoverReadyQueue.end(), setCoverComparator);
            SetCoverCandidate candidate = setCoverReadyQueue.back();
            setCoverReadyQueue.pop_back();
            setCoverLiveMemberCount -= candidate.memberCount;

            if (loadAssignedCluster(assignedCluster, candidate.representativeId) != DB_LOCAL_ID_INVALID) {
                continue;
            }

            // Drop already-assigned members, compacting the survivors in place
            // within the pool region [memberOffset, memberOffset + memberCount).
            DBLocalId *members = setCoverMemberPool.data() + candidate.memberOffset;
            size_t validCount = 0;
            for (size_t i = 0; i < candidate.memberCount; i++) {
                if (loadAssignedCluster(assignedCluster, members[i]) == DB_LOCAL_ID_INVALID) {
                    members[validCount++] = members[i];
                }
            }

            if (validCount <= 1) {
                continue;
            }

            if (validCount != candidate.memberCount) {
                candidate.memberCount = validCount;
                setCoverLiveMemberCount += validCount;
                setCoverReadyQueue.push_back(candidate);
                std::push_heap(setCoverReadyQueue.begin(), setCoverReadyQueue.end(), setCoverComparator);
                continue;
            }

            for (size_t i = 0; i < candidate.memberCount; i++) {
                storeAssignedCluster(assignedCluster, members[i], candidate.representativeId);
            }
        }

        // Compaction only touches consumer-private structures (setCoverReadyQueue,
        // setCoverMemberPool, setCoverLiveMemberCount), never shared state, so release
        // the mutex during the copy so worker threads can keep pushing results.
        lock.unlock();
        compactSetCoverMemberPool();
        lock.lock();

        if (allCalculationsDone &&
            reorderBufferedCount == 0 &&
            setCoverReadyQueue.empty()) {
            break;
        }
    }
}

void clusterThreadFuncGreedy(ClusterAssignment* assignedCluster) {
    while (true) {
        std::unique_lock<std::mutex> lock(clusterMutex);
        
        clusterCondition.wait(lock, [] {
            return reorderFilled[currentProcessPosition % reorderCapacity] != 0 ||
                   allCalculationsDone;
        });

        if (allCalculationsDone && reorderBufferedCount == 0) {
            break;
        }

        while (reorderFilled[currentProcessPosition % reorderCapacity] != 0) {
            const size_t slot = currentProcessPosition % reorderCapacity;
            ClusterResult result = std::move(reorderSlots[slot]);
            reorderFilled[slot] = 0;
            reorderBufferedCount--;
            currentProcessPosition++;
            reorderSpaceCondition.notify_all();

            if (loadAssignedCluster(assignedCluster, result.representativeId) != DB_LOCAL_ID_INVALID) {
                continue;
            }

            std::vector<DBLocalId> validMemberIds;
            validMemberIds.reserve(result.memberIds.size());
            for (DBLocalId memberId : result.memberIds) {
                if (loadAssignedCluster(assignedCluster, memberId) == DB_LOCAL_ID_INVALID) {
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
    }
}

int doAlign2clust(Parameters &par, DBWriter &resultWriter, DBReader<DBKeyType> &alnDbr, DBWriter *alnWriter) {
    DBReader<DBKeyType> *seqDbr = new DBReader<DBKeyType>(
        par.db1.c_str(), par.db1Index.c_str(), par.threads, 
        DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
    );
    seqDbr->open(DBReader<DBKeyType>::SORT_BY_LENGTH);
 
    DBReader<DBKeyType> *cluDbr = nullptr;
    DBReader<DBKeyType> *cluSeqDbr = nullptr;
    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false) {
        std::string cluIndex = par.filterCluDBFile + ".index";
        cluDbr = new DBReader<DBKeyType>(
            par.filterCluDBFile.c_str(), cluIndex.c_str(), par.threads, 
            DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
        );
        cluDbr->open(DBReader<DBKeyType>::LINEAR_ACCCESS);

        std::string cluSeqIndex = par.filterSeqDBFile + ".index";
        cluSeqDbr = new DBReader<DBKeyType>(
            par.filterSeqDBFile.c_str(), cluSeqIndex.c_str(), par.threads, 
            DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX
        );
            
        cluSeqDbr->open(DBReader<DBKeyType>::LINEAR_ACCCESS);
    } else if (par.filterCluDBFile.empty() != par.filterSeqDBFile.empty()) {
        Debug(Debug::ERROR)<< "Error: Both filterCluDBFile and filterSeqDBFile should be provided together.\n";
        EXIT(EXIT_FAILURE);
    }


    const size_t dbSize = seqDbr->getSize();

    BaseMatrix *subMat = new SubstitutionMatrix(
        par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0
    );
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
                                    
    float scorePerColThreshold = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    Debug(Debug::INFO) << "Score per column threshold for filtering: " << scorePerColThreshold << "\n";
    
    EvalueComputation evaluer(seqDbr->getAminoAcidDBSize(), subMat);
    int32_t xDrop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());
    
    ClusterAssignment *assignedCluster = new(std::nothrow) ClusterAssignment[dbSize];
    Util::checkAllocation(assignedCluster, "Can not allocate assignedCluster memory in Align2Clust");
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

    // Ring size = out-of-order window, sized from the memory budget (OOM-aware) and
    // capped by the result count. sequenceIdx runs over [0, endRange); every index
    // publishes exactly one result.
    const size_t align2clustResultCount = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();
    const size_t reorderCapacityChosen = computeReorderCapacity(par, dbSize, mode, align2clustResultCount);
    {
        std::lock_guard<std::mutex> lock(clusterMutex);
        reorderCapacity = reorderCapacityChosen;
        reorderSlots.clear();
        reorderSlots.resize(reorderCapacity);
        reorderFilled.assign(reorderCapacity, 0);
        reorderBufferedCount = 0;
        setCoverReadyQueue.clear();
        setCoverReadyQueue.shrink_to_fit();
        setCoverMemberPool.clear();
        setCoverMemberPool.shrink_to_fit();
        setCoverLiveMemberCount = 0;
        currentProcessPosition = 0;
        currentPrefSize = 0;
        allCalculationsDone = false;
    }

    std::thread clusterThread(clusterThreadFunc, assignedCluster);
    
    Timer timer;
    timer.reset();
    PrefInfo *prefRepSizePair = nullptr;
    
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
                const DBKeyType clusterId = seqDbr->getDbKey(i);
                const size_t alnId = alnDbr.getId(clusterId);
                const char *data = alnDbr.getData(alnId, thread_idx);
                const size_t dataSize = alnDbr.getEntryLen(alnId);
                prefRepSizePair[i].id = seqDbr->getId(clusterId);
                prefRepSizePair[i].size = (*data == '\0') ? 1 : Util::countLines(data, dataSize);
            }
        }
        SORT_PARALLEL(prefRepSizePair, prefRepSizePair + dbSize, PrefInfo::compareBySizeAndId);
    }

    timer.reset();

    size_t endRange = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();
    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
    Debug::Progress progress(endRange);
    size_t db_maxseqlen = (cluSeqDbr != nullptr)
        ? std::max(seqDbr->getMaxSeqLen(), cluSeqDbr->getMaxSeqLen())
        : seqDbr->getMaxSeqLen();
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
        std::vector<std::pair<DBKeyType, unsigned short>> targetsWithDiagonal;
        targetsWithDiagonal.reserve(1000);

        const bool includeAlignFiles = (alnWriter != nullptr);
        std::string alnResultBuffer;
        std::vector<char> alnLineBuffer;
        if (includeAlignFiles) {
            alnLineBuffer.resize(1024 + 32768 * 4);
        }

#pragma omp for schedule(dynamic, 1) nowait
        for (size_t i = 0; i < endRange; i++) {
            progress.updateProgress();
            ClusterResult clusterResult;
            clusterResult.sequenceIdx = i;
            targetsWithDiagonal.clear();
            if (includeAlignFiles) {
                alnResultBuffer.clear();
            }

            size_t representativeId;
            DBKeyType queryKey;

            if (mode == Parameters::SET_COVER) {
                representativeId = prefRepSizePair[i].id;
                queryKey = seqDbr->getDbKey(representativeId);
                clusterResult.prefSize = prefRepSizePair[i].size;   // precomputed in the prefix pass
            } else { // GREEDY || GREEDY_MEM
                queryKey = seqDbr->getDbKey(i);
                representativeId = seqDbr->getId(queryKey);
                clusterResult.prefSize = 0;                         // greedy has no currentPrefSize gate
            }
            clusterResult.representativeId = representativeId;

            // Representative already assigned to another cluster: this cluster is discarded
            // by the cluster thread anyway, so skip parsing and aligning it entirely. prefSize
            // is already set (precomputed for set-cover), so the currentPrefSize gate stays
            // correct.
            if (loadAssignedCluster(assignedCluster, representativeId) != DB_LOCAL_ID_INVALID) {
                pushClusterResult(std::move(clusterResult));
                continue;
            }

            const size_t alignmentId = alnDbr.getId(queryKey);
            char *alignmentData = alnDbr.getData(alignmentId, threadIdx);
            size_t queryId = representativeId;
            char *querySequence = seqDbr->getData(queryId, threadIdx);
            size_t queryLength = seqDbr->getSeqLen(queryId);
            query.mapSequence(queryId, queryKey, querySequence, queryLength);
            blockAligner.initQuery(&query);
            matcher.initQuery(&query);

            size_t prefSize = 0;
            while (*alignmentData != '\0') {
                hit_t hit = QueryMatcher::parsePrefilterHit(alignmentData);
                if (mode == Parameters::SET_COVER) {
                    targetsWithDiagonal.push_back(std::make_pair(hit.seqId, hit.diagonal));
                } else {
                    const size_t targetId = seqDbr->getId(hit.seqId);
                    if (loadAssignedCluster(assignedCluster, targetId) == DB_LOCAL_ID_INVALID) {
                            targetsWithDiagonal.push_back(std::make_pair(hit.seqId, hit.diagonal));
                    }
                }
                alignmentData = Util::skipLine(alignmentData);
                prefSize++;
            }
            clusterResult.prefSize = prefSize;   // exact parsed count for the aligned path

            for (size_t targetIdx = 0; targetIdx < targetsWithDiagonal.size(); targetIdx++) {
                // Representative assigned meanwhile: the cluster thread discards clusters
                // whose representative is assigned, so this result (partial or empty) is
                // never used; stop aligning the rest. Safe in set-cover too: prefSize is
                // already fully counted, so this is just the k-th-target form of the
                // (k=0) rep-skip above.
                if (loadAssignedCluster(assignedCluster, representativeId) != DB_LOCAL_ID_INVALID) {
                    break;
                }

                const DBKeyType targetKey = targetsWithDiagonal[targetIdx].first;
                const unsigned short diagonal = targetsWithDiagonal[targetIdx].second;
                const size_t targetId = seqDbr->getId(targetKey);

                const bool isIdentity = (queryKey == targetKey);
                if (isIdentity) {
                    clusterResult.memberIds.push_back(queryId);
                    if (includeAlignFiles) {
                        // Identity hit: no alignment needed. mmseqs forces coverage/seqId to
                        // 1.0 for identity (see Alignment.cpp), so emit a full-length self
                        // record directly instead of running Smith-Waterman.
                        std::string backtrace = par.addBacktrace ? std::string(query.L, 'M') : std::string();
                        Matcher::result_t selfResult(queryKey, query.L, 1.0f, 1.0f, 1.0f, 0.0,
                            query.L, 0, query.L - 1, query.L, 0, query.L - 1, query.L, backtrace);
                        appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), selfResult, par.addBacktrace);
                    }
                    continue;
                }

                // Skip the (expensive) alignment if the target was assigned meanwhile.
                // Safe in set-cover too: an assigned target is monotonic, so it would be
                // dropped by the cluster thread's re-evaluation anyway.
                if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) {
                    continue;
                }

                char *targetSequence = seqDbr->getData(targetId, threadIdx);
                size_t targetLength = seqDbr->getSeqLen(targetId);
                target.mapSequence(targetId, targetKey, targetSequence, targetLength);

                if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
                    continue;
                }

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
                    seqId = Util::computeSeqId(par.seqIdMode, identicalCount, query.L, target.L, ungappedAlignment.alnLen);
                }
                
                bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) continue;

                if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
                    if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) continue;
                    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                        // check all the member from filtering file
                        const size_t cluId = cluDbr->getId(targetKey);
                        char *cluData = cluDbr->getData(cluId, threadIdx);
                        const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                        size_t numClu = Util::countLines(cluData, cluDataSize);
                        bool allpass = true;
                        char buffer[1024];
                        if (numClu > 1) { // if not singleton
                            while (*cluData != '\0') {
                                Util::parseKey(cluData, buffer);

                                const DBKeyType elementKey = Util::fast_atoi<DBKeyType>(buffer);
                                if (elementKey == targetKey) {
                                    cluData = Util::skipLine(cluData);
                                    continue;
                                }
                                const size_t elementId = cluSeqDbr->getId(elementKey);
                                char *elementSequence = cluSeqDbr->getData(elementId, threadIdx);
                                size_t elementLength = cluSeqDbr->getSeqLen(elementId);
                                short concatedDiagonal = diagonal;
                                
                                // 1. ungapped alignment
                                element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
                                    allpass = false;
                                    break;
                                }
                                BlockAligner::UngappedAln_res elementUngappedAlignment = blockAligner.ungappedAlign(&element, concatedDiagonal);
                                
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
                                float elementSeqId = Util::computeSeqId(par.seqIdMode, elementIdenticalCount, query.L, elementLength, elementUngappedAlignment.alnLen);
                                bool elementHasSeqId = elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                                
                                if (!(elementHasAlnLen && elementHasCoverage && elementHasSeqId && elementHasEvalue)) {
                                    // 3. gapped alignment
                                    Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(concatedDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                                                        swMode, par.seqIdMode, false);
                                    if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                        allpass = false;
                                        break;
                                    }
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
                            ungappedAlignment.qStart, ungappedAlignment.qEnd, query.L,
                            ungappedAlignment.tStart, ungappedAlignment.tEnd, targetLength, backtrace);
                        appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), ungappedResult, par.addBacktrace);
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

                if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) continue;

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

                    // s_align gappedAlignment = blockAligner.align(&target, newQueryStartPos, newTargetStartPos, 
                                                                    //    gappedBacktrace, xDrop, par.covThr, par.covMode);
                    s_align gappedAlignment = blockAligner.bandedalign(&target, newQueryStartPos, newTargetStartPos, 
                                                                       gappedBacktrace, xDrop, par.covThr, par.covMode);
                    unsigned int gappedAlnLength = gappedBacktrace.size();
                    double gappedSeqId = Util::computeSeqId(par.seqIdMode, gappedAlignment.identicalAACnt, 
                                                           query.L, targetLength, gappedAlnLength);
                    Matcher::result_t result = Matcher::result_t(
                        targetKey, gappedAlignment.score1, gappedAlignment.qCov, gappedAlignment.tCov, 
                        gappedSeqId, gappedAlignment.evalue, gappedAlnLength,
                        gappedAlignment.qStartPos1, gappedAlignment.qEndPos1, query.L, 
                        gappedAlignment.dbStartPos1, gappedAlignment.dbEndPos1, targetLength, gappedBacktrace
                    );
                    if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, 
                                                par.alnLenThr, par.covMode, par.covThr)) {
                        if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) continue;
                        if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                            // check all the member from filtering file
                            const size_t cluId = cluDbr->getId(targetKey);
                            char *cluData = cluDbr->getData(cluId, threadIdx);
                            const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                            size_t numClu = Util::countLines(cluData, cluDataSize);
                            bool allpass = true;
                            char buffer[1024];
                            if (numClu > 1) { // if not singleton
                                while (*cluData != '\0') {
                                    Util::parseKey(cluData, buffer);
                                    const DBKeyType elementKey = Util::fast_atoi<DBKeyType>(buffer);
                                    if (elementKey == targetKey) {
                                        cluData = Util::skipLine(cluData);
                                        continue;
                                    }
                                    const size_t elementId = cluSeqDbr->getId(elementKey);
                                    char *elementSequence = cluSeqDbr->getData(elementId, threadIdx);
                                    size_t elementLength = cluSeqDbr->getSeqLen(elementId);
                                    // short concatedDiagonal = diagonal;
                                    short concatedDiagonal = 0;
                                    
                                    // 1. ungapped alignment
                                    element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                    if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
                                        allpass = false;
                                        break;
                                    }
                                    BlockAligner::UngappedAln_res elementUngappedAlignment = blockAligner.ungappedAlign(&element, concatedDiagonal);
                                    
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
                                    float elementSeqId = Util::computeSeqId(par.seqIdMode, elementIdenticalCount, query.L, elementLength, elementUngappedAlignment.alnLen);
                                    bool elementHasSeqId = elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                                    
                                    if (!(elementHasAlnLen && elementHasCoverage && elementHasSeqId && elementHasEvalue)) {
                                        // 3. gapped alignment
                                        Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(concatedDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                                                            swMode, par.seqIdMode, false);
                                        if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                            allpass = false;
                                            break;
                                        }
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
                        }
                        clusterResult.memberIds.push_back(targetId);
                    }
                }
            }

            if (includeAlignFiles) {
                alnWriter->writeData(alnResultBuffer.c_str(), alnResultBuffer.length(), queryKey, threadIdx);
            }
            pushClusterResult(std::move(clusterResult));
        }
    }

    {
        std::lock_guard<std::mutex> lock(clusterMutex);
        allCalculationsDone = true;
    }
    clusterCondition.notify_one();
    reorderSpaceCondition.notify_all();
    
    if (clusterThread.joinable()) {
        clusterThread.join(); 
    }

    for (size_t i = 0; i < dbSize; ++i) {
        if (loadAssignedCluster(assignedCluster, i) == DB_LOCAL_ID_INVALID) {
            storeAssignedCluster(assignedCluster, i, i);
        }
    }

    std::pair<DBKeyType, DBKeyType> *assignment = new std::pair<DBKeyType, DBKeyType>[dbSize];
    
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            const DBLocalId representativeId = loadAssignedCluster(assignedCluster, i);
            if (representativeId == DB_LOCAL_ID_INVALID) {
                Debug(Debug::ERROR) << "There must be an error: " << i 
                                    << " is not assigned to a cluster\n";
                continue;
            }

            assignment[i].first = seqDbr->getDbKey(representativeId);
            assignment[i].second = seqDbr->getDbKey(i);
        }
    }
    
    SORT_PARALLEL(assignment, assignment + dbSize);

    size_t clusterCount = (dbSize > 0) ? 1 : 0;
    for (size_t i = 1; i < dbSize; i++) {
        clusterCount += (assignment[i].first != assignment[i - 1].first);
    }

    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << clusterCount << "\n";
    
    writeData(&resultWriter, assignment, dbSize);
    
    delete[] assignedCluster;
    delete[] assignment;
    if (prefRepSizePair != nullptr) {
        delete[] prefRepSizePair;
    }
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
    alnDbr.open(DBReader<DBKeyType>::LINEAR_ACCCESS);
    int dbtype =  Parameters::DBTYPE_CLUSTER_RES;

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, dbtype);
    resultWriter.open();

    // Optional alignment-result output; path derived from the cluster DB (db3 + "_aln").
    DBWriter *alnWriter = nullptr;
    if (par.includeAlignFiles) {
        std::string alnDb = par.db3 + "_aln";
        std::string alnDbIndex = alnDb + ".index";
        alnWriter = new DBWriter(alnDb.c_str(), alnDbIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
        alnWriter->open();
    }

    int status = doAlign2clust(par, resultWriter, alnDbr, alnWriter);

    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";

    resultWriter.close();
    if (alnWriter != nullptr) {
        alnWriter->close();
        delete alnWriter;
    }
    alnDbr.close();

    return status;
}
