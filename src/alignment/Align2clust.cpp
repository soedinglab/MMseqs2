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
#include "MMseqsNode.h"
#include "FileUtil.h"
#include <atomic>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <ctime>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <sys/stat.h>

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
    DBLocalId len;
    size_t size;

    static bool compareBySizeAndId(const PrefInfo &first, const PrefInfo &second){
        if(first.size > second.size)
            return true;
        if(second.size > first.size)
            return false;
        if(first.len > second.len)
            return true;
        if(second.len > first.len)
            return false;
        if(first.id < second.id)
            return true;
        if(second.id < first.id)
            return false;
        return false;
    }
};

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
static std::vector<ClusterResult> reorderSlots;    // slot storage, size == reorderCapacity
static std::vector<unsigned char> reorderFilled;   // 1 == slot holds an unconsumed result
static size_t reorderCapacity = 0;                 // max out-of-order window
static size_t reorderBufferedCount = 0;            // number of filled slots
// Max-heap of candidates by member count (plain vector so the pool stays compactable).
static std::vector<SetCoverCandidate> setCoverCandidates;
static SetCoverComparator setCoverComparator;
// Shared backing store for all member ids referenced by setCoverCandidates.
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

// Env override for the reorder-buffer size; 0 (unset/invalid) means size from memory.
static size_t getReorderBufferLimitFromEnv() {
    const char *envValue = getenv("MMSEQS_ALIGN2CLUST_REORDER_LIMIT");
    if (envValue == nullptr || *envValue == '\0') {
        return 0;
    }

    char *end = nullptr;
    unsigned long long parsedValue = strtoull(envValue, &end, 10);
    if (end == envValue || *end != '\0' || parsedValue == 0) {
        Debug(Debug::WARNING) << "Ignoring invalid MMSEQS_ALIGN2CLUST_REORDER_LIMIT=" << envValue
                              << "; sizing the reorder buffer automatically\n";
        return 0;
    }
    return static_cast<size_t>(parsedValue);
}

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
    const size_t userCap = getReorderBufferLimitFromEnv();  // 0 == auto
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
        // Wait for this result's slot to free up. The next-in-order slot is always
        // free, so the producer the consumer is waiting on never blocks (deadlock-free).
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

static void writeClustering(DBWriter *dbWriter, const std::pair<DBKeyType, DBKeyType> * results, size_t dbSize) {
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

// the minimum ungapped score per column of a pair that is aligned with gaps
static float getScorePerColThreshold(const Parameters &par) {
    std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
    return parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.98);
}

struct AlignmentScoring {
    BaseMatrix *subMat;
    SubstitutionMatrix::FastMatrix fastMatrix;
    float scorePerColThreshold;
    EvalueComputation evaluer;

    AlignmentScoring(const Parameters &par, size_t dbResidues)
        : subMat(new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0)), fastMatrix(SubstitutionMatrix::createAsciiSubMat(*subMat)),
          scorePerColThreshold(getScorePerColThreshold(par)), evaluer(dbResidues, subMat) {}

    ~AlignmentScoring() {
        delete[] fastMatrix.matrix;
        delete[] fastMatrix.matrixData;
        delete subMat;
    }
};

// false if no diagonal can reach the required overlap
static bool canBeCoveredOnDiagonals(const Parameters &par, int queryLength, int targetLength, int &minDiag, int &maxDiag) {
    int diagBand = std::max(queryLength, targetLength) + targetLength;
    if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
        diagBand = static_cast<int>((1.0f - par.covThr) * std::max(queryLength, targetLength) + (1.0f - par.seqIdThr) * (queryLength + targetLength)) + 1;
    }
    const int minOverlap = static_cast<int>(std::ceil(Util::computeMinOverlap(Parameters::CLUST_LINEAR_DIAG_BAND_COV_SCALE * par.covThr, par.covMode,
                                                                            static_cast<float>(queryLength), static_cast<float>(targetLength))));
    minDiag = std::max(-diagBand, minOverlap - targetLength);
    maxDiag = std::min(diagBand, queryLength - minOverlap);
    return minDiag <= maxDiag;
}

// whether every other --filter-cludb-file member of the cluster of targetKey meets the clustering criteria against the query, from startDiagonal
static bool alignClusterMembers(const Parameters &par, BlockAligner &blockAligner, Matcher &matcher, bool &matcherInitialized, Sequence &query, Sequence &element,
                                DBReader<DBKeyType> *cluDbr, DBReader<DBKeyType> *cluSeqDbr, DBKeyType targetKey, short startDiagonal, unsigned int swMode,
                                unsigned int threadIdx, std::string *pendingMemberAln, char *alnLineBuffer) {
    const char *querySequence = query.getSeqData();
    const size_t cluId = cluDbr->getId(targetKey);
    char *cluData = cluDbr->getData(cluId, threadIdx);
    const size_t cluDataSize = cluDbr->getEntryLen(cluId);
    size_t numClu = Util::countLines(cluData, cluDataSize);
    char buffer[1024];
    if (pendingMemberAln != NULL) {
        pendingMemberAln->clear();
    }
    if (numClu <= 1) { // singleton
        return true;
    }
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
        short elementDiagonal = startDiagonal;

        // 1. ungapped alignment
        element.mapSequence(elementId, elementKey, elementSequence, elementLength);
        if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
            return false;
        }
        if (par.diagonalMode == Parameters::DIAGONAL_MODE_UNGAPPED) {
            int minDiag, maxDiag;
            if (canBeCoveredOnDiagonals(par, query.L, element.L, minDiag, maxDiag) == false) {
                return false;
            }
            elementDiagonal = blockAligner.bestDiagonal(&element, minDiag, maxDiag);
        }
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
        float elementSeqId = Util::computeSeqId(par.seqIdMode, elementIdenticalCount, query.L, elementLength, elementUngappedAlignment.alnLen);
        bool elementHasSeqId = elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());

        if (!(elementHasAlnLen && elementHasCoverage && elementHasSeqId && elementHasEvalue)) {
            // 3. gapped alignment
            if (matcherInitialized == false) {
                matcher.initQuery(&query);
                matcherInitialized = true;
            }
            Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(elementDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                                                swMode, par.seqIdMode, false);
            if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                return false;
            }
            // stage member alignment (flushed only if all members pass)
            if (pendingMemberAln != NULL) {
                appendAlignmentResult(*pendingMemberAln, alnLineBuffer, res_element, par.addBacktrace);
            }
        } else if (pendingMemberAln != NULL) {
            // member passed ungapped: gap-free (all 'M') record
            std::string elementBacktrace = par.addBacktrace ? std::string(elementUngappedAlignment.alnLen, 'M') : std::string();
            Matcher::result_t elementResult(elementKey, elementUngappedAlignment.score, elementUngappedAlignment.qcov,
                elementUngappedAlignment.tcov, elementSeqId, elementUngappedAlignment.eval, elementUngappedAlignment.alnLen,
                elementUngappedAlignment.qStart, elementUngappedAlignment.qEnd, query.L,
                elementUngappedAlignment.tStart, elementUngappedAlignment.tEnd, elementLength, elementBacktrace);
            appendAlignmentResult(*pendingMemberAln, alnLineBuffer, elementResult, par.addBacktrace);
        }
        cluData = Util::skipLine(cluData);
    }
    return true;
}

// whether the target meets the clustering criteria against the query of blockAligner, with its --filter-cludb-file members if cluDbr is given
static bool alignTarget(const Parameters &par, BlockAligner &blockAligner, bool &blockAlignerInitialized, Matcher *matcher, bool &matcherInitialized,
                        Sequence &query, Sequence &target, Sequence *element, DBReader<DBKeyType> *cluDbr, DBReader<DBKeyType> *cluSeqDbr,
                        float scorePerColThreshold, unsigned int swMode, unsigned short diagonal, unsigned int threadIdx,
                        const ClusterAssignment *assignedCluster, size_t targetId,
                        std::string *alnResultBuffer, std::string *pendingMemberAln, char *alnLineBuffer) {
    const char *querySequence = query.getSeqData();
    const char *targetSequence = target.getSeqData();
    if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
        return false;
    }

    if (par.diagonalMode == Parameters::DIAGONAL_MODE_UNGAPPED) {
        int minDiag, maxDiag;
        if (canBeCoveredOnDiagonals(par, query.L, target.L, minDiag, maxDiag) == false) {
            return false;
        }
        diagonal = blockAligner.bestDiagonal(&target, minDiag, maxDiag);
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
    const auto isAssigned = [&]() {
        return assignedCluster != NULL && loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID;
    };
    if (isAssigned()) {
        return false;
    }

    if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
        if (isAssigned()) {
            return false;
        }
        if (cluDbr != NULL && alignClusterMembers(par, blockAligner, *matcher, matcherInitialized, query, *element, cluDbr, cluSeqDbr,
                                                  target.getDbKey(), diagonal, swMode, threadIdx, pendingMemberAln, alnLineBuffer) == false) {
            return false;
        }
        if (alnResultBuffer != NULL) {
            std::string backtrace = par.addBacktrace ? std::string(ungappedAlignment.alnLen, 'M') : std::string();
            Matcher::result_t ungappedResult(target.getDbKey(), ungappedAlignment.score, ungappedAlignment.qcov,
                ungappedAlignment.tcov, seqId, ungappedAlignment.eval, ungappedAlignment.alnLen,
                ungappedAlignment.qStart, ungappedAlignment.qEnd, query.L,
                ungappedAlignment.tStart, ungappedAlignment.tEnd, target.L, backtrace);
            appendAlignmentResult(*alnResultBuffer, alnLineBuffer, ungappedResult, par.addBacktrace);
            // flush staged member alignments (empty unless filter gate ran)
            *alnResultBuffer += *pendingMemberAln;
        }
        return true;
    }

    float currentScorePerCol = static_cast<float>(ungappedAlignment.score) / static_cast<float>(ungappedAlignment.diagonalLen);
    if (currentScorePerCol < scorePerColThreshold) {
        return false;
    }

    int alignmentLength = ungappedAlignment.alnLen;
    int queryStartPos = ungappedAlignment.qStart;
    int targetStartPos = ungappedAlignment.tStart;
    int newQueryStartPos = queryStartPos;
    int newTargetStartPos = targetStartPos;

    if (queryStartPos == -1 || targetStartPos == -1 || alignmentLength < 3) {
        return false;
    }

    if (isAssigned()) {
        return false;
    }

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
    if (foundConsecutiveMatchSeed == false) {
        return false;
    }

    std::string gappedBacktrace;
    if (blockAlignerInitialized == false) {
        blockAligner.initQuery(&query);
        blockAlignerInitialized = true;
    }
    const int32_t xDrop = (MIN_SIZE * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid());
    s_align gappedAlignment = blockAligner.bandedalign(&target, newQueryStartPos, newTargetStartPos,
                                                       gappedBacktrace, xDrop, par.covThr, par.covMode);
    // bandedalign fails with evalue < 0 and an empty backtrace, whose length the default seq-id mode divides by
    if (gappedAlignment.evalue < 0 || gappedBacktrace.empty()) {
        return false;
    }
    unsigned int gappedAlnLength = gappedBacktrace.size();
    double gappedSeqId = Util::computeSeqId(par.seqIdMode, gappedAlignment.identicalAACnt,
                                           query.L, target.L, gappedAlnLength);
    Matcher::result_t result = Matcher::result_t(
        target.getDbKey(), gappedAlignment.score1, gappedAlignment.qCov, gappedAlignment.tCov,
        gappedSeqId, gappedAlignment.evalue, gappedAlnLength,
        gappedAlignment.qStartPos1, gappedAlignment.qEndPos1, query.L,
        gappedAlignment.dbStartPos1, gappedAlignment.dbEndPos1, target.L, gappedBacktrace
    );
    if (Alignment::checkCriteria(result, false, par.evalThr, par.seqIdThr,
                                 par.alnLenThr, par.covMode, par.covThr) == false) {
        return false;
    }
    if (isAssigned()) {
        return false;
    }
    if (cluDbr != NULL && alignClusterMembers(par, blockAligner, *matcher, matcherInitialized, query, *element, cluDbr, cluSeqDbr,
                                              target.getDbKey(), 0, swMode, threadIdx, pendingMemberAln, alnLineBuffer) == false) {
        return false;
    }
    if (alnResultBuffer != NULL) {
        appendAlignmentResult(*alnResultBuffer, alnLineBuffer, result, par.addBacktrace);
        // flush staged member alignments (empty unless filter gate ran)
        *alnResultBuffer += *pendingMemberAln;
    }
    return true;
}

static void (*clusterThreadFunc)(ClusterAssignment*) = nullptr;

void clusterThreadFuncSetcover(ClusterAssignment* assignedCluster) {
    while (true) {
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
                setCoverCandidates.push_back(candidate);
                std::push_heap(setCoverCandidates.begin(), setCoverCandidates.end(), setCoverComparator);
                continue;
            }

            for (size_t i = 0; i < candidate.memberCount; i++) {
                storeAssignedCluster(assignedCluster, members[i], candidate.representativeId);
            }
        }

        lock.unlock();
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
    seqDbr->open(DBReader<DBKeyType>::NOSORT);
 
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

    AlignmentScoring scoring(par, seqDbr->getAminoAcidDBSize());
    Debug(Debug::INFO) << "Score per column threshold for filtering: " << scoring.scorePerColThreshold << "\n";

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
        seqDbr->close();
        delete seqDbr;
        return EXIT_FAILURE;
    }

    const size_t align2clustResultCount = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();
    const size_t reorderCapacityChosen = computeReorderCapacity(par, dbSize, mode, align2clustResultCount);
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
    
    Timer timer;
    timer.reset();
    PrefInfo *prefRepSizePair = nullptr;
    DBLocalId *lengthOrder = nullptr;
    if (mode != Parameters::SET_COVER) {
        bool alreadyDescending = true;
        for (size_t i = 1; i < dbSize; i++) {
            if (seqDbr->getSeqLen(i) > seqDbr->getSeqLen(i - 1)) {
                alreadyDescending = false;
                break;
            }
        }
        if (alreadyDescending == false) {
            lengthOrder = new(std::nothrow) DBLocalId[dbSize];
            Util::checkAllocation(lengthOrder, "Can not allocate lengthOrder memory in align2clust");
#pragma omp parallel for schedule(static)
            for (size_t i = 0; i < dbSize; i++) {
                lengthOrder[i] = static_cast<DBLocalId>(i);
            }
            DBReader<DBKeyType> *order = seqDbr;
            SORT_PARALLEL(lengthOrder, lengthOrder + dbSize,
                          [order](DBLocalId first, DBLocalId second) {
                              const size_t firstLen = order->getSeqLen(first);
                              const size_t secondLen = order->getSeqLen(second);
                              if (firstLen != secondLen) {
                                  return firstLen > secondLen;
                              }
                              return first < second;
                          });
        }
        Debug(Debug::INFO) << "Length order: " << (lengthOrder == nullptr ? "db is already descending" : "built")
                           << " (" << timer.lap() << ")\n";
    }
    
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
                prefRepSizePair[i].len = static_cast<DBLocalId>(seqDbr->getSeqLen(prefRepSizePair[i].id));
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
        Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, db_maxseqlen, scoring.subMat, &scoring.evaluer,
                       par.compBiasCorrection, par.compBiasCorrectionScale,
                       par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
                       0.0, par.zdrop);
        Sequence query(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, scoring.subMat, 0, false, par.compBiasCorrection);
        Sequence target(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, scoring.subMat, 0, false, par.compBiasCorrection);
        Sequence element(db_maxseqlen, Parameters::DBTYPE_AMINO_ACIDS, scoring.subMat, 0, false, par.compBiasCorrection);
        BlockAligner blockAligner(Parameters::DBTYPE_AMINO_ACIDS, db_maxseqlen, scoring.subMat, &scoring.fastMatrix,
                                 &scoring.evaluer, par.compBiasCorrection, par.compBiasCorrectionScale,
                                 -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
        std::vector<std::pair<size_t, unsigned short>> targetsWithDiagonal;
        targetsWithDiagonal.reserve(1000);
        std::string querySeqBuffer;
        querySeqBuffer.reserve(db_maxseqlen + 1);

        const bool includeAlignFiles = (alnWriter != nullptr);
        std::string alnResultBuffer;
        // Staged member alignments; flushed only if the allpass-gate fully passes.
        std::string pendingMemberAln;
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
                representativeId = (lengthOrder != nullptr) ? lengthOrder[i] : static_cast<DBLocalId>(i);
                queryKey = seqDbr->getDbKey(representativeId);
                clusterResult.prefSize = 0;                         // greedy has no currentPrefSize gate
            }
            clusterResult.representativeId = representativeId;

            if (loadAssignedCluster(assignedCluster, representativeId) != DB_LOCAL_ID_INVALID) {
                pushClusterResult(std::move(clusterResult));
                continue;
            }

            const size_t alignmentId = alnDbr.getId(queryKey);
            char *alignmentData = alnDbr.getData(alignmentId, threadIdx);
            size_t prefSize = 0;
            while (*alignmentData != '\0') {
                hit_t hit = QueryMatcher::parsePrefilterHit(alignmentData);
                const size_t targetId = seqDbr->getId(hit.seqId);
                if (mode == Parameters::SET_COVER || loadAssignedCluster(assignedCluster, targetId) == DB_LOCAL_ID_INVALID) {
                    targetsWithDiagonal.push_back(std::make_pair(targetId, hit.diagonal));
                }
                alignmentData = Util::skipLine(alignmentData);
                prefSize++;
            }
            clusterResult.prefSize = prefSize;   // exact parsed count for the aligned path

            size_t queryId = representativeId;
            size_t queryLength = seqDbr->getSeqLen(queryId);
            // every entry holds the representative itself, so it is read only when another target is left
            if (targetsWithDiagonal.size() > 1) {
                querySeqBuffer.assign(seqDbr->getData(queryId, threadIdx), queryLength);
                query.mapSequence(queryId, queryKey, querySeqBuffer.c_str(), queryLength);
                blockAligner.setQuery(&query);
            }
            // the gapped aligners are only set up once an ungapped check fails
            bool blockAlignerInitialized = false;
            bool matcherInitialized = false;

            for (size_t targetIdx = 0; targetIdx < targetsWithDiagonal.size(); targetIdx++) {
                if (loadAssignedCluster(assignedCluster, representativeId) != DB_LOCAL_ID_INVALID) {
                    break;
                }

                const size_t targetId = targetsWithDiagonal[targetIdx].first;
                const DBKeyType targetKey = seqDbr->getDbKey(targetId);

                const bool isIdentity = (queryKey == targetKey);
                if (isIdentity) {
                    clusterResult.memberIds.push_back(queryId);
                    if (includeAlignFiles) {
                        std::string backtrace = par.addBacktrace ? std::string(queryLength, 'M') : std::string();
                        Matcher::result_t selfResult(queryKey, queryLength, 1.0f, 1.0f, 1.0f, 0.0,
                            queryLength, 0, queryLength - 1, queryLength, 0, queryLength - 1, queryLength, backtrace);
                        appendAlignmentResult(alnResultBuffer, alnLineBuffer.data(), selfResult, par.addBacktrace);
                    }
                    continue;
                }

                if (loadAssignedCluster(assignedCluster, targetId) != DB_LOCAL_ID_INVALID) {
                    continue;
                }

                char *targetSequence = seqDbr->getData(targetId, threadIdx);
                size_t targetLength = seqDbr->getSeqLen(targetId);
                target.mapSequence(targetId, targetKey, targetSequence, targetLength);
                if (alignTarget(par, blockAligner, blockAlignerInitialized, &matcher, matcherInitialized, query, target, &element, cluDbr, cluSeqDbr,
                                scoring.scorePerColThreshold, swMode, targetsWithDiagonal[targetIdx].second, threadIdx, assignedCluster, targetId,
                                includeAlignFiles ? &alnResultBuffer : NULL, includeAlignFiles ? &pendingMemberAln : NULL, alnLineBuffer.data())) {
                    clusterResult.memberIds.push_back(targetId);
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
    
    writeClustering(&resultWriter, assignment, dbSize);
    
    delete[] assignedCluster;
    delete[] assignment;
    if (prefRepSizePair != nullptr) {
        delete[] prefRepSizePair;
    }
    if (lengthOrder != nullptr) {
        delete[] lengthOrder;
    }
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

// --createdb-mode 3: a rep with its targets in a split: { rep, first target }, rep length, target count, then the rep residues and the target keys
struct RepTargets {
    ClusterPair pair;
    unsigned int repLength;
    unsigned int targetCount;

    // the size of the record at record, 0 if bytes do not hold all of it
    static size_t recordSize(const char *record, size_t bytes) {
        if (bytes < sizeof(RepTargets)) {
            return 0;
        }
        RepTargets header;
        memcpy(&header, record, sizeof(RepTargets));
        const size_t size = sizeof(RepTargets) + header.repLength + header.targetCount * sizeof(DBKeyType);
        return (bytes >= size) ? size : 0;
    }
};

// --createdb-mode 3: a row of <pref>.partitions, the entries of the keys from firstKey on in <pref>.<node>
struct PrefPartition {
    size_t firstKey;
    size_t prefOffset;
    size_t prefBytes;
    size_t cost;
    size_t node;
};

// --createdb-mode 3: a record that a node sends, where its targets start in the targets of the split and the node that aligns it
struct SendRecord {
    RepTargets header;
    size_t targetOffset;
    size_t node;

    static bool compareByNodeAndRep(const SendRecord &first, const SendRecord &second) {
        return (first.node != second.node) ? first.node < second.node : ClusterPair::compareByRep(first.header.pair, second.header.pair);
    }
};

// --createdb-mode 3: the targets are cut into blocks of blockKeys keys, which the nodes take in turn
static size_t getAlignNode(size_t blockKeys, size_t targetKey) {
    return (targetKey / blockKeys) % MMseqsNode::count;
}

// the pref partitions of all nodes are cut by cost into this many groups per node, which go to the nodes in turn
static const size_t PUSH_GROUPS_PER_NODE = 4;

// --createdb-mode 3: the nodes send the pref partitions in groups of equal cost in turn, a rep with its targets of a node in a split to that node
static std::string sendRepTargets(Parameters &par, DBReader<DBKeyType> &seqDbr, const std::vector<size_t> &nodeKeyOffsets, const std::vector<size_t> &splitKeyOffsets,
                                  size_t blockKeys, size_t memoryLimit, const std::string &repTargets) {
    const size_t dbSize = seqDbr.getSize();
    // the targets of a record: those from first on of its node (the targets of an entry are grouped by node) in the split of the first one, at most maxRecordTargets
    const size_t batchBytes = std::max<size_t>(1, memoryLimit / 16);
    const size_t maxRecordTargets = std::min<size_t>(UINT_MAX, std::max<size_t>(1, batchBytes / 2 / sizeof(DBKeyType)));
    const auto getRecordEnd = [&](const DBKeyType *targets, size_t first, size_t targetCount) {
        const size_t splitEnd = *std::upper_bound(splitKeyOffsets.begin(), splitKeyOffsets.end(), static_cast<size_t>(targets[first]));
        size_t nodeEnd = targetCount;
        if (MMseqsNode::count > 1) {
            const size_t node = getAlignNode(blockKeys, targets[first]);
            nodeEnd = std::partition_point(targets + first, targets + targetCount, [&](DBKeyType key) { return getAlignNode(blockKeys, key) == node; }) - targets;
        }
        return static_cast<size_t>(std::lower_bound(targets + first, targets + std::min(nodeEnd, first + maxRecordTargets), static_cast<DBKeyType>(splitEnd)) - targets);
    };

    const std::string partitionsFile = par.db2 + ".partitions";
    if (FileUtil::fileExists(partitionsFile.c_str()) == false) {
        Debug(Debug::ERROR) << "Cannot find " << partitionsFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    const std::string partitionsData = MMseqsNode::waitForFile(partitionsFile);
    std::vector<PrefPartition> partitions;
    size_t totalCost = 0;
    for (char *line = (char *) partitionsData.c_str(); *line != '\0'; line = Util::skipLine(line)) {
        const char *columns[4];
        if (Util::getWordsOfLine(line, columns, 4) < 4) {
            Debug(Debug::ERROR) << "Cannot read " << partitionsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        PrefPartition partition;
        partition.firstKey = Util::fast_atoi<size_t>(columns[0]);
        partition.prefOffset = Util::fast_atoi<size_t>(columns[1]);
        partition.prefBytes = Util::fast_atoi<size_t>(columns[2]);
        partition.cost = Util::fast_atoi<size_t>(columns[3]);
        partition.node = std::upper_bound(nodeKeyOffsets.begin(), nodeKeyOffsets.end(), partition.firstKey) - nodeKeyOffsets.begin() - 1;
        // a partition without targets has nothing to align
        if (partition.cost == 0) {
            continue;
        }
        if (partition.firstKey >= dbSize) {
            Debug(Debug::ERROR) << partitionsFile << " has a partition of key " << partition.firstKey << ", which is not in " << par.db1 << "\n";
            EXIT(EXIT_FAILURE);
        }
        partitions.push_back(partition);
        totalCost += partition.cost;
    }
    const size_t groupCount = PUSH_GROUPS_PER_NODE * MMseqsNode::count;
    const size_t groupCost = std::max<size_t>(1, (totalCost + groupCount - 1) / groupCount);
    std::vector<PrefPartition> nodePartitions;
    for (size_t i = 0, group = 0, costBefore = 0; i < partitions.size(); i++) {
        while (group + 1 < groupCount && costBefore >= (group + 1) * groupCost) {
            group++;
        }
        if (group % MMseqsNode::count == MMseqsNode::rank) {
            nodePartitions.push_back(partitions[i]);
        }
        costBefore += partitions[i].cost;
    }

    // a pref split with the next one and their targets, the reps of a batch, the next ones and the records of a batch take a 16th of the memory each
    const std::string repTargetsFile = repTargets + "." + SSTR(MMseqsNode::rank);
    FILE *repTargetsHandle = FileUtil::openAndDelete(repTargetsFile.c_str(), "w");
    std::string repTargetsIndexData;
    size_t repTargetsBytes = 0;
    std::vector<size_t> prefSplitOffsets;
    std::vector<DBKeyType> entryRepKeys;
    std::vector<size_t> entryTargetOffsets;
    std::vector<size_t> entryRecordOffsets;
    std::vector<size_t> entryRecordBytes;
    std::vector<DBKeyType> splitTargets;
    std::vector<size_t> batchEntryOffsets;
    std::vector<DBKeyType> repKeys;
    std::vector<DBKeyType> nextRepKeys;
    std::vector<SendRecord> records;
    std::vector<size_t> recordOffsets;
    std::vector<char> flushBuffer;
    std::vector<size_t> nodeOffsets(MMseqsNode::count + 1);
    DBReader<DBKeyType> *prefDb = NULL;
    for (size_t i = 0; i < nodePartitions.size(); i++) {
        const PrefPartition &partition = nodePartitions[i];
        const std::string prefFile = par.db2 + "." + SSTR(partition.node);
        if (i == 0 || partition.node != nodePartitions[i - 1].node) {
            if (prefDb != NULL) {
                prefDb->close();
                delete prefDb;
            }
            prefDb = new DBReader<DBKeyType>(prefFile.c_str(), (par.db2 + ".index." + SSTR(partition.node)).c_str(), par.threads,
                                             DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
            prefDb->open(DBReader<DBKeyType>::NOSORT);
        }
        // the entries of the partition in splits that end where an entry starts after prefSplitBytes more bytes
        const size_t prefEnd = partition.prefOffset + partition.prefBytes;
        const size_t partitionEndEntry = prefDb->findEntryByOffset(prefEnd);
        // a pref byte is read into two buffers and its targets, an entry also takes 64 bytes of offsets and keys
        const double entriesPerByte = static_cast<double>(partitionEndEntry - prefDb->findEntryByOffset(partition.prefOffset)) / std::max<size_t>(1, partition.prefBytes);
        const size_t prefSplitBytes = std::max<size_t>(1, static_cast<size_t>(static_cast<double>(memoryLimit / 16) / (3.0 + 64.0 * entriesPerByte)));
        prefSplitOffsets.clear();
        for (size_t offset = partition.prefOffset; offset < prefEnd; offset += prefSplitBytes) {
            const size_t entry = prefDb->findEntryByOffset(offset);
            if (prefSplitOffsets.empty() || prefSplitOffsets.back() != entry) {
                prefSplitOffsets.push_back(entry);
            }
        }
        prefSplitOffsets.push_back(partitionEndEntry);
        size_t nextRepKey = partition.firstKey;
        size_t partitionBytes = 0;
        for (size_t split = 0; split + 1 < prefSplitOffsets.size(); split++) {
            const size_t splitBegin = prefSplitOffsets[split];
            const size_t entryCount = prefSplitOffsets[split + 1] - splitBegin;
            prefDb->loadData(splitBegin, prefSplitOffsets[split + 1]);
            if (split + 2 < prefSplitOffsets.size()) {
                prefDb->prefetchData(prefSplitOffsets[split + 1], prefSplitOffsets[split + 2]);
            }
            entryRepKeys.resize(entryCount);
            entryTargetOffsets.assign(entryCount + 1, 0);
            entryRecordOffsets.assign(entryCount + 1, 0);
            entryRecordBytes.resize(entryCount);
#pragma omp parallel for schedule(dynamic, 100)
            for (size_t entry = 0; entry < entryCount; entry++) {
                const char *data = prefDb->getData(splitBegin + entry, 0);
                const size_t lines = std::count(data, data + prefDb->getEntryLen(splitBegin + entry), '\n');
                if (lines == 0) {
                    Debug(Debug::ERROR) << prefFile << " has an entry without its key\n";
                    EXIT(EXIT_FAILURE);
                }
                entryTargetOffsets[entry + 1] = lines - 1;
            }
            for (size_t entry = 0; entry < entryCount; entry++) {
                entryTargetOffsets[entry + 1] += entryTargetOffsets[entry];
                partitionBytes += prefDb->getEntryLen(splitBegin + entry);
            }
            splitTargets.resize(entryTargetOffsets[entryCount]);
#pragma omp parallel for schedule(dynamic, 100)
            for (size_t entry = 0; entry < entryCount; entry++) {
                char *data = prefDb->getData(splitBegin + entry, 0);
                const DBKeyType repKey = Util::fast_atoi<DBKeyType>(data);
                entryRepKeys[entry] = repKey;
                DBKeyType *targets = splitTargets.data() + entryTargetOffsets[entry];
                const size_t targetCount = entryTargetOffsets[entry + 1] - entryTargetOffsets[entry];
                data = Util::skipLine(data);
                for (size_t target = 0; target < targetCount; target++, data = Util::skipLine(data)) {
                    targets[target] = Util::fast_atoi<DBKeyType>(data);
                    if (targets[target] >= dbSize || targets[target] == repKey || (target > 0 && targets[target] <= targets[target - 1])) {
                        Debug(Debug::ERROR) << "Entry " << repKey << " of " << prefFile << " has the target " << targets[target] << " twice, out of order or not in " << par.db1 << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                // the targets of a node follow each other in key order, so a record holds the targets of one node and the rep goes to a node once a split
                if (MMseqsNode::count > 1) {
                    std::stable_sort(targets, targets + targetCount, [&](DBKeyType first, DBKeyType second) {
                        return getAlignNode(blockKeys, first) < getAlignNode(blockKeys, second);
                    });
                }
                const size_t repLength = (targetCount > 0) ? seqDbr.getSeqLen(repKey) : 0;
                entryRecordBytes[entry] = targetCount * sizeof(DBKeyType);
                for (size_t first = 0; first < targetCount; first = getRecordEnd(targets, first, targetCount)) {
                    entryRecordOffsets[entry + 1]++;
                    entryRecordBytes[entry] += sizeof(RepTargets) + repLength;
                }
            }
            // a batch takes entries while their reps and records fit, at least one
            batchEntryOffsets.assign(1, 0);
            for (size_t entry = 0, repBytes = 0, recordBytes = 0; entry < entryCount; entry++) {
                if (entryRepKeys[entry] < nextRepKey || entryRepKeys[entry] >= nodeKeyOffsets[partition.node + 1]) {
                    Debug(Debug::ERROR) << prefFile << " has the entry " << entryRepKeys[entry] << " out of key order or outside its partition\n";
                    EXIT(EXIT_FAILURE);
                }
                nextRepKey = entryRepKeys[entry] + 1;
                entryRecordOffsets[entry + 1] += entryRecordOffsets[entry];
                const size_t entryRepBytes = (entryTargetOffsets[entry + 1] > entryTargetOffsets[entry]) ? seqDbr.getEntryLen(entryRepKeys[entry]) : 0;
                if (entry > batchEntryOffsets.back() && (repBytes + entryRepBytes > batchBytes || recordBytes + entryRecordBytes[entry] > batchBytes)) {
                    batchEntryOffsets.push_back(entry);
                    repBytes = 0;
                    recordBytes = 0;
                }
                repBytes += entryRepBytes;
                recordBytes += entryRecordBytes[entry];
            }
            batchEntryOffsets.push_back(entryCount);
            const auto getRepKeys = [&](size_t batch, std::vector<DBKeyType> &keys) {
                keys.clear();
                for (size_t entry = batchEntryOffsets[batch]; entry < batchEntryOffsets[batch + 1]; entry++) {
                    if (entryTargetOffsets[entry + 1] > entryTargetOffsets[entry]) {
                        keys.push_back(entryRepKeys[entry]);
                    }
                }
            };
            getRepKeys(0, nextRepKeys);
            for (size_t batch = 0; batch + 1 < batchEntryOffsets.size(); batch++) {
                repKeys.swap(nextRepKeys);
                seqDbr.loadData(repKeys);
                if (batch + 2 < batchEntryOffsets.size()) {
                    getRepKeys(batch + 1, nextRepKeys);
                    seqDbr.prefetchData(nextRepKeys);
                }
                const size_t firstEntry = batchEntryOffsets[batch];
                const size_t endEntry = batchEntryOffsets[batch + 1];
                const size_t firstRecord = entryRecordOffsets[firstEntry];
                records.resize(entryRecordOffsets[endEntry] - firstRecord);
                if (records.empty()) {
                    continue;
                }
#pragma omp parallel for schedule(dynamic, 100)
                for (size_t entry = firstEntry; entry < endEntry; entry++) {
                    const DBKeyType *targets = splitTargets.data() + entryTargetOffsets[entry];
                    const size_t targetCount = entryTargetOffsets[entry + 1] - entryTargetOffsets[entry];
                    size_t record = entryRecordOffsets[entry] - firstRecord;
                    for (size_t first = 0, last; first < targetCount; first = last, record++) {
                        last = getRecordEnd(targets, first, targetCount);
                        RepTargets &header = records[record].header;
                        header.pair.repKey = entryRepKeys[entry];
                        header.pair.memberKey = targets[first];
                        header.repLength = static_cast<unsigned int>(seqDbr.getSeqLen(entryRepKeys[entry]));
                        header.targetCount = static_cast<unsigned int>(last - first);
                        records[record].targetOffset = entryTargetOffsets[entry] + first;
                        records[record].node = getAlignNode(blockKeys, targets[first]);
                    }
                }
                // the records of every node follow each other in rep order, in which the node aligns them
                SORT_PARALLEL(records.begin(), records.end(), SendRecord::compareByNodeAndRep);
                recordOffsets.resize(records.size() + 1);
                recordOffsets[0] = 0;
                for (size_t record = 0; record < records.size(); record++) {
                    recordOffsets[record + 1] = recordOffsets[record] + sizeof(RepTargets) + records[record].header.repLength + records[record].header.targetCount * sizeof(DBKeyType);
                }
                flushBuffer.resize(recordOffsets.back());
#pragma omp parallel
                {
                    unsigned int thread_idx = 0;
#ifdef OPENMP
                    thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 100)
                    for (size_t record = 0; record < records.size(); record++) {
                        const RepTargets &header = records[record].header;
                        char *data = flushBuffer.data() + recordOffsets[record];
                        memcpy(data, &header, sizeof(RepTargets));
                        memcpy(data + sizeof(RepTargets), seqDbr.getData(header.pair.repKey, thread_idx), header.repLength);
                        memcpy(data + sizeof(RepTargets) + header.repLength, splitTargets.data() + records[record].targetOffset, header.targetCount * sizeof(DBKeyType));
                    }
                }
                for (size_t node = 0; node < MMseqsNode::count; node++) {
                    const size_t nodeFirstRecord = std::partition_point(records.begin(), records.end(), [&](const SendRecord &record) {
                        return record.node < node;
                    }) - records.begin();
                    nodeOffsets[node] = recordOffsets[nodeFirstRecord];
                }
                nodeOffsets[MMseqsNode::count] = recordOffsets.back();
                if (fwrite(flushBuffer.data(), sizeof(char), flushBuffer.size(), repTargetsHandle) != flushBuffer.size()) {
                    Debug(Debug::ERROR) << "Cannot write " << repTargetsFile << "\n";
                    EXIT(EXIT_FAILURE);
                }
                MMseqsNode::addExchangeFlush(repTargetsIndexData, repTargetsBytes, nodeOffsets);
            }
        }
        if (partitionBytes != partition.prefBytes) {
            Debug(Debug::ERROR) << prefFile << " has an entry that ends after its partition\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (prefDb != NULL) {
        prefDb->close();
        delete prefDb;
    }
    return MMseqsNode::closeExchange(repTargetsHandle, repTargetsFile, repTargets + ".index", repTargetsIndexData);
}

// a node publishes its log after this many bytes or seconds since it published last
static const size_t PUBLISH_BYTES = 4 * 1024 * 1024;
static const time_t PUBLISH_SECONDS = 10;
// the records that a node aligns at once are cut into this many batches per thread, which the threads take in rep order
static const size_t ALIGN_BATCHES_PER_THREAD = 16;
// the master prefetches the assigned word of the pair this many pairs ahead as it takes pairs one by one
static const size_t PREFETCH_PAIRS = 16;
// the master clusters with at most this many threads, a quarter of its threads if fewer, besides the threads it aligns with
static const int CLUSTER_THREADS = 4;
// a block of targets: its bytes, the entries of every pref file sampled for the median target distance, and the blocks per node it must span
static const size_t BLOCK_BYTES = 16 * 1024 * 1024;
static const size_t BLOCK_SAMPLE_ENTRIES = 1024;
static const size_t BLOCKS_PER_DISTANCE = 8;

// --createdb-mode 3: a node aligns its records in rep order, a window at a time, logs accepted pairs, and skips reps and targets the master assigned
static void alignRepTargets(Parameters &par, DBReader<DBKeyType> &seqDbr, size_t blockKeys, size_t lagKeys, size_t memoryLimit, size_t representativeResidues,
                            const std::string &repTargetsIndex, const std::string &repTargets, const std::string &accepted, const std::string &assignedChanges) {
    const size_t dbSize = seqDbr.getSize();
    const size_t maxSeqLen = seqDbr.getMaxSeqLen();
    AlignmentScoring scoring(par, representativeResidues);
    std::vector<MMseqsNode::ExchangeBlock> repBlocks;
    std::vector<FILE *> repFiles;
    MMseqsNode::openExchange(repTargetsIndex, repTargets, repBlocks, repFiles);
    const std::string acceptedLog = accepted + "." + SSTR(MMseqsNode::rank);
    FILE *acceptedHandle = FileUtil::openAndDelete(acceptedLog.c_str(), "w");
    size_t acceptedBytes = 0;

    // the assigned bit ownKey(key) of ownWords for a target of this node, of a rep in ring word key / 64 modulo the ring size, as sent by the master
    const size_t blockCount = (dbSize + blockKeys - 1) / blockKeys;
    const size_t ownBlocks = (blockCount > MMseqsNode::rank) ? (blockCount - MMseqsNode::rank + MMseqsNode::count - 1) / MMseqsNode::count : 0;
    std::vector<uint64_t> ownWords(ownBlocks * (blockKeys / 64), 0);
    std::vector<uint64_t> ring(lagKeys / 64 + 1, 0);
    size_t ringFirstWord = 0;
    size_t walkedKey = 0;
    const auto getOwnKey = [&](size_t key) {
        return key / blockKeys / MMseqsNode::count * blockKeys + key % blockKeys;
    };
    const auto isAssigned = [&](size_t key) {
        if (getAlignNode(blockKeys, key) == MMseqsNode::rank) {
            const size_t ownKey = getOwnKey(key);
            return ((__atomic_load_n(&ownWords[ownKey / 64], __ATOMIC_RELAXED) >> (ownKey % 64)) & 1) != 0;
        }
        const size_t word = key / 64;
        const size_t firstWord = __atomic_load_n(&ringFirstWord, __ATOMIC_ACQUIRE);
        if (word < firstWord || word - firstWord >= ring.size()) {
            return false;
        }
        const uint64_t bits = __atomic_load_n(&ring[word % ring.size()], __ATOMIC_ACQUIRE);
        // the slot may hold a later word once the ring moved on meanwhile
        return __atomic_load_n(&ringFirstWord, __ATOMIC_ACQUIRE) == firstWord && ((bits >> (key % 64)) & 1) != 0;
    };
    const auto setRingKey = [&](size_t key) {
        const size_t word = key / 64;
        if (word >= ringFirstWord && word - ringFirstWord < ring.size()) {
            __atomic_fetch_or(&ring[word % ring.size()], 1ULL << (key % 64), __ATOMIC_RELEASE);
        }
    };
    // a step of the master appends its clustered key, first and count of entering ring words, ring key count, node key offsets, then those words and keys
    const size_t headerWords = 4 + MMseqsNode::count + 1;
    size_t changesBytes = 0;
    std::vector<uint64_t> changes;
    const auto readChanges = [&](size_t from, size_t words) {
        changes.resize(words);
        if (words > 0) {
            MMseqsNode::readLog(assignedChanges, from, from + words * sizeof(uint64_t), reinterpret_cast<char *>(changes.data()));
        }
    };
    const auto readAssignedChanges = [&]() {
        size_t publishedKey;
        size_t publishedBytes;
        if (MMseqsNode::readLogProgress(assignedChanges, publishedKey, publishedBytes) == false) {
            return;
        }
        while (changesBytes + headerWords * sizeof(uint64_t) <= publishedBytes) {
            readChanges(changesBytes, headerWords);
            const size_t stepWalkedKey = changes[0];
            const size_t sliceFirstWord = changes[1];
            const size_t sliceWords = changes[2];
            const size_t lagKeyCount = changes[3];
            const size_t nodeKeyOffset = changes[4 + MMseqsNode::rank];
            const size_t nodeKeyEnd = changes[5 + MMseqsNode::rank];
            const size_t commonStart = changesBytes + headerWords * sizeof(uint64_t);
            const size_t ownStart = commonStart + (sliceWords + lagKeyCount + nodeKeyOffset) * sizeof(uint64_t);
            changesBytes = commonStart + (sliceWords + lagKeyCount + changes[4 + MMseqsNode::count]) * sizeof(uint64_t);
            if (changesBytes > publishedBytes || stepWalkedKey < walkedKey || nodeKeyEnd < nodeKeyOffset) {
                Debug(Debug::ERROR) << "Cannot read " << assignedChanges << "\n";
                EXIT(EXIT_FAILURE);
            }
            // the words before the clustered key leave the ring, cleared before the release of ringFirstWord, and their slots take the words that enter it
            const size_t firstWord = stepWalkedKey / 64;
            for (size_t word = ringFirstWord; word < std::min(firstWord, ringFirstWord + ring.size()); word++) {
                __atomic_store_n(&ring[word % ring.size()], 0, __ATOMIC_RELAXED);
            }
            __atomic_store_n(&ringFirstWord, firstWord, __ATOMIC_RELEASE);
            walkedKey = stepWalkedKey;
            readChanges(commonStart, sliceWords + lagKeyCount);
            for (size_t i = 0; i < sliceWords; i++) {
                if (sliceFirstWord + i >= firstWord && sliceFirstWord + i - firstWord < ring.size()) {
                    __atomic_fetch_or(&ring[(sliceFirstWord + i) % ring.size()], changes[i], __ATOMIC_RELEASE);
                }
            }
            for (size_t i = sliceWords; i < changes.size(); i++) {
                setRingKey(changes[i]);
            }
            readChanges(ownStart, nodeKeyEnd - nodeKeyOffset);
            for (size_t i = 0; i < changes.size(); i++) {
                if (getAlignNode(blockKeys, changes[i]) != MMseqsNode::rank) {
                    Debug(Debug::ERROR) << "Cannot read " << assignedChanges << ": a key of another node\n";
                    EXIT(EXIT_FAILURE);
                }
                const size_t ownKey = getOwnKey(changes[i]);
                __atomic_fetch_or(&ownWords[ownKey / 64], 1ULL << (ownKey % 64), __ATOMIC_RELAXED);
            }
        }
    };

    // a chunk holds the records that every block buffers in rep order, read into the chunk that the window before the current one does not use
    std::vector<char> chunkRecords[2];
    std::vector<size_t> chunkRecordOffsets[2];
    std::vector<std::pair<ClusterPair, size_t> > chunkOrder[2];
    std::vector<size_t> chunkTargetBytes[2];
    size_t chunkPos[2] = { 0, 0 };
    size_t chunk = 0;
    const ClusterPair lastRecord = { DB_KEY_INVALID, DB_KEY_INVALID };
    // a window takes records while their targets fit into an eighth of the memory, at least one record, whose targets are in one split
    const auto readWindow = [&](size_t &windowChunk, size_t &windowBegin, size_t &windowEnd) {
        if (chunkPos[chunk] == chunkOrder[chunk].size()) {
            chunk = 1 - chunk;
            const size_t recordCount = MMseqsNode::readExchangeRange(repBlocks, memoryLimit / 16, RepTargets::recordSize, ClusterPair::compareByRep, lastRecord,
                                                                     chunkRecords[chunk], &chunkRecordOffsets[chunk], repTargets);
            chunkOrder[chunk].resize(recordCount);
            chunkTargetBytes[chunk].resize(recordCount);
#pragma omp parallel for schedule(dynamic, 100)
            for (size_t i = 0; i < recordCount; i++) {
                const char *record = chunkRecords[chunk].data() + chunkRecordOffsets[chunk][i];
                RepTargets header;
                memcpy(&header, record, sizeof(RepTargets));
                chunkOrder[chunk][i] = std::make_pair(header.pair, i);
                size_t targetBytes = 0;
                for (unsigned int j = 0; j < header.targetCount; j++) {
                    DBKeyType targetKey;
                    memcpy(&targetKey, record + sizeof(RepTargets) + header.repLength + j * sizeof(DBKeyType), sizeof(DBKeyType));
                    if (j == 0 && getAlignNode(blockKeys, targetKey) != MMseqsNode::rank) {
                        Debug(Debug::ERROR) << "Cannot read " << repTargets << ": a record has a target of another node\n";
                        EXIT(EXIT_FAILURE);
                    }
                    targetBytes += seqDbr.getEntryLen(targetKey);
                }
                chunkTargetBytes[chunk][i] = targetBytes;
            }
            SORT_PARALLEL(chunkOrder[chunk].begin(), chunkOrder[chunk].end(), [](const std::pair<ClusterPair, size_t> &first, const std::pair<ClusterPair, size_t> &second) {
                return ClusterPair::compareByRep(first.first, second.first);
            });
            chunkPos[chunk] = 0;
        }
        windowChunk = chunk;
        windowBegin = chunkPos[chunk];
        windowEnd = windowBegin;
        for (size_t windowBytes = 0; windowEnd < chunkOrder[chunk].size(); windowEnd++) {
            windowBytes += chunkTargetBytes[chunk][chunkOrder[chunk][windowEnd].second];
            if (windowEnd > windowBegin && windowBytes > memoryLimit / 8) {
                break;
            }
        }
        chunkPos[chunk] = windowEnd;
    };
    // the targets to read for a window, without those whose rep or target is assigned
    std::vector<std::vector<DBKeyType> > threadKeys(par.threads);
    const auto getWindowKeys = [&](size_t windowChunk, size_t windowBegin, size_t windowEnd, std::vector<DBKeyType> &keys) {
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            threadKeys[thread_idx].clear();
#pragma omp for schedule(dynamic, 100)
            for (size_t i = windowBegin; i < windowEnd; i++) {
                const char *record = chunkRecords[windowChunk].data() + chunkRecordOffsets[windowChunk][chunkOrder[windowChunk][i].second];
                RepTargets header;
                memcpy(&header, record, sizeof(RepTargets));
                if (isAssigned(header.pair.repKey)) {
                    continue;
                }
                for (unsigned int j = 0; j < header.targetCount; j++) {
                    DBKeyType targetKey;
                    memcpy(&targetKey, record + sizeof(RepTargets) + header.repLength + j * sizeof(DBKeyType), sizeof(DBKeyType));
                    if (isAssigned(targetKey) == false) {
                        threadKeys[thread_idx].push_back(targetKey);
                    }
                }
            }
        }
        keys.clear();
        for (int thread = 0; thread < par.threads; thread++) {
            keys.insert(keys.end(), threadKeys[thread].begin(), threadKeys[thread].end());
        }
        SORT_PARALLEL(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    };

    // the windows read bands of keys from their first key on that fit an eighth of the memory: whole on one node, the own blocks on several
    size_t bandBegin = 0;
    size_t bandEnd = 0;
    std::vector<DBKeyType> bandKeys;
    const auto readBand = [&](const std::vector<DBKeyType> &windowKeys, bool isPrefetch) {
        if (windowKeys.empty() || (windowKeys.front() >= bandBegin && windowKeys.back() < bandEnd)) {
            return;
        }
        const size_t begin = windowKeys.front();
        const size_t end = seqDbr.getSplitKeyOffsets(begin, dbSize, memoryLimit / 8 * MMseqsNode::count, 0)[1];
        const bool fitsBand = windowKeys.back() < end
                              && (MMseqsNode::count == 1 || (MMseqsNode::count - 1) * blockKeys * seqDbr.getEntryLen(end - 1) > DBReader<DBKeyType>::READ_MERGE_GAP);
        const bool readsRange = fitsBand && MMseqsNode::count == 1;
        const std::vector<DBKeyType> *keysToRead = &windowKeys;
        if (fitsBand && MMseqsNode::count > 1) {
            bandKeys.clear();
            for (size_t block = begin / blockKeys; block * blockKeys < end; block++) {
                if (getAlignNode(blockKeys, block * blockKeys) != MMseqsNode::rank) {
                    continue;
                }
                for (size_t key = std::max(begin, block * blockKeys); key < std::min(end, (block + 1) * blockKeys); key++) {
                    bandKeys.push_back(key);
                }
            }
            keysToRead = &bandKeys;
        }
        if (isPrefetch) {
            if (readsRange) {
                seqDbr.prefetchData(begin, end);
            } else {
                seqDbr.prefetchData(*keysToRead);
            }
            return;
        }
        bandBegin = begin;
        bandEnd = fitsBand ? end : begin;
        if (readsRange) {
            seqDbr.loadData(begin, end);
        } else {
            seqDbr.loadData(*keysToRead);
        }
    };

    size_t windowChunk;
    size_t windowBegin;
    size_t windowEnd;
    readWindow(windowChunk, windowBegin, windowEnd);
    const size_t firstRepKey = (windowEnd > windowBegin) ? static_cast<size_t>(chunkOrder[windowChunk][windowBegin].first.repKey) : dbSize;
    MMseqsNode::publishLog(acceptedHandle, acceptedLog, firstRepKey, 0);
    size_t publishedBytes = 0;
    time_t publishTime = time(NULL);
    std::mutex acceptedMutex;
    std::vector<DBKeyType> keys;
    std::vector<DBKeyType> nextKeys;
    std::vector<unsigned char> batchDone;
    std::vector<std::vector<ClusterPair> > batchPairs;
    getWindowKeys(windowChunk, windowBegin, windowEnd, keys);
    while (windowEnd > windowBegin) {
        // the reps of a window are checked in the ring, which the master fills to lagKeys after the clustered key
        const size_t windowRepKey = chunkOrder[windowChunk][windowBegin].first.repKey;
        if (windowRepKey >= walkedKey + lagKeys) {
            MMseqsNode::publishLog(acceptedHandle, acceptedLog, windowRepKey, acceptedBytes);
            publishedBytes = acceptedBytes;
            for (unsigned int polls = 0; windowRepKey >= walkedKey + lagKeys; polls++) {
                MMseqsNode::sleepBeforePoll(polls);
                readAssignedChanges();
            }
            publishTime = time(NULL);
        }
        readBand(keys, false);
        size_t nextChunk;
        size_t nextBegin;
        size_t nextEnd;
        readWindow(nextChunk, nextBegin, nextEnd);
        getWindowKeys(nextChunk, nextBegin, nextEnd, nextKeys);
        readBand(nextKeys, true);
        const size_t nextRepKey = (nextEnd > nextBegin) ? static_cast<size_t>(chunkOrder[nextChunk][nextBegin].first.repKey) : dbSize;
        const std::vector<std::pair<ClusterPair, size_t> > &order = chunkOrder[windowChunk];
        const size_t recordCount = windowEnd - windowBegin;
        const size_t batchCount = std::min(recordCount, par.threads * ALIGN_BATCHES_PER_THREAD);
        batchDone.assign(batchCount, 0);
        batchPairs.assign(batchCount, std::vector<ClusterPair>());
        size_t firstUnfinishedBatch = 0;
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence query(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, scoring.subMat, 0, false, par.compBiasCorrection);
            Sequence target(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, scoring.subMat, 0, false, par.compBiasCorrection);
            BlockAligner blockAligner(Parameters::DBTYPE_AMINO_ACIDS, maxSeqLen, scoring.subMat, &scoring.fastMatrix,
                                     &scoring.evaluer, par.compBiasCorrection, par.compBiasCorrectionScale,
                                     -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
            std::string querySeqBuffer;
            std::vector<ClusterPair> pairs;
#pragma omp for schedule(dynamic, 1)
            for (size_t batch = 0; batch < batchCount; batch++) {
                pairs.clear();
                for (size_t i = windowBegin + recordCount * batch / batchCount; i < windowBegin + recordCount * (batch + 1) / batchCount; i++) {
                    const char *record = chunkRecords[windowChunk].data() + chunkRecordOffsets[windowChunk][order[i].second];
                    RepTargets header;
                    memcpy(&header, record, sizeof(RepTargets));
                    if (isAssigned(header.pair.repKey)) {
                        continue;
                    }
                    querySeqBuffer.assign(record + sizeof(RepTargets), header.repLength);
                    query.mapSequence(header.pair.repKey, header.pair.repKey, querySeqBuffer.c_str(), header.repLength);
                    blockAligner.setQuery(&query);
                    bool blockAlignerInitialized = false;
                    bool matcherInitialized = false;
                    const char *targetKeys = record + sizeof(RepTargets) + header.repLength;
                    for (unsigned int j = 0; j < header.targetCount; j++) {
                        DBKeyType targetKey;
                        memcpy(&targetKey, targetKeys + j * sizeof(DBKeyType), sizeof(DBKeyType));
                        // a target of an assigned rep is not loaded, though the ring may have moved past the rep since
                        if (isAssigned(targetKey) || std::binary_search(keys.begin(), keys.end(), targetKey) == false) {
                            continue;
                        }
                        target.mapSequence(targetKey, targetKey, seqDbr.getData(targetKey, thread_idx), seqDbr.getSeqLen(targetKey));
                        if (alignTarget(par, blockAligner, blockAlignerInitialized, NULL, matcherInitialized, query, target, NULL, NULL, NULL,
                                        scoring.scorePerColThreshold, 0, 0, thread_idx, NULL, 0, NULL, NULL, NULL)) {
                            const ClusterPair pair = { header.pair.repKey, targetKey };
                            pairs.push_back(pair);
                        }
                    }
                }
                // the pairs of the batches go to the log in batch order, so the reps before the first unfinished batch are done and their pairs are in the log
                std::lock_guard<std::mutex> lock(acceptedMutex);
                batchPairs[batch].swap(pairs);
                batchDone[batch] = 1;
                for (; firstUnfinishedBatch < batchCount && batchDone[firstUnfinishedBatch]; firstUnfinishedBatch++) {
                    std::vector<ClusterPair> &donePairs = batchPairs[firstUnfinishedBatch];
                    if (fwrite(donePairs.data(), sizeof(ClusterPair), donePairs.size(), acceptedHandle) != donePairs.size()) {
                        Debug(Debug::ERROR) << "Cannot write " << acceptedLog << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    acceptedBytes += donePairs.size() * sizeof(ClusterPair);
                    std::vector<ClusterPair>().swap(donePairs);
                }
                const time_t now = time(NULL);
                if (acceptedBytes - publishedBytes >= PUBLISH_BYTES || now - publishTime >= PUBLISH_SECONDS) {
                    const size_t doneRepKey = (firstUnfinishedBatch < batchCount) ? static_cast<size_t>(order[windowBegin + recordCount * firstUnfinishedBatch / batchCount].first.repKey) : nextRepKey;
                    MMseqsNode::publishLog(acceptedHandle, acceptedLog, doneRepKey, acceptedBytes);
                    publishedBytes = acceptedBytes;
                    publishTime = now;
                    readAssignedChanges();
                }
            }
        }
        keys.swap(nextKeys);
        windowChunk = nextChunk;
        windowBegin = nextBegin;
        windowEnd = nextEnd;
    }
    // readExchangeRange stops at a block that is out of order, so every block must be read to its end
    for (size_t i = 0; i < repBlocks.size(); i++) {
        if (repBlocks[i].pos != repBlocks[i].end || repBlocks[i].bufferPos != repBlocks[i].bufferEnd) {
            Debug(Debug::ERROR) << "Cannot read " << repTargets << ": its records are not in rep order\n";
            EXIT(EXIT_FAILURE);
        }
    }
    MMseqsNode::publishLog(acceptedHandle, acceptedLog, dbSize, acceptedBytes);
    if (fclose(acceptedHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << acceptedLog << "\n";
        EXIT(EXIT_FAILURE);
    }
    MMseqsNode::closeExchangeFiles(repFiles, repTargets);
    seqDbr.close();
}

// an accepted pair of a node with its bit in the taken bits of the pairs that the master reads from its log in a step
struct AcceptedPair {
    ClusterPair pair;
    size_t takenBit;
    unsigned int node;

    static bool compareByRep(const AcceptedPair &first, const AcceptedPair &second) {
        return ClusterPair::compareByRep(first.pair, second.pair);
    }

    static bool compareByRepKey(const AcceptedPair &first, const AcceptedPair &second) {
        return first.pair.repKey < second.pair.repKey;
    }
};

// the taken bits of the pairs of a log in whole words; the bits of the last partial one wait in carry
struct TakenBits {
    FILE *file;
    size_t bytes;
    uint64_t carry;
    unsigned int carryBits;
};

static void appendTakenBits(TakenBits &taken, const uint64_t *bits, size_t bitCount, std::vector<uint64_t> &words, const std::string &fileName) {
    words.clear();
    for (size_t i = 0; i < bitCount; i += 64) {
        const unsigned int wordBits = std::min<size_t>(64, bitCount - i);
        taken.carry |= bits[i / 64] << taken.carryBits;
        if (taken.carryBits + wordBits >= 64) {
            words.push_back(taken.carry);
            taken.carry = (taken.carryBits > 0) ? bits[i / 64] >> (64 - taken.carryBits) : 0;
            taken.carryBits = taken.carryBits + wordBits - 64;
        } else {
            taken.carryBits += wordBits;
        }
    }
    if (fwrite(words.data(), sizeof(uint64_t), words.size(), taken.file) != words.size()) {
        Debug(Debug::ERROR) << "Cannot write " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    taken.bytes += words.size() * sizeof(uint64_t);
}

// --createdb-mode 3: the master follows the logs of accepted pairs and clusters them greedily in key order; returns the number of taken pairs
static size_t clusterAcceptedPairs(int threads, size_t dbSize, size_t memoryLimit, size_t blockKeys, size_t lagKeys,
                                   const std::string &accepted, const std::string &assignedChanges, const std::string &taken, const std::string &assignedKeys) {
    const size_t assignedWords = (dbSize + 63) / 64;
    uint64_t *assigned = new(std::nothrow) uint64_t[assignedWords];
    Util::checkAllocation(assigned, "Can not allocate assigned memory in align2clust");
    memset(assigned, 0, assignedWords * sizeof(uint64_t));
    std::vector<std::vector<DBKeyType> > threadAssignedKeys(threads);
    const auto isAssigned = [&](DBKeyType key) {
        return ((assigned[key / 64] >> (key % 64)) & 1) != 0;
    };
    const auto setAssigned = [&](DBKeyType key, std::vector<DBKeyType> &keys) {
        if ((__sync_fetch_and_or(&assigned[key / 64], 1ULL << (key % 64)) & (1ULL << (key % 64))) == 0) {
            keys.push_back(key);
        }
    };

    // the first rep that a node is not done with and the bytes of its log, once it has published them, and the bytes of it that are clustered
    std::vector<unsigned char> nodeHasProgress(MMseqsNode::count, 0);
    std::vector<size_t> nodeRepKeys(MMseqsNode::count, 0);
    std::vector<size_t> nodeLogBytes(MMseqsNode::count, 0);
    std::vector<size_t> nodeReadBytes(MMseqsNode::count, 0);
    std::vector<TakenBits> nodeTaken(MMseqsNode::count);
    std::vector<std::vector<uint64_t> > nodeTakenBits(MMseqsNode::count);
    std::vector<size_t> nodeTakenCounts(MMseqsNode::count, 0);
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        const TakenBits takenBits = { FileUtil::openAndDelete((taken + "." + SSTR(node)).c_str(), "w"), 0, 0, 0 };
        nodeTaken[node] = takenBits;
    }
    FILE *changesHandle = FileUtil::openAndDelete(assignedChanges.c_str(), "w");
    size_t changesBytes = 0;
    // a step clusters at most a 16th of the memory of pairs, each log read up to its share unless the rep it starts in spans more
    const size_t stepPairs = std::max<size_t>(1, memoryLimit / 16 / sizeof(AcceptedPair));
    const size_t nodePairs = std::max<size_t>(1, stepPairs / MMseqsNode::count);
    const size_t partBytes = std::max<size_t>(1, std::min(nodePairs * sizeof(ClusterPair), memoryLimit / 64) / sizeof(ClusterPair)) * sizeof(ClusterPair);
    std::vector<AcceptedPair> readPairs;
    std::vector<size_t> nodeFirstPairs(MMseqsNode::count + 1, 0);
    std::vector<char> logData;
    std::vector<size_t> repOffsets;
    std::vector<std::vector<size_t> > threadRepOffsets(threads);
    std::vector<unsigned char> pairIsFree;
    std::vector<size_t> repKeyOffsets;
    std::vector<std::pair<DBKeyType, size_t> > keyReps;
    std::vector<unsigned char> repSharesKey;
    std::vector<uint64_t> words;
    std::vector<uint64_t> record;
    std::vector<std::vector<size_t> > threadSectionCounts(threads, std::vector<size_t>(MMseqsNode::count + 1));
    std::vector<size_t> sectionOffsets(MMseqsNode::count + 2);
    size_t memberCount = 0;
    unsigned int polls = 0;
    for (size_t walkedKey = 0; walkedKey < dbSize;) {
        bool hasNewProgress = false;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            const std::string log = accepted + "." + SSTR(node);
            size_t repKey;
            size_t logBytes;
            if ((nodeHasProgress[node] && nodeRepKeys[node] == dbSize) || MMseqsNode::readLogProgress(log, repKey, logBytes) == false
                || (nodeHasProgress[node] && repKey == nodeRepKeys[node] && logBytes == nodeLogBytes[node])) {
                continue;
            }
            if (logBytes % sizeof(ClusterPair) != 0 || logBytes < nodeLogBytes[node] || repKey < nodeRepKeys[node]) {
                Debug(Debug::ERROR) << "Cannot read " << log << ".progress\n";
                EXIT(EXIT_FAILURE);
            }
            hasNewProgress = true;
            nodeHasProgress[node] = 1;
            nodeRepKeys[node] = repKey;
            nodeLogBytes[node] = logBytes;
        }
        size_t doneRepKey = dbSize;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            doneRepKey = std::min(doneRepKey, nodeHasProgress[node] ? nodeRepKeys[node] : 0);
        }
        if (doneRepKey <= walkedKey) {
            if (hasNewProgress == false) {
                MMseqsNode::sleepBeforePoll(polls++);
            }
            continue;
        }
        polls = 0;

        // the pairs of the reps before doneRepKey from every log; the step clusters the reps before the first one of a log that is not read in full
        readPairs.clear();
        size_t stepRepKey = doneRepKey;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            const std::string log = accepted + "." + SSTR(node);
            nodeFirstPairs[node] = readPairs.size();
            size_t nodeRepKey = doneRepKey;
            for (size_t from = nodeReadBytes[node]; from < nodeLogBytes[node] && nodeRepKey == doneRepKey;) {
                const size_t to = std::min(nodeLogBytes[node], from + partBytes);
                logData.resize(to - from);
                MMseqsNode::readLog(log, from, to, logData.data());
                for (size_t i = 0; from < to; i++, from += sizeof(ClusterPair)) {
                    AcceptedPair pair;
                    memcpy(&pair.pair, logData.data() + i * sizeof(ClusterPair), sizeof(ClusterPair));
                    pair.takenBit = readPairs.size() - nodeFirstPairs[node];
                    pair.node = static_cast<unsigned int>(node);
                    if (pair.pair.repKey >= doneRepKey) {
                        break;
                    }
                    if (pair.takenBit >= nodePairs && pair.pair.repKey != readPairs[nodeFirstPairs[node]].pair.repKey) {
                        nodeRepKey = pair.pair.repKey;
                        break;
                    }
                    readPairs.push_back(pair);
                }
                if (from < to) {
                    break;
                }
            }
            stepRepKey = std::min(stepRepKey, nodeRepKey);
        }
        nodeFirstPairs[MMseqsNode::count] = readPairs.size();
        // the pairs of every node before stepRepKey, which start its pairs, are clustered; the others are read again in a later step
        size_t pairCount = 0;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            const AcceptedPair firstNotDone = { { static_cast<DBKeyType>(stepRepKey), static_cast<DBKeyType>(stepRepKey) }, 0, 0 };
            const size_t nodeEnd = std::lower_bound(readPairs.begin() + nodeFirstPairs[node], readPairs.begin() + nodeFirstPairs[node + 1], firstNotDone, AcceptedPair::compareByRepKey) - readPairs.begin();
            nodeTakenCounts[node] = nodeEnd - nodeFirstPairs[node];
            nodeReadBytes[node] += nodeTakenCounts[node] * sizeof(ClusterPair);
            nodeTakenBits[node].assign((nodeTakenCounts[node] + 63) / 64, 0);
            if (pairCount != nodeFirstPairs[node]) {
                memmove(&readPairs[pairCount], &readPairs[nodeFirstPairs[node]], nodeTakenCounts[node] * sizeof(AcceptedPair));
            }
            pairCount += nodeTakenCounts[node];
        }
        SORT_PARALLEL(readPairs.begin(), readPairs.begin() + pairCount, AcceptedPair::compareByRep);
        const AcceptedPair *pairs = readPairs.data();
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            threadRepOffsets[thread_idx].clear();
#pragma omp for schedule(static)
            for (size_t i = 0; i < pairCount; i++) {
                if (i == 0 || pairs[i].pair.repKey != pairs[i - 1].pair.repKey) {
                    threadRepOffsets[thread_idx].push_back(i);
                }
            }
        }
        repOffsets.clear();
        for (int thread = 0; thread < threads; thread++) {
            repOffsets.insert(repOffsets.end(), threadRepOffsets[thread].begin(), threadRepOffsets[thread].end());
        }
        repOffsets.push_back(pairCount);
        const size_t repCount = repOffsets.size() - 1;

        // a pair whose rep or member is assigned before the step stays out, as bits are never cleared; a rep with free members may take them and itself
        pairIsFree.resize(pairCount);
        repKeyOffsets.assign(repCount + 1, 0);
#pragma omp parallel for schedule(dynamic, 1024)
        for (size_t rep = 0; rep < repCount; rep++) {
            const bool repIsFree = isAssigned(pairs[repOffsets[rep]].pair.repKey) == false;
            size_t freeCount = 0;
            for (size_t i = repOffsets[rep]; i < repOffsets[rep + 1]; i++) {
                pairIsFree[i] = repIsFree && isAssigned(pairs[i].pair.memberKey) == false;
                freeCount += pairIsFree[i];
            }
            repKeyOffsets[rep + 1] = (freeCount > 0) ? freeCount + 1 : 0;
        }
        for (size_t rep = 0; rep < repCount; rep++) {
            repKeyOffsets[rep + 1] += repKeyOffsets[rep];
        }
        keyReps.resize(repKeyOffsets[repCount]);
#pragma omp parallel for schedule(dynamic, 1024)
        for (size_t rep = 0; rep < repCount; rep++) {
            size_t key = repKeyOffsets[rep];
            if (key == repKeyOffsets[rep + 1]) {
                continue;
            }
            keyReps[key++] = std::make_pair(pairs[repOffsets[rep]].pair.repKey, rep);
            for (size_t i = repOffsets[rep]; i < repOffsets[rep + 1]; i++) {
                if (pairIsFree[i]) {
                    keyReps[key++] = std::make_pair(pairs[i].pair.memberKey, rep);
                }
            }
        }
        // every rep after the first of a key shares it with an earlier rep of the step
        SORT_PARALLEL(keyReps.begin(), keyReps.end());
        repSharesKey.assign(repCount, 0);
#pragma omp parallel for schedule(static)
        for (size_t i = 1; i < keyReps.size(); i++) {
            if (keyReps[i].first == keyReps[i - 1].first) {
                __atomic_store_n(&repSharesKey[keyReps[i].second], 1, __ATOMIC_RELAXED);
            }
        }
        // no earlier rep of the step takes a key of the others, so they take all their free members at once
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 1024) reduction(+:memberCount)
            for (size_t rep = 0; rep < repCount; rep++) {
                if (repKeyOffsets[rep + 1] == repKeyOffsets[rep] || repSharesKey[rep]) {
                    continue;
                }
                setAssigned(pairs[repOffsets[rep]].pair.repKey, threadAssignedKeys[thread_idx]);
                for (size_t i = repOffsets[rep]; i < repOffsets[rep + 1]; i++) {
                    if (pairIsFree[i]) {
                        setAssigned(pairs[i].pair.memberKey, threadAssignedKeys[thread_idx]);
                        __sync_fetch_and_or(&nodeTakenBits[pairs[i].node][pairs[i].takenBit / 64], 1ULL << (pairs[i].takenBit % 64));
                        memberCount++;
                    }
                }
            }
        }
        // the reps that share a key with an earlier one take in key order, and a later rep that does not share has no key of theirs
        for (size_t rep = 0; rep < repCount; rep++) {
            const DBKeyType repKey = pairs[repOffsets[rep]].pair.repKey;
            if (repKeyOffsets[rep + 1] == repKeyOffsets[rep] || repSharesKey[rep] == 0 || isAssigned(repKey)) {
                continue;
            }
            for (size_t i = repOffsets[rep]; i < repOffsets[rep + 1]; i++) {
                if (i + PREFETCH_PAIRS < pairCount) {
                    __builtin_prefetch(&assigned[pairs[i + PREFETCH_PAIRS].pair.memberKey / 64]);
                }
                const DBKeyType memberKey = pairs[i].pair.memberKey;
                if (pairIsFree[i] && isAssigned(memberKey) == false) {
                    setAssigned(repKey, threadAssignedKeys[0]);
                    setAssigned(memberKey, threadAssignedKeys[0]);
                    nodeTakenBits[pairs[i].node][pairs[i].takenBit / 64] |= 1ULL << (pairs[i].takenBit % 64);
                    memberCount++;
                }
            }
        }
        const size_t walkedBefore = walkedKey;
        walkedKey = stepRepKey;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            appendTakenBits(nodeTaken[node], nodeTakenBits[node].data(), nodeTakenCounts[node], words, taken + "." + SSTR(node));
        }

        // the keys of the step by block node, after those in the ring of every node, [walkedKey, walkedBefore + lagKeys), and the entering ring words
        const size_t sliceFirstWord = std::min(assignedWords, (walkedBefore + lagKeys) / 64);
        const size_t sliceWords = std::min(assignedWords, (walkedKey + lagKeys + 63) / 64) - sliceFirstWord;
#pragma omp parallel for schedule(static)
        for (int thread = 0; thread < threads; thread++) {
            std::vector<size_t> &sectionCounts = threadSectionCounts[thread];
            sectionCounts.assign(MMseqsNode::count + 1, 0);
            for (size_t i = 0; i < threadAssignedKeys[thread].size(); i++) {
                const DBKeyType key = threadAssignedKeys[thread][i];
                sectionCounts[getAlignNode(blockKeys, key) + 1]++;
                sectionCounts[0] += (key >= walkedKey && key < walkedBefore + lagKeys);
            }
        }
        // a count becomes where the keys of a thread start in its section, the ring section first
        sectionOffsets[0] = 0;
        for (size_t section = 0; section <= MMseqsNode::count; section++) {
            sectionOffsets[section + 1] = sectionOffsets[section];
            for (int thread = 0; thread < threads; thread++) {
                const size_t count = threadSectionCounts[thread][section];
                threadSectionCounts[thread][section] = sectionOffsets[section + 1];
                sectionOffsets[section + 1] += count;
            }
        }
        const size_t headerWords = 4 + MMseqsNode::count + 1;
        record.resize(headerWords + sliceWords + sectionOffsets[MMseqsNode::count + 1]);
        record[0] = walkedKey;
        record[1] = sliceFirstWord;
        record[2] = sliceWords;
        record[3] = sectionOffsets[1];
        for (size_t node = 0; node <= MMseqsNode::count; node++) {
            record[4 + node] = sectionOffsets[node + 1] - sectionOffsets[1];
        }
        memcpy(record.data() + headerWords, assigned + sliceFirstWord, sliceWords * sizeof(uint64_t));
#pragma omp parallel for schedule(static)
        for (int thread = 0; thread < threads; thread++) {
            std::vector<size_t> &sectionPositions = threadSectionCounts[thread];
            for (size_t i = 0; i < threadAssignedKeys[thread].size(); i++) {
                const DBKeyType key = threadAssignedKeys[thread][i];
                record[headerWords + sliceWords + sectionPositions[getAlignNode(blockKeys, key) + 1]++] = key;
                if (key >= walkedKey && key < walkedBefore + lagKeys) {
                    record[headerWords + sliceWords + sectionPositions[0]++] = key;
                }
            }
            threadAssignedKeys[thread].clear();
        }
        if (fwrite(record.data(), sizeof(uint64_t), record.size(), changesHandle) != record.size()) {
            Debug(Debug::ERROR) << "Cannot write " << assignedChanges << "\n";
            EXIT(EXIT_FAILURE);
        }
        changesBytes += record.size() * sizeof(uint64_t);
        MMseqsNode::publishLog(changesHandle, assignedChanges, walkedKey, changesBytes);
    }
    // every node is done, so it reads neither the log of changes nor writes its log of accepted pairs any more
    if (fclose(changesHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << assignedChanges << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::remove(assignedChanges.c_str());
    FileUtil::remove((assignedChanges + ".progress").c_str());
    // every node reads the assigned bits of its keys: a rep that took members is in its clusters, a member is in no cluster of its own
    FILE *assignedHandle = FileUtil::openAndDelete(assignedKeys.c_str(), "w");
    if (fwrite(assigned, sizeof(uint64_t), assignedWords, assignedHandle) != assignedWords) {
        Debug(Debug::ERROR) << "Cannot write " << assignedKeys << "\n";
        EXIT(EXIT_FAILURE);
    }
    MMseqsNode::publishLog(assignedHandle, assignedKeys, dbSize, assignedWords * sizeof(uint64_t));
    if (fclose(assignedHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << assignedKeys << "\n";
        EXIT(EXIT_FAILURE);
    }
    delete[] assigned;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        FileUtil::remove((accepted + "." + SSTR(node) + ".progress").c_str());
        const std::string takenFile = taken + "." + SSTR(node);
        if (nodeTaken[node].carryBits > 0 && fwrite(&nodeTaken[node].carry, sizeof(uint64_t), 1, nodeTaken[node].file) != 1) {
            Debug(Debug::ERROR) << "Cannot write " << takenFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        MMseqsNode::publishLog(nodeTaken[node].file, takenFile, dbSize, nodeTaken[node].bytes + ((nodeTaken[node].carryBits > 0) ? sizeof(uint64_t) : 0));
        if (fclose(nodeTaken[node].file) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << takenFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return memberCount;
}

// --createdb-mode 3: the nodes of the targets align the reps that the pref partitions send, and the master clusters the accepted pairs greedily
static int align2clustLinfinity(Parameters &par) {
    Timer timer;
    // the dbtype of an earlier result must not mark this one as done before it is
    if (MMseqsNode::isMaster() && FileUtil::fileExists((par.db3 + ".dbtype").c_str())) {
        FileUtil::remove((par.db3 + ".dbtype").c_str());
    }
    DBReader<DBKeyType> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
    seqDbr.open(DBReader<DBKeyType>::NOSORT);
    const std::vector<size_t> nodeKeyOffsets = seqDbr.getFileKeyOffsets(MMseqsNode::count);
    const size_t dbSize = seqDbr.getSize();
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    // the master clusters with half of its memory and a few more threads, which mostly wait for pairs, while it aligns like every node
    const size_t memoryLimit = Util::computeMemory(0) / (MMseqsNode::isMaster() ? 2 : 1);
    const int clusterThreads = std::min(CLUSTER_THREADS, std::max(1, par.threads / 4));

    // a clusthashfast subset: its representatives are the bits of every node and the sixth column of its index, the residues of the e-values
    const bool hasRepresentativeBits = FileUtil::fileExists((par.db1 + ".bits." + SSTR(MMseqsNode::rank)).c_str());
    size_t representatives = 0;
    size_t representativeResidues = 0;
    for (size_t range = 0; range < seqDbr.getIndexSize(); range++) {
        representatives += seqDbr.getRangeRepresentativeCount(range);
        representativeResidues += seqDbr.getRangeRepresentativeCount(range) * (index[range].length - 2);
    }

    // every node takes the same splits and ring of the smallest memory limit, and fails its memory checks together with the others
    const size_t smallestMemoryLimit = MMseqsNode::getSmallestMemoryLimit(par.db3 + ".memory", memoryLimit);
    const size_t assignedWords = (dbSize + 63) / 64;
    if (assignedWords * sizeof(uint64_t) > smallestMemoryLimit / 2) {
        Debug(Debug::ERROR) << "The master needs " << assignedWords * sizeof(uint64_t) << " bytes for a bit of each of the " << dbSize << " keys, more than half of the smallest memory limit of "
                            << smallestMemoryLimit << " bytes. Run on nodes with more memory\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t splitBytes = std::max<size_t>(1, smallestMemoryLimit / 8);
    const size_t lagKeys = std::max<size_t>(64, smallestMemoryLimit / 256 * 8 / 64 * 64);
    std::vector<size_t> splitKeyOffsets;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        const std::vector<size_t> nodeSplitKeyOffsets = seqDbr.getSplitKeyOffsets(nodeKeyOffsets[node], nodeKeyOffsets[node + 1], splitBytes, 0);
        splitKeyOffsets.insert(splitKeyOffsets.end(), nodeSplitKeyOffsets.begin(), nodeSplitKeyOffsets.end() - 1);
    }
    splitKeyOffsets.push_back(dbSize);

    // the master cuts the targets into blocks of keys, which the nodes take in turn (see BLOCK_BYTES)
    const std::string blocksFile = par.db3 + ".blocks";
    if (MMseqsNode::isMaster()) {
        struct stat sb;
        if (stat(seqDbr.getDataFileNames()[0].c_str(), &sb) != 0) {
            Debug(Debug::ERROR) << "Cannot read " << seqDbr.getDataFileNames()[0] << "\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t entryBytes = std::max<size_t>(1, seqDbr.getDataSize() / std::max<size_t>(1, dbSize));
        const size_t fileSystemKeys = (static_cast<size_t>(sb.st_blksize) + entryBytes - 1) / entryBytes;
        std::vector<size_t> distances;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            DBReader<DBKeyType> prefDb((par.db2 + "." + SSTR(node)).c_str(), (par.db2 + ".index." + SSTR(node)).c_str(), 1,
                                       DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
            prefDb.open(DBReader<DBKeyType>::NOSORT);
            const size_t entryCount = std::min(prefDb.getSize(), BLOCK_SAMPLE_ENTRIES);
            if (entryCount > 0) {
                prefDb.loadData(0, entryCount);
            }
            for (size_t entry = 0; entry < entryCount; entry++) {
                char *data = prefDb.getData(entry, 0);
                const char *dataEnd = data + prefDb.getEntryLen(entry);
                const DBKeyType repKey = Util::fast_atoi<DBKeyType>(data);
                for (data = Util::skipLine(data); data < dataEnd && *data != '\0'; data = Util::skipLine(data)) {
                    const DBKeyType targetKey = Util::fast_atoi<DBKeyType>(data);
                    distances.push_back((targetKey > repKey) ? targetKey - repKey : repKey - targetKey);
                }
            }
            prefDb.close();
        }
        size_t blockKeys = (std::max<size_t>(BLOCK_BYTES, sb.st_blksize) + entryBytes - 1) / entryBytes;
        if (distances.empty() == false) {
            std::nth_element(distances.begin(), distances.begin() + distances.size() / 2, distances.end());
            blockKeys = std::max(fileSystemKeys, std::min(blockKeys, distances[distances.size() / 2] / (BLOCKS_PER_DISTANCE * MMseqsNode::count)));
        }
        blockKeys = std::max<size_t>(64, (blockKeys + 63) / 64 * 64);
        MMseqsNode::publishRunFile(blocksFile, SSTR(blockKeys) + "\n");
    }
    const size_t blockKeys = Util::fast_atoi<size_t>(MMseqsNode::waitForRunFile(blocksFile).c_str());
    const size_t ownBlocks = ((dbSize + blockKeys - 1) / blockKeys + MMseqsNode::count - 1) / MMseqsNode::count;
    if (ownBlocks * (blockKeys / 64) * sizeof(uint64_t) > smallestMemoryLimit / 8) {
        Debug(Debug::ERROR) << "An aligning node needs " << ownBlocks * (blockKeys / 64) * sizeof(uint64_t) << " bytes for a bit of each of the " << ownBlocks * blockKeys
                            << " keys of its blocks, more than an eighth of the smallest memory limit of " << smallestMemoryLimit << " bytes. Run on more nodes\n";
        EXIT(EXIT_FAILURE);
    }

    // a log of an earlier run must not be read as one of this run, so its progress is removed before the exchange barrier
    const std::string accepted = par.db3 + "_swap_1";
    const std::string assignedChanges = par.db3 + "_assigned";
    const std::string taken = par.db3 + "_taken";
    const std::string assignedKeys = par.db3 + "_assigned.bits";
    std::vector<std::string> progressFiles;
    progressFiles.push_back(accepted + "." + SSTR(MMseqsNode::rank) + ".progress");
    progressFiles.push_back(taken + "." + SSTR(MMseqsNode::rank) + ".progress");
    if (MMseqsNode::isMaster()) {
        progressFiles.push_back(assignedChanges + ".progress");
        progressFiles.push_back(assignedKeys + ".progress");
    }
    for (size_t i = 0; i < progressFiles.size(); i++) {
        if (FileUtil::fileExists(progressFiles[i].c_str())) {
            FileUtil::remove(progressFiles[i].c_str());
        }
    }

    const std::string repTargets = par.db3 + "_swap_0";
    const std::string repTargetsIndex = sendRepTargets(par, seqDbr, nodeKeyOffsets, splitKeyOffsets, blockKeys, memoryLimit, repTargets);
    Debug(Debug::INFO) << "Time for sending the rep targets: " << timer.lap() << "\n";
    Timer phaseTimer;

    // the master clusters the accepted pairs while the nodes align
    size_t memberCount = 0;
    std::thread clusterThread;
    if (MMseqsNode::isMaster()) {
        clusterThread = std::thread([&]() {
#ifdef OPENMP
            omp_set_num_threads(clusterThreads);
#endif
            memberCount = clusterAcceptedPairs(clusterThreads, dbSize, memoryLimit, blockKeys, lagKeys, accepted, assignedChanges, taken, assignedKeys);
        });
    }
    alignRepTargets(par, seqDbr, blockKeys, lagKeys, memoryLimit, representativeResidues, repTargetsIndex, repTargets, accepted, assignedChanges);

    // every node sends the taken pairs of its log to the nodes of their reps once the master has clustered all
    size_t takenKey;
    size_t takenBytes;
    MMseqsNode::waitForLog(taken + "." + SSTR(MMseqsNode::rank), takenKey, takenBytes);
    std::vector<uint64_t> takenBits(takenBytes / sizeof(uint64_t));
    MMseqsNode::readLog(taken + "." + SSTR(MMseqsNode::rank), 0, takenBits.size() * sizeof(uint64_t), reinterpret_cast<char *>(takenBits.data()));
    FileUtil::remove((taken + "." + SSTR(MMseqsNode::rank)).c_str());
    FileUtil::remove((taken + "." + SSTR(MMseqsNode::rank) + ".progress").c_str());
    if (MMseqsNode::isMaster()) {
        clusterThread.join();
        Debug(Debug::INFO) << "Number of clusters: " << representatives - memberCount << "\n";
    }
    Debug(Debug::INFO) << "Time for aligning: " << phaseTimer.lap() << "\n";
    const size_t nodeFirstKey = nodeKeyOffsets[MMseqsNode::rank];
    const size_t nodeKeys = nodeKeyOffsets[MMseqsNode::rank + 1] - nodeFirstKey;
    const size_t nodeWords = (nodeKeys + 63) / 64;
    const std::string clusters = par.db3 + "_swap_2";
    MMseqsNode::ClusterPairWriter clustersWriter(clusters, std::max<size_t>(64, memoryLimit / 16 / sizeof(ClusterPair)), ClusterPair::compareByRep, nodeKeyOffsets);
    const std::string acceptedFile = accepted + "." + SSTR(MMseqsNode::rank);
    FILE *acceptedHandle = FileUtil::openFileOrDie(acceptedFile.c_str(), "r", true);
    // a chunk of pairs takes whole words of taken bits and fits into the writer
    const size_t chunkPairs = clustersWriter.pairCapacity / 64 * 64;
    std::vector<ClusterPair> acceptedPairs(chunkPairs);
    for (size_t pairCount, firstWord = 0; (pairCount = fread(acceptedPairs.data(), sizeof(ClusterPair), chunkPairs, acceptedHandle)) > 0; firstWord += chunkPairs / 64) {
        const size_t wordCount = (pairCount + 63) / 64;
        if (firstWord + wordCount > takenBits.size()) {
            Debug(Debug::ERROR) << "Cannot read " << taken << "." << MMseqsNode::rank << "\n";
            EXIT(EXIT_FAILURE);
        }
        const uint64_t *takenWords = takenBits.data() + firstWord;
        size_t takenCount = 0;
#pragma omp parallel for schedule(static) reduction(+:takenCount)
        for (size_t word = 0; word < wordCount; word++) {
            takenCount += __builtin_popcountll(takenWords[word]);
        }
        if (clustersWriter.pairCount + takenCount > clustersWriter.pairCapacity) {
            clustersWriter.flush();
        }
#pragma omp parallel for schedule(static)
        for (size_t word = 0; word < wordCount; word++) {
            uint64_t bits = takenWords[word];
            if (bits == 0) {
                continue;
            }
            size_t pairOffset = __sync_fetch_and_add(&clustersWriter.pairCount, __builtin_popcountll(bits));
            for (; bits != 0; bits &= bits - 1) {
                const ClusterPair &pair = acceptedPairs[word * 64 + __builtin_ctzll(bits)];
                clustersWriter.pairs[pairOffset++] = pair;
            }
        }
    }
    if (fclose(acceptedHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << acceptedFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::remove(acceptedFile.c_str());
    std::vector<uint64_t>().swap(takenBits);
    const std::string clustersIndex = clustersWriter.close();
    // every node has aligned the records of every node
    FileUtil::remove((repTargets + "." + SSTR(MMseqsNode::rank)).c_str());
    if (MMseqsNode::isMaster()) {
        FileUtil::remove((repTargets + ".index").c_str());
    }

    // every node writes the clusters of its reps and its unassigned representatives on their own; readKeyBits fills a word more than the keys
    size_t assignedKey;
    size_t assignedBytes;
    MMseqsNode::waitForLog(assignedKeys, assignedKey, assignedBytes);
    FILE *assignedHandle = FileUtil::openFileOrDie(assignedKeys.c_str(), "r", true);
    std::vector<uint64_t> assignedBits(nodeWords + 1, 0);
    MMseqsNode::readKeyBits(fileno(assignedHandle), assignedKeys, 0, 0, nodeFirstKey, nodeFirstKey + nodeKeys, nodeFirstKey, assignedBits.data());
    if (fclose(assignedHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << assignedKeys << "\n";
        EXIT(EXIT_FAILURE);
    }
    uint64_t *singletonBits = NULL;
    if (hasRepresentativeBits) {
        singletonBits = MMseqsNode::readRepresentativeBits(par.db1, nodeKeyOffsets, nodeFirstKey, nodeKeys);
    } else {
        singletonBits = new(std::nothrow) uint64_t[nodeWords];
        Util::checkAllocation(singletonBits, "Can not allocate singletonBits memory in align2clust");
        memset(singletonBits, 0xFF, nodeWords * sizeof(uint64_t));
    }
    for (size_t word = 0; word < nodeWords; word++) {
        singletonBits[word] &= ~assignedBits[word];
    }
    std::vector<uint64_t>().swap(assignedBits);
    MMseqsNode::writeClusterPairs(clustersIndex, clusters, par.db3, singletonBits, nodeFirstKey, nodeKeys);
    delete[] singletonBits;
    MMseqsNode::removeExchanges(par.db3, clusters, std::vector<std::string>(1, clusters));
    if (MMseqsNode::isMaster()) {
        FileUtil::remove((par.db3 + ".memory").c_str());
        FileUtil::remove(blocksFile.c_str());
        FileUtil::remove(assignedKeys.c_str());
        FileUtil::remove((assignedKeys + ".progress").c_str());
        DBWriter::writeDbtypeFile(par.db3.c_str(), DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_CLUSTER_RES, Parameters::DBTYPE_EXTENDED_LINFINITY), false);
    }
    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";
    return EXIT_SUCCESS;
}

int align2clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(FileUtil::parseDbType(par.db1.c_str())) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isLinfinity) {
        MMseqsNode::init();
        if (par.clusteringMode != Parameters::GREEDY && par.clusteringMode != Parameters::GREEDY_MEM) {
            Debug(Debug::ERROR) << "align2clust supports only --cluster-mode " << Parameters::GREEDY << " and " << Parameters::GREEDY_MEM
                                << " for --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.filterCluDBFile.empty() == false || par.filterSeqDBFile.empty() == false) {
            Debug(Debug::ERROR) << "--filter-cludb-file and --filter-seqdb-file are not supported for --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            EXIT(EXIT_FAILURE);
        }
        const int prefDbtype = FileUtil::parseDbType(par.db2.c_str());
        if (Parameters::isEqualDbtype(prefDbtype, Parameters::DBTYPE_PREFILTER_RES) == false
            || (DBReader<DBKeyType>::getExtendedDbtype(prefDbtype) & Parameters::DBTYPE_EXTENDED_LINFINITY) == 0) {
            Debug(Debug::ERROR) << par.db2 << " is not a --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " prefilter result\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.compressed) {
            Debug(Debug::WARNING) << "--compressed cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            Debug(Debug::WARNING) << "We recompute with --compressed 0\n";
            par.compressed = 0;
        }
        if (par.includeAlignFiles) {
            Debug(Debug::WARNING) << "--include-align-files cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            Debug(Debug::WARNING) << "We recompute with --include-align-files 0\n";
            par.includeAlignFiles = false;
        }
        par.diagonalMode = Parameters::DIAGONAL_MODE_UNGAPPED;
        return align2clustLinfinity(par);
    }

    Timer timer;
    timer.reset();
    
    DBReader<DBKeyType> alnDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    alnDbr.open(DBReader<DBKeyType>::LINEAR_ACCCESS);
    int dbtype =  Parameters::DBTYPE_CLUSTER_RES;

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, dbtype);
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

    resultWriter.close();
    if (alnWriter != nullptr) {
        alnWriter->close();
        delete alnWriter;
    }
    alnDbr.close();

    return status;
}
