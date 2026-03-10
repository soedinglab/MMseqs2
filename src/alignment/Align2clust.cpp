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

#ifdef OPENMP
#include <omp.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#endif

#define MAX_SIZE 4096 //change
#define MIN_SIZE 32

struct ClusterResult {
    unsigned int sequenceIdx;
    unsigned int representativeId;
    unsigned int prefSize;
    std::vector<unsigned int> memberIds;

    bool revalidated = false;
    
    bool operator<(const ClusterResult& other) const {
        if (memberIds.size() > other.memberIds.size()) {
            return true;
        }
        if (other.memberIds.size() > memberIds.size()) {
            return false;
        }
        if (representativeId < other.representativeId) {
            return true;
        }
        if (other.representativeId < representativeId) {
            return false;            
        }
        return true;
    }
};

struct PrefInfo {
    size_t id;
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

struct GreedyComparator {
    bool operator()(const ClusterResult& a, const ClusterResult& b) const {
        return a.sequenceIdx > b.sequenceIdx;
    }
};

struct SetCoverComparator {
    bool operator()(const ClusterResult& a, const ClusterResult& b) const {
        if (a.memberIds.size() < b.memberIds.size()) {
            return true;
        }
        if (b.memberIds.size() < a.memberIds.size()) {
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
static std::priority_queue<ClusterResult, std::vector<ClusterResult>, GreedyComparator> clusterResultQueue;
static std::priority_queue<ClusterResult, std::vector<ClusterResult>, SetCoverComparator> setCoverReadyQueue;

static std::unique_ptr<std::vector<ClusterResult>> elements;
static unsigned int currentProcessPosition = 0; 
static unsigned int currentPrefSize = 0; 
static bool allCalculationsDone = false;

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

static void writeData(DBWriter *dbWriter, const std::pair<unsigned int, unsigned int> * results, size_t dbSize) {
    std::string resultString;
    resultString.reserve(1024 * 1024 * 1024);
    char buffer[32];
    unsigned int previousRepresentativeKey = UINT_MAX;
    
    for (size_t i = 0; i < dbSize; i++) {
        unsigned int currentRepresentativeKey = results[i].first;
        
        if (previousRepresentativeKey != currentRepresentativeKey) {
            if (previousRepresentativeKey != UINT_MAX) {
                dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey);
            }
            resultString.clear();
            char *outPos = Itoa::u32toa_sse2(currentRepresentativeKey, buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }
        
        unsigned int memberKey = results[i].second;
        if (memberKey != currentRepresentativeKey) {
            char *outPos = Itoa::u32toa_sse2(memberKey, buffer);
            resultString.append(buffer, (outPos - buffer - 1));
            resultString.push_back('\n');
        }
        
        previousRepresentativeKey = currentRepresentativeKey;
    }
    
    if (previousRepresentativeKey != UINT_MAX) {
        dbWriter->writeData(resultString.c_str(), resultString.length(), previousRepresentativeKey);
    }
}

static void (*clusterThreadFunc)(unsigned int*) = nullptr;

void clusterThreadFuncSetcover(unsigned int* assignedCluster) {
    std::vector<unsigned int> valid;
    valid.reserve(1024);
    
    while (true) {
        std::unique_lock<std::mutex> lock(clusterMutex);
        
        clusterCondition.wait(lock, [] {
            return ((!clusterResultQueue.empty() &&
                     clusterResultQueue.top().sequenceIdx == currentProcessPosition) ||
                    allCalculationsDone);
        });
        
        // 1) clusterResultQueue → setCoverReadyQueue
        while (!clusterResultQueue.empty() &&
               clusterResultQueue.top().sequenceIdx == currentProcessPosition) {
            ClusterResult result = clusterResultQueue.top();
            clusterResultQueue.pop();
            currentProcessPosition++;
            currentPrefSize = result.prefSize;
            
            if (result.memberIds.size() > 1) {
                setCoverReadyQueue.push(std::move(result));
            }
        }
        
        // 2) handle setCoverReadyQueue 
        while (!setCoverReadyQueue.empty() &&
               (allCalculationsDone ||
                setCoverReadyQueue.top().memberIds.size() > currentPrefSize)) {
            
            ClusterResult res = setCoverReadyQueue.top();
            setCoverReadyQueue.pop();
            
            if (assignedCluster[res.representativeId] != UINT_MAX) {
                continue;
            }
            
            valid.clear();
            size_t originalSize = res.memberIds.size();
            
            for (unsigned int mem : res.memberIds) {
                if (assignedCluster[mem] == UINT_MAX) {
                    valid.push_back(mem);
                }
            }
            
            if (valid.size() <= 1) {
                continue;
            }
            
            if (valid.size() != originalSize) {
                res.memberIds.swap(valid);
                setCoverReadyQueue.push(std::move(res));
                continue;
            }
            
            for (unsigned int mem : res.memberIds) {
                assignedCluster[mem] = res.representativeId;
            }
        }
        
        if (allCalculationsDone &&
            clusterResultQueue.empty() &&
            setCoverReadyQueue.empty()) {
            break;
        }
    }
}

void clusterThreadFuncGreedy(unsigned int* assignedCluster) {
    while (true) {
        std::unique_lock<std::mutex> lock(clusterMutex);
        
        clusterCondition.wait(lock, [] { 
            return (!clusterResultQueue.empty() && 
                    clusterResultQueue.top().sequenceIdx == currentProcessPosition) 
                   || allCalculationsDone; 
        });

        if (allCalculationsDone && clusterResultQueue.empty()) {
            break;
        }

        while (!clusterResultQueue.empty() && 
               clusterResultQueue.top().sequenceIdx == currentProcessPosition) {
            ClusterResult result = clusterResultQueue.top();
            clusterResultQueue.pop();
            currentProcessPosition++;
            
            if (assignedCluster[result.representativeId] != UINT_MAX) {
                continue;  
            }
                        
            std::vector<unsigned int> validMemberIds;
            validMemberIds.reserve(result.memberIds.size());
            for (unsigned int memberId : result.memberIds) {
                if (assignedCluster[memberId] == UINT_MAX) {
                    validMemberIds.push_back(memberId);
                }
            }
            
            if (validMemberIds.size() <= 1) {
                continue;
            }
            
            for (unsigned int memberId : validMemberIds) {
                assignedCluster[memberId] = result.representativeId;
            }
        }
    }
}

int doAlign2clust(Parameters &par, DBWriter &resultWriter, DBReader<unsigned int> &alnDbr) {
    DBReader<unsigned int> *seqDbr = new DBReader<unsigned int>(
        par.db1.c_str(), par.db1Index.c_str(), par.threads, 
        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
    );
    seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);
 
    DBReader<unsigned int> *cluDbr = nullptr;
    DBReader<unsigned int> *cluSeqDbr = nullptr;
    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false) {
        std::string cluIndex = par.filterCluDBFile + ".index";
        cluDbr = new DBReader<unsigned int>(
            par.filterCluDBFile.c_str(), cluIndex.c_str(), par.threads, 
            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
        );
        cluDbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        std::string cluSeqIndex = par.filterSeqDBFile + ".index";
        cluSeqDbr = new DBReader<unsigned int>(
            par.filterSeqDBFile.c_str(), cluSeqIndex.c_str(), par.threads, 
            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX
        );
            
        cluSeqDbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
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
    
    unsigned int *assignedCluster = new(std::nothrow) unsigned int[dbSize];
    Util::checkAllocation(assignedCluster, "Can not allocate assignedCluster memory in Align2Clust");
    std::fill_n(assignedCluster, dbSize, UINT_MAX);

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
                const unsigned int clusterId = seqDbr->getDbKey(i);
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
#pragma omp parallel
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = (unsigned int) omp_get_thread_num();
#endif
        Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, seqDbr->getMaxSeqLen(), subMat, &evaluer, 
                       par.compBiasCorrection, par.compBiasCorrectionScale, 
                       par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 
                       0.0, par.zdrop);
        Sequence query(seqDbr->getMaxSeqLen(), Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        Sequence target(seqDbr->getMaxSeqLen(), Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        Sequence element(seqDbr->getMaxSeqLen(), Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false, par.compBiasCorrection);
        BlockAligner blockAligner(Parameters::DBTYPE_AMINO_ACIDS, seqDbr->getMaxSeqLen(), subMat, &fastMatrix, 
                                 &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, 
                                 -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid());
        char buffer[1024 + 32768 * 4];
        std::vector<std::pair<unsigned int, unsigned short>> targetsWithDiagonal;
        targetsWithDiagonal.reserve(1000);

#pragma omp for schedule(dynamic, 1) nowait
        for (size_t i = 0; i < endRange; i++) {
            progress.updateProgress();
            ClusterResult clusterResult;
            clusterResult.sequenceIdx = i;
            targetsWithDiagonal.clear();
            
            size_t representativeId;
            unsigned int queryKey;
            
            if (mode == Parameters::SET_COVER) {
                representativeId = prefRepSizePair[i].id;
                queryKey = seqDbr->getDbKey(representativeId);
            } else { // GREEDY || GREEDY_MEM
                queryKey = seqDbr->getDbKey(i);
                representativeId = seqDbr->getId(queryKey);
            }
            
            const size_t alignmentId = alnDbr.getId(queryKey);
            char *alignmentData = alnDbr.getData(alignmentId, threadIdx);
            clusterResult.representativeId = representativeId;
            size_t queryId = representativeId;
            char *querySequence = seqDbr->getData(queryId, threadIdx);
            size_t queryLength = seqDbr->getSeqLen(queryId);
            query.mapSequence(queryId, queryKey, querySequence, queryLength);
            blockAligner.initQuery(&query);
            matcher.initQuery(&query);
            
            if (assignedCluster[representativeId] != UINT_MAX) {
                {
                    std::lock_guard<std::mutex> lock(clusterMutex);
                    clusterResultQueue.push(std::move(clusterResult));
                }
                continue;
            }

            unsigned int prefSize = 0;
            while (*alignmentData != '\0') {
                hit_t hit = QueryMatcher::parsePrefilterHit(alignmentData);
                const size_t targetId = seqDbr->getId(hit.seqId);
                if (assignedCluster[targetId] == UINT_MAX) { 
                        targetsWithDiagonal.push_back(std::make_pair(hit.seqId, hit.diagonal));
                }
                alignmentData = Util::skipLine(alignmentData);
                prefSize++;
            }
            clusterResult.prefSize = prefSize;

            for (size_t targetIdx = 0; targetIdx < targetsWithDiagonal.size(); targetIdx++) {
                const unsigned int targetKey = targetsWithDiagonal[targetIdx].first;
                const unsigned short diagonal = targetsWithDiagonal[targetIdx].second;
                const size_t targetId = seqDbr->getId(targetKey);
                
                const bool isIdentity = (queryKey == targetKey);
                if (isIdentity) {
                    clusterResult.memberIds.push_back(queryId);
                    continue;
                }
                
                char *targetSequence = seqDbr->getData(targetId, threadIdx);
                size_t targetLength = seqDbr->getSeqLen(targetId);
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
                    seqId = Util::computeSeqId(par.seqIdMode, identicalCount, query.L, target.L, ungappedAlignment.alnLen);
                }
                
                char *bufferEnd = Itoa::i32toa_sse2(ungappedAlignment.alnLen, buffer);
                size_t bufferLen = bufferEnd - buffer;
                std::string backtrace = "";
                if (par.addBacktrace) {
                    backtrace = std::string(buffer, bufferLen - 1);
                    backtrace.push_back('M');
                }
                
                if (isIdentity) {
                    ungappedAlignment.qcov = 1.0f;
                    ungappedAlignment.tcov = 1.0f;
                    seqId = 1.0f;
                }
                
                bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                //ugly temporary gyuri
                if (assignedCluster[targetId] != UINT_MAX) continue;

                if (isIdentity || (hasAlnLen && hasCoverage && hasSeqId && hasEvalue)) {
                    Matcher::result_t result = Matcher::result_t(
                        targetKey, ungappedAlignment.bitScore, ungappedAlignment.qcov, ungappedAlignment.tcov, 
                        seqId, ungappedAlignment.eval, ungappedAlignment.alnLen,
                        ungappedAlignment.qStart, ungappedAlignment.qEnd, query.L, 
                        ungappedAlignment.tStart, ungappedAlignment.tEnd, target.L, backtrace
                    );
                    if (assignedCluster[targetId] != UINT_MAX) continue;
                    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                        // check all the member from filtering file
                        const size_t cluId = cluDbr->getId(targetKey);
                        char *cluData = cluDbr->getData(cluId, threadIdx);
                        const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                        unsigned int numClu = Util::countLines(cluData, cluDataSize);
                        bool allpass = true;
                        char buffer[1024];
                        if (numClu > 1) { // if not singleton
                            while (*cluData != '\0') {
                                Util::parseKey(cluData, buffer);

                                const unsigned int elementKey = (unsigned int) strtoul(buffer, NULL, 10);
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
                                    allpass = false;
                                    break;
                                    // float currentScorePerCol = static_cast<float>(elementUngappedAlignment.score) / static_cast<float>(elementUngappedAlignment.diagonalLen);
                                    // if (currentScorePerCol < scorePerColThreshold) {
                                    //     allpass = false;
                                    //     break;
                                    // } 
                                    // else {
                                    //     // 3. gapped alignment
                                    //     Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(concatedDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                    //                                         swMode, par.seqIdMode, false);
                                    //     if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                    //         allpass = false;
                                    //         break;
                                    //     }
                                    // }   
                                }
                                cluData = Util::skipLine(cluData);
                            }
                        }
                        if (allpass == false) {
                            continue;
                        }
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

                //ugly temporary gyuri
                if (assignedCluster[targetId] != UINT_MAX) continue;

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
                        if (assignedCluster[targetId] != UINT_MAX) continue;
                        if (par.filteringFile.empty()==false){
                            // check all the member from filtering file
                            const size_t cluId = cluDbr->getId(targetKey);
                            char *cluData = cluDbr->getData(cluId, threadIdx);
                            const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                            unsigned int numClu = Util::countLines(cluData, cluDataSize);
                            bool allpass = true;
                            char buffer[1024];
                            if (numClu > 1) { // if not singleton
                                while (*cluData != '\0') {
                                    Util::parseKey(cluData, buffer);
                                    const unsigned int elementKey = (unsigned int) strtoul(buffer, NULL, 10);
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
                                        allpass = false;
                                        break;
                                        // float currentScorePerCol = static_cast<float>(elementUngappedAlignment.score) / static_cast<float>(elementUngappedAlignment.diagonalLen);
                                        // if (currentScorePerCol < scorePerColThreshold) {
                                        //     allpass = false;
                                        //     break;
                                        // } 
                                        // else {
                                        //     // 3. gapped alignment
                                        //     Matcher::result_t res_element = matcher.getSWResult(&element, static_cast<int>(concatedDiagonal), false, par.covMode, par.covThr, par.evalThr,
                                        //                                         swMode, par.seqIdMode, false);
                                        //     if (Alignment::checkCriteria(res_element, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr) == false) {
                                        //         allpass = false;
                                        //         break;
                                        //     }
                                        // }
                                    }
                                    cluData = Util::skipLine(cluData);
                                }
                            }
                            if (allpass == false) {
                                continue;
                            }
                        }
                        clusterResult.memberIds.push_back(targetId);
                    }
                }
            }
            
            bool isNextPosition = false;
            {
                std::lock_guard<std::mutex> lock(clusterMutex);
                if (clusterResult.sequenceIdx == currentProcessPosition) {
                    isNextPosition = true;
                }
                clusterResultQueue.push(std::move(clusterResult));
            }
            
            if (isNextPosition) {
                clusterCondition.notify_one();
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(clusterMutex);
        allCalculationsDone = true;
    }
    clusterCondition.notify_one();
    
    if (clusterThread.joinable()) {
        clusterThread.join(); 
    }

    for (size_t i = 0; i < dbSize; ++i) {
        if (assignedCluster[i] == UINT_MAX) {
            assignedCluster[i] = i;
        }
    }
    
    std::pair<unsigned int, unsigned int> *assignment = new std::pair<unsigned int, unsigned int>[dbSize];
    
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            if (assignedCluster[i] == UINT_MAX) {
                Debug(Debug::ERROR) << "There must be an error: " << i 
                                    << " is not assigned to a cluster\n";
                continue;
            }

            assignment[i].first = seqDbr->getDbKey(assignedCluster[i]);
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
    
    return 0;
}

int align2clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    
    Timer timer;
    timer.reset();
    
    DBReader<unsigned int> alnDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, 
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    int dbtype =  Parameters::DBTYPE_CLUSTER_RES;

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, dbtype);
    resultWriter.open();

    int status = doAlign2clust(par, resultWriter, alnDbr);
    
    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";
    
    resultWriter.close();
    alnDbr.close();
    
    return status;
}