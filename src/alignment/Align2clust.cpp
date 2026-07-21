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
#include "Align2ClustReducer.h"
#include "Align2ClustChunking.h"
#include "Align2ClustCheckpoint.h"
#include "MMseqsMPI.h"
#include <atomic>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

#define MIN_SIZE 32

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

// Serialize a single alignment record and append it to the per-representative buffer.
static void appendAlignmentResult(std::string &alnResultBuffer, char *lineBuffer, const Matcher::result_t &result, bool addBacktrace) {
    size_t len = Matcher::resultToBuffer(lineBuffer, result, addBacktrace);
    alnResultBuffer.append(lineBuffer, len);
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

// Env override for the MPI/checkpoint chunk size; 0 (unset/invalid) means use the default.
// Deliberately not derived from rank or thread count: chunk boundaries must stay identical
// across topology changes so a run can be resumed with a different node count (see
// Align2ClustChunking and Align2ClustCheckpoint).
static size_t getAlign2ClustChunkSizeFromEnv() {
    const size_t defaultChunkSize = 50000;
    const char *envValue = getenv("MMSEQS_ALIGN2CLUST_CHUNK_SIZE");
    if (envValue == nullptr || *envValue == '\0') {
        return defaultChunkSize;
    }

    char *end = nullptr;
    unsigned long long parsedValue = strtoull(envValue, &end, 10);
    if (end == envValue || *end != '\0' || parsedValue == 0) {
        Debug(Debug::WARNING) << "Ignoring invalid MMSEQS_ALIGN2CLUST_CHUNK_SIZE=" << envValue
                              << "; using default chunk size " << defaultChunkSize << "\n";
        return defaultChunkSize;
    }
    return static_cast<size_t>(parsedValue);
}

// Size the reorder buffer from free memory: subtract the resident per-sequence arrays,
// keep 10% headroom, divide by the worst-case result size. Capped by
// MMSEQS_ALIGN2CLUST_REORDER_LIMIT and by the number of results produced.
static size_t computeReorderCapacity(const Parameters &par, size_t dbSize, int mode, size_t resultCount) {
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);

    size_t fixedMemory = dbSize * sizeof(std::atomic<DBLocalId>);
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

// resultWriter is nullptr on non-master MPI ranks: only rank 0 performs the global
// reduction and writes the final clustering result (see the useDistributedGeneration
// branch below and the Align2ClustReducer class comment for why this cannot itself
// be decentralized).
int doAlign2clust(Parameters &par, DBWriter *resultWriter, DBReader<DBKeyType> &alnDbr, DBWriter *alnWriter) {
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

    int mode = par.clusteringMode;

    Align2ClustReducer::Mode reducerMode;
    if (mode == Parameters::SET_COVER) {
        reducerMode = Align2ClustReducer::SET_COVER;
        Debug(Debug::INFO) << "Using SET_COVER clustering mode\n";
    } else if (mode == Parameters::GREEDY || mode == Parameters::GREEDY_MEM) {
        reducerMode = Align2ClustReducer::GREEDY;
        Debug(Debug::INFO) << "Using GREEDY clustering mode\n";
    } else {
        Debug(Debug::ERROR) << "MMseqs2 align2clust doesn't support clustering mode: " << mode << "\n";
        delete[] fastMatrix.matrix;
        delete[] fastMatrix.matrixData;
        delete subMat;
        seqDbr->close();
        delete seqDbr;
        return EXIT_FAILURE;
    }

    const size_t endRange = (mode == Parameters::SET_COVER) ? dbSize : alnDbr.getSize();

    // Ring size = out-of-order window, sized from the memory budget (OOM-aware) and
    // capped by the result count. sequenceIdx runs over [0, endRange); every index
    // publishes exactly one result.
    const size_t reorderCapacityChosen = computeReorderCapacity(par, dbSize, mode, endRange);

    int mpiRank = 0;
    int mpiNumProc = 1;
#ifdef HAVE_MPI
    mpiRank = MMseqsMPI::rank;
    mpiNumProc = MMseqsMPI::numProc;
#endif

    // Deterministic, topology-independent chunking of [0, endRange) for MPI-distributed
    // candidate generation and checkpoint/resume. Chunk boundaries depend only on
    // endRange and chunkSize -- never on rank or thread count -- so a run interrupted
    // on N nodes can always be resumed on M nodes (including M == 1); see
    // Align2ClustChunking and Align2ClustCheckpoint.
    const size_t chunkSize = getAlign2ClustChunkSizeFromEnv();
    std::vector<Align2ClustChunking::Chunk> chunks = Align2ClustChunking::computeChunks(endRange, chunkSize);
    Align2ClustCheckpoint checkpoint(par, par.db3, dbSize, endRange, mode, chunkSize);

    bool anyChunkAlreadyDone = false;
    for (const Align2ClustChunking::Chunk &chunk : chunks) {
        if (checkpoint.isChunkDone(chunk.index)) {
            anyChunkAlreadyDone = true;
            break;
        }
    }
    // Every rank evaluates mpiNumProc/anyChunkAlreadyDone identically (shared
    // filesystem, same CLI parameters), so this decision is consistent across ranks.
    const bool useDistributedGeneration = (mpiNumProc > 1) || anyChunkAlreadyDone;

    if (useDistributedGeneration && par.includeAlignFiles) {
        Debug(Debug::ERROR) << "align2clust: --include-align-files is not supported together with "
                            << "MPI-distributed candidate generation (numProc > 1) or a resumed "
                            << "checkpointed run. Re-run on a single rank with the checkpoint "
                            << "directory (" << checkpoint.getChunkDir() << ") removed, or without "
                            << "--include-align-files.\n";
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
        return EXIT_FAILURE;
    }

    Align2ClustReducer reducer(reducerMode, dbSize, reorderCapacityChosen);
    if (mpiRank == 0) {
        reducer.start();
    }

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

    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
    size_t db_maxseqlen = (cluSeqDbr != nullptr)
        ? std::max(seqDbr->getMaxSeqLen(), cluSeqDbr->getMaxSeqLen())
        : seqDbr->getMaxSeqLen();

    // Candidate generation for [rangeStart, rangeEnd), pushing every produced result to
    // `sink`. Parameterized so it can run once over the full [0, endRange) range against
    // the live `reducer` (single-process/fresh-run fast path, byte-identical to the
    // pre-MPI behavior) or once per not-yet-done checkpoint chunk against a fresh
    // distributed-collect-only reducer whose buffered results are written to a chunk
    // file (see useDistributedGeneration below).
    auto generateCandidates = [&](Align2ClustReducer &sink, size_t rangeStart, size_t rangeEnd) {
        Debug::Progress progress(rangeEnd - rangeStart);
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
        // Staged member alignments; flushed only if the allpass-gate fully passes.
        std::string pendingMemberAln;
        std::vector<char> alnLineBuffer;
        if (includeAlignFiles) {
            alnLineBuffer.resize(1024 + 32768 * 4);
        }

#pragma omp for schedule(dynamic, 1) nowait
        for (size_t i = rangeStart; i < rangeEnd; i++) {
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
            if (sink.isAssigned(representativeId)) {
                sink.pushResult(std::move(clusterResult));
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
                    if (!sink.isAssigned(targetId)) {
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
                if (sink.isAssigned(representativeId)) {
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
                if (sink.isAssigned(targetId)) {
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
                if (sink.isAssigned(targetId)) continue;

                if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
                    if (sink.isAssigned(targetId)) continue;
                    if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                        // check all the member from filtering file
                        const size_t cluId = cluDbr->getId(targetKey);
                        char *cluData = cluDbr->getData(cluId, threadIdx);
                        const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                        size_t numClu = Util::countLines(cluData, cluDataSize);
                        bool allpass = true;
                        char buffer[1024];
                        if (includeAlignFiles) {
                            pendingMemberAln.clear();
                        }
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
                                short elementDiagonal = diagonal;
                                
                                // 1. ungapped alignment
                                element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
                                    allpass = false;
                                    break;
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
                                        elementUngappedAlignment.qStart, elementUngappedAlignment.qEnd, query.L,
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
                            ungappedAlignment.qStart, ungappedAlignment.qEnd, query.L,
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

                if (sink.isAssigned(targetId)) continue;

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
                        if (sink.isAssigned(targetId)) continue;
                        if (par.filterCluDBFile.empty()== false && par.filterSeqDBFile.empty()== false){
                            // check all the member from filtering file
                            const size_t cluId = cluDbr->getId(targetKey);
                            char *cluData = cluDbr->getData(cluId, threadIdx);
                            const size_t cluDataSize = cluDbr->getEntryLen(cluId);
                            size_t numClu = Util::countLines(cluData, cluDataSize);
                            bool allpass = true;
                            char buffer[1024];
                            if (includeAlignFiles) {
                                pendingMemberAln.clear();
                            }
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
                                    short elementDiagonal = 0;
                                    
                                    // 1. ungapped alignment
                                    element.mapSequence(elementId, elementKey, elementSequence, elementLength);
                                    if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
                                        allpass = false;
                                        break;
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
                                            elementUngappedAlignment.qStart, elementUngappedAlignment.qEnd, query.L,
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
            sink.pushResult(std::move(clusterResult));
        }
    }
    };  // end generateCandidates lambda

    if (useDistributedGeneration == false) {
        // Fast path: numProc == 1 and no pre-existing checkpoint chunks for this
        // fingerprint. Behavior is unchanged from before Phase 3: one pass over the
        // full range straight into the live reducer, no chunk files written or read.
        generateCandidates(reducer, 0, endRange);
    } else {
        checkpoint.ensureChunkDirExists();
        for (const Align2ClustChunking::Chunk &chunk : chunks) {
            if (Align2ClustChunking::isChunkOwnedByRank(chunk.index, mpiRank, mpiNumProc) == false) {
                continue;
            }
            if (checkpoint.isChunkDone(chunk.index)) {
                continue;
            }
            // A worker rank has no live view of cross-rank assignment state, so it
            // generates into a fresh distributed-collect-only reducer (no consumer
            // thread, isAssigned() always false -- see the Align2ClustReducer class
            // comment) and simply buffers whatever candidates it produces. The
            // reducer that eventually runs on rank 0 re-validates every candidate
            // against the true assignment state at consumption time, so this
            // superset of candidates is always safe.
            Align2ClustReducer collector(reducerMode, dbSize, 1, /*distributedCollectOnly=*/true);
            collector.start();
            generateCandidates(collector, chunk.start, chunk.end);
            collector.finish();
            checkpoint.writeChunk(chunk.index, collector.getCollectedResults());
        }

#ifdef HAVE_MPI
        if (mpiNumProc > 1) {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif

        if (mpiRank == 0) {
            // Only rank 0 re-reads the chunks (in global chunk order, which matches
            // ascending sequenceIdx order since chunks are contiguous ranges) and
            // feeds them into the single live reducer that performs the ordered
            // set-cover/greedy reduction.
            for (const Align2ClustChunking::Chunk &chunk : chunks) {
                std::vector<ClusterResult> results = checkpoint.readChunk(chunk.index);
                for (ClusterResult &result : results) {
                    reducer.pushResult(std::move(result));
                }
            }
        }
    }

    if (mpiRank != 0) {
        // Non-master ranks only ever generate/checkpoint candidates; the global
        // reduction and final DB write happen exclusively on rank 0 (see the
        // Align2ClustReducer class comment for why this step cannot itself be
        // decentralized). resultWriter/alnWriter are nullptr here (see align2clust()).
        delete[] fastMatrix.matrix;
        delete[] fastMatrix.matrixData;
        delete subMat;
        if (prefRepSizePair != nullptr) {
            delete[] prefRepSizePair;
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
        return EXIT_SUCCESS;
    }

    reducer.finish();

    reducer.finalizeSingletons();

    std::pair<DBKeyType, DBKeyType> *assignment = new std::pair<DBKeyType, DBKeyType>[dbSize];
    
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            const DBLocalId representativeId = reducer.getAssignment(i);
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
    
    writeClustering(resultWriter, assignment, dbSize);

    if (useDistributedGeneration) {
        // Safe to remove only now that the final clustering result has been fully
        // written; an interrupted run must keep its checkpoint directory intact so it
        // can be resumed (with any node/rank count) rather than recomputing candidates.
        checkpoint.cleanup();
    }

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
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    
    Timer timer;
    timer.reset();
    
    DBReader<DBKeyType> alnDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    alnDbr.open(DBReader<DBKeyType>::LINEAR_ACCCESS);
    int dbtype =  Parameters::DBTYPE_CLUSTER_RES;

    // Optional alignment-result output; path derived from the cluster DB (db3 + "_aln").
    // Alignment output needs a backtrace (-a) and score+cov+seqid so member records carry a CIGAR.
    // Checked on every rank (identical parameters everywhere) so a bad flag combination
    // fails the same way on every rank rather than only on rank 0.
    if (par.includeAlignFiles) {
        const unsigned int effectiveSwMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
        if (par.addBacktrace == false || effectiveSwMode != Matcher::SCORE_COV_SEQID) {
            Debug(Debug::ERROR) << "Writing alignment files requires backtrace and score+cov+seqid alignment.\n"
                                << "Please re-run with '-a 1' and '--alignment-mode "
                                << Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID << "'.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // Only rank 0 opens the output DBs: on a shared filesystem every other rank would
    // otherwise try to create/write the same output paths concurrently. doAlign2clust
    // never dereferences resultWriter/alnWriter on non-master ranks (see its
    // "mpiRank != 0" early return).
    DBWriter *resultWriter = nullptr;
    DBWriter *alnWriter = nullptr;
    if (MMseqsMPI::isMaster()) {
        resultWriter = new DBWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, dbtype);
        resultWriter->open();

        if (par.includeAlignFiles) {
            std::string alnDb = par.db3 + "_aln";
            std::string alnDbIndex = alnDb + ".index";
            alnWriter = new DBWriter(alnDb.c_str(), alnDbIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
            alnWriter->open();
        }
    }

    int status = doAlign2clust(par, resultWriter, alnDbr, alnWriter);

    Debug(Debug::INFO) << "Time for run Align2Clust: " << timer.lap() << " sec\n";

    if (resultWriter != nullptr) {
        resultWriter->close();
        delete resultWriter;
    }
    if (alnWriter != nullptr) {
        alnWriter->close();
        delete alnWriter;
    }
    alnDbr.close();

    return status;
}
