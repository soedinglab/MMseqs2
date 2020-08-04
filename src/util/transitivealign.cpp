#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "Alignment.h"
#include "itoa.h"
#include "BacktraceTranslator.h"
#include "AlignmentSymmetry.h"
#include "DistanceCalculator.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif

int transitivealign(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> sequenceDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sequenceDbr.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        sequenceDbr.readMmapedDataInMemory();
    }
    const int querySeqType = sequenceDbr.getDbtype();

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    }

    DBReader<unsigned int> alnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);


    std::string tmpRes = par.db3+".tmp";
    std::string tmpResIndex = par.db3+".tmp.index";
    DBWriter resultWriter(tmpRes.c_str(), tmpResIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    EvalueComputation evaluer(sequenceDbr.getAminoAcidDBSize(), subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(alnReader.getSize()) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(alnReader.getSize() - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif

            Matcher matcher(querySeqType, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.gapOpen.aminoacids, par.gapExtend.aminoacids, par.zdrop);

//            Sequence query(par.maxSeqLen, targetSeqType, subMat, par.kmerSize, par.spacedKmer, par.compBiasCorrection);
//            Sequence target(par.maxSeqLen, targetSeqType, subMat, par.kmerSize, par.spacedKmer, par.compBiasCorrection);

            char * buffer = new char[1024 + 32768*4];
            BacktraceTranslator btTranslate;
            std::vector<Matcher::result_t> results;
            results.reserve(300);
            std::vector<Matcher::result_t> outputResults;
            outputResults.reserve(300);

#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();

                const unsigned int alnKey = alnReader.getDbKey(id);
                char *data = alnReader.getData(id, thread_idx);

                results.clear();
                Matcher::readAlignmentResults(results, data, false);
                resultWriter.writeStart(thread_idx);
                for (size_t entryIdx_i = 0; entryIdx_i < results.size(); entryIdx_i++) {
                    const unsigned int queryId = sequenceDbr.getId(results[entryIdx_i].dbKey);
                    const unsigned int queryKey = sequenceDbr.getDbKey(queryId);
                    // we need A->B->C to infer A->C
                    // in center start the oriontation is B->A
                    // so we need to swap the result A->B
                    Matcher::result_t swappedResult(results[entryIdx_i]);
                    Matcher::result_t::swapResult(swappedResult, evaluer, true);
                    char *querySeq = sequenceDbr.getData(queryId, thread_idx);

                    char * tmpBuff = Itoa::u32toa_sse2((uint32_t) queryKey, buffer);
                    *(tmpBuff-1) = '\t';
                    const unsigned int queryIdLen = tmpBuff - buffer;
                    if(queryKey == alnKey){
                        for (size_t aliId = 0; aliId < results.size(); aliId++) {
                            size_t len = Matcher::resultToBuffer(tmpBuff, results[aliId], true, true);
                            resultWriter.writeAdd(buffer, queryIdLen + len, thread_idx);
                        }
                        continue;
                    }

                    for (size_t entryIdx_j = 0; entryIdx_j < results.size(); entryIdx_j++) {
                        const unsigned int targetId = sequenceDbr.getId(results[entryIdx_j].dbKey);
                        char *targetSeq = sequenceDbr.getData(targetId, thread_idx);

                        if (Util::canBeCovered(par.covThr, par.covMode, swappedResult.qLen, results[entryIdx_j].dbLen) == false) {
                            continue;
                        }
                        const bool isIdentity = (queryId == targetId && par.includeIdentity) ? true : false;
                        Matcher::result_t result;
                        if(results[entryIdx_i].dbKey == results[entryIdx_j].dbKey ){
                            unsigned int score = DistanceCalculator::computeSubstitutionDistance(querySeq, targetSeq, results[entryIdx_i].dbLen, fastMatrix.matrix);
                            double bitScore = evaluer.computeBitScore(score);
                            double evalue   = evaluer.computeEvalue(score, results[entryIdx_j].dbLen);
                            result.dbKey = results[entryIdx_j].dbKey;
                            result.dbLen = results[entryIdx_j].dbLen;
                            result.score = bitScore;
                            result.qLen = results[entryIdx_j].dbLen;
                            result.dbEndPos = results[entryIdx_j].dbLen-1;
                            result.qEndPos  = results[entryIdx_j].dbLen-1;
                            result.dbStartPos = 0;
                            result.qStartPos  = 0;
                            result.eval    = evalue;
                            result.score    = bitScore;
                            result.seqId = 1.0f;
                            result.alnLength = results[entryIdx_j].dbLen;
                            result.backtrace = "";
                            result.backtrace.insert(0, result.alnLength, 'M');;
//                            result.backtrace.push_back('M');
                        }else{
                            btTranslate.translateResult(swappedResult, results[entryIdx_j], result);
                            Matcher::updateResultByRescoringBacktrace(querySeq, targetSeq, fastMatrix.matrix, evaluer, par.gapOpen.aminoacids, par.gapExtend.aminoacids, result);
                        }
                        // checkCriteria and Util::canBeCovered always work together
                        if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            outputResults.push_back(result);
                        }
                    }
                    SORT_SERIAL(outputResults.begin(), outputResults.end(), Matcher::compareHits);
                    for (size_t aliId = 0; aliId < outputResults.size(); aliId++) {
                        size_t len = Matcher::resultToBuffer(tmpBuff, outputResults[aliId], true, true);
                        resultWriter.writeAdd(buffer, queryIdLen + len, thread_idx);
                    }
                    outputResults.clear();
                }
                resultWriter.writeEnd(alnKey, thread_idx);
            }
            delete [] buffer;
        }
        alnReader.remapData();
    }
    resultWriter.close();
    alnReader.close();

    // logic to merge the results
    size_t maxTargetId = sequenceDbr.getLastKey();
    char * targetElementExists = new char[maxTargetId + 1];
    memset(targetElementExists, 0, sizeof(char) * (maxTargetId + 1));
#pragma omp parallel for
    for (size_t i = 0; i < sequenceDbr.getSize(); ++i) {
        unsigned int key = sequenceDbr.getDbKey(i);
        targetElementExists[key] = 1;
    }



    DBReader<unsigned int> resultDbr(tmpRes.c_str(), tmpResIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const size_t resultSize = resultDbr.getSize();
    Debug(Debug::INFO) << "Computing offsets.\n";
    size_t *targetElementSize = new size_t[maxTargetId + 2]; // extra element for offset + 1 index id
    memset(targetElementSize, 0, sizeof(size_t) * (maxTargetId + 2));
    Debug::Progress progress(resultSize);

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp  for schedule(dynamic, 100)
        for (size_t i = 0; i < resultSize; ++i) {
            progress.updateProgress();
            const unsigned int resultId = resultDbr.getDbKey(i);
            char queryKeyStr[1024];
            char *tmpBuff = Itoa::u32toa_sse2((uint32_t) resultId, queryKeyStr);
            *(tmpBuff) = '\0';
            char *data = resultDbr.getData(i, thread_idx);
            char dbKeyBuffer[255 + 1];
            while (*data != '\0') {
                Util::parseKey(data, dbKeyBuffer);
                size_t targetKeyLen = strlen(dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                char *nextLine = Util::skipLine(data);
                size_t lineLen = nextLine - data;
                lineLen -= (targetKeyLen + 1);
                __sync_fetch_and_add(&(targetElementSize[dbKey]), lineLen);
                data = nextLine;
            }
        }
    }

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    // compute splits
    std::vector<std::pair<unsigned int, size_t > > splits;
    std::vector<std::pair<std::string , std::string > > splitFileNames;
    size_t bytesToWrite = 0;
    for (size_t i = 0; i <= maxTargetId; i++) {
        bytesToWrite += targetElementSize[i];
        if (bytesToWrite > memoryLimit){
            splits.push_back(std::make_pair(i-1, bytesToWrite-targetElementSize[i]));
            bytesToWrite = targetElementSize[i];
        }
    }
    splits.push_back(std::make_pair(maxTargetId, bytesToWrite));
    AlignmentSymmetry::computeOffsetFromCounts(targetElementSize, maxTargetId + 1);

    std::string parOutDbStr(par.db3);
    std::string parOutDbIndexStr(par.db3Index);

    unsigned int prevDbKeyToWrite = 0;
    size_t prevBytesToWrite = 0;
    for (size_t split = 0; split < splits.size(); split++) {
        unsigned int dbKeyToWrite = splits[split].first;
        size_t bytesToWrite = splits[split].second;
        char *tmpData = new char[bytesToWrite];
        Util::checkAllocation(tmpData, "Could not allocate tmpData memory in doswap");
        Debug(Debug::INFO) << "\nReading results.\n";
#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
            for (size_t i = 0; i < resultSize; ++i) {
                progress.updateProgress();
                char *data = resultDbr.getData(i, thread_idx);
                unsigned int queryKey = resultDbr.getDbKey(i);
                char queryKeyStr[1024];
                char *tmpBuff = Itoa::u32toa_sse2((uint32_t) queryKey, queryKeyStr);
                *(tmpBuff) = '\0';
                char dbKeyBuffer[255 + 1];
                while (*data != '\0') {
                    Util::parseKey(data, dbKeyBuffer);
                    size_t targetKeyLen = strlen(dbKeyBuffer);
                    char *nextLine = Util::skipLine(data);
                    size_t oldLineLen = nextLine - data;
                    size_t newLineLen = oldLineLen;
                    newLineLen -= (targetKeyLen + 1);
                    //newLineLen += queryKeyLen;
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    // update offset but do not copy memory
                    size_t offset = __sync_fetch_and_add(&(targetElementSize[dbKey]), newLineLen) - prevBytesToWrite;
                    if (dbKey >= prevDbKeyToWrite && dbKey <= dbKeyToWrite) {
                        //memcpy(&tmpData[offset], queryKeyStr, queryKeyLen);
                        memcpy(&tmpData[offset], data + (targetKeyLen + 1), oldLineLen - (targetKeyLen + 1));
                    }
                    data = nextLine;
                }
            }
        }
        //revert offsets
        for (unsigned int i = maxTargetId + 1; i > 0; i--) {
            targetElementSize[i] = targetElementSize[i - 1];
        }
        targetElementSize[0] = 0;

        Debug(Debug::INFO) << "\nOutput database: " << par.db3 << "\n";
        std::string splitDbw = parOutDbStr + "_" + SSTR(split);
        std::pair<std::string, std::string> splitNamePair = (splits.size() > 1) ? std::make_pair(splitDbw,
                                                                                                 splitDbw + ".index") :
                                                            std::make_pair(par.db3, par.db3Index);
        splitFileNames.push_back(splitNamePair);

        DBWriter resultWriter(splitNamePair.first.c_str(), splitNamePair.second.c_str(), par.threads, par.compressed,
                              resultDbr.getDbtype());
        resultWriter.open();
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            // we are reusing this vector also for the prefiltering results
            // qcov is used for pScore because its the first float value
            // and alnLength for diagonal because its the first int value after
            std::vector<Matcher::result_t> curRes;
            std::string ss;
            ss.reserve(100000);

#pragma omp for schedule(dynamic, 100)
            for (size_t i = prevDbKeyToWrite; i <= dbKeyToWrite; ++i) {
                progress.updateProgress();

                char *data = &tmpData[targetElementSize[i] - prevBytesToWrite];
                size_t dataSize = targetElementSize[i + 1] - targetElementSize[i];

                if (dataSize > 0) {
                    resultWriter.writeData(data, dataSize, i, thread_idx);
                }

            }
        };
        Debug(Debug::INFO) << "\n";
        if (splits.size() > 1) {
            resultWriter.close(true);
        } else {
            resultWriter.close();
        }
        prevDbKeyToWrite = dbKeyToWrite + 1;
        prevBytesToWrite += bytesToWrite;
        delete[] tmpData;
    }

    DBReader<unsigned int>::removeDb(tmpRes);

    if(splits.size() > 1){
        DBWriter::mergeResults(parOutDbStr, parOutDbIndexStr, splitFileNames);
    }

    sequenceDbr.close();
    delete subMat;

    return EXIT_SUCCESS;
}
