#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "Parameters.h"
#include "AlignmentSymmetry.h"
#include "PrefilteringIndexReader.h"
#include "IndexReader.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif

int doswap(Parameters& par, bool isGeneralMode) {
    const char * parResultDb;
    const char * parResultDbIndex;
    const char * parOutDb;
    const char * parOutDbIndex;

    if (isGeneralMode) {
        parResultDb = par.db1.c_str();
        parResultDbIndex = par.db1Index.c_str();
        parOutDb = par.db2.c_str();
        parOutDbIndex = par.db2Index.c_str();

    } else {
        parResultDb = par.db3.c_str();
        parResultDbIndex = par.db3Index.c_str();
        parOutDb = par.db4.c_str();
        parOutDbIndex = par.db4Index.c_str();
    }
    std::string parResultDbStr(parResultDb);
    std::string parResultDbIndexStr(parResultDbIndex);

    std::string parOutDbStr(parOutDb);
    std::string parOutDbIndexStr(parOutDbIndex);

    BaseMatrix *subMat = NULL;
    EvalueComputation *evaluer = NULL;
    size_t aaResSize = 0;
    unsigned int maxTargetId = 0;
    char *targetElementExists = NULL;
    if (isGeneralMode) {
        DBReader<unsigned int> resultReader(parResultDb, parResultDbIndex, par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        resultReader.open(DBReader<unsigned int>::SORT_BY_OFFSET);
        //search for the maxTargetId (value of first column) in parallel
        Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char key[255];
#pragma omp for schedule(dynamic, 100) reduction(max:maxTargetId)
            for (size_t i = 0; i < resultReader.getSize(); ++i) {
                progress.updateProgress();
                char *data = resultReader.getData(i, thread_idx);
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int dbKey = std::strtoul(key, NULL, 10);
                    maxTargetId = std::max(maxTargetId, dbKey);
                    data = Util::skipLine(data);
                }
            }
        };
        resultReader.close();
    } else {
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        IndexReader query(par.db1, par.threads, IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        aaResSize = query.sequenceReader->getAminoAcidDBSize();

        IndexReader target(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        maxTargetId = target.sequenceReader->getLastKey();

        targetElementExists = new char[maxTargetId + 1];
        memset(targetElementExists, 0, sizeof(char) * (maxTargetId + 1));
#pragma omp parallel for
        for (size_t i = 0; i < target.sequenceReader->getSize(); ++i) {
            unsigned int key = target.sequenceReader->getDbKey(i);
            targetElementExists[key] = 1;
        }
        int gapOpen, gapExtend;
        if (Parameters::isEqualDbtype(target.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
            gapOpen = par.gapOpen.values.nucleotide();
            gapExtend =  par.gapExtend.values.nucleotide();
        } else {
            // keep score bias at 0.0 (improved ROC)
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
            gapOpen = par.gapOpen.values.aminoacid();
            gapExtend = par.gapExtend.values.aminoacid();
        }
        evaluer = new EvalueComputation(aaResSize, subMat, gapOpen, gapExtend);
    }

    DBReader<unsigned int> resultDbr(parResultDb, parResultDbIndex, par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultDbr.open(DBReader<unsigned int>::SORT_BY_OFFSET);

    const size_t resultSize = resultDbr.getSize();
    Debug(Debug::INFO) << "Computing offsets.\n";
    size_t *targetElementSize = new size_t[maxTargetId + 2]; // extra element for offset + 1 index id
    memset(targetElementSize, 0, sizeof(size_t) * (maxTargetId + 2));
    {
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
                size_t queryKeyLen = strlen(queryKeyStr);
                char *data = resultDbr.getData(i, thread_idx);
                char dbKeyBuffer[255 + 1];
                while (*data != '\0') {
                    Util::parseKey(data, dbKeyBuffer);
                    size_t targetKeyLen = strlen(dbKeyBuffer);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    char *nextLine = Util::skipLine(data);
                    size_t lineLen = nextLine - data;
                    lineLen -= targetKeyLen;
                    lineLen += queryKeyLen;
                    __sync_fetch_and_add(&(targetElementSize[dbKey]), lineLen);
                    data = nextLine;
                }
            }
        }
    }

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    size_t bytesForTargetElements = sizeof(size_t) * (maxTargetId + 2);
    memoryLimit = (memoryLimit > bytesForTargetElements) ? (memoryLimit - bytesForTargetElements) : 0;

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

    const char empty = '\0';

    unsigned int prevDbKeyToWrite = 0;
    size_t prevBytesToWrite = 0;
    for (size_t split = 0; split < splits.size(); split++) {
        unsigned int dbKeyToWrite = splits[split].first;
        size_t bytesToWrite = splits[split].second;
        char *tmpData = new(std::nothrow) char[bytesToWrite];
        Util::checkAllocation(tmpData, "Cannot allocate tmpData memory");
        Debug(Debug::INFO) << "\nReading results.\n";
        Debug::Progress progress(resultSize);
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
                size_t queryKeyLen = strlen(queryKeyStr);
                char dbKeyBuffer[255 + 1];
                while (*data != '\0') {
                    Util::parseKey(data, dbKeyBuffer);
                    size_t targetKeyLen = strlen(dbKeyBuffer);
                    char *nextLine = Util::skipLine(data);
                    size_t oldLineLen = nextLine - data;
                    size_t newLineLen = oldLineLen;
                    newLineLen -= targetKeyLen;
                    newLineLen += queryKeyLen;
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    // update offset but do not copy memory
                    size_t offset = __sync_fetch_and_add(&(targetElementSize[dbKey]), newLineLen) - prevBytesToWrite;
                    if(dbKey >= prevDbKeyToWrite && dbKey <=  dbKeyToWrite){
                        memcpy(&tmpData[offset], queryKeyStr, queryKeyLen);
                        memcpy(&tmpData[offset + queryKeyLen], data + targetKeyLen, oldLineLen - targetKeyLen);
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

        Debug(Debug::INFO) << "\nOutput database: " << parOutDbStr << "\n";
        bool isAlignmentResult = false;
        bool hasBacktrace = false;
        const char *entry[255];
        for (size_t i = 0; i < resultDbr.getSize(); i++){
            char *data = resultDbr.getData(i, 0);
            if (*data == '\0'){
                continue;
            }
            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            isAlignmentResult = columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT;
            hasBacktrace = columns >= Matcher::ALN_RES_WITH_BT_COL_CNT;
            break;
        }

        std::string splitDbw = parOutDbStr + "_" + SSTR(split);
        std::pair<std::string, std::string> splitNamePair = (splits.size() > 1) ? std::make_pair(splitDbw, splitDbw + ".index") :
                                                            std::make_pair(parOutDb, parOutDbIndex) ;
        splitFileNames.push_back(splitNamePair);
        Debug::Progress progress2(dbKeyToWrite - prevDbKeyToWrite  + 1);

        DBWriter resultWriter(splitNamePair.first.c_str(), splitNamePair.second.c_str(), par.threads, par.compressed, resultDbr.getDbtype());
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
            curRes.reserve(300);

            char buffer[1024 + 32768*4];
            std::string ss;
            ss.reserve(100000);

#pragma omp for schedule(dynamic, 100)
            for (size_t i = prevDbKeyToWrite; i <= dbKeyToWrite; ++i) {
                progress2.updateProgress();

                char *data = &tmpData[targetElementSize[i] - prevBytesToWrite];
                size_t dataSize = targetElementSize[i + 1] - targetElementSize[i];

                if (isGeneralMode) {
                    if (dataSize > 0) {
                        resultWriter.writeData(data, dataSize, i, thread_idx);
                    }
                    continue;
                }

                bool evalBreak = false;
                while (dataSize > 0) {
                    if (isAlignmentResult) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
                        Matcher::result_t::swapResult(res, *evaluer, hasBacktrace);
                        if (res.eval > par.evalThr) {
                            evalBreak = true;
                        } else {
                            curRes.emplace_back(res);
			}
                    } else {
                        hit_t hit = QueryMatcher::parsePrefilterHit(data);
                        hit.diagonal = static_cast<unsigned short>(static_cast<short>(hit.diagonal) * -1);
                        curRes.emplace_back(hit.seqId, hit.prefScore, 0, 0, 0, -static_cast<float>(hit.prefScore), hit.diagonal, 0, 0, 0, 0, 0, 0, "");
                    }
                    char *nextLine = Util::skipLine(data);
                    size_t lineLen = nextLine - data;
                    dataSize -= lineLen;
                    data = nextLine;
                }

                if (curRes.empty() == false) {
                    if (curRes.size() > 1) {
                        SORT_SERIAL(curRes.begin(), curRes.end(), Matcher::compareHits);
                    }

                    for (size_t j = 0; j < curRes.size(); j++) {
                        const Matcher::result_t &res = curRes[j];
                        if (isAlignmentResult) {
                            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                            ss.append(buffer, len);
                        } else {
                            hit_t hit;
                            hit.seqId = res.dbKey;
                            hit.prefScore = res.score;
                            hit.diagonal = res.alnLength;
                            size_t len = QueryMatcher::prefilterHitToBuffer(buffer, hit);
                            ss.append(buffer, len);
                        }
                    }

                    resultWriter.writeData(ss.c_str(), ss.size(), i, thread_idx);
                    ss.clear();

                    curRes.clear();
                } else if (evalBreak == true || targetElementExists[i] == 1) {
                    resultWriter.writeData(&empty, 0, i, thread_idx);
                }
            }
        };
        Debug(Debug::INFO) << "\n";
        if(splits.size() > 1){
            resultWriter.close(true);
        }else{
            resultWriter.close();
        }

        prevDbKeyToWrite = dbKeyToWrite + 1;
        prevBytesToWrite += bytesToWrite;
        delete[] tmpData;
    }
    if(splits.size() > 1){
        DBWriter::mergeResults(parOutDbStr, parOutDbIndexStr, splitFileNames);
    }

    if (evaluer != NULL) {
        delete evaluer;
    }

    if (subMat != NULL) {
        delete subMat;
    }

    resultDbr.close();
    if (targetElementExists != NULL) {
        delete[] targetElementExists;
    }
    delete[] targetElementSize;
    return EXIT_SUCCESS;
}

int swapdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    return doswap(par, true);
}

int swapresults(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    return doswap(par, false);
}
