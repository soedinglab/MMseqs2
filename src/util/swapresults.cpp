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

#ifdef OPENMP
#include <omp.h>
#endif

class IndexReader {
public:
    DBReader<unsigned int> *reader;
    IndexReader(const std::string &data) {
        std::string index = data;
        index.append(".index");

        bool isIndex = false;
        std::string indexData = PrefilteringIndexReader::searchForIndex(data);
        if (indexData != "") {
            Debug(Debug::INFO) << "Use index: " << indexData << "\n";
            std::string indexIndex = indexData;
            indexIndex.append(".index");
            indexReader = new DBReader<unsigned int>(indexData.c_str(), indexIndex.c_str());
            indexReader->open(DBReader<unsigned int>::NOSORT);
            if ((isIndex = PrefilteringIndexReader::checkIfIndexFile(indexReader)) == true) {
                reader = PrefilteringIndexReader::openNewReader(indexReader, false);
            } else {
                Debug(Debug::WARNING) << "Outdated index version. Please recompute it with 'createindex'!\n";
                indexReader->close();
                delete indexReader;
                indexReader = NULL;
            }
        }

        if (isIndex == false) {
            reader = new DBReader<unsigned int>(data.c_str(), index.c_str(), DBReader<unsigned int>::USE_INDEX);
            reader->open(DBReader<unsigned int>::NOSORT);
        }
    }

    ~IndexReader() {
        reader->close();
        delete reader;

        if (indexReader != NULL) {
            indexReader->close();
            delete indexReader;
        }
    }
private:
    DBReader<unsigned int> *indexReader = NULL;
};

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

    size_t aaResSize = 0;
    unsigned int maxTargetId = 0;
    char *targetElementExists = NULL;
    if (isGeneralMode) {
        Debug(Debug::INFO) << "Result database: " << parResultDbStr << "\n";
        DBReader<unsigned int> resultReader(parResultDb, parResultDbIndex);
        resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
        //search for the maxTargetId (value of first column) in parallel
#pragma omp parallel
        {
            char key[255];
#pragma omp for schedule(dynamic, 100) reduction(max:maxTargetId)
            for (size_t i = 0; i < resultReader.getSize(); ++i) {
                Debug::printProgress(i);
                char *data = resultReader.getData(i);
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
        Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
        IndexReader query(par.db1);
        aaResSize = query.reader->getAminoAcidDBSize();

        Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
        IndexReader target(par.db2);
        maxTargetId = target.reader->getLastKey();

        targetElementExists = new char[maxTargetId + 1];
        memset(targetElementExists, 0, sizeof(char) * (maxTargetId + 1));
#pragma omp parallel for
        for (size_t i = 0; i < target.reader->getSize(); ++i) {
            unsigned int key = target.reader->getDbKey(i);
            targetElementExists[key] = 1;
        }
    }
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    EvalueComputation evaluer(aaResSize, &subMat, par.gapOpen, par.gapExtend, true);

    Debug(Debug::INFO) << "Result database: " << parResultDbStr << "\n";
    DBReader<unsigned int> resultDbr(parResultDb, parResultDbIndex);
    resultDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const size_t resultSize = resultDbr.getSize();
    Debug(Debug::INFO) << "Computing offsets.\n";
    size_t *targetElementSize = new size_t[maxTargetId + 2]; // extra element for offset + 1 index id
    memset(targetElementSize, 0, sizeof(size_t) * (maxTargetId + 2));
#pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = 0; i < resultSize; ++i) {
        Debug::printProgress(i);
        const unsigned int resultId = resultDbr.getDbKey(i);
        char queryKeyStr[1024];
        char *tmpBuff = Itoa::u32toa_sse2((uint32_t) resultId, queryKeyStr);
        *(tmpBuff) = '\0';
        size_t queryKeyLen = strlen(queryKeyStr);
        char *data = resultDbr.getData(i);
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

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
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
        char *tmpData = new char[bytesToWrite];
        Util::checkAllocation(tmpData, "Could not allocate tmpData memory in doswap");
        Debug(Debug::INFO) << "\nReading results.\n";
#pragma omp parallel for schedule(dynamic, 10)
        for (size_t i = 0; i < resultSize; ++i) {
            Debug::printProgress(i);
            char *data = resultDbr.getData(i);
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
        //revert offsets
        for (unsigned int i = maxTargetId + 1; i > 0; i--) {
            targetElementSize[i] = targetElementSize[i - 1];
        }
        targetElementSize[0] = 0;

        Debug(Debug::INFO) << "\nOutput database: " << parOutDbStr << "\n";
        char *entry[255];
        bool isAlignmentResult = false;
        bool hasBacktrace = false;
        for (size_t i = 0; i < resultDbr.getSize(); i++){
            if (resultDbr.getSeqLens(i) <= 1){
                continue;
            }
            const size_t columns = Util::getWordsOfLine(resultDbr.getData(i), entry, 255);
            isAlignmentResult = columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT;
            hasBacktrace = columns >= Matcher::ALN_RES_WITH_BT_COL_CNT;
            break;
        }

        std::string splitDbw = parOutDbStr + "_" + SSTR(split);
        std::pair<std::string, std::string> splitNamePair = std::make_pair(splitDbw, splitDbw + ".index");
        splitFileNames.push_back(splitNamePair);

        DBWriter resultWriter(splitNamePair.first.c_str(), splitNamePair.second.c_str(), par.threads);
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
            char buffer[1024+32768];
            std::string ss;
            ss.reserve(100000);

#pragma omp for schedule(dynamic, 100)
            for (size_t i = prevDbKeyToWrite; i <= dbKeyToWrite; ++i) {
                Debug::printProgress(i);

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
                        double rawScore = evaluer.computeRawScoreFromBitScore(res.score);
                        res.eval = evaluer.computeEvalue(rawScore, res.dbLen);
                        if (res.eval > par.evalThr) {
                            evalBreak = true;
                            goto outer;
                        }
                        unsigned int qstart = res.qStartPos;
                        unsigned int qend = res.qEndPos;
                        unsigned int qLen = res.qLen;
                        res.qStartPos = res.dbStartPos;
                        res.qEndPos = res.dbEndPos;
                        res.qLen = res.dbLen;
                        res.dbStartPos = qstart;
                        res.dbEndPos = qend;
                        res.dbLen = qLen;
                        if (hasBacktrace) {
                            for (size_t j = 0; j < res.backtrace.size(); j++) {
                                if (res.backtrace.at(j) == 'I') {
                                    res.backtrace.at(j) = 'D';
                                } else if (res.backtrace.at(j) == 'D') {
                                    res.backtrace.at(j) = 'I';
                                }
                            }
                        }
                        curRes.emplace_back(res);
                    } else {
                        hit_t hit = QueryMatcher::parsePrefilterHit(data);
                        hit.diagonal = static_cast<unsigned short>(static_cast<short>(hit.diagonal) * -1);
                        curRes.emplace_back(hit.seqId, 0, hit.pScore, 0, 0, -hit.pScore, hit.diagonal, 0, 0, 0, 0, 0, 0, "");
                    }
                    outer:
                    char *nextLine = Util::skipLine(data);
                    size_t lineLen = nextLine - data;
                    dataSize -= lineLen;
                    data = nextLine;
                }

                if (curRes.empty() == false) {
                    if (curRes.size() > 1) {
                        std::sort(curRes.begin(), curRes.end(), Matcher::compareHits);
                    }

                    for (size_t j = 0; j < curRes.size(); j++) {
                        const Matcher::result_t &res = curRes[j];
                        if (isAlignmentResult) {
                            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                            ss.append(buffer, len);
                        } else {
                            hit_t hit;
                            hit.seqId = res.dbKey;
                            hit.pScore = res.qcov;
                            hit.diagonal = res.alnLength;
                            size_t len = QueryMatcher::prefilterHitToBuffer(buffer, hit);
                            ss.append(buffer, len);
                        }
                    }

                    resultWriter.writeData(ss.c_str(), ss.size(), i, thread_idx);
                    ss = "";

                    curRes.clear();
                } else if (evalBreak == true || targetElementExists[i] == 1) {
                    resultWriter.writeData(&empty, 0, i, thread_idx);
                }
            }
        };
        Debug(Debug::INFO) << "\n";
        resultWriter.close();

        prevDbKeyToWrite = dbKeyToWrite + 1;
        prevBytesToWrite += bytesToWrite;
        delete[] tmpData;
    }
    DBWriter::mergeResults(parOutDbStr, parOutDbIndexStr, splitFileNames);

    resultDbr.close();
    if (targetElementExists != NULL) {
        delete[] targetElementExists;
    }
    delete[] targetElementSize;
    return EXIT_SUCCESS;
}

int swapdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    return doswap(par, true);
}

int swapresults(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);
    return doswap(par, false);
}
