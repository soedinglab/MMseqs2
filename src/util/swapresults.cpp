#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "Parameters.h"
#include "MemoryMapped.h"

#include "omptl/omptl_algorithm"

#include <mutex>

#ifdef OPENMP
#include <omp.h>
#endif

struct AlignmentResultEntry {
    unsigned int key;
    double evalue;
    std::string result;

    AlignmentResultEntry(unsigned int key, double evalue, const std::string &result) :
            key(key), evalue(evalue), result(result) {}
};

struct compareKey {
    bool operator()(const AlignmentResultEntry &lhs,
                    const AlignmentResultEntry &rhs) const {
        return lhs.key < rhs.key;
    }
};

struct compareEval {
    bool operator()(const AlignmentResultEntry &lhs,
                    const AlignmentResultEntry &rhs) const {
        return lhs.evalue < rhs.evalue;
    }
};

std::pair<size_t,size_t> longestLine(const char* name) {
    MemoryMapped file(name, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char * buf = (char *) file.getData();

    size_t lines = 0;
    size_t maxLen = 0;

#pragma omp parallel reduction(+:lines) reduction(max:maxLen)
    {
        int thread_num = 0,num_threads = 1;
        size_t startChunk,endChunk;

#ifdef OPENMP
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
#endif
        startChunk =  thread_num * file.size() / num_threads;
        endChunk = (thread_num + 1) * file.size() / num_threads;

        size_t i = startChunk;
        while(i && i<file.size() && buf[i] != '\n') // wedge on line start
            i++;
        if(!i)
            i++; // jump the '\n'
        size_t charCount = 0;
        while(i<file.size() && i < endChunk)
        {

            if (buf[i] == '\n')
            {
                lines++;
                maxLen = std::max(maxLen,charCount);
                charCount = 0;
            } else if (buf[i] != '\0') {
                charCount++;
            }
            i++;
        }

        while(i<file.size()  && buf[i] != '\n')
        {
            if (buf[i] != '\0')
                charCount++;
            i++;
        }
        //lines++;
        maxLen = std::max(maxLen,charCount);

    }

    file.close();
    return std::make_pair(lines,maxLen);
}

void writeSwappedResults(const std::string &outputDB, const std::string &outputDBIndex,
                         std::vector<AlignmentResultEntry> *resMap) {
    int num_threads = 1;
#ifdef OPENMP
    num_threads = omp_get_num_threads();
#endif

    DBWriter resultWriter(outputDB.c_str(), outputDBIndex.c_str(), num_threads);
    resultWriter.open();

    size_t size = resMap->size();
#pragma omp parallel
    {
        int thread_num = 0;
        size_t start, end, orgStart, orgEnd;

#ifdef OPENMP
        thread_num = omp_get_thread_num();
#endif
        orgStart = thread_num * size / num_threads;
        orgEnd = (thread_num + 1) * size / num_threads;

        // Wedge the start and end pos to complete targetKey chunks
        start = orgStart;
        while (orgStart && start < size && resMap->at(orgStart).key == resMap->at(start).key) {
            start++;
        }

        end = orgEnd;
        while (end < size && resMap->at(orgEnd).key == resMap->at(end).key) {
            end++;
        }

        // When db is too small, only the first thread will deal with the swap
        if (size < static_cast<size_t>(num_threads)) {
            if (thread_num == 0) {
                end = size;
            } else {
                start = end;
            }
        }

        if (end - start) {
            unsigned int lastKey = resMap->at(start).key;

            std::vector<AlignmentResultEntry> curRes;
            for (size_t i = start; i < end; ++i) {
                // If we enter a new chunk, flush the results to file
                if (lastKey != resMap->at(i).key) {
                    omptl::sort(curRes.begin(), curRes.end(), compareEval());

                    std::ostringstream ss;
                    for (size_t j = 0; j < curRes.size(); j++) {
                        ss << curRes[j].result;
                    }

                    std::string result = ss.str();
                    resultWriter.writeData(result.c_str(), result.size(), lastKey, thread_num);
                    curRes.clear();
                }

                curRes.push_back(resMap->at(i));
                lastKey = resMap->at(i).key;
            }

            omptl::sort(curRes.begin(), curRes.end(), compareEval());

            std::ostringstream ss;
            for (size_t j = 0; j < curRes.size(); j++) {
                ss << curRes[j].result;
            }

            std::string result = ss.str();
            resultWriter.writeData(result.c_str(), result.size(), lastKey, thread_num);
        }
    }
    resultWriter.close();
}

void swapBt(std::string &bt) {
    for (size_t i = 0; i < bt.size(); i++) {
        if (bt.at(i) == 'I') {
            bt.at(i) = 'D';
        } else if (bt.at(i) == 'D') {
            bt.at(i) = 'I';
        }
    }
}

void swapResults(DBReader<unsigned int> &resultReader,
                    std::vector<AlignmentResultEntry> *resMap,
                    unsigned int targetKeyMin,
                    unsigned int targetKeyMax,
                    size_t aaResSize,
                    const char *scoringMatrix) {
    SubstitutionMatrix subMat(scoringMatrix, 2.0, 0.0);
    EvalueComputation evaluer(aaResSize, &subMat, Matcher::GAP_OPEN, Matcher::GAP_EXTEND, true);

    int thread_num = 0, num_threads = 1;
    size_t start, end;
    size_t count = 0;
    size_t size = resultReader.getSize();
    std::mutex lock;

#pragma omp parallel private(thread_num,num_threads,start,end)
    {
#ifdef OPENMP
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
#endif
        start = thread_num * size / num_threads;
        end = (thread_num + 1) * size / num_threads;

        for (size_t i = start; i < end; ++i) {
            unsigned int queryKey = resultReader.getDbKey(i);
            char *data = resultReader.getData(i);

            char *wordCnt[255];
            size_t cols = Util::getWordsOfLine(data, wordCnt, 254);
            if (cols >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                std::vector<Matcher::result_t> alnRes = Matcher::readAlignmentResults(data);

                for (size_t j = 0; j < alnRes.size(); j++) {
                    Matcher::result_t &res = alnRes[j];

                    unsigned int targetKey = res.dbKey;
                    if (targetKey >= targetKeyMin && targetKey < targetKeyMax) {
                        double rawScore = evaluer.computeBitScore(res.score);
//                        double rawScore = evaluer.computeRawScoreFromBitScore(res.score);
                        res.eval = evaluer.computeEvalue(rawScore, res.dbLen);
                        unsigned int qstart = res.qStartPos;
                        unsigned int qend = res.qEndPos;
                        unsigned int qLen = res.qLen;
                        res.qStartPos = res.dbStartPos;
                        res.qEndPos = res.dbEndPos;
                        res.qLen = res.dbLen;
                        res.dbStartPos = qstart;
                        res.dbEndPos = qend;
                        res.dbLen = qLen;

                        res.dbKey = queryKey;

                        bool addBacktrace = cols >= Matcher::ALN_RES_WITH_BT_COL_CNT;
                        if (addBacktrace) {
                            swapBt(res.backtrace);

                        }

                        std::string result = Matcher::resultToString(res, addBacktrace);

                        lock.lock();
                        (*resMap)[count].key = targetKey;
                        (*resMap)[count].evalue = res.eval;
                        (*resMap)[count].result += result;
                        count++;
                        lock.unlock();
                    }
                }
            } else // prefilter case
            {
                char* curData = data;
                while (*curData != '\0') {
                    hit_t hit = parsePrefilterHit(curData);
                
                    unsigned int targetKey = hit.seqId;
                    if(targetKey >= targetKeyMin && targetKey < targetKeyMax)
                    {
                        hit.seqId = queryKey;

                        float eval = exp(-hit.prefScore);

                        std::string result = prefilterHitToString(hit);
                        lock.lock();
                        (*resMap)[count].key = targetKey;
                        (*resMap)[count].evalue = eval;
                        (*resMap)[count].result += result;
                        count++;
                        lock.unlock();
                    }

                    curData = Util::skipLine(curData);
                }
            }
        }
    }

    resultReader.close();
}

int doSwapSort(Parameters &par, unsigned int procNumber = 0, unsigned int nbOfProc = 1) {
    AlignmentResultEntry nullEntry(0, 0, std::string());

    unsigned int targetKeyMin = 0;
    unsigned int targetKeyMax = static_cast<unsigned int>(-1);
    std::pair<size_t, size_t> searchFileLines = longestLine(par.db3.c_str());
    unsigned long numberOfEntries = searchFileLines.first;

    if (nbOfProc > 1) {
        DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str());
        targetReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
        unsigned int lastKey = targetReader.getLastKey();
        targetReader.close();

        targetKeyMin = procNumber * (lastKey + 1) / nbOfProc;
        targetKeyMax = (procNumber + 1) * (lastKey + 1) / nbOfProc;
    }

    std::vector<AlignmentResultEntry> resMap(numberOfEntries, nullEntry);

    // the results length should be <= the maximum original length + 10
    size_t lineLength = searchFileLines.second + 10;
    // the length difference stems from target id field
    for (size_t i = 0; i < numberOfEntries; i++) {
        resMap[i].result.reserve(lineLength);
    }


    size_t aaResSize = 0;
    {
        DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
        qdbr.open(DBReader<unsigned int>::NOSORT);
        aaResSize = qdbr.getAminoAcidDBSize();
        qdbr.close();
    }

    const char* scoringMatrix = par.scoringMatrixFile.c_str();

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    swapResults(resultReader, &resMap, targetKeyMin, targetKeyMax, aaResSize, scoringMatrix);
    resultReader.close();

    omptl::sort(resMap.begin(), resMap.end(), compareKey()); // sort by target id

   if (nbOfProc > 1)  {
       std::pair<std::string, std::string> tmpDb = Util::createTmpFileNames(par.db4, par.db4Index, procNumber);
       writeSwappedResults(tmpDb.first, tmpDb.second, &resMap);
   } else {
       writeSwappedResults(par.db4, par.db4Index, &resMap);

   }

    // In case of MPI parallelization, merge the partial results
    if (nbOfProc > 1) {
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (procNumber == 0) {
            std::vector<std::pair<std::string, std::string>> partialResFiles;
            for (unsigned int proc = 0; proc < nbOfProc; ++proc) {
                partialResFiles.push_back(Util::createTmpFileNames(par.db4, par.db4Index, proc));
            }

            DBWriter::mergeResults(par.db4, par.db4Index, partialResFiles);
        }
    }

    return EXIT_SUCCESS;
}

int swapresults(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#ifdef HAVE_MPI
    int status = doSwapSort(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doSwapSort(par);
#endif

    EXIT(status);
}
