#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "QueryMatcher.h"
#include "Parameters.h"

#include "omptl/omptl_algorithm"

#include <mutex>

#ifdef OPENMP
#include <omp.h>
#endif

// ((TargetKey,eVal),resultLine)
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

void writeSwappedResults(Parameters &par, std::vector<AlignmentResultEntry> *resMap,
                         unsigned int procNumber = 0, unsigned int nbOfProc = 1) {
    std::pair<std::string, std::string> outputDB =
            nbOfProc > 1 ? Util::createTmpFileNames(par.db4, par.db4Index, procNumber)
                         : std::make_pair(par.db4, par.db4Index);

    DBWriter resultWriter(outputDB.first.c_str(), outputDB.second.c_str(),
                          static_cast<unsigned int>(par.threads));
    resultWriter.open();

    size_t size = resMap->size();
#pragma omp parallel
    {
        int thread_num = 0, num_threads = 1;
        size_t start, end, orgStart, orgEnd;

#ifdef OPENMP
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
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

                    std::string result;
                    for (size_t j = 0; j < curRes.size(); j++) {
                        result.append(curRes[j].result);
                    }

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

void swapAlnResults(Parameters &par, std::vector<AlignmentResultEntry> *resMap,
                    unsigned int targetKeyMin, unsigned int targetKeyMax) {
    // evalue correction for the swapping
    double *kmnByLen = NULL;
    double lambda, logK, lambdaLog2, logKLog2;
    {
        DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
        qdbr.open(DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str());
        tdbr.open(DBReader<unsigned int>::NOSORT);

        SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix(subMat.getMatrixName(),
                                                                                      Matcher::GAP_OPEN,
                                                                                      Matcher::GAP_EXTEND, false);

        int seqLen = static_cast<int>(par.maxSeqLen);
        kmnByLen = new double[seqLen];
        for (int len = 0; len < seqLen; len++) {
            kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta,
                                                        qdbr.getAminoAcidDBSize(), qdbr.getSize());
        }

        lambda = stats.lambda;
        logK = log(stats.K);
        lambdaLog2 = lambda / log(2.0);
        logKLog2 = logK / log(2.0);

        qdbr.close();
        tdbr.close();
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    //std::cout<<sizeof(Matcher::result_t)<<std::endl; -> 96bytes

    int thread_num = 0, num_threads = 1;
    size_t start, end;
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

            //std::cout << 100.0*((float)i-start)/((float)end-start)<<std::endl;

            char *wordCnt[255];
            size_t cols = Util::getWordsOfLine(data, wordCnt, 254);
            if (cols >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                std::vector<Matcher::result_t> alnRes = Matcher::readAlignmentResults(data);

                for (size_t j = 0; j < alnRes.size(); j++) {
                    Matcher::result_t &res = alnRes[j];

                    unsigned int targetKey = res.dbKey;
                    if (targetKey >= targetKeyMin && targetKey < targetKeyMax) {
                        double rawScore = BlastScoreUtils::bitScoreToRawScore(res.score, lambdaLog2, logKLog2);
                        res.eval = BlastScoreUtils::computeEvalue(rawScore, kmnByLen[res.dbLen], lambda);
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
                        resMap->push_back(AlignmentResultEntry(targetKey, res.eval, result));
                        lock.unlock();
                    }
                }
            } else // prefilter case
            {
                char *curData = data;

                while (*curData != '\0') {
                    hit_t hit = parsePrefilterHit(curData);

                    unsigned int targetKey = hit.seqId;
                    if (targetKey >= targetKeyMin && targetKey < targetKeyMax) {
                        hit.seqId = queryKey;

                        float eval = exp(-hit.prefScore);

                        std::string result = prefilterHitToString(hit);
                        lock.lock();
                        resMap->push_back(AlignmentResultEntry(targetKey, eval, result));
                        lock.unlock();
                    }

                    curData = Util::skipLine(curData);
                }
            }
        }
    }

    delete[] kmnByLen;
    resultReader.close();
}


int doSwapSort(Parameters &par, unsigned int procNumber = 0, unsigned int nbOfProc = 1) {
    std::vector<AlignmentResultEntry> resMap;
    unsigned int targetKeyMin = 0, targetKeyMax = -1;

    if (nbOfProc > 1) {
        DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str());
        targetReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
        unsigned int lastKey = targetReader.getLastKey();
        targetReader.close();

        targetKeyMin = procNumber * (lastKey + 1) / nbOfProc;
        targetKeyMax = (procNumber + 1) * (lastKey + 1) / nbOfProc;
    }

    swapAlnResults(par, &resMap, targetKeyMin, targetKeyMax);

    omptl::sort(resMap.begin(), resMap.end(), compareKey()); // sort by target id

    writeSwappedResults(par, &resMap, procNumber, nbOfProc);

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

#ifdef HAVE_MPI
    int status = doSwapSort(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doSwapSort(par);
#endif

    EXIT(status);
}
