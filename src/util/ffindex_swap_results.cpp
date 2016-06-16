#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Alignment.h"
#include "Util.h"

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

typedef std::map<unsigned int, std::string *> SwapIt;

SwapIt readAllKeysIntoMap(const std::string &datafile) {
    SwapIt result;

    char dbKey[255 + 1];
    std::ifstream resultFile(datafile);
    std::string line;
    while (std::getline(resultFile, line)) {
        size_t i;
        for (i = 0; i < line.size() - 1; i++) {
            if (line[i] != '\0') {
                break;
            }
        }

        char *data = const_cast<char *>(line.c_str()) + i;

        Util::parseKey(data, dbKey);
        unsigned int key = static_cast<unsigned int>(strtoul(dbKey, NULL, 10));

        if (dbKey[0] == '\0')
            continue;

        SwapIt::iterator it = result.find(key);
        if (it == result.end()) {
            result[key] = new std::string();
        }
    }
    resultFile.close();

    return result;
}

void processSplit(DBReader<unsigned int> &dbr,
                  FILE *dataFile, const std::map<unsigned int, std::string *> &map,
                  size_t startIndex, size_t domainSize) {
    for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
        std::ostringstream ss;
        char dbKey[255 + 1];
        int c1;
        while ((c1 = fgetc(dataFile)) != EOF && c1 != (int) '\0') {
            ss << c1;
        }
        ss << '\0';
        std::string result = ss.str();

        std::string outerKey = SSTR(dbr.getDbKey(i));
        char *data = (char *) result.c_str();
        while (*data != '\0') {
            // extract key from results (ids must be always at the first position)
            Util::parseKey(data, dbKey);
            unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            std::string *entry = NULL;
            SwapIt::const_iterator it = map.find(key);
            if (it != map.end()) {
                entry = it->second;
            } else {
                entry = new std::string();
            }
            // write data to map
            entry->append(outerKey);
            // next db key
            char *endPosOfId = data + Util::skipNoneWhitespace(data);
            data = Util::skipLine(data);
            entry->append(endPosOfId, data);
        }
    }
}

int doSwap(Parameters &par, DBReader<unsigned int> &rdbr,
           const std::pair<std::string, std::string> &resultdb,
           const size_t startIndex, const size_t domainSize) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str());
    tdbr.open(DBReader<unsigned int>::NOSORT);

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix(subMat.getMatrixName(),
                                                                                  Matcher::GAP_OPEN,
                                                                                  Matcher::GAP_EXTEND);

    int seqLen = static_cast<int>(par.maxSeqLen);
    double *kmnByLen = new double[seqLen];
    for (int len = 0; len < seqLen; len++) {
        kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta,
                                                    qdbr.getAminoAcidDBSize(), qdbr.getSize());
    }
    const double logK = log(stats.K);
    const double lambdaLog2 = stats.lambda / log(2.0);
    const double logKLog2 = logK / log(2.0);

    // read all keys
    const SwapIt swapMap = readAllKeysIntoMap(rdbr.getDataFileName());

    // create and open db write
    DBWriter splitWriter(resultdb.first.c_str(), resultdb.second.c_str(),
                         static_cast<unsigned int>(par.threads));
    splitWriter.open();

    FILE *dataFile = fopen(rdbr.getDataFileName(), "r");
    processSplit(rdbr, dataFile, swapMap, startIndex, domainSize);
    fclose(dataFile);

    Debug(Debug::INFO) << "Start to swap results. Write to " << resultdb.first << ".\n";

    // Write sorted results and delete memory
    char *wordCnt[255];
#pragma omp parallel for schedule(dynamic, 100)
    for (SwapIt::const_iterator it = swapMap.begin(); it != swapMap.end(); ++it) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned int id = it->first;
        std::string result = *it->second;
        if (result.size() > 1) {
            size_t cols = Util::getWordsOfLine((char *) result.c_str(), wordCnt, 254);
            if (Matcher::ALN_RES_WITH_OUT_BT_COL_CNT >= cols) {
                std::vector<Matcher::result_t> alnRes = Matcher::readAlignmentResults((char *) result.c_str());
                for (size_t i = 0; i < alnRes.size(); i++) {
                    Matcher::result_t &res = alnRes[i];
                    double rawScore = BlastScoreUtils::bitScoreToRawScore(res.score, lambdaLog2, logKLog2);
                    res.eval = BlastScoreUtils::computeEvalue(rawScore, kmnByLen[res.dbLen], stats.lambda);
                    unsigned int qstart = res.qStartPos;
                    unsigned int qend = res.qEndPos;
                    unsigned int qLen = res.qLen;
                    res.qStartPos = res.dbStartPos;
                    res.qEndPos = res.dbEndPos;
                    res.qLen = res.dbLen;
                    res.dbStartPos = qstart;
                    res.dbEndPos = qend;
                    res.dbLen = qLen;
                }

                std::stable_sort(alnRes.begin(), alnRes.end(), Matcher::compareHits);

                std::ostringstream swResultsSs;
                for (size_t i = 0; i < alnRes.size(); i++) {
                    bool addBacktrace = Matcher::ALN_RES_WITH_BT_COL_CNT == cols;
                    swResultsSs << Matcher::resultToString(alnRes[i], addBacktrace);
                }
                result = swResultsSs.str();
            }
        }

        splitWriter.write(result.c_str(), result.size(), SSTR(id).c_str(), thread_idx);

        // remove just the value (string *)
        delete it->second;
    }

    Debug(Debug::INFO) << "Done.\n";
    splitWriter.close();

    qdbr.close();
    tdbr.close();
    delete[] kmnByLen;

    return EXIT_SUCCESS;
}

int doSwap(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomain(reader.getSize(), mpiRank, mpiNumProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db4, par.db4Index, mpiRank);

    int status = doSwap(par, reader, tmpOutput, dbFrom, dbSize);
    reader.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // master reduces results
    if (mpiRank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (unsigned int proc = 0; proc < mpiNumProc; ++proc) {
            std::pair<std::string, std::string> names = Util::createTmpFileNames(par.db4, par.db4Index, proc);
            splitFiles.push_back(std::make_pair(names.first, names.first + ".index"));
        }
        Alignment::mergeAndRemoveTmpDatabases(par.db4, par.db4 + ".index", splitFiles);
    }

    return status;
}

int doSwap(Parameters &par) {
    size_t resultSize;

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    resultSize = reader.getSize();

    int status = doSwap(par, reader, std::make_pair(par.db4, par.db4Index), 0, resultSize);
    reader.close();

    return status;
}

int swapresults(int argc, const char *argv[]) {
    MMseqsMPI::init(argc, argv);

    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n");
    usage.append("USAGE: <queryDb> <targetDb> <ffindexDB> <fastaDB>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.swapresults, 4);

#ifdef HAVE_MPI
    int status = doSwap(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doSwap(par);
#endif

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return status;
}
