//
// Created by mad on 10/21/15.
//
#include <string>
#include <vector>
#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"

int rescorediagonal(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    DBReader<unsigned int> *qdbr = NULL;
    DBReader<unsigned int> *tdbr = NULL;
    Debug(Debug::WARNING) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    qdbr->readMmapedDataInMemory();
    Debug(Debug::WARNING) << "Target  file: " << par.db2 << "\n";
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), (par.db2 + ".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }
    Debug(Debug::WARNING) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_res(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    dbr_res.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::WARNING) << "Result database: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBWriter::BINARY_MODE);
    resultWriter.open();
    const size_t flushSize = 1000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbr_res.getSize()) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbr_res.getSize() - (i * flushSize), flushSize);

#pragma omp parallel for schedule(dynamic, 10)
        for (size_t id = start; id < (start + bucketSize); id++) {
            Debug::printProgress(id);
            char buffer[100];
            std::string prefResultsOutString;
            prefResultsOutString.reserve(1000000);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char *data = dbr_res.getData(id);
            unsigned int queryId = qdbr->getId(dbr_res.getDbKey(id));
            char *querySeq = qdbr->getData(queryId);
            unsigned int queryLen = qdbr->getSeqLens(queryId) - 2; // - 2 because of /0/n
            std::vector<hit_t> results = parsePrefilterHits(data);
            for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                unsigned int targetId = tdbr->getId(results[entryIdx].seqId);
                unsigned int targetLen = tdbr->getSeqLens(targetId) - 2; // - 2 because of /0/n
                short diagonal = results[entryIdx].diagonal;
                unsigned short distanceToDiagonal = abs(diagonal);
                unsigned int diagonalLen = 0;
                unsigned int distance = 0;
                if (diagonal >= 0 && distanceToDiagonal < queryLen) {
                    diagonalLen = std::min(targetLen, queryLen - distanceToDiagonal);
                    distance = DistanceCalculator::computeHammingDistance(querySeq + distanceToDiagonal,
                                                                          tdbr->getData(targetId), diagonalLen);
                } else if (diagonal < 0 && distanceToDiagonal < targetLen) {
                    diagonalLen = std::min(targetLen - distanceToDiagonal, queryLen);
                    distance = DistanceCalculator::computeHammingDistance(querySeq,
                                                                          tdbr->getData(targetId) + distanceToDiagonal,
                                                                          diagonalLen);
                }

                float seqId = (static_cast<float>(diagonalLen) - static_cast<float>(distance)) /
                              static_cast<float>(diagonalLen);
                float targetCov = static_cast<float>(diagonalLen) / static_cast<float>(targetLen);
                float queryCov = static_cast<float>(diagonalLen) / static_cast<float>(queryLen);
                if (targetCov >= (par.targetCovThr - std::numeric_limits<float>::epsilon()) &&
                    seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon())) {
                    int len = snprintf(buffer, 100, "%s\t%d\t%.2f\t%.2f\t%.2f\n", SSTR(results[entryIdx].seqId).c_str(),
                                       diagonal, seqId, queryCov, targetCov);
                    prefResultsOutString.append(buffer, len);
                }
            }
            // write prefiltering results string to ffindex database
            const size_t prefResultsLength = prefResultsOutString.length();
            char *prefResultsOutData = (char *) prefResultsOutString.c_str();
            resultWriter.writeData(prefResultsOutData, prefResultsLength, SSTR(qdbr->getDbKey(queryId)).c_str(),
                                   thread_idx);
        }
        dbr_res.remapData();
    }
    Debug(Debug::WARNING) << "Done." << "\n";
    dbr_res.close();
    resultWriter.close();
    qdbr->close();
    delete qdbr;
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }
    return 0;
}


