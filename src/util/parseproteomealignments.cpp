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

#ifdef OPENMP
#include <omp.h>
#endif

static std::vector<unsigned int> indexSources(DBReader<DBKeyType> &dbr, std::vector<size_t> &sourceToIdx, std::vector<DBKeyType> &sourceIds) {
    std::vector<unsigned int> numEntries(dbr.getSourceSize(), 0);
#pragma omp parallel
    {
        std::vector<unsigned int> localNumEntries(numEntries.size(), 0);
#pragma omp for schedule(static)
        for (size_t i = 0; i < dbr.getSize(); i++) {
            localNumEntries[static_cast<size_t>(dbr.getLookupFileNumber(dbr.getLookupIdByKey(dbr.getDbKey(i))))]++;
        }
#pragma omp critical
        {
            for (size_t source = 0; source < numEntries.size(); source++) {
                numEntries[source] += localNumEntries[source];
            }
        }
    }
    std::vector<unsigned int> sourceNumEntries;
    sourceToIdx.assign(numEntries.size(), SIZE_MAX);
    for (size_t source = 0; source < numEntries.size(); source++) {
        if (numEntries[source] > 0) {
            sourceToIdx[source] = sourceIds.size();
            sourceIds.push_back(static_cast<DBKeyType>(source));
            sourceNumEntries.push_back(numEntries[source]);
        }
    }
    return sourceNumEntries;
}

int parseproteomealignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // Open DBs
    DBReader<DBKeyType> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_LOOKUP|DBReader<DBKeyType>::USE_SOURCE);
    qdbr.open(DBReader<DBKeyType>::NOSORT);

    DBReader<DBKeyType> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_LOOKUP|DBReader<DBKeyType>::USE_SOURCE);
    tdbr.open(DBReader<DBKeyType>::NOSORT);

    std::vector<size_t> qSourceToIdx;
    std::vector<DBKeyType> qSourceIds;
    const std::vector<unsigned int> qSourceNumEntries = indexSources(qdbr, qSourceToIdx, qSourceIds);
    std::vector<size_t> tSourceToIdx;
    std::vector<DBKeyType> tSourceIds;
    indexSources(tdbr, tSourceToIdx, tSourceIds);
    const size_t qdbSourceSize = qSourceIds.size();
    const size_t tdbSourceSize = tSourceIds.size();

    DBReader<DBKeyType> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_INDEX);
    alndbr.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    int proteomeDBType = DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, proteomeDBType);
    resultWriter.open();

    unsigned int *scores = new(std::nothrow) unsigned int[tdbSourceSize * qdbSourceSize];
    Util::checkAllocation(scores, "Can not allocate scores memory in parseproteomealignments");
    std::fill(scores, scores + (tdbSourceSize * qdbSourceSize), 0);
    unsigned int **scoreLookupTable = new unsigned int*[qdbSourceSize]; 
    for (size_t i = 0; i < qdbSourceSize; ++i) { scoreLookupTable[i] = scores + i * tdbSourceSize; }
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(alndbr.getSize()) / static_cast<double>(flushSize)));

    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(alndbr.getSize() - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char buffer[1024 + 32768*4];
            std::vector<unsigned int> localMatchResults(tdbSourceSize, 0);

#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();

                const DBKeyType queryDbKey = alndbr.getDbKey(id);
                // const unsigned int qId = qdbr.getId(queryDbKey);
                size_t qLookupId = qdbr.getLookupIdByKey(queryDbKey);
                const size_t querySourceIdx = qSourceToIdx[qdbr.getLookupFileNumber(qLookupId)];
                char *data = alndbr.getData(id, thread_idx);
                // localMatchResults.clear();
                std::fill(localMatchResults.begin(), localMatchResults.end(), 0);
                while (*data != '\0') {
                    Util::parseKey(data, buffer);
                    const DBKeyType targetDbKey = Util::fast_atoi<DBKeyType>(buffer);
                    // const unsigned int tId = tdbr.getId(targetDbKey);  
                    size_t tLookupId = tdbr.getLookupIdByKey(targetDbKey);
                    const size_t targetSourceIdx = tSourceToIdx[tdbr.getLookupFileNumber(tLookupId)];
                    if (localMatchResults[targetSourceIdx] == 0) {
                        localMatchResults[targetSourceIdx] = 1;
                    }
                    data = Util::skipLine(data);
                }
                for (size_t t = 0; t < tdbSourceSize; ++t) {
                    if (localMatchResults[t] != 0) {
#pragma omp atomic
                        scoreLookupTable[querySourceIdx][t] += localMatchResults[t];
                    }
                }
            }
        }
    }

    std::string alnResultsOutString;
    alnResultsOutString.reserve(1024*1024);
    char buffer[1024];
    for (size_t q = 0; q < qdbSourceSize; ++q) {
        for (size_t t = 0; t < tdbSourceSize; ++t) {
            char* tmpBuffer = Itoa::u64toa_sse2(static_cast<uint64_t>(tSourceIds[t]), buffer);
            *(tmpBuffer - 1) = '\t';
            float score = static_cast<float>(scoreLookupTable[q][t]) / static_cast<float>(qSourceNumEntries[q]);
            tmpBuffer = Util::fastSeqIdToBuffer(score, tmpBuffer);
            *(tmpBuffer - 1) = '\n';
            alnResultsOutString.append(buffer, tmpBuffer - buffer);
        }
        resultWriter.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), qSourceIds[q], 0);
        alnResultsOutString.clear();
    }
    delete[] scores;
    delete[] scoreLookupTable;
    
    alndbr.close();
    qdbr.remapData();
    qdbr.close();
    tdbr.remapData();
    tdbr.close();
    resultWriter.close();

    return EXIT_SUCCESS;
}
