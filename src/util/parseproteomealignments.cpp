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

int parseproteomealignments(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // Open DBs
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_LOOKUP|DBReader<unsigned int>::USE_SOURCE);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    size_t qdbSourceSize = qdbr.getSourceSize();
   
    DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_LOOKUP|DBReader<unsigned int>::USE_SOURCE);
    tdbr.open(DBReader<unsigned int>::NOSORT);
    size_t tdbSourceSize = tdbr.getSourceSize();

    DBReader<unsigned int> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int proteomeDBType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, proteomeDBType);
    resultWriter.open();

    unsigned int *scores = new unsigned int[tdbSourceSize * qdbSourceSize]; 
    std::fill(scores, scores + (tdbSourceSize * qdbSourceSize), 0);
    Util::checkAllocation(scores, "Can not allocate scores memory in parseproteomealignments");
    unsigned int **scoreLookupTable = new unsigned int*[qdbSourceSize]; 
    for (size_t i = 0; i < qdbSourceSize; ++i) { scoreLookupTable[i] = scores + i * tdbSourceSize; }
    unsigned int * qSourceIdtoNumEntries = new unsigned int[qdbSourceSize];
    Util::checkAllocation(qSourceIdtoNumEntries, "Can not allocate qSourceIdtoNumEntries memory in parseproteomealignments");
    memset(qSourceIdtoNumEntries, 0, qdbSourceSize * sizeof(unsigned int));

    for (size_t i = 0; i < qdbr.getSize(); i++) {
        unsigned int dbKey = qdbr.getDbKey(i);
        size_t lookupId = qdbr.getLookupIdByKey(dbKey);
        const unsigned int proteomeSourceId = qdbr.getLookupFileNumber(lookupId);
        qSourceIdtoNumEntries[proteomeSourceId]++;
    }

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

                const unsigned int queryDbKey = alndbr.getDbKey(id);
                // const unsigned int qId = qdbr.getId(queryDbKey);
                size_t qLookupId = qdbr.getLookupIdByKey(queryDbKey);
                const unsigned int querySourceId = qdbr.getLookupFileNumber(qLookupId);
                char *data = alndbr.getData(id, thread_idx);
                // localMatchResults.clear();
                std::fill(localMatchResults.begin(), localMatchResults.end(), 0);
                while (*data != '\0') {
                    Util::parseKey(data, buffer);
                    const unsigned int targetDbKey = (unsigned int) strtoul(buffer, NULL, 10);
                    // const unsigned int tId = tdbr.getId(targetDbKey);  
                    size_t tLookupId = tdbr.getLookupIdByKey(targetDbKey);
                    const unsigned int targetSourceId = tdbr.getLookupFileNumber(tLookupId);
                    if (localMatchResults[targetSourceId] == 0) {
                        localMatchResults[targetSourceId] = 1;
                    }
                    data = Util::skipLine(data);
                }
                for (size_t t = 0; t < tdbSourceSize; ++t) {
                    if (localMatchResults[t] != 0) {
#pragma omp atomic
                        scoreLookupTable[querySourceId][t] += localMatchResults[t];
                    }
                }
            }
        }

        std::string alnResultsOutString;
        alnResultsOutString.reserve(1024*1024);
        char buffer[1024];

        for (size_t q = 0; q < qdbSourceSize; ++q) {
            for (size_t t = 0; t < tdbSourceSize; ++t) {
                char* tmpBuffer = Itoa::i32toa_sse2(t, buffer);
                *(tmpBuffer - 1) = '\t';
                float score = static_cast<float>(scoreLookupTable[q][t]) / static_cast<float>(qSourceIdtoNumEntries[q]);
                tmpBuffer = Util::fastSeqIdToBuffer(score, tmpBuffer);
                *(tmpBuffer - 1) = '\n';
                alnResultsOutString.append(buffer, tmpBuffer - buffer);
            }
            resultWriter.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), q, 0);
            alnResultsOutString.clear();
        }
    }
    delete[] scores;
    delete[] scoreLookupTable;
    delete[] qSourceIdtoNumEntries;
    
    alndbr.close();
    qdbr.remapData();
    qdbr.close();
    tdbr.remapData();
    tdbr.close();
    resultWriter.close();

    return EXIT_SUCCESS;
}


