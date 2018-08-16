#include <string>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2repseq(const Parameters &par, DBReader<unsigned int> &resultReader,
                   const std::string &outDb, const size_t dbFrom, const size_t dbSize) {

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str());
    qDbr.open(DBReader<unsigned int>::NOSORT);

    std::string outIndex(outDb);
    outIndex.append(".index");

    DBWriter resultWriter(outDb.c_str(), outIndex.c_str(), par.threads);
    resultWriter.open();

    Debug(Debug::INFO) << "Start computing representative sequences.\n";
#pragma omp parallel for schedule(dynamic, 100)
    for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
        Debug::printProgress(id);
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        char *results = resultReader.getData(id);
        if (*results == '\0') {
            Debug(Debug::WARNING) << "Empty result for entry " << id << "!\n";
            continue;
        }

        char dbKey[255];
        Util::parseKey(results, dbKey);
        const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
        const size_t edgeId = qDbr.getId(key);

        std::string result(qDbr.getData(edgeId));
        unsigned int queryKey = resultReader.getDbKey(id);
        resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
    }

    resultWriter.close();
    qDbr.close();

    Debug(Debug::INFO) << "\nDone.\n";
    return EXIT_SUCCESS;
}

int result2repseq(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int status;

    size_t resultSize = resultReader.getSize();
    status = result2repseq(par, resultReader, par.db3, 0, resultSize);

    resultReader.close();

    return status;
}
