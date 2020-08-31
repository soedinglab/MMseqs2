#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2repseq(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> seqReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, seqReader.getDbtype());
    resultWriter.open();
    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        char dbKey[255];
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < resultReader.getSize(); ++id) {
            progress.updateProgress();

            char *results = resultReader.getData(id, thread_idx);
            if (*results == '\0') {
                continue;
            }

            Util::parseKey(results, dbKey);
            const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            const size_t edgeId = seqReader.getId(key);
            resultWriter.writeData(seqReader.getData(edgeId, thread_idx), seqReader.getEntryLen(edgeId) - 1, resultReader.getDbKey(id), thread_idx);
        }
    }
    resultWriter.close(true);
    resultReader.close();
    seqReader.close();
    DBReader<unsigned int>::softlinkDb(par.db1, par.db3, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}
