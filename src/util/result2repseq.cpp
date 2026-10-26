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

    DBReader<DBKeyType> seqReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    seqReader.open(DBReader<DBKeyType>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        seqReader.readMmapedDataInMemory();
    }

    DBReader<DBKeyType> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    resultReader.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

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
            const DBKeyType key = Util::fast_atoi<DBKeyType>(dbKey);
            const size_t edgeId = seqReader.getId(key);
            if (edgeId == DB_ENTRY_NOT_FOUND) {
                Debug(Debug::ERROR) << "Sequence " << key << " does not exist in sequence database.\n";
                EXIT(EXIT_FAILURE);
            }
            resultWriter.writeData(seqReader.getData(edgeId, thread_idx), seqReader.getEntryLen(edgeId) - 1, resultReader.getDbKey(id), thread_idx);
        }
    }
    resultWriter.close(true);
    resultReader.close();
    seqReader.close();
    DBReader<DBKeyType>::softlinkDb(par.db1, par.db3, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}
