#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "IndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

int mergeresultsbyset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, true, 0);

    DBReader<unsigned int> setReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    setReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

//    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
//    resultReader.open(DBReader<unsigned int>::NOSORT);


    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    unsigned int databaseType = IndexReader::USER_SELECT;
    int dbtype = FileUtil::parseDbType(par.db2.c_str());
    if(Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_INDEX_DB)||
       Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_AMINO_ACIDS)||
       Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES)||
       Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)){
        databaseType = IndexReader::ALIGNMENTS;
    }
    IndexReader resultReader(par.db2, par.threads,
                              databaseType,
                             (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    int dbType = resultReader.sequenceReader->getDbtype();
    dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, dbType);
    dbw.open();
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);
        char dbKey[255];
#pragma omp for schedule(static)
        for (size_t i = 0; i < setReader.getSize(); ++i) {
            char *data = setReader.getData(i, thread_idx);
            // go through the results in the cluster and add them to one entry
            while (*data != '\0'){
                Util::parseKey(data, dbKey);
                unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
                size_t id = resultReader.sequenceReader->getId(key);
                if (id == UINT_MAX) {
                    Debug(Debug::ERROR) << "Invalid key " << key << " in entry " << i << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                buffer.append(resultReader.sequenceReader->getData(id, thread_idx));
                data = Util::skipLine(data);
            }
            dbw.writeData(buffer.c_str(), buffer.length(), setReader.getDbKey(i), thread_idx);
            buffer.clear();
        }
    }
    dbw.close();
    setReader.close();

    return EXIT_SUCCESS;
}
