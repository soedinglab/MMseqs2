#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"
#include "Util.h"
#include "omptl/omptl_algorithm"

#include <list>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

int mergeorfcontigs(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, true);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str());
    resultReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> contigReader(par.db2.c_str(), par.db2Index.c_str());
    contigReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    dbw.open();
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(1024);
        char dbKey[255];
#pragma omp for schedule(static)
        for (size_t i = 0; i < contigReader.getSize(); ++i) {
            char *data = contigReader.getData(i);
            // go through the results in the cluster and add them to one entry
            while (*data != '\0'){
                Util::parseKey(data, dbKey);
                unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
                buffer.append(resultReader.getDataByDBKey(key));
                data = Util::skipLine(data);
            }
            dbw.writeData(buffer.c_str(), buffer.length(), contigReader.getDbKey(i), thread_idx);
        }
        buffer.clear();
    };
    dbw.close();
    contigReader.close();
    resultReader.close();

    return EXIT_SUCCESS;
}
