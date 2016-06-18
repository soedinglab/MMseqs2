#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int count(int argn, const char **argv) {
    std::string usage;
    usage.append("Count number of lines in every entry.\n");
    usage.append("USAGE: <ffindexInDB> <ffindexOutDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.count, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads));
    writer.open();

    size_t entries = reader.getSize();
#pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < entries; ++i) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        char* data = reader.getData(i);
        size_t length = reader.getSeqLens(i) - 1;

        size_t lines = 0;
        for(size_t j = 0; j < length; ++j) {
            char c = data[j];
            if(c == '\n') {
                lines++;
            }
        }

        unsigned int id = reader.getDbKey(i);
        std::string result = SSTR(lines);
        result.append("\n");
        writer.write(result.c_str(), result.length(), SSTR(id).c_str(), thread_idx);
    }

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}
