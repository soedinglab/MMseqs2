#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "CompressedA3M.h"

#ifdef OPENMP
#include <omp.h>
#endif

int convertca3m(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);


    DBReader<std::string> reader((par.db1 + "_ca3m.ffdata").c_str(), (par.db1 + "_ca3m.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<std::string>::NOSORT);

    DBReader<unsigned int> sequences((par.db1 + "_sequence.ffdata").c_str(), (par.db1 + "_sequence.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    sequences.open(DBReader<unsigned int>::SORT_BY_LINE);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CA3M_DB);
    writer.open();

    Debug::Progress progress(reader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::vector<Matcher::result_t> results;
        results.reserve(1000);

        char buffer[1024 + 32768*4];

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            results.clear();

            unsigned int key;
            CompressedA3M::extractMatcherResults(key, results, reader.getData(i, thread_idx), reader.getEntryLen(i), sequences, true);

            writer.writeStart(thread_idx);
            for (size_t j = 0; j < results.size(); j++) {
                size_t len = Matcher::resultToBuffer(buffer, results[j], true);
                writer.writeAdd(buffer, len, thread_idx);
            }

            writer.writeEnd(key, thread_idx);
        }
    }
    writer.close();
    sequences.close();
    reader.close();

    return EXIT_SUCCESS;
}


