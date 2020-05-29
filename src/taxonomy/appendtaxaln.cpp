#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int appendtaxaln(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // open tax assignments per sequence
    DBReader<unsigned int> taxReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    taxReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> alnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_TAXONOMICAL_RESULT);
    writer.open();


    Debug::Progress progress(taxReader.getSize());

    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        // per thread variables
        const char *entry[2048];
        std::string lineToWrite;
        lineToWrite.reserve(4096);

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < taxReader.getSize(); ++i) {
            progress.updateProgress();
            unsigned int seqKey = taxReader.getDbKey(i);
            char *results = taxReader.getData(i, thread_idx);

            lineToWrite = "";
            // each seqToTax is a single line so the loop will happen once
            while (*results != '\0') {
                Util::getWordsOfLine(results, entry, 255);
                unsigned int taxid = Util::fast_atoi<unsigned int>(entry[0]);
                lineToWrite += std::to_string(taxid);
                lineToWrite += "\t";

                char *seqToTargetAlnData = alnReader.getDataByDBKey(seqKey, thread_idx);
                // we care about only the first line in the alignment entry (toher do not refer to q vs t)
                while (*seqToTargetAlnData != '\n' && *seqToTargetAlnData != '\0') {
                    lineToWrite.push_back(*seqToTargetAlnData);
                    seqToTargetAlnData++; 
                }
                lineToWrite += "\n";
                results = Util::skipLine(results);
            }
            writer.writeData(lineToWrite.c_str(), lineToWrite.size(), seqKey, thread_idx);
            lineToWrite.clear();

        }
    };
    Debug(Debug::INFO) << "\n";

    writer.close();
    taxReader.close();
    alnReader.close();

    return EXIT_SUCCESS;
}
