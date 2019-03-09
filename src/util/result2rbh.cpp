#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2rbh(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbw.open();
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        const char *entry[255];

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            Debug::printProgress(id);

            unsigned int AdbID = resultReader.getDbKey(id);
            char *results = resultReader.getData(id, thread_idx);
            char *startEntryAtoB = results;
            char *endEntryAtoB = results;
            size_t entryAtoBlength = 0;

            int bestAtoBbitScore = 0; // initialize - no match of current A to any B
            unsigned int bestB = 0;

            while (*results != '\0') {
                Util::getWordsOfLine(results, entry, 255);
                unsigned int BdbID = Util::fast_atoi<int>(entry[0]);
                int currAlnScore = Util::fast_atoi<int>(entry[1]);

                results = Util::skipLine(results);

                // initialize with first line: this line is A-->B best hit:
                if (bestAtoBbitScore == 0) {
                    bestAtoBbitScore = currAlnScore;
                    bestB = BdbID;
                    endEntryAtoB = results;
                }
                // after the first line the bitscore can only decrease:
                else if (bestAtoBbitScore < currAlnScore) {
                    Debug(Debug::ERROR) << "The merged results are assumed to be sorted by decreasing bitscore.\n";
                    EXIT(EXIT_FAILURE);
                }
                // go over B-->A best hits:
                else if (currAlnScore == bestAtoBbitScore) {
                    if (BdbID == bestB) {
                        entryAtoBlength = (endEntryAtoB - startEntryAtoB);
                        break;
                        // no need to check anymore...
                    }
                } else {
                    break;
                    // worse bitscore - no need to check anymore...
                }
            }

            dbw.writeData(startEntryAtoB, entryAtoBlength, AdbID, thread_idx);
        }
    }
    dbw.close(true);
    resultReader.close();

    return EXIT_SUCCESS;
}
