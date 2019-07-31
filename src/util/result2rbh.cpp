#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2rbh(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbw.open();
    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        const char *entry[255];
        std::string rbhBsString;
        rbhBsString.reserve(100000);

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();

            unsigned int AdbID = resultReader.getDbKey(id);
            char *results = resultReader.getData(id, thread_idx);
            int bestAtoBbitScore = 0; // initialize - no match of current A to any B

            char *startRbhB = results;
            char *endRbhB = results;
            size_t RbhBlength  = 0;

            // go over A->B direction (getting bestAtoBbitScore):
            while (*results != '\0') {
                Util::getWordsOfLine(results, entry, 255);
                int currAlnScore = Util::fast_atoi<int>(entry[1]);
                if (bestAtoBbitScore == 0) {
                    // first iteration is an A->B line - update bestAtoBbitScore
                    bestAtoBbitScore = currAlnScore;
                    results = Util::skipLine(results);
                } else {
                    // this is a B->A line - the bitscore can only decrease:
                    if (bestAtoBbitScore < currAlnScore) {
                        Debug(Debug::ERROR) << "The merged results are assumed to be sorted by decreasing bitscore.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    if (currAlnScore < bestAtoBbitScore) {
                        // worse bitscore - no need to check anymore...
                        break;
                    }
                    // current B->A has bestAtoBbitScore - we retain this line for writing!
                    startRbhB = results;
                    results = Util::skipLine(results); 
                    endRbhB = results;
                    RbhBlength = (endRbhB - startRbhB);
                    rbhBsString.append(startRbhB, RbhBlength);
                }
            }

            // write entry for A (write null byte);
            dbw.writeData(rbhBsString.c_str(), rbhBsString.length(), AdbID, thread_idx);
            rbhBsString.clear();
        }
    }
    dbw.close(true);
    resultReader.close();

    return EXIT_SUCCESS;
}
