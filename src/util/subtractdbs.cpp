//
// Created by mad on 2/6/16.
//

#include <climits>
#include <list>
#include <vector>
#include <sys/time.h>
#include <Matcher.h>
#include "DBReader.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Util.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

void dosubstractresult(std::string leftDb, std::string rightDb, std::string outDb,
                       size_t maxLineLength, double evalThreshold, int threads)
{
    Debug(Debug::INFO) << "Remove " << rightDb << " ids from "<< leftDb << "\n";
    DBReader<unsigned int> leftDbr(leftDb.c_str(), (leftDb + std::string(".index")).c_str());
    leftDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> rightDbr(rightDb.c_str(), (rightDb + std::string(".index")).c_str());
    rightDbr.open(DBReader<unsigned int>::NOSORT);
    DBWriter writer(outDb.c_str(), (outDb + std::string(".index")).c_str(), threads);
    writer.open();
    const size_t LINE_BUFFER_SIZE = 1000000;
#pragma omp parallel
    {
        char * lineBuffer = new char[LINE_BUFFER_SIZE];
        char * key = new char[255];
        std::string minusResultsOutString;
        minusResultsOutString.reserve(maxLineLength);
#pragma omp  for schedule(static)
        for (size_t id = 0; id < leftDbr.getSize(); id++) {
            std::map<unsigned int, bool> elementLookup;
            int thread_idx = 0;
            const char *leftData = leftDbr.getData(id);
            unsigned int leftDbKey = leftDbr.getDbKey(id);
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // fill element id look up with left side elementLookup
            {
                char *data = (char *) leftData;
                char *entry[255];
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int dbKey = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    // its an aln result (parse e-value)
                    if (columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                        evalue = strtod(entry[3], NULL);
                    }
                    if(evalue <= evalThreshold){
                        elementLookup[dbKey] = true;
                    }
                    data = Util::skipLine(data);
                }
            }
            // get all data for the leftDbkey from rightDbr
            // check if right ids are in elementsId
            char *data = rightDbr.getDataByDBKey(leftDbKey);
            char *entry[255];

            if (data != NULL) {
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int element = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    if (columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                        evalue = strtod(entry[3], NULL);
                    }
                    if(evalue <= evalThreshold) {
                        elementLookup[element] = false;
                    }
                    data = Util::skipLine(data);
                }
            }
            // write only elementLookup that are not found in rightDbr (id != UINT_MAX)
            {
                char *data = (char *) leftData;
                size_t dataLength = leftDbr.getSeqLens(id);

                while (*data != '\0') {
                    Util::parseKey(data, key);
                    if (!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
                        Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                        data = Util::skipLine(data);
                        continue;
                    }
                    unsigned int elementIdx = std::strtoul(key, NULL, 10);
                    if (elementLookup[elementIdx]) {
                        minusResultsOutString.append(lineBuffer);
                        minusResultsOutString.append("\n");
                    }
                    data = Util::skipLine(data);
                }
            }

            // create merged string
            if (minusResultsOutString.length() >= maxLineLength) {
                Debug(Debug::ERROR) << "ERROR: Buffer overflow at id: " << leftDbr.getDbKey(id) <<
                " during the merging.\n";
                Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << maxLineLength << " < " <<
                minusResultsOutString.length() <<
                ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
                continue; // read next id
            }
            // write result
            char *mergeResultsOutData = (char *) minusResultsOutString.c_str();
            writer.writeData(mergeResultsOutData, minusResultsOutString.length(), SSTR(leftDbKey).c_str(), thread_idx);
            minusResultsOutString.clear();
        }
        delete [] lineBuffer;
        delete [] key;
    } //end of OMP
    leftDbr.close();
    rightDbr.close();
    // close all reader
    writer.close();
    Debug(Debug::INFO) << "Stored results in " << outDb << "\n";
}

int subtractdbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    struct timeval start, end;
    gettimeofday(&start, NULL);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    dosubstractresult(par.db1, par.db2, par.db3, 1000000, par.evalProfile, par.threads);
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for profile substracting: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return 0;
}
