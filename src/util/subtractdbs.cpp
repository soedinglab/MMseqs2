//
// Created by mad on 2/6/16.
//

#include <climits>
#include <list>
#include <vector>
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
    Debug(Debug::INFO) << "Remove " << rightDb << " ids from " << leftDb << "\n";
    DBReader<unsigned int> leftDbr(leftDb.c_str(), (leftDb + std::string(".index")).c_str());
    leftDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> rightDbr(rightDb.c_str(), (rightDb + std::string(".index")).c_str());
    rightDbr.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output databse: " << outDb << "\n";
    DBWriter writer(outDb.c_str(), (outDb + std::string(".index")).c_str(), threads);
    writer.open();
    const size_t LINE_BUFFER_SIZE = 1000000;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        char * lineBuffer = new char[LINE_BUFFER_SIZE];
        char * key = new char[255];
        std::string minusResultsOutString;
        minusResultsOutString.reserve(maxLineLength);
#pragma omp  for schedule(dynamic, 10)
        for (size_t id = 0; id < leftDbr.getSize(); id++) {
            std::map<unsigned int, bool> elementLookup;
            const char *leftData = leftDbr.getData(id);
            unsigned int leftDbKey = leftDbr.getDbKey(id);

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
                while (*data != '\0') {
                    char *start = data;
                    data = Util::skipLine(data);
                    Util::parseKey(start, key);
                    unsigned int elementIdx = std::strtoul(key, NULL, 10);
                    if (elementLookup[elementIdx]) {
                        minusResultsOutString.append(start, data - start);
                    }
                }
            }
            
            // write result
            char *mergeResultsOutData = (char *) minusResultsOutString.c_str();
            writer.writeData(mergeResultsOutData, minusResultsOutString.length(), leftDbKey, thread_idx);
            minusResultsOutString.clear();
        }
        delete [] lineBuffer;
        delete [] key;
    }
    writer.close();

    leftDbr.close();
    rightDbr.close();
}

int subtractdbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    dosubstractresult(par.db1, par.db2, par.db3, 1000000, par.evalProfile, par.threads);
    return EXIT_SUCCESS;
}
