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

#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

void dosubstractresult(std::string leftDb, std::string rightDb, std::string outDb,
                       size_t maxLineLength, double evalThreshold, int threads, int compressed)
{
    Debug(Debug::INFO) << "Remove " << rightDb << " ids from " << leftDb << "\n";
    DBReader<unsigned int> leftDbr(leftDb.c_str(), (leftDb + std::string(".index")).c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    leftDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> rightDbr(rightDb.c_str(), (rightDb + std::string(".index")).c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    rightDbr.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output databse: " << outDb << "\n";
    DBWriter writer(outDb.c_str(), (outDb + std::string(".index")).c_str(), threads, compressed, leftDbr.getDbtype());
    writer.open();
    const size_t LINE_BUFFER_SIZE = 1000000;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char *entry[255];
        char * lineBuffer = new char[LINE_BUFFER_SIZE];
        char * key = new char[255];
        std::string minusResultsOutString;
        minusResultsOutString.reserve(maxLineLength);

#pragma omp  for schedule(dynamic, 10)
        for (size_t id = 0; id < leftDbr.getSize(); id++) {
            std::map<unsigned int, bool> elementLookup;
            const char *leftData = leftDbr.getData(id, thread_idx);
            unsigned int leftDbKey = leftDbr.getDbKey(id);

            // fill element id look up with left side elementLookup
            {
                char *data = (char *) leftData;
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int dbKey = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    // its an aln result (parse e-value)
                    if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
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
            char *data = rightDbr.getDataByDBKey(leftDbKey, thread_idx);

            if (data != NULL) {
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int element = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
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
    par.parseParameters(argc, argv, command, true, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    dosubstractresult(par.db1, par.db2, par.db3, 1000000, par.evalProfile, par.threads, par.compressed);
    return EXIT_SUCCESS;
}
