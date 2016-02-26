//
// Created by mad on 2/6/16.
//

#include <climits>
#include <list>
#include <vector>
#include "DBReader.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Util.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

void dosubstractresult(std::string leftDb, std::string rightDb, std::string outDb,
                       size_t maxLineLength, int threads)
{
    Debug(Debug::INFO) << "Remove " << rightDb << " ids from "<< leftDb << "\n";
    DBReader<unsigned int> leftDbr(leftDb.c_str(), (leftDb + std::string(".index")).c_str());
    leftDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> rightDbr(rightDb.c_str(), (rightDb + std::string(".index")).c_str());
    rightDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter writer(outDb.c_str(), (outDb + std::string(".index")).c_str(), threads);
    writer.open();
    const size_t LINE_BUFFER_SIZE = 1000000;
    const size_t MAX_ELEMENT_SIZE = 1000000;
#pragma omp parallel
    {
        unsigned int * elementsId = new unsigned int[MAX_ELEMENT_SIZE];
        char * lineBuffer = new char[LINE_BUFFER_SIZE];
        char * key = new char[255];
        std::string minusResultsOutString;
        minusResultsOutString.reserve(maxLineLength);
#pragma omp  for schedule(static)
        for (size_t id = 0; id < leftDbr.getSize(); id++) {
            int thread_idx = 0;
            const char *leftData = leftDbr.getData(id);
            unsigned int leftDbKey = leftDbr.getDbKey(id);
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            // fill element id look up with left side elements
            size_t elementCount = 0;
            {
                char *data = (char *) leftData;
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    elementsId[elementCount] = std::strtoul(key, NULL, 10);
                    data = Util::skipLine(data);
                    elementCount++;
                    if(elementCount >= MAX_ELEMENT_SIZE){
                        Debug(Debug::ERROR) << "Set with id="<< leftDbKey <<" has more elements than " << MAX_ELEMENT_SIZE << " increase MAX_ELEMENT_SIZE!\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
            }
            // get all data for the leftDbkey from rightDbr
            // check if right ids are in elementsId
            char *data = rightDbr.getDataByDBKey(leftDbKey);
            if (data != NULL) {
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int element = std::strtoul(key, NULL, 10);
                    for (size_t i = 0; i < elementCount; i++) {
                        if (elementsId[i] == element) {
                            elementsId[i] = UINT_MAX;
                            break;
                        }
                    }
                    data = Util::skipLine(data);
                }
            }
            // write only elements that are not found in rightDbr (id != UINT_MAX)
            {
                char *data = (char *) leftData;
                size_t dataLength = leftDbr.getSeqLens(id);

                size_t elementIdx = 0;
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    if (!Util::getLine(data, dataLength, lineBuffer, LINE_BUFFER_SIZE)) {
                        Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                        continue;
                    }
                    if (elementsId[elementIdx] == std::strtoul(key, NULL, 10)) {
                        minusResultsOutString.append(lineBuffer);
                        minusResultsOutString.append("\n");
                    }
                    elementIdx++;
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
            writer.write(mergeResultsOutData, minusResultsOutString.length(), SSTR(leftDbKey).c_str(), thread_idx);
            minusResultsOutString.clear();
        }
        delete [] elementsId;
        delete [] lineBuffer;
        delete [] key;
    } //end of OMP
    leftDbr.close();
    rightDbr.close();
    // close all reader
    writer.close();
    Debug(Debug::INFO) << "Stored results in " << outDb << "\n";
}

int substractresult(int argc,const char **argv)
{

    std::string usage;
    usage.append("Removes all entries with same ID from resultDbLeft contained in resultDbRight ( out = left - right) \n");
    usage.append("USAGE: <resultDbLeft> <resultDbRight> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.substractresult, 3);
    dosubstractresult(par.db1, par.db2, par.db3, 1000000, par.threads);

    return 0;
}
