#include <sstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createseqfiledb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> clusters(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    clusters.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> bodies(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    bodies.open(DBReader<unsigned int>::NOSORT);
    bodies.readMmapedDataInMemory();

    DBReader<unsigned int> headers(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headers.open(DBReader<unsigned int>::NOSORT);
    headers.readMmapedDataInMemory();

    DBWriter msaOut(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    msaOut.open();

    const size_t numClusters = clusters.getSize();
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < numClusters; ++i){
            std::string resultStr;
            char* data = clusters.getData(i, thread_idx);

            size_t entries = Util::countLines(data, clusters.getEntryLen(i) - 1);
            if (entries < (unsigned int) par.minSequences || entries > (unsigned int) par.maxSequences) {
                continue;
            }

            std::string entry;
            std::istringstream clusterEntries(data);
            size_t entries_num = 0;
            char dbKey[255 + 1];
            resultStr.clear();
            while (std::getline(clusterEntries, entry)) {
                entries_num++;
                Util::parseKey((char*)entry.c_str(), dbKey);
                const unsigned int entryId = (unsigned int) strtoul(dbKey, NULL, 10);

                char* header = headers.getDataByDBKey(entryId, thread_idx);
                if (header == NULL) {
                    Debug(Debug::WARNING) << "Entry " << entry << " does not contain a header!" << "\n";
                    continue;
                }

                char* body = bodies.getDataByDBKey(entryId, thread_idx);
                if (body == NULL) {
                    Debug(Debug::WARNING) << "Entry " << entry << " does not contain a sequence!" << "\n";
                    continue;
                }
                size_t lineLen = Util::skipLine(header) - header;
                std::string headerStr(header, lineLen);
                lineLen = Util::skipLine(body) - body;
                std::string bodyStr(body, lineLen);

                if (entries_num == 1 && par.hhFormat) {
                    std::string consensusHeader(headerStr);
                    resultStr.push_back('#');
                    resultStr.append(headerStr);
                    resultStr.push_back('>');
                    resultStr.append(Util::removeAfterFirstSpace(consensusHeader));
                    resultStr.append("_consensus\n");
                    resultStr.append(bodyStr);
                    resultStr.push_back('>');
                    resultStr.append(headerStr);
                    resultStr.append(bodyStr);
                } else {
                    resultStr.push_back('>');
                    resultStr.append(headerStr);
                    resultStr.append(bodyStr);
                }
            }

            unsigned int key = clusters.getDbKey(i);
            msaOut.writeData(resultStr.c_str(), resultStr.length(), key, thread_idx);
        }
    };

    msaOut.close();
    headers.close();
    bodies.close();
    clusters.close();

    return EXIT_SUCCESS;
}
