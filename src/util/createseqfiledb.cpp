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
    par.parseParameters(argc, argv, command, 3);

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

    unsigned int* dataLengths = clusters.getSeqLens();
    const size_t numClusters = clusters.getSize();
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < numClusters; ++i){
            std::ostringstream fastaStream;
            char* data = clusters.getData(i, thread_idx);

            size_t entries = Util::countLines(data, dataLengths[i] - 1);
            if (entries < (unsigned int) par.minSequences || entries > (unsigned int) par.maxSequences) {
                continue;
            }

            std::string entry;
            std::istringstream clusterEntries(data);
            size_t entries_num = 0;
            char dbKey[255 + 1];
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

                if (entries_num == 1 && par.hhFormat) {
                    std::string consensusHeader(header);
                    fastaStream << "#" << header
                                << ">" << Util::removeAfterFirstSpace(consensusHeader) << "_consensus\n" << body
                                << ">" << header << body;
                } else {
                    fastaStream << ">" << header << body;
                }
            }

            std::string fasta = fastaStream.str();
            unsigned int key = clusters.getDbKey(i);
            msaOut.writeData(fasta.c_str(), fasta.length(), key, thread_idx);
        }
    };

    msaOut.close();
    headers.close();
    bodies.close();
    clusters.close();

    return EXIT_SUCCESS;
}
