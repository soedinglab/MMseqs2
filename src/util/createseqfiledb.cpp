#include <sstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"

#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createseqfiledb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> clusters(par.db2.c_str(), par.db2Index.c_str());
    clusters.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> bodies(par.db1.c_str(), par.db1Index.c_str());
    bodies.open(DBReader<unsigned int>::NOSORT);
    bodies.readMmapedDataInMemory();

    std::string headerFilename(par.db1);
    headerFilename.append("_h");

    std::string headerIndexFilename(par.db1);
    headerIndexFilename.append("_h.index");

    DBReader<unsigned int> headers(headerFilename.c_str(), headerIndexFilename.c_str());
    headers.open(DBReader<unsigned int>::NOSORT);
    headers.readMmapedDataInMemory();

    DBWriter msaOut(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads));
    msaOut.open();

    unsigned int* dataLengths = clusters.getSeqLens();

    size_t entries = clusters.getSize();
#pragma omp for schedule(dynamic, 100)
	for (size_t i = 0; i < entries; ++i){
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
		std::ostringstream fastaStream;
        char* data = clusters.getData(i);

        size_t entries = Util::countLines(data, dataLengths[i] - 1);
        if (entries < (unsigned int) par.minSequences || entries > (unsigned int) par.maxSequences) {
            continue;
        }

        std::string entry;
        std::istringstream clusterEntries(data);
        size_t entries_num = 0;
		while (std::getline(clusterEntries, entry)) {
            entries_num++;

            char* rest;
            unsigned int entryId = strtoul(entry.c_str(), &rest, 10);
            if ((rest != entry.c_str() && *rest != '\0') || errno == ERANGE) {
                Debug(Debug::WARNING) << "Invalid entry in line " << entries_num << "!\n";
                continue;
            }

			char* header = headers.getDataByDBKey(entryId);
            if (header == NULL) {
                Debug(Debug::WARNING) << "Entry " << entry << " does not contain a header!" << "\n";
                continue;
            }

            char* body = bodies.getDataByDBKey(entryId);
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
        std::string key = SSTR(clusters.getDbKey(i));
        msaOut.writeData(fasta.c_str(), fasta.length(), key.c_str(), thread_idx);
    }

	msaOut.close();
    headers.close();
    bodies.close();
    clusters.close();

    return EXIT_SUCCESS;
}
