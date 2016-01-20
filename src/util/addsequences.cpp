#include <sstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"

#include "Util.h"

int addsequences(int argc, const char **argv)
{
    std::string usage("Adds sequences in fasta format to an mmseqs clustering.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de).\n");
    usage.append("USAGE: <clusteredDB> <fastaInDB> <fastaOut>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.addSequences, 3);

    DBReader<unsigned int> clusters(par.db1.c_str(), par.db1Index.c_str());
    clusters.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> bodies(par.db2.c_str(), par.db2Index.c_str());
    bodies.open(DBReader<unsigned int>::NOSORT);

    std::string headerFilename(par.db2);
    headerFilename.append("_h");

    std::string headerIndexFilename(par.db2);
    headerIndexFilename.append("_h.index");

    DBReader<unsigned int> headers(headerFilename.c_str(), headerIndexFilename.c_str());
    headers.open(DBReader<unsigned int>::NOSORT);

    DBWriter msaOut(par.db3.c_str(), par.db3Index.c_str());
    msaOut.open();

    unsigned int* dataLengths = clusters.getSeqLens();

	for (size_t i = 0; i < clusters.getSize(); i++){
		std::ostringstream fastaStream;

        char* data = clusters.getData(i);

        if(par.minSequences > 1) {
            size_t entries = Util::count_lines(data, dataLengths[i] - 1);
            if(entries < par.minSequences) {
                continue;
            }
        }

        std::string entry;
        std::istringstream clusterEntries(data);
		while (std::getline(clusterEntries, entry)) {
            unsigned int entryId = strtoul(entry.c_str(), NULL, 10);
			char* header = headers.getDataByDBKey(entryId);
            if(header == NULL) {
                Debug(Debug::WARNING) << "Entry " << entry << " does not contain a header!" << "\n";
                continue;
            }

            char* body = bodies.getDataByDBKey(entryId);
            if(body == NULL) {
                Debug(Debug::WARNING) << "Entry " << entry << " does not contain a sequence!" << "\n";
                continue;
            }

            fastaStream << ">" << header << body;
		}

        std::string fasta = fastaStream.str();
        std::string key = SSTR(clusters.getDbKey(i));
        msaOut.write(fasta.c_str(), fasta.length(), key.c_str());
    }

	msaOut.close();
    headers.close();
    bodies.close();
    clusters.close();

    return EXIT_SUCCESS;
}
