#include <sstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"

#include "Util.h"

int clusteringtofastadb (int argc, const char **argv)
{
    std::string usage("Convert a mmseqs ffindex clustering to an clustering fasta format.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de).\n");
    usage.append("USAGE:  <clusteredDB> <fastaHeaderInDB> <fastaBodyInDB> <fastaOut>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 4);

    DBReader<unsigned int> clusters(par.db1.c_str(), par.db1Index.c_str());
    clusters.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> headers(par.db2.c_str(), par.db2Index.c_str());
    headers.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> bodies(par.db3.c_str(), par.db3Index.c_str());
    bodies.open(DBReader<unsigned int>::NOSORT);

    DBWriter msaOut(par.db4.c_str(), par.db4Index.c_str());
    msaOut.open();

	for (size_t i = 0; i < clusters.getSize(); i++){
		std::ostringstream fastaStream;

        std::string entry;
        std::istringstream clusterEntries(clusters.getData(i));
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
    bodies.close();
    headers.close();
    clusters.close();

    return EXIT_SUCCESS;
}
