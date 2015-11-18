#include <sstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

#include "Util.h"

void printUsageCusteringToFasta(){
    Debug(Debug::ERROR) << "\nConvert a mmseqs ffindex clustering to an clustering fasta format.\n";
    Debug(Debug::ERROR) << "Written by Milot Mirdita (milot@mirdita.de), Martin Steinegger (Martin.Steinegger@campus.lmu.de) and Maria Hauser (mhauser@genzentrum.lmu.de).\n\n";
    Debug(Debug::ERROR) << "USAGE:  <clusteredDB> <fastaHeaderInDB> <fastaBodyInDB> <msaOutDB>\n";
}

int clusteringtofastadb (int argc, const char **argv)
{
    if (argc < 4){
        printUsageCusteringToFasta();
        return EXIT_FAILURE;
    }

    std::string clusteredDB(argv[0]);
    std::string fastaHeaderInDB(argv[1]);
    std::string fastaBodyInDB(argv[2]);
    std::string msaOutDB(argv[3]);

    std::string clusteredDBIndex(std::string(clusteredDB).append(".index"));
    std::string fastaHeaderInDBIndex(std::string(fastaHeaderInDB).append(".index"));
    std::string fastaBodyInDBIndex(std::string(fastaBodyInDB).append(".index"));
    std::string msaOutDBIndex(std::string(msaOutDB).append(".index"));

    Debug(Debug::WARNING) << "Clustered database file: " << clusteredDB << "\n";
    Debug(Debug::WARNING) << "Fasta header input file: " << fastaHeaderInDB << "\n";
    Debug(Debug::WARNING) << "Fasta body input file: " << fastaBodyInDB << "\n";

    DBReader clusters(clusteredDB.c_str(), clusteredDBIndex.c_str());
    clusters.open(DBReader::NOSORT);

    DBReader headers(fastaHeaderInDB.c_str(), fastaHeaderInDBIndex.c_str());
    headers.open(DBReader::NOSORT);

    DBReader bodies(fastaBodyInDB.c_str(), fastaBodyInDBIndex.c_str());
    bodies.open(DBReader::NOSORT);

    DBWriter msaOut(msaOutDB.c_str(), msaOutDBIndex.c_str());
    msaOut.open();

    Debug(Debug::WARNING) << "Start writing results to " << msaOutDB << "\n";

	for (size_t i = 0; i < clusters.getSize(); i++){
		std::ostringstream fastaStream;

        std::string entry;
        std::istringstream clusterEntries(clusters.getData(i));
		while (std::getline(clusterEntries, entry)) {
            const char* cEntry = entry.c_str();
			char* header = headers.getDataByDBKey(cEntry);
            char* body   =  bodies.getDataByDBKey(cEntry);
            if(header == NULL || body == NULL) {
                Debug(Debug::WARNING) << "Entry (" << entry << ") is incomplete!" << "\n";
                continue;
            }
            fastaStream << ">" << header << body;
		}

        std::string fasta = fastaStream.str();
        std::string key = clusters.getDbKey(i);
        msaOut.write(fasta.c_str(), fasta.length(), key.c_str());
    }

	msaOut.close();
    bodies.close();
    headers.close();
    clusters.close();

    return EXIT_SUCCESS;
}
