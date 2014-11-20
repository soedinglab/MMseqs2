#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>

#include <sys/stat.h>


#include "DBReader.h"
#include "Debug.h"

void printUsageCusteringToFasta(){
    std::string usage("\nConvert a mmseqs ffindex clustering to an clustering fasta format.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de) & Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE:  <clusteredDB> <fastaHeaderInDB> <fastaBodyInDB> <msaOutDB>\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv,
               std::string* clusteredDB,
               std::string* fastaHeaderInDB,
               std::string* fastaBodyInDB,
               std::string* msaOutDB){
	if (argc < 5){
        printUsageCusteringToFasta();
        exit(EXIT_FAILURE);
    }

    clusteredDB->assign(argv[1]);
    fastaHeaderInDB->assign(argv[2]);
    fastaBodyInDB->assign(argv[3]);
    msaOutDB->assign(argv[4]);
}

FILE* openFileOrDie(const char * fileName, const char * mode) {
	struct stat st;
	FILE* file;
    if(stat(fileName, &st) == 0) { errno = EEXIST; perror(fileName); exit(EXIT_FAILURE); }
    file = fopen(fileName, mode);
    if(file == NULL) { perror(fileName); exit(EXIT_FAILURE); }
	return file;
}

int clusteringtofastadb (int argc, const char **argv)
{
    
    std::string clusteredDB = "";
    std::string fastaHeaderInDB = "";
    std::string fastaBodyInDB = "";
    std::string msaOutDB = "";

    parseArgs(argc, argv, &clusteredDB, &fastaHeaderInDB, &fastaBodyInDB, &msaOutDB);
    Debug(Debug::WARNING) << "Clustered database file: " << clusteredDB << "\n";
    DBReader clusters(clusteredDB.c_str(), std::string(clusteredDB+".index").c_str());
    clusters.open(DBReader::NOSORT);
    
	Debug(Debug::WARNING) << "Fasta header input file: " << fastaHeaderInDB << "\n";
    DBReader headers(fastaHeaderInDB.c_str(), std::string(fastaHeaderInDB+".index").c_str());
    headers.open(DBReader::NOSORT);
    
	Debug(Debug::WARNING) << "Fasta body input file: " << fastaBodyInDB << "\n";
    DBReader bodies(fastaBodyInDB.c_str(), std::string(fastaBodyInDB+".index").c_str());
    bodies.open(DBReader::NOSORT);
    
	std::string msaOutIndex = std::string(msaOutDB + ".index");

	FILE* msaData  = openFileOrDie(msaOutDB.c_str(), "w");
	FILE* msaIndex = openFileOrDie(msaOutIndex.c_str(), "w+");

	char header_start[] = {'>'};
    char newline[] = {'\n'};
    Debug(Debug::WARNING) << "Start writing file to " << msaOutDB << "\n";
    
	size_t offset = 0;
	for (size_t i = 0; i < clusters.getSize(); i++){
        char* clusterKey     = clusters.getDbKey(i);
		std::istringstream clusterEntries(clusters.getData(i));
		std::string entry;
		std::string fasta;
		while (std::getline(clusterEntries, entry)) {
	    	char* cEntry = const_cast<char *>(entry.c_str());
			char* header = headers.getDataByDBKey(cEntry);
        	char* body   =  bodies.getDataByDBKey(cEntry);
			fasta += "> " + std::string(header) + std::string(body);
		}

		ffindex_insert_memory(msaData, msaIndex, &offset, const_cast<char *>(fasta.c_str()), fasta.length(), clusterKey);
    }
    Debug(Debug::WARNING) << "Done." << "\n";

	fclose(msaData);
	fclose(msaIndex);
    bodies.close();
    headers.close();
    clusters.close();

    return 0;
}
