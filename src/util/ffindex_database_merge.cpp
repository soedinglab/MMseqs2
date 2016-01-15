#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include <vector>
#include <utility>      // std::pair

void printUsageFFindexMergeDb(){
    std::string usage("\nMerge multiple ffindex files based on similar id into one file. \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <ffindexDB1> <ffindexDB2>\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv, 
	       std::string* ffindexSeqDB, 
	       std::string* ffindexOutDB, 
	       std::vector<std::pair<std::string,std::string> > * files){
    if (argc < 2){
        printUsageFFindexMergeDb();
        exit(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    ffindexOutDB->assign(argv[2]);

    int i = 3;
    while (i < argc){
    	files->push_back( std::make_pair<std::string, std::string>(std::string(argv[i]), std::string(argv[i])+".index" ) );
    	i++;
    }
}



int mergeffindex (int argc, const char * argv[])
{

    std::string seqDB = "";
    std::string outDB = ""; 
    std::vector<std::pair<std::string, std::string> > filenames;
    parseArgs(argc, argv, &seqDB, &outDB, &filenames); 

    DBReader<unsigned int> qdbr(seqDB.c_str(), std::string(seqDB+".index").c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str());
    writer.open();
    writer.mergeFiles(&qdbr, filenames, 1000000);
    writer.close();
    qdbr.close();
    return 0;
}
