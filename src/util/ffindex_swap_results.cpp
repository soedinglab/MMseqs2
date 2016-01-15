#include "DBReader.h"
#include "DBWriter.h"

#include "Debug.h"
#include <cstddef>
#include <stdio.h>
#include "Util.h"

void printUsageFFindexSwapResults(){
    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <ffindexDB> <fastaDB> [ffindexHeaderDB]\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexSeqDB, std::string* fastaOutDB) {
    if (argc < 2) {
        printUsageFFindexSwapResults();
        EXIT(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    fastaOutDB->assign(argv[2]);
}

int swapresults (int argc, const char * argv[])
{
    std::string ffindexResDB = "";
    std::string outDB = "";
    size_t splitSize = 1;
    parseArgs(argc, argv, &ffindexResDB, &outDB);
    Debug(Debug::WARNING) << "FFindex input file is " << ffindexResDB << "\n";
    DBWriter writer(outDB.c_str(), (outDB + ".index").c_str());
    writer.swapResults(ffindexResDB, splitSize);
    writer.close();
    return 0;
}
