#include <iostream>
#include <unistd.h>
#include <string>
#include <signal.h>
#include <execinfo.h>

#include "TimeTest.h"

#ifdef OPENMP
#include <omp.h>
#endif

void printUsage(){

    std::string usage("\nRuns time benchmark on a database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_pref ffindexDBBase outputFile [opts]\n"
            "-m              \t[file]\tAmino acid substitution matrix file.\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexDBBase, std::string* logFile, std::string* scoringMatrixFile, size_t* maxSeqLen){
    if (argc < 3){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexDBBase->assign(argv[1]);
    logFile->assign(argv[2]);

    int i = 3;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (++i < argc){
                scoringMatrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--max-seq-len") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }

    if (strcmp (scoringMatrixFile->c_str(), "") == 0){
        printUsage();
        std::cerr << "\nPlease provide a scoring matrix file. You can find scoring matrix files in $INSTALLDIR/data/.\n";
        exit(EXIT_FAILURE);
    }
}



int main (int argc, const char * argv[])
{
    size_t maxSeqLen = 50000;

    std::string db = "";
    std::string scoringMatrixFile = "";
    std::string logFile = "";

    parseArgs(argc, argv, &db, &logFile, &scoringMatrixFile, &maxSeqLen);

    std::string dbIndex = db + ".index";

    TimeTest* tt = new TimeTest(db,
            dbIndex,
            db,
            dbIndex,
            scoringMatrixFile,
            maxSeqLen,
            logFile);

    tt->runTimeTest();
    
    delete tt;
}
