
#include <iostream>
#include <time.h>
#include "Clustering.h"
#include "CommandDeclarations.h"


void printUsageClustering(){

    std::string usage("\nCalculates clustering of a sequence database based on Smith Waterman alignment scores with set cover algorithm.\n"
            "ATTENTION: ffindex sequence database should contain all the sequences participating in the clustering (queries and targets).");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: mmseqs_clu ffindexSeqDB ffindexAlnResultsDB ffindexOutDB [opts]\n"
             "-g              \t[file]\tgreedy clustering by sequence length (default: set cover clustering algorithm).\n"
             "-s              \t[float]\tMinimum sequence identity  (default = 0.0)\n"
             "-c              \t\tCheck clusters (default = off)\n"
             "-v              \t[int]\tVerbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO (default=3).\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexAlnDBBase, std::string* ffindexOutDBBase, std::string* ffindexSeqDBBase, int* clusteringMode, float* seqIdThr, int* verbosity, int* validateClustering){
    if (argc < 3){
        printUsageClustering();
        exit(EXIT_FAILURE);
    }
    ffindexSeqDBBase->assign(argv[1]);
    ffindexAlnDBBase->assign(argv[2]);
    ffindexOutDBBase->assign(argv[3]);

    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-g") == 0){
                *clusteringMode = Clustering::GREEDY;
                i++;
        }
        if (strcmp(argv[i], "-c") == 0){
            *validateClustering = 1;
            i++;
        }
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *seqIdThr = atof(argv[i]);
                i++;
            }
            else {
                printUsageClustering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-v") == 0){
            if (++i < argc){
                *verbosity = atoi(argv[i]);
                if (*verbosity < 0 || *verbosity > 3){
                    Debug(Debug::ERROR) << "Wrong value for verbosity, please choose one in the range [0:3].\n";
                    exit(1);
                }
                i++;
            }
            else {
                printUsageClustering();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsageClustering();
            Debug(Debug::ERROR) << "Wrong argument: " << argv[i-1] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

int cluster(int argc, const char ** argv)
{
    int verbosity = Debug::INFO;

    std::string seqDB = "";
    std::string alnDB = "";
    std::string outDB = "";
    int clusteringMode = Clustering::SET_COVER;
    float seqIdThr = 0.0;
    int validateClustering = 0;

    parseArgs(argc, argv, &alnDB, &outDB, &seqDB, &clusteringMode, &seqIdThr, &verbosity, &validateClustering);

    Debug::setDebugLevel(verbosity);

    Clustering* clu = new Clustering(seqDB, seqDB + ".index", alnDB, alnDB + ".index", outDB, outDB + ".index", seqIdThr, validateClustering);

    clu->run(clusteringMode);
    return 0;
}

