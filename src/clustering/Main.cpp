
#include <iostream>
#include <time.h>

#include "Clustering.h"

void printUsage(){

    std::string usage("\nCalculates clustering of a sequence database based on Smith Waterman alignment scores with set cover algorithm.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: mmseqs_clu ffindexSeqDB ffindexAlnResultsDB ffindexOutDB [opts]\n"
             "-g              \t[file]\tgreedy clustering by sequence length (default: set cover clustering algorithm).\n"
             "-s              \t[float]\tMinimum sequence identity  (default = 0.0)\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexAlnDBBase, std::string* ffindexOutDBBase, std::string* ffindexSeqDBBase, int* clusteringMode, float* seqIdThr){
    if (argc < 3){
        printUsage();
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
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *seqIdThr = atof(argv[i]);
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
            std::cerr << "Wrong argument: " << argv[i-1] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

int main(int argc, const char * argv[])
{

    std::string seqDB = "";
    std::string alnDB = "";
    std::string outDB = "";
    int clusteringMode = Clustering::SET_COVER;
    float seqIdThr = 0.0;

    parseArgs(argc, argv, &alnDB, &outDB, &seqDB, &clusteringMode, &seqIdThr);

    Clustering* clu = new Clustering(seqDB, seqDB + ".index", alnDB, alnDB + ".index", outDB, outDB + ".index", seqIdThr);

    clu->run(clusteringMode);
    
}

