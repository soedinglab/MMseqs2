#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"


#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::ifstream
#include <vector>
#include <utility>      // std::pair

void printUsage(){
    std::string usage("\nMerge multiple ffindex files based on simular id into one file. \n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: ffindex_database_merge ffindexQueryDB ffindexOutputDB ffindexOutDB ffindexFILES*\n");
    usage.append("-s\t[int]\tcreate splitted index for n nodes\n");
    usage.append("-k              \t[int]\tk-mer size in the range [4:7] (default=6).\n");
    usage.append("-a              \t[int]\tAmino acid alphabet size (default=21).\n");
    usage.append("--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n");
    usage.append("--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n");


    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv,
               std::string* ffindexSeqDB,
               std::string* ffindexOutDB,
               int *kmerSize, int *alphabetSize,
               int *maxSeqLen, int * skip, int * splitt){
    if (argc < 2){
        printUsage();
        exit(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    ffindexOutDB->assign(argv[2]);

    
    int i = 3;
    while (i < argc){
        if (strcmp(argv[i], "-k") == 0){
            if (++i < argc){
                *splitt = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-k") == 0){
        if (++i < argc){
            *kmerSize = atoi(argv[i]);
            if (*kmerSize < 4 || *kmerSize > 7){
                Debug(Debug::ERROR) << "Please choose k in the range [4:7].\n";
                EXIT(EXIT_FAILURE);
            }
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    else if (strcmp(argv[i], "-a") == 0){
        if (++i < argc){
            *alphabetSize = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    else if (strcmp(argv[i], "--max-seq-len") == 0){
        if (++i < argc){
            *maxSeqLen = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }        else if (strcmp(argv[i], "--skip") == 0){
        if (++i < argc){
            *skip = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
            Debug(Debug::ERROR) << "WRONG PARAMETER";
            printUsage();
            exit(EXIT_FAILURE);
            break;
        }
    }
}



int main (int argc, const char * argv[])
{
    
    std::string seqDB = "";
    std::string outDB = "";

    int splitt = 1;
    int kmerSize =  6;
    int alphabetSize = 21;
    int maxSeqLen = 50000;
    int skip = 0;
    std::string scoringMatrixFile = "/Users/mad/Documents/workspace/mmseqs/data/blosum62.out";
    parseArgs(argc, argv, &seqDB, &outDB, &kmerSize, &alphabetSize, &maxSeqLen, &skip, &splitt);
    DBReader dbr(seqDB.c_str(), std::string(seqDB+".index").c_str());
    dbr.open(DBReader::SORT);

    
    BaseMatrix* subMat = Prefiltering::getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0f);
    Sequence seq(maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, subMat);
    
    PrefilteringIndexReader::createIndexFile(outDB, &dbr, &seq, splitt, alphabetSize, kmerSize, skip );

    
    // write code
    dbr.close();
    
    
    delete subMat;
    return 0;
}
