#include <iostream>
#include <time.h>
#include <unistd.h>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include <execinfo.h>

#include "Alignment.h"

void printUsage(){

    std::string usage("\nCalculates Smith-Waterman alignment scores.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: kClust2_pref ffindexSequenceDBBase ffindexPrefilteringDBBase ffindexOutDBBase [opts]\n"
            "-m\t[file]\tAmino acid substitution matrix file.\n"
            "-e\t[float]\tMaximum e-value (default=0.01).\n"
            "-c\t[float]\tMinimum alignment coverage (default=0.8).\n"
            "-s\t[int]\tMaximum sequence length (default=50000).\n"
            "-r\t[int]\tMaximum result alignment number per query sequence (default=10).\n"
            "-n\t\tNucleotide sequences input.\n");
    std::cout << usage;
}

void parseArgs(int argc, char** argv, std::string* seqDB, std::string* prefDB, std::string* matrixFile, std::string* outDB, double* evalThr, double* covThr, int* maxSeqLen, int* maxAlnNum, int* seqType){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    seqDB->assign(argv[1]);
    prefDB->assign(argv[2]);
    outDB->assign(argv[3]);
    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (*seqType == Sequence::NUCLEOTIDES){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            if (++i < argc){
                matrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -m\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-e") == 0){
            if (++i < argc){
                *evalThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -e\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-c") == 0){
            if (++i < argc){
                *covThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -c\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -s\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-r") == 0){
            if (++i < argc){
                *maxAlnNum = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -a\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-n") == 0){
            if (strcmp(matrixFile->c_str(), "") != 0){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            *seqType = Sequence::NUCLEOTIDES;
            i++;
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

bool compareHits (Matcher::result_t first, Matcher::result_t second){
    if (first.score > second.score)
        return true;
    return false;
}

int main(int argc, char **argv){

    std::string seqDB = "";
    std::string prefDB = ""; 
    std::string matrixFile = ""; 
    std::string outDB = "";

    double evalThr = 0.001;
    double covThr = 0.8;
    int maxSeqLen = 50000;
    int maxAlnNum = 10;
    int seqType = Sequence::AMINO_ACIDS;

    parseArgs(argc, argv, &seqDB, &prefDB, &matrixFile, &outDB, &evalThr, &covThr, &maxSeqLen, &maxAlnNum, &seqType);

    std::cout << "max. evalue:\t" << evalThr 
        << "\nmin. sequence coverage:\t" << covThr 
        << "\nmax. sequence length:\t" << maxSeqLen 
        << "\nmax. alignment number per query:\t" << maxAlnNum << "\n\n";

    Alignment* aln = new Alignment(seqDB, prefDB, outDB, matrixFile, evalThr, covThr, maxSeqLen, seqType);

    std::cout << "Calculation of Smith-Waterman alignments...\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    aln->run(maxAlnNum);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for alignments calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return 0;
}


