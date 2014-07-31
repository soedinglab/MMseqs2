#include <iostream>
#include <time.h>
#include <unistd.h>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include <execinfo.h>

#include "Alignment.h"

void printUsage(){

    std::string usage("\nCalculates Smith-Waterman alignment scores between all sequences in the query database and the sequences of the target database which passed the prefiltering.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_pref <queryDB> <targetDB> <prefResultsDB> <outDB> [opts]\n"
            "-e          \t[float]\tMaximum e-value (default=0.01).\n"
            "-c          \t[float]\tMinimum alignment coverage (default=0.8).\n"
            "-cpu        \t[int]\tNumber of cores used for the computation (default=all cores).\n"
            "--max-seq-len\t[int]\tMaximum sequence length (default=50000).\n"
            "--max-seqs\t[int]\tMaximum alignment results per query sequence (default=300).\n"
            "--max-rejected\t[int]\tMaximum rejected alignments before alignment calculation for a query is aborted. (default=INT_MAX)\n"
            "--nucleotides\t\tNucleotide sequences input.\n"
            "--sub-mat  \t[file]\tAmino acid substitution matrix file.\n"
            "-v         \t[int]\tVerbosity level: 0=NOTHING, 1=ERROR, 2=WARNING, 3=INFO (default=3).\n");
    Debug(Debug::INFO) << usage;
}

void parseArgs(int argc, char** argv, std::string* qseqDB, std::string* tseqDB, std::string* prefDB, std::string* matrixFile, std::string* outDB, double* evalThr, double* covThr, int* maxSeqLen, int* maxAlnNum, int* seqType, int* verbosity, int* maxRejected, int* threads){
    if (argc < 5){
        printUsage();
        exit(EXIT_FAILURE);
    }

    qseqDB->assign(argv[1]);
    tseqDB->assign(argv[2]);
    prefDB->assign(argv[3]);
    outDB->assign(argv[4]);
    int i = 5;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (*seqType == Sequence::NUCLEOTIDES){
                Debug(Debug::ERROR) << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            if (++i < argc){
                matrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--sub-mat") == 0){
            if (*seqType == Sequence::NUCLEOTIDES){
                Debug(Debug::ERROR) << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            if (++i < argc){
                matrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
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
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
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
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
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
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--max-rejected") == 0){
            if (++i < argc){
                *maxRejected = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--max-seqs") == 0){
            if (++i < argc){
                *maxAlnNum = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--nucleotides") == 0){
            *seqType = Sequence::NUCLEOTIDES;
            i++;
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
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-cpu") == 0){
            if (++i < argc){
                *threads = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }

    if (strcmp (matrixFile->c_str(), "") == 0){
        printUsage();
        Debug(Debug::ERROR) << "\nPlease provide a scoring matrix file. You can find scoring matrix files in $INSTALLDIR/data/.\n";
        exit(EXIT_FAILURE);
    }
}

bool compareHits (Matcher::result_t first, Matcher::result_t second){
    if (first.score > second.score)
        return true;
    return false;
}

// this is needed because with GCC4.7 omp_get_num_threads() returns just 1.
int omp_thread_count() {
    int n = 0;
#pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

int main(int argc, char **argv){

    int verbosity = Debug::INFO;
    int threads = 1;
#ifdef OPENMP
    threads = omp_thread_count();
#endif

    std::string qseqDB = "";
    std::string tseqDB = "";
    std::string prefDB = ""; 
    std::string outDB = "";

    double evalThr = 0.001;
    double covThr = 0.8;
    int maxSeqLen = 50000;
    int maxAlnNum = 300;
    int maxRejected = INT_MAX;
    int seqType = Sequence::AMINO_ACIDS;

    // get the path of the scoring matrix
    char* mmdir = getenv ("MMDIR");
    if (mmdir == 0){
        std::cerr << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        exit(1);
    }
    std::string matrixFile(mmdir);
    matrixFile = matrixFile + "/data/blosum62.out";

    Debug::setDebugLevel(Debug::INFO);

    Debug(Debug::WARNING) << "Program call:\n";
    for (int i = 0; i < argc; i++)
        Debug(Debug::WARNING) << argv[i] << " ";
    Debug(Debug::WARNING) << "\n\n";

    parseArgs(argc, argv, &qseqDB, &tseqDB, &prefDB, &matrixFile, &outDB, &evalThr, &covThr, &maxSeqLen, &maxAlnNum, &seqType, &verbosity, &maxRejected, &threads);

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

    Debug::setDebugLevel(verbosity);

    Debug(Debug::WARNING) 
        << "max. evalue:                       \t" << evalThr 
        << "\nmin. sequence coverage:          \t" << covThr 
        << "\nmax. sequence length:            \t" << maxSeqLen 
        << "\nmax. alignment results per query:\t" << maxAlnNum 
        << "\nmax rejected sequences per query:\t";
    if (maxRejected == INT_MAX)
        Debug(Debug::WARNING) << "off\n\n";
    else
        Debug(Debug::WARNING) << maxRejected << "\n\n";

    std::string qseqDBIndex = qseqDB + ".index";
    std::string tseqDBIndex = tseqDB + ".index";
    std::string prefDBIndex = prefDB + ".index";
    std::string outDBIndex = outDB+ ".index";

    Debug(Debug::WARNING) << "Init data structures...\n";
    Alignment* aln = new Alignment(qseqDB, qseqDBIndex, tseqDB, tseqDBIndex, prefDB, prefDBIndex, outDB, outDBIndex, matrixFile, evalThr, covThr, maxSeqLen, seqType);

    Debug(Debug::WARNING) << "Calculation of Smith-Waterman alignments.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    aln->run(maxAlnNum, maxRejected);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for alignments calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    delete aln;

    return 0;
}


