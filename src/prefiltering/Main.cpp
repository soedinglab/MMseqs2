#include <iostream>
#include <unistd.h>
#include <string>
#include <signal.h>
#include <execinfo.h>

#include "Prefiltering.h"

#ifdef OPENMP
#include <omp.h>
#endif

void mmseqs_debug_catch_signal(int sig_num)
{
  if(sig_num == SIGILL)
  {
    fprintf(stderr, "Your CPU does not support all the latest features that this version of mmseqs makes use of!\n"
                    "Please run on a newer machine.");
    exit(sig_num);
  }
  else
  {
    fprintf (stderr, "\n\n-------------------------8<-----------------------\nExiting on error!\n");
    fprintf (stderr, "Signal %d received\n",sig_num);
    perror("ERROR (can be bogus)");

/*    fprintf(stderr, "globalId: %d\n", globalId);
    fprintf(stderr, "getIndexTableExited: %d\n", getIndexTableExited);
    fprintf(stderr, "globalPos: %d\n", globalPos);*/

    fprintf(stderr, "Backtrace:");
    void *buffer[30];
    int nptrs = backtrace(buffer, 30);
    backtrace_symbols_fd(buffer, nptrs, 2);
    fprintf (stderr, "------------------------->8-----------------------\n\n"
        "Send the binary program that caused this error and the coredump (ls core.*).\n"
        "Or send the backtrace:"
        "\n$ gdb -ex=bt --batch PROGRAMM_NAME CORE_FILE\n"
        "If there is no core file, enable coredumps in your shell and run again:\n"
        "$ ulimit -c unlimited\n\n");
  }

  exit(1);
}

void mmseqs_cuticle_init()
{
  struct sigaction handler;
  handler.sa_handler = mmseqs_debug_catch_signal;
  sigemptyset(&handler.sa_mask);
  handler.sa_flags = 0;

  sigaction(SIGFPE, &handler, NULL);
  sigaction(SIGSEGV, &handler, NULL);
  sigaction(SIGBUS, &handler, NULL);
  sigaction(SIGABRT, &handler, NULL);
}

void printUsage(){

    std::string usage("\nCalculates similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de) & Martin Steinegger (Martin.Steinegger@campus.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_pref ffindexQueryDBBase ffindexTargetDBBase ffindexOutDBBase [opts]\n"
            "-m              \t[file]\tAmino acid substitution matrix file.\n"
            "-s              \t[float]\tSensitivity in the range [2:9] (default=7.2)\n"
            "-k              \t[int]\tk-mer size (default=6).\n"
            "-a              \t[int]\tAmino acid alphabet size (default=21).\n"
            "--z-score-thr   \t[float]\tZ-score threshold [default: 300.0]\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n"
            "--nucleotides   \t\tNucleotide sequences input.\n"
            "--max-res-num   \t[int]\tMaximum result sequences per query (default=100)\n"
            "--aa-bias-corr  \t\tLocal amino acid composition bias correction.\n"
            "--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n"); 
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* ffindexOutDBBase, std::string* scoringMatrixFile, float* sens, int* kmerSize, int* alphabetSize, float* zscoreThr, size_t* maxSeqLen, int* seqType, size_t* maxResListLen, bool* aaBiasCorrection, int* skip){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexQueryDBBase->assign(argv[1]);
    ffindexTargetDBBase->assign(argv[2]);
    ffindexOutDBBase->assign(argv[3]);

    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (*seqType == Sequence::NUCLEOTIDES){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
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
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *sens = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-k") == 0){
            if (++i < argc){
                *kmerSize = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-a") == 0){
            if (++i < argc){
                *alphabetSize = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--z-score-thr") == 0){
            if (++i < argc){
                *zscoreThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided" << argv[i-1] << "\n";
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
         else if (strcmp(argv[i], "--nucleotides") == 0){
            if (strcmp(scoringMatrixFile->c_str(), "") != 0){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }                                                           
            *seqType = Sequence::NUCLEOTIDES;
            i++;
        }
       else if (strcmp(argv[i], "--max-res-num") == 0){
            if (++i < argc){
                *maxResListLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--aa-bias-corr") == 0){
            *aaBiasCorrection = true;
            i++;
        }
        else if (strcmp(argv[i], "--skip") == 0){
            if (++i < argc){
                *skip = atoi(argv[i]);
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
    mmseqs_cuticle_init();

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int kmerSize =  6;
    int alphabetSize = 21;
    size_t maxSeqLen = 50000;
    size_t maxResListLen = 100;
    float sensitivity = 7.2f;
    int skip = 0;
    int seqType = Sequence::AMINO_ACIDS;
    bool aaBiasCorrection = false;
    float zscoreThr = 50.0f;

    std::string queryDB = "";
    std::string targetDB = "";
    std::string outDB = "";
    std::string scoringMatrixFile = "";

    parseArgs(argc, argv, &queryDB, &targetDB, &outDB, &scoringMatrixFile, &sensitivity, &kmerSize, &alphabetSize, &zscoreThr, &maxSeqLen, &seqType, &maxResListLen, &aaBiasCorrection, &skip);

    if (seqType == Sequence::NUCLEOTIDES)
        alphabetSize = 5;

    std::cout << "k-mer size: " << kmerSize << "\n";
    std::cout << "Alphabet size: " << alphabetSize << "\n";
    std::cout << "Sensitivity: " << sensitivity << "\n";
    std::cout << "Z-score threshold: " << zscoreThr << "\n";

    std::string queryDBIndex = queryDB + ".index";
    std::string targetDBIndex = targetDB + ".index";
    std::string outDBIndex = outDB + ".index";

    Prefiltering* pref = new Prefiltering(queryDB, queryDBIndex, targetDB, targetDBIndex, outDB, outDBIndex, scoringMatrixFile, sensitivity, kmerSize, alphabetSize, zscoreThr, maxSeqLen, seqType, aaBiasCorrection, skip);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    std::cout << "Starting prefiltering scores calculation.\n";
    pref->run(maxResListLen);

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return 0;
}
