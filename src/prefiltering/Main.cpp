#include <iostream>
#include <time.h>
#include <unistd.h>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include <execinfo.h>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/SubstitutionMatrix.h"
#include "../commons/Sequence.h"
#include "../commons/NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

int globalId;

void kClust2_debug_catch_signal(int sig_num)
{
  if(sig_num == SIGILL)
  {
    fprintf(stderr, "Your CPU does not support all the latest features that this version of kClust2 makes use of!\n"
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

void kClust2_cuticle_init()
{
  struct sigaction handler;
  handler.sa_handler = kClust2_debug_catch_signal;
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
    usage.append("USAGE: kClust2_pref ffindexQueryDBBase ffindexTargetDBBase ffindexOutDBBase [opts]\n"
            "-m              \t[file]\tAmino acid substitution matrix file.\n"
            "-s              \t[float]\tSensitivity in the range [2:9] (default=7.2)\n"
            "-k              \t[int]\tk-mer size (default=6).\n"
            "-a              \t[int]\tAmino acid alphabet size (default=21).\n"
            "--z-score-thr   \t[float]\tZ-score threshold [default: 300.0]\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n"
            "--nucleotides   \t\tNucleotide sequences input.\n"
            "--max-res-num   \t[int]\tMaximum result sequences per query (default=100)\n"
            "--thr-calc-method\t[int]\tMethod for the prefiltering threshold calculation.\n\t\t1: derived from the scores distribution for each query sequence, \n\t\t2: derived from the score distribution of the reverted sequence\n\t\t(default=1)\n"
            "--aa-bias-corr  \t\tLocal amino acid composition bias correction.\n"
            "--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n"); 
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* ffindexOutDBBase, std::string* scoringMatrixFile, float* sens, int* kmerSize, int* alphabetSize, float* zscoreThr, size_t* maxSeqLen, int* seqType, size_t* maxResListLen, int* thresholdCalcMethod, bool* aaBiasCorrection, int* skip){
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
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided" << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--thr-calc-method") == 0){
            if (++i < argc){
                *thresholdCalcMethod = atoi(argv[i]);
                if (*thresholdCalcMethod < 1 && *thresholdCalcMethod > 2){
                    std::cerr << "Wrong threshold calculation method.\n";
                    exit(EXIT_FAILURE);
                }
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
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
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

BaseMatrix* getSubstitutionMatrix(std::string scoringMatrixFile, int alphabetSize){
    std::cout << "Substitution matrices...\n";
    BaseMatrix* subMat;
    if (alphabetSize < 21){
        SubstitutionMatrix* sMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), 8.0);
        subMat = new ReducedMatrix(sMat->probMatrix, alphabetSize);
    }
    else
        subMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), 8.0);

    return subMat;
}


IndexTable* getIndexTable (int alphabetSize, int kmerSize, DBReader* dbr, Sequence* seq, int dbSize, int skip){

    std::cout << "Index table: counting k-mers...\n";
    // fill and init the index table
    IndexTable* indexTable = new IndexTable(alphabetSize, kmerSize, skip);

    for (int id = 0; id < dbSize; id++){
        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(id, dbr->getDbKey(id), seqData);
        indexTable->addKmerCount(seq);
    }

    std::cout << "Index table init...\n";
    indexTable->init();

    std::cout << "Index table: fill...\n";
    for (int id = 0; id < dbSize; id++){
        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";

        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(id, dbr->getDbKey(id), seqData);
        indexTable->addSequence(seq);
    }

/*    std::cout << "Index table: removing duplicate entries...\n";
    indexTable->removeDuplicateEntries();
*/
    return indexTable;
}

/* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
 * K-mer similarity threshold is set to meet a certain DB match probability.
 * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
 */
std::pair<short,double> setKmerThreshold (DBReader* dbr, Sequence** seqs, BaseMatrix* subMat, ExtendedSubstitutionMatrix* _2merSubMatrix, ExtendedSubstitutionMatrix* _3merSubMatrix, int kmerSize, int maxSeqLen, double targetKmerMatchProb, double toleratedDeviation){

    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif

    size_t targetDbSize = dbr->getSize();
    if (targetDbSize > 100000)
        targetDbSize = 100000;
    IndexTable* indexTable = getIndexTable(subMat->alphabetSize, kmerSize, dbr, seqs[0], targetDbSize, 0); 

    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

    int targetSeqLenSum = 0;
    for (size_t i = 0; i < targetDbSize; i++)
        targetSeqLenSum += dbr->getSeqLens()[i];

    // generate a small random sequence set for testing 
    int querySetSize = dbr->getSize();
    if (querySetSize > 1000)
        querySetSize = 1000;
    int* querySeqs = new int[querySetSize];
    srand(1);
    for (int i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % dbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability 
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb
    short kmerThrMin = 20;
    short kmerThrMax = 200;
    short kmerThrMid;

    size_t dbMatchesSum;
    size_t querySeqLenSum;
    size_t dbMatchesExp_pc;
    // 1000 * 350 * 100000 * 350
    size_t lenSum_pc = 12250000000000;

    double kmersPerPos = 0.0;
    double kmerMatchProb;

    double kmerMatchProbMax = targetKmerMatchProb + (toleratedDeviation * targetKmerMatchProb);
    double kmerMatchProbMin = targetKmerMatchProb - (toleratedDeviation * targetKmerMatchProb);
    std::cout << "Searching for a k-mer threshold with a k-mer match probability per query position within range [" << kmerMatchProbMin << ":" << kmerMatchProbMax << "].\n";

    // adjust k-mer list length threshold
    while (kmerThrMax >= kmerThrMin){
        dbMatchesSum = 0;
        querySeqLenSum = 0;

        kmerThrMid = kmerThrMin + (kmerThrMax - kmerThrMin)*3/4;

        std::cout << "k-mer threshold range: [" << kmerThrMin  << ":" << kmerThrMax << "], trying threshold " << kmerThrMid << "\n";
        // determine k-mer match probability for kmerThrMid
#pragma omp parallel for schedule(static) 
        for (int i = 0; i < threads; i++){
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // set a current k-mer list length threshold and a high prefitlering threshold (we don't need the prefiltering results in this test run)
            matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix, indexTable, dbr->getSeqLens(), kmerThrMid, 1.0, kmerSize, dbr->getSize(), false, maxSeqLen, 500.0);
        }

#pragma omp parallel for schedule(dynamic, 10) reduction (+: dbMatchesSum, querySeqLenSum, kmersPerPos)
        for (int i = 0; i < querySetSize; i++){
            int id = querySeqs[i];

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char* seqData = dbr->getData(id);
            seqs[thread_idx]->mapSequence(id, dbr->getDbKey(id), seqData);

            matchers[thread_idx]->matchQuery(seqs[thread_idx]);

            kmersPerPos += seqs[thread_idx]->stats->kmersPerPos;
            dbMatchesSum += seqs[thread_idx]->stats->dbMatches;
            querySeqLenSum += seqs[thread_idx]->L;

        }

        kmersPerPos /= (double)querySetSize;

        // add pseudo-counts
        dbMatchesExp_pc = (size_t)(((double)lenSum_pc) * kmersPerPos * pow((1.0/((double)(subMat->alphabetSize-1))), kmerSize));

        // match probability with pseudocounts
        kmerMatchProb = ((double)dbMatchesSum + dbMatchesExp_pc) / ((double) (querySeqLenSum * targetSeqLenSum + lenSum_pc)); 

        for (int j = 0; j < threads; j++){
            delete matchers[j];
        }

        std::cout << "k-mer match probability: " << kmerMatchProb << "\n";
        if (kmerMatchProb < kmerMatchProbMin)
            kmerThrMax = kmerThrMid - 1;
        else if (kmerMatchProb > kmerMatchProbMax)
            kmerThrMin = kmerThrMid + 1;
        else if (kmerMatchProb > kmerMatchProbMin && kmerMatchProb < kmerMatchProbMax){
            // delete data structures used before returning
            delete[] querySeqs;
            delete[] matchers;
            delete indexTable;
            return std::pair<short, double> (kmerThrMid, kmerMatchProb);
        }
    }
    delete[] querySeqs;
    delete[] matchers;
    delete indexTable;
    return std::pair<short, double> (0, 0.0);
}


int main (int argc, const char * argv[])
{
    kClust2_cuticle_init();

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int kmerSize =  6;
    int alphabetSize = 21;
    size_t maxSeqLen = 50000;
    size_t maxResListLen = 100;
    size_t BUFFER_SIZE = 1000000;
    float sensitivity = 7.2f;
    int skip = 0;
    int seqType = Sequence::AMINO_ACIDS;
    int thresholdCalcMethod = 1;
    bool aaBiasCorrection = false;
    float zscoreThr = 300.0f;

    std::string queryDB = "";
    std::string queryDBIndex = "";
    std::string targetDB = "";
    std::string targetDBIndex = "";
    std::string outDB = "";
    std::string outDBIndex = "";
    std::string scoringMatrixFile = "";

    parseArgs(argc, argv, &queryDB, &targetDB, &outDB, &scoringMatrixFile, &sensitivity, &kmerSize, &alphabetSize, &zscoreThr, &maxSeqLen, &seqType, &maxResListLen, &thresholdCalcMethod, &aaBiasCorrection, &skip);

    if (seqType == Sequence::NUCLEOTIDES)
        alphabetSize = 5;

    std::cout << "k-mer size: " << kmerSize << "\n";
    std::cout << "Alphabet size: " << alphabetSize << "\n";
    std::cout << "Sensitivity: " << sensitivity << "\n";
    std::cout << "Z-score threshold: " << zscoreThr << "\n";

    queryDBIndex = queryDB + ".index";
    targetDBIndex = targetDB + ".index";
    outDBIndex = outDB + ".index";

    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    std::cout << "Using " << threads << " threads.\n";
#endif
    std::cout << "\n";

    std::cout << "Init data structures...\n";
    DBReader* qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader::NOSORT);

    DBReader* tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open(DBReader::SORT);

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    std::cout << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";
    std::cout << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

    // init the substitution matrices
    BaseMatrix* subMat;
    if (seqType == Sequence::AMINO_ACIDS)
        subMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize);
    else
        subMat = new NucleotideMatrix();

    ExtendedSubstitutionMatrix* _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    ExtendedSubstitutionMatrix* _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

    // init all thread-specific data structures 
    Sequence** seqs = new Sequence*[threads];
    std::list<int>** reslens = new std::list<int>*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType);
        reslens[thread_idx] = new std::list<int>();
    }

    int queryDBSize = qdbr->getSize();
    char** outBuffers = new char*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

    // set the k-mer similarity threshold
    float p_match =  pow(2.0, sensitivity) * 1.0e-08 * (float) (skip + 1); // old target value: 1.5e-06, reached with sens = 7.2 approximately
    std::cout << "\nAdjusting k-mer similarity threshold within +-10% deviation from the target k-mer match probability (target probability = " << p_match << ")...\n";
    std::pair<short, double> ret = setKmerThreshold (tdbr, seqs, subMat, _2merSubMatrix, _3merSubMatrix, kmerSize, maxSeqLen, p_match, 0.1);
    short kmerThr = ret.first; //174; //ret.first; // 103;
    double kmerMatchProb =  ret.second; //1.16383e-09; //ret.second; // 1.57506e-06;
    if (kmerThr == 0.0){
        std::cout << "Could not set the probability within +-10% deviation. Trying +-15% deviation.\n";
        ret = setKmerThreshold (tdbr, seqs, subMat, _2merSubMatrix, _3merSubMatrix, kmerSize, maxSeqLen, p_match, 0.15);
        kmerThr = ret.first;
        kmerMatchProb = ret.second;
    }
    if (kmerThr == 0.0){
        std::cout << "ERROR: Could not adjust the k-mer list length threshold!\n";
        std::cout << "Please report this error to the developers Maria Hauser mhauser@genzentrum.lmu.de and Martin Steinegger Martin.Steinegger@campus.lmu.de\n";
        std::cout << "In the meantime, try to change your parameters k, a, and/or sensitivity.\n";
        exit(1);
    }
    std::cout << "... done.\n";

    std::cout << "k-mer similarity threshold: " << kmerThr << "\n"; 
    std::cout << "k-mer match probability: " << kmerMatchProb << "\n\n";

    // initialise the index table and the matcher structures for the database
    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType);
    IndexTable* indexTable = getIndexTable(alphabetSize, kmerSize, tdbr, seq, tdbr->getSize(), skip);
    delete seq;

    std::cout << "Initializing data structures...";
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif 
        matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix, indexTable, tdbr->getSeqLens(), kmerThr, kmerMatchProb, kmerSize, tdbr->getSize(), aaBiasCorrection, maxSeqLen, zscoreThr);
    }
    std::cout << "... done.\n";

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    // calculate prefiltering scores for each sequence in the database
    std::cout << "Calculating prefiltering scores!\n";
    size_t kmersPerPos = 0.0;
    size_t dbMatches = 0;
    size_t resSize = 0;
    int empty = 0;

#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, empty, dbMatches)
    for (int id = 0; id < queryDBSize; id++){

        if (id % 1000000 == 0 && id > 0){
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
            fflush(stdout);
        }
        else if (id % 10000 == 0 && id > 0) {
            std::cout << ".";
            fflush(stdout);
        }

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get query sequence
        char* seqData = qdbr->getData(id);
        seqs[thread_idx]->mapSequence(id, qdbr->getDbKey(id), seqData);

        // calculate prefitlering results
        std::list<hit_t>* prefResults;
        if (thresholdCalcMethod == 1)
            prefResults = matchers[thread_idx]->matchQuery(seqs[thread_idx]);
        else
            prefResults = matchers[thread_idx]->matchQueryRevSeq(seqs[thread_idx]);
       
        // write prefiltering results to a string
        std::stringstream prefResultsOut;
        size_t l = 0;
        for (std::list<hit_t>::iterator iter = prefResults->begin(); iter != prefResults->end(); iter++){
            if (iter->seqId >= tdbr->getSize()){
                std::cout << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id)<< " -> " << iter->seqId << "\t" << iter->prefScore << "\n";
            }
            prefResultsOut << tdbr->getDbKey(iter->seqId) << "\t" << iter->prefScore << "\t" << iter->eval << "\n";
            l++;
            // maximum allowed result list length is reached
            if (l == maxResListLen)
                break;
        }

        // write prefiltering results string to ffindex database
        std::string prefResultsOutString = prefResultsOut.str();
        const char* prefResultsOutData = prefResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(prefResultsOutData)){
            std::cerr << "Tried to process the prefiltering list for the query " << qdbr->getDbKey(id) << " , the length of the list = " << prefResults->size() << "\n";
            std::cerr << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffers[thread_idx], prefResultsOutData, prefResultsOutString.length()*sizeof(char));
        dbw->write(outBuffers[thread_idx], prefResultsOutString.length(), qdbr->getDbKey(id), thread_idx);

        // update statistics counters
        if (prefResults->size() == 0)
            empty++;
        kmersPerPos += (size_t) seqs[thread_idx]->stats->kmersPerPos;
        dbMatches += seqs[thread_idx]->stats->dbMatches;
        resSize += prefResults->size();
        reslens[thread_idx]->push_back(prefResults->size());

    }
    std::cout << "\n\n";

    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }

    // calculate and print statistics
    kmersPerPos = (long double) kmersPerPos / queryDBSize;
    size_t dbMatchesPerSeq = dbMatches/(size_t)queryDBSize;
    size_t prefPassedPerSeq = resSize/(size_t)queryDBSize;
    std::cout << kmersPerPos << " k-mers per position.\n";
    std::cout << dbMatchesPerSeq << " DB matches per sequence.\n";
    std::cout << prefPassedPerSeq << " sequences passed prefiltering per query sequence";
    if (prefPassedPerSeq > 100)
        std::cout << " (ATTENTION: max. " << maxResListLen << " best scoring sequences were written to the output prefiltering database).\n";
    else
        std::cout << ".\n";

    int mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
    std::advance(it, mid);
    std::cout << "Median result list size: " << *it << "\n";
    std::cout << empty << " sequences with 0 size result lists.\n";

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    // merge output ffindex databases
    std::cout << "Merging the results...\n";
    qdbr->close();
    if (queryDB.compare(targetDB) != 0)
        tdbr->close();
    dbw->close();
    std::cout << "done.\n";

    for (int i = 0; i < threads; i++){
        delete seqs[i];
        delete[] outBuffers[i];
        delete matchers[i];
        delete reslens[i];
    }
    delete[] seqs;
    delete[] outBuffers;
    delete[] matchers;
    delete[] reslens;

    delete indexTable;
    delete subMat;
    delete _2merSubMatrix;
    delete _3merSubMatrix;

    return 0;
}
