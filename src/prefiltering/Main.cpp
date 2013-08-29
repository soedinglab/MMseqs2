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
#include "ExtendedSubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

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

  abort();
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
    usage.append("USAGE: kClust2_pref ffindexQueryDBBase ffindexTargetDBBase scoringMatrixFile ffindexOutDBBase [opts]\n"
            "-t\t[float]\tPrefiltering threshold (minimum half bits per query position).\n"
            "-k\t[int]\tk-mer size (default=6).\n"
            "-a\t[int]\tAmino acid alphabet size (default=21).\n"
            "-m\t[int]\tMaximum sequence length (default=50000).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* scoringMatrixFile, std::string* ffindexOutDBBase, int* kmerSize, int* alphabetSize, size_t* maxSeqLen){
    if (argc < 5){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexQueryDBBase->assign(argv[1]);
    ffindexTargetDBBase->assign(argv[2]);
    scoringMatrixFile->assign(argv[3]);
    ffindexOutDBBase->assign(argv[4]);
    int i = 5;
    while (i < argc){
        if (strcmp(argv[i], "-k") == 0){
            if (++i < argc){
                *kmerSize = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -k\n";
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
                std::cerr << "No value provided for -a\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-m") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -m\n";
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
        SubstitutionMatrix* sMat = new SubstitutionMatrix (scoringMatrixFile.c_str());
        subMat = new ReducedMatrix(sMat->probMatrix, alphabetSize);
    }
    else
        subMat = new SubstitutionMatrix (scoringMatrixFile.c_str());
    //BaseMatrix::print(subMat->subMatrix, subMat->int2aa, subMat->alphabetSize);

    return subMat;
}


IndexTable* getIndexTable (int alphabetSize, int kmerSize, DBReader* dbr, Sequence* seq, int dbSize){

    std::cout << "Index table: counting k-mers...\n";
    // fill and init the index table
    IndexTable* indexTable = new IndexTable(alphabetSize, kmerSize);

    for (int id = 0; id < dbSize; id++){
        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
        seq->id = id;
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(seqData);
        indexTable->addKmerCount(seq);
    }

    std::cout << "Index table init...\n";
    indexTable->init();

    std::cout << "Index table: fill...\n";
    for (int id = 0; id < dbSize; id++){
        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";

        seq->id = id;
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(seqData);
        indexTable->addSequence(seq);
    }

    return indexTable;
}

/* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
 * K-mer similarity threshold is set to meet a certain DB match probability.
 * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
 */
short setKmerThreshold (DBReader* dbr, Sequence** seqs, BaseMatrix* subMat, int kmerSize, float targetKmerMatchProb, float toleratedDeviation){

    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif

    size_t targetDbSize = dbr->getSize();
    if (targetDbSize > 100000)
        targetDbSize = 100000;
    IndexTable* indexTable = getIndexTable(subMat->alphabetSize, kmerSize, dbr, seqs[0], targetDbSize); 

    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

    ExtendedSubstitutionMatrix* _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    ExtendedSubstitutionMatrix* _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

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
    float kmerMatchProb;

    float kmerMatchProbMax = targetKmerMatchProb + (toleratedDeviation * targetKmerMatchProb);
    float kmerMatchProbMin = targetKmerMatchProb - (toleratedDeviation * targetKmerMatchProb);
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
            matchers[thread_idx] = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, dbr->getSeqLens(), kmerThrMid, kmerSize, dbr->getSize(), subMat->alphabetSize);
        }

#pragma omp parallel for schedule(static, 10) reduction (+: dbMatchesSum, querySeqLenSum)
        for (int i = 0; i < querySetSize; i++){
            int id = querySeqs[i];

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char* seqData = dbr->getData(id);
            seqs[thread_idx]->id = id;
            seqs[thread_idx]->mapSequence(seqData);

            matchers[thread_idx]->matchQuery(seqs[thread_idx]);

            dbMatchesSum += seqs[thread_idx]->stats->dbMatches;
            querySeqLenSum += seqs[thread_idx]->L;
        }

        kmerMatchProb = ((float)dbMatchesSum / (float) querySeqLenSum)/targetSeqLenSum; 
        std::cout << "\tMatch probability: " << kmerMatchProb << "\n"; 

        if (kmerMatchProb < kmerMatchProbMin)
            kmerThrMax = kmerThrMid - 1;
        else if (kmerMatchProb > kmerMatchProbMax)
            kmerThrMin = kmerThrMid + 1;
        else if (kmerMatchProb > kmerMatchProbMin && kmerMatchProb < kmerMatchProbMax){
            for (int j = 0; j < threads; j++){
                delete matchers[j];
            }
            delete[] querySeqs;
            delete[] matchers;
            delete indexTable;
            delete _2merSubMatrix;
            delete _3merSubMatrix;
            return kmerThrMid;
        }
    }
    return 0.0;
}


int main (int argc, const char * argv[])
{
    kClust2_cuticle_init();

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int kmerSize =  6;
    int alphabetSize = 21;
    size_t maxSeqLen = 50000;
    size_t BUFFER_SIZE = 10000000;

    std::string queryDB = "";
    std::string queryDBIndex = "";
    std::string targetDB = "";
    std::string targetDBIndex = "";
    std::string outDB = "";
    std::string outDBIndex = "";
    std::string scoringMatrixFile = "";

    parseArgs(argc, argv, &queryDB, &targetDB, &scoringMatrixFile, &outDB, &kmerSize, &alphabetSize, &maxSeqLen);

    std::cout << "Query database: " << queryDB << "\n";
    std::cout << "Target database: " << targetDB << "\n";
    std::cout << "k-mer size: " << kmerSize << "\n";
    std::cout << "Alphabet size: " << alphabetSize << "\n";

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
    qdbr->open();

    DBReader* tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open();

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    // init the substitution matrices
    BaseMatrix* subMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize);

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
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);
        reslens[thread_idx] = new std::list<int>();
    }

    int queryDBSize = qdbr->getSize();
    char** outBuffers = new char*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

    // set the k-mer similarity threshold
    std::cout << "\nAdjusting k-mer similarity threshold...\n";
    short kmerThr = setKmerThreshold (tdbr, seqs, subMat, kmerSize, 1.5e-06, 0.1);
    if (kmerThr == 0.0)
        kmerThr = setKmerThreshold (tdbr, seqs, subMat, kmerSize, 1.5e-06, 0.15);
    if (kmerThr == 0.0){
        std::cout << "ERROR: Could not adjust the k-mer list length threshold!\n";
        std::cout << "Please report this error to the developers Maria Hauser mhauser@genzentrum.lmu.de and Martin Steinegger Martin.Steinegger@campus.lmu.de\n";
        exit(1);
    }
    std::cout << "k-mer similarity threshold: " << kmerThr << "\n\n"; 

    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);
    IndexTable* indexTable = getIndexTable(alphabetSize, kmerSize, tdbr, seq, tdbr->getSize());
    delete seq;
 
    std::cout << "Initializing data structures...";
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif 
        matchers[thread_idx] = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, tdbr->getSeqLens(), kmerThr, kmerSize, tdbr->getSize(), subMat->alphabetSize);
    }
    std::cout << "done!\n";

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    // calculate prefiltering scores for each sequence in the database
    std::cout << "Calculating prefiltering scores!\n";
    long double kmersPerPos = 0.0;
    size_t dbMatches = 0;
    size_t resSize = 0;
    int empty = 0;

#pragma omp parallel for schedule(static, 10) reduction (+: kmersPerPos, resSize, empty, dbMatches)

    for (int id = 0; id < queryDBSize; id++){

        if (id % 1000000 == 0 && id > 0){
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
            size_t prefPassedPerSeq = resSize/id;
            std::cout << "\t" << prefPassedPerSeq << " sequences passed prefiltering per query sequence.\n";
        }

        //std::cout << qdbr->getDbKey(id) << "\n";
        std::list<hit_t>* prefResults;
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        char* seqData = qdbr->getData(id);
        seqs[thread_idx]->id = id;
        seqs[thread_idx]->mapSequence(seqData);

        prefResults = matchers[thread_idx]->matchQuery(seqs[thread_idx]);

        if (prefResults->size() == 0)
            empty++;
        
        if (prefResults->size() > 100000){
            std::cerr << "The prefiltering list for the query " << qdbr->getDbKey(id) << " is too long, the length of the list = " << prefResults->size() << ", maximum allowed is 100 000.\n";
        }

        std::stringstream prefResultsOut;
        int l = 0;
        for (std::list<hit_t>::iterator iter = prefResults->begin(); iter != prefResults->end(); iter++){
            prefResultsOut << tdbr->getDbKey(iter->seqId) << "\t" << iter->prefScore << "\t" << iter->eval << "\n";
            l++;
            // maximum allowed result list length is 100000
            if (l == 100000)
                break;
        }

        std::string prefResultsOutString = prefResultsOut.str();
        const char* prefResultsOutData = prefResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(prefResultsOutData)){
            std::cerr << "Tried to process the prefiltering list for the query " << qdbr->getDbKey(id) << " , the length of the list = " << prefResults->size() << "\n";
            std::cerr << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffers[thread_idx], prefResultsOutData, prefResultsOutString.length()*sizeof(char));
        dbw->write(outBuffers[thread_idx], prefResultsOutString.length(), qdbr->getDbKey(id), thread_idx);

        kmersPerPos += seqs[thread_idx]->stats->kmersPerPos;
        dbMatches += seqs[thread_idx]->stats->dbMatches;
        resSize += prefResults->size();
        reslens[thread_idx]->push_back(prefResults->size());
    }

    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }

    //    matcher->printStats();
    kmersPerPos /= queryDBSize;
    size_t dbMatchesPerSeq = dbMatches/(size_t)queryDBSize;
    size_t prefPassedPerSeq = resSize/(size_t)queryDBSize;
    std::cout << kmersPerPos << " k-mers per position.\n";
    std::cout << dbMatchesPerSeq << " DB matches per sequence.\n";
    std::cout << prefPassedPerSeq << " sequences passed prefiltering per query sequence.\n";

    int mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
    std::advance(it, mid);
    std::cout << "Median result list size: " << *it << "\n";
/*    std::ofstream myfile;
    myfile.open ("/net/cluster/user/maria/test/ressizes.dat");
    for (it = reslens[0]->begin(); it != reslens[0]->end(); it++)
        myfile << *it << "\n";
    myfile.close();
*/
    std::cout << empty << " sequences with zero size result lists.\n";

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "Time for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

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
