#include <iostream>
#include <time.h>
#include <unistd.h>
#include <string>

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

void printUsage(){

    std::string usage("\nCalculates similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de) & Martin Steinegger (Martin.Steinegger@campus.lmu.de)\n\n");
    usage.append("USAGE: kClust2_pref ffindexQueryDBBase ffindexTargetDBBase scoringMatrixFile ffindexOutDBBase [opts]\n"
            "-t\t[float]\tPrefiltering threshold (minimum half bits per query position).\n"
            "-k\t[int]\tk-mer size (default=6).\n"
            "-a\t[int]\tAmino acid alphabet size (default=21).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* scoringMatrixFile, std::string* ffindexOutDBBase, int* kmerSize, int* alphabetSize){
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


IndexTable* getIndexTable (int alphabetSize, int kmerSize, DBReader* dbr, Sequence* seq){

    int dbSize = dbr->getSize();

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

    std::cout << "Index table: removing duplicate entries...\n";
    indexTable->removeDuplicateEntries();

    return indexTable;
}

short getKmerThreshold (DBReader* dbr, IndexTable* indexTable, Sequence** seqs, BaseMatrix* subMat, int kmerSize, float targetKmerMatchProb){

    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif

    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

    ExtendedSubstitutionMatrix* _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    ExtendedSubstitutionMatrix* _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

    int seqLenSum = 0;
    for (size_t i = 0; i < dbr->getSize(); i++)
        seqLenSum += dbr->seqLens[i];

    // generate a small random sequence set for testing 
    int* testSeqs = new int[1000];
    for (int i = 0; i < 1000; i++){
        testSeqs[i] = rand() % dbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability 
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb
    short kmerThrMin = 20;
    short kmerThrMax = 160;
    short kmerThrMid;

    int dbMatchesSum;
    int querySeqLens;
    float kmerMatchProb;

    float kmerMatchProbMax = targetKmerMatchProb + (0.1 * targetKmerMatchProb);
    float kmerMatchProbMin = targetKmerMatchProb - (0.1 * targetKmerMatchProb);
    std::cout << "Searching for a k-mer threshold with a k-mer match probability per query position within range [" << kmerMatchProbMin << ":" << kmerMatchProbMax << "].\n";
    
    // adjust k-mer list length threshold
    while (kmerThrMax >= kmerThrMin){
        dbMatchesSum = 0;
        querySeqLens = 0;

        kmerThrMid = kmerThrMin + (kmerThrMax - kmerThrMin)/2;

        // determine k-mer match probability for kmerThrMid
#pragma omp parallel for schedule(static) 
        for (int i = 0; i < threads; i++){
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // set a current k-mer list length threshold and a high prefitlering threshold (we don't need the prefiltering results in this test run)
            matchers[thread_idx] = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, dbr->seqLens, kmerThrMid, kmerSize, dbr->getSize(), subMat->alphabetSize);
        }

#pragma omp parallel for schedule(static, 10) reduction (+: dbMatchesSum, querySeqLens)
        for (int i = 0; i < 1000; i++){
            int id = testSeqs[i];
        
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char* seqData = dbr->getData(id);
            seqs[thread_idx]->id = id;
            seqs[thread_idx]->mapSequence(seqData);

            matchers[thread_idx]->matchQuery(seqs[thread_idx]);
            
            dbMatchesSum += seqs[thread_idx]->stats->dbMatches;
            querySeqLens += seqs[thread_idx]->L;
        }

        kmerMatchProb = ((float)dbMatchesSum / (float) querySeqLens)/seqLenSum; 
        std::cout << "k-mer threshold range: [" << kmerThrMin  << ":" << kmerThrMax << "]\n";
        std::cout << "\tTried threshold " << kmerThrMid << ", match probability: " << kmerMatchProb << "\n"; 
        
        if (kmerMatchProb < kmerMatchProbMin)
            kmerThrMax = kmerThrMid - 1;
        else if (kmerMatchProb > kmerMatchProbMax)
            kmerThrMin = kmerThrMid + 1;
        else if (kmerMatchProb > kmerMatchProbMin && kmerMatchProb < kmerMatchProbMax){
            delete[] testSeqs;
            return kmerThrMid;
        }
    }
    std::cout << "ERROR: Could not adjust the k-mer list length threshold!\n";
    std::cout << "Aborted at k-mer threshold " << kmerThrMid << " and match probability " << kmerMatchProb << "\n";
    std::cout << "Tolerated range of match probability would be [" << kmerMatchProbMin  << ":" << kmerMatchProbMax << "]\n";
    std::cout << "Please report this error to the developers Maria Hauser mhauser@genzentrum.lmu.de and Martin Steinegger Martin.Steinegger@campus.lmu.de\n";
    exit(1);
}



int main (int argc, const char * argv[])
{
    int kmerSize =  6;
    int alphabetSize = 21;
    size_t maxSeqLen = 40000;
    size_t BUFFER_SIZE = 10000000;

    clock_t c = clock();
    std::string queryDB = "";
    std::string queryDBIndex = "";
    std::string targetDB = "";
    std::string targetDBIndex = "";
    std::string outDB = "";
    std::string outDBIndex = "";
    std::string scoringMatrixFile = "";

    parseArgs(argc, argv, &queryDB, &targetDB, &scoringMatrixFile, &outDB, &kmerSize, &alphabetSize);
    
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

    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);

    IndexTable* indexTable = getIndexTable(alphabetSize, kmerSize, tdbr, seq);

    delete seq;
    
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
    short kmerThr = getKmerThreshold (tdbr, indexTable, seqs, subMat, kmerSize, 1.5e-06);
    std::cout << "k-mer similarity threshold: " << kmerThr << "\n"; 

    std::cout << "Initializing data structures...";
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif 
        matchers[thread_idx] = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, tdbr->seqLens, kmerThr, kmerSize, tdbr->getSize(), subMat->alphabetSize);
    }
    std::cout << "done!";

    c = clock() -c ;
    int sec = (int)((float)c/CLOCKS_PER_SEC);
    std::cout << "Time for init: " << sec/60 << " m " << (sec % 60) << "s\n";
    c = clock();

    // calculate prefiltering scores for each sequence in the database
    std::cout << "Calculating prefiltering scores!\n";
    long double kmersPerPos = 0.0;
    size_t dbMatches = 0;
    size_t resSize = 0;
    int empty = 0;

#pragma omp parallel for schedule(static, 10) reduction (+: kmersPerPos, resSize, empty, dbMatches)

    for (int id = 0; id < queryDBSize; id++){

        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
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

        std::stringstream prefResultsOut;
        for (std::list<hit_t>::iterator iter = prefResults->begin(); iter != prefResults->end(); iter++){
            //            std::cout << tdbr->getDbKey(iter->seqId);
            //            std::cout << "\tscore: " << iter->prefScore << "\n";
            prefResultsOut << tdbr->getDbKey(iter->seqId) << "\t" << iter->prefScore << "\t" << iter->eval << "\n";
        }

        std::string prefResultsOutString = prefResultsOut.str();
        const char* prefResultsOutData = prefResultsOutString.c_str();
        if (BUFFER_SIZE < prefResultsOutString.length()){
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
    
    int mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
    std::advance(it, mid);
    std::cout << "Median result list size: " << *it << "\n";
    std::ofstream myfile;
    myfile.open ("~/test/ressizes.dat");
    for (it = reslens[0]->begin(); it != reslens[0]->end(); it++)
        myfile << *it << "\n";
    myfile.close();

    //    matcher->printStats();
    kmersPerPos /= queryDBSize;
    size_t dbMatchesPerSeq = dbMatches/(size_t)queryDBSize;
    long double prefPassedPerSeq = (long double)resSize/(long double)queryDBSize;
    std::cout << kmersPerPos << " k-mers per position.\n";
    std::cout << dbMatchesPerSeq << " DB matches per sequence.\n";
    std::cout << prefPassedPerSeq << " sequences passed prefiltering per query sequence.\n";
    std::cout << empty << " sequences with zero size result lists.\n";

    c = clock() - c;
    sec = (int)((float)c/CLOCKS_PER_SEC);
    std::cout << "Time for the scores calculation: " << sec/60 << " m " << sec%60 << "s\n";

    std::cout << "Merging the results...\n";
    qdbr->close();
    if (queryDB.compare(targetDB) != 0)
        tdbr->close();
    dbw->close();
    std::cout << "done.\n";

    return 0;
}
