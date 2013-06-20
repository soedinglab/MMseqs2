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
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de) & Marting Steinegger (Martin.Steinegger@campus.lmu.de)\n\n");
    usage.append("USAGE: kClust2_pref ffindexQueryDBBase ffindexTargetDBBase scoringMatrixFile ffindexOutDBBase [opts]\n"
            "-t\t[float]\tPrefiltering threshold (minimum half bits per query position).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* scoringMatrixFile, std::string* ffindexOutDBBase, float* prefThr){
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
        if (strcmp(argv[i], "-t") == 0){
            if (++i < argc){
                *prefThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -t\n";
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

int main (int argc, const char * argv[])
{
    float prefThr = 0.55f;
    int kmerSize =  6;
    short kmerThr = 27;
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
    
    parseArgs(argc, argv, &queryDB, &targetDB, &scoringMatrixFile, &outDB, &prefThr);

    std::cout << "prefThr: " << prefThr << "\n";

    queryDBIndex = queryDB + ".index";
    targetDBIndex = targetDB + ".index";
    outDBIndex = outDB + ".index";

     int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    std::cout << "Using " << threads << " threads.\n";
#endif

   std::cout << "Init data structures:\n";
    DBReader* qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open();

    DBReader* tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open();

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    // init the substitution matrices

    SubstitutionMatrix* subMat = new SubstitutionMatrix (scoringMatrixFile.c_str());

//    ReducedMatrix* subMat = new ReducedMatrix(sMat->probMatrix, 13);
//    BaseMatrix::print(subMat->subMatrix, subMat->int2aa, subMat->alphabetSize);
    
    int targetDBSize = tdbr->getSize();
    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);
    
    std::cout << "Index table: counting k-mers...\n";
    // fill and init the index table
    IndexTable* indexTable = new IndexTable(subMat->alphabetSize, kmerSize);
    for (int id = 0; id < targetDBSize; id++){
        if (id % 1000000 == 0)
            std::cout << id << "\n";
        seq->id = id;
        char* seqData = tdbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(seqData);
        indexTable->addKmerCount(seq);
    }
    std::cout << "Index table init...\n";
    indexTable->init();

    std::cout << "Index table: fill...\n";
    for (int id = 0; id < targetDBSize; id++){

        if (id % 1000000 == 0)
            std::cout << id << "\n";

        seq->id = id;
        char* seqData = tdbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(seqData);
        indexTable->addSequence(seq);
    }

    delete subMat;
    delete seq;


    std::cout << "Index table: removing duplicate entries...\n";
    indexTable->removeDuplicateEntries();

    std::cout << "Substitution matrices...\n";
    
    subMat = new SubstitutionMatrix (scoringMatrixFile.c_str());
    ExtendedSubstitutionMatrix* _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    ExtendedSubstitutionMatrix* _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];
    Sequence** seqs = new Sequence*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);
        matchers[thread_idx] = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, tdbr->seqLens, kmerThr, prefThr, kmerSize, targetDBSize, subMat->alphabetSize); 
    }

    c = clock() -c ;
    int sec = (int)((float)c/CLOCKS_PER_SEC);
    std::cout << "Time for init: " << sec/60 << " m " << (sec % 60) << "s\n";
    c = clock();

    // calculate prefiltering scores for each sequence in the database
    std::cout << "Calculating prefiltering scores!\n";
    int queryDBSize = qdbr->getSize();
    if (queryDBSize > 100)
        queryDBSize = 100;
    char** outBuffers = new char*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

    double kmersPerPos = 0.0;
    size_t dbMatches = 0;
    int resSize = 0;
    int empty = 0;

#pragma omp parallel for schedule(static, 10) reduction (+: kmersPerPos, dbMatches, resSize, empty)
    for (int id = 0; id < queryDBSize; id++){
        std::list<hit_t>* prefResults;
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        char* seqData = qdbr->getData(id);
        seqs[thread_idx]->id = id;
        seqs[thread_idx]->mapSequence(seqData);
//        if (id > 0 && id % 10 == 0)
//            std::cout << id << "\n";
    //        std::cout << "Sequence: " << qdbr->getDbKey(id) << " (length: " << seqs[thread_idx]->L << ")\n";
//        std::cout << seqData << "\n";
//        seqs[thread_idx]->print();

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
            std::cerr << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            exit(1);
        }
        memcpy(outBuffers[thread_idx], prefResultsOutData, prefResultsOutString.length()*sizeof(char));
        dbw->write(outBuffers[thread_idx], prefResultsOutString.length(), qdbr->getDbKey(id), thread_idx);

        kmersPerPos += seqs[thread_idx]->stats->kmersPerPos;
        dbMatches += seqs[thread_idx]->stats->dbMatches;
        resSize += prefResults->size();
    }
//    matcher->printStats();
    kmersPerPos /= queryDBSize;
    size_t dbMatchesPerSeq = dbMatches/(size_t)queryDBSize;
    double prefPassedPerSeq = (double)resSize/(double)queryDBSize;
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
