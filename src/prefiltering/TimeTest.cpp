#include "TimeTest.h"

TimeTest::TimeTest(std::string queryDB,
        std::string queryDBIndex,
        std::string targetDB,
        std::string targetDBIndex,
        std::string scoringMatrixFile,
        int alphabetSize,
        size_t maxSeqLen,
        std::string logFile){

    std::cout << "Init data structures...\n";
    BUFFER_SIZE = 1000000;

    this->logFile = logFile;

    this->maxSeqLen = maxSeqLen;

    this->alphabetSize = alphabetSize;

    this->threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    std::cout << "Using " << threads << " threads.\n";
#endif
    std::cout << "\n";

    this->qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader::NOSORT);

    this->tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open(DBReader::SORT);

    std::cout << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";
    std::cout << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

    subMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), 8.0);

    _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

    this->seqs = new Sequence*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS);
    }

    std::cout << "Init done.\n\n";
}

TimeTest::~TimeTest(){

    for (int i = 0; i < threads; i++){
        delete seqs[i];
    }
    delete[] seqs;

    delete subMat;
    delete _2merSubMatrix;
    delete _3merSubMatrix;
}

void TimeTest::runTimeTest (){

    std::ofstream logFileStream;
    logFileStream.open(logFile.c_str());

    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

    int targetSeqLenSum = 0;
    for (size_t i = 0; i < tdbr->getSize(); i++)
        targetSeqLenSum += tdbr->getSeqLens()[i];

    // generate a small random sequence set for testing 
    int querySetSize = tdbr->getSize();
    if (querySetSize > 1000)
        querySetSize = 1000;

    int* querySeqs = new int[querySetSize];
    srand(1);
    size_t querySeqLenSum = 0;
    for (int i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % tdbr->getSize();
        querySeqLenSum += tdbr->getSeqLens()[querySeqs[i]];
    }

    short kmerThrPerPosMin = 10;
    short kmerThrPerPosMax = 20;

    // adjust k-mer list length threshold
    for (int kmerSize = 4; kmerSize <= 7; kmerSize++){
        short kmerThrMin = kmerThrPerPosMin * kmerSize;
        int kmerThrMax = kmerThrPerPosMax * kmerSize;

        std::cout << "------------------  k = " << kmerSize << " -----------------------------\n";
        IndexTable* indexTable = Prefiltering::getIndexTable(tdbr, seqs[0], alphabetSize, kmerSize, 0, tdbr->getSize());
        for (short kmerThr = kmerThrMin; kmerThr < kmerThrMax; kmerThr += kmerSize){
            std::cout << "k = " << kmerSize << "\n";
            std::cout << "k-mer threshold = " << kmerThr << "\n";

            size_t dbMatchesSum = 0;

            double kmersPerPos = 0.0;
            double kmerMatchProb = 0.0;

            // determine k-mer match probability and k-mer list length for kmerThr
#pragma omp parallel for schedule(static) 
            for (int i = 0; i < threads; i++){
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                // set a current k-mer list length threshold and a high prefitlering threshold (we don't need the prefiltering results in this test run)
                matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix, indexTable, tdbr->getSeqLens(), kmerThr, 1.0, kmerSize, tdbr->getSize(), false, maxSeqLen, 500.0);
            }

            struct timeval start, end;
            gettimeofday(&start, NULL);

#pragma omp parallel for schedule(dynamic, 10) reduction (+: dbMatchesSum, kmersPerPos)
            for (int i = 0; i < querySetSize; i++){
                int id = querySeqs[i];

                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                char* seqData = tdbr->getData(id);
                seqs[thread_idx]->mapSequence(id, tdbr->getDbKey(id), seqData);

                matchers[thread_idx]->matchQuery(seqs[thread_idx]);

                kmersPerPos += seqs[thread_idx]->stats->kmersPerPos;
                dbMatchesSum += seqs[thread_idx]->stats->dbMatches;
            }

            kmersPerPos /= (double)querySetSize;
            kmerMatchProb = ((double)dbMatchesSum) / ((double) (querySeqLenSum * targetSeqLenSum));

            std::cout << "kmerPerPos: " << kmersPerPos << "\n";
            std::cout << "k-mer match probability: " << kmerMatchProb << "\n";

            gettimeofday(&end, NULL);
            int sec = end.tv_sec - start.tv_sec;
            std::cout << "Time: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";

            logFileStream << kmersPerPos << "\t" << kmerMatchProb << "\t" << sec << "\n";

            for (int j = 0; j < threads; j++){
                delete matchers[j];
            }
        }
        delete indexTable;
    }

    logFileStream.close();
    delete[] querySeqs;
    delete[] matchers;
}
 
