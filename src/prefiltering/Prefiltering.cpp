#include "Prefiltering.h"

Prefiltering::Prefiltering(std::string queryDB,
        std::string targetDB,
        std::string outDB,
        std::string scoringMatrixFile,
        float sensitivity,
        int kmerSize,
        int alphabetSize,
        float zscoreThr,
        size_t maxSeqLen,
        int seqType,
        bool aaBiasCorrection,
        int skip){


    BUFFER_SIZE = 1000000;

    std::string queryDBIndex = queryDB + ".index";
    std::string targetDBIndex = targetDB + ".index";
    std::string outDBIndex = outDB + ".index";

    this->kmerSize = kmerSize;
    this->alphabetSize = alphabetSize;

    threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    std::cout << "Using " << threads << " threads.\n";
#endif
    std::cout << "\n";

    std::cout << "Init data structures...\n";
    this->qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader::NOSORT);

    this->tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open(DBReader::SORT);

    this->dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    std::cout << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";
    std::cout << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

    // init the substitution matrices
    if (seqType == Sequence::AMINO_ACIDS)
        subMat = getSubstitutionMatrix(scoringMatrixFile, 8.0);
    else
        subMat = new NucleotideMatrix();

    _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

    // init all thread-specific data structures 
    this->seqs = new Sequence*[threads];
    this->reslens = new std::list<int>*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType);
        reslens[thread_idx] = new std::list<int>();
    }

    outBuffers = new char*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

    // set the k-mer similarity threshold
    float p_match =  pow(2.0, sensitivity) * 1.0e-08 * (float) (skip + 1); // old target value: 1.5e-06, reached with sens = 7.2 approximately
    std::cout << "\nAdjusting k-mer similarity threshold within +-10% deviation from the target k-mer match probability (target probability = " << p_match << ")...\n";
    std::pair<short, double> ret = setKmerThreshold (qdbr, maxSeqLen, p_match, 0.1);
    short kmerThr = ret.first; //174; //ret.first; // 103;
    double kmerMatchProb = ret.second; //1.16383e-09; //ret.second; // 1.57506e-06;
    if (kmerThr == 0.0){
        std::cout << "Could not set the probability within +-10% deviation. Trying +-15% deviation.\n";
        ret = setKmerThreshold (qdbr, maxSeqLen, p_match, 0.15);
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
    this->indexTable = getIndexTable(tdbr, seq, tdbr->getSize(), skip);
    delete seq;

    std::cout << "Initializing data structures...";
    matchers = new QueryTemplateMatcher*[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix, indexTable, tdbr->getSeqLens(), kmerThr, kmerMatchProb, kmerSize, tdbr->getSize(), aaBiasCorrection, maxSeqLen, zscoreThr);
    }
    std::cout << "... done.\n";


}

Prefiltering::~Prefiltering(){
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
}

void Prefiltering::run(size_t maxResListLen){

    double kmersPerPos = 0.0;
    int dbMatches = 0;

    int empty = 0;
    int resSize = 0;

    int queryDBSize = qdbr->getSize();
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
        prefResults = matchers[thread_idx]->matchQuery(seqs[thread_idx]);

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

    // merge output ffindex databases
    std::cout << "Merging the results...\n";
    qdbr->close();
    if (strcmp(qdbr->getIndexFileName(), tdbr->getIndexFileName()) != 0)
        tdbr->close();
    dbw->close();
    std::cout << "done.\n";
}


BaseMatrix* Prefiltering::getSubstitutionMatrix(std::string scoringMatrixFile, float bitFactor){
    std::cout << "Substitution matrices...\n";
    BaseMatrix* subMat;
    if (alphabetSize < 21){
        SubstitutionMatrix* sMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), bitFactor);
        subMat = new ReducedMatrix(sMat->probMatrix, alphabetSize);
    }
    else
        subMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), bitFactor);

    return subMat;
}


IndexTable* Prefiltering::getIndexTable (DBReader* dbr, Sequence* seq, int dbSize, int skip){

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

    std::cout << "Index table: init...\n";
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

    std::cout << "Index table: removing duplicate entries...\n";
    indexTable->removeDuplicateEntries();

    return indexTable;
}

std::pair<short,double> Prefiltering::setKmerThreshold (DBReader* dbr, int maxSeqLen, double targetKmerMatchProb, double toleratedDeviation){

    size_t targetDbSize = dbr->getSize();
    if (targetDbSize > 100000)
        targetDbSize = 100000;
    IndexTable* indexTable = getIndexTable(dbr, seqs[0], targetDbSize);

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
