#include "Prefiltering.h"
#include "../commons/Util.h"

Prefiltering::Prefiltering(std::string queryDB,
        std::string queryDBIndex,
        std::string targetDB,
        std::string targetDBIndex,
        std::string outDB,
        std::string outDBIndex,
        std::string scoringMatrixFile,
        float sensitivity,
        int kmerSize,
        int alphabetSize,
        float zscoreThr,
        size_t maxSeqLen,
        int seqType,
        bool aaBiasCorrection,
        int splitSize,
        int skip):    outDB(outDB),
    outDBIndex(outDBIndex),
    kmerSize(kmerSize),
    alphabetSize(alphabetSize),
    zscoreThr(zscoreThr),
    maxSeqLen(maxSeqLen),
    seqType(seqType),
    aaBiasCorrection(aaBiasCorrection),
    splitSize(splitSize),
    skip(skip)
{

    this->threads = 1;
#ifdef OPENMP
    this->threads = omp_get_max_threads();
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif
    Debug(Debug::INFO) << "\n";

    this->qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader::NOSORT);

    this->tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open(DBReader::SORT);

    if (this->splitSize == 0)
        this->splitSize = tdbr->getSize();

    std::string outDBTmp = outDB + "_tmp";
    std::string outDBIndexTmp = outDBIndex.c_str()+std::string("_tmp");

    tmpDbw = new DBWriter(outDBTmp.c_str(), outDBIndexTmp.c_str(), threads);
    tmpDbw->open();

    Debug(Debug::INFO) << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";
    Debug(Debug::INFO) << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

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
    Debug(Debug::INFO) << "\nAdjusting k-mer similarity threshold within +-10% deviation from the reference time value, sensitivity = " << sensitivity << ")...\n";
    std::pair<short, double> ret = setKmerThreshold (tdbr, sensitivity, 0.1);
    this->kmerThr = ret.first;
    this->kmerMatchProb = ret.second;

    Debug(Debug::WARNING) << "k-mer similarity threshold: " << kmerThr << "\n";
    Debug(Debug::WARNING) << "k-mer match probability: " << kmerMatchProb << "\n\n";

    // initialise the index table and the matcher structures for the database
    // Init for next split
    this->matchers = new QueryTemplateMatcher*[threads];
}

Prefiltering::~Prefiltering(){
    for (int i = 0; i < threads; i++){
        delete seqs[i];
        delete[] outBuffers[i];
        delete reslens[i];
    }
    delete[] seqs;
    delete[] outBuffers;
    delete[] reslens;

    delete subMat;
    delete _2merSubMatrix;
    delete _3merSubMatrix;
}

void Prefiltering::run(size_t maxResListLen){

    size_t kmersPerPos = 0;
    size_t dbMatches = 0;

    size_t resSize = 0;
    size_t realResSize = 0;

    size_t queryDBSize = qdbr->getSize();
    int splitCount = 0;
    int* notEmpty = new int[queryDBSize];
    memset(notEmpty, 0, queryDBSize*sizeof(int));

    // splits template database into chunks
    int step = 0;
    for(unsigned int splitStart = 0; splitStart < tdbr->getSize(); splitStart += splitSize ){
        splitCount++;
        std::string idSuffix;
        std::stringstream idSuffixStream;
        if (splitStart==0) {
            idSuffixStream << "";
        }else {
            idSuffixStream << "." << splitCount;
        }
        idSuffix = idSuffixStream.str();


        Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType);
        this->indexTable = getIndexTable(tdbr, seq, alphabetSize, kmerSize, splitStart, splitStart + splitSize , skip);
        delete seq;
        int stepCnt = (tdbr->getSize() + splitSize - 1) / splitSize;
        Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step " << ++step << " of " << stepCnt <<  ")\n";

        struct timeval start, end;
        gettimeofday(&start, NULL);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < this->threads; i++){
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix,
                    indexTable, tdbr->getSeqLens(), kmerThr,
                    kmerMatchProb, kmerSize, tdbr->getSize(),
                    aaBiasCorrection, maxSeqLen, zscoreThr);
        }

#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, realResSize, dbMatches)
        for (size_t id = 0; id < queryDBSize; id++){

            Log::printProgress(id);

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // get query sequence
            char* seqData = qdbr->getData(id);
            seqs[thread_idx]->mapSequence(id, qdbr->getDbKey(id), seqData);

            // calculate prefitlering results
            std::pair<hit_t *, size_t> prefResults = matchers[thread_idx]->matchQuery(seqs[thread_idx], tdbr->getId(seqs[thread_idx]->getDbKey()));

            const size_t resultSize = prefResults.second;
            if(writePrefilterOutput(thread_idx, idSuffix, id, maxResListLen, prefResults) != 0)
                continue; // couldnt write result - too many results

            // update statistics counters
            if (resultSize != 0)
                notEmpty[id] = 1;
            kmersPerPos += (size_t) seqs[thread_idx]->stats->kmersPerPos;
            dbMatches += seqs[thread_idx]->stats->dbMatches;
            resSize += resultSize;
            realResSize += std::min(resultSize, maxResListLen);
            reslens[thread_idx]->push_back(resultSize);
        } // step end
        if (queryDBSize > 1000)
            Debug(Debug::INFO) << "\n";
        Debug(Debug::INFO) << "\n";

        for (int j = 0; j < threads; j++){
            delete matchers[j];
        }
        delete indexTable;

        gettimeofday(&end, NULL);
        int sec = end.tv_sec - start.tv_sec;
        Debug(Debug::WARNING) << "\nTime for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    } // prefiltering scores calculation one split end
    int empty = 0;
    for (unsigned int i = 0; i < qdbr->getSize(); i++){
        if (notEmpty[i] == 0){
            empty++;
        }
    }

    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }
    // correction because of x splits
    kmersPerPos = kmersPerPos / splitCount;
    // print statistics
    this->printStatistics(queryDBSize, kmersPerPos, resSize, realResSize, dbMatches, empty, maxResListLen, reslens[0]);

    // close reader to reduce memory
    qdbr->close();
    if (strcmp(qdbr->getIndexFileName(), tdbr->getIndexFileName()) != 0)
        tdbr->close();

    // merge output ffindex databases
    Debug(Debug::INFO) << "Merging the results...\n";
    tmpDbw->close(); // sorts the index
   
    std::cout << "Threads merged...\n";
    DBReader tmpReader(tmpDbw->getDataFileName(), tmpDbw->getIndexFileName());
    tmpReader.open(DBReader::SORT);

    this->dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1);
    dbw->open();

    for (size_t id = 0; id < queryDBSize; id++){
        std::stringstream mergeResultsOut;
        for(int split = 0; split < splitCount; split++)
            mergeResultsOut << tmpReader.getData(id+(split*queryDBSize));

        std::string mergeResultsOutString = mergeResultsOut.str();
        const char* mergeResultsOutData = mergeResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(mergeResultsOutData)){
            Debug(Debug::ERROR) << "ERROR: Buffer overflow during the merging.\n";
            exit(EXIT_FAILURE);
        }
        memcpy(outBuffers[0], mergeResultsOutData, mergeResultsOutString.length()*sizeof(char));
        dbw->write(outBuffers[0], mergeResultsOutString.length(),  tmpReader.getDbKey(id), 0);

    }
    std::cout << "Closing tmp reader...\n";
    tmpReader.close();
    remove(tmpDbw->getDataFileName());
    remove(tmpDbw->getIndexFileName());
    std::cout << "Closing dbw...\n";
    dbw->close();
    std::cout << "done.\n";

}

// write prefiltering to ffindex database
int Prefiltering::writePrefilterOutput( int thread_idx, std::string idSuffix, size_t id,
        size_t maxResListLen, std::pair<hit_t *, size_t> prefResults){
    // write prefiltering results to a string
    std::stringstream prefResultsOut;
    size_t l = 0;
    hit_t * resultVector = prefResults.first;
    const size_t resultSize = prefResults.second;

    for (size_t i = 0; i < resultSize; i++){
        hit_t * res = resultVector + i;

        if (res->seqId >= tdbr->getSize()) {
            Debug(Debug::ERROR) << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id)<< " -> " << res->seqId << "\t" << res->prefScore << "\n";
            continue;
        }
        prefResultsOut << tdbr->getDbKey(res->seqId) << "\t" << res->zScore << "\t" << res->prefScore << "\n";
        l++;
        // maximum allowed result list length is reached
        if (l == maxResListLen)
            break;
    }
    // write prefiltering results string to ffindex database
    std::string prefResultsOutString = prefResultsOut.str();
    const char* prefResultsOutData = prefResultsOutString.c_str();
    if (BUFFER_SIZE < strlen(prefResultsOutData)){
        Debug(Debug::ERROR) << "Tried to process the prefiltering list for the query " << qdbr->getDbKey(id) << " , the length of the list = " << resultSize << "\n";
        Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
        return -1;
    }
    memcpy(outBuffers[thread_idx], prefResultsOutData, prefResultsOutString.length()*sizeof(char));
    std::stringstream keyStream;
    keyStream << qdbr->getDbKey(id) << idSuffix;
    tmpDbw->write(outBuffers[thread_idx], prefResultsOutString.length(), (char *) keyStream.str().c_str() , thread_idx);
    return 0;

}


void Prefiltering::printStatistics(size_t queryDBSize, size_t kmersPerPos,
        size_t resSize,  size_t realResSize,   size_t dbMatches,
        int empty, size_t maxResListLen,
        std::list<int>* reslens){
    size_t dbMatchesPerSeq = dbMatches/queryDBSize;
    size_t prefPassedPerSeq = resSize/queryDBSize;
    size_t prefRealPassedPerSeq = realResSize/queryDBSize;
    Debug(Debug::WARNING) << kmersPerPos/queryDBSize << " k-mers per position.\n";
    Debug(Debug::WARNING) << dbMatchesPerSeq << " DB matches per sequence.\n";
    Debug(Debug::WARNING) << prefPassedPerSeq << " sequences passed prefiltering per query sequence";
    if (prefPassedPerSeq > maxResListLen)
        Debug(Debug::WARNING) << " (ATTENTION: max. " << maxResListLen << " best scoring sequences were written to the output prefiltering database).\n";
    else
        Debug(Debug::WARNING) << ".\n";
    Debug(Debug::WARNING) << prefRealPassedPerSeq << " sequences per query sequence were really written, after restricting maximum list length to " << maxResListLen << "\n";

    int mid = reslens->size() / 2;
    std::list<int>::iterator it = reslens->begin();
    std::advance(it, mid);
    Debug(Debug::INFO) << "Median result list size: " << *it << "\n";
    Debug(Debug::INFO) << empty << " sequences with 0 size result lists.\n";
}

BaseMatrix* Prefiltering::getSubstitutionMatrix(std::string scoringMatrixFile, float bitFactor){
    Debug(Debug::INFO) << "Substitution matrices...\n";
    BaseMatrix* subMat;
    if (alphabetSize < 21){
        SubstitutionMatrix* sMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), bitFactor);
        subMat = new ReducedMatrix(sMat->probMatrix, alphabetSize);
    }
    else
        subMat = new SubstitutionMatrix (scoringMatrixFile.c_str(), bitFactor);

    return subMat;
}


IndexTable* Prefiltering::getIndexTable (DBReader* dbr, Sequence* seq, int alphabetSize,
        int kmerSize, size_t dbFrom, size_t dbTo, int skip){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    IndexTable* indexTable = new IndexTable(alphabetSize, kmerSize, skip);
    dbTo=std::min(dbTo,dbr->getSize());
    for (unsigned int id = dbFrom; id < dbTo; id++){
        Log::printProgress(id-dbFrom);
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(id, dbr->getDbKey(id), seqData);
        indexTable->addKmerCount(seq);
    }

    if ((dbTo-dbFrom) > 10000)
        Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";
    indexTable->init();

    Debug(Debug::INFO) << "Index table: fill...\n";
    for (unsigned int id = dbFrom; id < dbTo; id++){
        Log::printProgress(id-dbFrom);
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(id, dbr->getDbKey(id), seqData);
        indexTable->addSequence(seq);
    }

    if ((dbTo-dbFrom) > 10000)
        Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: removing duplicate entries...\n";
    indexTable->removeDuplicateEntries();
    Debug(Debug::INFO) << "Index table init done.\n\n";

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for index table init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    return indexTable;
}

std::pair<short,double> Prefiltering::setKmerThreshold (DBReader* dbr, double sensitivity, double toleratedDeviation){

    size_t targetDbSize = std::min( dbr->getSize(), (size_t) 100000);
    IndexTable* indexTable = getIndexTable(dbr, seqs[0], alphabetSize, kmerSize, 0, targetDbSize);

    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

    int targetSeqLenSum = 0;
    for (size_t i = 0; i < targetDbSize; i++)
        targetSeqLenSum += dbr->getSeqLens()[i];

    // generate a small random sequence set for testing 
    int querySetSize = std::min ( dbr->getSize(), (size_t)1000);

    int* querySeqs = new int[querySetSize];
    srand(1);
    for (int i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % dbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability 
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb
    short kmerThrMin = 3 * kmerSize;
    short kmerThrMax = 80 * kmerSize;
    short kmerThrMid;

    size_t dbMatchesSum;
    size_t querySeqLenSum;
    size_t dbMatchesExp_pc;
    // 1000 * 350 * 100000 * 350
    size_t lenSum_pc = 12250000000000;

    double kmersPerPos = 0.0;
    double kmerMatchProb;

    // parameters for searching
    // fitted function: Time ~ alpha * kmer_list_len + beta * kmer_match_prob + gamma
    double alpha;
    double beta;
    double gamma;

    // the parameters of the fitted function depend on k
    if (kmerSize == 4){ 
        alpha = 6.974347e-01; // 6.717981e-01;
        beta = 6.954641e+05; // 6.990462e+05;
        gamma = 1.194005; // 1.718601;
    }
    else if (kmerSize == 5){ 
        alpha = 2.133863e-01; // 2.013548e-01;
        beta = 7.612418e+05; // 7.781889e+05;
        gamma = 1.959421; // 1.997792;
    }
    else if (kmerSize == 6){ 
        alpha = 1.141648e-01; // 1.114936e-01
        beta = 9.033168e+05; // 9.331253e+05;
        gamma = 1.411142; // 1.416222;
    }
    else if (kmerSize == 7){ 
        alpha = 6.356035e-02; //7.123599e-02; //6.438574e-02; // 6.530289e-02;
        beta = 3.432066e+06; //3.148479e+06; //3.480680e+06; // 3.243035e+06;
        gamma = 1.238231; //1.304421; // 1.753651; //1.137125;
    }
    else{
        Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
        exit(EXIT_FAILURE);
    }

    // Run using k=6, a=21, with k-mer similarity threshold 103
    // k-mer list length was 117.6, k-mer match probability 1.12735e-06
    // time reference value for these settings is 15.6
    // Since we want to represent the time value in the form base^sensitivity with base=2.0, it yields sensitivity = 3.96 (~4.0)
    double base = 2.0;
    double timevalMax = pow(base, sensitivity) * (1.0 + toleratedDeviation);
    double timevalMin = pow(base, sensitivity) * (1.0 - toleratedDeviation);

    // in case the time value cannot be set within the threshold boundaries, the best value that could be reached will be returned with a warning
    double timevalBest = 0.0;
    short kmerThrBest = 0;
    double kmerMatchProbBest = 0.0;

    // adjust k-mer list length threshold
    while (kmerThrMax >= kmerThrMin){
        dbMatchesSum = 0;
        querySeqLenSum = 0;

        kmerThrMid = kmerThrMin + (kmerThrMax - kmerThrMin)*3/4;

        Debug(Debug::INFO) << "k-mer threshold range: [" << kmerThrMin  << ":" << kmerThrMax << "], trying threshold " << kmerThrMid << "\n";
        // determine k-mer match probability for kmerThrMid
#pragma omp parallel for schedule(static) 
        for (int i = 0; i < threads; i++){
            int thread_idx = 0;

#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // set a current k-mer list length threshold and a high prefitlering threshold (we don't need the prefiltering results in this test run)
            matchers[thread_idx] = new QueryTemplateMatcher(subMat, _2merSubMatrix, _3merSubMatrix, indexTable, dbr->getSeqLens(), kmerThrMid, 1.0, kmerSize, dbr->getSize(), aaBiasCorrection, maxSeqLen, 500.0);
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

            matchers[thread_idx]->matchQuery(seqs[thread_idx], UINT_MAX);

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

        // check the parameters
        double timeval = alpha * kmersPerPos + beta * kmerMatchProb + gamma;
        Debug(Debug::INFO) << "\tk-mers per position = " << kmersPerPos << ", k-mer match probability: " << kmerMatchProb << "\n";
        Debug(Debug::INFO) << "\ttime value = " << timeval << ", allowed range: [" << timevalMin << ":" << timevalMax << "]\n";
        if (timeval < timevalMin){
            if ((timevalMin - timeval) < (timevalMin - timevalBest) || (timevalMin - timeval) < (timevalBest - timevalMax)){
                // save new best values
                timevalBest = timeval;
                kmerThrBest = kmerThrMid;
                kmerMatchProbBest = kmerMatchProb;
            }
            kmerThrMax = kmerThrMid - 1;
        }
        else if (timeval > timevalMax){
            if ((timeval - timevalMax) < (timevalMin - timevalBest) || (timeval - timevalMax) < (timevalBest - timevalMax)){
                // save new best values
                timevalBest = timeval;
                kmerThrBest = kmerThrMid;
                kmerMatchProbBest = kmerMatchProb;
            }
            kmerThrMin = kmerThrMid + 1;
        }
        else if (timeval >= timevalMin && timeval <= timevalMax){
            // delete data structures used before returning
            delete[] querySeqs;
            delete[] matchers;
            delete indexTable;
            Debug(Debug::WARNING) << "\nk-mer threshold set, yielding sensitivity " << (log(timeval)/log(base)) << "\n\n";
            return std::pair<short, double> (kmerThrMid, kmerMatchProb);
        }
    }
    delete[] querySeqs;
    delete[] matchers;
    delete indexTable;

    Debug(Debug::WARNING) << "\nCould not set the k-mer threshold to meet the time value. Using the best value obtained so far, yielding sensitivity = " << (log(timevalBest)/log(base)) << "\n\n";
    return std::pair<short, double> (kmerThrBest, kmerMatchProbBest);
}
