#include "Prefiltering.h"


Prefiltering::Prefiltering(std::string queryDB,
        std::string queryDBIndex,
        std::string targetDB,
        std::string targetDBIndex,
        std::string outDB,
        std::string outDBIndex,
        std::string scoringMatrixFile,
        float sensitivity,
        int kmerSize,
        int maxResListLen,
        int alphabetSize,
        float zscoreThr,
        size_t maxSeqLen,
        int seqType,
        bool aaBiasCorrection,
        int splitSize,
        int skip):    outDB(outDB),
    outDBIndex(outDBIndex),
    kmerSize(kmerSize),
    maxResListLen(maxResListLen),
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

    
    DBWriter::errorIfFileExist(outDB.c_str());
    DBWriter::errorIfFileExist(outDBIndex.c_str());
    
    if (splitSize == 0)
        splitSize = tdbr->getSize();

    
    Debug(Debug::INFO) << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";
    Debug(Debug::INFO) << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

    // init the substitution matrices
    if (seqType == Sequence::NUCLEOTIDES)
        subMat = new NucleotideMatrix();
    else
        subMat = getSubstitutionMatrix(scoringMatrixFile, 8.0);


    _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);

    // init all thread-specific data structures 
    this->seqs = new Sequence*[threads];
    this->reslens = new std::list<int>*[threads];
    this->notEmpty = new int[qdbr->getSize()];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType, subMat);
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
        delete matchers[i];
        delete reslens[i];
    }
    delete[] seqs;
    delete[] outBuffers;
    delete[] matchers;
    delete[] reslens;
    delete notEmpty;

    delete indexTable;
    delete subMat;
    delete _2merSubMatrix;
    delete _3merSubMatrix;
}


void Prefiltering::run(){
    // splits template database into x sequence steps
    int step = 0;
    int stepCnt = (tdbr->getSize() + splitSize - 1) / splitSize;
    std::vector<std::pair<std::string, std::string> > splitFiles;
    for(unsigned int splitStart = 0; splitStart < tdbr->getSize(); splitStart += splitSize ){
        step++;
        Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step " << step << " of " << stepCnt <<  ")\n";
        std::pair<std::string, std::string> filenamePair = createTmpFileNames(outDB,outDBIndex,step);
        splitFiles.push_back(filenamePair);
        
        this->run (splitStart, splitSize,
                   filenamePair.first.c_str(),
                   filenamePair.second.c_str() );
        
        this->printStatistics();
    } // prefiltering scores calculation end
    
    // merge output ffindex databases
    this->mergeOutput(splitFiles);
    // remove temp databases
    this->removeDatabaes(splitFiles);
    // close reader to reduce memory
    this->closeReader();
}

std::pair<std::string, std::string> Prefiltering::createTmpFileNames(std::string db, std::string dbindex, int numb){
    std::string splitSuffix = "_tmp_" + SSTR(numb);
    std::string dataFile  = db + splitSuffix;
    std::string indexFile = dbindex + splitSuffix;
    return std::make_pair(dataFile, indexFile);
}

void Prefiltering::run(int mpi_rank, int mpi_num_procs){
    int splitStart, splitSize;
    
    Util::decompose_domain(tdbr->getSize(), mpi_rank,
                     mpi_num_procs, &splitStart,
                     &splitSize);
    
    std::pair<std::string, std::string> filenamePair = createTmpFileNames(outDB, outDBIndex, mpi_rank);

    this->run (splitStart, splitSize,
               filenamePair.first.c_str(),
               filenamePair.second.c_str());
    this->printStatistics();
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(mpi_rank == 0){ // master reduces results
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for(int procs = 0; procs < mpi_num_procs; procs++){
            splitFiles.push_back(createTmpFileNames(outDB, outDBIndex, procs));
        }
        // merge output ffindex databases
        this->mergeOutput(splitFiles);
        // remove temp databases
        this->removeDatabaes(splitFiles);
        // close reader to reduce memory
        this->closeReader();
    } else {
        // close reader to reduce memory
        this->closeReader();
    }
}





void Prefiltering::run (size_t dbFrom,size_t dbSize,
                        std::string resultDB, std::string resultDBIndex){
    
    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), threads);
    tmpDbw.open();
    size_t queryDBSize = qdbr->getSize();
    
    memset(notEmpty, 0, queryDBSize*sizeof(int)); // init notEmpty
    
    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, seqType, subMat);
    this->indexTable = getIndexTable(tdbr, seq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize , skip);
    delete seq;

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

    int kmersPerPos = 0;
    int dbMatches = 0;
    int resSize = 0;
#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, dbMatches)
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
        std::list<hit_t>* prefResults = matchers[thread_idx]->matchQuery(seqs[thread_idx]);
        // add id if exist in targetDb
        addIdIfExistInTargetDb(thread_idx, prefResults);
        
        // write
        if(writePrefilterOutput(&tmpDbw, thread_idx, id, prefResults)!=0)
            continue; // couldnt write result because of to much results
        
        // update statistics counters
        if (prefResults->size() != 0)
            notEmpty[id] = 1;
        kmersPerPos += (size_t) seqs[thread_idx]->stats->kmersPerPos;
        dbMatches += seqs[thread_idx]->stats->dbMatches;
        resSize += prefResults->size();
        reslens[thread_idx]->push_back(prefResults->size());
    } // step end

    this->kmersPerPos = kmersPerPos;
    this->dbMatches = dbMatches;
    this->resSize = resSize;
    if (queryDBSize > 1000)
        Debug(Debug::INFO) << "\n";
    Debug(Debug::WARNING) << "\n";
    
    for (int j = 0; j < threads; j++){
        delete matchers[j];
    }
    delete indexTable;
    
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nTime for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    
    tmpDbw.close(); // sorts the index

}

void Prefiltering::closeReader(){
    qdbr->close();
    if (strcmp(qdbr->getIndexFileName(), tdbr->getIndexFileName()) != 0)
        tdbr->close();
}

void Prefiltering::mergeOutput(std::vector<std::pair<std::string, std::string> > filenames){
    Debug(Debug::INFO) << "Merging the results...\n\n";
    const size_t file_count = filenames.size();

    DBWriter dbw(outDB.c_str(), outDBIndex.c_str(), 1);
    dbw.open();
    // open DBReader
    DBReader * filesToMerge [file_count];
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file] = new DBReader(filenames[file].first.c_str(),
                                          filenames[file].second.c_str());
        filesToMerge[file]->open(DBReader::NOSORT);
    }
    for (size_t id = 0; id < qdbr->getSize(); id++){
        std::stringstream mergeResultsOut;
        // get all data for the id from all files
        for(size_t file = 0; file < file_count; file++){
            mergeResultsOut << filesToMerge[file]->getData(id);
        }
        // create merged string
        std::string mergeResultsOutString = mergeResultsOut.str();
        if (BUFFER_SIZE < mergeResultsOutString.length()){
            Debug(Debug::ERROR) << "ERROR: Buffer overflow at id: " << qdbr->getDbKey(id) << " during the merging.\n";
            Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << mergeResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            continue; // read next id
        }
        // write result
        const char* mergeResultsOutData = mergeResultsOutString.c_str();
        memcpy(outBuffers[0], mergeResultsOutData, mergeResultsOutString.length() * sizeof(char));
        dbw.write(outBuffers[0], mergeResultsOutString.length(), qdbr->getDbKey(id), 0);
    }
    // close all reader
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file]->close();
        delete filesToMerge[file];
    }
    // sort result list
    dbw.close();
}

void Prefiltering::removeDatabaes(std::vector<std::pair<std::string, std::string> > filenames) {
    for (int i = 0; i < filenames.size(); i++) {
        remove(filenames[i].first.c_str());
        remove(filenames[i].second.c_str());
    }
}

void Prefiltering::addIdIfExistInTargetDb(int thread_idx, std::list<hit_t>* prefResults) {
    const size_t identityId = tdbr->getId(seqs[thread_idx]->getDbKey());
    if (identityId != UINT_MAX){
        std::list<hit_t>::iterator res = std::find_if(prefResults->begin(),
                                                      prefResults->end(),
                                                      HitIdEqual(identityId));
        if(res == prefResults->end()) {
            hit_t hit = {identityId, 0, 0};
            prefResults->push_front(hit);
        }else{
            hit_t identityHit = *res;
            prefResults->erase (res);
            prefResults->push_front(identityHit);
        }
    }
}

// write prefiltering to ffindex database
int Prefiltering::writePrefilterOutput(DBWriter * dbWriter, int thread_idx, size_t id, std::list<hit_t>* prefResults){
    // write prefiltering results to a string
    std::stringstream prefResultsOut;
    size_t l = 0;
    for (std::list<hit_t>::iterator iter = prefResults->begin(); iter != prefResults->end(); iter++){
        if (iter->seqId >= tdbr->getSize()){
            Debug(Debug::INFO) << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id)<< " -> " << iter->seqId << "\t" << iter->prefScore << "\n";
        }
        prefResultsOut << tdbr->getDbKey(iter->seqId) << "\t" << iter->zScore << "\t" << iter->prefScore << "\n";
        l++;
        // maximum allowed result list length is reached
        if (l == this->maxResListLen)
            break;
    }
    // write prefiltering results string to ffindex database
    std::string prefResultsOutString = prefResultsOut.str();
    const char* prefResultsOutData = prefResultsOutString.c_str();
    if (BUFFER_SIZE < strlen(prefResultsOutData)){
        Debug(Debug::ERROR) << "Tried to process the prefiltering list for the query " << qdbr->getDbKey(id) << " , the length of the list = " << prefResults->size() << "\n";
        Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
        return -1;
    }
    memcpy(outBuffers[thread_idx], prefResultsOutData, prefResultsOutString.length()*sizeof(char));
    dbWriter->write(outBuffers[thread_idx], prefResultsOutString.length(), qdbr->getDbKey(id), thread_idx);
    return 0;

}


void Prefiltering::printStatistics(){
    
    size_t queryDBSize = qdbr->getSize();
    int empty = 0;
    for (unsigned int i = 0; i < qdbr->getSize(); i++){
        if (notEmpty[i] == 0){
            //Debug(Debug::INFO) << "No prefiltering results for id " << i << ", " << qdbr->getDbKey(i) << ", len = " << strlen(qdbr->getData(i)) << "\n";
            empty++;
        }
    }
    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }
    
    size_t dbMatchesPerSeq = dbMatches/queryDBSize;
    size_t prefPassedPerSeq = resSize/queryDBSize;
    Debug(Debug::INFO) << kmersPerPos/queryDBSize << " k-mers per position.\n";
    Debug(Debug::INFO) << dbMatchesPerSeq << " DB matches per sequence.\n";
    Debug(Debug::INFO) << prefPassedPerSeq << " sequences passed prefiltering per query sequence";
    if (prefPassedPerSeq > maxResListLen)
        Debug(Debug::INFO) << " (ATTENTION: max. " << maxResListLen << " best scoring sequences were written to the output prefiltering database).\n";
    else
        Debug(Debug::INFO) << ".\n";

    int mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
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
        alpha = 1.141648e-01; // 1.114936e-01;
        beta = 9.033168e+05; // 9.331253e+05;
        gamma = 1.411142; // 1.416222;
    }
    else if (kmerSize == 7){ 
        alpha = 7.123599e-02; //6.438574e-02; // 6.530289e-02;
        beta = 3.148479e+06; //3.480680e+06; // 3.243035e+06;
        gamma = 1.304421; // 1.753651; //1.137125;
    }
    else{
        Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
        EXIT(EXIT_FAILURE);
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
