#include "Prefiltering.h"
#include "PrefilteringIndexReader.h"
#include "../commons/Util.h"
#include "IndexTableGlobal.h"
#include "IndexTableLocal.h"
#include "QueryTemplateMatcherGlobal.h"
#include "QueryTemplateMatcherLocal.h"
#include "QueryTemplateMatcher.h"


Prefiltering::Prefiltering(std::string queryDB,
                           std::string queryDBIndex,
                           std::string targetDB,
                           std::string targetDBIndex,
                           std::string outDB,
                           std::string outDBIndex,
                           Parameters par):
    outDB(outDB),
    outDBIndex(outDBIndex),
    kmerSize(par.kmerSize),
    spacedKmer(par.spacedKmer),
    sensitivity(par.sensitivity),
    maxResListLen(par.maxResListLen),
    alphabetSize(par.alphabetSize),
    zscoreThr(par.zscoreThr),
    maxSeqLen(par.maxSeqLen),
    querySeqType(par.querySeqType),
    targetSeqType(par.targetSeqType),
    aaBiasCorrection(par.compBiasCorrection),
    split(par.split),
    skip(par.skip),
    isLocal(par.localSearch)
{
    if(this->split == 0 )
        this->split = 1;

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
    templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tdbr);
    if(templateDBIsIndex == true){ // exchange reader with old ffindex reader
        this->tidxdbr = this->tdbr;
        this->tdbr = PrefilteringIndexReader::openNewReader(tdbr);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
        this->kmerSize     = data.kmerSize;
        this->alphabetSize = data.alphabetSize;
        this->skip         = data.skip;
        this->split        = data.split;
        this->spacedKmer   = (data.spacedKmer == 1) ? true : false;
        this->isLocal      = (data.local == 1) ? true : false;
    }
    
    DBWriter::errorIfFileExist(outDB.c_str());
    DBWriter::errorIfFileExist(outDBIndex.c_str());
    
    Debug(Debug::INFO) << "Query database: " << par.db1 << "(size=" << qdbr->getSize() << ")\n";
    Debug(Debug::INFO) << "Target database: " << par.db2 << "(size=" << tdbr->getSize() << ")\n";

    // init the substitution matrices
    switch (querySeqType) {
        case Sequence::NUCLEOTIDES:
            subMat = new NucleotideMatrix();
            _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
            _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);
            break;
        case Sequence::AMINO_ACIDS:
            subMat = getSubstitutionMatrix(par.scoringMatrixFile, alphabetSize, 8.0);
            _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
            _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);
            break;
        case Sequence::HMM_PROFILE:
            subMat = getSubstitutionMatrix(par.scoringMatrixFile, alphabetSize, 8.0); // needed for Background distrubutions
            _2merSubMatrix = NULL;
            _3merSubMatrix = NULL;
            break;
    }
    
    // init all thread-specific data structures 
    this->qseq = new Sequence*[threads];
    this->reslens = new std::list<int>*[threads];
    this->notEmpty = new int[qdbr->getSize()];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif 
        qseq[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, querySeqType, kmerSize, spacedKmer);
        reslens[thread_idx] = new std::list<int>();
    }

}

Prefiltering::~Prefiltering(){
    for (int i = 0; i < threads; i++){
        delete qseq[i];
        reslens[i]->clear();
        delete reslens[i];
    }
    delete[] qseq;
    delete[] reslens;
    delete[] notEmpty;
    delete qdbr;
    delete tdbr;
    if(templateDBIsIndex == true)
        delete tidxdbr;
    delete subMat;
    if(_2merSubMatrix != NULL)
        delete _2merSubMatrix;
    if(_3merSubMatrix != NULL)
        delete _3merSubMatrix;
}


void Prefiltering::run(){
    // splits template database into x sequence steps
    unsigned int splitCounter =  this->split;
    std::vector<std::pair<std::string, std::string> > splitFiles;
    for(unsigned int split = 0; split < splitCounter; split++){
        std::pair<std::string, std::string> filenamePair = createTmpFileNames(outDB,outDBIndex,split);
        splitFiles.push_back(filenamePair);

        this->run (split, splitCounter,
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

void Prefiltering::mergeOutput(std::vector<std::pair<std::string, std::string> > filenames){
    DBWriter writer(outDB.c_str(), outDBIndex.c_str());
    writer.open();
    writer.mergeFiles(qdbr, filenames, BUFFER_SIZE);
    writer.close();
}

std::pair<std::string, std::string> Prefiltering::createTmpFileNames(std::string db, std::string dbindex, int numb){
    std::string splitSuffix = std::string("_tmp_") + SSTR(numb);
    std::string dataFile  = db + splitSuffix;
    std::string indexFile = dbindex + splitSuffix;
    return std::make_pair(dataFile, indexFile);
}

void Prefiltering::run(int mpi_rank, int mpi_num_procs){

    std::pair<std::string, std::string> filenamePair = createTmpFileNames(outDB, outDBIndex, mpi_rank);

    this->run (mpi_rank, mpi_num_procs,
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


QueryTemplateMatcher** Prefiltering::createQueryTemplateMatcher(BaseMatrix *m,
                                                                IndexTable *indexTable,
                                                                unsigned int *seqLens,
                                                                short kmerThr,
                                                                double kmerMatchProb,
                                                                int kmerSize,
                                                                int dbSize,
                                                                bool aaBiasCorrection,
                                                                int maxSeqLen,
                                                                float zscoreThr){
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher*[threads];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < this->threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        if (isLocal) {
            matchers[thread_idx] = new QueryTemplateMatcherLocal(m, indexTable, seqLens, kmerThr,
                                                                  kmerMatchProb, kmerSize, dbSize,
                                                                  aaBiasCorrection, maxSeqLen, zscoreThr);

        }else{
            matchers[thread_idx] = new QueryTemplateMatcherGlobal(m, indexTable, seqLens, kmerThr,
                                                                  kmerMatchProb, kmerSize, dbSize,
                                                                  aaBiasCorrection, maxSeqLen, zscoreThr);
        }
        if(querySeqType == Sequence::HMM_PROFILE){
            matchers[thread_idx]->setProfileMatrix(qseq[thread_idx]->profile_matrix);
        }else {
            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix->scoreMatrix, _2merSubMatrix->scoreMatrix );
        }
    }
    return matchers;
}

IndexTable * Prefiltering::getIndexTable(int split, int splitCount){
    if(templateDBIsIndex == true ){
        return PrefilteringIndexReader::generateIndexTable(tidxdbr,split);
    }else{
        size_t dbFrom = 0;
        size_t dbSize = 0;
        Util::decomposeDomainByAminoaAcid(tdbr->getAminoAcidDBSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                          split, splitCount, &dbFrom, &dbSize);
        Sequence tseq(maxSeqLen, subMat->aa2int, subMat->int2aa, targetSeqType, kmerSize, spacedKmer);
        return generateIndexTable(tdbr, &tseq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize, isLocal, skip);
    }
}

void Prefiltering::run (size_t split, size_t splitCount,
                        std::string resultDB, std::string resultDBIndex){

    Debug(Debug::WARNING) << "Process prefiltering step " << split << " of " << splitCount <<  "\n\n";

    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), threads);
    tmpDbw.open();
    size_t queryDBSize = qdbr->getSize();
    
    memset(notEmpty, 0, queryDBSize * sizeof(int)); // init notEmpty
    
    IndexTable * indexTable = getIndexTable(split, splitCount);

    // set the k-mer similarity threshold
    Debug(Debug::INFO) << "\nAdjusting k-mer similarity threshold within +-10% deviation from the reference time value, sensitivity = " << sensitivity << ")...\n";
    std::pair<short, double> ret = setKmerThreshold (indexTable, qdbr, tdbr, sensitivity, 0.1);
    //std::pair<short, double> ret = std::pair<short, double>(105, 8.18064e-05);
    //std::pair<short, double> ret = std::pair<short, double>(88, 8.18064e-05);
    //std::pair<short, double> ret = std::pair<short, double>(80, 8.18064e-05);
    this->kmerThr = ret.first;
    this->kmerMatchProb = ret.second;

    Debug(Debug::WARNING) << "k-mer similarity threshold: " << kmerThr << "\n";
    Debug(Debug::WARNING) << "k-mer match probability: " << kmerMatchProb << "\n\n";

    
    struct timeval start, end;

    gettimeofday(&start, NULL);
    QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat,indexTable, tdbr->getSeqLens(), kmerThr,
                             kmerMatchProb, kmerSize, tdbr->getSize(),
                             aaBiasCorrection, maxSeqLen, zscoreThr);

    size_t kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step "<< split << " of " << splitCount << ")\n";

#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches)
    for (size_t id = 0; id < queryDBSize; id++){ 
        Log::printProgress(id);
        
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get query sequence
        char* seqData = qdbr->getData(id);
        qseq[thread_idx]->mapSequence(id, qdbr->getDbKey(id), seqData);
        
        // calculate prefitlering results
        std::pair<hit_t *, size_t> prefResults = matchers[thread_idx]->matchQuery(qseq[thread_idx],
                                                                                  tdbr->getId(qseq[thread_idx]->getDbKey()));
        const size_t resultSize = prefResults.second;
        // write
        if(writePrefilterOutput(&tmpDbw, thread_idx, id, prefResults) != 0)
            continue; // couldnt write result because of too much results
        
        // update statistics counters
        if (resultSize != 0)
            notEmpty[id] = 1;

        kmersPerPos += (size_t) matchers[thread_idx]->getStatistics()->kmersPerPos;
        dbMatches += matchers[thread_idx]->getStatistics()->dbMatches;
        doubleMatches += matchers[thread_idx]->getStatistics()->doubleMatches;
        resSize += resultSize;
        realResSize += std::min(resultSize, maxResListLen);
        reslens[thread_idx]->push_back(resultSize);
    } // step end

    this->kmersPerPos = kmersPerPos;
    this->dbMatches = dbMatches;
    this->resSize = resSize;
    this->realResSize = realResSize;
    this->doubleMatches = doubleMatches;

    if (queryDBSize > 1000)
        Debug(Debug::INFO) << "\n";
    Debug(Debug::WARNING) << "\n";
    
    for (int j = 0; j < threads; j++){
        delete matchers[j];
    }
    delete[] matchers;
    delete indexTable;
    
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nTime for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    
    tmpDbw.close(); // sorts the index

}

void Prefiltering::closeReader(){
    qdbr->close();
    //if (strcmp(qdbr->getIndexFileName(), tdbr->getIndexFileName()) != 0)
    tdbr->close();
    if(templateDBIsIndex)
        tidxdbr->close();
}

void Prefiltering::removeDatabaes(std::vector<std::pair<std::string, std::string> > filenames) {
    for (size_t i = 0; i < filenames.size(); i++) {
        remove(filenames[i].first.c_str());
        remove(filenames[i].second.c_str());
    }
}

// write prefiltering to ffindex database
int Prefiltering::writePrefilterOutput(DBWriter * dbWriter, int thread_idx, size_t id, std::pair<hit_t *,size_t> prefResults){
    // write prefiltering results to a string
    size_t l = 0;
    hit_t * resultVector = prefResults.first;
    const size_t resultSize = prefResults.second;
    std::string prefResultsOutString;
    prefResultsOutString.reserve(BUFFER_SIZE);
    char buffer [100];

    for (size_t i = 0; i < resultSize; i++){
        hit_t * res = resultVector + i;

        if (res->seqId >= tdbr->getSize()) {
            Debug(Debug::INFO) << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id)<< " -> " << res->seqId << "\t" << res->prefScore << "\n";
        }
        const int len = snprintf(buffer,100,"%s\t%.4f\t%d\n",tdbr->getDbKey(res->seqId), res->zScore, res->prefScore);
        prefResultsOutString.append( buffer, len );
        l++;
        // maximum allowed result list length is reached
        if (l == this->maxResListLen)
            break;
    }
    // write prefiltering results string to ffindex database
    const size_t prefResultsLength = prefResultsOutString.length();
    if (BUFFER_SIZE < prefResultsLength){
        Debug(Debug::ERROR) << "Tried to process the prefiltering list for the query " << qdbr->getDbKey(id) << " , the length of the list = " << resultSize << "\n";
        Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << BUFFER_SIZE << " < " << prefResultsLength << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
        return -1;
    }
    char* prefResultsOutData = (char *) prefResultsOutString.c_str();
    dbWriter->write(prefResultsOutData, prefResultsLength, qdbr->getDbKey(id), thread_idx);
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
//    size_t prefRealPassedPerSeq = realResSize/queryDBSize;
    size_t doubleMatcherPerQuerySeq = doubleMatches/queryDBSize;
    Debug(Debug::INFO) << kmersPerPos/queryDBSize << " k-mers per position.\n";
    Debug(Debug::INFO) << dbMatchesPerSeq << " DB matches per sequence.\n";
    if(isLocal){
         Debug(Debug::INFO) << doubleMatcherPerQuerySeq << " Double diagonal matches per sequence.\n";
    }
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

BaseMatrix* Prefiltering::getSubstitutionMatrix(std::string scoringMatrixFile, int alphabetSize, float bitFactor){
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


void Prefiltering::countKmersForIndexTable (DBReader* dbr, Sequence* seq,
                                                   IndexTable* indexTable,
                                                   size_t dbFrom, size_t dbTo){
    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    dbTo=std::min(dbTo,dbr->getSize());
    for (unsigned int id = dbFrom; id < dbTo; id++){
        Log::printProgress(id - dbFrom);
        char* seqData = dbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(id, dbr->getDbKey(id), seqData);
        indexTable->addKmerCount(seq);
    }
}


void Prefiltering::fillDatabase(DBReader* dbr, Sequence* seq, IndexTable * indexTable,
                                size_t dbFrom, size_t dbTo)
{
    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";
    indexTable->initMemory();
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
    indexTable->revertPointer();
    indexTable->removeDuplicateEntries();
    Debug(Debug::INFO) << "Index table init done.\n\n";

}


IndexTable* Prefiltering::generateIndexTable (DBReader* dbr, Sequence* seq, int alphabetSize,
        int kmerSize, size_t dbFrom, size_t dbTo, bool isLocal, int skip){

    struct timeval start, end;
    gettimeofday(&start, NULL);
    IndexTable * indexTable;
	if (isLocal) {
        indexTable = new IndexTableLocal(alphabetSize, kmerSize, skip);
    } else{
        indexTable = new IndexTableGlobal(alphabetSize, kmerSize, skip);

    }

    countKmersForIndexTable(dbr, seq, indexTable, dbFrom, dbTo);

    Debug(Debug::INFO) << "\n";
    
    fillDatabase(dbr, seq, indexTable, dbFrom, dbTo);
    
    gettimeofday(&end, NULL);
    
    indexTable->printStatisitic(seq->int2aa);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for index table init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    return indexTable;
}

std::pair<short,double> Prefiltering::setKmerThreshold (IndexTable * indexTable, DBReader* qdbr, DBReader* tdbr,
                                                        double sensitivity, double toleratedDeviation){

    size_t targetDbSize = indexTable->getSize();
    size_t targetSeqLenSum = 0;
    for (size_t i = 0; i < targetDbSize; i++)
        targetSeqLenSum += tdbr->getSeqLens()[i];

    // generate a small random sequence set for testing 
    size_t querySetSize = std::min ( tdbr->getSize(), (size_t) 1000);
    int* querySeqs = new int[querySetSize];
    srand(1);
    for (size_t i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % tdbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability 
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb
    size_t kmerThrMin = 3 * kmerSize;
    size_t kmerThrMax = 80 * kmerSize;
    size_t kmerThrMid;

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
        EXIT(EXIT_FAILURE);
    }
    //alpha = 0; //7.123599e-02; //6.438574e-02; // 6.530289e-02;
    //beta = 1; //3.148479e+06; //3.480680e+06; // 3.243035e+06;
    //gamma = 0; //1.304421; // 1.753651; //1.137125;


    // Run using k=6, a=21, with k-mer similarity threshold 103
    // k-mer list length was 117.6, k-mer match probability 1.12735e-06
    // time reference value for these settings is 15.6
    // Since we want to represent the time value in the form base^sensitivity with base=2.0, it yields sensitivity = 3.96 (~4.0)
    double base = 2.0;
    double timevalMax = pow(base, sensitivity) * (1.0 + toleratedDeviation);
    double timevalMin = pow(base, sensitivity) * (1.0 - toleratedDeviation);

    // in case the time value cannot be set within the threshold boundaries, the best value that could be reached will be returned with a warning
    double timevalBest = 0.0;
    size_t kmerThrBest = 0;
    double kmerMatchProbBest = 0.0;

    // adjust k-mer list length threshold
    while (kmerThrMax >= kmerThrMin){
        dbMatchesSum = 0;
        querySeqLenSum = 0;

        kmerThrMid = kmerThrMin + (kmerThrMax - kmerThrMin)*3/4;

        Debug(Debug::INFO) << "k-mer threshold range: [" << kmerThrMin  << ":" << kmerThrMax << "], trying threshold " << kmerThrMid << "\n";
        // determine k-mer match probability for kmerThrMid
        QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable,
                                 tdbr->getSeqLens(), kmerThrMid, 1.0, kmerSize, tdbr->getSize(),
                                 aaBiasCorrection, maxSeqLen, 500.0);

#pragma omp parallel for schedule(dynamic, 10) reduction (+: dbMatchesSum, querySeqLenSum, kmersPerPos)
        for (size_t i = 0; i < querySetSize; i++){
            int id = querySeqs[i];

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char* seqData = tdbr->getData(id);
            qseq[thread_idx]->mapSequence(id, tdbr->getDbKey(id), seqData);

            matchers[thread_idx]->matchQuery(qseq[thread_idx],UINT_MAX);

            kmersPerPos += matchers[thread_idx]->getStatistics()->kmersPerPos;
            dbMatchesSum += matchers[thread_idx]->getStatistics()->dbMatches;
            querySeqLenSum += qseq[thread_idx]->L;
        }

        kmersPerPos /= (double)querySetSize;

        // add pseudo-counts
        dbMatchesExp_pc = (size_t)(((double)lenSum_pc) * kmersPerPos * pow((1.0/((double)(subMat->alphabetSize-1))), kmerSize));

        // match probability with pseudocounts
        kmerMatchProb = ((double)dbMatchesSum + dbMatchesExp_pc) / ((double) (querySeqLenSum * targetSeqLenSum + lenSum_pc));

        for (int j = 0; j < threads; j++){
            delete matchers[j];
        }
        delete[] matchers;

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
            Debug(Debug::WARNING) << "\nk-mer threshold set, yielding sensitivity " << (log(timeval)/log(base)) << "\n\n";
            return std::pair<short, double> (kmerThrMid, kmerMatchProb);
        }
    }
    delete[] querySeqs;

    Debug(Debug::WARNING) << "\nCould not set the k-mer threshold to meet the time value. Using the best value obtained so far, yielding sensitivity = " << (log(timevalBest)/log(base)) << "\n\n";
    return std::pair<short, double> (kmerThrBest, kmerMatchProbBest);
}
