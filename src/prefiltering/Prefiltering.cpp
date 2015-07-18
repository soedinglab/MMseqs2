#include <cstddef>

#include "Prefiltering.h"
#include "PrefilteringIndexReader.h"
#include "Util.h"
#include "IndexTableGlobal.h"
#include "IndexTableLocal.h"
#include "QueryTemplateMatcherGlobal.h"
#include "QueryTemplateMatcherExactMatch.h"
#include "QueryTemplateMatcherLocal.h"
#include "QueryTemplateMatcher.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

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
        kmerScore(par.kmerScore),
        spacedKmer(par.spacedKmer),
        sensitivity(par.sensitivity),
        maxResListLen(par.maxResListLen),
        alphabetSize(par.alphabetSize),
        zscoreThr(par.zscoreThr),
        maxSeqLen(par.maxSeqLen),
        querySeqType(par.querySeqType),
        targetSeqType(par.targetSeqType),
        aaBiasCorrection(par.compBiasCorrection),
        fastMode(par.fastMode),
        split(par.split),
        skip(par.skip),
        searchMode(par.searchMode)
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
        this->searchMode   = (data.local == 1) ? ((par.searchMode >= 1) ? par.searchMode : Parameters::SEARCH_LOCAL)
                                               : Parameters::SEARCH_GLOBAL;
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

    } // prefiltering scores calculation end

    // merge output ffindex databases
    if(splitCounter > 1){
        this->mergeOutput(splitFiles);
    }else{
        std::rename(splitFiles[0].first.c_str(),  outDB.c_str());
        std::rename(splitFiles[0].second.c_str(), outDBIndex.c_str());
    }
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


QueryTemplateMatcher ** Prefiltering::createQueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable,
                                                                 unsigned int *seqLens, short kmerThr,
                                                                 double kmerMatchProb, int kmerSize,
                                                                 size_t effectiveKmerSize, size_t dbSize,
                                                                 bool aaBiasCorrection, bool fastMode,
                                                                 unsigned int maxSeqLen, float zscoreThr, int searchMode,
                                                                 size_t maxHitsPerQuery) {
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher *[threads];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < this->threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        if (searchMode == Parameters::SEARCH_LOCAL) {
            matchers[thread_idx] = new QueryTemplateMatcherLocal(m, indexTable, seqLens, kmerThr,
                                                                 kmerMatchProb, kmerSize, effectiveKmerSize, dbSize,
                                                                 fastMode, maxSeqLen, maxHitsPerQuery);
        } else if(searchMode == Parameters::SEARCH_LOCAL_FAST) {
            matchers[thread_idx] = new QueryTemplateMatcherExactMatch(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                                      kmerSize, dbSize, maxSeqLen, maxHitsPerQuery);
        } else {
            matchers[thread_idx] = new QueryTemplateMatcherGlobal(m, indexTable, seqLens, kmerThr,
                                                                  kmerMatchProb, kmerSize, effectiveKmerSize, dbSize,
                                                                  aaBiasCorrection, maxSeqLen, zscoreThr);
        }
        if(querySeqType == Sequence::HMM_PROFILE){
            matchers[thread_idx]->setProfileMatrix(qseq[thread_idx]->profile_matrix);
        } else {
            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix->scoreMatrix, _2merSubMatrix->scoreMatrix );
        }
    }
    return matchers;
}
IndexTable * Prefiltering::getIndexTable(int split, int splitCount){
    if(templateDBIsIndex == true ){
        return PrefilteringIndexReader::generateIndexTable(tidxdbr, split);
    }else{
        size_t dbFrom = 0;
        size_t dbSize = 0;
        Util::decomposeDomainByAminoaAcid(tdbr->getAminoAcidDBSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                          split, splitCount, &dbFrom, &dbSize);
        Sequence tseq(maxSeqLen, subMat->aa2int, subMat->int2aa, targetSeqType, kmerSize, spacedKmer);
        return generateIndexTable(tdbr, &tseq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize, searchMode, skip);
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
    std::pair<short, double> calibration;
    calibration = setKmerThreshold(indexTable, qdbr, tdbr, sensitivity, 0.1, kmerScore);
    //std::pair<short, double> ret = std::pair<short, double>(105, 8.18064e-05);
    //std::pair<short, double> ret = std::pair<short, double>(88, 8.18064e-05);
    //std::pair<short, double> ret = std::pair<short, double>(80, 8.18064e-05);
    this->kmerThr = calibration.first;
    this->kmerMatchProb = calibration.second;

    Debug(Debug::WARNING) << "k-mer similarity threshold: " << kmerThr << "\n";
    Debug(Debug::WARNING) << "k-mer match probability: " << kmerMatchProb << "\n\n";


    struct timeval start, end;

    gettimeofday(&start, NULL);
    QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable, tdbr->getSeqLens(), kmerThr,
                                                                  kmerMatchProb, kmerSize,
                                                                  qseq[0]->getEffectiveKmerSize(), tdbr->getSize(),
                                                                  aaBiasCorrection, fastMode, maxSeqLen, zscoreThr,
                                                                  searchMode, maxResListLen);

    size_t kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step "<< split << " of " << splitCount << ")\n";

#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
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
        querySeqLenSum += qseq[thread_idx]->L;
        diagonalOverflow += matchers[thread_idx]->getStatistics()->diagonalOverflow;
        resSize += resultSize;
        realResSize += std::min(resultSize, maxResListLen);
        reslens[thread_idx]->push_back(resultSize);
    } // step end

    statistics_t stats(kmersPerPos/queryDBSize, dbMatches/queryDBSize, doubleMatches/queryDBSize,
                       querySeqLenSum, diagonalOverflow, resSize/queryDBSize);

    this->printStatistics(stats);


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
    // needed to speed up merge later one
    // this sorts this datafile according to the index file

    if(splitCount > 1){
        DBReader tmpDbr(tmpDbw.getDataFileName(), tmpDbw.getIndexFileName());
        DBWriter tmpDbw2(resultDB.c_str(), resultDBIndex.c_str(), threads);
        tmpDbw2.open();
        tmpDbw2.sortDatafileByIdOrder(&tmpDbr);
        tmpDbr.close();
        tmpDbw2.close();
    }
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


void Prefiltering::printStatistics(statistics_t &stats){

    size_t empty = 0;
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

//    size_t prefRealPassedPerSeq = realResSize/queryDBSize;
    Debug(Debug::INFO) << stats.kmersPerPos << " k-mers per position.\n";
    Debug(Debug::INFO) << stats.dbMatches << " DB matches per sequence.\n";
    if(searchMode){
        Debug(Debug::INFO) << stats.doubleMatches << " Double diagonal matches per sequence.\n";
        Debug(Debug::INFO) << stats.diagonalOverflow << " Overflows .\n";
    }
    Debug(Debug::INFO) << stats.resultsPassedPrefPerSeq << " sequences passed prefiltering per query sequence";
    if (stats.resultsPassedPrefPerSeq > maxResListLen)
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
        subMat = new ReducedMatrix(sMat->probMatrix, sMat->subMatrixPseudoCounts, alphabetSize, bitFactor);
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
        //TODO - dbFrom?!?
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
                                              int kmerSize, size_t dbFrom, size_t dbTo, int searchMode, int skip){

    struct timeval start, end;
    gettimeofday(&start, NULL);
    IndexTable * indexTable;
    if (searchMode == Parameters::SEARCH_LOCAL || searchMode == Parameters::SEARCH_LOCAL_FAST) {
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

std::pair<short, double> Prefiltering::setKmerThreshold(IndexTable *indexTable, DBReader *qdbr, DBReader *tdbr, float sensitivity,
                                                        double toleratedDeviation, const int kmerScore) {

    size_t targetDbSize = indexTable->getSize();
    size_t targetSeqLenSum = 0;
    for (size_t i = 0; i < targetDbSize; i++){
        targetSeqLenSum += (tdbr->getSeqLens()[i] - qseq[0]->getEffectiveKmerSize());
    }
    // generate a small random sequence set for testing
    size_t querySetSize = std::min ( tdbr->getSize(), (size_t) 1000);
    unsigned int * querySeqs = new unsigned int[querySetSize];
    srand(1);
    for (size_t i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % qdbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb
    size_t kmerThrBest;

    size_t dbMatchesExp_pc;
    // 1000 * 350 * 100000 * 350
    size_t lenSum_pc = 12250000000000;

    double kmerMatchProb;

    const unsigned int sens =  (int) ceil(sensitivity);

    kmerThrBest =  kmerScore;
    if(kmerThrBest == INT_MAX){
        if (kmerSize == 6){
            switch(sens){
                case 1:
                    kmerThrBest = 125;
                    break;
                case 2:
                    kmerThrBest = 115;
                    break;
                case 3:
                    kmerThrBest = 110;
                    break;
                case 4:
                    kmerThrBest = 100;
                    break;
                case 5:
                    kmerThrBest = 95;
                    break;
                case 6:
                    kmerThrBest = 90;
                    break;
                case 7:
                    kmerThrBest = 85;
                    break;
                case 8:
                    kmerThrBest = 80;
                    break;
                case 9:
                    kmerThrBest = 70;
                    break;
                default:
                    kmerThrBest = 100;
                    break;
            }

        }
        else if (kmerSize == 7){
            switch(sens){
                case 1:
                    kmerThrBest = 130;
                    break;
                case 2:
                    kmerThrBest = 120;
                    break;
                case 3:
                    kmerThrBest = 110;
                    break;
                case 4:
                    kmerThrBest = 100;
                    break;
                case 5:
                    kmerThrBest = 95;
                    break;
                case 6:
                    kmerThrBest = 90;
                    break;
                case 7:
                    kmerThrBest = 85;
                    break;
                case 8:
                    kmerThrBest = 80;
                    break;
                case 9:
                    kmerThrBest = 75;
                    break;
                default:
                    kmerThrBest = 100;
                    break;
            }
        }
        else{
            Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    //alpha = 0; //7.123599e-02; //6.438574e-02; // 6.530289e-02;


    Debug(Debug::INFO) << "k-mer threshold threshold: " << kmerThrBest << "\n";
    statistics_t stats;
    stats = computeStatisticForKmerThreshold(indexTable, querySetSize, querySeqs, searchMode, kmerThrBest);
    // match probability with pseudocounts
    // add pseudo-counts (1/20^6 * kmerPerPos * length of pseudo counts)
    dbMatchesExp_pc = (size_t)(((double)lenSum_pc) * stats.kmersPerPos * pow((1.0/((double)(subMat->alphabetSize-1))), kmerSize));
    // match probability with pseudocounts
    kmerMatchProb = ((double)stats.dbMatches + dbMatchesExp_pc) / ((double) (stats.querySeqLen * targetSeqLenSum + lenSum_pc));
    // compute match prob for local match
    if(searchMode == true){
        kmerMatchProb = ((double) stats.doubleMatches) / ((double) (stats.querySeqLen * targetSeqLenSum));
        kmerMatchProb /= 256;
    }
    kmerMatchProb = std::max(kmerMatchProb, std::numeric_limits<double>::min());
    Debug(Debug::INFO) << "\tk-mers per position = " << stats.kmersPerPos << ", k-mer match probability: " << kmerMatchProb << "\n";


    delete[] querySeqs;

    return std::pair<short, double> (kmerThrBest, kmerMatchProb);
}

statistics_t Prefiltering::computeStatisticForKmerThreshold(IndexTable *indexTable, size_t querySetSize,
                                                            unsigned int *querySeqsIds, bool reverseQuery, size_t kmerThrMid) {
    // determine k-mer match probability for kmerThrMid
    QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable, tdbr->getSeqLens(), kmerThrMid,
                                                                  1.0, kmerSize, qseq[0]->getEffectiveKmerSize(),
                                                                  tdbr->getSize(), aaBiasCorrection, false, maxSeqLen,
                                                                  500.0, searchMode, LONG_MAX);
    size_t dbMatchesSum = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t diagonalOverflow = 0;
    size_t resultsPassedPref = 0;
    double kmersPerPos = 0.0;

#pragma omp parallel for schedule(dynamic, 10) reduction (+: dbMatchesSum, doubleMatches, kmersPerPos, querySeqLenSum, diagonalOverflow, resultsPassedPref)
    for (size_t i = 0; i < querySetSize; i++){
        size_t id = querySeqsIds[i];

        int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
        char* seqData = qdbr->getData(id);
        qseq[thread_idx]->mapSequence(id, qdbr->getDbKey(id), seqData);
        if(reverseQuery == true){
            qseq[thread_idx]->reverse();
        }
        matchers[thread_idx]->matchQuery(qseq[thread_idx], UINT_MAX);

        kmersPerPos    += matchers[thread_idx]->getStatistics()->kmersPerPos;
        dbMatchesSum   += matchers[thread_idx]->getStatistics()->dbMatches;
        querySeqLenSum += matchers[thread_idx]->getStatistics()->querySeqLen;
        doubleMatches  += matchers[thread_idx]->getStatistics()->doubleMatches;
        diagonalOverflow += matchers[thread_idx]->getStatistics()->diagonalOverflow;
        resultsPassedPref += matchers[thread_idx]->getStatistics()->resultsPassedPrefPerSeq;
    }
    // clean up memory
    for (int j = 0; j < threads; j++){
        delete matchers[j];
    }
    delete[] matchers;

    return statistics_t(kmersPerPos / (double)querySetSize, dbMatchesSum, doubleMatches,
                        querySeqLenSum, diagonalOverflow, resultsPassedPref/ querySetSize);
}
