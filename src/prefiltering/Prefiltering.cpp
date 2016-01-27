#include "Prefiltering.h"
#include "PrefilteringIndexReader.h"
#include "Util.h"
#include "IndexTable.h"
#include "QueryTemplateLocalFast.h"
#include "QueryTemplateMatcherLocal.h"
#include "QueryTemplateMatcher.h"
#include "QueryScore.h"
#include <regex.h>

#include <cstddef>

#ifdef OPENMP
#include <omp.h>
#endif

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
        maxSeqLen(par.maxSeqLen),
        querySeqType(par.querySeqType),
        targetSeqType(par.targetSeqType),
        diagonalScoring(par.diagonalScoring),
        minDiagScoreThr(par.minDiagScoreThr),
        aaBiasCorrection(par.compBiasCorrection),
        split(par.split),
        splitMode(par.splitMode),
        searchMode(par.searchMode)
{
    if(this->split == 0 )
        this->split = 1;

    this->threads = par.threads;
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif
    Debug(Debug::INFO) << "\n";
    DBWriter::errorIfFileExist(outDB.c_str());
    DBWriter::errorIfFileExist(outDBIndex.c_str());
    this->qdbr = new DBReader<unsigned int>(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);

    this->tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open(DBReader<unsigned int>::NOSORT);
    std::string check(targetDB);
    check.replace(check.find(queryDB), queryDB.length(), "");
    regex_t regex;
    regcomp(&regex, "^\\.s?k[5-7]$", REG_EXTENDED | REG_NEWLINE);
    int nomatch = regexec(&regex, check.c_str(), 0, NULL, 0);
    regfree(&regex);
    // if no match found or two matches found (we want exactly one match)
    sameQTDB = (queryDB.compare(targetDB) == 0 || (nomatch == false) );

    templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tdbr);
    if(templateDBIsIndex == true){ // exchange reader with old ffindex reader
        this->tidxdbr = this->tdbr;
        this->tdbr = PrefilteringIndexReader::openNewReader(tdbr);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
        this->kmerSize     = data.kmerSize;
        this->alphabetSize = data.alphabetSize;
        this->split        = data.split;
        this->spacedKmer   = (data.spacedKmer == 1) ? true : false;
        this->searchMode   = (data.local == 1) ? ((par.searchMode >= 1) ? par.searchMode : Parameters::SEARCH_LOCAL)
                                               : Parameters::SEARCH_GLOBAL;
    }

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
            subMat = getSubstitutionMatrix(par.scoringMatrixFile, alphabetSize, 8.0, false);
            this->alphabetSize = subMat->alphabetSize;

            _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
            _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);
            break;
        case Sequence::HMM_PROFILE:
            subMat = getSubstitutionMatrix(par.scoringMatrixFile, alphabetSize, 8.0, false); // needed for Background distrubutions
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
    if(splitCounter > 1 && splitMode == Parameters::TARGET_DB_SPLIT) {
        maxResListLen = (maxResListLen / splitCounter)  + 1;
    }
    for(unsigned int split = 0; split < splitCounter; split++){
        std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(outDB,outDBIndex,split);
        splitFiles.push_back(filenamePair);

        this->run(split, splitCounter, splitMode, filenamePair.first.c_str(),
                  filenamePair.second.c_str());

    } // prefiltering scores calculation end

    // merge output ffindex databases
    if(splitCounter > 1){
        mergeFiles(splitFiles, splitMode);
    }else{
        std::rename(splitFiles[0].first.c_str(),  outDB.c_str());
        std::rename(splitFiles[0].second.c_str(), outDBIndex.c_str());
    }
//    // remove temp databases
//    this->removeDatabaes(splitFiles);
    // close reader to reduce memory
    this->closeReader();
}

void Prefiltering::mergeOutput(std::vector<std::pair<std::string, std::string> > filenames){

    struct timeval start, end;

    gettimeofday(&start, NULL);

    if(filenames.size() < 2){
        std::rename(filenames[0].first.c_str(),  outDB.c_str());
        std::rename(filenames[0].second.c_str(), outDBIndex.c_str());
        Debug(Debug::INFO) << "No mergeing needed.\n";
        return;
    }

    std::list<std::pair<std::string, std::string>> files(filenames.begin(), filenames.end());
    size_t mergeStep = 0;
    while(true){
        std::pair<std::string, std::string> file1 = files.front();
        files.pop_front();
        std::pair<std::string, std::string> file2 = files.front();
        files.pop_front();
        std::pair<std::string, std::string> out   = std::make_pair((outDB + "_merge_"+ SSTR(mergeStep)).c_str(), (outDBIndex + "_merge_"+ SSTR(mergeStep)).c_str());
        DBWriter writer(out.first.c_str(), out.second.c_str(), threads);
        writer.open();
        writer.mergeFilePair(file1.first.c_str(), file1.second.c_str(), file2.first.c_str(), file2.second.c_str());
        // remove split
        remove(file1.first.c_str()); remove(file1.second.c_str());
        remove(file2.first.c_str()); remove(file2.second.c_str());
        writer.close();
        // push back the current merge to result to the end
        files.push_back(out);
        mergeStep++;
        if(files.size() == 1 )
            break;
    }
    std::pair<std::string, std::string> out = files.front();
    Debug(Debug::WARNING) << out.first  << " " << out.second << "\n";

    std::rename(out.first.c_str(),  outDB.c_str());
    std::rename(out.second.c_str(), outDBIndex.c_str());
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nTime for mergeing results: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

}

void Prefiltering::run(int mpi_rank, int mpi_num_procs) {

    std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(outDB, outDBIndex, mpi_rank);
    if(splitMode == Parameters::TARGET_DB_SPLIT){
        maxResListLen = (maxResListLen / mpi_num_procs) + 1;
    }
    this->run(mpi_rank, mpi_num_procs, splitMode, filenamePair.first.c_str(),
              filenamePair.second.c_str());
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(mpi_rank == 0){ // master reduces results
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for(int procs = 0; procs < mpi_num_procs; procs++){
            splitFiles.push_back(Util::createTmpFileNames(outDB, outDBIndex, procs));
        }
        // merge output ffindex databases
        mergeFiles(splitFiles, splitMode);
        // remove temp databases
        // this->removeDatabaes(splitFiles);
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
                                                                 bool aaBiasCorrection, bool diagonalScoring,
                                                                 unsigned int maxSeqLen, int searchMode,
                                                                 size_t maxHitsPerQuery) {
    QueryTemplateMatcher** matchers = new QueryTemplateMatcher *[threads];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < this->threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        if (searchMode == Parameters::SEARCH_LOCAL) {
            matchers[thread_idx] = new QueryTemplateMatcherLocal(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                                 kmerSize, effectiveKmerSize, dbSize, maxSeqLen,
                                                                 maxHitsPerQuery);
        } else if(searchMode == Parameters::SEARCH_LOCAL_FAST) {
            matchers[thread_idx] = new QueryTemplateLocalFast(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                              kmerSize, dbSize, maxSeqLen, effectiveKmerSize,
                                                              maxHitsPerQuery, aaBiasCorrection, diagonalScoring, minDiagScoreThr);
        } else {
            Debug(Debug::ERROR) << "Seach mode is not valid.\n";
            EXIT(EXIT_FAILURE);
        }
        if(querySeqType == Sequence::HMM_PROFILE){
            matchers[thread_idx]->setProfileMatrix(qseq[thread_idx]->profile_matrix);
        } else {
            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix->scoreMatrix, _2merSubMatrix->scoreMatrix );
        }
    }
    return matchers;
}

IndexTable * Prefiltering::getIndexTable(int split, size_t dbFrom, size_t dbSize) {
    if(templateDBIsIndex == true ){
        return PrefilteringIndexReader::generateIndexTable(tidxdbr, split, diagonalScoring);
    }else{
        Sequence tseq(maxSeqLen, subMat->aa2int, subMat->int2aa, targetSeqType, kmerSize, spacedKmer);
        return generateIndexTable(tdbr, &tseq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize, searchMode, diagonalScoring);
    }
}

void Prefiltering::run(size_t split, size_t splitCount, int splitMode, std::string resultDB, std::string resultDBIndex) {

    Debug(Debug::WARNING) << "Process prefiltering step " << split << " of " << splitCount <<  "\n\n";

    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), threads);
    tmpDbw.open();

    size_t dbFrom = 0;
    size_t dbSize =  tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();
    IndexTable * indexTable = NULL;
    if(splitMode == Parameters::TARGET_DB_SPLIT){
        Util::decomposeDomainByAminoaAcid(tdbr->getAminoAcidDBSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                          split, splitCount, &dbFrom, &dbSize);
        indexTable = getIndexTable(split, dbFrom, dbSize);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        Util::decomposeDomainByAminoaAcid(qdbr->getAminoAcidDBSize(), qdbr->getSeqLens(), qdbr->getSize(),
                                          split, splitCount, &queryFrom, &querySize);
        indexTable = getIndexTable(0, dbFrom, dbSize); // create the whole index table
    } else{
        Debug(Debug::ERROR) << "Wrong split mode. This should not happen. Please contact developer.\n";
        EXIT(EXIT_FAILURE);
    }
    // create index table based on split parameter
    // run small query sample against the index table to calibrate p-match
    std::pair<short, double> calibration;


    const int kmerScore = getKmerThreshold(this->sensitivity, this->kmerScore);
    if(diagonalScoring == true){
        calibration = std::pair<short, double>(kmerScore, 0.0);
    }else{
        calibration = setKmerThreshold(indexTable, qdbr, tdbr, kmerScore);
    }
    //std::pair<short, double> ret = std::pair<short, double>(105, 8.18064e-05);
    this->kmerThr = calibration.first;
    this->kmerMatchProb = calibration.second;

    Debug(Debug::WARNING) << "k-mer similarity threshold: " << kmerThr << "\n";
    Debug(Debug::WARNING) << "k-mer match probability: " << kmerMatchProb << "\n\n";

    struct timeval start, end;
    gettimeofday(&start, NULL);
    QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable,
                                                                  tdbr->getSeqLens() + dbFrom, // offset for split mode
                                                                  kmerThr, kmerMatchProb, kmerSize,
                                                                  qseq[0]->getEffectiveKmerSize(), dbSize,
                                                                  aaBiasCorrection, diagonalScoring,
                                                                  maxSeqLen, searchMode, maxResListLen);
    size_t kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    size_t totalQueryDBSize = querySize;
    memset(notEmpty, 0, qdbr->getSize() * sizeof(int)); // init notEmpty

    Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step "<< split << " of " << splitCount << ")\n";
    Debug(Debug::WARNING) << "Query db start  "<< queryFrom << " to " << queryFrom + querySize << "\n";
    Debug(Debug::WARNING) << "Target db start  "<< dbFrom << " to " << dbFrom + dbSize << "\n";

#pragma omp parallel for schedule(dynamic, 100) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
    for (size_t id = queryFrom; id < queryFrom + querySize; id++){
        Log::printProgress(id);

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get query sequence
        char* seqData = qdbr->getData(id);
        unsigned int qKey = qdbr->getDbKey(id);
        qseq[thread_idx]->mapSequence(id, qKey, seqData);
        // only the corresponding split should include the id (hack for the hack)
        unsigned int targetSeqId = UINT_MAX;
        if(id >= dbFrom && id < (dbFrom + dbSize) && sameQTDB){
            targetSeqId = tdbr->getId(qseq[thread_idx]->getDbKey());
            if(targetSeqId != UINT_MAX){
                targetSeqId = targetSeqId - dbFrom;
            }
        }
        // calculate prefiltering results
        std::pair<hit_t *, size_t> prefResults = matchers[thread_idx]->matchQuery(qseq[thread_idx], targetSeqId);
        const size_t resultSize = prefResults.second;
        // write
        if(writePrefilterOutput(&tmpDbw, thread_idx, id, prefResults, dbFrom, diagonalScoring) != 0)
            continue; // could not write result because of too much results

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

    statistics_t stats(kmersPerPos/ totalQueryDBSize, dbMatches/ totalQueryDBSize, doubleMatches/ totalQueryDBSize,
                       querySeqLenSum, diagonalOverflow, resSize/ totalQueryDBSize);
    size_t empty = 0;
    for (size_t id = queryFrom; id < queryFrom + querySize; id++){
        if (notEmpty[id] == 0){
            empty++;
        }
    }
    this->printStatistics(stats, empty);

    if (totalQueryDBSize > 1000)
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
    // sorts this datafile according to the index file
    if(splitMode == Parameters::TARGET_DB_SPLIT) {
        if (splitCount > 1 && splitMode == Parameters::TARGET_DB_SPLIT) {
            DBReader<unsigned int> tmpDbr(tmpDbw.getDataFileName(), tmpDbw.getIndexFileName());
            tmpDbr.open(DBReader<unsigned int>::NOSORT);
            DBWriter tmpDbw2((resultDB + "_tmp").c_str(), (resultDBIndex + "_tmp").c_str(), threads);
            tmpDbw2.open();
            tmpDbw2.sortDatafileByIdOrder(&tmpDbr);
            tmpDbr.close();
            tmpDbw2.close();
            remove(resultDB.c_str());
            remove(resultDBIndex.c_str());
            std::rename((resultDB + "_tmp").c_str(), resultDB.c_str());
            std::rename((resultDBIndex + "_tmp").c_str(), resultDBIndex.c_str());
        }
    }
}

void Prefiltering::closeReader(){
    qdbr->close();
    //if (strcmp(qdbr->getIndexFileName(), tdbr->getIndexFileName()) != 0)
    tdbr->close();
    if(templateDBIsIndex)
        tidxdbr->close();
}

// write prefiltering to ffindex database
int Prefiltering::writePrefilterOutput(DBWriter *dbWriter, int thread_idx, size_t id,
                                       std::pair<hit_t *, size_t> prefResults,
                                       size_t seqIdOffset, bool diagonalScoring) {
    // write prefiltering results to a string
    size_t l = 0;
    hit_t * resultVector = prefResults.first;
    const size_t resultSize = prefResults.second;
    std::string prefResultsOutString;
    prefResultsOutString.reserve(BUFFER_SIZE);
    char buffer [100];

    for (size_t i = 0; i < resultSize; i++){
        hit_t * res = resultVector + i;
        size_t targetSeqId = res->seqId + seqIdOffset;
        if (targetSeqId >= tdbr->getSize()) {
            Debug(Debug::INFO) << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id)<< " -> " << targetSeqId << "\t" << res->prefScore << "\n";
        }
        int len;
        if(diagonalScoring == true){
            len = snprintf(buffer, 100, "%s\t%d\t%d\n", SSTR(tdbr->getDbKey(targetSeqId)).c_str(),
                           res->prefScore, res->diagonal);
        }else {
            len = snprintf(buffer, 100, "%s\t%.4f\t%d\n", SSTR(tdbr->getDbKey(targetSeqId)).c_str(),
                           res->pScore, res->prefScore);
        }
        prefResultsOutString.append( buffer, len );
        l++;
        // maximum allowed result list length is reached
        if (l >= this->maxResListLen)
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
    dbWriter->write(prefResultsOutData, prefResultsLength, SSTR(qdbr->getDbKey(id)).c_str(), thread_idx);
    return 0;

}


void Prefiltering::printStatistics(statistics_t &stats, size_t empty) {
    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }

//    size_t prefRealPassedPerSeq = realResSize/queryDBSize;
    Debug(Debug::INFO) << "\n" << stats.kmersPerPos << " k-mers per position.\n";
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

BaseMatrix * Prefiltering:: getSubstitutionMatrix(std::string scoringMatrixFile, int alphabetSize, float bitFactor,
                                                  bool ignoreX) {
    Debug(Debug::INFO) << "Substitution matrices...\n";
    BaseMatrix* subMat;
    if (alphabetSize < 21){
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, alphabetSize, bitFactor);
    }else if(ignoreX == true){
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2);
        subMat = new SubstitutionMatrixWithoutX(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.subMatrix, bitFactor);
        subMat->print(subMat->subMatrix, subMat->int2aa, subMat->alphabetSize);
    }else{
        subMat = new SubstitutionMatrix(scoringMatrixFile.c_str(), bitFactor, -0.2);
    }
    return subMat;
}


void Prefiltering::fillDatabase(DBReader<unsigned int>* dbr, Sequence* seq, IndexTable * indexTable,
                                size_t dbFrom, size_t dbTo)
{

    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    size_t aaCount = 0;
    dbTo=std::min(dbTo,dbr->getSize());
    for (unsigned int id = dbFrom; id < dbTo; id++){
        Log::printProgress(id - dbFrom);
        char* seqData = dbr->getData(id);
        unsigned int qKey = dbr->getDbKey(id);
        seq->mapSequence(id - dbFrom, qKey, seqData);
        aaCount += seq->L;
        indexTable->addKmerCount(seq);
    }
    Debug(Debug::INFO) << "\n";

    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";
    size_t tableEntriesNum = 0;
    for(size_t i = 0; i < indexTable->getTableSize(); i++){
        tableEntriesNum += (size_t) indexTable->getTable(i);
    }
    indexTable->initMemory(dbTo - dbFrom, tableEntriesNum, aaCount);
    indexTable->init();

    Debug(Debug::INFO) << "Index table: fill...\n";
    for (unsigned int id = dbFrom; id < dbTo; id++){
        Log::printProgress(id - dbFrom);
        char* seqData = dbr->getData(id);
        //TODO - dbFrom?!?
        unsigned int qKey = dbr->getDbKey(id);
        seq->mapSequence(id - dbFrom, qKey, seqData);
        indexTable->addSequence(seq);
    }

    if ((dbTo-dbFrom) > 10000)
        Debug(Debug::INFO) << "\n";

    Debug(Debug::INFO) << "Index table: removing duplicate entries...\n";
    indexTable->revertPointer();
    Debug(Debug::INFO) << "Index table init done.\n\n";

}

IndexTable * Prefiltering::generateIndexTable(DBReader<unsigned int>*dbr, Sequence *seq, int alphabetSize, int kmerSize,
                                              size_t dbFrom, size_t dbTo, int searchMode, bool diagonalScoring) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    IndexTable * indexTable;
    if (searchMode == Parameters::SEARCH_LOCAL || searchMode == Parameters::SEARCH_LOCAL_FAST) {
        indexTable = new IndexTable(alphabetSize, kmerSize, diagonalScoring);
    } else{
        Debug(Debug::ERROR) << "Seach mode is not valid.\n";
        EXIT(EXIT_FAILURE);
    }


    fillDatabase(dbr, seq, indexTable, dbFrom, dbTo);

    gettimeofday(&end, NULL);

    indexTable->printStatisitic(seq->int2aa);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for index table init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    return indexTable;
}

std::pair<short, double> Prefiltering::setKmerThreshold(IndexTable *indexTable, DBReader<unsigned int> *qdbr,
                                                        DBReader<unsigned int> *tdbr, const int kmerScore) {
    size_t targetDbSize = indexTable->getSize();
    size_t targetSeqLenSum = 0;
    for (size_t i = 0; i < targetDbSize; i++){
        targetSeqLenSum += (tdbr->getSeqLens(i) - qseq[0]->getEffectiveKmerSize());
    }
    // generate a small random sequence set for testing
    size_t querySetSize = std::min ( qdbr->getSize(), (size_t) 1000);
    unsigned int * querySeqs = new unsigned int[querySetSize];
    srand(1);
    for (size_t i = 0; i < querySetSize; i++){
        querySeqs[i] = rand() % qdbr->getSize();
    }

    // do a binary search through the k-mer list length threshold space to adjust the k-mer list length threshold in order to get a match probability
    // for a list of k-mers at one query position as close as possible to targetKmerMatchProb

    size_t dbMatchesExp_pc;
    // 1000 * 350 * 100000 * 350
    size_t lenSum_pc = 12250000000000;

    double kmerMatchProb;


    statistics_t stats;
    stats = computeStatisticForKmerThreshold(indexTable, querySetSize, querySeqs, searchMode, kmerScore);
    // match probability with pseudocounts
    // add pseudo-counts (1/20^6 * kmerPerPos * length of pseudo counts)
    dbMatchesExp_pc = (size_t)(((double)lenSum_pc) * stats.kmersPerPos * pow((1.0/((double)(subMat->alphabetSize-1))), kmerSize));
    // match probability with pseudocounts
    kmerMatchProb = ((double)stats.dbMatches + dbMatchesExp_pc) / ((double) (stats.querySeqLen * targetSeqLenSum + lenSum_pc));
    // compute match prob for local match
    if(searchMode == Parameters::SEARCH_LOCAL || searchMode == Parameters::SEARCH_LOCAL_FAST){
        kmerMatchProb = ((double) stats.doubleMatches) / ((double) (stats.querySeqLen * targetSeqLenSum));
        kmerMatchProb /= 256;
    }
    kmerMatchProb = std::max(kmerMatchProb, std::numeric_limits<double>::min());
    Debug(Debug::INFO) << "\tk-mers per position = " << stats.kmersPerPos << ", k-mer match probability: " << kmerMatchProb << "\n";

    delete[] querySeqs;

    return std::pair<short, double> (kmerScore, kmerMatchProb);
}

statistics_t Prefiltering::computeStatisticForKmerThreshold(IndexTable *indexTable, size_t querySetSize,
                                                            unsigned int *querySeqsIds, bool reverseQuery, size_t kmerThrMid) {
    // determine k-mer match probability for kmerThrMid
    QueryTemplateMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable, tdbr->getSeqLens(), kmerThrMid,
                                                                  1.0, kmerSize, qseq[0]->getEffectiveKmerSize(),
                                                                  indexTable->getSize(), aaBiasCorrection, false, maxSeqLen,
                                                                  searchMode, LONG_MAX);
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
        unsigned int qKey = qdbr->getDbKey(id);
        qseq[thread_idx]->mapSequence(id, qKey, seqData);
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

void Prefiltering::mergeFiles(std::vector<std::pair<std::string, std::string>> splitFiles, int mode) {
    if(mode == Parameters::TARGET_DB_SPLIT){
        this->mergeOutput(splitFiles);
    }else if (mode == Parameters::QUERY_DB_SPLIT){
        const char * datafilesNames[splitFiles.size()];
        const char * indexFilesNames[splitFiles.size()];
        for(size_t i = 0; i < splitFiles.size(); i++){
            datafilesNames[i] = splitFiles[i].first.c_str();
            indexFilesNames[i] = splitFiles[i].second.c_str();
        }
        DBWriter::mergeFFindexFile(outDB.c_str(), outDBIndex.c_str(), "w", datafilesNames, indexFilesNames, splitFiles.size() );
    }
}

const int Prefiltering::getKmerThreshold(const float sensitivity, const int score) {
    const unsigned int sens =  sensitivity;

    int kmerThrBest = kmerScore;
    if(kmerThrBest == INT_MAX){
        if (kmerSize == 5){
            kmerThrBest = 123.75 - (sens * 8.75);
        } else if (kmerSize == 6){
            kmerThrBest = 133.75 - (sens * 8.75);
        } else if (kmerSize == 7){
            kmerThrBest = 154.75 - (sens * 9.75);
        }
        else{
            Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return kmerThrBest;
}
