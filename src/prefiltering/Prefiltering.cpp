#include "Prefiltering.h"
#include "PrefilteringIndexReader.h"
#include "QueryMatcher.h"
#include "NucleotideMatrix.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrixWithoutX.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"

#include <regex.h>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

Prefiltering::Prefiltering(const std::string& queryDB,
                           const std::string& queryDBIndex,
                           const std::string& targetDB,
                           const std::string& targetDBIndex,
                           const std::string& outDB,
                           const std::string& outDBIndex,
                           Parameters &par):
        outDB(outDB),
        outDBIndex(outDBIndex),
        kmerSize(par.kmerSize),
        kmerScore(par.kmerScore),
        spacedKmer(par.spacedKmer),
        sensitivity(par.sensitivity),
        resListOffset(par.resListOffset),
        maxResListLen(par.maxResListLen),
        alphabetSize(par.alphabetSize),
        maxSeqLen(par.maxSeqLen),
        querySeqType(par.querySeqType),
        targetSeqType(par.targetSeqType),
        diagonalScoring(par.diagonalScoring),
        minDiagScoreThr(par.minDiagScoreThr),
        aaBiasCorrection(par.compBiasCorrection),
        split(par.split),
        splitMode(par.splitMode)
{
    this->threads = par.threads;
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#else
    this->threads = 1;
#endif
    Debug(Debug::INFO) << "\n";
//    FileUtil::errorIfFileExist(outDB.c_str());
//    FileUtil::errorIfFileExist(outDBIndex.c_str());
    this->qdbr = new DBReader<unsigned int>(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    //  check if when qdb and tdb have the same name an index extention exists
    std::string check(targetDB);
    size_t pos = check.find(queryDB);
    int nomatch = true;
    if(pos == 0){
        check.replace(0, queryDB.length(), "");
        regex_t regex;
        regcomp(&regex, "^\\.s?k[5-7]$", REG_EXTENDED | REG_NEWLINE);
        nomatch = regexec(&regex, check.c_str(), 0, NULL, 0);
        regfree(&regex);
    }
    // if no match found or two matches found (we want exactly one match)
    sameQTDB = (queryDB.compare(targetDB) == 0 || (nomatch == false) );
    includeIdentical = par.includeIdentity;
    std::string indexDB = searchForIndex(targetDB);

    //TODO optimize this. Dont read twice the target index. This seems to be slow
    if(indexDB != ""){
        Debug(Debug::INFO) << "Use index  " << indexDB << "\n";

        this->tdbr = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str());
    }else{
        Debug(Debug::INFO) << "Cound not find precomputed index. Compute index.\n";
        this->tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str());
    }
    tdbr->open(DBReader<unsigned int>::NOSORT);
    templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tdbr);
    if(templateDBIsIndex == true){ // exchange reader with old ffindex reader
        this->tidxdbr = this->tdbr;
        this->tdbr = PrefilteringIndexReader::openNewReader(tdbr);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
        this->kmerSize     = data.kmerSize;
        this->alphabetSize = data.alphabetSize;
        this->split        = data.split;
        this->spacedKmer   = (data.spacedKmer == 1) ? true : false;
    }
    if(templateDBIsIndex == false && sameQTDB == true){
        this->tdbr->close();
        delete tdbr;
        this->tdbr = qdbr;
    }else if(templateDBIsIndex == false){
        this->tdbr->close();
        delete tdbr;
        this->tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str());
        tdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }

    Debug(Debug::INFO) << "Query database: " << par.db1 << "(size=" << qdbr->getSize() << ")\n";
    Debug(Debug::INFO) << "Target database: " << par.db2 << "(size=" << tdbr->getSize() << ")\n";

    const size_t totalMemoryInByte =  Util::getTotalSystemMemory();
    size_t neededSize = estimateMemoryConsumption(1,
                                                  tdbr->getSize(), tdbr->getAminoAcidDBSize(), alphabetSize,
                                                  kmerSize == 0 ? // if auto detect kmerSize
                                                  IndexTable::computeKmerSize(tdbr->getAminoAcidDBSize()) : kmerSize,
                                                  threads);
    if(neededSize > 0.9 * totalMemoryInByte){ // memory is not enough to compute everything at once
        std::pair<int, int> splitingSetting = optimizeSplit(totalMemoryInByte, tdbr, alphabetSize, kmerSize, threads);
        if(splitingSetting.second == -1){
            Debug(Debug::ERROR) << "Can not fit databased into " << totalMemoryInByte <<" byte. Please use a computer with more main memory.\n";
            EXIT(EXIT_FAILURE);
        }
        if(par.split == Parameters::AUTO_SPLIT_DETECTION && templateDBIsIndex == false){
            this->split = splitingSetting.second;
            if(kmerSize == 0){ // set k-mer based on aa size in database
                // if we have less than 10Mio * 335 amino acids use 6mers
                kmerSize = splitingSetting.first;
            }
        }
        if(splitMode == Parameters::DETECT_BEST_DB_SPLIT) {
            splitMode = Parameters::TARGET_DB_SPLIT;
        }
    } else { // memory is  enough to compute everything with split setting
        if(kmerSize == 0){
            const int tmpSplit = (par.split > 1) ? par.split : 1;
            size_t aaSize = tdbr->getAminoAcidDBSize() / tmpSplit;
            kmerSize = IndexTable::computeKmerSize(aaSize);
        }

        if(this->split == Parameters::AUTO_SPLIT_DETECTION)
            this->split = 1;
        if(splitMode == Parameters::DETECT_BEST_DB_SPLIT){
            if(templateDBIsIndex == true && this->split > 1){
                splitMode = Parameters::TARGET_DB_SPLIT;
            }else{
#ifdef HAVE_MPI
                splitMode = Parameters::QUERY_DB_SPLIT;
#else
                splitMode = Parameters::TARGET_DB_SPLIT;
#endif
            }
        }
    }
    Debug(Debug::INFO) << "Use kmer size " << kmerSize << " and split " << this->split << " using split mode " << this->splitMode <<"\n";
    neededSize = estimateMemoryConsumption((splitMode == Parameters::TARGET_DB_SPLIT) ? split : 1, tdbr->getSize(), tdbr->getAminoAcidDBSize(), alphabetSize, kmerSize,
                                                        threads);
    //Debug(Debug::INFO) << "Split target databases into " << split << " parts because of memory constraint.\n";
    Debug(Debug::INFO) << "Needed memory (" << neededSize << " byte) of total memory (" << totalMemoryInByte << " byte)\n";
    if(neededSize > 0.9 * totalMemoryInByte){
        Debug(Debug::WARNING) << "WARNING: MMseqs processes needs more main memory than available."
                                 "Increase the size of --split or set it to 0 to automatic optimize target database split.\n";
        if(templateDBIsIndex == true){
            Debug(Debug::WARNING) << "WARNING: Split has to be computed by createindex if precomputed index is used.\n";
        }
    }


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
    this->notEmpty = new char[qdbr->getSize()];

#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        qseq[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, querySeqType, kmerSize, spacedKmer, aaBiasCorrection);
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

    if(sameQTDB == false){
        delete tdbr;
    }
    if(templateDBIsIndex == true){
        delete tidxdbr;
    }
    delete subMat;
    if(_2merSubMatrix != NULL)
        delete _2merSubMatrix;
    if(_3merSubMatrix != NULL)
        delete _3merSubMatrix;
}

void Prefiltering::run(size_t fromSplit, size_t splits){
    // splits template database into x sequence steps
    unsigned int splitCounter = this->split;
    if(splitCounter > 1 && splitMode == Parameters::TARGET_DB_SPLIT) {
        maxResListLen = (maxResListLen / splitCounter)  + 1;
    }
    std::vector<std::pair<std::string, std::string> > splitFiles;
    for(unsigned int split = fromSplit; split < (fromSplit+splits); split++){
        std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(outDB,outDBIndex,split);
        splitFiles.push_back(filenamePair);
        this->run(split, splitCounter, splitMode, filenamePair.first.c_str(),
                  filenamePair.second.c_str());
    } // prefiltering scores calculation end
    // merge output ffindex databases
    if(splitCounter > 1){
        mergeFiles(splitFiles, splitMode, outDB, outDBIndex);
    }else{
        std::rename(splitFiles[0].first.c_str(),  outDB.c_str());
        std::rename(splitFiles[0].second.c_str(), outDBIndex.c_str());
    }
    // close reader to reduce memory
    this->closeReader();
}

void Prefiltering::mergeOutput(const std::string& outDB, const std::string& outDBIndex,
                               const std::vector<std::pair<std::string, std::string>>& filenames){
    struct timeval start, end;
    gettimeofday(&start, NULL);
    if(filenames.size() < 2){
        std::rename(filenames[0].first.c_str(),  outDB.c_str());
        std::rename(filenames[0].second.c_str(), outDBIndex.c_str());
        Debug(Debug::INFO) << "No merging needed.\n";
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
        DBWriter writer(out.first.c_str(), out.second.c_str(), 1);
        writer.open(1024 * 1024 * 1024); // 1 GB buffer
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
    Debug(Debug::WARNING) << "\nTime for merging results: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
}

QueryMatcher ** Prefiltering::createQueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable,
                                                         unsigned int *seqLens, short kmerThr,
                                                         double kmerMatchProb, int kmerSize,
                                                         size_t effectiveKmerSize, size_t dbSize,
                                                         bool aaBiasCorrection, bool diagonalScoring,
                                                         unsigned int maxSeqLen, size_t maxHitsPerQuery) {
    QueryMatcher** matchers = new QueryMatcher *[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < this->threads; i++){
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        matchers[thread_idx] = new QueryMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                kmerSize, dbSize, maxSeqLen, effectiveKmerSize,
                                                maxHitsPerQuery, aaBiasCorrection, diagonalScoring, minDiagScoreThr);
        if(querySeqType == Sequence::HMM_PROFILE){
            matchers[thread_idx]->setProfileMatrix(qseq[thread_idx]->profile_matrix);
        } else {
            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix->scoreMatrix, _2merSubMatrix->scoreMatrix );
        }
    }
    return matchers;
}

IndexTable * Prefiltering::getIndexTable(int split, size_t dbFrom, size_t dbSize, int threads) {
    if(templateDBIsIndex == true ){
        return PrefilteringIndexReader::generateIndexTable(tidxdbr, split, diagonalScoring);
    }else{
        Sequence tseq(maxSeqLen, subMat->aa2int, subMat->int2aa, targetSeqType, kmerSize, spacedKmer, aaBiasCorrection);
        return generateIndexTable(tdbr, &tseq, subMat, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize, diagonalScoring, threads);
    }
}

void Prefiltering::run(size_t split, size_t splitCount, int splitMode, const std::string& resultDB, const std::string& resultDBIndex) {

    Debug(Debug::WARNING) << "Process prefiltering step " << split << " of " << splitCount <<  "\n\n";
    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), threads);
    tmpDbw.open();
    size_t dbFrom = 0;
    size_t dbSize =  tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();
    IndexTable * indexTable = NULL;
    if(splitMode == Parameters::TARGET_DB_SPLIT){
        Util::decomposeDomainByAminoAcid(tdbr->getAminoAcidDBSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                         split, splitCount, &dbFrom, &dbSize);
        //TODO fix this what if we have 10 chunks but only 4 servers (please fix me)
        indexTable = getIndexTable(split, dbFrom, dbSize, threads);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        Util::decomposeDomainByAminoAcid(qdbr->getAminoAcidDBSize(), qdbr->getSeqLens(), qdbr->getSize(),
                                         split, splitCount, &queryFrom, &querySize);
        indexTable = getIndexTable(0, dbFrom, dbSize, threads); // create the whole index table
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
    QueryMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable,
                                                          tdbr->getSeqLens() + dbFrom, // offset for split mode
                                                          kmerThr, kmerMatchProb, kmerSize,
                                                          qseq[0]->getEffectiveKmerSize(), dbSize,
                                                          aaBiasCorrection, diagonalScoring,
                                                          maxSeqLen, maxResListLen);
    size_t kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    size_t totalQueryDBSize = querySize;
    memset(notEmpty, 0, qdbr->getSize() * sizeof(char)); // init notEmpty

    Debug(Debug::WARNING) << "Starting prefiltering scores calculation (step "<< split << " of " << splitCount << ")\n";
    Debug(Debug::WARNING) << "Query db start  "<< queryFrom << " to " << queryFrom + querySize << "\n";
    Debug(Debug::WARNING) << "Target db start  "<< dbFrom << " to " << dbFrom + dbSize << "\n";

#pragma omp parallel for schedule(dynamic, 10) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
    for (size_t id = queryFrom; id < queryFrom + querySize; id++){
        Debug::printProgress(id);

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
        if(id >= dbFrom && id < (dbFrom + dbSize) && (sameQTDB || includeIdentical) ){
            targetSeqId = tdbr->getId(qseq[thread_idx]->getDbKey());
            if(targetSeqId != UINT_MAX){
                targetSeqId = targetSeqId - dbFrom;
            }
        }
        // calculate prefiltering results
        std::pair<hit_t *, size_t> prefResults = matchers[thread_idx]->matchQuery(qseq[thread_idx], targetSeqId);
        const size_t resultSize = prefResults.second;
        // write
        if(writePrefilterOutput(&tmpDbw, thread_idx, id, prefResults, dbFrom, diagonalScoring, resListOffset) != 0)
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

    // sort by ids
    // needed to speed up merge later one
    // sorts this datafile according to the index file
    if(splitMode == Parameters::TARGET_DB_SPLIT) {
        if (splitCount > 1 && splitMode == Parameters::TARGET_DB_SPLIT) {
            DBReader<unsigned int> resultReader(tmpDbw.getDataFileName(), tmpDbw.getIndexFileName());
            resultReader.open(DBReader<unsigned int>::NOSORT);
            DBWriter resultWriter((resultDB + "_tmp").c_str(), (resultDBIndex + "_tmp").c_str(), threads);
            resultWriter.open();
            resultWriter.sortDatafileByIdOrder(resultReader);
            resultReader.close();
            resultWriter.close();
            remove(resultDB.c_str());
            remove(resultDBIndex.c_str());
            std::rename((resultDB + "_tmp").c_str(), resultDB.c_str());
            std::rename((resultDBIndex + "_tmp").c_str(), resultDBIndex.c_str());
        }
    }
}

void Prefiltering::closeReader(){
    qdbr->close();
    if(sameQTDB == false || templateDBIsIndex == true){
        tdbr->close();
    }
    if(templateDBIsIndex){
        tidxdbr->close();
    }
}

// write prefiltering to ffindex database
int Prefiltering::writePrefilterOutput(DBWriter *dbWriter, int thread_idx, size_t id,
                                       const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                                       bool diagonalScoring, size_t resultOffsetPos) {
    // write prefiltering results to a string
    size_t l = 0;
    hit_t * resultVector = prefResults.first + resultOffsetPos;
    const size_t resultSize = (prefResults.second < resultOffsetPos) ? 0 : prefResults.second - resultOffsetPos;
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
    char* prefResultsOutData = (char *) prefResultsOutString.c_str();
    dbWriter->writeData(prefResultsOutData, prefResultsLength, SSTR(qdbr->getDbKey(id)).c_str(), thread_idx);
    return 0;
}

void Prefiltering::printStatistics(const statistics_t &stats, size_t empty) {
    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (int i = 1; i < threads; i++){
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }
    Debug(Debug::INFO) << "\n" << stats.kmersPerPos << " k-mers per position.\n";
    Debug(Debug::INFO) << stats.dbMatches << " DB matches per sequence.\n";
    Debug(Debug::INFO) << stats.diagonalOverflow << " Overflows .\n";
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

BaseMatrix * Prefiltering:: getSubstitutionMatrix(const std::string& scoringMatrixFile, int alphabetSize, float bitFactor,
                                                  bool ignoreX) {
    Debug(Debug::INFO) << "Substitution matrices...\n";
    BaseMatrix* subMat;
    if (alphabetSize < 21){
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, alphabetSize, bitFactor);
    }else if(ignoreX == true){
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2);
        subMat = new SubstitutionMatrixWithoutX(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.subMatrix, bitFactor);
    }else{
        subMat = new SubstitutionMatrix(scoringMatrixFile.c_str(), bitFactor, -0.2);
    }
    return subMat;
}


void Prefiltering::fillDatabase(DBReader<unsigned int>* dbr, Sequence* seq, IndexTable * indexTable,
                                BaseMatrix *subMat, size_t dbFrom, size_t dbTo, bool diagonalScoring, int threads)
{
    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    size_t aaCount = 0;
    dbTo=std::min(dbTo,dbr->getSize());
    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    size_t tableSize = 0;

    size_t dbSize = dbTo - dbFrom;
    size_t * sequenceOffSet = new size_t[dbSize];
    size_t aaDbSize = 0;
    sequenceOffSet[0] = 0;
    for (unsigned int id = dbFrom; id < dbTo; id++){
        int seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 2, 0);
        aaDbSize += seqLen; // remove /n and /0
        size_t idFromNull = (id - dbFrom);
        if(id < dbTo - 1){
            sequenceOffSet[idFromNull + 1] =  sequenceOffSet[idFromNull] + seqLen;
        }
        if(Util::overlappingKmers(seqLen, seq->getEffectiveKmerSize() > 0)){
            tableSize += 1;
        }
    }
    SequenceLookup * sequenceLookup = new SequenceLookup(dbSize, aaDbSize);
#pragma omp parallel
    {
        Indexer idxer(subMat->alphabetSize, seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
#pragma omp for schedule(dynamic, 100) reduction(+:aaCount, totalKmerCount, maskedResidues)
        for (unsigned int id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);
            char *seqData = dbr->getData(id);
            unsigned int qKey = dbr->getDbKey(id);
            s.mapSequence(id - dbFrom, qKey, seqData);

            maskedResidues += Util::maskLowComplexity(subMat, &s, s.L, 12, 3,
                                                      indexTable->getAlphabetSize(), seq->aa2int[(unsigned char) 'X']);

            aaCount += s.L;
            totalKmerCount += indexTable->addKmerCount(&s, &idxer);
            sequenceLookup->addSequence(&s, sequenceOffSet[id-dbFrom]);
        }
    }
    delete [] sequenceOffSet;
    dbr->remapData();
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: Masked residues: " << maskedResidues << "\n";
    //TODO find smart way to remove extrem k-mers without harming huge protein families
//    size_t lowSelectiveResidues = 0;
//    const float dbSize = static_cast<float>(dbTo - dbFrom);
//    for(size_t kmerIdx = 0; kmerIdx < indexTable->getTableSize(); kmerIdx++){
//        size_t res = (size_t) indexTable->getTable(kmerIdx);
//        float selectivityOfKmer = (static_cast<float>(res)/dbSize);
//        if(selectivityOfKmer > 0.005){
//            indexTable->getTable()[kmerIdx] = 0;
//            lowSelectiveResidues += res;
//        }
//    }
//    Debug(Debug::INFO) << "Index table: Remove "<< lowSelectiveResidues <<" none selective residues\n";
//    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";
    size_t tableEntriesNum = 0;
    for(size_t i = 0; i < indexTable->getTableSize(); i++){
        tableEntriesNum += (size_t) indexTable->getTable(i);
    }

#pragma omp parallel
    {
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
        Indexer idxer(subMat->alphabetSize, seq->getKmerSize());
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        size_t threadFrom, threadSize;
        char **table = (indexTable->getTable());
        Util::decomposeDomainSizet(tableEntriesNum, (size_t *) table, indexTable->getTableSize(),
                                   thread_idx, threads, &threadFrom, &threadSize);
//        std::stringstream stream;
//        stream << thread_idx << "\t" << threadFrom << "\t" << threadSize;
//        std::cout << stream.str() << std::endl;

#pragma omp barrier
        if(thread_idx == 0){
            if(diagonalScoring == true){
                indexTable->initMemory(dbTo - dbFrom, tableEntriesNum, aaCount, sequenceLookup, tableSize);
            }else{
                indexTable->initMemory(dbTo - dbFrom, tableEntriesNum, aaCount, NULL, tableSize);
            }
            indexTable->init();
            Debug(Debug::INFO) << "Index table: fill...\n";
        }
#pragma omp barrier
        for (unsigned int id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            if(thread_idx == 0) {
                Debug::printProgress(id - dbFrom);
            }
            //char *seqData = dbr->getData(id);
            //TODO - dbFrom?!?
            unsigned int qKey = dbr->getDbKey(id);
            //seq->mapSequence(id - dbFrom, qKey, seqData);
            s.mapSequence(id - dbFrom, qKey, sequenceLookup->getSequence(id-dbFrom));
//            Util::maskLowComplexity(subMat, seq, seq->L, 12, 3,
//                                    indexTable->getAlphabetSize(), seq->aa2int['X']);
            indexTable->addSequence(&s, &idxer, threadFrom, threadSize);
        }
    }
    if(diagonalScoring == false){
        delete sequenceLookup;
    }
    if ((dbTo-dbFrom) > 10000)
        Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: removing duplicate entries...\n";
    indexTable->revertPointer();
    Debug(Debug::INFO) << "Index table init done.\n\n";

}

IndexTable * Prefiltering::generateIndexTable(DBReader<unsigned int>*dbr, Sequence *seq, BaseMatrix * subMat, int alphabetSize, int kmerSize,
                                              size_t dbFrom, size_t dbTo, bool diagonalScoring, int threads) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    IndexTable * indexTable;
    indexTable = new IndexTable(alphabetSize, kmerSize);
    fillDatabase(dbr, seq, indexTable, subMat, dbFrom, dbTo, diagonalScoring, threads);
    gettimeofday(&end, NULL);
    indexTable->printStatisitic(seq->int2aa);
    int sec = end.tv_sec - start.tv_sec;
    dbr->remapData();
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

    statistics_t stats;
    stats = computeStatisticForKmerThreshold(indexTable, querySetSize, querySeqs, true, kmerScore);

    // compute match prob for local match
    double kmerMatchProb = ((double) stats.doubleMatches) / ((double) (stats.querySeqLen * targetSeqLenSum));
    kmerMatchProb /= 256;

    kmerMatchProb = std::max(kmerMatchProb, std::numeric_limits<double>::min());
    Debug(Debug::INFO) << "\tk-mers per position = " << stats.kmersPerPos << ", k-mer match probability: " << kmerMatchProb << "\n";
    delete[] querySeqs;
    return std::pair<short, double> (kmerScore, kmerMatchProb);
}

statistics_t Prefiltering::computeStatisticForKmerThreshold(IndexTable *indexTable, size_t querySetSize,
                                                            unsigned int *querySeqsIds, bool reverseQuery, size_t kmerThrMid) {
    // determine k-mer match probability for kmerThrMid
    QueryMatcher ** matchers = createQueryTemplateMatcher(subMat, indexTable, tdbr->getSeqLens(), kmerThrMid,
                                                          1.0, kmerSize, qseq[0]->getEffectiveKmerSize(),
                                                          indexTable->getSize(), aaBiasCorrection, false, maxSeqLen,
                                                          LONG_MAX);
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

void Prefiltering::mergeFiles(const std::vector<std::pair<std::string, std::string>>& splitFiles, int mode, std::string outDB, std::string outDBIndex) {
    if(mode == Parameters::TARGET_DB_SPLIT){
        mergeOutput(outDB, outDBIndex, splitFiles);
    }else if (mode == Parameters::QUERY_DB_SPLIT){
        const char ** datafilesNames = new const char *[splitFiles.size()];
        const char ** indexFilesNames = new const char *[splitFiles.size()];
        for(size_t i = 0; i < splitFiles.size(); i++){
            datafilesNames[i]  = splitFiles[i].first.c_str();
            indexFilesNames[i] = splitFiles[i].second.c_str();
        }
        DBWriter::mergeResults(outDB.c_str(), outDBIndex.c_str(), datafilesNames, indexFilesNames, splitFiles.size());
        delete [] datafilesNames;
        delete [] indexFilesNames;
    }
}

int Prefiltering::getKmerThreshold(const float sensitivity, const int score) {
    const unsigned int sens =  sensitivity;
    int kmerThrBest = kmerScore;
    if(kmerThrBest == INT_MAX){
        if (kmerSize == 5){
            kmerThrBest = 123.75 - (sens * 8.75);
        } else if (kmerSize == 6){
            kmerThrBest = 138.75 - (sens * 8.75);
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

size_t Prefiltering::estimateMemoryConsumption(int split, size_t dbSize, size_t resSize, int alphabetSize, int kmerSize,
                                               int threads) {
    // for each residue in the database we need 7 byte
    size_t dbSizeSplit = (dbSize) / split;
    size_t residueSize = (resSize/split * 7);
    // 21^7 * pointer size is needed for the index
    size_t indexTableSize   =  pow(alphabetSize, kmerSize) * sizeof(char *);
    // memory needed for the threads
    // This memory is an approx. for Countint32Array and QueryTemplateLocalFast
    size_t threadSize =  threads * (
                 (dbSizeSplit * 2 * sizeof(IndexEntryLocal)) // databaseHits in QueryMatcher
               + (dbSizeSplit * sizeof(CounterResult)) // databaseHits in QueryMatcher
               + (QueryMatcher::MAX_RES_LIST_LEN * sizeof(hit_t))
               + (dbSizeSplit * 2 * sizeof(CounterResult) * 2) // BINS * binSize, (binSize = dbSize * 2 / BINS)
                                                               // 2 is a security factor the size can increase during run
               );
    // extended matrix
    size_t extenededMatrix = sizeof(std::pair<short,unsigned int>) * pow(pow(alphabetSize,3),2);
    extenededMatrix += sizeof(std::pair<short,unsigned int>) * pow(pow(alphabetSize,2),2);
    // some memory needed to keep the index, ....
    size_t background =  dbSize * 22;
    return residueSize + indexTableSize + threadSize + background + extenededMatrix;
}

std::string Prefiltering::searchForIndex(const std::string &pathToDB) {
    for(size_t spaced = 0; spaced < 2; spaced++) {
        for (size_t k = 5; k <= 7; k++) {
            std::string outIndexName(pathToDB); // db.sk6
            std::string s = (spaced == true) ? "s" : "";
            outIndexName.append(".").append(s).append("k").append(SSTR(k));
            if (FileUtil::fileExists(outIndexName.c_str()) == true) {
                return outIndexName;
            }
        }
    }
    return "";
}

std::pair<int, int> Prefiltering::optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr,
                                                int alphabetSize, int externalKmerSize, int threads) {
    for(int split = 1; split < 100; split++ ){
        for(int kmerSize = 6; kmerSize <= 7; kmerSize++){
            if(kmerSize==externalKmerSize || externalKmerSize == 0){ // 0: set k-mer based on aa size in database
                size_t aaUpperBoundForKmerSize = IndexTable::getUpperBoundAACountForKmerSize(kmerSize);
                if((tdbr->getAminoAcidDBSize() / split) < aaUpperBoundForKmerSize){
                    size_t neededSize = estimateMemoryConsumption(split, tdbr->getSize(), tdbr->getAminoAcidDBSize(), alphabetSize,
                                                                  kmerSize, threads);
                    if(neededSize < 0.9 * totalMemoryInByte){
                        return std::make_pair(kmerSize, split);
                    }
                }
            }
        }
    }
    return std::make_pair(-1, -1);
}
