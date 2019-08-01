#include "Prefiltering.h"
#include "NucleotideMatrix.h"
#include "ReducedMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrixProfileStates.h"
#include "PatternCompiler.h"
#include "FileUtil.h"
#include "IndexBuilder.h"
#include "Timer.h"
#include "ByteParser.h"
#include "Parameters.h"

namespace prefilter {
#include "ExpOpt3_8_polished.cs32.lib.h"
}

#ifdef OPENMP
#include <omp.h>
#endif

Prefiltering::Prefiltering(const std::string &queryDB,
                           const std::string &queryDBIndex,
                           const std::string &targetDB,
                           const std::string &targetDBIndex,
                           int querySeqType, int targetSeqType_,
                           const Parameters &par) :
        queryDB(queryDB),
        queryDBIndex(queryDBIndex),
        targetDB(targetDB),
        targetDBIndex(targetDBIndex),
        _2merSubMatrix(NULL),
        _3merSubMatrix(NULL),
        splits(par.split),
        kmerSize(par.kmerSize),
        spacedKmerPattern(par.spacedKmerPattern),
        localTmp(par.localTmp),
        spacedKmer(par.spacedKmer != 0),
        alphabetSize(par.alphabetSize),
        maskMode(par.maskMode),
        maskLowerCaseMode(par.maskLowerCaseMode),
        splitMode(par.splitMode),
        scoringMatrixFile(par.scoringMatrixFile),
        seedScoringMatrixFile(par.seedScoringMatrixFile),
        targetSeqType(targetSeqType_),
        maxResListLen(par.maxResListLen),
        prevMaxResListLengths(par.prevMaxResListLengths),
        kmerScore(par.kmerScore),
        sensitivity(par.sensitivity),
        maxSeqLen(par.maxSeqLen),
        querySeqType(querySeqType),
        diagonalScoring(par.diagonalScoring),
        minDiagScoreThr(static_cast<unsigned int>(par.minDiagScoreThr)),
        aaBiasCorrection(par.compBiasCorrection != 0),
        covThr(par.covThr), covMode(par.covMode), includeIdentical(par.includeIdentity),
        preloadMode(par.preloadMode),
        threads(static_cast<unsigned int>(par.threads)), compressed(par.compressed) {
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif
    sameQTDB = isSameQTDB();

    // init the substitution matrices
    switch (querySeqType & 0x7FFFFFFF) {
        case Parameters::DBTYPE_NUCLEOTIDES:
            kmerSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 1.0, false, true);
            ungappedSubMat = kmerSubMat;
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        case Parameters::DBTYPE_AMINO_ACIDS:
            kmerSubMat = getSubstitutionMatrix(seedScoringMatrixFile, alphabetSize, 8.0, false, false);
            ungappedSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 2.0, false, false);
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        case Parameters::DBTYPE_HMM_PROFILE:
            // needed for Background distributions
            kmerSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0, false, false);
            ungappedSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 2.0, false, false);
            break;
        case Parameters::DBTYPE_PROFILE_STATE_PROFILE:
            kmerSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0, true, false);
            ungappedSubMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 2.0, false, false);
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        default:
            Debug(Debug::ERROR) << "Query sequence type not implemented!\n";
            EXIT(EXIT_FAILURE);
    }

    if (Parameters::isEqualDbtype(FileUtil::parseDbType(targetDB.c_str()), Parameters::DBTYPE_INDEX_DB)) {
        int dataMode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;
        if(preloadMode == Parameters::PRELOAD_MODE_AUTO){
            if(sensitivity > 6.0){
                preloadMode = Parameters::PRELOAD_MODE_FREAD;
            }else{
                preloadMode = Parameters::PRELOAD_MODE_MMAP_TOUCH;
            }
        }
        if (preloadMode == Parameters::PRELOAD_MODE_FREAD) {
            dataMode |= DBReader<unsigned int>::USE_FREAD;
        }
        tdbr = new DBReader<unsigned int>(targetDB.c_str(), (targetDB + ".index").c_str(), threads, dataMode);
        tdbr->open(DBReader<unsigned int>::NOSORT);

        templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tdbr);
        if (templateDBIsIndex == true) {
            // exchange reader with old reader
            tidxdbr = tdbr;
            bool touch = false;
            if (preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH) {
                touch = true;
                tidxdbr->readMmapedDataInMemory();
            }
            tdbr = PrefilteringIndexReader::openNewReader(tdbr, PrefilteringIndexReader::DBR1DATA, PrefilteringIndexReader::DBR1INDEX, false, threads, touch, touch);
            PrefilteringIndexReader::printSummary(tidxdbr);
            PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
            for(size_t i = 0; i < par.prefilter.size(); i++){
                if(par.prefilter[i]->wasSet && par.prefilter[i]->uniqid == par.PARAM_K.uniqid){
                    if(kmerSize != 0 && data.kmerSize != kmerSize){
                        Debug(Debug::WARNING) << "Index was created with -k " << data.kmerSize << " but the prefilter was called with -k " << kmerSize << "!\n";
                        Debug(Debug::WARNING) << "Search with -k " <<  data.kmerSize << "\n";
                    }
                }
                if(par.prefilter[i]->wasSet && par.prefilter[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid){
                    if(data.alphabetSize != alphabetSize){
                        Debug(Debug::WARNING) << "Index was created with --alph-size  " << data.alphabetSize << " but the prefilter was called with --alph-size " << alphabetSize << "!\n";
                        Debug(Debug::WARNING) << "Current search will use --alph-size " <<  data.alphabetSize  << "\n";
                    }
                }
                if(par.prefilter[i]->wasSet && par.prefilter[i]->uniqid == par.PARAM_SPACED_KMER_MODE.uniqid){
                    if(data.spacedKmer != spacedKmer){
                        Debug(Debug::WARNING) << "Index was created with --spaced-kmer-mode " << data.spacedKmer << " but the prefilter was called with --spaced-kmer-mode " << spacedKmer << "!\n";
                        Debug(Debug::WARNING) <<  "Current search will use  --spaced-kmer-mode " <<  data.spacedKmer << "\n";
                    }
                }
                if(par.prefilter[i]->wasSet && par.prefilter[i]->uniqid == par.PARAM_NO_COMP_BIAS_CORR.uniqid){
                    if(data.compBiasCorr != aaBiasCorrection && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)){
                        Debug(Debug::WARNING) << "Index was created with --comp-bias-corr " << data.compBiasCorr  <<" please recreate index with --comp-bias-corr " << aaBiasCorrection << "!\n";
                        Debug(Debug::WARNING) <<  "Current search will use --comp-bias-corr " <<  data.compBiasCorr << "\n";
                    }
                }
            }

            kmerSize = data.kmerSize;
            alphabetSize = data.alphabetSize;
            targetSeqType = data.seqType;
            spacedKmer   = (data.spacedKmer == 1) ? true : false;
            maxSeqLen = data.maxSeqLength;
            aaBiasCorrection = data.compBiasCorr;

            if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) &&
                Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
                Debug(Debug::ERROR) << "Query-profiles cannot be searched against a target-profile database!\n";
                EXIT(EXIT_FAILURE);
            }

            splits = 1;
            spacedKmer = data.spacedKmer != 0;
            spacedKmerPattern = PrefilteringIndexReader::getSpacedPattern(tidxdbr);
            seedScoringMatrixFile = ScoreMatrixFile(PrefilteringIndexReader::getSubstitutionMatrix(tidxdbr));
        } else {
            Debug(Debug::ERROR) << "Outdated index version. Please recompute it with 'createindex'!\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);

        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tdbr->readMmapedDataInMemory();
            tdbr->mlock();
        }

        templateDBIsIndex = false;
    }

    // restrict amount of allocated memory if all results are requested
    // INT_MAX would allocate 72GB RAM per thread for no reason
    maxResListLen = std::min(tdbr->getSize(), maxResListLen);

    // investigate if it makes sense to mask the profile consensus sequence
    if (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        maskMode = 0;
    }

    takeOnlyBestKmer = (par.exactKmerMatching==1) ||
                       (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(querySeqType,Parameters::DBTYPE_AMINO_ACIDS)) ||
                       (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) && Parameters::isEqualDbtype(querySeqType,Parameters::DBTYPE_NUCLEOTIDES));

    // memoryLimit in bytes
    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = par.splitMemoryLimit;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }

    if (templateDBIsIndex == false && sameQTDB == true) {
        qdbr = tdbr;
    } else {
        qdbr = new DBReader<unsigned int>(queryDB.c_str(), queryDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }
    Debug(Debug::INFO) << "Query database size: " << qdbr->getSize() << " type: "<< Parameters::getDbTypeName(querySeqType) << "\n";

    setupSplit(*tdbr, alphabetSize - 1, querySeqType,
               threads, templateDBIsIndex, memoryLimit, qdbr->getSize(),
               maxResListLen, kmerSize, splits, splitMode);

    if(Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) == false){
        const bool isProfileSearch = Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) ||
                                     Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE);
        kmerThr = getKmerThreshold(sensitivity, isProfileSearch, kmerScore, kmerSize);
    }else {
        kmerThr = 0;
    }

    Debug(Debug::INFO) << "Target database size: " << tdbr->getSize() << " type: " <<Parameters::getDbTypeName(targetSeqType) << "\n";

    if (splitMode == Parameters::QUERY_DB_SPLIT) {
        // create the whole index table
        getIndexTable(0, 0, tdbr->getSize());
    } else if (splitMode == Parameters::TARGET_DB_SPLIT) {
        sequenceLookup = NULL;
        indexTable = NULL;
    } else {
        Debug(Debug::ERROR) << "Invalid split mode: " << splitMode << "\n";
        EXIT(EXIT_FAILURE);
    }
}

Prefiltering::~Prefiltering() {
    if (sameQTDB == false) {
        qdbr->close();
        delete qdbr;
    }

    if (indexTable != NULL) {
        delete indexTable;
    }

    if (sequenceLookup != NULL) {
        delete sequenceLookup;
    }

    tdbr->close();
    delete tdbr;

    if (templateDBIsIndex == true) {
        tidxdbr->close();
        delete tidxdbr;
    }
    if(kmerSubMat != ungappedSubMat){
        delete ungappedSubMat;
    }
    delete kmerSubMat;


    if (_2merSubMatrix != NULL && templateDBIsIndex == false) {
        ScoreMatrix::cleanup(_2merSubMatrix);
    }
    if (_3merSubMatrix != NULL && templateDBIsIndex == false) {
        ScoreMatrix::cleanup(_3merSubMatrix);
    }
}

void Prefiltering::reopenTargetDb() {
    if (templateDBIsIndex == true) {
        tidxdbr->close();
        delete tidxdbr;
        tidxdbr = NULL;
    }

    tdbr->close();
    delete tdbr;

    Debug(Debug::INFO) << "Index table not compatible with chosen settings. Compute index.\n";
    tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    tdbr->open(DBReader<unsigned int>::NOSORT);

    if (preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tdbr->readMmapedDataInMemory();
        tdbr->mlock();
    }

    templateDBIsIndex = false;
}

void Prefiltering::setupSplit(DBReader<unsigned int>& tdbr, const int alphabetSize, const unsigned int querySeqTyp, const int threads,
                              const bool templateDBIsIndex, const size_t memoryLimit, const size_t qDbSize,
                               size_t &maxResListLen, int &kmerSize, int &split, int &splitMode) {
    size_t memoryNeeded = estimateMemoryConsumption(1,
                                                  tdbr.getSize(), tdbr.getAminoAcidDBSize(), maxResListLen, alphabetSize,
                                                  kmerSize == 0 ? // if auto detect kmerSize
                                                  IndexTable::computeKmerSize(tdbr.getAminoAcidDBSize()) : kmerSize, querySeqTyp,
                                                  threads);
    
    int optimalSplitMode = Parameters::TARGET_DB_SPLIT;
    if (memoryNeeded > 0.9 * memoryLimit) {
        if (splitMode == Parameters::QUERY_DB_SPLIT) {
            Debug(Debug::ERROR) << "--split-mode was set to query-split (" << Parameters::QUERY_DB_SPLIT << ") but memory limit requires target-split." << 
                                   " Please use a computer with more main memory or run with default --split-mode setting.\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
#ifdef HAVE_MPI
        // TODO Currently we dont support split target indeces, we need MPI TARGET_DB_SPLIT when we have split index support again
        // if (templateDBIsIndex ) {
        //     optimalSplitMode = Parameters::TARGET_DB_SPLIT;
        // } else ...
        // if enough memory + MPI -> optimal mode is query
        optimalSplitMode = Parameters::QUERY_DB_SPLIT;
#endif
    }

    // user split mode is legal and respected, only set this if we are in automatic detection
    if (splitMode == Parameters::DETECT_BEST_DB_SPLIT) {
        splitMode = optimalSplitMode;
    }

    // ideally we always run without splitting
    int minimalNumSplits = 1;
    // get minimal number of splits in case of target split
    // we EXITed already in query split mode
    if (memoryNeeded > 0.9 * memoryLimit) {
        // memory is not enough to compute everything at once
        //TODO add PROFILE_STATE (just 6-mers)
        std::pair<int, int> splitSettings = Prefiltering::optimizeSplit(memoryLimit, &tdbr, alphabetSize, kmerSize, querySeqTyp, threads);
        if (splitSettings.second == -1) {
            Debug(Debug::ERROR) << "Cannot fit databased into " << ByteParser::format(memoryLimit) << ". Please use a computer with more main memory.\n";
            EXIT(EXIT_FAILURE);
        }
        if (kmerSize == 0) {
            // set k-mer based on aa size in database
            // if we have less than 10Mio * 335 amino acids use 6mers
            kmerSize = splitSettings.first;
        }
        minimalNumSplits = splitSettings.second;
    } 

    int optimalNumSplits = minimalNumSplits;

    size_t sizeOfDbToSplit = tdbr.getSize();
    if (splitMode == Parameters::QUERY_DB_SPLIT) {
        sizeOfDbToSplit = qDbSize;
    }
#ifdef HAVE_MPI
    optimalNumSplits = std::max(MMseqsMPI::numProc, optimalNumSplits);
#endif
    optimalNumSplits = std::min((int)sizeOfDbToSplit, optimalNumSplits);

    // set the final number of splits
    if (split == Parameters::AUTO_SPLIT_DETECTION) {
        split = optimalNumSplits;
    } 
    
    // templateDBIsIndex = false when called from indexdb
    if ((split < minimalNumSplits) && (templateDBIsIndex)) {
        Debug(Debug::WARNING) << "split was set to " << split << " but at least " << minimalNumSplits << " are required. Please run with default paramerters\n";
    } else if (split > (int)sizeOfDbToSplit) {
        Debug(Debug::ERROR) << "split was set to " << split << " but the db to split has only " << sizeOfDbToSplit << " sequences. Please run with default paramerters\n";
        EXIT(EXIT_FAILURE);
    }

    if (kmerSize == 0) {
        size_t aaSize = tdbr.getAminoAcidDBSize() / std::max(split, 1);
        kmerSize = IndexTable::computeKmerSize(aaSize);
    }

    // in TARGET_DB_SPLIT we have to reduce the number of prefilter hits can produce,
    // so that the merged database does not contain more than maxResListLen
    if (splitMode == Parameters::TARGET_DB_SPLIT && split > 1) {
        size_t fourTimesStdDeviation = 4 * sqrt(static_cast<double>(maxResListLen) / static_cast<double>(split));
        maxResListLen = std::max(static_cast<size_t>(1), (maxResListLen / split) + fourTimesStdDeviation);
    }

    if (split > 1) {
        Debug(Debug::INFO) << Parameters::getSplitModeName(splitMode) << " split mode. Searching through " << split << " splits\n";
    }

    size_t memoryNeededPerSplit = estimateMemoryConsumption((splitMode == Parameters::TARGET_DB_SPLIT) ? split : 1, tdbr.getSize(),
                                             tdbr.getAminoAcidDBSize(), maxResListLen, alphabetSize, kmerSize, querySeqTyp, threads);
    Debug(Debug::INFO) << "Estimated memory consumption: " << ByteParser::format(memoryNeededPerSplit) << "\n";
    if (memoryNeededPerSplit > 0.9 * memoryLimit) {
        Debug(Debug::WARNING) << "Process needs more than " << ByteParser::format(memoryLimit) << " main memory.\n" <<
                                 "Increase the size of --split or set it to 0 to automatically optimize target database split.\n";
        if (templateDBIsIndex == true) {
            Debug(Debug::WARNING) << "Computed index is too large. Avoid using the index.\n";
        }
    }
}

void Prefiltering::mergeOutput(const std::string &outDB, const std::string &outDBIndex,
                               const std::vector<std::pair<std::string, std::string>> &filenames) {
    Timer timer;
    if (filenames.size() < 2) {
        DBReader<unsigned int>::moveDb(filenames[0].first, outDB);
        Debug(Debug::INFO) << "No merging needed.\n";
        return;
    }

    const std::pair<std::string, std::string> tmpDb = Util::databaseNames((outDB + "_merged"));
    DBWriter writer(tmpDb.first.c_str(), tmpDb.second.c_str(), 1, compressed, Parameters::DBTYPE_PREFILTER_RES);
    writer.open(1024 * 1024 * 1024); // 1 GB buffer
    writer.mergeFilePair(filenames);
    writer.close();
    for (size_t i = 0; i < filenames.size(); i++) {
        // remove split
        DBReader<unsigned int>::removeDb(filenames[i].first);
    }
    // sort merged entries by evalue
    DBReader<unsigned int> dbr(tmpDb.first.c_str(), tmpDb.second.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter dbw(outDB.c_str(), outDBIndex.c_str(), threads, compressed, Parameters::DBTYPE_PREFILTER_RES);
    dbw.open(1024 * 1024 * 1024);
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string result;
        result.reserve(BUFFER_SIZE);
        char buffer[100];
// dynamic schedule causes mis-order of the input files during merge. Needs to be fixed there
#pragma omp for schedule(static)
        for (size_t id = 0; id < dbr.getSize(); id++) {
            unsigned int dbKey = dbr.getDbKey(id);
            char *data = dbr.getData(id, thread_idx);
            std::vector<hit_t> hits = QueryMatcher::parsePrefilterHits(data);
            if (hits.size() > 1) {
                std::sort(hits.begin(), hits.end(), hit_t::compareHitsByScoreAndId);
            }
            for(size_t hit_id = 0; hit_id < hits.size(); hit_id++){
                int len = QueryMatcher::prefilterHitToBuffer(buffer, hits[hit_id]);
                result.append(buffer, len);
            }
            dbw.writeData(result.c_str(), result.size(), dbKey, thread_idx);
            result.clear();
        }
    }

    // TODO: we close with "true" because multiple calls to mergeOutput call mergeFilePair that expects two file (merged input)
    dbw.close(true);
    dbr.close();

    DBReader<unsigned int>::removeDb(tmpDb.first);

    Debug(Debug::INFO) << "\nTime for merging results: " << timer.lap() << "\n";
}


ScoreMatrix *Prefiltering::getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize) {
    // profile only uses the 2mer, 3mer matrix
    if (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        return NULL;
    }

    ScoreMatrix *result = NULL;
    if (templateDBIsIndex == true) {
        switch(kmerSize) {
            case 2:
                result = PrefilteringIndexReader::get2MerScoreMatrix(tidxdbr, false);
                break;
            case 3:
                result = PrefilteringIndexReader::get3MerScoreMatrix(tidxdbr, false);
                break;
            default:
                break;
        }
    }

    if (result != NULL) {
        return result;
    }
    return ExtendedSubstitutionMatrix::calcScoreMatrix(matrix, kmerSize);

}

// TODO reimplement split index feature
void Prefiltering::getIndexTable(int /*split*/, size_t dbFrom, size_t dbSize) {
    if (templateDBIsIndex == true) {
        indexTable = PrefilteringIndexReader::generateIndexTable(tidxdbr, false);
        //TODO check masked mode here
        sequenceLookup = PrefilteringIndexReader::getSequenceLookup(tidxdbr, false);
    } else {
        Timer timer;

        Sequence tseq(maxSeqLen, targetSeqType, kmerSubMat, kmerSize, spacedKmer, aaBiasCorrection, true, spacedKmerPattern);
        int localKmerThr = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) ||
                Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE) ||
                Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES) ||
                            (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) == false && takeOnlyBestKmer == true) ) ? 0 : kmerThr;

        // remove X or N for seeding
        int adjustAlphabetSize = (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) ||
                                  Parameters::isEqualDbtype(targetSeqType,Parameters::DBTYPE_AMINO_ACIDS))
                           ? alphabetSize -1 : alphabetSize;
        indexTable = new IndexTable(adjustAlphabetSize, kmerSize, false);
        SequenceLookup **maskedLookup   = maskMode == 1 || maskLowerCaseMode == 1 ? &sequenceLookup : NULL;
        SequenceLookup **unmaskedLookup = maskMode == 0 ? &sequenceLookup : NULL;

        Debug(Debug::INFO) << "Index table k-mer threshold: " << localKmerThr << " at k-mer size " << kmerSize << " \n";
        IndexBuilder::fillDatabase(indexTable, maskedLookup, unmaskedLookup, *kmerSubMat,  &tseq, tdbr, dbFrom, dbFrom + dbSize, localKmerThr, maskMode, maskLowerCaseMode);

        // sequenceLookup has to be temporarily present to speed up masking
        // afterwards its not needed anymore in DIAG_SCORE_OFF
        if (diagonalScoring == Parameters::DIAG_SCORE_OFF) {
            delete sequenceLookup;
            sequenceLookup = NULL;
        }

        indexTable->printStatistics(kmerSubMat->int2aa);
        tdbr->remapData();
        Debug(Debug::INFO) << "Time for index table init: " << timer.lap() << "\n";
    }

    // init the substitution matrices
    switch (querySeqType  & 0x7FFFFFFF) {
        case Parameters::DBTYPE_AMINO_ACIDS:
            // Do not add X
            kmerSubMat->alphabetSize = kmerSubMat->alphabetSize - 1;
            _2merSubMatrix = getScoreMatrix(*kmerSubMat, 2);
            _3merSubMatrix = getScoreMatrix(*kmerSubMat, 3);
            kmerSubMat->alphabetSize = alphabetSize;
            break;
        case Parameters::DBTYPE_HMM_PROFILE:
        case Parameters::DBTYPE_PROFILE_STATE_PROFILE:
        case Parameters::DBTYPE_NUCLEOTIDES:
        default:
            if (_2merSubMatrix != NULL && templateDBIsIndex == false) {
                delete _2merSubMatrix;
            }

            if (_3merSubMatrix != NULL && templateDBIsIndex == false) {
                delete _3merSubMatrix;
            }
            _2merSubMatrix = NULL;
            _3merSubMatrix = NULL;
            break;
    }
}

bool Prefiltering::isSameQTDB() {
    //  check if when qdb and tdb have the same name an index extension exists
    std::string check(targetDB);
    size_t pos = check.find(queryDB);
    int match = false;
    if (pos == 0) {
        check.replace(0, queryDB.length(), "");
        // TODO name changed to .idx
        PatternCompiler regex("^\\.s?k[5-7]$");
        match = regex.isMatch(check.c_str());
    }
    // if no match found or two matches found (we want exactly one match)
    return (queryDB.compare(targetDB) == 0 || (match == true));
}

void Prefiltering::runAllSplits(const std::string &resultDB, const std::string &resultDBIndex) {
    runSplits(resultDB, resultDBIndex, 0, splits, false);
}

#ifdef HAVE_MPI
void Prefiltering::runMpiSplits(const std::string &resultDB, const std::string &resultDBIndex, const std::string &localTmpPath) {
    if(compressed == true && splitMode == Parameters::TARGET_DB_SPLIT){
            Debug(Debug::ERROR) << "The output of the prefilter cannot be compressed during target split mode. Please remove --compress.\n";
            EXIT(EXIT_FAILURE);
    }

    // if split size is great than nodes than we have to
    // distribute all splits equally over all nodes
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(int i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }

    size_t fromSplit = 0;
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }

    size_t splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    // setting names in case of localTmp path
    std::string procTmpResultDB = localTmpPath;
    std::string procTmpResultDBIndex = localTmpPath;
    if (localTmpPath == "") {
        procTmpResultDB = resultDB;
        procTmpResultDBIndex = resultDBIndex;
    } else {
        procTmpResultDB = procTmpResultDB + "/" + FileUtil::baseName(resultDB);
        procTmpResultDBIndex = procTmpResultDBIndex + "/" + FileUtil::baseName(resultDBIndex);

        if (FileUtil::directoryExists(localTmpPath.c_str()) == false) {
            Debug(Debug::INFO) << "Local tmp dir " << localTmpPath << " does not exist or is not a directory\n";
            if (FileUtil::makeDir(localTmpPath.c_str()) == false) {
                Debug(Debug::ERROR) << "Cannot create local tmp dir " << localTmpPath << "\n";
                EXIT(EXIT_FAILURE);
            } else {
                Debug(Debug::INFO) << "Created local tmp dir " << localTmpPath << "\n";
            }
        }
    }

    std::pair<std::string, std::string> result = Util::createTmpFileNames(procTmpResultDB, procTmpResultDBIndex, MMseqsMPI::rank);
    bool merge = (splitMode == Parameters::QUERY_DB_SPLIT);

    int hasResult = runSplits(result.first, result.second, fromSplit, splitCount, merge) == true ? 1 : 0;

    if (localTmpPath != "") {
        std::pair<std::string, std::string> resultShared = Util::createTmpFileNames(resultDB, resultDBIndex, MMseqsMPI::rank);
        // moveDb takes care if file doesn't exist
        DBReader<unsigned int>::moveDb(result.first, resultShared.first);
    }
    

    int *results = NULL;
    if (MMseqsMPI::isMaster()) {
        results = new int[MMseqsMPI::numProc]();
    }

    MPI_Gather(&hasResult, 1, MPI_INT, results, 1, MPI_INT, MMseqsMPI::MASTER, MPI_COMM_WORLD);
    if (MMseqsMPI::isMaster()) {
        // gather does not write the result of the master into the array
        results[MMseqsMPI::MASTER] = hasResult;

        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int i = 0; i < MMseqsMPI::numProc; ++i) {
            if (results[i] == 1) {
                std::pair<std::string, std::string> resultOfRanki = Util::createTmpFileNames(resultDB, resultDBIndex, i);
                splitFiles.push_back(std::make_pair(resultOfRanki.first, resultOfRanki.second));
            }
        }

        if (splitFiles.size() > 0) {
            // merge output databases
            mergeFiles(resultDB, resultDBIndex, splitFiles);
        } else {
            Debug(Debug::ERROR) << "Aborting. No results were computed!\n";
            EXIT(EXIT_FAILURE);
        }

        delete [] results;
    }

}
#endif

int Prefiltering::runSplits(const std::string &resultDB, const std::string &resultDBIndex,
                            size_t fromSplit, size_t splitProcessCount, bool merge) {
    if (fromSplit + splitProcessCount > static_cast<size_t>(splits)) {
        Debug(Debug::ERROR) << "Start split " << fromSplit << " plus split count " << splitProcessCount << " cannot be larger than splits " << splits << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t freeSpace = FileUtil::getFreeSpace(FileUtil::dirName(resultDB).c_str());
    size_t estimatedHDDMemory = estimateHDDMemoryConsumption(qdbr->getSize(), maxResListLen);
    if (freeSpace < estimatedHDDMemory) {
        std::string freeSpaceToPrint = ByteParser::format(freeSpace);
        std::string estimatedHDDMemoryToPrint = ByteParser::format(estimatedHDDMemory);
        Debug(Debug::WARNING) << "Hard disk might not have enough free space (" << freeSpaceToPrint << " left)."
                              << "The prefilter result might need up to " << estimatedHDDMemoryToPrint << ".\n";
    }

    bool hasResult = false;
    if (splitProcessCount > 1) {
        if(compressed == true && splitMode == Parameters::TARGET_DB_SPLIT){
            Debug(Debug::ERROR) << "The output of the prefilter cannot be compressed during target split mode. Please remove --compress.\n";
            EXIT(EXIT_FAILURE);
        }
        // splits template database into x sequence steps
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for (size_t i = fromSplit; i < (fromSplit + splitProcessCount); i++) {
            std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(resultDB, resultDBIndex, i);
            if (runSplit(filenamePair.first.c_str(), filenamePair.second.c_str(), i, merge)) {
                splitFiles.push_back(filenamePair);

            }
        }
        if (splitFiles.size() > 0) {
            mergeFiles(resultDB, resultDBIndex, splitFiles);
            hasResult = true;
        }
    } else if (splitProcessCount == 1) {
        if (runSplit(resultDB.c_str(), resultDBIndex.c_str(), fromSplit, merge)) {
            hasResult = true;
        }
    }

    return hasResult;
}

bool Prefiltering::runSplit(const std::string &resultDB, const std::string &resultDBIndex, size_t split, bool merge) {
    Debug(Debug::INFO) << "Process prefiltering step " << (split + 1) << " of " << splits << "\n\n";

    size_t dbFrom = 0;
    size_t dbSize = tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();
    
    size_t resListOffset = 0;
    std::vector<std::string> prevMaxVals = Util::split(prevMaxResListLengths, ",");
    for (size_t i = 0; i < prevMaxVals.size(); ++i) {
        size_t value = strtoull(prevMaxVals[i].c_str(), NULL, 10);
        if (splitMode == Parameters::TARGET_DB_SPLIT) {
            size_t fourTimesStdDeviation = splits > 1 ? static_cast<size_t>(4 * sqrt(static_cast<double>(value) / static_cast<double>(splits))) : 0;
            resListOffset += std::max(static_cast<size_t>(1), (value / splits) + fourTimesStdDeviation);
        } else {
            resListOffset += value;
        }
    }

    // create index table based on split parameter
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        Util::decomposeDomainByAminoAcid(tdbr->getDataSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                         split, splits, &dbFrom, &dbSize);
        if (dbSize == 0) {
            return false;
        }

        if (indexTable != NULL) {
            delete indexTable;
            indexTable = NULL;
        }

        if (sequenceLookup != NULL) {
            delete sequenceLookup;
            sequenceLookup = NULL;
        }

        getIndexTable(split, dbFrom, dbSize);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        Util::decomposeDomainByAminoAcid(qdbr->getDataSize(), qdbr->getSeqLens(), qdbr->getSize(), split, splits, &queryFrom, &querySize);
        if (querySize == 0) {
            return false;
        }
    }

    Debug(Debug::INFO) << "k-mer similarity threshold: " << kmerThr << "\n";

    double kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    size_t totalQueryDBSize = querySize;

    unsigned int localThreads = 1;
#ifdef OPENMP
    localThreads = std::min((unsigned int)threads, (unsigned int)querySize);
#endif

    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), localThreads, compressed, Parameters::DBTYPE_PREFILTER_RES);
    tmpDbw.open();

    // init all thread-specific data structures
    char *notEmpty = new char[querySize];
    memset(notEmpty, 0, querySize * sizeof(char)); // init notEmpty

    std::list<int> **reslens = new std::list<int> *[localThreads];
    for (unsigned int i = 0; i < localThreads; ++i) {
        reslens[i] = new std::list<int>();
    }

    Debug(Debug::INFO) << "Starting prefiltering scores calculation (step " << (split + 1) << " of " << splits << ")\n";
    Debug(Debug::INFO) << "Query db start " << (queryFrom + 1) << " to " << queryFrom + querySize << "\n";
    Debug(Debug::INFO) << "Target db start " << (dbFrom + 1) << " to " << dbFrom + dbSize << "\n";
    Debug::Progress progress(querySize);

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(maxSeqLen, querySeqType, kmerSubMat, kmerSize, spacedKmer, aaBiasCorrection, true, spacedKmerPattern);

        QueryMatcher matcher(indexTable, sequenceLookup, kmerSubMat,  ungappedSubMat,
                            kmerThr, kmerSize, dbSize, maxSeqLen, maxResListLen, aaBiasCorrection,
                            diagonalScoring, minDiagScoreThr, takeOnlyBestKmer, resListOffset);

        if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
            matcher.setProfileMatrix(seq.profile_matrix);
        } else {
            matcher.setSubstitutionMatrix(_3merSubMatrix, _2merSubMatrix);
        }

#pragma omp for schedule(dynamic, 2) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
        for (size_t id = queryFrom; id < queryFrom + querySize; id++) {
            progress.updateProgress();
            // get query sequence
            char *seqData = qdbr->getData(id, thread_idx);
            unsigned int qKey = qdbr->getDbKey(id);
            seq.mapSequence(id, qKey, seqData);
            size_t targetSeqId = UINT_MAX;
            if (sameQTDB || includeIdentical) {
                targetSeqId = tdbr->getId(seq.getDbKey());
                // only the corresponding split should include the id (hack for the hack)
                if (targetSeqId >= dbFrom && targetSeqId < (dbFrom + dbSize) && targetSeqId != UINT_MAX) {
                    targetSeqId = targetSeqId - dbFrom;
                    if(targetSeqId > tdbr->getSize()){
                        Debug(Debug::ERROR) << "targetSeqId: " << targetSeqId << " > target database size: "  << tdbr->getSize() <<  "\n";
                        EXIT(EXIT_FAILURE);
                    }
                }else{
                    targetSeqId = UINT_MAX;
                }
            }
            // calculate prefiltering results
            std::pair<hit_t *, size_t> prefResults = matcher.matchQuery(&seq, targetSeqId);
            size_t resultSize = prefResults.second;
            // write
            writePrefilterOutput(&tmpDbw, thread_idx, id, prefResults, dbFrom);

            // update statistics counters
            if (resultSize != 0) {
                notEmpty[id - queryFrom] = 1;
            }

            kmersPerPos += matcher.getStatistics()->kmersPerPos;
            dbMatches += matcher.getStatistics()->dbMatches;
            doubleMatches += matcher.getStatistics()->doubleMatches;
            querySeqLenSum += seq.L;
            diagonalOverflow += matcher.getStatistics()->diagonalOverflow;
            resSize += resultSize;
            realResSize += std::min(resultSize, maxResListLen);
            reslens[thread_idx]->emplace_back(resultSize);
        } // step end
    }

    if (Debug::debugLevel >= Debug::INFO) {
        statistics_t stats(kmersPerPos / static_cast<double>(totalQueryDBSize),
                           dbMatches / totalQueryDBSize,
                           doubleMatches / totalQueryDBSize,
                           querySeqLenSum, diagonalOverflow,
                           resSize / totalQueryDBSize);

        size_t empty = 0;
        for (size_t id = 0; id < querySize; id++) {
            if (notEmpty[id] == 0) {
                empty++;
            }
        }

        printStatistics(stats, reslens, localThreads, empty, maxResListLen);
    }

    if (splitMode == Parameters::TARGET_DB_SPLIT && splits == 1) {
#ifdef HAVE_MPI
        // if a mpi rank processed a single split, it must have it merged before all ranks can be united
        tmpDbw.close(true);
#else
        tmpDbw.close(merge);
#endif
    } else {
        tmpDbw.close(merge);
    }

    // sort by ids
    // needed to speed up merge later on
    // sorts this datafile according to the index file
    if (splitMode == Parameters::TARGET_DB_SPLIT && splits > 1) {
        // delete indexTable to free memory:
        if (indexTable != NULL) {
            delete indexTable;
            indexTable = NULL;
        }
        DBReader<unsigned int> resultReader(tmpDbw.getDataFileName(), tmpDbw.getIndexFileName(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        resultReader.open(DBReader<unsigned int>::NOSORT);
        resultReader.readMmapedDataInMemory();
        
        const std::pair<std::string, std::string> tempDb = Util::databaseNames((resultDB + "_tmp"));
        DBWriter resultWriter(tempDb.first.c_str(), tempDb.second.c_str(), localThreads, compressed, Parameters::DBTYPE_PREFILTER_RES);
        resultWriter.open();
        resultWriter.sortDatafileByIdOrder(resultReader);
        resultWriter.close(true);
        resultReader.close();
        DBReader<unsigned int>::removeDb(resultDB);
        DBReader<unsigned int>::moveDb(tempDb.first, resultDB);
    }

    for (unsigned int i = 0; i < localThreads; i++) {
        reslens[i]->clear();
        delete reslens[i];
    }
    delete[] reslens;
    delete[] notEmpty;

    return true;
}

void Prefiltering::writePrefilterOutput(DBWriter *dbWriter, unsigned int thread_idx, size_t id,
                                        const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset) {
    hit_t *resultVector = prefResults.first;
    std::string prefResultsOutString;
    prefResultsOutString.reserve(BUFFER_SIZE);
    char buffer[100];
    for (size_t i = 0; i < prefResults.second; i++) {
        hit_t *res = resultVector + i;
        // correct the 0 indexed sequence id again to its real identifier
        size_t targetSeqId = res->seqId + seqIdOffset;
        // replace id with key
        res->seqId = tdbr->getDbKey(targetSeqId);
        if (targetSeqId >= tdbr->getSize()) {
            Debug(Debug::WARNING) << "Wrong prefiltering result for query: " << qdbr->getDbKey(id) << " -> " << targetSeqId
                                  << "\t" << res->prefScore << "\n";
        }

        // TODO: check if this should happen when diagonalScoring == Parameters::DIAG_SCORE_OFF
        if (covThr > 0.0 && (covMode == Parameters::COV_MODE_BIDIRECTIONAL || covMode == Parameters::COV_MODE_QUERY)) {
            float queryLength = static_cast<float>(qdbr->getSeqLens(id));
            float targetLength = static_cast<float>(tdbr->getSeqLens(targetSeqId));
            if (Util::canBeCovered(covThr, covMode, queryLength, targetLength) == false) {
                continue;
            }
        }

        // write prefiltering results to a string
        int len = QueryMatcher::prefilterHitToBuffer(buffer, *res);
        // TODO: error handling for len
        prefResultsOutString.append(buffer, len);
    }
    dbWriter->writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), qdbr->getDbKey(id), thread_idx);
}

void Prefiltering::printStatistics(const statistics_t &stats, std::list<int> **reslens,
                                   unsigned int resLensSize, size_t empty, size_t maxResults) {
    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (unsigned int i = 1; i < resLensSize; i++) {
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }
    Debug(Debug::INFO) << "\n" << stats.kmersPerPos << " k-mers per position\n";
    Debug(Debug::INFO) << stats.dbMatches << " DB matches per sequence\n";
    Debug(Debug::INFO) << stats.diagonalOverflow << " overflows\n";
    Debug(Debug::INFO) << stats.resultsPassedPrefPerSeq << " sequences passed prefiltering per query sequence";
    if (stats.resultsPassedPrefPerSeq > maxResults)
        Debug(Debug::WARNING) << " (ATTENTION: max. " << maxResults
                           << " best scoring sequences were written to the output prefiltering database)\n";
    else
        Debug(Debug::INFO) << "\n";
    size_t mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
    std::advance(it, mid);
    Debug(Debug::INFO) << *it << " median result list length\n";
    Debug(Debug::INFO) << empty << " sequences with 0 size result lists\n";
}


BaseMatrix *Prefiltering::getSubstitutionMatrix(const ScoreMatrixFile &scoringMatrixFile, size_t alphabetSize, float bitFactor, bool profileState, bool isNucl) {
    BaseMatrix *subMat;

    if (isNucl){
        subMat = new NucleotideMatrix(scoringMatrixFile.nucleotides, bitFactor, 0.0);
    } else if (alphabetSize < 21) {
        SubstitutionMatrix sMat(scoringMatrixFile.aminoacids, bitFactor, -0.2f);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2int, sMat.int2aa, sMat.alphabetSize, alphabetSize, bitFactor);
    }else if(profileState == true){
        SubstitutionMatrix sMat(scoringMatrixFile.aminoacids, bitFactor, -0.2f);
        subMat = new SubstitutionMatrixProfileStates(sMat.matrixName, sMat.probMatrix, sMat.pBack,
                                                     sMat.subMatrixPseudoCounts, bitFactor, 0.0, 8);
    } else {
        subMat = new SubstitutionMatrix(scoringMatrixFile.aminoacids, bitFactor, -0.2f);
    }
    return subMat;
}

void Prefiltering::mergeFiles(const std::string &outDB, const std::string &outDBIndex,
                              const std::vector<std::pair<std::string, std::string>> &splitFiles) {
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        mergeOutput(outDB, outDBIndex, splitFiles);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        DBWriter::mergeResults(outDB, outDBIndex, splitFiles);
    }
}

int Prefiltering::getKmerThreshold(const float sensitivity, const bool isProfile, const int kmerScore, const int kmerSize) {
    double kmerThrBest = kmerScore;
    if (kmerScore == INT_MAX) {
        if(isProfile){
            if (kmerSize == 5) {
                float base = 140.75;
                kmerThrBest = base - (sensitivity * 8.75);
            } else if (kmerSize == 6) {
                float base = 155.75;
                kmerThrBest = base - (sensitivity * 8.75);
            } else if (kmerSize == 7) {
                float base = 171.75;
                kmerThrBest = base - (sensitivity * 9.75);
            } else {
                Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
                EXIT(EXIT_FAILURE);
            }
        }else{
            if (kmerSize == 5) {
                float base = 160.75;
                kmerThrBest = base - (sensitivity * 12.75);
            } else if (kmerSize == 6) {
                float base = 163.2;
                kmerThrBest = base - (sensitivity * 8.917);
            } else if (kmerSize == 7) {
                float base = 186.15;
                kmerThrBest = base - (sensitivity * 11.22);
            } else {
                Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
    return static_cast<int>(kmerThrBest);
}

size_t Prefiltering::estimateMemoryConsumption(int split, size_t dbSize, size_t resSize,
                                               size_t maxResListLen,
                                               int alphabetSize, int kmerSize, unsigned int querySeqType,
                                               int threads) {
    // for each residue in the database we need 7 byte
    size_t dbSizeSplit = (dbSize) / split;
    size_t residueSize = (resSize / split * 7);
    // 21^7 * pointer size is needed for the index
    size_t indexTableSize = static_cast<size_t>(pow(alphabetSize, kmerSize)) * sizeof(size_t *);
    // memory needed for the threads
    // This memory is an approx. for Countint32Array and QueryTemplateLocalFast
    size_t threadSize = threads * (
            (dbSizeSplit * 2 * sizeof(IndexEntryLocal)) // databaseHits in QueryMatcher
            + (dbSizeSplit * sizeof(CounterResult)) // databaseHits in QueryMatcher
            + (maxResListLen * sizeof(hit_t))
            + (dbSizeSplit * 2 * sizeof(CounterResult) * 2) // BINS * binSize, (binSize = dbSize * 2 / BINS)
            // 2 is a security factor the size can increase during run
    );

    // extended matrix
    size_t extendedMatrix = 0;
    if(Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)){
        extendedMatrix = sizeof(std::pair<short, unsigned int>) * static_cast<size_t>(pow(pow(alphabetSize, 3), 2));
        extendedMatrix += sizeof(std::pair<short, unsigned int>) * pow(pow(alphabetSize, 2), 2);
    }
    // some memory needed to keep the index, ....
    size_t background = dbSize * 22;
    // return result in bytes
    return residueSize + indexTableSize + threadSize + background + extendedMatrix;
}

size_t Prefiltering::estimateHDDMemoryConsumption(size_t dbSize, size_t maxResListLen) {
    // 21 bytes is roughly the size of an entry
    // 2x because the merge doubles the hdd demand
    return 2 * (21 * dbSize * maxResListLen);
}

std::pair<int, int> Prefiltering::optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr,
                                                int alphabetSize, int externalKmerSize, unsigned int querySeqType, unsigned int threads) {
    for (int optSplit = 1; optSplit < 100; optSplit++) {
        for (int optKmerSize = 6; optKmerSize <= 7; optKmerSize++) {
            if (optKmerSize == externalKmerSize || externalKmerSize == 0) { // 0: set k-mer based on aa size in database
                size_t aaUpperBoundForKmerSize = IndexTable::getUpperBoundAACountForKmerSize(optKmerSize);
                if ((tdbr->getAminoAcidDBSize() / optSplit) < aaUpperBoundForKmerSize) {
                    size_t neededSize = estimateMemoryConsumption(optSplit, tdbr->getSize(), tdbr->getAminoAcidDBSize(),
                                                                  0, alphabetSize, optKmerSize, querySeqType, threads);
                    if (neededSize < 0.9 * totalMemoryInByte) {
                        return std::make_pair(optKmerSize, optSplit);
                    }
                }
            }
        }
    }

    return std::make_pair(-1, -1);
}


