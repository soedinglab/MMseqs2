#include "Prefiltering.h"
#include "NucleotideMatrix.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrixProfileStates.h"
#include "DBWriter.h"
#include "QueryMatcherTaxonomyHook.h"

#include "PatternCompiler.h"
#include "FileUtil.h"
#include "IndexBuilder.h"
#include "Timer.h"
#include "ByteParser.h"
#include "Parameters.h"
#include "MemoryMapped.h"
#include "FastSort.h"
#include <sys/mman.h>

#ifdef OPENMP
#include <omp.h>
#endif

Prefiltering::Prefiltering(const std::string &queryDB,
                           const std::string &queryDBIndex,
                           const std::string &targetDB,
                           const std::string &targetDBIndex,
                           int querySeqType, int targetSeqType,
                           const Parameters &par) :
        queryDB(queryDB),
        queryDBIndex(queryDBIndex),
        targetDB(targetDB),
        targetDBIndex(targetDBIndex),
        splits(par.split),
        kmerSize(par.kmerSize),
        spacedKmerPattern(par.spacedKmerPattern),
        localTmp(par.localTmp),
        spacedKmer(par.spacedKmer != 0),
        maskMode(par.maskMode),
        maskLowerCaseMode(par.maskLowerCaseMode),
        maskProb(par.maskProb),
        splitMode(par.splitMode),
        scoringMatrixFile(par.scoringMatrixFile),
        seedScoringMatrixFile(par.seedScoringMatrixFile),
        targetSeqType(targetSeqType),
        targetSearchMode(par.targetSearchMode),
        maxResListLen(par.maxResListLen),
        sensitivity(par.sensitivity),
        maxSeqLen(par.maxSeqLen),
        querySeqType(querySeqType),
        diagonalScoring(par.diagonalScoring),
        minDiagScoreThr(static_cast<unsigned int>(par.minDiagScoreThr)),
        aaBiasCorrection(par.compBiasCorrection != 0),
        aaBiasCorrectionScale(par.compBiasCorrectionScale),
        covThr(par.covThr), covMode(par.covMode), includeIdentical(par.includeIdentity),
        preloadMode(par.preloadMode),
        threads(static_cast<unsigned int>(par.threads)),
        compressed(par.compressed) {
    sameQTDB = isSameQTDB();

    // init the substitution matrices
    switch (querySeqType & Parameters::DBTYPE_MASK) {
        case Parameters::DBTYPE_NUCLEOTIDES:
            kmerSubMat = getSubstitutionMatrix(scoringMatrixFile, par.alphabetSize, 1.0, false, true);
            ungappedSubMat = kmerSubMat;
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        case Parameters::DBTYPE_AMINO_ACIDS:
            kmerSubMat = getSubstitutionMatrix(seedScoringMatrixFile, par.alphabetSize, 8.0, false, false);
            ungappedSubMat = getSubstitutionMatrix(scoringMatrixFile, par.alphabetSize, 2.0, false, false);
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        case Parameters::DBTYPE_HMM_PROFILE:
            // needed for Background distributions
            kmerSubMat = getSubstitutionMatrix(scoringMatrixFile, par.alphabetSize, 8.0, false, false);
            ungappedSubMat = getSubstitutionMatrix(scoringMatrixFile, par.alphabetSize, 2.0, false, false);
            alphabetSize = kmerSubMat->alphabetSize;
            break;
        default:
            Debug(Debug::ERROR) << "Query sequence type not implemented!\n";
            EXIT(EXIT_FAILURE);
    }

    if (Parameters::isEqualDbtype(FileUtil::parseDbType(targetDB.c_str()), Parameters::DBTYPE_INDEX_DB)) {
        if (preloadMode == Parameters::PRELOAD_MODE_AUTO) {
            if (sensitivity > 6.0) {
                preloadMode = Parameters::PRELOAD_MODE_FREAD;
            } else {
                preloadMode = Parameters::PRELOAD_MODE_MMAP_TOUCH;
            }
        }

        tidxdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tidxdbr->open(DBReader<unsigned int>::NOSORT);

        templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tidxdbr);
        if (templateDBIsIndex == true) {
            tdbr = PrefilteringIndexReader::openNewReader(tidxdbr, PrefilteringIndexReader::DBR1DATA, PrefilteringIndexReader::DBR1INDEX, false, threads, false, false);
            PrefilteringIndexReader::printSummary(tidxdbr);
            PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
            for (size_t i = 0; i < par.prefilter.size(); i++) {
                if (par.prefilter[i]->wasSet == false) {
                    continue;
                }
                if(par.prefilter[i]->uniqid == par.PARAM_K.uniqid) {
                    if (kmerSize != 0 && data.kmerSize != kmerSize) {
                        Debug(Debug::WARNING) << "Index was created with -k " << data.kmerSize << " but the prefilter was called with -k " << kmerSize << "!\n";
                        Debug(Debug::WARNING) << "Search with -k " << data.kmerSize << "\n";
                    }
                }
                if(par.prefilter[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid) {
                    if (data.alphabetSize != alphabetSize) {
                        Debug(Debug::WARNING) << "Index was created with --alph-size  " << data.alphabetSize << " but the prefilter was called with --alph-size " << alphabetSize << "!\n";
                        Debug(Debug::WARNING) << "Current search will use --alph-size " << data.alphabetSize << "\n";
                    }
                }
                if(par.prefilter[i]->uniqid == par.PARAM_SPACED_KMER_MODE.uniqid) {
                    if (data.spacedKmer != spacedKmer) {
                        Debug(Debug::WARNING) << "Index was created with --spaced-kmer-mode " << data.spacedKmer << " but the prefilter was called with --spaced-kmer-mode " << spacedKmer << "!\n";
                        Debug(Debug::WARNING) << "Current search will use  --spaced-kmer-mode " << data.spacedKmer << "\n";
                    }
                }
                if(par.prefilter[i]->uniqid == par.PARAM_NO_COMP_BIAS_CORR.uniqid) {
                    if (data.compBiasCorr != aaBiasCorrection && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
                        Debug(Debug::WARNING) << "Index was created with --comp-bias-corr " << data.compBiasCorr << " please recreate index with --comp-bias-corr " << aaBiasCorrection << "!\n";
                        Debug(Debug::WARNING) << "Current search will use --comp-bias-corr " << data.compBiasCorr << "\n";
                    }
                }
                if(par.prefilter[i]->uniqid == par.PARAM_SPLIT.uniqid) {
                    if (splitMode == Parameters::TARGET_DB_SPLIT && data.splits != splits) {
                        Debug(Debug::WARNING) << "Index was created with --splits " << data.splits << " please recreate index with --splits " << splits << "!\n";
                        Debug(Debug::WARNING) << "Current search will use --splits " << data.splits << "\n";
                    }
                }
            }

            kmerSize = data.kmerSize;
            alphabetSize = data.alphabetSize;
            targetSeqType = data.seqType;
            spacedKmer = data.spacedKmer == 1 ? true : false;
            // the query database could have longer sequences than the target database, do not cut them short
            maxSeqLen = std::max(maxSeqLen, (size_t)data.maxSeqLength);
            aaBiasCorrection = data.compBiasCorr;

            if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) &&
                Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
                Debug(Debug::ERROR) << "Query-profiles cannot be searched against a target-profile database!\n";
                EXIT(EXIT_FAILURE);
            }

            splits = data.splits;
            if (data.splits > 1) {
                splitMode = Parameters::TARGET_DB_SPLIT;
            }
            spacedKmer = data.spacedKmer != 0;
            spacedKmerPattern = PrefilteringIndexReader::getSpacedPattern(tidxdbr);
            seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(PrefilteringIndexReader::getSubstitutionMatrix(tidxdbr));
        } else {
            Debug(Debug::ERROR) << "Outdated index version. Please recompute it with 'createindex'!\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
        templateDBIsIndex = false;
    }

    // restrict amount of allocated memory if all results are requested
    // INT_MAX would allocate 72GB RAM per thread for no reason
    maxResListLen = std::min(tdbr->getSize(), maxResListLen);

    // investigate if it makes sense to mask the profile consensus sequence
    if (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        maskMode = 0;
    }

    takeOnlyBestKmer = (par.exactKmerMatching==1) ||
                       (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(querySeqType,Parameters::DBTYPE_AMINO_ACIDS)) ||
                       (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) && Parameters::isEqualDbtype(querySeqType,Parameters::DBTYPE_NUCLEOTIDES)) ||
                       (targetSearchMode == 1);

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    if (templateDBIsIndex == false && sameQTDB == true) {
        qdbr = tdbr;
    } else {
        qdbr = new DBReader<unsigned int>(queryDB.c_str(), queryDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }
    Debug(Debug::INFO) << "Query database size: " << qdbr->getSize() << " type: " << Parameters::getDbTypeName(querySeqType) << "\n";

    setupSplit(*tdbr, alphabetSize - 1, querySeqType,
               threads, templateDBIsIndex, memoryLimit, qdbr->getSize(),
               maxResListLen, kmerSize, splits, splitMode);

    if(Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) == false){
        const bool isProfileSearch = Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) ||
                                     Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE);
        const bool queryCPC = DBReader<unsigned int>::getExtendedDbtype(querySeqType) & Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS;
        const bool targetCPC = DBReader<unsigned int>::getExtendedDbtype(targetSeqType) & Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS;
        const bool contextPseudoCnts = queryCPC || targetCPC;
        kmerThr = getKmerThreshold(sensitivity, isProfileSearch, contextPseudoCnts, par.kmerScore.values, kmerSize);
    }else {
        kmerThr = 0;
    }

    Debug(Debug::INFO) << "Target database size: " << tdbr->getSize() << " type: " <<Parameters::getDbTypeName(targetSeqType) << "\n";

    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)) {
        kmerSubMat->alphabetSize = kmerSubMat->alphabetSize - 1;
        _2merSubMatrix = getScoreMatrix(*kmerSubMat, 2);
        _3merSubMatrix = getScoreMatrix(*kmerSubMat, 3);
        kmerSubMat->alphabetSize = alphabetSize;
    }

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



    if (par.taxonList.length() > 0) {
        taxonomyHook = new QueryMatcherTaxonomyHook(targetDB, tdbr, par.taxonList);
    } else {
        taxonomyHook = NULL;
    }
}

Prefiltering::~Prefiltering() {
    if (taxonomyHook != NULL) {
        delete taxonomyHook;
    }

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

    if (templateDBIsIndex == false || preloadMode == Parameters::PRELOAD_MODE_FREAD) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(_3merSubMatrix);
        ExtendedSubstitutionMatrix::freeScoreMatrix(_2merSubMatrix);
    }

    if (kmerSubMat != ungappedSubMat) {
        delete ungappedSubMat;
    }
    delete kmerSubMat;
}

void Prefiltering::setupSplit(DBReader<unsigned int>& tdbr, const int alphabetSize, const unsigned int querySeqTyp, const int threads,
                              const bool templateDBIsIndex, const size_t memoryLimit, const size_t qDbSize,
                              size_t &maxResListLen, int &kmerSize, int &split, int &splitMode) {
    size_t memoryNeeded = estimateMemoryConsumption(1, tdbr.getSize(), tdbr.getAminoAcidDBSize(), maxResListLen, alphabetSize,
                                                    kmerSize == 0 ? // if auto detect kmerSize
                                                    IndexTable::computeKmerSize(tdbr.getAminoAcidDBSize()) : kmerSize, querySeqTyp, threads);

    int optimalSplitMode = Parameters::TARGET_DB_SPLIT;
    if (memoryNeeded > 0.9 * memoryLimit) {
        if (splitMode == Parameters::QUERY_DB_SPLIT) {
            Debug(Debug::ERROR) << "--split-mode was set to query-split (" << Parameters::QUERY_DB_SPLIT << ") but memory limit requires target-split." <<
                                " Please use a computer with more main memory or run with default --split-mode setting.\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
#ifdef HAVE_MPI
        if (templateDBIsIndex) {
             optimalSplitMode = Parameters::TARGET_DB_SPLIT;
         } else {
             optimalSplitMode = Parameters::QUERY_DB_SPLIT;
         }
#else
        optimalSplitMode = Parameters::QUERY_DB_SPLIT;
#endif
    }

    // user split mode is legal and respected, only set this if we are in automatic detection
    if (splitMode == Parameters::DETECT_BEST_DB_SPLIT) {
        splitMode = optimalSplitMode;
    }

    // ideally we always run without splitting
    size_t minimalNumSplits = 1;
    // get minimal number of splits in case of target split
    // we EXITed already in query split mode
    if (memoryNeeded > 0.9 * memoryLimit) {
        // memory is not enough to compute everything at once
        //TODO add PROFILE_STATE (just 6-mers)
        std::pair<int, int> splitSettings = Prefiltering::optimizeSplit(memoryLimit, &tdbr, alphabetSize, kmerSize, querySeqTyp, threads);
        if (splitSettings.second == -1) {
            Debug(Debug::ERROR) << "Cannot fit databases into " << ByteParser::format(memoryLimit) << ". Please use a computer with more main memory.\n";
            EXIT(EXIT_FAILURE);
        }
        if (kmerSize == 0) {
            // set k-mer based on aa size in database
            // if we have less than 10Mio * 335 amino acids use 6mers
            kmerSize = splitSettings.first;
        }
        minimalNumSplits = splitSettings.second;
    }

    size_t optimalNumSplits = minimalNumSplits;
    size_t sizeOfDbToSplit = tdbr.getSize();
    if (splitMode == Parameters::QUERY_DB_SPLIT) {
        sizeOfDbToSplit = qDbSize;
    }
#ifdef HAVE_MPI
    optimalNumSplits = std::max(static_cast<size_t>(std::max(MMseqsMPI::numProc, 1)), optimalNumSplits);
#endif
    optimalNumSplits = std::min(sizeOfDbToSplit, optimalNumSplits);

    // set the final number of splits
    if (split == 0) {
        if(optimalNumSplits > INT_MAX){
            Debug(Debug::ERROR) << "optimalNumSplits is greater INT_MAX\n";
            EXIT(EXIT_FAILURE);
        }
        split = optimalNumSplits;
    }

    // templateDBIsIndex = false when called from indexdb
    if ((static_cast<size_t>(split) < minimalNumSplits) && (templateDBIsIndex)) {
        Debug(Debug::WARNING) << "split was set to " << split << " but at least " << minimalNumSplits << " are required. Please run with default paramerters\n";
    } else if (static_cast<size_t>(split) > sizeOfDbToSplit) {
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

void Prefiltering::mergeTargetSplits(const std::string &outDB, const std::string &outDBIndex, const std::vector<std::pair<std::string, std::string>> &fileNames, unsigned int threads) {
    // we assume that the hits are in the same order
    const size_t splits = fileNames.size();

    if (splits < 2) {
        DBReader<unsigned int>::moveDb(fileNames[0].first, outDB);
        Debug(Debug::INFO) << "No merging needed.\n";
        return;
    }

    Timer timer;
    Debug(Debug::INFO) << "Merging " << splits << " target splits to " << FileUtil::baseName(outDB) << "\n";
    DBReader<unsigned int> reader1(fileNames[0].first.c_str(), fileNames[0].second.c_str(), 1, DBReader<unsigned int>::USE_INDEX);
    reader1.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::Index *index1 = reader1.getIndex();

    size_t totalSize = 0;
    for (size_t id = 0; id < reader1.getSize(); id++) {
        totalSize += index1[id].length;
    }
    for (size_t i = 1; i < splits; ++i) {
        DBReader<unsigned int> reader2(fileNames[i].first.c_str(), fileNames[i].second.c_str(), 1, DBReader<unsigned int>::USE_INDEX);
        reader2.open(DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int>::Index *index2 = reader2.getIndex();
        size_t currOffset = 0;
        for (size_t id = 0; id < reader1.getSize(); id++) {
            // add length for file1 and file2 and subtract -1 for one null byte
            size_t seqLen = index1[id].length + index2[id].length - 1;
            totalSize += index2[id].length - 1;
            index1[id].length = seqLen;
            index1[id].offset = currOffset;
            currOffset += seqLen;
        }
        reader2.close();
    }
    reader1.setDataSize(totalSize);

    FILE ** files = new FILE*[fileNames.size()];
    char ** dataFile = new char*[fileNames.size()];
    size_t * dataFileSize = new size_t[fileNames.size()];
    size_t globalIdOffset = 0;
    for (size_t i = 0; i < splits; ++i) {
        files[i] = FileUtil::openFileOrDie(fileNames[i].first.c_str(), "r", true);
        dataFile[i] = static_cast<char*>(FileUtil::mmapFile(files[i], &dataFileSize[i]));
#ifdef HAVE_POSIX_MADVISE
        if (dataFileSize[i] > 0 && posix_madvise (dataFile[i], dataFileSize[i], POSIX_MADV_SEQUENTIAL) != 0){
            Debug(Debug::ERROR) << "posix_madvise returned an error " << fileNames[i].first << "\n";
        }
#endif

    }
    Debug(Debug::INFO) << "Preparing offsets for merging: " << timer.lap() << "\n";
    // merge target splits data files and sort the hits at the same time
    // TODO: compressed?
    DBWriter writer(outDB.c_str(), outDBIndex.c_str(), threads, 0, Parameters::DBTYPE_PREFILTER_RES);
    writer.open();

    Debug::Progress progress(reader1.getSize());
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string result;
        result.reserve(1024);
        std::vector<hit_t> hits;
        hits.reserve(300);
        char buffer[1024];
        size_t * currentDataFileOffset = new size_t[splits];
        memset(currentDataFileOffset, 0, sizeof(size_t)*splits);
        size_t currentId = __sync_fetch_and_add(&(globalIdOffset), 1);
        size_t prevId = 0;
        while(currentId < reader1.getSize()){
            progress.updateProgress();
            for(size_t file = 0; file < splits; file++){
                size_t tmpId = prevId;
                size_t pos;
                for(pos = currentDataFileOffset[file]; pos < dataFileSize[file] && tmpId != currentId; pos++){
                    tmpId += (dataFile[file][pos] == '\0');
                    currentDataFileOffset[file] = pos;
                }
                currentDataFileOffset[file] = pos;
                QueryMatcher::parsePrefilterHits(&dataFile[file][pos], hits);
            }
            if (hits.size() > 1) {
                SORT_SERIAL(hits.begin(), hits.end(), hit_t::compareHitsByScoreAndId);
            }
            for (size_t i = 0; i < hits.size(); ++i) {
                int len = QueryMatcher::prefilterHitToBuffer(buffer, hits[i]);
                result.append(buffer, len);
            }
            writer.writeData(result.c_str(), result.size(), reader1.getDbKey(currentId), thread_idx);
            hits.clear();
            result.clear();
            prevId = currentId;
            currentId = __sync_fetch_and_add(&(globalIdOffset), 1);
        }

            delete[] currentDataFileOffset;
    }
    writer.close();
    reader1.close();

    for (size_t i = 0; i < splits; ++i) {
        DBReader<unsigned int>::removeDb(fileNames[i].first);
        FileUtil::munmapData(dataFile[i], dataFileSize[i]);
        if (fclose(files[i]) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << fileNames[i].first << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    delete [] dataFile;
    delete [] dataFileSize;
    delete [] files;

    Debug(Debug::INFO) << "Time for merging target splits: " << timer.lap() << "\n";
}


ScoreMatrix Prefiltering::getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize) {
    if (templateDBIsIndex == true) {
        switch(kmerSize) {
            case 2:
                return PrefilteringIndexReader::get2MerScoreMatrix(tidxdbr, preloadMode);
            case 3:
                return PrefilteringIndexReader::get3MerScoreMatrix(tidxdbr, preloadMode);
            default:
                break;
        }
    }
    return ExtendedSubstitutionMatrix::calcScoreMatrix(matrix, kmerSize);

}

void Prefiltering::getIndexTable(int split, size_t dbFrom, size_t dbSize) {
    if (templateDBIsIndex == true) {
        indexTable = PrefilteringIndexReader::getIndexTable(split, tidxdbr, preloadMode);
        // only the ungapped alignment needs the sequence lookup, we can save quite some memory here
        if (diagonalScoring) {
            sequenceLookup = PrefilteringIndexReader::getSequenceLookup(split, tidxdbr, preloadMode);
        }
    } else {
        Timer timer;

        Sequence tseq(maxSeqLen, targetSeqType, kmerSubMat, kmerSize, spacedKmer, aaBiasCorrection, true, spacedKmerPattern);
        int localKmerThr = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) ||
                            Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES) ||
                            (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE) == false && targetSearchMode == 0 && takeOnlyBestKmer == true) ) ? 0 : kmerThr;

        // remove X or N for seeding
        int adjustAlphabetSize = (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES) ||
                                  Parameters::isEqualDbtype(targetSeqType,Parameters::DBTYPE_AMINO_ACIDS))
                                 ? alphabetSize -1 : alphabetSize;
        indexTable = new IndexTable(adjustAlphabetSize, kmerSize, false);
        SequenceLookup **unmaskedLookup = maskMode == 0 && maskLowerCaseMode == 0 ? &sequenceLookup : NULL;
        SequenceLookup **maskedLookup   = maskMode == 1 || maskLowerCaseMode == 1 ? &sequenceLookup : NULL;

        Debug(Debug::INFO) << "Index table k-mer threshold: " << localKmerThr << " at k-mer size " << kmerSize << " \n";
        IndexBuilder::fillDatabase(indexTable, maskedLookup, unmaskedLookup, *kmerSubMat,
                                   _3merSubMatrix, _2merSubMatrix,
                                   &tseq, tdbr, dbFrom, dbFrom + dbSize,
                                   localKmerThr, maskMode, maskLowerCaseMode, maskProb, targetSearchMode);

        // sequenceLookup has to be temporarily present to speed up masking
        // afterwards its not needed anymore without diagonal scoring
        if (diagonalScoring == false) {
            delete sequenceLookup;
            sequenceLookup = NULL;
        }

        indexTable->printStatistics(kmerSubMat->num2aa);
        tdbr->remapData();
        Debug(Debug::INFO) << "Time for index table init: " << timer.lap() << "\n";
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
void Prefiltering::runMpiSplits(const std::string &resultDB, const std::string &resultDBIndex, const std::string &localTmpPath, const int runRandomId) {
    if(compressed == true && splitMode == Parameters::TARGET_DB_SPLIT){
            Debug(Debug::WARNING) << "The output of the prefilter cannot be compressed during target split mode. "
                                     "Prefilter result will not be compressed.\n";
            compressed = false;
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

    std::pair<std::string, std::string> result = Util::createTmpFileNames(procTmpResultDB, procTmpResultDBIndex, MMseqsMPI::rank + runRandomId);
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
            mergePrefilterSplits(resultDB, resultDBIndex, splitFiles);
        } else {
            DBWriter writer(resultDB.c_str(), resultDBIndex.c_str(), 1, compressed, Parameters::DBTYPE_PREFILTER_RES);
            writer.open();
            writer.close();
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
            Debug(Debug::WARNING) << "The output of the prefilter cannot be compressed during target split mode. "
                                     "Prefilter result will not be compressed.\n";
            compressed = false;
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
            mergePrefilterSplits(resultDB, resultDBIndex, splitFiles);
            if (splitFiles.size() > 1) {
                DBReader<unsigned int> resultReader(resultDB.c_str(), resultDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
                resultReader.open(DBReader<unsigned int>::NOSORT);
                resultReader.readMmapedDataInMemory();
                const std::pair<std::string, std::string> tempDb = Util::databaseNames(resultDB + "_tmp");
                DBWriter resultWriter(tempDb.first.c_str(), tempDb.second.c_str(), threads, compressed, Parameters::DBTYPE_PREFILTER_RES);
                resultWriter.open();
                resultWriter.sortDatafileByIdOrder(resultReader);
                resultWriter.close(true);
                resultReader.close();
                DBReader<unsigned int>::removeDb(resultDB);
                DBReader<unsigned int>::moveDb(tempDb.first, resultDB);
            }
            hasResult = true;
        }
    } else if (splitProcessCount == 1) {
        if (runSplit(resultDB.c_str(), resultDBIndex.c_str(), fromSplit, merge)) {
            hasResult = true;
        }
    } else if (splitProcessCount == 0) {
        DBWriter writer(resultDB.c_str(), resultDBIndex.c_str(), 1, compressed, Parameters::DBTYPE_PREFILTER_RES);
        writer.open();
        writer.close();
        hasResult = false;
    }

    return hasResult;
}

bool Prefiltering::runSplit(const std::string &resultDB, const std::string &resultDBIndex, size_t split, bool merge) {
    Debug(Debug::INFO) << "Process prefiltering step " << (split + 1) << " of " << splits << "\n\n";

    size_t dbFrom = 0;
    size_t dbSize = tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();

    // create index table based on split parameter
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        tdbr->decomposeDomainByAminoAcid(split, splits, &dbFrom, &dbSize);
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
        qdbr->decomposeDomainByAminoAcid(split, splits, &queryFrom, &querySize);
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

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)threads, querySize), (size_t)1);
#endif

    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), localThreads, compressed, Parameters::DBTYPE_PREFILTER_RES);
    tmpDbw.open();

    // init all thread-specific data structures
    char *notEmpty = new char[querySize];
    memset(notEmpty, 0, querySize * sizeof(char)); // init notEmpty

    std::list<int> **reslens = new std::list<int> *[localThreads];
    for (size_t i = 0; i < localThreads; ++i) {
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
        Sequence seq(qdbr->getMaxSeqLen(), querySeqType, kmerSubMat, kmerSize, spacedKmer, aaBiasCorrection, true, spacedKmerPattern);
        QueryMatcher matcher(indexTable, sequenceLookup, kmerSubMat,  ungappedSubMat,
                             kmerThr, kmerSize, dbSize, std::max(tdbr->getMaxSeqLen(),qdbr->getMaxSeqLen()), maxResListLen, aaBiasCorrection, aaBiasCorrectionScale,
                             diagonalScoring, minDiagScoreThr, takeOnlyBestKmer, targetSeqType==Parameters::DBTYPE_NUCLEOTIDES);

        if (seq.profile_matrix != NULL) {
            matcher.setProfileMatrix(seq.profile_matrix);
        } else if (_3merSubMatrix.isValid() && _2merSubMatrix.isValid()) {
            matcher.setSubstitutionMatrix(&_3merSubMatrix, &_2merSubMatrix);
        } else {
            matcher.setSubstitutionMatrix(NULL, NULL);
        }

        if (taxonomyHook != NULL) {
            matcher.setQueryMatcherHook(taxonomyHook);
        }

        char buffer[128];
        std::string result;
        result.reserve(1000000);

#pragma omp for schedule(dynamic, 1) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
        for (size_t id = queryFrom; id < queryFrom + querySize; id++) {
            progress.updateProgress();
            // get query sequence
            char *seqData = qdbr->getData(id, thread_idx);
            unsigned int qKey = qdbr->getDbKey(id);
            seq.mapSequence(id, qKey, seqData, qdbr->getSeqLen(id));
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
            if (taxonomyHook != NULL) {
                taxonomyHook->setDbFrom(dbFrom);
            }
            std::pair<hit_t *, size_t> prefResults = matcher.matchQuery(&seq, targetSeqId, targetSeqType==Parameters::DBTYPE_NUCLEOTIDES);
            size_t resultSize = prefResults.second;
            const float queryLength = static_cast<float>(qdbr->getSeqLen(id));
            for (size_t i = 0; i < resultSize; i++) {
                hit_t *res = prefResults.first + i;
                // correct the 0 indexed sequence id again to its real identifier
                size_t targetSeqId1 = res->seqId + dbFrom;
                // replace id with key
                res->seqId = tdbr->getDbKey(targetSeqId1);
                if (UNLIKELY(targetSeqId1 >= tdbr->getSize())) {
                    Debug(Debug::WARNING) << "Wrong prefiltering result for query: " << qdbr->getDbKey(id) << " -> " << targetSeqId1 << "\t" << res->prefScore << "\n";
                }

                // TODO: check if this should happen when diagonalScoring == false
                if (covThr > 0.0 && (covMode == Parameters::COV_MODE_BIDIRECTIONAL
                                               || covMode == Parameters::COV_MODE_QUERY
                                               || covMode == Parameters::COV_MODE_LENGTH_SHORTER )) {
                    const float targetLength = static_cast<float>(tdbr->getSeqLen(targetSeqId1));
                    if (Util::canBeCovered(covThr, covMode, queryLength, targetLength) == false) {
                        continue;
                    }
                }

                // write prefiltering results to a string
                int len = QueryMatcher::prefilterHitToBuffer(buffer, *res);
                result.append(buffer, len);
            }
            tmpDbw.writeData(result.c_str(), result.length(), qKey, thread_idx);
            result.clear();

            // update statistics counters
            if (resultSize != 0) {
                notEmpty[id - queryFrom] = 1;
            }

            if (Debug::debugLevel >= Debug::INFO) {
                kmersPerPos += matcher.getStatistics()->kmersPerPos;
                dbMatches += matcher.getStatistics()->dbMatches;
                doubleMatches += matcher.getStatistics()->doubleMatches;
                querySeqLenSum += seq.L;
                diagonalOverflow += matcher.getStatistics()->diagonalOverflow;
                resSize += resultSize;
                realResSize += std::min(resultSize, maxResListLen);
                reslens[thread_idx]->emplace_back(resultSize);
            }
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
        // free memory early since the merge might need quite a bit of memory
        if (indexTable != NULL) {
            delete indexTable;
            indexTable = NULL;
        }
        if (sequenceLookup != NULL) {
            delete sequenceLookup;
            sequenceLookup = NULL;
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

    for (size_t i = 0; i < localThreads; i++) {
        reslens[i]->clear();
        delete reslens[i];
    }
    delete[] reslens;
    delete[] notEmpty;

    return true;
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


BaseMatrix *Prefiltering::getSubstitutionMatrix(const MultiParam<NuclAA<std::string>> &scoringMatrixFile, MultiParam<NuclAA<int>> alphabetSize, float bitFactor, bool profileState, bool isNucl) {
    BaseMatrix *subMat;

    if (isNucl){
        subMat = new NucleotideMatrix(scoringMatrixFile.values.nucleotide().c_str(), bitFactor, 0.0);
    } else if (alphabetSize.values.aminoacid() < 21) {
        SubstitutionMatrix sMat(scoringMatrixFile.values.aminoacid().c_str(), bitFactor, -0.2f);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, alphabetSize.values.aminoacid(), bitFactor);
    }else if(profileState == true){
        SubstitutionMatrix sMat(scoringMatrixFile.values.aminoacid().c_str(), bitFactor, -0.2f);
        subMat = new SubstitutionMatrixProfileStates(sMat.matrixName, sMat.probMatrix, sMat.pBack,
                                                     sMat.subMatrixPseudoCounts, bitFactor, 0.0, 8);
    } else {
        subMat = new SubstitutionMatrix(scoringMatrixFile.values.aminoacid().c_str(), bitFactor, -0.2f);
    }
    return subMat;
}

void Prefiltering::mergePrefilterSplits(const std::string &outDB, const std::string &outDBIndex,
                              const std::vector<std::pair<std::string, std::string>> &splitFiles) {
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        mergeTargetSplits(outDB, outDBIndex, splitFiles, threads);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        DBWriter::mergeResults(outDB, outDBIndex, splitFiles);
    }
}

int Prefiltering::getKmerThreshold(const float sensitivity, const bool isProfile, const bool hasContextPseudoCnts, const SeqProf<int> kmerScore, const int kmerSize) {
    if (isProfile == true && kmerScore.profile() != INT_MAX) {
        return kmerScore.profile();
    }
    if (isProfile == false && kmerScore.sequence() != INT_MAX) {
        return kmerScore.sequence();
    }
    float kmerThrBest = FLT_MAX;
    int paramType = isProfile ? Parameters::DBTYPE_HMM_PROFILE : Parameters::DBTYPE_AMINO_ACIDS;
    for(size_t i = 0; i < externalThreshold.size(); i++){
        if(kmerSize == externalThreshold[i].kmerSize && externalThreshold[i].sequenceType == paramType){
            return static_cast<int>(externalThreshold[i].base - (externalThreshold[i].sensPerStep * sensitivity));
        }
    }
    if (isProfile == true) {
        if (hasContextPseudoCnts == true) {
            if (kmerSize == 5) {
                float base = 97.75;
                kmerThrBest = base - (sensitivity * 8.75);
            } else if (kmerSize == 6) {
                float base = 132.75;
                kmerThrBest = base - (sensitivity * 8.75);
            } else if (kmerSize == 7) {
                float base = 158.75;
                kmerThrBest = base - (sensitivity * 9.75);
            } else {
                Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            if (kmerSize == 5) {
                float base = 108.8;
                kmerThrBest = base - (sensitivity * 4.7);
            } else if (kmerSize == 6) {
                float base = 134.35;
                kmerThrBest = base - (sensitivity * 6.15);
            } else if (kmerSize == 7) {
                float base = 149.15;
                kmerThrBest = base - (sensitivity * 6.85);
            } else {
                Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid\n";
                EXIT(EXIT_FAILURE);
            }
        }
    } else {
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
            Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid\n";
            EXIT(EXIT_FAILURE);
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
    size_t indexTableSize = static_cast<size_t>(pow(alphabetSize, kmerSize)) * sizeof(size_t);
    // memory needed for the threads
    // This memory is an approx. for Countint32Array and QueryTemplateLocalFast
    size_t threadSize = threads * (
            (dbSizeSplit * 2 * sizeof(IndexEntryLocal)) // databaseHits in QueryMatcher
            + (dbSizeSplit * 1.5 * sizeof(CounterResult)) // databaseHits in QueryMatcher
            // 1.5 is a security factor
            + (maxResListLen * sizeof(hit_t))
            + (dbSizeSplit * 2 * sizeof(CounterResult) * 2) // BINS * binSize, (binSize = dbSize * 2 / BINS)
              // 2 is a security factor the size can increase during run
    );
    size_t dbReaderSize = dbSize * (sizeof(DBReader<unsigned int>::Index) + sizeof(unsigned int)); // DB index size

    // extended matrix
    size_t extendedMatrix = 0;
    if(Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)){
        extendedMatrix = sizeof(std::pair<short, unsigned int>) * static_cast<size_t>(pow(pow(alphabetSize, 3), 2));
        extendedMatrix += sizeof(std::pair<short, unsigned int>) * pow(pow(alphabetSize, 2), 2);
    }
    // some memory needed to keep the index, ....
    size_t background = dbSize * 22;
    // return result in bytes
    return residueSize + indexTableSize + threadSize + background + extendedMatrix + dbReaderSize;
}

size_t Prefiltering::estimateHDDMemoryConsumption(size_t dbSize, size_t maxResListLen) {
    // 21 bytes is roughly the size of an entry
    // 2x because the merge doubles the hdd demand
    return 2 * (21 * dbSize * maxResListLen);
}

std::pair<int, int> Prefiltering::optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr,
                                                int alphabetSize, int externalKmerSize, unsigned int querySeqType, unsigned int threads) {

    int startKmerSize = (externalKmerSize == 0) ? 6 : externalKmerSize;
    int endKmerSize   = (externalKmerSize == 0) ? 7 : externalKmerSize;

    if(Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        startKmerSize = (externalKmerSize == 0) ? 14 : externalKmerSize;
        endKmerSize   = (externalKmerSize == 0) ? 15 : externalKmerSize;
    }

    for (int optKmerSize = endKmerSize; optKmerSize >= startKmerSize ; optKmerSize--) {
        size_t aaUpperBoundForKmerSize = (SIZE_MAX - 1);
        if(externalKmerSize == 0){
            if(Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
                aaUpperBoundForKmerSize = IndexTable::getUpperBoundNucCountForKmerSize(optKmerSize);
            }else{
                aaUpperBoundForKmerSize = IndexTable::getUpperBoundAACountForKmerSize(optKmerSize);
            }
        }
        for (int optSplit = 1; optSplit < 1000; optSplit++) {
            if ((tdbr->getAminoAcidDBSize() / optSplit) < aaUpperBoundForKmerSize) {
                size_t neededSize = estimateMemoryConsumption(optSplit, tdbr->getSize(),
                                                              tdbr->getAminoAcidDBSize(),
                                                              0, alphabetSize, optKmerSize, querySeqType,
                                                              threads);
                if (neededSize < 0.9 * totalMemoryInByte) {
                    return std::make_pair(optKmerSize, optSplit);
                }
            }
        }
    }

    return std::make_pair(-1, -1);
}



