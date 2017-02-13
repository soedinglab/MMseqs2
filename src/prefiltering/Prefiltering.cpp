#include "Prefiltering.h"
#include "NucleotideMatrix.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrixWithoutX.h"
#include "ExtendedSubstitutionMatrix.h"

#include <regex.h>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

Prefiltering::Prefiltering(const std::string &targetDB,
                           const std::string &targetDBIndex,
                           const Parameters &par) :
        split(par.split), targetDB(targetDB), targetDBIndex(targetDBIndex),
        kmerSize(par.kmerSize), spacedKmer(par.spacedKmer != 0), alphabetSize(par.alphabetSize),
        maskResidues(par.maskResidues != 0), splitMode(par.splitMode), scoringMatrixFile(par.scoringMatrixFile),
        maxResListLen(par.maxResListLen), kmerScore(par.kmerScore),
        sensitivity(par.sensitivity), resListOffset(par.resListOffset), maxSeqLen(par.maxSeqLen),
        querySeqType(par.querySeqType), targetSeqType(par.targetSeqType),
        diagonalScoring(par.diagonalScoring != 0), minDiagScoreThr(static_cast<unsigned int>(par.minDiagScoreThr)),
        aaBiasCorrection(par.compBiasCorrection != 0), covThr(par.covThr), includeIdentical(par.includeIdentity),
        threads(static_cast<unsigned int>(par.threads)) {
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif

    std::string indexDB = PrefilteringIndexReader::searchForIndex(targetDB);
    if (indexDB != "") {
        Debug(Debug::INFO) << "Use index  " << indexDB << "\n";
        tdbr = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tdbr);
        if (templateDBIsIndex == true) {
            // exchange reader with old ffindex reader
            tidxdbr = tdbr;
            tdbr = PrefilteringIndexReader::openNewReader(tdbr);
            PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tidxdbr);
            kmerSize = data.kmerSize;
            alphabetSize = data.alphabetSize;
            maskResidues = data.maskResidues != 0;
            split = data.split;
            spacedKmer = data.spacedKmer != 0;
            scoringMatrixFile = PrefilteringIndexReader::getSubstitutionMatrixName(tidxdbr);
        } else {
            Debug(Debug::ERROR) << "Outdated index version. Please recompute it with 'createindex'!\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        Debug(Debug::INFO) << "Could not find precomputed index. Compute index.\n";
        tdbr = new DBReader<unsigned int>(targetDB.c_str(), targetDBIndex.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);

        templateDBIsIndex = false;
    }

    if (par.noPreload == false) {
        tdbr->readMmapedDataInMemory();
        tdbr->mlock();
    }

    Debug(Debug::INFO) << "Target database: " << targetDB << "(size=" << tdbr->getSize() << ")\n";

    setupSplit(*tdbr, alphabetSize, threads, templateDBIsIndex, &kmerSize, &split, &splitMode);

    if (split > 1 && splitMode == Parameters::TARGET_DB_SPLIT) {
        maxResListLen = (maxResListLen / split) + 1;
    }

    // init the substitution matrices
    switch (querySeqType) {
        case Sequence::NUCLEOTIDES:
            subMat = new NucleotideMatrix();
            _2merSubMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
            _3merSubMatrix = getScoreMatrix(*subMat, 3);
            break;
        case Sequence::AMINO_ACIDS:
            subMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0, false);
            alphabetSize = subMat->alphabetSize;
            _2merSubMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
            _3merSubMatrix = getScoreMatrix(*subMat, 3);
            break;
        case Sequence::HMM_PROFILE:
            // needed for Background distributions
            subMat = getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0, false);
            _2merSubMatrix = NULL;
            _3merSubMatrix = NULL;
            break;
        default:
            Debug(Debug::ERROR) << "Query sequence type not implemented!\n";
            EXIT(EXIT_FAILURE);
    }

    if (splitMode == Parameters::QUERY_DB_SPLIT) {
        // create the whole index table
        const int kmerThr = getKmerThreshold(sensitivity);
        indexTable = getIndexTable(0, 0, tdbr->getSize(), kmerThr, threads);
    } else {
        indexTable = NULL;
    }
}

Prefiltering::~Prefiltering() {
    if (indexTable != NULL) {
        delete indexTable;
    }

    tdbr->close();
    delete tdbr;

    if (templateDBIsIndex == true) {
        tidxdbr->close();
        delete tidxdbr;
    }

    delete subMat;
    if (_2merSubMatrix != NULL) {
        ScoreMatrix::cleanup(_2merSubMatrix);
    }
    if (_3merSubMatrix != NULL && templateDBIsIndex == false) {
        ScoreMatrix::cleanup(_3merSubMatrix);
    }
}

void Prefiltering::setupSplit(DBReader<unsigned int>& dbr, const int alphabetSize, const int threads, const bool templateDBIsIndex, int *kmerSize, int *split, int *splitMode) {
    const size_t totalMemoryInByte = Util::getTotalSystemMemory();
    size_t neededSize = estimateMemoryConsumption(1,
                                                  dbr.getSize(), dbr.getAminoAcidDBSize(), alphabetSize,
                                                  *kmerSize == 0 ? // if auto detect kmerSize
                                                  IndexTable::computeKmerSize(dbr.getAminoAcidDBSize()) : *kmerSize,
                                                  threads);
    if (neededSize > 0.9 * totalMemoryInByte) { // memory is not enough to compute everything at once
        std::pair<int, int> splitingSetting = Prefiltering::optimizeSplit(totalMemoryInByte, &dbr, alphabetSize, *kmerSize, threads);
        if (splitingSetting.second == -1) {
            Debug(Debug::ERROR) << "Can not fit databased into " << totalMemoryInByte
                                << " byte. Please use a computer with more main memory.\n";
            EXIT(EXIT_FAILURE);
        }
        if (*split == Parameters::AUTO_SPLIT_DETECTION && templateDBIsIndex == false) {
            *split = splitingSetting.second;
            if (kmerSize == 0) {
                // set k-mer based on aa size in database
                // if we have less than 10Mio * 335 amino acids use 6mers
                *kmerSize = splitingSetting.first;
            }
        }
        if (*splitMode == Parameters::DETECT_BEST_DB_SPLIT) {
            *splitMode = Parameters::TARGET_DB_SPLIT;
        }
    } else { // memory is  enough to compute everything with split setting
        if (*kmerSize == 0) {
            const int tmpSplit = (*split > 1) ? *split : 1;
            size_t aaSize = dbr.getAminoAcidDBSize() / tmpSplit;
            *kmerSize = IndexTable::computeKmerSize(aaSize);
        }

        if (*split == Parameters::AUTO_SPLIT_DETECTION) {
            *split = 1;
        }

        if (*splitMode == Parameters::DETECT_BEST_DB_SPLIT) {
            if (templateDBIsIndex == true && *split > 1) {
                *splitMode = Parameters::TARGET_DB_SPLIT;
            } else {
#ifdef HAVE_MPI
                *splitMode = Parameters::QUERY_DB_SPLIT;
#else
                *splitMode = Parameters::TARGET_DB_SPLIT;
#endif
            }
        }
    }

    Debug(Debug::INFO) << "Use kmer size " << *kmerSize << " and split "
                       << *split << " using split mode " << *splitMode << "\n";
    neededSize = estimateMemoryConsumption((*splitMode == Parameters::TARGET_DB_SPLIT) ? *split : 1, dbr.getSize(),
                                           dbr.getAminoAcidDBSize(), alphabetSize, *kmerSize, threads);
    Debug(Debug::INFO) << "Needed memory (" << neededSize << " byte) of total memory (" << totalMemoryInByte
                       << " byte)\n";
    if (neededSize > 0.9 * totalMemoryInByte) {
        Debug(Debug::WARNING) << "WARNING: MMseqs processes needs more main memory than available."
                "Increase the size of --split or set it to 0 to automatic optimize target database split.\n";
        if (templateDBIsIndex == true) {
            Debug(Debug::WARNING) << "WARNING: Split has to be computed by createindex if precomputed index is used.\n";
        }
    }
}

void Prefiltering::mergeOutput(const std::string &outDB, const std::string &outDBIndex,
                               const std::vector<std::pair<std::string, std::string>> &filenames) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    if (filenames.size() < 2) {
        std::rename(filenames[0].first.c_str(), outDB.c_str());
        std::rename(filenames[0].second.c_str(), outDBIndex.c_str());
        Debug(Debug::INFO) << "No merging needed.\n";
        return;
    }
    std::list<std::pair<std::string, std::string>> files(filenames.begin(), filenames.end());
    size_t mergeStep = 0;
    while (true) {
        std::pair<std::string, std::string> file1 = files.front();
        files.pop_front();
        std::pair<std::string, std::string> file2 = files.front();
        files.pop_front();
        std::pair<std::string, std::string> out = std::make_pair((outDB + "_merge_" + SSTR(mergeStep)),
                                                                 (outDBIndex + "_merge_" + SSTR(mergeStep)));
        DBWriter writer(out.first.c_str(), out.second.c_str(), 1);
        writer.open(1024 * 1024 * 1024); // 1 GB buffer
        writer.mergeFilePair(file1.first.c_str(), file1.second.c_str(), file2.first.c_str(), file2.second.c_str());
        // remove split
        int error = 0;
        error += remove(file1.first.c_str()); error += remove(file1.second.c_str());
        error += remove(file2.first.c_str()); error += remove(file2.second.c_str());
        if(error != 0){
            Debug(Debug::ERROR) << "Error while deleting files in mergeOutput!\n";
            EXIT(EXIT_FAILURE);
        }
        writer.close();
        // push back the current merge to result to the end
        files.push_back(out);
        mergeStep++;
        if (files.size() == 1) {
            break;
        }
    }
    const std::pair<std::string, std::string> &out = files.front();
    Debug(Debug::INFO) << out.first << " " << out.second << "\n";

    std::rename(out.first.c_str(), outDB.c_str());
    std::rename(out.second.c_str(), outDBIndex.c_str());
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for merging results: " << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
}

QueryMatcher **Prefiltering::createQueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable, Sequence** qseq,
                                                        unsigned int *seqLens, short kmerThr,
                                                        double kmerMatchProb, int kmerSize,
                                                        unsigned int effectiveKmerSize, size_t dbSize,
                                                        bool aaBiasCorrection, bool diagonalScoring,
                                                        unsigned int maxSeqLen, size_t maxHitsPerQuery) {
    QueryMatcher **matchers = new QueryMatcher *[threads];
#pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < threads; i++) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        matchers[thread_idx] = new QueryMatcher(m, indexTable, seqLens, kmerThr, kmerMatchProb,
                                                kmerSize, dbSize, maxSeqLen, effectiveKmerSize,
                                                maxHitsPerQuery, aaBiasCorrection, diagonalScoring, minDiagScoreThr);
        if (querySeqType == Sequence::HMM_PROFILE) {
            matchers[thread_idx]->setProfileMatrix(qseq[thread_idx]->profile_matrix);
        } else {
            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix, _2merSubMatrix);
        }
    }
    return matchers;
}

ScoreMatrix *Prefiltering::getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize) {
    if (templateDBIsIndex == true) {
        switch(kmerSize) {
//            case 2:
//                return PrefilteringIndexReader::get2MerScoreMatrix(tidxdbr);
            case 3:
                return PrefilteringIndexReader::get3MerScoreMatrix(tidxdbr);
            default:
                Debug(Debug::ERROR) << "Invalid k-mer score matrix!\n";
                EXIT(EXIT_FAILURE);
        }
    } else {
        return ExtendedSubstitutionMatrix::calcScoreMatrix(matrix, kmerSize);
    }
}

IndexTable *Prefiltering::getIndexTable(int split, size_t dbFrom, size_t dbSize, int kmerThr, unsigned int threads) {
    if (templateDBIsIndex == true) {
        return PrefilteringIndexReader::generateIndexTable(tidxdbr, split, diagonalScoring);
    } else {
        Sequence tseq(maxSeqLen, subMat->aa2int, subMat->int2aa, targetSeqType, kmerSize, spacedKmer, aaBiasCorrection);
        return generateIndexTable(tdbr, &tseq, subMat, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize,
                                  diagonalScoring, maskResidues, kmerThr, threads);
    }
}
void Prefiltering::runAllSplits(const std::string &queryDB, const std::string &queryDBIndex,
                                const std::string &resultDB, const std::string &resultDBIndex) {
    runSplits(queryDB, queryDBIndex, resultDB, resultDBIndex, 0, split);
}

void Prefiltering::runSplits(const std::string &queryDB, const std::string &queryDBIndex,
                             const std::string &resultDB, const std::string &resultDBIndex,
                             size_t fromSplit, size_t splits) {
//    FileUtil::errorIfFileExist(outDB.c_str());
//    FileUtil::errorIfFileExist(outDBIndex.c_str());

    //  check if when qdb and tdb have the same name an index extension exists
    std::string check(targetDB);
    size_t pos = check.find(queryDB);
    int nomatch = true;
    if (pos == 0) {
        check.replace(0, queryDB.length(), "");
        regex_t regex;
        regcomp(&regex, "^\\.s?k[5-7]$", REG_EXTENDED | REG_NEWLINE);
        nomatch = regexec(&regex, check.c_str(), 0, NULL, 0);
        regfree(&regex);
    }
    // if no match found or two matches found (we want exactly one match)
    bool sameQTDB = (queryDB.compare(targetDB) == 0 || (nomatch == false));
    DBReader<unsigned int> *qdbr;
    if (templateDBIsIndex == false && sameQTDB == true) {
        qdbr = tdbr;
    } else {
        qdbr = new DBReader<unsigned int>(queryDB.c_str(), queryDBIndex.c_str());
        qdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }
    Debug(Debug::INFO) << "Query database: " << queryDB << "(size=" << qdbr->getSize() << ")\n";

    // splits template database into x sequence steps
    std::vector<std::pair<std::string, std::string> > splitFiles;
    for (size_t s = fromSplit; s < (fromSplit + splits); s++) {
        std::pair<std::string, std::string> filenamePair = Util::createTmpFileNames(resultDB, resultDBIndex, s);
        splitFiles.push_back(filenamePair);
        runSplit(qdbr, filenamePair.first.c_str(), filenamePair.second.c_str(), s, split, sameQTDB);
    }
    // prefiltering scores calculation end

    // merge output ffindex databases
    if (split > 1) {
        mergeFiles(splitFiles, resultDB, resultDBIndex);
    } else {
        std::rename(splitFiles[0].first.c_str(), resultDB.c_str());
        std::rename(splitFiles[0].second.c_str(), resultDBIndex.c_str());
    }

    if (sameQTDB == false) {
        qdbr->close();
        delete qdbr;
    }
}

void Prefiltering::runSplit(DBReader<unsigned int>* qdbr, const std::string &resultDB, const std::string &resultDBIndex,
                            size_t split, size_t splitCount, bool sameQTDB) {

    // init all thread-specific data structures
    Sequence **qseq = new Sequence *[threads];
    std::list<int> **reslens = new std::list<int> *[threads];
    char *notEmpty = new char[qdbr->getSize()];

#pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < threads; i++) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        qseq[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa,
                                        querySeqType, kmerSize, spacedKmer, aaBiasCorrection);
        reslens[thread_idx] = new std::list<int>();
    }

    Debug(Debug::INFO) << "Process prefiltering step " << split << " of " << splitCount << "\n\n";
    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), threads);
    tmpDbw.open();
    size_t dbFrom = 0;
    size_t dbSize = tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();

    const int kmerThr = getKmerThreshold(sensitivity);
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        Util::decomposeDomainByAminoAcid(tdbr->getAminoAcidDBSize(), tdbr->getSeqLens(), tdbr->getSize(),
                                         split, splitCount, &dbFrom, &dbSize);
        //TODO fix this what if we have 10 chunks but only 4 servers (please fix me)
        if (indexTable != NULL) {
            delete indexTable;
        }
        indexTable = getIndexTable(split, dbFrom, dbSize, kmerThr, threads);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        Util::decomposeDomainByAminoAcid(qdbr->getAminoAcidDBSize(), qdbr->getSeqLens(), qdbr->getSize(),
                                         split, splitCount, &queryFrom, &querySize);
    }
    // create index table based on split parameter
    // run small query sample against the index table to calibrate p-match

    const double kmerMatchProb = (diagonalScoring == true) ? 0.0f
                                                           : setKmerThreshold(indexTable, qdbr, qseq, tdbr,
                                                                              kmerThr, qseq[0]->getEffectiveKmerSize());

    Debug(Debug::INFO) << "k-mer similarity threshold: " << kmerThr << "\n";
    Debug(Debug::INFO) << "k-mer match probability: " << kmerMatchProb << "\n\n";

    struct timeval start, end;
    gettimeofday(&start, NULL);
    QueryMatcher **matchers = createQueryTemplateMatcher(subMat, indexTable, qseq,
                                                         tdbr->getSeqLens() + dbFrom, // offset for split mode
                                                         kmerThr, kmerMatchProb, kmerSize,
                                                         qseq[0]->getEffectiveKmerSize(), dbSize,
                                                         aaBiasCorrection, diagonalScoring, maxSeqLen, maxResListLen);
    size_t kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    size_t totalQueryDBSize = querySize;
    memset(notEmpty, 0, qdbr->getSize() * sizeof(char)); // init notEmpty

    Debug(Debug::INFO) << "Starting prefiltering scores calculation (step " << split << " of " << splitCount << ")\n";
    Debug(Debug::INFO) << "Query db start  " << queryFrom << " to " << queryFrom + querySize << "\n";
    Debug(Debug::INFO) << "Target db start  " << dbFrom << " to " << dbFrom + dbSize << "\n";

#pragma omp parallel for schedule(dynamic, 10) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
    for (size_t id = queryFrom; id < queryFrom + querySize; id++) {
        Debug::printProgress(id);

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        // get query sequence
        char *seqData = qdbr->getData(id);
        unsigned int qKey = qdbr->getDbKey(id);
        qseq[thread_idx]->mapSequence(id, qKey, seqData);
        // only the corresponding split should include the id (hack for the hack)
        size_t targetSeqId = UINT_MAX;
        if (id >= dbFrom && id < (dbFrom + dbSize) && (sameQTDB || includeIdentical)) {
            targetSeqId = tdbr->getId(qseq[thread_idx]->getDbKey());
            if (targetSeqId != UINT_MAX) {
                targetSeqId = targetSeqId - dbFrom;
            }
        }
        // calculate prefiltering results
        std::pair<hit_t *, size_t> prefResults = matchers[thread_idx]->matchQuery(qseq[thread_idx], targetSeqId);
        size_t resultSize = prefResults.second;
        // write
        writePrefilterOutput(qdbr, &tmpDbw, thread_idx, id, prefResults, dbFrom, diagonalScoring, resListOffset);

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
        reslens[thread_idx]->emplace_back(resultSize);
    } // step end
    statistics_t stats(kmersPerPos / totalQueryDBSize, dbMatches / totalQueryDBSize, doubleMatches / totalQueryDBSize,
                       querySeqLenSum, diagonalOverflow, resSize / totalQueryDBSize);
    size_t empty = 0;
    for (size_t id = queryFrom; id < queryFrom + querySize; id++) {
        if (notEmpty[id] == 0) {
            empty++;
        }
    }
    printStatistics(stats, reslens, empty);

    if (totalQueryDBSize > 1000) {
        Debug(Debug::INFO) << "\n";
    }
    Debug(Debug::INFO) << "\n";

    for (unsigned int j = 0; j < threads; j++) {
        delete matchers[j];
    }
    delete[] matchers;

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for prefiltering scores calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60)
                       << " m " << (sec % 60) << "s\n";
    tmpDbw.close(); // sorts the index

    // sort by ids
    // needed to speed up merge later one
    // sorts this datafile according to the index file
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

    for (unsigned int i = 0; i < threads; i++) {
        delete qseq[i];
        reslens[i]->clear();
        delete reslens[i];
    }
    delete[] qseq;
    delete[] reslens;
    delete[] notEmpty;
}

// write prefiltering to ffindex database
void Prefiltering::writePrefilterOutput(DBReader<unsigned int> *qdbr, DBWriter *dbWriter, unsigned int thread_idx, size_t id,
                                        const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                                        bool diagonalScoring, size_t resultOffsetPos) {
    // write prefiltering results to a string
    size_t l = 0;
    hit_t *resultVector = prefResults.first + resultOffsetPos;
    const size_t resultSize = (prefResults.second < resultOffsetPos) ? 0 : prefResults.second - resultOffsetPos;
    std::string prefResultsOutString;
    prefResultsOutString.reserve(BUFFER_SIZE);
    char buffer[100];
    for (size_t i = 0; i < resultSize; i++) {
        hit_t *res = resultVector + i;
        size_t targetSeqId = res->seqId + seqIdOffset;
        if (targetSeqId >= tdbr->getSize()) {
            Debug(Debug::INFO) << "Wrong prefiltering result: Query: " << qdbr->getDbKey(id) << " -> " << targetSeqId
                               << "\t" << res->prefScore << "\n";
        }
        if (covThr > 0.0) {
            // check if the sequences could pass the coverage threshold
            float queryLength = static_cast<float>(qdbr->getSeqLens(id));
            float targetLength = static_cast<float>(tdbr->getSeqLens(targetSeqId));
            if ((queryLength / targetLength < covThr) || (targetLength / queryLength < covThr))
                continue;

        }

        int len;
        if (diagonalScoring == true) {
            len = snprintf(buffer, 100, "%s\t%.2e\t%d\n", SSTR(tdbr->getDbKey(targetSeqId)).c_str(),
                           res->pScore, (short) res->diagonal);
        } else {
            len = snprintf(buffer, 100, "%s\t%.4f\t%d\n", SSTR(tdbr->getDbKey(targetSeqId)).c_str(),
                           res->pScore, res->prefScore);
        }
        // TODO: error handling for len
        prefResultsOutString.append(buffer, len);
        l++;
        // maximum allowed result list length is reached
        if (l >= maxResListLen)
            break;
    }
    // write prefiltering results string to ffindex database
    const size_t prefResultsLength = prefResultsOutString.length();
    char *prefResultsOutData = (char *) prefResultsOutString.c_str();
    dbWriter->writeData(prefResultsOutData, prefResultsLength, SSTR(qdbr->getDbKey(id)).c_str(), thread_idx);
}

void Prefiltering::printStatistics(const statistics_t &stats, std::list<int> **reslens, size_t empty) {
    // sort and merge the result list lengths (for median calculation)
    reslens[0]->sort();
    for (unsigned int i = 1; i < threads; i++) {
        reslens[i]->sort();
        reslens[0]->merge(*reslens[i]);
    }
    Debug(Debug::INFO) << "\n" << stats.kmersPerPos << " k-mers per position.\n";
    Debug(Debug::INFO) << stats.dbMatches << " DB matches per sequence.\n";
    Debug(Debug::INFO) << stats.diagonalOverflow << " Overflows.\n";
    Debug(Debug::INFO) << stats.resultsPassedPrefPerSeq << " sequences passed prefiltering per query sequence";
    if (stats.resultsPassedPrefPerSeq > maxResListLen)
        Debug(Debug::INFO) << " (ATTENTION: max. " << maxResListLen
                           << " best scoring sequences were written to the output prefiltering database).\n";
    else
        Debug(Debug::INFO) << ".\n";
    size_t mid = reslens[0]->size() / 2;
    std::list<int>::iterator it = reslens[0]->begin();
    std::advance(it, mid);
    Debug(Debug::INFO) << "Median result list size: " << *it << "\n";
    Debug(Debug::INFO) << empty << " sequences with 0 size result lists.\n";
}

BaseMatrix *Prefiltering::getSubstitutionMatrix(const std::string &scoringMatrixFile, size_t alphabetSize,
                                                float bitFactor, bool ignoreX) {
    Debug(Debug::INFO) << "Substitution matrices...\n";
    BaseMatrix *subMat;
    if (alphabetSize < 21) {
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2f);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, alphabetSize, bitFactor);
    } else if (ignoreX == true) {
        SubstitutionMatrix sMat(scoringMatrixFile.c_str(), bitFactor, -0.2f);
        subMat = new SubstitutionMatrixWithoutX(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.subMatrix, bitFactor);
    } else {
        subMat = new SubstitutionMatrix(scoringMatrixFile.c_str(), bitFactor, -0.2f);
    }
    return subMat;
}


void Prefiltering::fillDatabase(DBReader<unsigned int> *dbr, Sequence *seq, IndexTable *indexTable,
                                BaseMatrix *subMat, size_t dbFrom, size_t dbTo,
                                bool diagonalScoring, bool maskResiudes, int kmerThr, unsigned int threads) {
    Debug(Debug::INFO) << "Index table: counting k-mers...\n";
    // fill and init the index table
    size_t aaCount = 0;
    dbTo = std::min(dbTo, dbr->getSize());
    size_t maskedResidues = 0;
    size_t totalKmerCount = 0;
    size_t tableSize = 0;

    size_t dbSize = dbTo - dbFrom;
    size_t *sequenceOffSet = new size_t[dbSize];
    char *idScoreLookup = new char[subMat->alphabetSize];
    for (int aa = 0; aa < subMat->alphabetSize; aa++){
        short score = subMat->subMatrix[aa][aa];
        if (score > CHAR_MAX || score < CHAR_MIN) {
            Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
        }
        idScoreLookup[aa] = (char) score;
    }

    size_t aaDbSize = 0;
    sequenceOffSet[0] = 0;
    for (size_t id = dbFrom; id < dbTo; id++) {
        int seqLen = std::max(static_cast<int>(dbr->getSeqLens(id)) - 2, 0);
        aaDbSize += seqLen; // remove /n and /0
        size_t idFromNull = (id - dbFrom);
        if (id < dbTo - 1) {
            sequenceOffSet[idFromNull + 1] = sequenceOffSet[idFromNull] + seqLen;
        }
        if (Util::overlappingKmers(seqLen, seq->getEffectiveKmerSize() > 0)) {
            tableSize += 1;
        }
    }

    SequenceLookup *sequenceLookup = new SequenceLookup(dbSize, aaDbSize);
#pragma omp parallel
    {
        Indexer idxer(static_cast<unsigned int>(subMat->alphabetSize), seq->getKmerSize());
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
        unsigned int * buffer = new unsigned int[seq->getMaxLen()];
#pragma omp for schedule(dynamic, 100) reduction(+:aaCount, totalKmerCount, maskedResidues)
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            Debug::printProgress(id - dbFrom);
            char *seqData = dbr->getData(id);
            unsigned int qKey = dbr->getDbKey(id);
            s.mapSequence(id - dbFrom, qKey, seqData);
            maskedResidues += Util::maskLowComplexity(subMat, &s, s.L, 12, 3,
                                                      indexTable->getAlphabetSize(), seq->aa2int[(unsigned char) 'X'],
                                                      true, true, true, true);

            aaCount += s.L;
            totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
            sequenceLookup->addSequence(&s, id - dbFrom, sequenceOffSet[id - dbFrom]);
        }
        delete [] buffer;
    }
    delete[] sequenceOffSet;
    dbr->remapData();
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Index table: Masked residues: " << maskedResidues << "\n";
    if(totalKmerCount == 0) {
        Debug(Debug::ERROR) << "No k-mer could be extracted for the database " << dbr->getDataFileName() << ".\n"
                            << "Maybe the sequences length is less than 14 residues";
        Debug(Debug::ERROR) << "\n";
        if (maskedResidues == true){
            Debug(Debug::ERROR) <<  " or contains only low complexity regions.";
            Debug(Debug::ERROR) << "Use --do-not-mask to deactivate low complexity filter.\n";
        }
        EXIT(EXIT_FAILURE);
    }
    //TODO find smart way to remove extrem k-mers without harming huge protein families
//    size_t lowSelectiveResidues = 0;
//    const float dbSize = static_cast<float>(dbTo - dbFrom);
//    for(size_t kmerIdx = 0; kmerIdx < indexTable->getTableSize(); kmerIdx++){
//        size_t res = (size_t) indexTable->getOffset(kmerIdx);
//        float selectivityOfKmer = (static_cast<float>(res)/dbSize);
//        if(selectivityOfKmer > 0.005){
//            indexTable->getOffset()[kmerIdx] = 0;
//            lowSelectiveResidues += res;
//        }
//    }
//    Debug(Debug::INFO) << "Index table: Remove "<< lowSelectiveResidues <<" none selective residues\n";
//    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";
    size_t tableEntriesNum = 0;
    for (size_t i = 0; i < indexTable->getTableSize(); i++) {
        tableEntriesNum += (size_t) indexTable->getOffset(i);
    }

#pragma omp parallel
    {
        Sequence s(seq->getMaxLen(), seq->aa2int, seq->int2aa,
                   seq->getSeqType(), seq->getKmerSize(), seq->isSpaced(), false);
        Indexer idxer(static_cast<unsigned int>(subMat->alphabetSize), seq->getKmerSize());
        unsigned int thread_idx = 0;
        IndexEntryLocalTmp * buffer = new IndexEntryLocalTmp[seq->getMaxLen()];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        size_t threadFrom, threadSize;
        size_t *offsets = indexTable->getOffsets();
        Util::decomposeDomainSize(tableEntriesNum, offsets, indexTable->getTableSize() + 1,
                                  thread_idx, threads, &threadFrom, &threadSize);

//        Debug(Debug::WARNING) << thread_idx << "\t" << threadFrom << "\t" << threadSize << "\n";

#pragma omp barrier
        if (thread_idx == 0) {
            if (diagonalScoring == true) {
                indexTable->initMemory(tableEntriesNum, sequenceLookup, tableSize);
            } else {
                indexTable->initMemory(tableEntriesNum, NULL, tableSize);
            }
            indexTable->init();
            Debug(Debug::INFO) << "Index table: fill...\n";
        }
#pragma omp barrier
        for (size_t id = dbFrom; id < dbTo; id++) {
            s.resetCurrPos();
            if (thread_idx == 0) {
                Debug::printProgress(id - dbFrom);
            }
            //char *seqData = dbr->getData(id);
            //TODO - dbFrom?!?
            unsigned int qKey = dbr->getDbKey(id);
            //seq->mapSequence(id - dbFrom, qKey, seqData);
            s.mapSequence(id - dbFrom, qKey, sequenceLookup->getSequence(id - dbFrom));
//            Util::maskLowComplexity(subMat, seq, seq->L, 12, 3,
//                                    indexTable->getAlphabetSize(), seq->aa2int['X']);
            indexTable->addSequence(&s, &idxer, buffer, threadFrom, threadSize, kmerThr, idScoreLookup);
        }
        delete [] buffer;
    }
    if (diagonalScoring == false) {
        delete sequenceLookup;
    }
    delete [] idScoreLookup;
    if ((dbTo-dbFrom) > 10000) {
        Debug(Debug::INFO) << "\n";
    }
    Debug(Debug::INFO) << "Index table: removing duplicate entries...\n";
    indexTable->revertPointer();
    Debug(Debug::INFO) << "Index table init done.\n\n";

}

IndexTable* Prefiltering::generateIndexTable(DBReader<unsigned int> *dbr, Sequence *seq, BaseMatrix *subMat,
                                             int alphabetSize, int kmerSize, size_t dbFrom, size_t dbTo,
                                             bool diagonalScoring, bool maskResidues, int kmerThr, unsigned int threads) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    IndexTable *indexTable = new IndexTable(alphabetSize, kmerSize, false);
    fillDatabase(dbr, seq, indexTable, subMat, dbFrom, dbTo, diagonalScoring, maskResidues, kmerThr, threads);

    gettimeofday(&end, NULL);
    indexTable->printStatistics(seq->int2aa);
    time_t sec = end.tv_sec - start.tv_sec;
    dbr->remapData();
    Debug(Debug::INFO) << "Time for index table init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m "
                       << (sec % 60) << "s\n\n\n";
    return indexTable;
}

double Prefiltering::setKmerThreshold(IndexTable *indexTable, DBReader<unsigned int> *qdbr, Sequence** qseq,
                                      DBReader<unsigned int> *tdbr, int kmerThr, unsigned int effectiveKmerSize) {
    size_t targetDbSize = indexTable->getSize();
    size_t targetSeqLenSum = 0;

    for (size_t i = 0; i < targetDbSize; i++) {
        targetSeqLenSum += (tdbr->getSeqLens(i) - effectiveKmerSize);
    }

    // generate a small random sequence set for testing
    size_t querySetSize = std::min(qdbr->getSize(), (size_t) 1000);
    unsigned int *querySeqs = new unsigned int[querySetSize];
    srand(1);
    for (size_t i = 0; i < querySetSize; i++) {
        querySeqs[i] = rand() % qdbr->getSize();
    }

    statistics_t stats = computeStatisticForKmerThreshold(qdbr, qseq, indexTable, querySetSize, querySeqs, true, kmerThr);

    // compute match prob for local match
    double kmerMatchProb = ((double) stats.doubleMatches) / ((double) (stats.querySeqLen * targetSeqLenSum));
    kmerMatchProb /= 256;

    kmerMatchProb = std::max(kmerMatchProb, std::numeric_limits<double>::min());
    Debug(Debug::INFO) << "\tk-mers per position = " << stats.kmersPerPos << ", k-mer match probability: "
                       << kmerMatchProb << "\n";
    delete[] querySeqs;
    return kmerMatchProb;
}

statistics_t Prefiltering::computeStatisticForKmerThreshold(DBReader<unsigned int> *qdbr, Sequence** qseq,
                                                            IndexTable *indexTable, size_t querySetSize,
                                                            unsigned int *querySeqsIds, bool reverseQuery,
                                                            int kmerThrMid) {
    // determine k-mer match probability for kmerThrMid
    QueryMatcher **matchers = createQueryTemplateMatcher(subMat, indexTable, qseq, tdbr->getSeqLens(), kmerThrMid,
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
    for (size_t i = 0; i < querySetSize; i++) {
        size_t id = querySeqsIds[i];

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        char *seqData = qdbr->getData(id);
        unsigned int qKey = qdbr->getDbKey(id);
        qseq[thread_idx]->mapSequence(id, qKey, seqData);
        if (reverseQuery == true) {
            qseq[thread_idx]->reverse();
        }
        matchers[thread_idx]->matchQuery(qseq[thread_idx], UINT_MAX);
        kmersPerPos += matchers[thread_idx]->getStatistics()->kmersPerPos;
        dbMatchesSum += matchers[thread_idx]->getStatistics()->dbMatches;
        querySeqLenSum += matchers[thread_idx]->getStatistics()->querySeqLen;
        doubleMatches += matchers[thread_idx]->getStatistics()->doubleMatches;
        diagonalOverflow += matchers[thread_idx]->getStatistics()->diagonalOverflow;
        resultsPassedPref += matchers[thread_idx]->getStatistics()->resultsPassedPrefPerSeq;
    }
    // clean up memory
    for (unsigned int j = 0; j < threads; j++) {
        delete matchers[j];
    }
    delete[] matchers;

    return statistics_t(kmersPerPos / (double) querySetSize, dbMatchesSum, doubleMatches,
                        querySeqLenSum, diagonalOverflow, resultsPassedPref / querySetSize);
}

void Prefiltering::mergeFiles(const std::vector<std::pair<std::string, std::string>> &splitFiles,
                              const std::string &outDB, const std::string &outDBIndex) {
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        mergeOutput(outDB, outDBIndex, splitFiles);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        const char **datafilesNames = new const char *[splitFiles.size()];
        const char **indexFilesNames = new const char *[splitFiles.size()];
        for (size_t i = 0; i < splitFiles.size(); i++) {
            datafilesNames[i] = splitFiles[i].first.c_str();
            indexFilesNames[i] = splitFiles[i].second.c_str();
        }
        DBWriter::mergeResults(outDB.c_str(), outDBIndex.c_str(), datafilesNames, indexFilesNames, splitFiles.size());
        delete[] datafilesNames;
        delete[] indexFilesNames;
    } else {
        Debug(Debug::ERROR) << "Invalid split mode: " << splitMode << "\n";
        EXIT(EXIT_FAILURE);
    }
}

int Prefiltering::getKmerThreshold(const float sensitivity) {
    double kmerThrBest = kmerScore;
    if (kmerThrBest == INT_MAX) {
        if (kmerSize == 5) {
            kmerThrBest = 123.75 - (sensitivity * 8.75);
        } else if (kmerSize == 6) {
            kmerThrBest = 138.75 - (sensitivity * 8.75);
        } else if (kmerSize == 7) {
            kmerThrBest = 154.75 - (sensitivity * 9.75);
        } else {
            Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return static_cast<int>(kmerThrBest);
}

size_t Prefiltering::estimateMemoryConsumption(int split, size_t dbSize, size_t resSize,
                                               int alphabetSize, int kmerSize, unsigned int threads) {
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
            + (QueryMatcher::MAX_RES_LIST_LEN * sizeof(hit_t))
            + (dbSizeSplit * 2 * sizeof(CounterResult) * 2) // BINS * binSize, (binSize = dbSize * 2 / BINS)
            // 2 is a security factor the size can increase during run
    );
    // extended matrix
    size_t extendedMatrix = sizeof(std::pair<short, unsigned int>) * static_cast<size_t>(pow(pow(alphabetSize, 3), 2));
    extendedMatrix += sizeof(std::pair<short, unsigned int>) * pow(pow(alphabetSize, 2), 2);
    // some memory needed to keep the index, ....
    size_t background = dbSize * 22;
    return residueSize + indexTableSize + threadSize + background + extendedMatrix;
}

std::pair<int, int> Prefiltering::optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr,
                                                int alphabetSize, int externalKmerSize, unsigned int threads) {
    for (int split = 1; split < 100; split++) {
        for (int kmerSize = 6; kmerSize <= 7; kmerSize++) {
            if (kmerSize == externalKmerSize || externalKmerSize == 0) { // 0: set k-mer based on aa size in database
                size_t aaUpperBoundForKmerSize = IndexTable::getUpperBoundAACountForKmerSize(kmerSize);
                if ((tdbr->getAminoAcidDBSize() / split) < aaUpperBoundForKmerSize) {
                    size_t neededSize = estimateMemoryConsumption(split, tdbr->getSize(), tdbr->getAminoAcidDBSize(),
                                                                  alphabetSize,
                                                                  kmerSize, threads);
                    if (neededSize < 0.9 * totalMemoryInByte) {
                        return std::make_pair(kmerSize, split);
                    }
                }
            }
        }
    }

    return std::make_pair(-1, -1);
}
