#include <SubstitutionMatrixProfileStates.h>
#include <QueryMatcher.h>
#include "Alignment.h"
#include "Util.h"
#include "Debug.h"

#include "Matcher.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

Alignment::Alignment(const std::string &querySeqDB, const std::string &querySeqDBIndex,
                     const std::string &targetSeqDB, const std::string &targetSeqDBIndex,
                     const std::string &prefDB, const std::string &prefDBIndex,
                     const std::string &outDB, const std::string &outDBIndex,
                     const Parameters &par) :

        covThr(par.covThr), covMode(par.covMode), seqIdMode(par.seqIdMode), evalThr(par.evalThr), seqIdThr(par.seqIdThr),
        includeIdentity(par.includeIdentity), addBacktrace(par.addBacktrace), realign(par.realign), scoreBias(par.scoreBias),
        threads(static_cast<unsigned int>(par.threads)), outDB(outDB), outDBIndex(outDBIndex),
        maxSeqLen(par.maxSeqLen), compBiasCorrection(par.compBiasCorrection), altAlignment(par.altAlignment), qdbr(NULL), qSeqLookup(NULL),
        tdbr(NULL), tidxdbr(NULL), tSeqLookup(NULL), templateDBIsIndex(false) {


    unsigned int alignmentMode = par.alignmentMode;
    if (alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED) {
        Debug(Debug::ERROR) << "Use rescorediagonal for ungapped alignment mode.\n";
        EXIT(EXIT_FAILURE);
    }

    if (addBacktrace == true) {
        alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }

    if (realign == true) {
        alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        if (addBacktrace == false) {
            Debug(Debug::WARNING) << "Turn on backtrace for realign.\n";
            addBacktrace = true;
        }
    }

    if (altAlignment > 0) {
        if(querySeqType==Sequence::NUCLEOTIDES){
            Debug(Debug::ERROR) << "Alternative alignments are not supported for nucleotides.\n";
            EXIT(EXIT_FAILURE);
        }
//        if(realign==true){
//            Debug(Debug::ERROR) << "Alternative alignments do not supported realignment.\n";
//            EXIT(EXIT_FAILURE);
//        }
        alignmentMode = (alignmentMode > Parameters::ALIGNMENT_MODE_SCORE_COV) ? alignmentMode : Parameters::ALIGNMENT_MODE_SCORE_COV;
    }

    initSWMode(alignmentMode);

    std::string scoringMatrixFile = par.scoringMatrixFile;
    std::string indexDB = PrefilteringIndexReader::searchForIndex(targetSeqDB);
    if (indexDB.length() > 0) {
        Debug(Debug::INFO) << "Use index  " << indexDB << "\n";

        tidxdbr = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str());
        tidxdbr->open(DBReader<unsigned int>::NOSORT);

        templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tidxdbr);
        if (templateDBIsIndex == true) {
            tSeqLookup = PrefilteringIndexReader::getUnmaskedSequenceLookup(tidxdbr, par.noPreload == false);
            if (tSeqLookup == NULL) {
                Debug(Debug::WARNING) << "No unmasked index available. Falling back to sequence database.\n";
                templateDBIsIndex = false;
            } else {
                PrefilteringIndexReader::printSummary(tidxdbr);
                PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(tidxdbr);
                targetSeqType = meta.seqType;
                tdbr = PrefilteringIndexReader::openNewReader(tidxdbr, par.noPreload == false);
                scoringMatrixFile = PrefilteringIndexReader::getSubstitutionMatrixName(tidxdbr);
            }
        }

        if (templateDBIsIndex == false) {
            tidxdbr->close();
            delete tidxdbr;
            tidxdbr = NULL;
        }
    }

    if (templateDBIsIndex == false) {
        tdbr = new DBReader<unsigned int>(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.noPreload == false) {
            tdbr->readMmapedDataInMemory();
            tdbr->mlock();
        }
    }

    sameQTDB = (targetSeqDB.compare(querySeqDB) == 0);
    if (sameQTDB == true) {
        qdbr = tdbr;
        qSeqLookup = tSeqLookup;
        querySeqType = targetSeqType;
    } else {
        // open the sequence, prefiltering and output databases
        qdbr = new DBReader<unsigned int>(querySeqDB.c_str(), querySeqDBIndex.c_str());
        qdbr->open(DBReader<unsigned int>::NOSORT);

        //size_t freeSpace =  FileUtil::getFreeSpace(FileUtil::dirName(outDB).c_str());
        //size_t estimatedHDDMemory = estimateHDDMemoryConsumption(qdbr->getSize(),
        //                                                         std::min(static_cast<size_t>(par.maxAccept), par.maxResListLen));
        //if (freeSpace < estimatedHDDMemory){
        //    Debug(Debug::ERROR) << "Hard disk has not enough space (" << freeSpace << " bytes left) "
        //                        << "to store " << estimatedHDDMemory << " bytes of results.\n"
        //                        << "Please free disk space and start MMseqs again.\n";
        //    EXIT(EXIT_FAILURE);
        //}

        qdbr->readMmapedDataInMemory();
        qdbr->mlock();
        querySeqType = qdbr->getDbtype();
    }

    if (qdbr->getSize() <= threads) {
        threads = qdbr->getSize();
    }

#ifdef OPENMP
    omp_set_num_threads(threads);
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif

    if (templateDBIsIndex == false) {
        querySeqType = qdbr->getDbtype();
        targetSeqType = tdbr->getDbtype();
    }
    if (querySeqType == -1 || targetSeqType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if (querySeqType == Sequence::HMM_PROFILE && targetSeqType == Sequence::HMM_PROFILE) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if (querySeqType != Sequence::HMM_PROFILE && targetSeqType == Sequence::PROFILE_STATE_SEQ) {
        Debug(Debug::ERROR) << "The query has to be a profile when using a target profile state database.\n";
        EXIT(EXIT_FAILURE);
    } else if (querySeqType == Sequence::HMM_PROFILE && targetSeqType == Sequence::PROFILE_STATE_SEQ) {
        querySeqType = Sequence::PROFILE_STATE_PROFILE;
    }
    Debug(Debug::INFO) << "Query database type: " << DBReader<unsigned int>::getDbTypeName(querySeqType) << "\n";
    Debug(Debug::INFO) << "Target database type: " << DBReader<unsigned int>::getDbTypeName(targetSeqType) << "\n";

    prefdbr = new DBReader<unsigned int>(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);

    if (querySeqType == Sequence::NUCLEOTIDES) {
        m = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, scoreBias);
        gapOpen = 7;
        gapExtend = 1;
    } else if (querySeqType == Sequence::PROFILE_STATE_PROFILE){
        SubstitutionMatrix s(par.scoringMatrixFile.c_str(), 2.0, scoreBias);
        this->m = new SubstitutionMatrixProfileStates(s.matrixName, s.probMatrix, s.pBack, s.subMatrixPseudoCounts, 2.0, scoreBias, 255);
        gapOpen = par.gapOpen;
        gapExtend = par.gapExtend;
    } else {
        // keep score bias at 0.0 (improved ROC)
        m = new SubstitutionMatrix(scoringMatrixFile.c_str(), 2.0, scoreBias);
        gapOpen = par.gapOpen;
        gapExtend = par.gapExtend;
    }

    if (realign == true) {
        realign_m = new SubstitutionMatrix(scoringMatrixFile.c_str(), 2.0, scoreBias-0.2f);
    } else {
        realign_m = NULL;
    }
}

void Alignment::initSWMode(unsigned int alignmentMode) {
    switch (alignmentMode) {
        case Parameters::ALIGNMENT_MODE_FAST_AUTO:
            if(covThr > 0.0 && seqIdThr == 0.0) {
                swMode = Matcher::SCORE_COV; // fast
            } else if(covThr > 0.0  && seqIdThr > 0.0) { // if seq id is needed
                swMode = Matcher::SCORE_COV_SEQID; // slowest
            } else {
                swMode = Matcher::SCORE_ONLY;
            }
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV:
            swMode = Matcher::SCORE_COV; // fast
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
            swMode = Matcher::SCORE_COV_SEQID; // slowest
            break;
        default:
            swMode = Matcher::SCORE_ONLY;
            break;
    }

    // print out mode and check for errors
    switch (swMode) {
        case Matcher::SCORE_ONLY:
            Debug(Debug::INFO) << "Compute score only.\n";
            break;
        case Matcher::SCORE_COV:
            Debug(Debug::INFO) << "Compute score and coverage.\n";
            break;
        case Matcher::SCORE_COV_SEQID:
            Debug(Debug::INFO) << "Compute score, coverage and sequence id.\n";
            break;
        default:
            Debug(Debug::ERROR) << "Wrong swMode mode.\n";
            EXIT(EXIT_FAILURE);
    }
}

Alignment::~Alignment() {
    if (realign == true) {
        delete realign_m;
    }
    delete m;

    tdbr->close();
    delete tdbr;

    if (templateDBIsIndex == true) {
        delete tSeqLookup;
        tidxdbr->close();
        delete tidxdbr;
    }

    if (sameQTDB == false) {
        qdbr->close();
        delete qdbr;
    }

    prefdbr->close();
    delete prefdbr;
}

void Alignment::run(const unsigned int mpiRank, const unsigned int mpiNumProc,
                    const unsigned int maxAlnNum, const unsigned int maxRejected) {

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(prefdbr->getAminoAcidDBSize(), prefdbr->getSeqLens(),
                                     prefdbr->getSize(), mpiRank, mpiNumProc, &dbFrom, &dbSize);

    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << (dbFrom + dbSize) << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDB, outDBIndex, mpiRank);
    run(tmpOutput.first, tmpOutput.second, dbFrom, dbSize, maxAlnNum, maxRejected);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for (unsigned int proc = 0; proc < mpiNumProc; proc++) {
            splitFiles.push_back(Util::createTmpFileNames(outDB, outDBIndex, proc));
        }

        // merge output databases
        DBWriter::mergeResults(outDB, outDBIndex, splitFiles);
    }
}

void Alignment::run(const unsigned int maxAlnNum, const unsigned int maxRejected) {
    run(outDB, outDBIndex, 0, prefdbr->getSize(), maxAlnNum, maxRejected);
}

void Alignment::run(const std::string &outDB, const std::string &outDBIndex,
                    const size_t dbFrom, const size_t dbSize,
                    const unsigned int maxAlnNum, const unsigned int maxRejected) {
    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;

    DBWriter dbw(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw.open();

    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), this->m, gapOpen, gapExtend, true);
    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 1000000;
    if(totalMemory > prefdbr->getDataSize()){
        flushSize = dbSize;
    }
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string alnResultsOutString;
        alnResultsOutString.reserve(1024*1024);
        char buffer[1024+32768];
        Sequence qSeq(maxSeqLen, querySeqType, m, 0, false, compBiasCorrection);
        Sequence dbSeq(maxSeqLen, targetSeqType, m, 0, false, compBiasCorrection);
        Matcher matcher(querySeqType, maxSeqLen, m, &evaluer, compBiasCorrection, gapOpen, gapExtend);
        Matcher *realigner = NULL;
        if (realign ==  true) {
            realigner = new Matcher(querySeqType, maxSeqLen, realign_m, &evaluer, compBiasCorrection, gapOpen, gapExtend);
        }

        size_t iterations = static_cast<size_t>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));
        for (size_t i = 0; i < iterations; i++) {
            size_t start = dbFrom + (i * flushSize);
            size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);

#pragma omp for schedule(dynamic, 5) reduction(+: alignmentsNum, totalPassedNum)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);

                // get the prefiltering list
                char *data = prefdbr->getData(id);
                unsigned int queryDbKey = prefdbr->getDbKey(id);
                setQuerySequence(qSeq, id, queryDbKey);

                matcher.initQuery(&qSeq);
                // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list
                std::vector<Matcher::result_t> swResults;
                size_t passedNum = 0;
                unsigned int rejected = 0;

                while (*data != '\0' && passedNum < maxAlnNum && rejected < maxRejected) {
                    // DB key of the db sequence
                    char dbKeyBuffer[255 + 1];
                    char * words[10];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);

                    size_t elements = Util::getWordsOfLine(data, words, 10);
                    int diagonal = INT_MAX;
                    // Prefilter result (need to make this better)
                    if(elements == 3){
                        hit_t hit = QueryMatcher::parsePrefilterHit(data);
                        diagonal = hit.diagonal;
                    }

                    setTargetSequence(dbSeq, dbKey);
                    // check if the sequences could pass the coverage threshold
                    if(Util::canBeCovered(covThr, covMode, static_cast<float>(qSeq.L), static_cast<float>(dbSeq.L)) == false )
                    {
                        rejected++;
                        data = Util::skipLine(data);
                        continue;
                    }
                    const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB)) ? true : false;

                    // calculate Smith-Waterman alignment
                    Matcher::result_t res = matcher.getSWResult(&dbSeq, diagonal, covMode, covThr, evalThr, swMode, seqIdMode, isIdentity);
                    alignmentsNum++;

                    //set coverage and seqid if identity
                    if (isIdentity) {
                        res.qcov = 1.0f;
                        res.dbcov = 1.0f;
                        res.seqId = 1.0f;
                    }
                    if(checkCriteriaAndAddHitToList(res, isIdentity, swResults)){
                        passedNum++;
                        totalPassedNum++;
                        rejected = 0;
                    }else{
                        rejected++;
                    }

                    data = Util::skipLine(data);
                }
                if(altAlignment > 0 && realign == false ){
                    computeAlternativeAlignment(queryDbKey, dbSeq, swResults, matcher, evalThr, swMode);
                }

                // write the results
                std::sort(swResults.begin(), swResults.end(), Matcher::compareHits);
                if (realign == true) {
                    realigner->initQuery(&qSeq);
                    for (size_t result = 0; result < swResults.size(); result++) {
                        setTargetSequence(dbSeq, swResults[result].dbKey);
                        const bool isIdentity = (queryDbKey == swResults[result].dbKey && (includeIdentity || sameQTDB)) ? true : false;
                        Matcher::result_t res = realigner->getSWResult(&dbSeq, INT_MAX, covMode, covThr, FLT_MAX,
                                                                       Matcher::SCORE_COV_SEQID, seqIdMode, isIdentity);
                        swResults[result].backtrace  = res.backtrace;
                        swResults[result].qStartPos  = res.qStartPos;
                        swResults[result].qEndPos    = res.qEndPos;
                        swResults[result].dbStartPos = res.dbStartPos;
                        swResults[result].dbEndPos   = res.dbEndPos;
                        swResults[result].alnLength  = res.alnLength;
                        swResults[result].seqId      = res.seqId;
                        swResults[result].qcov       = res.qcov;
                        swResults[result].dbcov      = res.dbcov;
                    }
                    if(altAlignment> 0 ){
                        computeAlternativeAlignment(queryDbKey, dbSeq, swResults, matcher, FLT_MAX, Matcher::SCORE_COV_SEQID);
                    }
                }

                // put the contents of the swResults list into ffindex DB
                for (size_t result = 0; result < swResults.size(); result++) {
                    size_t len = Matcher::resultToBuffer(buffer, swResults[result], addBacktrace);
                    alnResultsOutString.append(buffer, len);
                }
                dbw.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), qSeq.getDbKey(), thread_idx);
                alnResultsOutString.clear();
            }

#pragma omp barrier
            if (thread_idx == 0) {
                prefdbr->remapData();
            }
#pragma omp barrier
        }

        if (realign == true) {
            delete realigner;
        }
    }

    dbw.close();

    Debug(Debug::INFO) << "\nAll sequences processed.\n\n";
    Debug(Debug::INFO) << alignmentsNum << " alignments calculated.\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds ("
                       << ((float) totalPassedNum / (float) alignmentsNum) << " of overall calculated).\n";

    size_t hits = totalPassedNum / dbSize;
    size_t hits_rest = totalPassedNum % dbSize;
    float hits_f = ((float) hits) + ((float) hits_rest) / (float) dbSize;
    Debug(Debug::INFO) << hits_f << " hits per query sequence.\n";
}

inline void Alignment::setQuerySequence(Sequence &seq, size_t id, unsigned int key) {
    if (qSeqLookup != NULL) {
        std::pair<const unsigned char*, const unsigned int> sequence = qSeqLookup->getSequence(id);
        seq.mapSequence(id, key, sequence);
    } else {
        // map the query sequence
        char *querySeqData = qdbr->getDataByDBKey(key);
        if (querySeqData == NULL) {
#pragma omp critical
            {
                Debug(Debug::ERROR) << "ERROR: Query sequence " << key
                                    << " is required in the prefiltering, "
                                    << "but is not contained in the query sequence database!\n"
                                    << "Please check your database.\n";
                EXIT(EXIT_FAILURE);
            }
        }

        seq.mapSequence(id, key, querySeqData);
    }
}

inline void Alignment::setTargetSequence(Sequence &seq, unsigned int key) {
    if (tSeqLookup != NULL) {
        size_t id = tdbr->getId(key);
        std::pair<const unsigned char*, const unsigned int> sequence = tSeqLookup->getSequence(id);
        seq.mapSequence(id, key, sequence);
    } else {
        char *dbSeqData = tdbr->getDataByDBKey(key);
        if (dbSeqData == NULL) {
#pragma omp critical
            {
                Debug(Debug::ERROR) << "ERROR: Sequence " << key
                                    << " is required in the prefiltering,"
                                    << "but is not contained in the target sequence database!\n"
                                    << "Please check your database.\n";
                EXIT(EXIT_FAILURE);
            }
        }
        seq.mapSequence(static_cast<size_t>(-1), key, dbSeqData);
    }
}


size_t Alignment::estimateHDDMemoryConsumption(int dbSize, int maxSeqs) {
    return 2 * (dbSize * maxSeqs * 21 * 1.75);
}


bool Alignment::checkCriteriaAndAddHitToList(Matcher::result_t &res, bool isIdentity, std::vector<Matcher::result_t> &swHits){
    const bool evalOk = (res.eval <= evalThr); // -e
    const bool seqIdOK = (res.seqId >= seqIdThr); // --min-seq-id
    const bool covOK = Util::hasCoverage(covThr, covMode, res.qcov, res.dbcov);
    // check first if it is identity
    if (isIdentity
        ||
        // general accaptance criteria
        ( evalOk   &&
          seqIdOK  &&
          covOK
        ))
    {
        swHits.push_back(res);
        return true;
    } else {
        return false;
    }
}

void Alignment::computeAlternativeAlignment(unsigned int queryDbKey, Sequence &dbSeq,
                                            std::vector<Matcher::result_t> &swResults,
                                            Matcher &matcher, float evalThr, int swMode) {
    int xIndex = m->aa2int[static_cast<int>('X')];
    size_t firstItResSize = swResults.size();
    for(size_t i = 0; i < firstItResSize; i++) {
        const bool isIdentity = (queryDbKey == swResults[i].dbKey && (includeIdentity || sameQTDB))
                                ? true : false;
        if (isIdentity == true) {
            continue;
        }
        setTargetSequence(dbSeq, swResults[i].dbKey);
        for (int pos = swResults[i].dbStartPos; pos < swResults[i].dbEndPos; ++pos) {
            dbSeq.int_sequence[pos] = xIndex;
        }
        bool nextAlignment = true;
        for (int altAli = 0; altAli < altAlignment && nextAlignment; altAli++) {
            Matcher::result_t res = matcher.getSWResult(&dbSeq, INT_MAX, covMode, covThr, evalThr, swMode,
                                                        seqIdMode, isIdentity);
            nextAlignment = checkCriteriaAndAddHitToList(res, isIdentity, swResults);
            if (nextAlignment == true) {
                for (int pos = res.dbStartPos; pos < res.dbEndPos; pos++) {
                    dbSeq.int_sequence[pos] = xIndex;
                }
            }
        }
    }
}
