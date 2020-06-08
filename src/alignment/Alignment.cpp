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
#include "LinsearchIndexReader.h"
#include "IndexReader.h"
#include "Parameters.h"


#ifdef OPENMP
#include <omp.h>
#endif

Alignment::Alignment(const std::string &querySeqDB,
                     const std::string &targetSeqDB,
                     const std::string &prefDB, const std::string &prefDBIndex,
                     const std::string &outDB, const std::string &outDBIndex,
                     const Parameters &par) :

        covThr(par.covThr), canCovThr(par.covThr), covMode(par.covMode), seqIdMode(par.seqIdMode), evalThr(par.evalThr), seqIdThr(par.seqIdThr),
        alnLenThr(par.alnLenThr), includeIdentity(par.includeIdentity), addBacktrace(par.addBacktrace), realign(par.realign), scoreBias(par.scoreBias),
        threads(static_cast<unsigned int>(par.threads)), compressed(par.compressed), outDB(outDB), outDBIndex(outDBIndex),
        maxSeqLen(par.maxSeqLen), compBiasCorrection(par.compBiasCorrection), altAlignment(par.altAlignment), qdbr(NULL), qDbrIdx(NULL),
        tdbr(NULL), tDbrIdx(NULL) {


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
        realignCov = par.covThr;
        covThr = 0.0;
        if (addBacktrace == false) {
            Debug(Debug::WARNING) << "Turn on backtrace for realign.\n";
            addBacktrace = true;
        }
    }

    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    tDbrIdx = new IndexReader(targetSeqDB, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    tdbr = tDbrIdx->sequenceReader;
    targetSeqType = tdbr->getDbtype();
    sameQTDB = (targetSeqDB.compare(querySeqDB) == 0);
    if (sameQTDB == true) {
        qDbrIdx = tDbrIdx;
        qdbr = tdbr;
        querySeqType = targetSeqType;
    } else {
        // open the sequence, prefiltering and output databases
        qDbrIdx = new IndexReader(par.db1, par.threads,  IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0 );
        qdbr = qDbrIdx->sequenceReader;
        querySeqType = qdbr->getDbtype();
    }

    if (altAlignment > 0) {
        if(Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)){
            Debug(Debug::ERROR) << "Alternative alignments are not supported for nucleotides.\n";
            EXIT(EXIT_FAILURE);
        }
//        if(realign==true){
//            Debug(Debug::ERROR) << "Alternative alignments do not supported realignment.\n";
//            EXIT(EXIT_FAILURE);
//        }
        alignmentMode = (alignmentMode > Parameters::ALIGNMENT_MODE_SCORE_COV) ? alignmentMode : Parameters::ALIGNMENT_MODE_SCORE_COV;
    }
    swMode = initSWMode(alignmentMode, par.covThr, par.seqIdThr);

    if (par.wrappedScoring)
    {
        maxSeqLen = maxSeqLen * 2;
        if(!Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)){
            Debug(Debug::ERROR) << "Wrapped scoring is only supported for nucleotides.\n";
            EXIT(EXIT_FAILURE);
        }

        if (realign == true) {
            Debug(Debug::ERROR) << "Alternative alignments do not supported wrapped scoring.\n";
            EXIT(EXIT_FAILURE);
        }
    }


    //qdbr->readMmapedDataInMemory();
    // make sure to touch target after query, so if there is not enough memory for the query, at least the targets
    // might have had enough space left to be residung in the page cache
    if (sameQTDB == false && tDbrIdx == NULL && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tdbr->readMmapedDataInMemory();
    }

    if (qdbr->getSize() <= threads) {
        threads = qdbr->getSize();
    }

    if (querySeqType == -1 || targetSeqType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) == false && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        Debug(Debug::ERROR) << "The query has to be a profile when using a target profile state database.\n";
        EXIT(EXIT_FAILURE);
    } else if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        querySeqType = Parameters::DBTYPE_PROFILE_STATE_PROFILE;
    }
    Debug(Debug::INFO) << "Query database size: "  << qdbr->getSize() << " type: " << Parameters::getDbTypeName(querySeqType) << "\n";
    Debug(Debug::INFO) << "Target database size: " << tdbr->getSize() << " type: " << Parameters::getDbTypeName(targetSeqType) << "\n";

    prefdbr = new DBReader<unsigned int>(prefDB.c_str(), prefDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    prefdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    reversePrefilterResult = (Parameters::isEqualDbtype(prefdbr->getDbtype(), Parameters::DBTYPE_PREFILTER_REV_RES));

    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        m = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, scoreBias);
        gapOpen = par.gapOpen.nucleotides;
        gapExtend = par.gapExtend.nucleotides;
        zdrop = par.zdrop;
    } else if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)){
        SubstitutionMatrix s(par.scoringMatrixFile.aminoacids, 2.0, scoreBias);
        this->m = new SubstitutionMatrixProfileStates(s.matrixName, s.probMatrix, s.pBack, s.subMatrixPseudoCounts, 2.0, scoreBias, 219);
        gapOpen = par.gapOpen.aminoacids;
        gapExtend = par.gapExtend.aminoacids;
    } else {
        // keep score bias at 0.0 (improved ROC)
        m = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, scoreBias);
        gapOpen = par.gapOpen.aminoacids;
        gapExtend = par.gapExtend.aminoacids;
    }

    if (realign == true) {
        if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
            realign_m = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, scoreBias-0.2f);
        } else {
            realign_m = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, scoreBias-0.2f);
        }
    } else {
        realign_m = NULL;
    }
}

unsigned int Alignment::initSWMode(unsigned int alignmentMode, float covThr, float seqIdThr) {
    unsigned int swMode = Matcher::SCORE_ONLY;
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
            Debug(Debug::INFO) << "Compute score only\n";
            break;
        case Matcher::SCORE_COV:
            Debug(Debug::INFO) << "Compute score and coverage\n";
            break;
        case Matcher::SCORE_COV_SEQID:
            Debug(Debug::INFO) << "Compute score, coverage and sequence identity\n";
            break;
        default:
            Debug(Debug::ERROR) << "Wrong swMode mode\n";
            EXIT(EXIT_FAILURE);
    }
    return swMode;
}

Alignment::~Alignment() {
    if (realign == true) {
        delete realign_m;
    }
    delete m;

    if (tDbrIdx != NULL) {
        delete tDbrIdx;
    }else{
        tdbr->close();
        delete tdbr;
    }

    if (sameQTDB == false) {
        if(qDbrIdx != NULL){
            delete qDbrIdx;
        }else{
            qdbr->close();
            delete qdbr;
        }
    }

    prefdbr->close();
    delete prefdbr;
}

void Alignment::run(const unsigned int mpiRank, const unsigned int mpiNumProc,
                    const unsigned int maxAlnNum, const unsigned int maxRejected, bool wrappedScoring) {

    size_t dbFrom = 0;
    size_t dbSize = 0;
    prefdbr->decomposeDomainByAminoAcid( mpiRank, mpiNumProc, &dbFrom, &dbSize);

    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << (dbFrom + dbSize) << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDB, outDBIndex, mpiRank);
    run(tmpOutput.first, tmpOutput.second, dbFrom, dbSize, maxAlnNum, maxRejected, true, wrappedScoring);

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

void Alignment::run(const unsigned int maxAlnNum, const unsigned int maxRejected, bool wrappedScoring) {
    run(outDB, outDBIndex, 0, prefdbr->getSize(), maxAlnNum, maxRejected, false, wrappedScoring);
}

void Alignment::run(const std::string &outDB, const std::string &outDBIndex,
                    const size_t dbFrom, const size_t dbSize,
                    const unsigned int maxAlnNum, const unsigned int maxRejected, bool merge, bool wrappedScoring) {
    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;
    DBWriter dbw(outDB.c_str(), outDBIndex.c_str(), threads, compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    // handle no alignment case early, below would divide by 0 otherwise
    if (dbSize == 0) {
        dbw.close(merge);
        return;
    }

    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), this->m, gapOpen, gapExtend);
    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 1000000;
    if(totalMemory > prefdbr->getTotalDataSize()){
        flushSize = dbSize;
    }

    size_t iterations = static_cast<size_t>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = dbFrom + (i * flushSize);
        size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);

#pragma omp parallel num_threads(threads)
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
            Matcher matcher(querySeqType,
                                (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ? maxSeqLen : std::max(tdbr->getMaxSeqLen(), qdbr->getMaxSeqLen()),
                                 m, &evaluer, compBiasCorrection, gapOpen, gapExtend, zdrop);
            Matcher *realigner = NULL;
            if (realign ==  true && wrappedScoring == false) {
                realigner = new Matcher(querySeqType,
                                       (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ? maxSeqLen : std::max(tdbr->getMaxSeqLen(), qdbr->getMaxSeqLen()),
                                       realign_m, &evaluer, compBiasCorrection, gapOpen, gapExtend, zdrop);
            }

            std::vector<Matcher::result_t> swResults;
            swResults.reserve(300);
            std::vector<Matcher::result_t> swRealignResults;
            swRealignResults.reserve(300);
            std::vector<hit_t> shortResults;
            shortResults.reserve(300);

#pragma omp for schedule(dynamic, 5) reduction(+: alignmentsNum, totalPassedNum)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();

                // get the prefiltering list
                char *data = prefdbr->getData(id, thread_idx);
                unsigned int queryDbKey = prefdbr->getDbKey(id);
                size_t origQueryLen = 0;
                std::string queryToWrap;
                // only load query data if data != \0
                if(*data != '\0'){
                    size_t qId = qdbr->getId(queryDbKey);
                    char *querySeqData = qdbr->getData(qId, thread_idx);
                    if (querySeqData == NULL) {
                        Debug(Debug::ERROR) << "Query sequence " << queryDbKey
                                            << " is required in the prefiltering, but is not contained in the query sequence database.\nPlease check your database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    size_t queryLen = qdbr->getSeqLen(qId);
                    origQueryLen = queryLen;
                    if (wrappedScoring) {
                        queryToWrap = std::string(querySeqData,queryLen);
                        queryToWrap = queryToWrap + queryToWrap;
                        querySeqData = (char*)(queryToWrap).c_str();
                        queryLen = origQueryLen*2;
                    }

                    qSeq.mapSequence(qId, queryDbKey, querySeqData, queryLen);
                    matcher.initQuery(&qSeq);
                }

                // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list
                size_t passedNum = 0;
                unsigned int rejected = 0;
                while (*data != '\0' && passedNum < maxAlnNum && rejected < maxRejected) {
                    // DB key of the db sequence
                    char dbKeyBuffer[255 + 1];
                    const char* words[10];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);

                    size_t elements = Util::getWordsOfLine(data, words, 10);
                    short diagonal = 0;
                    bool isReverse = false;
                    // Prefilter result (need to make this better)
                    if(elements == 3){
                        hit_t hit = QueryMatcher::parsePrefilterHit(data);
                        isReverse = reversePrefilterResult && (hit.prefScore < 0);
                        diagonal = static_cast<short>(hit.diagonal);
                    }
                    size_t dbId = tdbr->getId(dbKey);
                    char *dbSeqData = tdbr->getData(dbId, thread_idx);

                    if (dbSeqData == NULL) {
                        Debug(Debug::ERROR) << "Sequence " << dbKey <<" is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    dbSeq.mapSequence(dbId, dbKey, dbSeqData, tdbr->getSeqLen(dbId));
                    // check if the sequences could pass the coverage threshold
                    if(Util::canBeCovered(canCovThr, covMode, static_cast<float>(origQueryLen), static_cast<float>(dbSeq.L)) == false) {
                        rejected++;
                        data = Util::skipLine(data);
                        continue;
                    }
                    const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB)) ? true : false;

                    // calculate Smith-Waterman alignment
                    Matcher::result_t res = matcher.getSWResult(&dbSeq, static_cast<int>(diagonal), isReverse, covMode, covThr, evalThr, swMode, seqIdMode, isIdentity, wrappedScoring);
                    alignmentsNum++;

                    //set coverage and seqid if identity
                    if (isIdentity) {
                        res.qcov = 1.0f;
                        res.dbcov = 1.0f;
                        res.seqId = 1.0f;
                    }
                    if(checkCriteria(res, isIdentity, evalThr, seqIdThr, alnLenThr, covMode, covThr)){

                        swResults.emplace_back(res);
                        passedNum++;
                        totalPassedNum++;
                        rejected = 0;
                    }else{
                        rejected++;
                    }

                    data = Util::skipLine(data);
                }
                if(altAlignment > 0 && realign == false && wrappedScoring == false){
                    computeAlternativeAlignment(queryDbKey, dbSeq, swResults, matcher, evalThr, swMode, thread_idx);
                }

                if(wrappedScoring && shortResults.size() > 1)
                    std::sort(shortResults.begin(), shortResults.end(), hit_t::compareHitsByScoreAndId);

                // write the results
                if(swResults.size() > 1)
                    std::sort(swResults.begin(), swResults.end(), Matcher::compareHits);
                if (realign == true) {
                    realigner->initQuery(&qSeq);
                    for (size_t result = 0; result < swResults.size(); result++) {
                        size_t dbId = tdbr->getId(swResults[result].dbKey);
                        char *dbSeqData = tdbr->getData(dbId, thread_idx);
                        if (dbSeqData == NULL) {
                            Debug(Debug::ERROR) << "Sequence " << swResults[result].dbKey <<" is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
                            EXIT(EXIT_FAILURE);
                        }
                        dbSeq.mapSequence(static_cast<size_t>(-1), swResults[result].dbKey, dbSeqData,
                                          tdbr->getSeqLen(dbId));
                        const bool isIdentity = (queryDbKey == swResults[result].dbKey && (includeIdentity || sameQTDB)) ? true : false;
                        Matcher::result_t res = realigner->getSWResult(&dbSeq, INT_MAX, false, covMode, covThr, FLT_MAX,
                                                                       Matcher::SCORE_COV_SEQID, seqIdMode, isIdentity);
                        const bool covOK = Util::hasCoverage(realignCov, covMode, res.qcov, res.dbcov);
                        if(covOK == true|| isIdentity){
                            swResults[result].backtrace  = res.backtrace;
                            swResults[result].qStartPos  = res.qStartPos;
                            swResults[result].qEndPos    = res.qEndPos;
                            swResults[result].dbStartPos = res.dbStartPos;
                            swResults[result].dbEndPos   = res.dbEndPos;
                            swResults[result].alnLength  = res.alnLength;
                            swResults[result].seqId      = res.seqId;
                            swResults[result].qcov       = res.qcov;
                            swResults[result].dbcov      = res.dbcov;
                            swRealignResults.push_back(swResults[result]);
                        }
                    }
                    swResults = swRealignResults;
                    if(altAlignment > 0){
                        computeAlternativeAlignment(queryDbKey, dbSeq, swResults, matcher, FLT_MAX, Matcher::SCORE_COV_SEQID, thread_idx);
                    }
                }

                // put the contents of the swResults list into a result DB
                for (size_t result = 0; result < swResults.size(); result++) {
                    size_t len = Matcher::resultToBuffer(buffer, swResults[result], addBacktrace);
                    alnResultsOutString.append(buffer, len);
                }

                for (size_t result = 0; result < shortResults.size(); result++) {
                    size_t len = snprintf(buffer, 100, "%u\t%d\t%d\n", shortResults[result].seqId, shortResults[result].prefScore,
                                          shortResults[result].diagonal);
                    alnResultsOutString.append(buffer, len);
                }

                dbw.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), queryDbKey, thread_idx);
                alnResultsOutString.clear();
                swResults.clear();
                swRealignResults.clear();
                shortResults.clear();
            }
            if (realign == true) {
                delete realigner;
            }
#pragma omp barrier
            if (thread_idx == 0) {
                prefdbr->remapData();
            }
#pragma omp barrier
        }


    }

    dbw.close(merge);

    Debug(Debug::INFO) << "\n" << alignmentsNum << " alignments calculated.\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds ("
                       << ((float) totalPassedNum / (float) alignmentsNum) << " of overall calculated).\n";

    size_t hits = totalPassedNum / dbSize;
    size_t hits_rest = totalPassedNum % dbSize;
    float hits_f = ((float) hits) + ((float) hits_rest) / (float) dbSize;
    Debug(Debug::INFO) << hits_f << " hits per query sequence.\n";
}

size_t Alignment::estimateHDDMemoryConsumption(int dbSize, int maxSeqs) {
    return 2 * (dbSize * maxSeqs * 21 * 1.75);
}


bool Alignment::checkCriteria(Matcher::result_t &res, bool isIdentity, double evalThr, double seqIdThr, int alnLenThr, int covMode, float covThr) {
    const bool evalOk = (res.eval <= evalThr); // -e
    const bool seqIdOK = (res.seqId >= seqIdThr); // --min-seq-id
    const bool covOK = Util::hasCoverage(covThr, covMode, res.qcov, res.dbcov);
    const bool alnLenOK = Util::hasAlignmentLength(alnLenThr, res.alnLength);

    // check first if it is identity
    if (isIdentity
        ||
        // general accaptance criteria
        ( evalOk   &&
          seqIdOK  &&
          covOK    &&
          alnLenOK
        ))
    {
        return true;
    } else {
        return false;
    }
}

void Alignment::computeAlternativeAlignment(unsigned int queryDbKey, Sequence &dbSeq,
                                            std::vector<Matcher::result_t> &swResults,
                                            Matcher &matcher, float evalThr, int swMode, int thread_idx) {
    unsigned char xIndex = m->aa2num[static_cast<int>('X')];
    size_t firstItResSize = swResults.size();
    for(size_t i = 0; i < firstItResSize; i++) {
        const bool isIdentity = (queryDbKey == swResults[i].dbKey && (includeIdentity || sameQTDB))
                                ? true : false;
        if (isIdentity == true) {
            continue;
        }
        size_t dbId = tdbr->getId(swResults[i].dbKey);
        char *dbSeqData = tdbr->getData(dbId, thread_idx);
        if (dbSeqData == NULL) {
            Debug(Debug::ERROR) << "Sequence " << swResults[i].dbKey <<" is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
            EXIT(EXIT_FAILURE);
        }

        dbSeq.mapSequence(dbId, swResults[i].dbKey, dbSeqData, tdbr->getSeqLen(dbId));
        for (int pos = swResults[i].dbStartPos; pos < swResults[i].dbEndPos; ++pos) {
            dbSeq.numSequence[pos] = xIndex;
        }
        bool nextAlignment = true;
        for (int altAli = 0; altAli < altAlignment && nextAlignment; altAli++) {
            Matcher::result_t res = matcher.getSWResult(&dbSeq, INT_MAX, false, covMode, covThr, evalThr, swMode,
                                                        seqIdMode, isIdentity);
            nextAlignment = checkCriteria(res, isIdentity, evalThr, seqIdThr, alnLenThr, covMode, covThr);
            if (nextAlignment == true) {
                swResults.emplace_back(res);
                for (int pos = res.dbStartPos; pos < res.dbEndPos; pos++) {
                    dbSeq.numSequence[pos] = xIndex;
                }
            }
        }
    }
}
