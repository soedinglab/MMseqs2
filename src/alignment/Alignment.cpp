#include "Alignment.h"
#include "Util.h"
#include "Debug.h"

#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "PrefilteringIndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

Alignment::Alignment(const std::string &querySeqDB, const std::string &querySeqDBIndex,
                     const std::string &targetSeqDB, const std::string &targetSeqDBIndex,
                     const std::string &prefDB, const std::string &prefDBIndex,
                     const std::string &outDB, const std::string &outDBIndex,
                     const Parameters &par) :
        covThr(par.covThr), targetCovThr(par.targetCovThr), evalThr(par.evalThr), seqIdThr(par.seqIdThr),
        includeIdentity(par.includeIdentity), addBacktrace(par.addBacktrace), realign(par.realign),
        threads(static_cast<unsigned int>(par.threads)), outDB(outDB), outDBIndex(outDBIndex),
        maxSeqLen(par.maxSeqLen), querySeqType(par.querySeqType), targetSeqType(par.targetSeqType),
        compBiasCorrection(par.compBiasCorrection), qseqdbr(NULL), qSeqLookup(NULL),
        tseqdbr(NULL), tidxdbr(NULL), tSeqLookup(NULL), templateDBIsIndex(false), earlyExit(par.earlyExit) {
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif

    int alignmentMode = par.alignmentMode;
    if (addBacktrace == true) {
        alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }

    if (realign == true) {
        alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        if (addBacktrace == false) {
            Debug(Debug::ERROR) << "Realign can only be used in combination with -a.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    swMode = Matcher::SCORE_ONLY;
    switch (alignmentMode) {
        case Parameters::ALIGNMENT_MODE_FAST_AUTO:
            if((covThr > 0.0 || targetCovThr > 0.0) && seqIdThr == 0.0) {
                swMode = Matcher::SCORE_COV; // fast
            } else if((covThr > 0.0 || targetCovThr > 0.0) && seqIdThr > 0.0) { // if seq id is needed
                swMode = Matcher::SCORE_COV_SEQID; // slowest
            }
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV:
            swMode = Matcher::SCORE_COV; // fast
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
            swMode = Matcher::SCORE_COV_SEQID; // fast
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

    std::string indexDB = PrefilteringIndexReader::searchForIndex(targetSeqDB);
    //TODO optimize this. Dont read twice the target index. This seems to be slow
    if (indexDB != "") {
        Debug(Debug::INFO) << "Use index  " << indexDB << "\n";
        tseqdbr = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str());
    } else {
        tseqdbr = new DBReader<unsigned int>(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
    }
    tseqdbr->open(DBReader<unsigned int>::NOSORT);

    std::string scoringMatrixFile = par.scoringMatrixFile;
    bool reopenTargetDb = false;
    templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tseqdbr);
    if (templateDBIsIndex == true) { // exchange reader with old ffindex reader
        tidxdbr = tseqdbr;
        PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(tseqdbr);
        if (meta.split != 1) {
            Debug(Debug::WARNING) << "Can not use split index for alignment.\n";
            reopenTargetDb = true;

        } else {
            tSeqLookup = PrefilteringIndexReader::getSequenceLookup(tidxdbr, 0);
            tseqdbr = PrefilteringIndexReader::openNewReader(tidxdbr);
            scoringMatrixFile = PrefilteringIndexReader::getSubstitutionMatrixName(tidxdbr);
        }
    } else {
        reopenTargetDb = true;
    }

    if (reopenTargetDb) {
        tseqdbr->close();
        delete tseqdbr;
        tidxdbr = NULL;

        tseqdbr = new DBReader<unsigned int>(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
        tseqdbr->open(DBReader<unsigned int>::NOSORT);
    }

    if (par.noPreload == false) {
        tseqdbr->readMmapedDataInMemory();
        tseqdbr->mlock();
    }

    sameQTDB = (targetSeqDB.compare(querySeqDB) == 0);
    if (sameQTDB == true) {
        qseqdbr = tseqdbr;
        qSeqLookup = tSeqLookup;
    } else {
        // open the sequence, prefiltering and output databases
        qseqdbr = new DBReader<unsigned int>(querySeqDB.c_str(), querySeqDBIndex.c_str());
        qseqdbr->open(DBReader<unsigned int>::NOSORT);
        qseqdbr->readMmapedDataInMemory();
        qseqdbr->mlock();
    }

    prefdbr = new DBReader<unsigned int>(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);

    if (querySeqType == Sequence::NUCLEOTIDES) {
        m = new NucleotideMatrix();
    } else {
        // keep score bias at 0.0 (improved ROC)
        m = new SubstitutionMatrix(scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    if (realign == true) {
        realign_m = new SubstitutionMatrix(scoringMatrixFile.c_str(), 2.0, -0.2f);
    } else {
        realign_m = NULL;
    }
}

Alignment::~Alignment() {
    if (realign == true) {
        delete realign_m;
    }
    delete m;

    tseqdbr->close();
    delete tseqdbr;

    if (templateDBIsIndex == true) {
        tidxdbr->close();
        delete tidxdbr;

        delete tSeqLookup;
    }

    if (sameQTDB == false) {
        qseqdbr->close();
        delete qseqdbr;
    } else {
        qSeqLookup = NULL;
    }

    prefdbr->close();
    delete prefdbr;
}

void Alignment::run(const unsigned int mpiRank, const unsigned int mpiNumProc,
                    const unsigned int maxAlnNum, const unsigned int maxRejected) {

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(getQueryDbSize(), qseqdbr->getSeqLens(), getQueryDbEntries(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << (dbFrom + dbSize) << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDB, outDBIndex, mpiRank);
    run(tmpOutput.first, tmpOutput.second, dbFrom, dbSize, maxAlnNum, maxRejected);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (mpiRank == 0) {
        // master reduces results
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
    const size_t flushSize = 1000000;
    size_t iterations = static_cast<size_t>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence qSeq(maxSeqLen, m->aa2int, m->int2aa, querySeqType, 0, false, compBiasCorrection);
        Sequence dbSeq(maxSeqLen, m->aa2int, m->int2aa, targetSeqType, 0, false, compBiasCorrection);
        Matcher matcher(maxSeqLen, m, getTargetDbSize(), getTargetDbEntries(), compBiasCorrection);
        Matcher *realigner = NULL;
        if (realign) {
            realigner = new Matcher(maxSeqLen, realign_m, getTargetDbSize(), getTargetDbEntries(), compBiasCorrection);
        }

        for (size_t i = 0; i < iterations; i++) {
            size_t start = dbFrom + (i * flushSize);
            size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);

            #pragma omp for schedule(dynamic, 100) reduction(+: alignmentsNum, totalPassedNum)
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
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    // sequence are identical if qID == dbID  (needed to cluster really short sequences)
                    const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB)) ? true : false;

                    setTargetSequence(dbSeq, dbKey);
                    // check if the sequences could pass the coverage threshold
                    if ((((float) qSeq.L) / ((float) dbSeq.L) < covThr) ||
                        (((float) dbSeq.L) / ((float) qSeq.L) < covThr)) {
                        rejected++;
                        data = Util::skipLine(data);
                        continue;
                    }
                    // calculate Smith-Waterman alignment
                    Matcher::result_t res = matcher.getSWResult(&dbSeq, getTargetDbEntries(), evalThr, swMode);
                    alignmentsNum++;
                    //set coverage and seqid if identity
                    if (isIdentity) {
                        res.qcov = 1.0f;
                        res.dbcov = 1.0f;
                        res.seqId = 1.0f;
                    }

                    const bool evalOk = (res.eval <= evalThr); // -e
                    const bool seqIdOK = (res.seqId >= seqIdThr); // --min-seq-id
                    const bool covOK = (res.qcov >= covThr && res.dbcov >= covThr); //-c
                    const bool targetCovOK = (res.dbcov >= targetCovThr); // --target-cov (-c = 0.0 if --targetCov)
                    // check first if it is identity
                    if (isIdentity
                        ||
                        // general accaptance criteria
                        ( evalOk   &&
                          seqIdOK  &&
                          covOK    &&
                          targetCovOK
                        ))
                    {

                        swResults.push_back(res);
                        passedNum++;
                        totalPassedNum++;
                        rejected = 0;
                    } else {
                        rejected++;
                    }
                    data = Util::skipLine(data);
                }

                // write the results
                std::sort(swResults.begin(), swResults.end(), Matcher::compareHits);
                if (realign == true) {
                    realigner->initQuery(&qSeq);
                    for (size_t result = 0; result < swResults.size(); result++) {
                        setTargetSequence(dbSeq, swResults[result].dbKey);
                        Matcher::result_t res = realigner->getSWResult(&dbSeq, getTargetDbEntries(), 0.0,
                                                                       Matcher::SCORE_COV_SEQID);
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
                }

                // put the contents of the swResults list into ffindex DB
                std::stringstream swResultsString;
                for (size_t result = 0; result < swResults.size(); result++) {
                    swResultsString << Matcher::resultToString(swResults[result], addBacktrace);
                }

                std::string swResultString = swResultsString.str();
                const char *swResultsData = swResultString.c_str();
                dbw.writeData(swResultsData, swResultString.length(), qSeq.getDbKey(), thread_idx);
                swResults.clear();
            }
            prefdbr->remapData();
        }
#ifndef HAVE_MPI
        if (earlyExit) {
            #pragma omp barrier
            if(thread_idx == 0) {
                dbw.close();
                Debug(Debug::INFO) << "Done. Exiting early now.\n";
            }
            #pragma omp barrier

            _Exit(EXIT_SUCCESS);
        }
#endif

        if (realigner) {
            delete realigner;
        }
    }
    Debug(Debug::INFO) << "\n";
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

inline size_t Alignment::getQueryDbEntries() const {
    if (qSeqLookup != NULL) {
        return qSeqLookup->getSequenceCount();
    } else {
        return qseqdbr->getSize();
    }
}

inline size_t Alignment::getQueryDbSize() const {
    if (qSeqLookup != NULL) {
        return qSeqLookup->getDataSize();
    } else {
        return qseqdbr->getAminoAcidDBSize();
    }
}

inline void Alignment::setQuerySequence(Sequence &seq, size_t id, unsigned int key) {
    if (qSeqLookup != NULL) {
        std::pair<const unsigned char*, const unsigned int> sequence = qSeqLookup->getSequence(id);
        seq.mapSequence(id, key, sequence);
    } else {
        // map the query sequence
        char *querySeqData = qseqdbr->getDataByDBKey(key);
        if (querySeqData == NULL) {
#pragma omp critical
            {
                Debug(Debug::ERROR) << "ERROR: Query sequence " << key
                                    << " is required in the prefiltering, but is not contained in the query sequence database!\n"
                                    << "Please check your database.\n";
                EXIT(EXIT_FAILURE);
            }
        }

        seq.mapSequence(id, key, querySeqData);
    }
}

inline size_t Alignment::getTargetDbEntries() const {
    if (tSeqLookup != NULL) {
        return tSeqLookup->getSequenceCount();
    } else {

        return tseqdbr->getSize();
    }
}

inline size_t Alignment::getTargetDbSize() const {
    if (tSeqLookup != NULL) {
        return tSeqLookup->getDataSize();
    } else {
        return tseqdbr->getAminoAcidDBSize();
    }
}

inline void Alignment::setTargetSequence(Sequence &seq, unsigned int key) {
    if (tSeqLookup != NULL) {
        size_t id = tseqdbr->getId(key);
        std::pair<const unsigned char*, const unsigned int> sequence = tSeqLookup->getSequence(id);
        seq.mapSequence(id, key, sequence);
    } else {
        char *dbSeqData = tseqdbr->getDataByDBKey(key);
        if (dbSeqData == NULL) {
#pragma omp critical
            {
                Debug(Debug::ERROR) << "ERROR: Sequence " << key
                                    << " is required in the prefiltering, but is not contained in the target sequence database!\n"
                                    <<
                                    "Please check your database.\n";
                EXIT(EXIT_FAILURE);
            }
        }
        //char *maskedDbSeq = seg[thread_idx]->maskseq(dbSeqData);
        seq.mapSequence(static_cast<size_t>(-1), key, dbSeqData);
    }
}
