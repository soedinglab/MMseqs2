#include "Alignment.h"
#include "Util.h"
#include "Debug.h"

#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"

#ifdef OPENMP
#include <omp.h>
#endif

Alignment::Alignment(const std::string &querySeqDB, const std::string &querySeqDBIndex,
                     const std::string &targetSeqDB, const std::string &targetSeqDBIndex,
                     const std::string &prefDB, const std::string &prefDBIndex,
                     const std::string &outDB, const std::string &outDBIndex,
                     const Parameters &par) :
        covThr(par.covThr), evalThr(par.evalThr), seqIdThr(par.seqIdThr),
        includeIdentity(par.includeIdentity), addBacktrace(par.addBacktrace), realign(par.realign),
        threads(static_cast<unsigned int>(par.threads)), outDB(outDB), outDBIndex(outDBIndex) {
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

    switch (alignmentMode) {
        case Parameters::ALIGNMENT_MODE_FAST_AUTO:
            if (covThr > 0.0 && seqIdThr == 0.0) {
                swMode = Matcher::SCORE_COV; // fast
            } else if (covThr > 0.0 && seqIdThr > 0.0) { // if seq id is needed
                swMode = Matcher::SCORE_COV_SEQID; // slowest
            }
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV:
            swMode = Matcher::SCORE_COV; // fast
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
            swMode = Matcher::SCORE_COV_SEQID; // fast
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

    if (par.querySeqType == Sequence::NUCLEOTIDES) {
        m = new NucleotideMatrix();
    } else {
        // keep score bias to 0.0 (improved ROC)
        m = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    qSeqs = new Sequence *[threads];
    dbSeqs = new Sequence *[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++) {
        qSeqs[i] = new Sequence(par.maxSeqLen, m->aa2int, m->int2aa,
                                par.querySeqType, 0, false, par.compBiasCorrection);
        dbSeqs[i] = new Sequence(par.maxSeqLen, m->aa2int, m->int2aa,
                                 par.targetSeqType, 0, false, par.compBiasCorrection);
    }

    // open the sequence, prefiltering and output databases
    qseqdbr = new DBReader<unsigned int>(querySeqDB.c_str(), querySeqDBIndex.c_str());
    qseqdbr->open(DBReader<unsigned int>::NOSORT);
    if (par.noPreload != false) {
        qseqdbr->readMmapedDataInMemory();
    }

    sameQTDB = (querySeqDB.compare(targetSeqDB) == 0);
    if (sameQTDB == true) {
        tseqdbr = qseqdbr;
    } else {
        tseqdbr = new DBReader<unsigned int>(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
        tseqdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.noPreload != false) {
            tseqdbr->readMmapedDataInMemory();
        }
    }

    prefdbr = new DBReader<unsigned int>(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);

    matchers = new Matcher *[threads];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++) {
        matchers[i] = new Matcher(par.maxSeqLen, m, tseqdbr->getAminoAcidDBSize(),
                                  tseqdbr->getSize(), par.compBiasCorrection);
    }

    if (realign == true) {
        realign_m = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, -0.2);
        realigner = new Matcher *[threads];
#pragma omp parallel for schedule(static)
        for (int i = 0; i < threads; i++) {
            realigner[i] = new Matcher(par.maxSeqLen, realign_m, tseqdbr->getAminoAcidDBSize(),
                                       tseqdbr->getSize(), par.compBiasCorrection);
        }
    } else {
        realign_m = NULL;
        realigner = NULL;
    }

    dbKeys = new unsigned int[threads];
}

Alignment::~Alignment() {
    if (realign == true) {
        for (int i = 0; i < threads; i++) {
            delete realigner[i];
        }
        delete[] realigner;
        delete realign_m;
    }

    for (int i = 0; i < threads; i++) {
        delete qSeqs[i];
        delete dbSeqs[i];
        delete matchers[i];
    }

    delete[] qSeqs;
    delete[] dbSeqs;
    delete[] matchers;
    delete[] dbKeys;
    delete m;
    delete qseqdbr;

    if (sameQTDB == false) {
        delete tseqdbr;
    }
    delete prefdbr;
}

void Alignment::run(const unsigned int mpiRank, const unsigned int mpiNumProc,
                    const unsigned int maxAlnNum, const unsigned int maxRejected) {

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(qseqdbr->getAminoAcidDBSize(), qseqdbr->getSeqLens(), qseqdbr->getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << (dbFrom + dbSize) << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDB, outDBIndex, mpiRank);
    run(tmpOutput.first, tmpOutput.second, dbFrom, dbSize, maxAlnNum, maxRejected);

    // close reader to reduce memory
    closeReaders();

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
        mergeAndRemoveTmpDatabases(outDB, outDBIndex, splitFiles);
    }
}


void Alignment::closeReaders() {
    qseqdbr->close();
    if (sameQTDB == false) {
        tseqdbr->close();
    }
    prefdbr->close();
}

void Alignment::run(const unsigned int maxAlnNum, const unsigned int maxRejected) {
    run(outDB.c_str(), outDBIndex.c_str(), 0, prefdbr->getSize(), maxAlnNum, maxRejected);
    closeReaders();
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
    for (size_t i = 0; i < iterations; i++) {
        size_t start = dbFrom + (i * flushSize);
        size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);
# pragma omp parallel for schedule(dynamic, 100) reduction (+: alignmentsNum, totalPassedNum)
        for (size_t id = start; id < (start + bucketSize); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            // get the prefiltering list
            char *data = prefdbr->getData(id);
            unsigned int queryDbKey = prefdbr->getDbKey(id);
            // map the query sequence
            char *querySeqData = qseqdbr->getDataByDBKey(queryDbKey);
            if (querySeqData == NULL) {
                #pragma omp critical
                {
                    Debug(Debug::ERROR) << "ERROR: Query sequence " << queryDbKey
                                        << " is required in the prefiltering, but is not contained in the query sequence database!\n"
                                        << "Please check your database.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            qSeqs[thread_idx]->mapSequence(id, queryDbKey, querySeqData);
            matchers[thread_idx]->initQuery(qSeqs[thread_idx]);
            // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list
            std::vector<Matcher::result_t> swResults;
            size_t passedNum = 0;
            unsigned int rejected = 0;
            while (*data != '\0' && passedNum < maxAlnNum && rejected < maxRejected) {
                // DB key of the db sequence
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                dbKeys[thread_idx] = dbKey;
                // sequence are identical if qID == dbID  (needed to cluster really short sequences)
                const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB)) ? true : false;
                // map the database sequence
                char *dbSeqData = tseqdbr->getDataByDBKey(dbKeys[thread_idx]);
                if (dbSeqData == NULL) {
# pragma omp critical
                    {
                        Debug(Debug::ERROR) << "ERROR: Sequence " << dbKeys[thread_idx]
                                            << " is required in the prefiltering, but is not contained in the target sequence database!\n"
                                            <<
                                            "Please check your database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                //char *maskedDbSeq = seg[thread_idx]->maskseq(dbSeqData);
                dbSeqs[thread_idx]->mapSequence(-1, dbKeys[thread_idx], dbSeqData);
                // check if the sequences could pass the coverage threshold
                if ((((float) qSeqs[thread_idx]->L) / ((float) dbSeqs[thread_idx]->L) < covThr) ||
                    (((float) dbSeqs[thread_idx]->L) / ((float) qSeqs[thread_idx]->L) < covThr)) {
                    rejected++;
                    data = Util::skipLine(data);
                    continue;
                }
                // calculate Smith-Waterman alignment
                Matcher::result_t res = matchers[thread_idx]->getSWResult(dbSeqs[thread_idx], tseqdbr->getSize(),
                                                                          evalThr, swMode);
                alignmentsNum++;
                //set coverage and seqid if identity
                if (isIdentity) {
                    res.qcov = 1.0f;
                    res.dbcov = 1.0f;
                    res.seqId = 1.0f;
                }

                // check first if it is identity
                if (isIdentity ||
                    // general acceptance criteria
                    (res.eval <= evalThr && res.seqId >= seqIdThr && res.qcov >= covThr && res.dbcov >= covThr)) {
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
                realigner[thread_idx]->initQuery(qSeqs[thread_idx]);
                for (size_t result = 0; result < swResults.size(); result++) {
                    char *dbSeqData = tseqdbr->getDataByDBKey(swResults[result].dbKey);
                    dbSeqs[thread_idx]->mapSequence(-1, swResults[result].dbKey, dbSeqData);
                    Matcher::result_t res = realigner[thread_idx]->getSWResult(dbSeqs[thread_idx],
                                                                               tseqdbr->getSize(), 0.0,
                                                                               Matcher::SCORE_COV_SEQID);
                    swResults[result].backtrace = res.backtrace;
                    swResults[result].qStartPos = res.qStartPos;
                    swResults[result].qEndPos = res.qEndPos;
                    swResults[result].dbStartPos = res.dbStartPos;
                    swResults[result].dbEndPos = res.dbEndPos;
                }
            }
            // put the contents of the swResults list into ffindex DB
            std::stringstream swResultsString;
            for (size_t result = 0; result < swResults.size(); result++) {
                swResultsString << Matcher::resultToString(swResults[result], addBacktrace);
            }

            std::string swResultString = swResultsString.str();
            const char *swResultsData = swResultString.c_str();
            dbw.writeData(swResultsData, swResultString.length(), SSTR(qSeqs[thread_idx]->getDbKey()).c_str(),
                          thread_idx);
            swResults.clear();
        }
        prefdbr->remapData();
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

void Alignment::mergeAndRemoveTmpDatabases(const std::string &out, const std::string &outIndex,
                                           const std::vector<std::pair<std::string, std::string >> &files) {
    const char **datafilesNames = new const char *[files.size()];
    const char **indexFilesNames = new const char *[files.size()];
    for (size_t i = 0; i < files.size(); i++) {
        datafilesNames[i] = files[i].first.c_str();
        indexFilesNames[i] = files[i].second.c_str();
    }
    DBWriter::mergeResults(out.c_str(), outIndex.c_str(), datafilesNames, indexFilesNames, files.size());
    delete[] datafilesNames;
    delete[] indexFilesNames;
}
