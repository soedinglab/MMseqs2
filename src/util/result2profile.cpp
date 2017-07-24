// Computes PSSMs from clustering or alignment result
// MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <algorithm>
#include <utility>

#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2profile(DBReader<unsigned int> &qDbr, Parameters &par, const std::string &outpath,
                   const size_t dbFrom, const size_t dbSize) {
    int localThreads = par.threads;
    if (static_cast<int>(qDbr.getSize()) <= par.threads) {
        localThreads = static_cast<int>(qDbr.getSize());
    }

#ifdef OPENMP
    omp_set_num_threads(localThreads);
#endif

    DBReader<unsigned int> *tDbr = NULL;
    DBReader<unsigned int> *tidxdbr = NULL;
    SequenceLookup *tSeqLookup = NULL;
    bool templateDBIsIndex = false;
    std::string scoringMatrixFile = par.scoringMatrixFile;

    unsigned int maxSequenceLength = 0;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        std::string indexDB = PrefilteringIndexReader::searchForIndex(par.db2);
        if (indexDB.length() > 0) {
            Debug(Debug::INFO) << "Use index  " << indexDB << "\n";

            tidxdbr = new DBReader<unsigned int>(indexDB.c_str(), (indexDB + ".index").c_str());
            tidxdbr->open(DBReader<unsigned int>::NOSORT);

            templateDBIsIndex = PrefilteringIndexReader::checkIfIndexFile(tidxdbr);
            if (templateDBIsIndex == true) {
                PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(tidxdbr);
                if (meta.maskMode == 1) {
                    Debug(Debug::WARNING) << "Cannot use masked index for profiles.\n";
                    templateDBIsIndex = false;
                } else {
                    PrefilteringIndexReader::printSummary(tidxdbr);
                    if (meta.maskMode == 0) {
                        tSeqLookup = PrefilteringIndexReader::getSequenceLookup(tidxdbr, par.noPreload == false);
                    } else if (meta.maskMode == 2) {
                        tSeqLookup = PrefilteringIndexReader::getUnmaskedSequenceLookup(tidxdbr, par.noPreload == false);
                    }
                    tDbr = PrefilteringIndexReader::openNewReader(tidxdbr, par.noPreload == false);
                    scoringMatrixFile = PrefilteringIndexReader::getSubstitutionMatrixName(tidxdbr);
                }
            }

            if (templateDBIsIndex == false) {
                tidxdbr->close();
                delete tidxdbr;
            }
        }

        if (templateDBIsIndex == false) {
            tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
            tDbr->open(DBReader<unsigned int>::NOSORT);
            if (par.noPreload == false) {
                tDbr->readMmapedDataInMemory();
                tDbr->mlock();
            }
        }

        unsigned int *lengths = qDbr.getSeqLens();
        for (size_t i = 0; i < qDbr.getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }

        lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }
    } else {
        tDbr = &qDbr;

        unsigned int *lengths = qDbr.getSeqLens();
        for (size_t i = 0; i < qDbr.getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(outpath.c_str(), std::string(outpath + ".index").c_str(), localThreads, DBWriter::BINARY_MODE);
    resultWriter.open();

    DBWriter *consensusWriter = NULL;
    if (!par.omitConsensus) {
        consensusWriter = new DBWriter(std::string(outpath + "_consensus").c_str(),
                     std::string(outpath + "_consensus.index").c_str(), localThreads);
        consensusWriter->open();
    }

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(scoringMatrixFile.c_str(), 2.0f, -0.2f);

    Debug(Debug::INFO) << "Start computing profiles.\n";
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, Matcher::GAP_OPEN, Matcher::GAP_EXTEND,
                                        true);
#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection);
        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, maxSetSize, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat);

        int sequenceType = Sequence::AMINO_ACIDS;
        if (par.queryProfile == true) {
            sequenceType = Sequence::HMM_PROFILE;
        }

        Sequence centerSequence(maxSequenceLength, subMat.aa2int, subMat.int2aa,
                                sequenceType, 0, false, par.compBiasCorrection);

#pragma omp for schedule(static)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif

            // Get the sequence from the queryDB
            unsigned int queryKey = resultReader.getDbKey(id);
            char *seqData = qDbr.getDataByDBKey(queryKey);
            if (seqData == NULL) {
                Debug(Debug::WARNING) << "Empty sequence " << id << ". Skipping.\n";
                continue;
            }
            centerSequence.mapSequence(0, queryKey, seqData);

            char *results = resultReader.getData(id);
            std::vector<Matcher::result_t> alnResults;
            std::vector<Sequence *> seqSet;
            char dbKey[255 + 1];
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                char *entry[255];
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // just add sequences if eval < thr. and if key is not the same as the query in case of sameDatabase
                if ((key != queryKey || sameDatabase == false)) {
                    if (columns > Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                        if (!(res.dbKey == centerSequence.getDbKey() && sameDatabase == true)) {
                            alnResults.push_back(res);
                        }
                    }

                    const size_t edgeId = tDbr->getId(key);
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLens(edgeId), subMat.aa2int, subMat.int2aa,
                                                          Sequence::AMINO_ACIDS, 0, false, false);

                    if (tSeqLookup != NULL) {
                        std::pair<const unsigned char*, const unsigned int> sequence = tSeqLookup->getSequence(edgeId);
                        edgeSequence->mapSequence(0, key, sequence);
                    } else {
                        char *dbSeqData = tDbr->getData(edgeId);
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
                        edgeSequence->mapSequence(0, key, dbSeqData);
                    }

                    seqSet.push_back(edgeSequence);
                }
                results = Util::skipLine(results);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = (alnResults.size() == seqSet.size())
                                               ? aligner.computeMSA(&centerSequence, seqSet, alnResults, true)
                                               : aligner.computeMSA(&centerSequence, seqSet, true);
            //MultipleAlignment::print(res, &subMat);

            size_t filteredSetSize = res.setSize;
            if (par.filterMsa == 1) {
                filter.filter(res.setSize, res.centerLength, static_cast<int>(par.cov * 100),
                              static_cast<int>(par.qid * 100), par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff,
                              (const char **) res.msaSequence, &filteredSetSize);
            }

            std::pair<const char *, std::string> pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength,
                                                                                         (const char **) res.msaSequence,
                                                                                         par.wg);
            char *data = (char *) pssmRes.first;
            size_t dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

            for (size_t i = 0; i < dataSize; i++) {
                // Avoid a null byte result
                data[i] = data[i] ^ 0x80;
            }

            //pssm.printProfile(res.centerLength);
            //calculator.printPSSM(res.centerLength);

            resultWriter.writeData(data, dataSize, queryKey, thread_idx);

            if (consensusWriter != NULL) {
                std::string consensusStr = pssmRes.second;
                consensusStr.push_back('\n');
                consensusWriter->writeData(consensusStr.c_str(), consensusStr.length(), queryKey, thread_idx);
            }

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }
    }

    // cleanup
    if (consensusWriter != NULL) {
        consensusWriter->close();
        delete consensusWriter;
    }
    resultWriter.close();

#ifndef HAVE_MPI
    if (par.earlyExit) {
        Debug(Debug::INFO) << "\nDone. Exiting early now.\n";
        _Exit(EXIT_SUCCESS);
    }
#endif

    resultReader.close();

    if (!sameDatabase) {
        tDbr->close();
        delete tDbr;
    }

    if (templateDBIsIndex == true) {
        delete tSeqLookup;
        tidxdbr->close();
        delete tidxdbr;
    }

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int result2profile(DBReader<unsigned int> &qDbr, Parameters &par,
                   const size_t dbFrom, const size_t dbSize) {
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";

    const std::string &outname = par.db4;
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outname, "", MMseqsMPI::rank);
    int status = result2profile(qDbr, par, tmpOutput.first, dbFrom, dbSize);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // master reduces results
    if (MMseqsMPI::isMaster() == 0) {
        const char *ext[2] = {"", "_consensus"};
        for (size_t i = 0; i < 2; i++) {
            std::vector<std::pair<std::string, std::string> > splitFiles;
            for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
                std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outname, "", procs);
                splitFiles.push_back(std::make_pair(tmpFile.first + ext[i], tmpFile.first + ext[i] + ".index"));

            }
            // merge output ffindex databases
            DBWriter::mergeResults(outname + ext[i], outname + ext[i] + ".index", splitFiles);
        }
    }

    return status;
}

int result2profile(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // default for result2profile filter MSA
    par.filterMsa = 1;
    par.parseParameters(argc, argv, command, 4);

    struct timeval start, end;
    gettimeofday(&start, NULL);


    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str());
    qDbr.open(DBReader<unsigned int>::NOSORT);

#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(qDbr.getAminoAcidDBSize(), qDbr.getSeqLens(), qDbr.getSize(),
                                     MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);

    int status = result2profile(qDbr, par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = result2profile(qDbr, par, par.db4, 0, qDbr.getSize());
#endif

    qDbr.close();

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m "
                          << (sec % 60) << "s\n";

    return status;
}
