#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "tantan.h"
#include "IndexReader.h"

#include <algorithm>
#include <utility>

#ifdef OPENMP
#include <omp.h>
#endif

int result2profile(DBReader<unsigned int> &resultReader, Parameters &par, const std::string &outDb, const std::string &outIndex, const size_t dbFrom, const size_t dbSize) {
    int localThreads = par.threads;
    if (static_cast<int>(resultReader.getSize()) <= par.threads) {
        localThreads = static_cast<int>(resultReader.getSize());
    }

    DBReader<unsigned int> *tDbr = NULL;
    IndexReader *tDbrIdx = NULL;
    bool templateDBIsIndex = false;

    int targetSeqType = -1;
    int targetDbtype = FileUtil::parseDbType(par.db2.c_str());
    if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        tDbrIdx = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tDbr = tDbrIdx->sequenceReader;
        templateDBIsIndex = true;
        targetSeqType = tDbr->getDbtype();
    }

    if (templateDBIsIndex == false) {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
        targetSeqType = tDbr->getDbtype();
    }

    DBReader<unsigned int> *qDbr = NULL;
    unsigned int maxSequenceLength = 0;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        qDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qDbr->readMmapedDataInMemory();
        }

        for (size_t i = 0; i < qDbr->getSize(); i++) {
            maxSequenceLength = std::max(static_cast<unsigned int>(qDbr->getSeqLen(i)), maxSequenceLength);
        }

    } else {
        qDbr = tDbr;
    }

    for (size_t i = 0; i < tDbr->getSize(); i++) {
        maxSequenceLength = std::max(static_cast<unsigned int>(tDbr->getSeqLen(i)), maxSequenceLength);
    }

    // qDbr->readMmapedDataInMemory();
    // make sure to touch target after query, so if there is not enough memory for the query, at least the targets
    // might have had enough space left to be residung in the page cache
    if (sameDatabase == false && templateDBIsIndex == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tDbr->readMmapedDataInMemory();
    }

    DBWriter resultWriter(outDb.c_str(), outIndex.c_str(), localThreads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultWriter.open();

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, -0.2f);
    ProbabilityMatrix probMatrix(subMat);
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);

    if (qDbr->getDbtype() == -1 || targetSeqType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(qDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }

    Debug(Debug::INFO) << "Query database size: " << qDbr->getSize() << " type: " << qDbr->getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << Parameters::getDbTypeName(targetSeqType) << "\n";

    const bool isFiltering = par.filterMsa != 0;
    int xAmioAcid = subMat.aa2num[static_cast<int>('X')];
    Debug::Progress progress(dbSize);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(qDbr->getDbtype(), maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, maxSetSize, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        Sequence centerSequence(maxSequenceLength, qDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        std::string result;
        result.reserve((maxSequenceLength + 1) * Sequence::PROFILE_READIN_SIZE);
        char *charSequence = new char[maxSequenceLength];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<Sequence *> seqSet;
        seqSet.reserve(300);

        char dbKey[255];
        const char *entry[255];
#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr->getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::ERROR) << "Sequence " << queryKey << " is not contained in the query sequence database\n";
                EXIT(EXIT_FAILURE);
            }
            centerSequence.mapSequence(0, queryKey, qDbr->getData(queryId, thread_idx), qDbr->getSeqLen(queryId));

            char *data = resultReader.getData(id, thread_idx);
            while (*data != '\0') {
                Util::parseKey(data, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                // in the same database case, we have the query repeated
                if ((key == queryKey && sameDatabase == true)) {
                    data = Util::skipLine(data);
                    continue;
                }

                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                float evalue = 0.0;
                if (columns >= 4) {
                    evalue = strtod(entry[3], NULL);
                }
                bool hasInclusionEval = (evalue < par.evalProfile);
                if (hasInclusionEval && columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Matcher::result_t res = Matcher::parseAlignmentRecord(data);
                    alnResults.push_back(res);
                }
                if (hasInclusionEval) {
                    const size_t edgeId = tDbr->getId(key);
                    if (edgeId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Sequence " << queryKey << " is not contained in the target sequence database\n";
                        EXIT(EXIT_FAILURE);
                    }
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLen(edgeId), targetSeqType, &subMat, 0, false, false);
                    edgeSequence->mapSequence(0, key, tDbr->getData(edgeId, thread_idx), tDbr->getSeqLen(edgeId));
                    seqSet.push_back(edgeSequence);
                }
                data = Util::skipLine(data);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = (alnResults.size() == seqSet.size())
                                               ? aligner.computeMSA(&centerSequence, seqSet, alnResults, true)
                                               : aligner.computeMSA(&centerSequence, seqSet, true);
            //MultipleAlignment::print(res, &subMat);
            alnResults.clear();

            size_t filteredSetSize = res.setSize;
            if (isFiltering) {
                filteredSetSize = filter.filter(res, static_cast<int>(par.covMSAThr * 100),
                              static_cast<int>(par.qid * 100), par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);
            }
            //MultipleAlignment::print(res, &subMat);

            for (size_t pos = 0; pos < res.centerLength; pos++) {
                if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                    Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength, (const char **) res.msaSequence, par.wg);
            if (par.maskProfile == true) {
                for (int i = 0; i < centerSequence.L; ++i) {
                    charSequence[i] = (unsigned char ) centerSequence.numSequence[i];
                }

                tantan::maskSequences(charSequence, charSequence + centerSequence.L,
                                      50 /*options.maxCycleLength*/,
                                      probMatrix.probMatrixPointers,
                                      0.005 /*options.repeatProb*/,
                                      0.05 /*options.repeatEndProb*/,
                                      0.9 /*options.repeatOffsetProbDecay*/,
                                      0, 0,
                                      0.9 /*options.minMaskProb*/,
                                      probMatrix.hardMaskTable);

                for (size_t pos = 0; pos < res.centerLength; pos++) {
                    if (charSequence[pos] == xAmioAcid) {
                        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                            pssmRes.prob[pos * Sequence::PROFILE_AA_SIZE + aa] = subMat.pBack[aa] * 0.5;
                        }
                        pssmRes.consensus[pos] = 'X';
                    }
                }
            }

            for (size_t pos = 0; pos < res.centerLength; pos++) {
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    result.push_back(Sequence::scoreMask(pssmRes.prob[pos * Sequence::PROFILE_AA_SIZE + aa]));
                }
                // write query, consensus sequence and neffM
                result.push_back(static_cast<unsigned char>(centerSequence.numSequence[pos]));
                result.push_back(static_cast<unsigned char>(subMat.aa2num[static_cast<int>(pssmRes.consensus[pos])]));
                unsigned char neff = MathUtil::convertNeffToChar(pssmRes.neffM[pos]);
                result.push_back(neff);
            }
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx);
            result.clear();

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                delete *it;
            }
            seqSet.clear();
        }
        delete[] charSequence;
    }
    resultWriter.close(true);

    if (!sameDatabase) {
        qDbr->close();
        delete qDbr;
    }
    if (tDbrIdx == NULL) {
        tDbr->close();
        delete tDbr;
    } else {
        delete tDbrIdx;
    }

    return EXIT_SUCCESS;
}

int result2profile(DBReader<unsigned int> &resultReader, Parameters &par, const size_t dbFrom, const size_t dbSize) {
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";

    const std::string outname = par.db4;
    const std::string outnameIndex = par.db4Index;

    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outname, outnameIndex, MMseqsMPI::rank);
    int status = result2profile(resultReader, par, tmpOutput.first, tmpOutput.second, dbFrom, dbSize);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outname, outnameIndex, procs);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));

        }
        DBWriter::mergeResults(outname, outnameIndex, splitFiles);
    }

    return status;
}

int result2profile(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // default for result2profile filter MSA
    par.filterMsa = 1;
    // no pseudo counts
    par.pca = 0.0;
    par.parseParameters(argc, argv, command, true, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    par.printParameters(command.cmd, argc, argv, *command.params);

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;
    resultReader.decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);

    int status = result2profile(resultReader, par, dbFrom, dbSize);
#else
    int status = result2profile(resultReader, par, par.db4, par.db4Index, 0, resultReader.getSize());
#endif

    if (MMseqsMPI::isMaster()) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    }

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    resultReader.close();

    return status;
}
