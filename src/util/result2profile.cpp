#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "PSSMMasker.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "tantan.h"
#include "IndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2profile(int argc, const char **argv, const Command &command, bool returnAlnRes) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // default for result2profile to filter MSA
    par.filterMsa = 1;
    par.pca = 0.0;
    if (returnAlnRes) {
        par.PARAM_FILTER_MAX_SEQ_ID.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_QID.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_QSC.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_COV.removeCategory(MMseqsParameter::COMMAND_EXPERT);
        par.PARAM_FILTER_NDIFF.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile || returnAlnRes) ? par.evalThr : par.evalProfile;
    par.printParameters(command.cmd, argc, argv, *command.params);

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t dbFrom = 0;
    size_t dbSize = 0;
#ifdef HAVE_MPI
    resultReader.decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db4, par.db4Index, MMseqsMPI::rank);
#else
    dbSize = resultReader.getSize();
    std::pair<std::string, std::string> tmpOutput = std::make_pair(par.db4, par.db4Index);
#endif

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
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        qDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qDbr->readMmapedDataInMemory();
        }
    } else {
        qDbr = tDbr;
    }
    const unsigned int maxSequenceLength = std::max(tDbr->getMaxSeqLen(), qDbr->getMaxSeqLen());

    // qDbr->readMmapedDataInMemory();
    // make sure to touch target after query, so if there is not enough memory for the query, at least the targets
    // might have had enough space left to be residung in the page cache
    if (sameDatabase == false && templateDBIsIndex == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tDbr->readMmapedDataInMemory();
    }

    int type = Parameters::DBTYPE_HMM_PROFILE;
    if (returnAlnRes) {
        type = Parameters::DBTYPE_ALIGNMENT_RES;
    }
    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), localThreads, par.compressed, type);
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
    if (Parameters::isEqualDbtype(qDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }

    Debug(Debug::INFO) << "Query database size: " << qDbr->getSize() << " type: " << qDbr->getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << Parameters::getDbTypeName(targetSeqType) << "\n";

    const bool isFiltering = par.filterMsa != 0 || returnAlnRes;
    Debug::Progress progress(dbSize - dbFrom);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(qDbr->getDbtype(), maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        MultipleAlignment aligner(maxSequenceLength, &subMat);
        PSSMCalculator calculator(&subMat, maxSequenceLength, maxSetSize, par.pca, par.pcb);
        PSSMMasker masker(maxSequenceLength, probMatrix, subMat);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        Sequence centerSequence(maxSequenceLength, qDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        Sequence edgeSequence(maxSequenceLength, targetSeqType, &subMat, 0, false, false);

        char dbKey[255];
        const char *entry[255];
        char buffer[2048];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<std::vector<unsigned char>> seqSet;
        seqSet.reserve(300);

        std::string result;
        result.reserve((maxSequenceLength + 1) * Sequence::PROFILE_READIN_SIZE);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr->getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query sequence " << queryKey << "\n";
                continue;
            }
            centerSequence.mapSequence(queryId, queryKey, qDbr->getData(queryId, thread_idx), qDbr->getSeqLen(queryId));

            bool isQueryInit = false;
            char *data = resultReader.getData(id, thread_idx);
            while (*data != '\0') {
                Util::parseKey(data, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                // in the same database case, we have the query repeated
                if (key == queryKey && sameDatabase == true) {
                    data = Util::skipLine(data);
                    continue;
                }

                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                float evalue = 0.0;
                if (columns >= 4) {
                    evalue = strtod(entry[3], NULL);
                }

                if (evalue < par.evalProfile) {
                    const size_t edgeId = tDbr->getId(key);
                    if (edgeId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Sequence " << key << " does not exist in target sequence database\n";
                        EXIT(EXIT_FAILURE);
                    }
                    edgeSequence.mapSequence(edgeId, key, tDbr->getData(edgeId, thread_idx), tDbr->getSeqLen(edgeId));
                    seqSet.emplace_back(std::vector<unsigned char>(edgeSequence.numSequence, edgeSequence.numSequence + edgeSequence.L));

                    if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                        alnResults.emplace_back(Matcher::parseAlignmentRecord(data));
                    } else {
                        // Recompute if not all the backtraces are present
                        if (isQueryInit == false) {
                            matcher.initQuery(&centerSequence);
                            isQueryInit = true;
                        }
                        alnResults.emplace_back(matcher.getSWResult(&edgeSequence, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
                    }
                }
                data = Util::skipLine(data);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = aligner.computeMSA(&centerSequence, seqSet, alnResults, true);
            if (returnAlnRes == false) {
                alnResults.clear();
            }
            size_t filteredSetSize = isFiltering == false ? res.setSize
                    : filter.filter(res, alnResults, (int)(par.covMSAThr * 100), (int)(par.qid * 100), par.qsc, (int)(par.filterMaxSeqId * 100), par.Ndiff);
            //MultipleAlignment::print(res, &subMat);

            if (returnAlnRes) {
                // do not count query
                for (size_t i = 0; i < (filteredSetSize - 1); ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, alnResults[i], true);
                    result.append(buffer, len);
                }
                alnResults.clear();
            } else {
                for (size_t pos = 0; pos < res.centerLength; pos++) {
                    if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                        Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }

                PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength, (const char **) res.msaSequence, par.wg);
                if (par.maskProfile == true) {
                    masker.mask(centerSequence, pssmRes);
                }
                pssmRes.toBuffer(centerSequence, subMat, result);
            }
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
            result.clear();

            MultipleAlignment::deleteMSA(&res);
            seqSet.clear();
        }
    }
    resultWriter.close(returnAlnRes == false);
    resultReader.close();

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

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db4, par.db4Index, procs);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));

        }
        DBWriter::mergeResults(par.db4, par.db4Index, splitFiles);
    }
#endif

    if (MMseqsMPI::isMaster()) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    }

    return EXIT_SUCCESS;
}

int result2profile(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command, false);
}

int filterresult(int argc, const char **argv, const Command &command) {
    return result2profile(argc, argv, command, true);
}

