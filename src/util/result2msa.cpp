#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "DBConcat.h"
#include "HeaderSummarizer.h"
#include "CompressedA3M.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2msa(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // do not filter by default
    par.filterMsa = 0;
    par.pca = 0.0;
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> queryHeaderReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    queryHeaderReader.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        qDbr.readMmapedDataInMemory();
        queryHeaderReader.readMmapedDataInMemory();
    }

    DBReader<unsigned int> *tDbr = &qDbr;
    DBReader<unsigned int> *targetHeaderReader = &queryHeaderReader;
    unsigned int maxSequenceLength = qDbr.getMaxSeqLen();

    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
        maxSequenceLength = std::max(qDbr.getMaxSeqLen(), tDbr->getMaxSeqLen());
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tDbr->readMmapedDataInMemory();
        }

        if (par.compressMSA == false) {
            targetHeaderReader = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            targetHeaderReader->open(DBReader<unsigned int>::NOSORT);

            if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
                targetHeaderReader->readMmapedDataInMemory();
            }
        }
    }

    DBConcat *seqConcat = NULL;
    DBReader<unsigned int> *refReader = NULL;
    std::string outDb = par.db4;
    std::string outIndex = par.db4Index;
    if (par.compressMSA) {
        std::string refData = outDb + "_sequence.ffdata";
        std::string refIndex = outDb + "_sequence.ffindex";
        // Use only 1 thread for concat to ensure the same order
        seqConcat = new DBConcat(par.db1, par.db1Index, par.db2, par.db2Index, refData, refIndex, 1, MMseqsMPI::isMaster());
        DBConcat hdrConcat(par.hdr1, par.hdr1Index, par.hdr2, par.hdr2Index, outDb + "_header.ffdata", outDb + "_header.ffindex", 1, MMseqsMPI::isMaster(), false, false, false, 2);

#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // When exporting in ca3m, we need to access with SORT_BY_LINE
        // mode in order to keep track of the original line number in the index file.
        refReader = new DBReader<unsigned int>(refData.c_str(), refIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
        refReader->open(DBReader<unsigned int>::SORT_BY_LINE);

        outDb = par.db4 + "_ca3m.ffdata";
        outIndex = par.db4 + "_ca3m.ffindex";
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t dbFrom = 0;
    size_t dbSize = 0;
#ifdef HAVE_MPI
    resultReader.decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDb, outIndex, MMseqsMPI::rank);
#else
    dbSize = resultReader.getSize();
    std::pair<std::string, std::string> tmpOutput = std::make_pair(outDb, outIndex);
#endif

   int localThreads = par.threads;
    if (static_cast<int>(resultReader.getSize()) <= par.threads) {
        localThreads = static_cast<int>(resultReader.getSize());
    }

    size_t mode = par.compressed;
    int type = Parameters::DBTYPE_MSA_DB;
    if (par.compressMSA) {
        mode |= Parameters::WRITER_LEXICOGRAPHIC_MODE;
        type = Parameters::DBTYPE_CA3M_DB;
    }
    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), localThreads, mode, type);
    resultWriter.open();

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, -0.2f);
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    if (qDbr.getDbtype() == -1 || tDbr->getDbtype() == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(qDbr.getDbtype(), Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }
    Debug(Debug::INFO) << "Query database size: " << qDbr.getSize() << " type: " << qDbr.getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << tDbr->getDbTypeName() << "\n";

    const bool isFiltering = par.filterMsa != 0;
    Debug::Progress progress(dbSize - dbFrom);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(qDbr.getDbtype(), maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, maxSetSize, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
        UniprotHeaderSummarizer summarizer;
        Sequence centerSequence(maxSequenceLength, qDbr.getDbtype(), &subMat, 0, false, par.compBiasCorrection);

        // which sequences where kept after filtering
        bool *kept = new bool[maxSetSize];
        for (size_t i = 0; i < maxSetSize; ++i) {
            kept[i] = 1;
        }

        char dbKey[255];
        const char *entry[255];

        std::vector<std::string> headers;
        headers.reserve(300);

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<Sequence *> seqSet;
        seqSet.reserve(300);

        std::string result;
        result.reserve(300 * 1024);

#pragma omp  for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query sequence " << queryKey << "\n";
                continue;
            }
            centerSequence.mapSequence(queryId, queryKey, qDbr.getData(queryId, thread_idx), qDbr.getSeqLen(queryId));

            // TODO: Do we still need this?
//            if (centerSequence.L) {
//                // remove last in it is a *
//                if(centerSequence.numSequence[centerSequence.L-1] == 20) {
//                    centerSequence.L--;
//                }
//            }

            size_t centerHeaderId = queryHeaderReader.getId(queryKey);
            if (centerHeaderId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query header " << queryKey << "\n";
                continue;
            }
            char* centerSequenceHeader = queryHeaderReader.getData(centerHeaderId, thread_idx);
            size_t centerHeaderLength = queryHeaderReader.getEntryLen(centerHeaderId) - 1;

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
                if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    alnResults.push_back(Matcher::parseAlignmentRecord(data));
                }

                const size_t edgeId = tDbr->getId(key);
                if (edgeId == UINT_MAX) {
                    Debug(Debug::ERROR) << "Sequence " << key << " does not exist in target sequence database\n";
                    EXIT(EXIT_FAILURE);
                }
                Sequence *edgeSequence = new Sequence(tDbr->getSeqLen(edgeId), tDbr->getDbtype(), &subMat, 0, false, false);
                edgeSequence->mapSequence(edgeId, key, tDbr->getData(edgeId, thread_idx), tDbr->getSeqLen(edgeId));
                seqSet.push_back(edgeSequence);

                data = Util::skipLine(data);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = (alnResults.size() == seqSet.size())
                                               ? aligner.computeMSA(&centerSequence, seqSet, alnResults, !par.allowDeletion)
                                               : aligner.computeMSA(&centerSequence, seqSet, !par.allowDeletion);
            //MultipleAlignment::print(res, &subMat);
            alnResults.clear();

            if (par.compressMSA == false) {
                if (isFiltering) {
                    filter.filter(res.setSize, res.centerLength, static_cast<int>(par.covMSAThr * 100), static_cast<int>(par.qid * 100), par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, (const char**)res.msaSequence, false);
                    filter.getKept(kept, res.setSize);
                }
                if (par.summarizeHeader) {
                    // gather headers for summary
                    for (size_t i = 0; i < res.setSize; i++) {
                        if (i == 0) {
                            headers.emplace_back(centerSequenceHeader, centerHeaderLength);
                        } else if (kept[i] == true) {
                            unsigned int id = seqSet[i - 1]->getId();
                            char *header = targetHeaderReader->getData(id, thread_idx);
                            size_t length = targetHeaderReader->getEntryLen(id) - 1;
                            headers.emplace_back(header, length);
                        }
                    }
                    result.append(1, '#');
                    result.append(par.summaryPrefix);
                    result.append(1, '-');
                    result.append(SSTR(queryKey));
                    result.append(1, '|');
                    result.append(summarizer.summarize(headers));
                    result.append(1, '\n');
                    headers.clear();
                }

                size_t start = 0;
                if (par.skipQuery == true) {
                    start = 1;
                }
                for (size_t i = start; i < res.setSize; i++) {
                    if (kept[i] == false) {
                        continue;
                    }

                    unsigned int key;
                    char *header;
                    size_t length;
                    if (i == 0) {
                        key = queryKey;
                        header = centerSequenceHeader;
                        length = centerHeaderLength - 1;
                    } else {
                        key = seqSet[i - 1]->getDbKey();
                        size_t id = seqSet[i - 1]->getId();
                        header = targetHeaderReader->getData(id, thread_idx);
                        length = targetHeaderReader->getEntryLen(id) - 1;
                    }

                    if (par.addInternalId) {
                        result.append(1, '#');
                        result.append(SSTR(key));
                        result.append(1, '\n');
                    }

                    result.append(1, '>');
                    result.append(header, length);
                    // need to allow insertion in the centerSequence
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        result.append(1, ((aa < MultipleAlignment::NAA) ? subMat.num2aa[(int) aa] : '-'));
                    }
                    result.append(1, '\n');
                }
            } else {
                size_t filteredSetSize = res.setSize;
                if (isFiltering) {
                    filteredSetSize = filter.filter(res, static_cast<int>(par.covMSAThr * 100), static_cast<int>(par.qid * 100), par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);
                }
                if (par.omitConsensus == false) {
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                            Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                            EXIT(EXIT_FAILURE);
                        }
                    }

                    PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength, (const char **) res.msaSequence, par.wg);
                    result.append(">consensus_");
                    result.append(centerSequenceHeader, centerHeaderLength);
                    result.append(pssmRes.consensus);
                    result.append("\n;");
                } else {
                    result.append(1, '>');
                    result.append(centerSequenceHeader, centerHeaderLength);
                    // Retrieve the master sequence
                    for (int pos = 0; pos < centerSequence.L; pos++) {
                        result.append(1, subMat.num2aa[centerSequence.numSequence[pos]]);
                    }
                    result.append("\n;");
                }

                Matcher::result_t queryAln;
                unsigned int newQueryKey = seqConcat->dbAKeyMap(queryKey);
                queryAln.qStartPos = 0;
                queryAln.dbStartPos = 0;
                queryAln.backtrace = std::string(centerSequence.L, 'M'); // only matches
                CompressedA3M::hitToBuffer(refReader->getId(newQueryKey), queryAln, result);
                for (size_t i = 0; i < res.alignmentResults.size(); ++i) {
                    unsigned int key = res.alignmentResults[i].dbKey;
                    unsigned int targetKey = seqConcat->dbBKeyMap(key);
                    unsigned int targetId = refReader->getId(targetKey);
                    CompressedA3M::hitToBuffer(targetId, res.alignmentResults[i], result);
                }
            }
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
            result.clear();

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                delete *it;
            }
            seqSet.clear();
        }

        delete[] kept;
    }
    resultWriter.close(true);
    resultReader.close();
    queryHeaderReader.close();
    qDbr.close();
    if (!sameDatabase) {
        if (par.compressMSA == false) {
            targetHeaderReader->close();
            delete targetHeaderReader;
        }
        tDbr->close();
        delete tDbr;
    }

    if (refReader != NULL) {
        refReader->close();
        delete refReader;
    }
    if (seqConcat != NULL) {
        delete seqConcat;
    }

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outDb, outIndex, procs);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));

        }
        DBWriter::mergeResults(outDb, outIndex, splitFiles, par.compressMSA);
    }
#endif
    return EXIT_SUCCESS;
}

