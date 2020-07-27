// Computes MSAs from clustering or alignment result

#include <string>
#include <vector>
#include <sstream>

#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBConcat.h"
#include "DBWriter.h"
#include "HeaderSummarizer.h"
#include "CompressedA3M.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2msa(Parameters &par, const std::string &resultData, const std::string &resultIndex,
               const size_t dbFrom, const size_t dbSize, DBConcat *referenceDBr = NULL) {
    if (par.compressMSA && referenceDBr == NULL) {
        Debug(Debug::ERROR) << "Need a sequence and header database for ca3m output!\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);
    if(par.preloadMode != Parameters::PRELOAD_MODE_MMAP){
        qDbr.readMmapedDataInMemory();
    }
    DBReader<unsigned int> queryHeaderReader(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    // NOSORT because the index should be in the same order as resultReader
    queryHeaderReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = &qDbr;
    DBReader<unsigned int> *tempateHeaderReader = &queryHeaderReader;

    size_t maxSequenceLength = 0;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);

        for (size_t i = 0; i < qDbr.getSize(); i++) {
            maxSequenceLength = std::max(qDbr.getEntryLen(i), maxSequenceLength);
        }

        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max(tDbr->getEntryLen(i), maxSequenceLength);
        }

        tempateHeaderReader = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tempateHeaderReader->open(DBReader<unsigned int>::NOSORT);
    } else {
        for (size_t i = 0; i < qDbr.getSize(); i++) {
        maxSequenceLength = std::max(qDbr.getEntryLen(i), maxSequenceLength);
        }
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t mode = par.compressed;
    int type = Parameters::DBTYPE_MSA_DB;
    if (par.compressMSA) {
        mode |= Parameters::WRITER_LEXICOGRAPHIC_MODE;
        type = Parameters::DBTYPE_CA3M_DB;
    }
    DBWriter resultWriter(resultData.c_str(), resultIndex.c_str(), par.threads, mode, type);
    resultWriter.open();

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, -0.2f);

    Debug(Debug::INFO) << "Start computing "
                       << (par.compressMSA ? "compressed" : "") << " multiple sequence alignments.\n";
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    if (qDbr.getDbtype() == -1 || tDbr->getDbtype() == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    if (Parameters::isEqualDbtype(qDbr.getDbtype(), Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE)){
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database.\n";
        EXIT(EXIT_FAILURE);
    }
    Debug(Debug::INFO) << "Query database size: "  << qDbr.getSize() << " type: " << qDbr.getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << tDbr->getDbTypeName() << "\n";

    const bool isFiltering = par.filterMsa != 0;
    Debug::Progress progress(dbSize-dbFrom);

#pragma omp parallel
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

        const char *entry[255];

#pragma omp  for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            // Get the sequence from the queryDB
            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            if (queryId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query sequence " << queryKey << ".\n";
                continue;
            }
            centerSequence.mapSequence(id, queryKey, qDbr.getData(queryId, thread_idx), qDbr.getSeqLen(queryId));

            // TODO: Do we still need this?
            if (centerSequence.L)
            {
                if(centerSequence.numSequence[centerSequence.L-1] == 20) // remove last in it is a *
                {
                    centerSequence.L--;
                }
            }
            char *centerSequenceHeader = queryHeaderReader.getDataByDBKey(queryKey, 0);

            char *results = resultReader.getData(id, thread_idx);
            std::vector<Matcher::result_t> alnResults;
            std::vector<Sequence *> seqSet;
            while (*results != '\0') {
                char dbKey[255 + 1];
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                // in the same database case, we have the query repeated
                if ((key == queryKey && sameDatabase == true)) {
                    results = Util::skipLine(results);
                    continue;
                }

                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                    alnResults.push_back(res);
                }

                const size_t edgeId = tDbr->getId(key);
                Sequence *edgeSequence = new Sequence(tDbr->getSeqLen(edgeId), tDbr->getDbtype(), &subMat, 0, false, false);


                char *dbSeqData = tDbr->getData(edgeId, thread_idx);
                if (dbSeqData == NULL) {
#pragma omp critical
                    {
                        Debug(Debug::ERROR) << "Sequence " << key << " is required in the prefiltering,"
                                            << "but is not contained in the target sequence database!\n"
                                            << "Please check your database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                edgeSequence->mapSequence(0, key, dbSeqData, tDbr->getSeqLen(edgeId));
                seqSet.push_back(edgeSequence);

                results = Util::skipLine(results);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = (alnResults.size() == seqSet.size())
                                               ? aligner.computeMSA(&centerSequence, seqSet, alnResults, !par.allowDeletion)
                                               : aligner.computeMSA(&centerSequence, seqSet, !par.allowDeletion);
            //MultipleAlignment::print(res, &subMat);

            alnResults = res.alignmentResults;
            size_t filteredSetSize = res.setSize;
            if (isFiltering) {
                filteredSetSize = filter.filter(res, static_cast<int>(par.covMSAThr * 100),
                              static_cast<int>(par.qid * 100), par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);
                filter.getKept(kept, res.setSize);
            }

            if (!par.compressMSA) {
                std::ostringstream msa;
                if (par.summarizeHeader) {
                    // gather headers for summary
                    std::vector<std::string> headers;
                    for (size_t i = 0; i < res.setSize; i++) {
                        if (i == 0) {
                            headers.push_back(std::string(centerSequenceHeader));
                        } else if (kept[i] == true) {
                            headers.push_back(tempateHeaderReader->getData(seqSet[i - 1]->getId(), thread_idx));
                        }
                    }

                    std::string summary = summarizer.summarize(headers);
                    msa << "#" << par.summaryPrefix.c_str() << "-" << queryKey << "|" << summary.c_str()
                        << "\n";
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
                    if (i == 0) {
                        key = queryKey;
                        header = centerSequenceHeader;
                    } else {
                        key = seqSet[i - 1]->getDbKey();
                        header = tempateHeaderReader->getDataByDBKey(key, thread_idx);
                    }
                    if (par.addInternalId) {
                        msa << "#" << key << "\n";
                    }

                    msa << ">" << header;

                    // need to allow insertion in the centerSequence
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        msa << ((aa < MultipleAlignment::NAA) ? subMat.num2aa[(int) aa] : '-');
                    }

                    msa << "\n";
                }

                std::string result = msa.str();
                resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
            } else {
                // Put the query sequence (master sequence) first in the alignment
                Matcher::result_t firstSequence;
                firstSequence.dbKey = queryKey;
                firstSequence.qStartPos = 0;
                firstSequence.dbStartPos = 0;
                firstSequence.backtrace = std::string(centerSequence.L, 'M'); // only matches

                alnResults.insert(alnResults.begin(), firstSequence);

                std::ostringstream msa;
                if (par.omitConsensus == false) {
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                            Debug(Debug::ERROR) <<  "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                            EXIT(EXIT_FAILURE);
                        }
                    }

                    PSSMCalculator::Profile pssmRes =
                            calculator.computePSSMFromMSA(filteredSetSize, res.centerLength,
                                                          (const char **) res.msaSequence, par.wg);
                    msa << ">consensus_" << queryHeaderReader.getDataByDBKey(queryKey, thread_idx) << pssmRes.consensus << "\n;";
                } else {
                    std::ostringstream centerSeqStr;
                    // Retrieve the master sequence
                    for (int pos = 0; pos < centerSequence.L; pos++) {
                        centerSeqStr << subMat.num2aa[centerSequence.numSequence[pos]];
                    }
                    msa << ">" << queryHeaderReader.getDataByDBKey(queryKey,  thread_idx) << centerSeqStr.str() << "\n;";
                }

                msa << CompressedA3M::fromAlignmentResult(alnResults, *referenceDBr);

                std::string result = msa.str();
                resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
            }

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }

        delete[] kept;
    }

    // cleanup
    resultWriter.close(true);
    resultReader.close();
    queryHeaderReader.close();
    qDbr.close();

    if (!sameDatabase) {
        tempateHeaderReader->close();
        delete tempateHeaderReader;
        tDbr->close();
        delete tDbr;
    }



    return EXIT_SUCCESS;
}

int result2msa(Parameters &par) {
    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    resultReader->open(DBReader<unsigned int>::NOSORT);
    size_t resultSize = resultReader->getSize();
    resultReader->close();
    delete resultReader;

    std::string outDb = par.db4;
    std::string outIndex = par.db4Index;
    DBConcat *referenceDBr = NULL;
    if (par.compressMSA) {
        std::string referenceSeqName(outDb);
        std::string referenceSeqIndexName(outDb);
        referenceSeqName.append("_sequence.ffdata");
        referenceSeqIndexName.append("_sequence.ffindex");

        // Use only 1 thread for concat to ensure the same order as the later header concat
        referenceDBr = new DBConcat(par.db1, par.db1Index, par.db2, par.db2Index,
                                    referenceSeqName, referenceSeqIndexName, 1);
        // When exporting in ca3m,
        // we need to have an access in SORT_BY_LINE
        // mode in order to keep track of the original
        // line number in the index file.
        referenceDBr->open(DBReader<unsigned int>::SORT_BY_LINE);

        std::string referenceHeadersName(outDb);
        std::string referenceHeadersIndexName(outDb);

        referenceHeadersName.append("_header.ffdata");
        referenceHeadersIndexName.append("_header.ffindex");

        // Use only 1 thread for concat to ensure the same order as the former sequence concat
        DBConcat referenceHeadersDBr(par.hdr1, par.hdr1Index, par.hdr2, par.hdr2Index,
                                     referenceHeadersName, referenceHeadersIndexName, 1);

        outDb.append("_ca3m.ffdata");
        outIndex = par.db4;
        outIndex.append("_ca3m.ffindex");
    }

    int status = result2msa(par, outDb, outIndex, 0, resultSize, referenceDBr);

    if (referenceDBr != NULL) {
        referenceDBr->close();
        delete referenceDBr;
    }
    return status;
}

int result2msa(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr->open(DBReader<unsigned int>::NOSORT);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    qDbr->decomposeDomainByAminoAcid(mpiRank, mpiNumProc, &dbFrom, &dbSize);
    qDbr->close();
    delete qDbr;


    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";

    DBConcat *referenceDBr = NULL;
    std::string outDb = par.db4;
    std::string outIndex = par.db4Index;
    if (par.compressMSA) {
        std::string referenceSeqName(outDb);
        std::string referenceSeqIndexName(outDb);
        referenceSeqName.append("_sequence.ffdata");
        referenceSeqIndexName.append("_sequence.ffindex");

        // Use only 1 thread for concat to ensure the same order as the later header concat
        referenceDBr = new DBConcat(par.db1, par.db1Index, par.db2, par.db2Index,
                                    referenceSeqName, referenceSeqIndexName, 1, MMseqsMPI::isMaster());

#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // When exporting in ca3m,
        // we need to have an access in SORT_BY_LINE
        // mode in order to keep track of the original
        // line number in the index file.
        referenceDBr->open(DBReader<unsigned int>::SORT_BY_LINE);

        if (MMseqsMPI::isMaster()) {
            std::string referenceHeadersName(outDb);
            std::string referenceHeadersIndexName(outDb);

            referenceHeadersName.append("_header.ffdata");
            referenceHeadersIndexName.append("_header.ffindex");

            // Use only 1 thread for concat to ensure the same order as the former sequence concat
            DBConcat referenceHeadersDBr(par.hdr1, par.hdr1Index, par.hdr2, par.hdr2Index,
                                         referenceHeadersName, referenceHeadersIndexName, 1);
        }

        outDb.append("_ca3m.ffdata");
        outIndex = par.db4;
        outIndex.append("_ca3m.ffindex");
    }

    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDb, outIndex, mpiRank);
    int status = result2msa(par, tmpOutput.first, tmpOutput.second, dbFrom, dbSize, referenceDBr);

    // close reader to reduce memory
    if (referenceDBr != NULL) {
        referenceDBr->close();
        delete referenceDBr;
    }

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for (unsigned int procs = 0; procs < mpiNumProc; procs++) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outDb, outIndex, procs);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));

        }
        DBWriter::mergeResults(outDb, outIndex, splitFiles, par.compressMSA);
    }

    return status;
}

int result2msa(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    // do not filter as default
    par.filterMsa = 0;
    par.pca = 0.0;
    par.parseParameters(argc, argv, command, true, 0, 0);

#ifdef HAVE_MPI
    int status = result2msa(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = result2msa(par);
#endif

    return status;
}
