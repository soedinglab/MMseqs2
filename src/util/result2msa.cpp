#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "IndexReader.h"
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
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
    std::vector<int> qid_vec;
    for (size_t qid_idx = 0; qid_idx < qid_str_vec.size(); qid_idx++) {
        float qid_float = strtod(qid_str_vec[qid_idx].c_str(), NULL);
        qid_vec.push_back(static_cast<int>(qid_float*100));
    }
    std::sort(qid_vec.begin(), qid_vec.end());

    const bool isCA3M = par.msaFormatMode == Parameters::FORMAT_MSA_CA3M || par.msaFormatMode == Parameters::FORMAT_MSA_CA3M_CONSENSUS;
    const bool shouldWriteNullByte = par.msaFormatMode != Parameters::FORMAT_MSA_STOCKHOLM_FLAT;

    DBReader<unsigned int> *tDbr = NULL;
    IndexReader *tDbrIdx = NULL;
    DBReader<unsigned int> *targetHeaderReader = NULL;
    IndexReader *targetHeaderReaderIdx = NULL;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;

    if (Parameters::isEqualDbtype(FileUtil::parseDbType(par.db2.c_str()), Parameters::DBTYPE_INDEX_DB)) {
        if (isCA3M == true) {
            Debug(Debug::ERROR) << "Cannot use result2msa with indexed target database for CA3M output\n";
            return EXIT_FAILURE;
        }
        uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        tDbrIdx = new IndexReader(par.db2, par.threads,
                                  extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                                  (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tDbr = tDbrIdx->sequenceReader;
        targetHeaderReaderIdx = new IndexReader(par.db2, par.threads,
                                                extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC ? IndexReader::SRC_HEADERS : IndexReader::HEADERS,
                                                (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        targetHeaderReader = targetHeaderReaderIdx->sequenceReader;
    } else {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tDbr->readMmapedDataInMemory();
        }
        if (isCA3M == false || sameDatabase) {
            targetHeaderReader = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            targetHeaderReader->open(DBReader<unsigned int>::NOSORT);
            if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
                targetHeaderReader->readMmapedDataInMemory();
            }
        }
    }

    DBReader<unsigned int> *qDbr = NULL;
    DBReader<unsigned int> *queryHeaderReader = NULL;
    if (!sameDatabase) {
        qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        qDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qDbr->readMmapedDataInMemory();
        }
        queryHeaderReader = new DBReader<unsigned int>(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        queryHeaderReader->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            queryHeaderReader->readMmapedDataInMemory();
        }
    } else {
        qDbr = tDbr;
        queryHeaderReader = targetHeaderReader;
    }
    const unsigned int maxSequenceLength = std::max(tDbr->getMaxSeqLen(), qDbr->getMaxSeqLen());

    DBConcat *seqConcat = NULL;
    DBReader<unsigned int> *refReader = NULL;
    std::string outDb = par.db4;
    std::string outIndex = par.db4Index;
    if (isCA3M == true) {
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

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, resultReader.getSize()), (size_t)1);
#endif

    size_t mode = par.compressed;
    int type = Parameters::DBTYPE_MSA_DB;
    if (isCA3M == true) {
        mode |= Parameters::WRITER_LEXICOGRAPHIC_MODE;
        type = Parameters::DBTYPE_CA3M_DB;
    } else if (par.msaFormatMode == Parameters::FORMAT_MSA_STOCKHOLM_FLAT) {
        type = Parameters::DBTYPE_OMIT_FILE;
    }
    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), localThreads, mode, type);
    resultWriter.open();

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
    if (qDbr->getDbtype() == -1 || tDbr->getDbtype() == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(qDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }
    Debug(Debug::INFO) << "Query database size: " << qDbr->getSize() << " type: " << qDbr->getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << tDbr->getDbTypeName() << "\n";

    const bool isFiltering = par.filterMsa != 0;
    Debug::Progress progress(dbSize - dbFrom);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(qDbr->getDbtype(), tDbr->getDbtype(), maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection,
                        par.compBiasCorrectionScale, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
        MultipleAlignment aligner(maxSequenceLength, &subMat);
        PSSMCalculator calculator(
            &subMat, maxSequenceLength, maxSetSize, par.pcmode, par.pca, par.pcb
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid()
            , par.gapPseudoCount
#endif
        );
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
        UniprotHeaderSummarizer summarizer;
        Sequence centerSequence(maxSequenceLength, qDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        Sequence edgeSequence(maxSequenceLength, tDbr->getDbtype(), &subMat, 0, false, false);

        // which sequences where kept after filtering
        bool *kept = new bool[maxSetSize];
        for (size_t i = 0; i < maxSetSize; ++i) {
            kept[i] = 1;
        }

        char dbKey[255];
        const char *entry[255];
        std::string accession;

        std::vector<std::string> headers;
        headers.reserve(300);

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<std::vector<unsigned char>> seqSet;
        seqSet.reserve(300);

        std::vector<unsigned int> seqKeys;
        seqKeys.reserve(300);

        std::string result;
        result.reserve(300 * 1024);
        char buffer[1024 + 32768*4];

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
            // TODO: Do we still need this?
//            if (centerSequence.L) {
//                // remove last in it is a *
//                if(centerSequence.numSequence[centerSequence.L-1] == 20) {
//                    centerSequence.L--;
//                }
//            }

            size_t centerHeaderId = queryHeaderReader->getId(queryKey);
            if (centerHeaderId == UINT_MAX) {
                Debug(Debug::WARNING) << "Invalid query header " << queryKey << "\n";
                continue;
            }

            char *centerSequenceHeader = queryHeaderReader->getData(centerHeaderId, thread_idx);
            size_t centerHeaderLength = queryHeaderReader->getEntryLen(centerHeaderId) - 1;

            if (par.msaFormatMode == Parameters::FORMAT_MSA_STOCKHOLM_FLAT) {
                accession = Util::parseFastaHeader(centerSequenceHeader);
            }

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

                const size_t edgeId = tDbr->getId(key);
                if (edgeId == UINT_MAX) {
                    Debug(Debug::ERROR) << "Sequence " << key << " does not exist in target sequence database\n";
                    EXIT(EXIT_FAILURE);
                }
                edgeSequence.mapSequence(edgeId, key, tDbr->getData(edgeId, thread_idx), tDbr->getSeqLen(edgeId));
                seqSet.emplace_back(std::vector<unsigned char>(edgeSequence.numSequence, edgeSequence.numSequence + edgeSequence.L));
                seqKeys.emplace_back(key);

                const size_t columns = Util::getWordsOfLine(data, entry, 255);
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
                data = Util::skipLine(data);
            }

            MultipleAlignment::MSAResult res = aligner.computeMSA(&centerSequence, seqSet, alnResults, !par.allowDeletion);
            //MultipleAlignment::print(res, &subMat);

            if (par.msaFormatMode == Parameters::FORMAT_MSA_FASTADB || par.msaFormatMode == Parameters::FORMAT_MSA_FASTADB_SUMMARY) {
                if (isFiltering) {
                    filter.filter(res.setSize, res.centerLength, static_cast<int>(par.covMSAThr * 100), qid_vec, par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable, (const char **) res.msaSequence, false);
                    filter.getKept(kept, res.setSize);
                }
                if (par.msaFormatMode == Parameters::FORMAT_MSA_FASTADB_SUMMARY) {
                    // gather headers for summary
                    for (size_t i = 0; i < res.setSize; i++) {
                        if (i == 0) {
                            headers.emplace_back(centerSequenceHeader, centerHeaderLength);
                        } else if (kept[i] == true) {
                            unsigned int key = seqKeys[i - 1];
                            size_t id = targetHeaderReader->getId(key);
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

                    char *header;
                    size_t length;
                    if (i == 0) {
                        header = centerSequenceHeader;
                        length = centerHeaderLength;
                    } else {
                        unsigned int key = seqKeys[i - 1];
                        size_t id = targetHeaderReader->getId(key);
                        header = targetHeaderReader->getData(id, thread_idx);
                        length = targetHeaderReader->getEntryLen(id) - 1;
                    }
                    bool isOnlyGap = true;
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        if (aa != MultipleAlignment::GAP) {
                            isOnlyGap = false;
                            break;
                        }
                    }

                    result.append(1, '>');
                    if(isOnlyGap) {
                        result.append("DUMMY\n");
                    }else{
                        result.append(header, length);
                    }
                    // need to allow insertion in the centerSequence
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        result.append(1, ((aa < MultipleAlignment::GAP) ? subMat.num2aa[(int) aa] : '-'));
                    }
                    result.append(1, '\n');
                }
            } else if (par.msaFormatMode == Parameters::FORMAT_MSA_STOCKHOLM_FLAT) {
                if (isFiltering) {
                    filter.filter(res.setSize, res.centerLength, static_cast<int>(par.covMSAThr * 100), qid_vec, par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable, (const char **) res.msaSequence, false);
                    filter.getKept(kept, res.setSize);
                }

                result.append("# STOCKHOLM 1.0\n");
                size_t start = 0;
                if (par.skipQuery == true) {
                    start = 1;
                    result.append("#=GF ID ");
                    result.append(Util::parseFastaHeader(centerSequenceHeader));
                    result.append(1, '\n');
                }
                for (size_t i = start; i < res.setSize; i++) {
                    if (kept[i] == false) {
                        continue;
                    }

                    const char *header;
                    bool isOnlyGap = true;
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        if (aa != MultipleAlignment::GAP) {
                            isOnlyGap = false;
                            break;
                        }
                    }


                    if (i == 0) {
                        if(isOnlyGap) {
                            header = "DUMMY";
                        }else{
                            header = centerSequenceHeader;
                        }
                    } else {
                        if(isOnlyGap) {
                            header = "DUMMY";
                        }else {
                            unsigned int key = seqKeys[i - 1];
                            size_t id = targetHeaderReader->getId(key);
                            header = targetHeaderReader->getData(id, thread_idx);
                        }
                    }
                    accession = Util::parseFastaHeader(header);

                    result.append(accession);
                    result.append(1, ' ');
                    // need to allow insertion in the centerSequence
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        result.append(1, ((aa < MultipleAlignment::GAP) ? subMat.num2aa[(int) aa] : '-'));
                    }
                    result.append(1, '\n');
                }
                result.append("//\n");
            } else if (par.msaFormatMode == Parameters::FORMAT_MSA_A3M || par.msaFormatMode == Parameters::FORMAT_MSA_A3M_ALN_INFO) {
                if (isFiltering) {
                    filter.filter(res.setSize, res.centerLength, static_cast<int>(par.covMSAThr * 100), qid_vec, par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable, (const char **) res.msaSequence, false);
                    filter.getKept(kept, res.setSize);
                }

                size_t start = (par.skipQuery == true) ? 1 : 0;
                for (size_t i = start; i < res.setSize; i++) {
                    if (kept[i] == false) {
                        continue;
                    }

                    result.push_back('>');
                    bool isOnlyGap = true;
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        char aa = res.msaSequence[i][pos];
                        if (aa != MultipleAlignment::GAP) {
                            isOnlyGap = false;
                            break;
                        }
                    }


                    if (i == 0) {
                        if(isOnlyGap){
                            result.append("DUMMY");
                        }else{
                            result.append(Util::parseFastaHeader(centerSequenceHeader));
                        }
                    } else {
                        unsigned int key = seqKeys[i - 1];
                        size_t id = targetHeaderReader->getId(key);
                        if(isOnlyGap){
                            result.append("DUMMY");
                        }else {
                            result.append(Util::parseFastaHeader(targetHeaderReader->getData(id, thread_idx)));
                        }
                        if (par.msaFormatMode == Parameters::FORMAT_MSA_A3M_ALN_INFO) {
                            size_t len = Matcher::resultToBuffer(buffer, alnResults[i - 1], false);
                            char* data = buffer;
                            data += Util::skipNoneWhitespace(data);
                            result.append(data, len - (data - buffer) - 1);
                        }
                    }
                    result.push_back('\n');

                    // need to allow insertion in the centerSequence
                    if(i == 0){
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            char aa = res.msaSequence[i][pos];
                            result.append(1, ((aa < MultipleAlignment::GAP) ? subMat.num2aa[(int) aa] : '-'));
                        }
                        result.append(1, '\n');
                    }else{
                        const std::vector<unsigned char> & seq = seqSet[i-1];
                        int seqStartPos = alnResults[i-1].dbStartPos;
                        size_t seqPos = 0;
                        const std::string & bt = alnResults[i-1].backtrace;
                        size_t btPos = 0;

                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            char aa = res.msaSequence[i][pos];

                            if(aa>=MultipleAlignment::GAP){
                                result.push_back('-');
                            }else if(aa<MultipleAlignment::GAP){
                                result.push_back( subMat.num2aa[(int) aa]);
                                btPos++;
                                seqPos++;
                            }
                            // skip insert
                            while(btPos < bt.size() && bt[btPos] == 'I') { btPos++;}

                            // add lower case deletions
                            while(btPos < bt.size() && bt[btPos] == 'D') {
                                result.push_back(tolower(subMat.num2aa[seq[seqStartPos+seqPos]]));
                                btPos++;
                                seqPos++;
                            }
                        }
                        result.append(1, '\n');
                    }
                }
            } else if (isCA3M == true) {
                size_t filteredSetSize = res.setSize;
                if (isFiltering) {
                    filteredSetSize = filter.filter(res, alnResults, static_cast<int>(par.covMSAThr * 100), qid_vec, par.qsc, static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff, par.filterMinEnable);
                }
                if (par.formatAlignmentMode == Parameters::FORMAT_MSA_CA3M_CONSENSUS) {
                    for (size_t pos = 0; pos < res.centerLength; pos++) {
                        if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                            Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                            EXIT(EXIT_FAILURE);
                        }
                    }

                    PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(
                        filteredSetSize, res.centerLength, (const char **) res.msaSequence,
#ifdef GAP_POS_SCORING
                        alnResults,
#endif
                        par.wg);
                    result.append(">consensus_");
                    result.append(centerSequenceHeader, centerHeaderLength);
                    for (int pos = 0; pos < centerSequence.L; pos++) {
                        result.push_back(subMat.num2aa[pssmRes.consensus[pos]]);
                    }
                    result.append("\n;");
                } else {
                    result.append(1, '>');
                    result.append(centerSequenceHeader, centerHeaderLength);
                    // Retrieve the master sequence
                    for (int pos = 0; pos < centerSequence.L; pos++) {
                        result.push_back(subMat.num2aa[centerSequence.numSequence[pos]]);
                    }
                    result.append("\n;");
                }
                Matcher::result_t queryAln;
                unsigned int newQueryKey = seqConcat->dbAKeyMap(queryKey);
                queryAln.qStartPos = 0;
                queryAln.dbStartPos = 0;
                queryAln.backtrace = std::string(centerSequence.L, 'M'); // only matches
                CompressedA3M::hitToBuffer(refReader->getId(newQueryKey), queryAln, result);
                for (size_t i = 0; i < alnResults.size(); ++i) {
                    unsigned int key = alnResults[i].dbKey;
                    unsigned int targetKey = seqConcat->dbBKeyMap(key);
                    unsigned int targetId = refReader->getId(targetKey);
                    CompressedA3M::hitToBuffer(targetId, alnResults[i], result);
                }
            }
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx, shouldWriteNullByte);
            result.clear();

            MultipleAlignment::deleteMSA(&res);
            seqSet.clear();
            seqKeys.clear();
            alnResults.clear();
        }

        delete[] kept;
    }
    resultWriter.close(true);
    if (shouldWriteNullByte == false) {
        FileUtil::remove(resultWriter.getIndexFileName());
    }
    resultReader.close();

    if (!sameDatabase) {
        qDbr->close();
        delete qDbr;
        queryHeaderReader->close();
        delete queryHeaderReader;
    }
    if (tDbrIdx == NULL) {
        tDbr->close();
        delete tDbr;
        if (targetHeaderReader != NULL) {
            targetHeaderReader->close();
            delete targetHeaderReader;
        }
    } else {
        delete tDbrIdx;
        delete targetHeaderReaderIdx;
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
        DBWriter::mergeResults(outDb, outIndex, splitFiles, isCA3M);
    }
#endif
    return EXIT_SUCCESS;
}

