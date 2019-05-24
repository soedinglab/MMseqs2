// Computes PSSMs from clustering or alignment result
// MMseqs just stores the position specific score in 1 byte
#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <tantan.h>
#include "IndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

int result2profile(DBReader<unsigned int> &resultReader, Parameters &par, const std::string &outpath,
                   const size_t dbFrom, const size_t dbSize) {
    int localThreads = par.threads;
    if (static_cast<int>(resultReader.getSize()) <= par.threads) {
        localThreads = static_cast<int>(resultReader.getSize());
    }

    DBReader<unsigned int> *tDbr = NULL;
    IndexReader * tDbrIdx = NULL;
    bool templateDBIsIndex = false;

    int targetSeqType = -1;
    int targetDbtype = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
        bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
        tDbrIdx = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        tDbr = tDbrIdx->sequenceReader;
        tDbr->getDbtype();
        templateDBIsIndex = true;
        targetSeqType = tDbr->getDbtype();
    }

    if (templateDBIsIndex == false) {
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
        targetSeqType = tDbr->getDbtype();
    }

    DBReader<unsigned int> *qDbr = NULL;
    unsigned int maxSequenceLength = 0;
    const bool sameDatabase = (par.db1.compare(par.db2) == 0) ? true : false;
    if (!sameDatabase) {
        qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qDbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qDbr->readMmapedDataInMemory();
        }

        unsigned int *lengths = qDbr->getSeqLens();
        for (size_t i = 0; i < qDbr->getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }

        lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }
    } else {
        qDbr = tDbr;

        unsigned int *lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }
    }

    // qDbr->readMmapedDataInMemory();
    // make sure to touch target after query, so if there is not enough memory for the query, at least the targets
    // might have had enough space left to be residung in the page cache
    if (sameDatabase == false && templateDBIsIndex == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tDbr->readMmapedDataInMemory();
    }

    DBWriter resultWriter(outpath.c_str(), std::string(outpath + ".index").c_str(), localThreads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultWriter.open();

    DBWriter *consensusWriter = NULL;
    if (!par.omitConsensus) {
        consensusWriter = new DBWriter(std::string(outpath + "_consensus").c_str(),
                                       std::string(outpath + "_consensus.index").c_str(), localThreads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
        consensusWriter->open();
    }

    // + 1 for query
    size_t maxSetSize = resultReader.maxCount('\n') + 1;

    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    ProbabilityMatrix probMatrix(subMat);

    Debug(Debug::INFO) << "Start computing profiles\n";
    EvalueComputation evalueComputation(tDbr->getAminoAcidDBSize(), &subMat, par.gapOpen, par.gapExtend);

    if (qDbr->getDbtype() == -1 || targetSeqType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(qDbr->getDbtype(), Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE)){
        Debug(Debug::ERROR) << "Only the query OR the target database can be a profile database\n";
        return EXIT_FAILURE;
    }

    Debug(Debug::INFO) << "Query database size: "  << qDbr->getSize() << " type: " << qDbr->getDbTypeName() << "\n";
    Debug(Debug::INFO) << "Target database size: " << tDbr->getSize() << " type: " << DBReader<unsigned int>::getDbTypeName(targetSeqType) << "\n";


    const bool isFiltering = par.filterMsa != 0;
    int xAmioAcid = subMat.aa2int[(int)'X'];
    Debug::Progress progress(dbSize);

#pragma omp parallel num_threads(localThreads)
    {
        Matcher matcher(qDbr->getDbtype(), maxSequenceLength, &subMat, &evalueComputation, par.compBiasCorrection, par.gapOpen, par.gapExtend);
        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, maxSetSize, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat, par.gapOpen, par.gapExtend);
        Sequence centerSequence(maxSequenceLength, qDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        std::string result;
        result.reserve(par.maxSeqLen * Sequence::PROFILE_READIN_SIZE * sizeof(char));
        char *charSequence = new char[maxSequenceLength];

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        const char *entry[255];

#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            // Get the sequence from the queryDB
            unsigned int queryKey = resultReader.getDbKey(id);

            size_t queryId = qDbr->getId(queryKey);

            char *dbSeqData = qDbr->getData(queryId, thread_idx);
            if (dbSeqData == NULL) {
#pragma omp critical
                {
                    Debug(Debug::ERROR) << "Sequence " << queryKey << " is required in the database,"
                                        << "but is not contained in the query sequence database!\n"
                                        << "Please check your database.\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            centerSequence.mapSequence(0, queryKey, dbSeqData);

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
                float evalue = 0.0;
                if(columns >= 4){
                    evalue = strtod (entry[3], NULL);
                }
                bool hasInclusionEval  = (evalue < par.evalProfile);
                if (hasInclusionEval && columns > Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                    Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                    alnResults.push_back(res);
                }
                if(hasInclusionEval){
                    const size_t edgeId = tDbr->getId(key);
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLens(edgeId),
                                                          targetSeqType, &subMat, 0, false, false);

                    char *dbSeqData = tDbr->getData(edgeId, thread_idx);
                    if (dbSeqData == NULL) {
#pragma omp critical
                        {
                            Debug(Debug::ERROR) << "Sequence " << key << " is required in the database,"
                                                << "but is not contained in the target sequence database!\n"
                                                << "Please check your database.\n";
                            EXIT(EXIT_FAILURE);
                        }
                    }
                    edgeSequence->mapSequence(0, key, dbSeqData);

                    seqSet.push_back(edgeSequence);
                }
                results = Util::skipLine(results);
            }

            // Recompute if not all the backtraces are present
            MultipleAlignment::MSAResult res = (alnResults.size() == seqSet.size())
                                               ? aligner.computeMSA(&centerSequence, seqSet, alnResults, true)
                                               : aligner.computeMSA(&centerSequence, seqSet, true);


//            MultipleAlignment::print(res, &subMat);

            size_t filteredSetSize = res.setSize;
            if (isFiltering) {
                filter.filter(res.setSize, res.centerLength, static_cast<int>(par.covMSAThr * 100),
                              static_cast<int>(par.qid * 100), par.qsc,
                              static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff,
                              (const char **) res.msaSequence, &filteredSetSize);
                filter.shuffleSequences((const char **) res.msaSequence, res.setSize);
            }

//                        MultipleAlignment::print(res, &subMat);


            for (size_t pos = 0; pos < res.centerLength; pos++) {
                if (res.msaSequence[0][pos] == MultipleAlignment::GAP) {
                    Debug(Debug::ERROR) <<  "Error in computePSSMFromMSA. First sequence of MSA is not allowed to contain gaps.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(filteredSetSize, res.centerLength,
                                                                            (const char **) res.msaSequence,
                                                                            par.wg);

            if(par.maskProfile == true){
                for (int i = 0; i < centerSequence.L; ++i) {
                    charSequence[i] = (char) centerSequence.int_sequence[i];
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

                for(size_t pos = 0; pos < res.centerLength; pos++) {
                    if(charSequence[pos] == xAmioAcid) {
                        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                            pssmRes.prob[pos*Sequence::PROFILE_AA_SIZE + aa] = subMat.pBack[aa]*0.5;
                        }
                        pssmRes.consensus[pos]='X';
                    }
                }
            }



            for(size_t pos = 0; pos <  res.centerLength; pos++){
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    result.push_back(Sequence::scoreMask(pssmRes.prob[pos*Sequence::PROFILE_AA_SIZE + aa]));
                }
                // write query, consensus sequence and neffM
                result.push_back(static_cast<unsigned char>(centerSequence.int_sequence[pos]));
                result.push_back(static_cast<unsigned char>(subMat.aa2int[static_cast<int>(pssmRes.consensus[pos])]));
                unsigned char neff = MathUtil::convertNeffToChar(pssmRes.neffM[pos]);
                result.push_back(neff);
            }


            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx);
            result.clear();
            if (consensusWriter != NULL) {
                std::string consensusStr = pssmRes.consensus;
                consensusStr.push_back('\n');
                consensusWriter->writeData(consensusStr.c_str(), consensusStr.length(), queryKey, thread_idx);
            }
            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }
        delete [] charSequence;
    }

    // cleanup
    if (consensusWriter != NULL) {
        consensusWriter->close(true);
        delete consensusWriter;
    }
    resultWriter.close(true);

    if (!sameDatabase) {
        qDbr->close();
        delete qDbr;
    }
    if(tDbrIdx==NULL) {
        tDbr->close();
        delete tDbr;
    }else{
        delete tDbrIdx;
    }



    return EXIT_SUCCESS;
}

int result2profile(DBReader<unsigned int> &resultReader, Parameters &par,
                   const size_t dbFrom, const size_t dbSize) {
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom + dbSize << "\n";

    const std::string outname = par.db4;
    const std::string outnameIndex = par.db4Index;

    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outname, outnameIndex, MMseqsMPI::rank);
    int status = result2profile(resultReader, par, tmpOutput.first, dbFrom, dbSize);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // master reduces results
    if (MMseqsMPI::isMaster()) {
        const char *ext[2] = {"", "_consensus"};
        for (size_t i = 0; i < 2; i++) {
            std::vector<std::pair<std::string, std::string> > splitFiles;
            for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
                std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outname, outnameIndex, procs);
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
    // no pseudo counts
    par.pca = 0.0;
    par.parseParameters(argc, argv, command, 4, false);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(resultReader.getDataSize(), resultReader.getSeqLens(), resultReader.getSize(),
                                     MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);

    int status = result2profile(resultReader, par, dbFrom, dbSize);
#else
    int status = result2profile(resultReader, par, par.db4, 0, resultReader.getSize());
#endif

    if (MMseqsMPI::isMaster()) {
        FileUtil::symlinkAbs(par.hdr1, par.hdr4);
        FileUtil::symlinkAbs(par.hdr1Index, par.hdr4Index);

        if (par.omitConsensus == false) {
            FileUtil::symlinkAbs(par.hdr1, par.db4 + "_consensus_h");
            FileUtil::symlinkAbs(par.hdr1Index, par.db4 + "_consensus_h.index");
        }
    }

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    resultReader.close();

    return status;
}
