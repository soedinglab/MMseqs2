// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

#include "Alignment.h"
#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBConcat.h"
#include "HeaderSummarizer.h"
#include "CompressedA3M.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

enum outputmode {
    MSA = 0,
    PSSM,
    CA3M,
    REPSEQ
};

size_t findMaxSetSize(DBReader<unsigned int> *reader) {
    //Find the max set size
    size_t maxSetSize = 0;
    for (size_t i = 0; i < reader->getSize(); i++) {
        char *data = reader->getData(i);
        size_t entry = 0;
        size_t position = 0;
        while (data[position] != '\0') {
            if (data[position] == '\n') {
                entry++;
            }
            position++;
        }
        maxSetSize = std::max(maxSetSize, entry);
    }
    return maxSetSize + 1;
}

MultipleAlignment::MSAResult computeAlignment(MultipleAlignment &aligner, Sequence *centerSequence,
                                              std::vector<Sequence *> seqSet,
                                              std::vector<Matcher::result_t> alnResults,
                                              bool allowDeletion, bool sameDatabase) {

    if (alnResults.size() > 0) {
        std::vector<Matcher::result_t> alnWithoutIdentity;
        for (size_t i = 0; i < alnResults.size(); i++) {
            if (alnResults[i].dbKey == centerSequence->getDbKey() && sameDatabase == true) { ;
            } else {
                alnWithoutIdentity.push_back(alnResults[i]);
            }
        }

        return aligner.computeMSA(centerSequence, seqSet, alnWithoutIdentity, !allowDeletion);
    } else {
        return aligner.computeMSA(centerSequence, seqSet, !allowDeletion);
    }
}


int result2outputmode(Parameters &par,const std::string &outpath,
                      const size_t dbFrom, const size_t dbSize,
                      const int mode, DBConcat* referenceDBr = NULL) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    if (mode == CA3M && referenceDBr == NULL) {
        Debug(Debug::ERROR) << "Need a sequence and header database for ca3m output!\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    qDbr->open(DBReader<unsigned int>::NOSORT);

    std::string headerNameQuery(par.db1);
    headerNameQuery.append("_h");

    std::string headerIndexNameQuery(par.db1);
    headerIndexNameQuery.append("_h.index");

    std::string headerNameTarget(par.db2);
    headerNameTarget.append("_h");

    std::string headerIndexNameTarget(par.db2);
    headerIndexNameTarget.append("_h.index");

    bool firstSeqRepr = par.db1.compare(par.db2) == 0; 
    if (firstSeqRepr)
        Debug(Debug::INFO) << "Using the first target sequence as center sequence for making each alignment.\n";
    
    DBReader<unsigned int> *queryHeaderReader = new DBReader<unsigned int>(headerNameQuery.c_str(),
                                                                           headerIndexNameQuery.c_str());
    // NOSORT because the index should be in the same order as resultReader
    queryHeaderReader->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = qDbr;
    DBReader<unsigned int> *tempateHeaderReader = queryHeaderReader;

    unsigned int maxSequenceLength = 0;
    unsigned int *lengths = qDbr->getSeqLens();
    for (size_t i = 0; i < qDbr->getSize(); i++) {
        maxSequenceLength = std::max(lengths[i], maxSequenceLength);
    }

    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tDbr->open(DBReader<unsigned int>::SORT_BY_LINE);

        lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max(lengths[i], maxSequenceLength);
        }

        tempateHeaderReader = new DBReader<unsigned int>(headerNameTarget.c_str(), headerIndexNameTarget.c_str());
        tempateHeaderReader->open(DBReader<unsigned int>::NOSORT);
    }

    std::string referenceName(outpath);
    std::string referenceIndexName(outpath);
    referenceIndexName.append(".index");

    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);


//    FileUtil::errorIfFileExist(par.db4.c_str());
//    FileUtil::errorIfFileExist(par.db4Index.c_str());

    DBWriter resultWriter(referenceName.c_str(), referenceIndexName.c_str(), par.threads, DBWriter::BINARY_MODE);
    resultWriter.open();

    DBWriter *concensusWriter = NULL;
    if (mode == PSSM) {
        concensusWriter = new DBWriter(std::string(outpath + "_consensus").c_str(),
                                       std::string(outpath + "_consensus.index").c_str(),
                                       par.threads, DBWriter::ASCII_MODE);
        concensusWriter->open();
    }

    size_t maxSetSize = findMaxSetSize(resultReader);
    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);

    Debug(Debug::INFO) << "Start computing ";
    switch(mode) {
        case MSA:
            Debug(Debug::INFO) << "multiple sequence alignments";
            break;
        case CA3M:
            Debug(Debug::INFO) << "compressed multiple sequence alignments";
            break;
        case REPSEQ:
            Debug(Debug::INFO) << "representative sequences";
            break;
        case PSSM:
        default:
            Debug(Debug::INFO) << "profiles";
            break;
    }
    Debug(Debug::INFO) << ".\n";

#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &subMat, tDbr->getAminoAcidDBSize(), tDbr->getSize(),
                        par.compBiasCorrection);

        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat);
        UniprotHeaderSummarizer summarizer;

        int sequenceType = Sequence::AMINO_ACIDS;
        if (par.profile == true) {
            sequenceType = Sequence::HMM_PROFILE;
        }

        Sequence *centerSequence = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa,
                                                sequenceType, 0, false, par.compBiasCorrection);

#pragma omp for schedule(static)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char *results = resultReader->getData(id);
            char dbKey[255 + 1];
            unsigned int centerSequenceKey;
            char *centerSequenceHeader;

            // check if already aligned results exists
            std::vector<Matcher::result_t> alnResults;

            // Get the sequence from the queryDB
            unsigned int queryKey = resultReader->getDbKey(id);
            char *seqData = qDbr->getDataByDBKey(queryKey);
            
            if (seqData != NULL && !firstSeqRepr)
            {
                centerSequence->mapSequence(0, queryKey, seqData);
                centerSequenceKey = queryKey;
                centerSequenceHeader = queryHeaderReader->getDataByDBKey(centerSequenceKey);
            }

            std::vector<Sequence *> seqSet;
            std::string *reprSeq = NULL;
            
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                double evalue = 0.0;
                char *entry[255];
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // its an aln result
                if (columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                    evalue = strtod(entry[3], NULL);
                }

                if(reprSeq == NULL)
                {
                    const size_t edgeId = tDbr->getId(key);
                    char *dbSeqData = tDbr->getData(edgeId);
                    if (firstSeqRepr)
                    {
                        centerSequence->mapSequence(0, key, dbSeqData);
                        centerSequenceKey = key;
                        centerSequenceHeader = tempateHeaderReader->getDataByDBKey(centerSequenceKey);
                    }
                    reprSeq = new std::string(dbSeqData);
                }


                // just add sequences if eval < thr. and if key is not the same as the query in case of sameDatabase
                if (evalue <= par.evalProfile && (key != queryKey || sameDatabase == false)) {
                    if (columns > Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                        alnResults.push_back(res);
                    }
                    const size_t edgeId = tDbr->getId(key);
                    char *dbSeqData = tDbr->getData(edgeId);
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLens(edgeId), subMat.aa2int, subMat.int2aa,
                                                          Sequence::AMINO_ACIDS, 0, false, false);
                    edgeSequence->mapSequence(0, key, dbSeqData);
                    seqSet.push_back(edgeSequence);
                }
                results = Util::skipLine(results);
            }

            // Recompute the backtrace if the center seq has to be the first seq
            // or if not all the backtraces are present
            if (firstSeqRepr || alnResults.size() != seqSet.size()) 
            {
                //Debug(Debug::INFO) << "BT info is missing... will recompute alignments.\n";
                alnResults.clear(); // will force to recompute alignments together with bt information
            }


            MultipleAlignment::MSAResult res = computeAlignment(aligner, centerSequence, seqSet,
                                                                alnResults, par.allowDeletion, sameDatabase);
            //MultipleAlignment::print(res, &subMat);

            alnResults = res.alignmentResults;

            MsaFilter::MsaFilterResult filterRes = filter.filter((const char **) res.msaSequence, res.setSize,
                                                                 res.centerLength, static_cast<int>(par.cov * 100),
                                                                 static_cast<int>(par.qid * 100), par.qsc,
                                                                 static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);
                                                                 
                                                                 
            char *data;
            size_t dataSize;
            std::ostringstream msa;
            std::string result;
            switch (mode) {
                case MSA: {
                    if (par.summarizeHeader) {
                        // gather headers for summary
                        std::vector<std::string> headers;
                        for (size_t i = 0; i < filterRes.setSize; i++) {
                            if (i == 0) {
                                headers.push_back(centerSequenceHeader);
                            } else {
                                headers.push_back(tempateHeaderReader->getDataByDBKey(seqSet[i]->getDbKey()));
                            }
                        }

                        std::string summary = summarizer.summarize(headers);
                        msa << "#" << par.summaryPrefix.c_str() << "-" << centerSequenceKey << "|" << summary.c_str() << "\n";
                    }

                    // TODO : the first sequence in the MSA seems to be overwritten by the query seq
                    for (size_t i = 0; i < filterRes.setSize; i++) {
                        unsigned int key;
                        char* header;
                        if (i == 0) {
                            key = centerSequenceKey;
                            header = centerSequenceHeader;
                        } else {
                            key = seqSet[i-1]->getDbKey();
                            header = tempateHeaderReader->getDataByDBKey(key);
                        }
                        if (par.addInternalId) {
                            msa << "#" << key << "\n";
                        }

                        msa << ">" << header;

                        // need to allow insertion in the centerSequence
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            char aa = filterRes.filteredMsaSequence[i][pos];
                            msa << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
                        }

                        msa << "\n";
                    }
                    result = msa.str();

                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;
                case REPSEQ: {
                    if (reprSeq != NULL) {
                        msa << *reprSeq;
                        delete reprSeq;
                    } else {
                        //msa << '\0';
                    }
                    result = msa.str();
                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;

                case CA3M: {
                    // Here the backtrace information should be present in the alnResults[i].backtrace for all i
                    std::string consensusStr;
                    std::vector<Matcher::result_t> filteredAln; // alignment information for the sequences that passed the filtering step

                    
                    if (par.useConsensus)
                    {
                        std::pair<const char*, std::string> pssmRes = calculator.computePSSMFromMSA(filterRes.setSize, res.centerLength,
                                                                                                    filterRes.filteredMsaSequence, par.wg);
                        consensusStr = pssmRes.second;
                    } else { // use query sequence as consensus sequence
                        std::ostringstream centerSeqStr;
                        // Retrieve the master sequence
                        for (int pos = 0; pos < centerSequence->L; pos++) {
                            centerSeqStr << subMat.int2aa[centerSequence->int_sequence[pos]];
                        }
                        consensusStr = centerSeqStr.str();

                    }

                    // Put the query sequence (master sequence) first in the alignment
                    Matcher::result_t firstSequence;
                    firstSequence.dbKey = centerSequenceKey;
                    firstSequence.qStartPos = 0;
                    firstSequence.dbStartPos = 0;
                    firstSequence.backtrace = std::string(centerSequence->L, 'M'); // only matches
                    filteredAln.push_back(firstSequence);

                    // Filtering out alnResults using the same filter as previously
                    // defined in object filterRes (with sequence coverage, etc.).
                    for (size_t i = 0; i < alnResults.size(); i++) {
                        // +1 : index 0 corresponds to the query sequence
                        if (filterRes.keep[i + 1]) {
                            filteredAln.push_back(alnResults.at(i));
                        }
                    }

                    // Write the consensus sequence
                    msa << ">consensus_" << queryHeaderReader->getDataByDBKey(queryKey) << consensusStr.c_str() << "\n;";

                    msa << CompressedA3M::fromAlignmentResult(filteredAln, *referenceDBr);

                    result = msa.str();
                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;
                case PSSM: {
                    std::pair<const char*, std::string> pssmRes = calculator.computePSSMFromMSA(filterRes.setSize, res.centerLength,
                                                                                                filterRes.filteredMsaSequence, par.wg);
                    data = (char*)pssmRes.first;
                    std::string consensusStr = pssmRes.second;
                    dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

                    consensusStr.push_back('\n');
                    concensusWriter->writeData(consensusStr.c_str(), consensusStr.length(), SSTR(queryKey).c_str(),
                                               thread_idx);

                    for (size_t i = 0; i < dataSize; i++) {
                        // Avoid a null byte result
                        data[i] = data[i] ^ 0x80;
                    }

                    //pssm.printProfile(res.centerLength);
                    //calculator.printPSSM(res.centerLength);
                }
                    break;
                default:
                    EXIT(EXIT_FAILURE);
            }

            resultWriter.writeData(data, dataSize, SSTR(queryKey).c_str(), thread_idx);

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }
        delete centerSequence;
    }

    // cleanup
    resultWriter.close();
    if (mode == PSSM) {
        concensusWriter->close();
        delete concensusWriter;
    }

    resultReader->close();
    delete resultReader;

    if (!sameDatabase) {
        tempateHeaderReader->close();
        tDbr->close();
        delete tempateHeaderReader;
        delete tDbr;
    }

    queryHeaderReader->close();
    qDbr->close();
    delete queryHeaderReader;
    delete qDbr;

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int result2outputmode(Parameters &par, int mode) {
  

    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    resultReader->open(DBReader<unsigned int>::NOSORT);
    size_t resultSize = resultReader->getSize();
    resultReader->close();
    delete resultReader;

    std::string outname = par.db4;
    DBConcat *referenceDBr = NULL;
    if (mode == CA3M) {
        std::string referenceSeqName(outname);
        std::string referenceSeqIndexName(outname);
        referenceSeqName.append("_sequence.ffdata");
        referenceSeqIndexName.append("_sequence.ffindex");

        // Use only 1 thread for concat to ensure the same order as the later header concat
        referenceDBr = new DBConcat(par.db1.c_str(), par.db1Index.c_str(), par.db2.c_str(), par.db2Index.c_str(),
                                    referenceSeqName.c_str(), referenceSeqIndexName.c_str(), 1);
        referenceDBr->concat();
        // When exporting in ca3m,
        // we need to have an access in SORT_BY_LINE
        // mode in order to keep track of the original
        // line number in the index file.
        referenceDBr->open(DBReader<unsigned int>::SORT_BY_LINE);

        std::string headerQuery(par.db1);
        headerQuery.append("_h");
        std::pair<std::string, std::string> query = Util::databaseNames(headerQuery);

        std::string headerTarget(par.db2);
        headerTarget.append("_h");
        std::pair<std::string, std::string> target = Util::databaseNames(headerTarget);

        std::string referenceHeadersName(outname);
        std::string referenceHeadersIndexName(outname);

        referenceHeadersName.append("_header.ffdata");
        referenceHeadersIndexName.append("_header.ffindex");

        // Use only 1 thread for concat to ensure the same order as the former sequence concat
        DBConcat referenceHeadersDBr(query.first.c_str(), query.second.c_str(),
                                     target.first.c_str(), target.second.c_str(),
                                     referenceHeadersName.c_str(), referenceHeadersIndexName.c_str(), 1);
        referenceHeadersDBr.concat();

        outname.append("_ca3m");
    }

    int status = result2outputmode(par, outname, 0, resultSize, mode, referenceDBr);

    
    if(referenceDBr != NULL) {
        referenceDBr->close();
        delete referenceDBr;
    }
    return status;
}

int result2outputmode(Parameters &par, int mode, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    qDbr->open(DBReader<unsigned int>::NOSORT);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(qDbr->getAminoAcidDBSize(), qDbr->getSeqLens(), qDbr->getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    qDbr->close();
    delete qDbr;

     
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom+dbSize << "\n";

    DBConcat *referenceDBr = NULL;
    std::string outname = par.db4;
    if (mode == CA3M) {
        std::string referenceSeqName(outname);
        std::string referenceSeqIndexName(outname);
        referenceSeqName.append("_sequence.ffdata");
        referenceSeqIndexName.append("_sequence.ffindex");

        // Use only 1 thread for concat to ensure the same order as the later header concat
        referenceDBr = new DBConcat(par.db1.c_str(), par.db1Index.c_str(), par.db2.c_str(), par.db2Index.c_str(),
                                    referenceSeqName.c_str(), referenceSeqIndexName.c_str(), 1);
        if(mpiRank == 0) {
            referenceDBr->concat();
        } else {
            referenceDBr->concat(false);
        }

#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // When exporting in ca3m,
        // we need to have an access in SORT_BY_LINE
        // mode in order to keep track of the original
        // line number in the index file.
        referenceDBr->open(DBReader<unsigned int>::SORT_BY_LINE);

        if(mpiRank == 0) {
            std::string headerQuery(par.db1);
            headerQuery.append("_h");
            std::pair<std::string, std::string> query = Util::databaseNames(headerQuery);

            std::string headerTarget(par.db2);
            headerTarget.append("_h");
            std::pair<std::string, std::string> target = Util::databaseNames(headerTarget);

            std::string referenceHeadersName(outname);
            std::string referenceHeadersIndexName(outname);

            referenceHeadersName.append("_header.ffdata");
            referenceHeadersIndexName.append("_header.ffindex");

            // Use only 1 thread for concat to ensure the same order as the former sequence concat
            DBConcat referenceHeadersDBr(query.first.c_str(), query.second.c_str(),
                                         target.first.c_str(), target.second.c_str(),
                                         referenceHeadersName.c_str(), referenceHeadersIndexName.c_str(), 1);
            referenceHeadersDBr.concat();
        }

        outname.append("_ca3m");
    }

    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outname, "", mpiRank);
    int status = result2outputmode(par, tmpOutput.first, dbFrom, dbSize, mode, referenceDBr);
    

    // close reader to reduce memory
    if(referenceDBr != NULL){
        referenceDBr->close();
        delete referenceDBr;
    }

    
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // master reduces results
    if(mpiRank == 0) {
        const char * ext[2] = {"", "_consensus"};
        size_t extCount = 1;
        if(mode == PSSM){
            extCount = 2;
        }
        for(size_t i = 0; i < extCount; i++){
            std::vector<std::pair<std::string, std::string> > splitFiles;
            for(unsigned int procs = 0; procs < mpiNumProc; procs++){
                std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outname, "", procs);
                splitFiles.push_back(std::make_pair(tmpFile.first + ext[i],  tmpFile.first + ext[i] + ".index"));

            }
            // merge output ffindex databases
            Alignment::mergeAndRemoveTmpDatabases(outname + ext[i], outname + ext[i] + ".index", splitFiles);
        }
    }

    return status;
}

int result2profile(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::WARNING) << "Compute profile.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

#ifdef HAVE_MPI
    int retCode = result2outputmode(par, PSSM, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int retCode = result2outputmode(par, PSSM);
#endif

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return retCode;
}

int result2msa(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    MMseqsMPI::init(argc, argv);

    struct timeval start, end;
    gettimeofday(&start, NULL);
    int retCode;

    outputmode mode;
    if (par.compressMSA) {
        mode = CA3M;
    } else if(par.onlyRepSeq) {
        mode = REPSEQ;
    } else {
        mode = MSA;
    }


    /*
    // Can only use the first sequence as representative sequence with a clustering result
    // otherwise the key would point to a sequence in the query DB instead of the target DB
    // TODO : can fix that by passing the firstSeqRepr to CompressedA3M::fromAlignmentResult
    if (par.firstSeqRepr && par.db1.compare(par.db2) != 0)
    {
        Debug(Debug::ERROR) << "Use first sequence as representative only with clustering results (<queryDB> and <targetDB> should be the same).\n";
        return -1;
    }*/
    
#ifdef HAVE_MPI
    retCode = result2outputmode(par, mode, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    retCode = result2outputmode(par, mode);
#endif

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return retCode;
}
