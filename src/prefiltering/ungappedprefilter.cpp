//
// Created by Martin Steinegger on 17.09.18.
//

#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "QueryMatcher.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrixProfileStates.h"

#ifdef OPENMP
#include <omp.h>
#endif


int doRescorealldiagonal(Parameters &par, DBReader<unsigned int> &qdbr, DBWriter &resultWriter, size_t dbStart, size_t dbSize) {
    int querySeqType = qdbr.getDbtype();
    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    DBReader<unsigned int> *tdbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.noPreload == false) {
            tdbr->readMmapedDataInMemory();
        }
    }
    const int targetSeqType = tdbr->getDbtype();
    if (querySeqType == Sequence::HMM_PROFILE && targetSeqType == Sequence::PROFILE_STATE_SEQ) {
        querySeqType = Sequence::PROFILE_STATE_PROFILE;
    }

    BaseMatrix *subMat;
    EvalueComputation * evaluer;
    int8_t * tinySubMat;
    if (querySeqType == Sequence::NUCLEOTIDES) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen, par.gapExtend, false);
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    } else if(targetSeqType == Sequence::PROFILE_STATE_SEQ ){
        SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), &sMat, par.gapOpen, par.gapExtend, false);
        subMat = new SubstitutionMatrixProfileStates(sMat.matrixName, sMat.probMatrix, sMat.pBack,
                                                     sMat.subMatrixPseudoCounts, 2.0, 0.0, 219);
        tinySubMat = new int8_t[sMat.alphabetSize*sMat.alphabetSize];
        for (int i = 0; i < sMat.alphabetSize; i++) {
            for (int j = 0; j < sMat.alphabetSize; j++) {
                tinySubMat[i*sMat.alphabetSize + j] = sMat.subMatrix[i][j];
            }
        }
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen, par.gapExtend, false);
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    }



#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        char buffer[1024+32768];
        std::vector<hit_t> shortResults;
        shortResults.reserve(std::max(static_cast<size_t >(1), tdbr->getSize()/5));
        Sequence qSeq(par.maxSeqLen, querySeqType, subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, targetSeqType, subMat, 0, false, par.compBiasCorrection);
        SmithWaterman aligner(par.maxSeqLen, subMat->alphabetSize, par.compBiasCorrection);

        std::string resultBuffer;
        resultBuffer.reserve(262144);
#pragma omp for schedule(dynamic, 1)
        for (size_t id = dbStart; id < (dbStart+dbSize); id++) {
            Debug::printProgress(id);
            char *querySeqData = qdbr.getData(id);
            size_t queryKey = qdbr.getDbKey(id);
            qSeq.mapSequence(id, queryKey, querySeqData);
//            qSeq.printProfileStatePSSM();
            if(qSeq.getSeqType() == Sequence::HMM_PROFILE || qSeq.getSeqType() == Sequence::PROFILE_STATE_PROFILE){
                aligner.ssw_init(&qSeq, qSeq.getAlignmentProfile(), subMat, subMat->alphabetSize, 0);
            }else{
                aligner.ssw_init(&qSeq, tinySubMat, subMat, subMat->alphabetSize, 0);
            }

            for (size_t tId = 0; tId < tdbr->getSize(); tId++) {
                unsigned int targetKey = tdbr->getDbKey(tId);
                const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameDB))? true : false;
                char * targetSeq = tdbr->getData(tId);
                tSeq.mapSequence(tId, targetKey, targetSeq);
//                tSeq.print();
                float queryLength = qSeq.L;
                float targetLength = tSeq.L;
                if(Util::canBeCovered(par.covThr, par.covMode, queryLength, targetLength)==false){
                    continue;
                }

                int score = aligner.ungapped_alignment(tSeq.int_sequence, tSeq.L);
                bool hasDiagScore = (score > par.minDiagScoreThr);
                double evalue = evaluer->computeEvalue(score, qSeq.L);
                bool hasEvalue = (evalue <= par.evalThr);
                // --filter-hits
                if (isIdentity || (hasDiagScore && hasEvalue)) {
                    hit_t hit;
                    hit.seqId = targetKey;
                    hit.pScore = score;
                    hit.prefScore = evalue;
                    hit.diagonal = 0;
                    shortResults.emplace_back(hit);
                }
            }

            std::sort(shortResults.begin(), shortResults.end(), hit_t::compareHitsByEValueAndId);

            for (size_t i = 0; i < shortResults.size(); ++i) {
                size_t len = QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
                resultBuffer.append(buffer, len);
            }

            resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            shortResults.clear();
        }
    }

    Debug(Debug::INFO) << "\nDone.\n";
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    delete [] tinySubMat;
    delete subMat;
    delete evaluer;
    return 0;
}

int ungappedprefilter(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);
    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    if (par.noPreload == false) {
        qdbr.readMmapedDataInMemory();
    }


#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;

    Util::decomposeDomainByAminoAcid(qdbr.getAminoAcidDBSize(), qdbr.getSeqLens(), qdbr.getSize(),
                                     MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db3, par.db3Index, MMseqsMPI::rank);
    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), par.threads);
    resultWriter.open();
    int status = doRescorealldiagonal(par, qdbr, resultWriter, dbFrom, dbSize);
    resultWriter.close();

    MPI_Barrier(MPI_COMM_WORLD);
    if(MMseqsMPI::rank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for(unsigned int proc = 0; proc < MMseqsMPI::numProc; ++proc){
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db3, par.db3Index, proc);
            splitFiles.push_back(std::make_pair(tmpFile.first,  tmpFile.second));
        }
        DBWriter::mergeResults(par.db3, par.db3Index, splitFiles);
    }
#else
    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    resultWriter.open();
    int status = doRescorealldiagonal(par, qdbr, resultWriter, 0, qdbr.getSize());
    resultWriter.close();
#endif
    return status;
}


