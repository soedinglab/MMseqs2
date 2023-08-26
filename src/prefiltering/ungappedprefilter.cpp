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
#include "FastSort.h"
#include "SubstitutionMatrixProfileStates.h"
#include "IndexReader.h"
#include "QueryMatcherTaxonomyHook.h"
#ifdef OPENMP
#include <omp.h>
#endif

int prefilterInternal(int argc, const char **argv, const Command &command, int mode) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();
    bool sameDB = (par.db2.compare(par.db1) == 0);
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader tDbrIdx(par.db2, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    IndexReader * qDbrIdx = NULL;
    DBReader<unsigned int> * qdbr = NULL;
    DBReader<unsigned int> * tdbr = tDbrIdx.sequenceReader;
    const int targetSeqType = tdbr->getDbtype();
    int querySeqType;
    if (sameDB == true) {
        qDbrIdx = &tDbrIdx;
        qdbr = tdbr;
        querySeqType = targetSeqType;
    } else {
        // open the sequence, prefiltering and output databases
        qDbrIdx = new IndexReader(par.db1, par.threads,  IndexReader::SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0);
        qdbr = qDbrIdx->sequenceReader;
        querySeqType = qdbr->getDbtype();
    }

    SequenceLookup * sequenceLookup = NULL;
    if(Parameters::isEqualDbtype(tDbrIdx.getDbtype(), Parameters::DBTYPE_INDEX_DB)){
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(tDbrIdx.index);
        if(data.splits == 1){
            sequenceLookup = PrefilteringIndexReader::getSequenceLookup(0, tDbrIdx.index, par.preloadMode);
        }
    }
    BaseMatrix *subMat;
    EvalueComputation * evaluer;
    int8_t * tinySubMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat);
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        evaluer = new EvalueComputation(tdbr->getAminoAcidDBSize(), subMat);
        tinySubMat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
        for (int i = 0; i < subMat->alphabetSize; i++) {
            for (int j = 0; j < subMat->alphabetSize; j++) {
                tinySubMat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
            }
        }
    }


    QueryMatcherTaxonomyHook * taxonomyHook = NULL;
    if(par.PARAM_TAXON_LIST.wasSet){
        taxonomyHook = new QueryMatcherTaxonomyHook(par.db2, tdbr, par.taxonList);
    }

    Debug::Progress progress(qdbr->getSize());
    std::vector<hit_t> shortResults;
    shortResults.reserve(tdbr->getSize()/2);

#ifdef OPENMP
    omp_set_nested(1);
#endif

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        char buffer[1024+32768];
        std::vector<hit_t> threadShortResults;
        Sequence qSeq(par.maxSeqLen, querySeqType, subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, targetSeqType, subMat, 0, false, par.compBiasCorrection);
        SmithWaterman aligner(par.maxSeqLen, subMat->alphabetSize,
                              par.compBiasCorrection, par.compBiasCorrectionScale, targetSeqType);

        std::string resultBuffer;
        resultBuffer.reserve(262144);
        for (size_t id = 0; id < qdbr->getSize(); id++) {
            char *querySeqData = qdbr->getData(id, thread_idx);
            size_t queryKey = qdbr->getDbKey(id);
            unsigned int querySeqLen = qdbr->getSeqLen(id);

            qSeq.mapSequence(id, queryKey, querySeqData, querySeqLen);
//            qSeq.printProfileStatePSSM();
            if(Parameters::isEqualDbtype(qSeq.getSeqType(), Parameters::DBTYPE_HMM_PROFILE) ){
                aligner.ssw_init(&qSeq, qSeq.getAlignmentProfile(), subMat);
            }else{
                aligner.ssw_init(&qSeq, tinySubMat, subMat);
            }
#pragma omp for schedule(static) nowait
            for (size_t tId = 0; tId < tdbr->getSize(); tId++) {
                unsigned int targetKey = tdbr->getDbKey(tId);
                if(taxonomyHook != NULL){
                    TaxID currTax = taxonomyHook->taxonomyMapping->lookup(targetKey);
                    if (taxonomyHook->expression->isAncestor(currTax) == false) {
                        continue;
                    }
                }

                const bool isIdentity = (queryKey == targetKey && (par.includeIdentity || sameDB))? true : false;
                if(sequenceLookup == NULL){
                    char * targetSeq = tdbr->getData(tId, thread_idx);
                    unsigned int targetSeqLen = tdbr->getSeqLen(tId);
                    tSeq.mapSequence(tId, targetKey, targetSeq, targetSeqLen);
                }else{
                    tSeq.mapSequence(tId, targetKey, sequenceLookup->getSequence(tId));
                }
                float queryLength = qSeq.L;
                float targetLength = tSeq.L;
                if(Util::canBeCovered(par.covThr, par.covMode, queryLength, targetLength)==false){
                    continue;
                }

                int score;
                if (mode == 0) {
                    score = aligner.ungapped_alignment(tSeq.numSequence, tSeq.L);
                } else {
                    std::string backtrace;
                    s_align res;
                    if (isIdentity) {
                        res = aligner.scoreIdentical(
                            tSeq.numSequence, tSeq.L, evaluer, Matcher::SCORE_ONLY, backtrace
                        );
                    } else {
                        res = aligner.ssw_align(
                            tSeq.numSequence,
                            tSeq.numConsensusSequence,
                            tSeq.getAlignmentProfile(),
                            tSeq.L,
                            backtrace,
                            par.gapOpen.values.aminoacid(),
                            par.gapExtend.values.aminoacid(),
                            Matcher::SCORE_ONLY,
                            par.evalThr,
                            evaluer,
                            par.covMode,
                            par.covThr,
                            par.correlationScoreWeight,
                            qSeq.L / 2,
                            tId
                        );
                    }
                    score = res.score1;
                }
                bool hasDiagScore = (score > par.minDiagScoreThr);
                double evalue = 0.0;
                // check if evalThr != inf
                if (par.evalThr < std::numeric_limits<double>::max()) {
                    evalue = evaluer->computeEvalue(score, qSeq.L);
                }
                bool hasEvalue = (evalue <= par.evalThr);
                // --filter-hits
                if (isIdentity || (hasDiagScore && hasEvalue)) {
                    hit_t hit;
                    hit.seqId = targetKey;
                    hit.prefScore = score;
                    hit.diagonal = 0;
                    threadShortResults.emplace_back(hit);
                }
            }
#pragma omp critical
            {
                shortResults.insert(shortResults.end(), threadShortResults.begin(), threadShortResults.end());
                threadShortResults.clear();
            }
#pragma omp barrier
#pragma omp master
            {
                SORT_PARALLEL(shortResults.begin(), shortResults.end(), hit_t::compareHitsByScoreAndId);
                size_t maxSeqs = std::min(par.maxResListLen, shortResults.size());
                for (size_t i = 0; i < maxSeqs; ++i) {
                    size_t len = QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
                    resultBuffer.append(buffer, len);
                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, 0);
                resultBuffer.clear();
                shortResults.clear();
                progress.updateProgress();
            }
#pragma omp barrier
        }
    }


    if(taxonomyHook != NULL){
        delete taxonomyHook;
    }

    if(sameDB == false){
        delete qDbrIdx;
    }

    delete [] tinySubMat;
    delete subMat;
    delete evaluer;

    resultWriter.close();
    return 0;
}

int ungappedprefilter(int argc, const char **argv, const Command &command) {
    return prefilterInternal(argc, argv, command, 0);
}

int gappedprefilter(int argc, const char **argv, const Command &command) {
    return prefilterInternal(argc, argv, command, 1);
}
