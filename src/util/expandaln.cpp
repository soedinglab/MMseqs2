#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "BacktraceTranslator.h"
#include "EvalueComputation.h"
#include "CompressedA3M.h"
#include "Sequence.h"
#include "Alignment.h"
#include "SubstitutionMatrix.h"
#include "FastSort.h"
#include <cassert>

#ifdef OPENMP
#include <omp.h>
#endif

void rescoreResultByBacktrace(Matcher::result_t &result, Sequence &qSeq, Sequence &tSeq, SubstitutionMatrix &subMat, float *compositionBias, EvalueComputation &evaluer, int gapOpen, int gapExtend, int seqIdMode) {
    size_t qPos = result.qStartPos;
    size_t tPos = result.dbStartPos;
    int score = 0;
    char lastState = '\0';
    int identities = 0;
    const bool isQueryProf = Parameters::isEqualDbtype(qSeq.getSeqType(), Parameters::DBTYPE_HMM_PROFILE);
    const bool isTargetProf = Parameters::isEqualDbtype(tSeq.getSeqType(), Parameters::DBTYPE_HMM_PROFILE);
//    for(int i = result.qStartPos; i < result.qEndPos; i++){
//        printf("%c",subMat.num2aa[qSeq.sequence[i]]);
//    }
//    Debug(Debug::INFO) << "\n";
//    for(int i = result.dbStartPos; i < result.dbEndPos; i++){
//        printf("%c",subMat.num2aa[tSeq.sequence[i]]);
//    }
//    Debug(Debug::INFO) << "\n";
    for (size_t i = 0; i < result.backtrace.size(); ++i) {
        char state = result.backtrace[i];
        if (state == 'M') {
            if (isTargetProf) {
                score += tSeq.profile_for_alignment[qSeq.numSequence[qPos] * tSeq.L + tPos]  + static_cast<short>((compositionBias[i] < 0.0)? compositionBias[i] - 0.5: compositionBias[i] + 0.5);;
            } else if (isQueryProf) {
                score += qSeq.profile_for_alignment[tSeq.numSequence[tPos] * qSeq.L + qPos];
            } else {
                score += subMat.subMatrix[qSeq.numSequence[qPos]][tSeq.numSequence[tPos]] + static_cast<short>((compositionBias[i] < 0.0)? compositionBias[i] - 0.5: compositionBias[i] + 0.5);
            }
            identities += qSeq.numSequence[qPos] == tSeq.numSequence[tPos] ? 1 : 0;
            qPos++;
            tPos++;
        } else if (state == 'I') {
            if (lastState == 'I') {
                // TODO no std::max(0, gapExtend)?
                score -= gapExtend;
            } else {
                score -= gapOpen;
            }
            tPos++;
        } else if (state == 'D') {
            if (lastState == 'D') {
                score -= gapExtend;
            } else {
                score -= gapOpen;
            }
            qPos++;
        }
        lastState = state;
    }
    result.eval = evaluer.computeEvalue(score, qSeq.L);
    result.score = static_cast<int>(evaluer.computeBitScore(score)+0.5);
    result.seqId = Util::computeSeqId(seqIdMode, identities, qSeq.L, tSeq.L, result.backtrace.size());
}

static bool compareHitsByKeyEvalScore(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.dbKey < second.dbKey)
        return true;
    if (second.dbKey < first.dbKey)
        return false;
    if (first.eval < second.eval)
        return true;
    if (second.eval < first.eval)
        return false;
    if (first.score > second.score)
        return true;
    if (second.score > first.score)
        return false;

    return false;
}

int expandaln(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    queryReader.open(DBReader<unsigned int>::NOSORT);
    const int queryDbType = queryReader.getDbtype();
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        queryReader.readMmapedDataInMemory();
    }

    DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    targetReader.open(DBReader<unsigned int>::NOSORT);
    const int targetDbType = targetReader.getDbtype();
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        targetReader.readMmapedDataInMemory();
    }

    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-profile is currently not supported.\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> *resultReader = NULL;
    DBReader<unsigned int> *ca3mSequenceReader = NULL;
    if (FileUtil::fileExists((par.db3 + "_ca3m.ffdata").c_str())) {

        resultReader = new DBReader<unsigned int>((par.db3 + "_ca3m.ffdata").c_str(), (par.db3 + "_ca3m.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        ca3mSequenceReader = new DBReader<unsigned int>((par.db3 + "_sequence.ffdata").c_str(), (par.db3 + "_sequence.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        ca3mSequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    } else {

        resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }

    DBReader<unsigned int> expansionReader(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    expansionReader.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        expansionReader.readMmapedDataInMemory();
    }

    Debug(Debug::INFO) << "Output database: " << par.db5 << "\n";
    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    BacktraceTranslator translator;
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);
    EvalueComputation evaluer(targetReader.getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    Debug::Progress progress(resultReader->getSize());

    Debug(Debug::INFO) << "Computing expanded alignment result...\n";
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence qSeq(par.maxSeqLen, queryDbType, &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, targetDbType, &subMat, 0, false, false);
        float *compositionBias = new float[par.maxSeqLen + 1]();

        std::vector<Matcher::result_t> expanded;
        expanded.reserve(300);

        std::vector<Matcher::result_t> results;
        results.reserve(1000);

        char buffer[1024];

        Matcher::result_t resultAC;
        resultAC.backtrace.reserve(par.maxSeqLen + 1);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader->getSize(); ++i) {
            progress.updateProgress();
            unsigned int queryKey = resultReader->getDbKey(i);

            size_t querySeqId = queryReader.getId(queryKey);
            qSeq.mapSequence(querySeqId, queryKey, queryReader.getData(querySeqId, thread_idx),
                             queryReader.getSeqLen(querySeqId));

            if(par.compBiasCorrection == true && Parameters::isEqualDbtype(queryDbType,Parameters::DBTYPE_AMINO_ACIDS)){
                SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, qSeq.numSequence, qSeq.L, compositionBias);
            }

            char *data = resultReader->getData(i, thread_idx);
            while (*data != '\0') {
                Matcher::result_t resultAB = Matcher::parseAlignmentRecord(data, false);

                if (resultAB.backtrace.empty()) {
                    Debug(Debug::ERROR) << "Alignment must contain a backtrace.\n";
                    EXIT(EXIT_FAILURE);
                }
//                Matcher::resultToBuffer(buffer, resultAB, true, true);
//                Debug(Debug::INFO) << buffer;

                unsigned int targetKey = resultAB.dbKey;
                size_t targetId = expansionReader.getId(targetKey);

                if (ca3mSequenceReader != NULL) {
                    unsigned int key;
                    CompressedA3M::extractMatcherResults(key, expanded, expansionReader.getData(targetId, thread_idx),
                                                         expansionReader.getEntryLen(targetId), *ca3mSequenceReader, false);
                } else {
                    Matcher::readAlignmentResults(expanded, expansionReader.getData(targetId, thread_idx), false);
                }
                for (size_t k = 0; k < expanded.size(); ++k) {
                    Matcher::result_t &resultBC = expanded[k];

                    if (resultBC.backtrace.empty()) {
                        Debug(Debug::ERROR) << "Alignment must contain a backtrace.\n";
                        EXIT(EXIT_FAILURE);
                    }
//                    Matcher::resultToBuffer(buffer, resultBC, true, true);
//                    Debug(Debug::INFO) << buffer;

                    translator.translateResult(resultAB, resultBC, resultAC);
                    if (resultAC.backtrace.empty()) {
                        continue;
                    }

                    if (Util::canBeCovered(par.covThr, par.covMode, resultAC.qLen, resultAC.dbLen) == false) {
                        continue;
                    }
                    size_t bcTargetSeqId = targetReader.getId(resultBC.dbKey);
                    tSeq.mapSequence(bcTargetSeqId, targetKey, targetReader.getData(bcTargetSeqId, thread_idx),
                                     targetReader.getSeqLen(bcTargetSeqId));

                    rescoreResultByBacktrace(resultAC, qSeq, tSeq, subMat, compositionBias,
                                             evaluer, par.gapOpen.aminoacids, par.gapExtend.aminoacids, par.seqIdMode);

                    if (Alignment::checkCriteria(resultAC, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        results.emplace_back(resultAC);
                    }
                }
                expanded.clear();
                data = Util::skipLine(data);
            }

            std::vector<Matcher::result_t> *finalResults = &results;
            if (par.expansionMode == 1) {
                // keep only the best hit to same target
                SORT_SERIAL(results.begin(), results.end(), compareHitsByKeyEvalScore);
                ssize_t lastKey = -1;
                for (size_t j = 0; j < results.size(); ++j) {
                    const Matcher::result_t& res = results[j];
                    if (res.dbKey != lastKey) {
                        expanded.emplace_back(res);
                    }
                    lastKey = res.dbKey;
                }
                finalResults = &expanded;
            }

            SORT_SERIAL(finalResults->begin(), finalResults->end(), Matcher::compareHits);

            writer.writeStart(thread_idx);
            for (size_t j = 0; j < finalResults->size(); ++j) {
                size_t len = Matcher::resultToBuffer(buffer, (*finalResults)[j], true, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(queryKey, thread_idx);
            expanded.clear();
            results.clear();
        }

        delete[] compositionBias;
    }
    Debug(Debug::INFO) << "\n";
    writer.close();
    expansionReader.close();
    resultReader->close();
    delete resultReader;
    if (ca3mSequenceReader != NULL) {
        ca3mSequenceReader->close();
        delete ca3mSequenceReader;
    }
    targetReader.close();
    queryReader.close();

    return EXIT_SUCCESS;
}
