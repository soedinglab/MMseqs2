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

void rescoreResultByBacktrace(Matcher::result_t &result, Sequence &qSeq, Sequence &tSeq, SubstitutionMatrix &subMat, float *compositionBias, int gapOpen, int gapExtend) {
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
                score += tSeq.profile_for_alignment[qSeq.numSequence[qPos] * tSeq.L + tPos]  + static_cast<short>((compositionBias[i] < 0.0) ? (compositionBias[i] - 0.5) : (compositionBias[i] + 0.5));
            } else if (isQueryProf) {
                score += qSeq.profile_for_alignment[tSeq.numSequence[tPos] * qSeq.L + qPos];
            } else {
                score += subMat.subMatrix[qSeq.numSequence[qPos]][tSeq.numSequence[tPos]] + static_cast<short>((compositionBias[i] < 0.0) ? (compositionBias[i] - 0.5) : (compositionBias[i] + 0.5));
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
            qPos++;
        } else if (state == 'D') {
            if (lastState == 'D') {
                score -= gapExtend;
            } else {
                score -= gapOpen;
            }
            tPos++;
        }
        lastState = state;
    }
    result.score = score;
//    result.eval = evaluer.computeEvalue(score, qSeq.L);
//    result.score = static_cast<int>(evaluer.computeBitScore(score)+0.5);
//    result.seqId = Util::computeSeqId(seqIdMode, identities, qSeq.L, tSeq.L, result.backtrace.size());
}

static bool compareHitsByKeyScore(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.dbKey < second.dbKey)
        return true;
    if (second.dbKey < first.dbKey)
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

    DBReader<unsigned int> aReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    aReader.open(DBReader<unsigned int>::NOSORT);
    const int aSeqDbType = aReader.getDbtype();
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        aReader.readMmapedDataInMemory();
    }

    bool isCa3m = false;
    DBReader<unsigned int> *resultAbReader = NULL;
    DBReader<unsigned int> *cReader = NULL;
    if (FileUtil::fileExists((par.db3 + "_ca3m.ffdata").c_str())) {
        resultAbReader = new DBReader<unsigned int>((par.db3 + "_ca3m.ffdata").c_str(), (par.db3 + "_ca3m.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        resultAbReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        cReader = new DBReader<unsigned int>((par.db3 + "_sequence.ffdata").c_str(), (par.db3 + "_sequence.ffindex").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        cReader->open(DBReader<unsigned int>::SORT_BY_LINE);
        isCa3m = true;
    } else {
        resultAbReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        resultAbReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        cReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        cReader->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            cReader->readMmapedDataInMemory();
        }
    }

    const int cSeqDbType = cReader->getDbtype();
    if (Parameters::isEqualDbtype(aSeqDbType, Parameters::DBTYPE_HMM_PROFILE) && Parameters::isEqualDbtype(cSeqDbType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-profile is currently not supported\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> resultBcReader(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultBcReader.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        resultBcReader.readMmapedDataInMemory();
    }

    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    BacktraceTranslator translator;
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);
    EvalueComputation evaluer(cReader->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    Debug::Progress progress(resultAbReader->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence aSeq(par.maxSeqLen, aSeqDbType, &subMat, 0, false, par.compBiasCorrection);
        Sequence cSeq(par.maxSeqLen, cSeqDbType, &subMat, 0, false, false);

        size_t compBufferSize = (par.maxSeqLen + 1) * sizeof(float);
        float *compositionBias = (float*)malloc(compBufferSize);
        memset(compositionBias, 0, compBufferSize);

        char buffer[1024 + 32000];

        std::vector<Matcher::result_t> resultsBc;
        resultsBc.reserve(300);

        Matcher::result_t resultAc;
        resultAc.backtrace.reserve(par.maxSeqLen + 1);
        Matcher::result_t currBestAc;

        std::vector<Matcher::result_t> resultsAc;
        resultsAc.reserve(1000);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultAbReader->getSize(); ++i) {
            progress.updateProgress();

            unsigned int resAbKey = resultAbReader->getDbKey(i);

            size_t aSeqId = aReader.getId(resAbKey);
            aSeq.mapSequence(aSeqId, resAbKey, aReader.getData(aSeqId, thread_idx), aReader.getSeqLen(aSeqId));

            if (par.compBiasCorrection == true && Parameters::isEqualDbtype(aSeqDbType, Parameters::DBTYPE_AMINO_ACIDS)) {
                if ((size_t)aSeq.L >= compBufferSize) {
                    compBufferSize = (size_t)aSeq.L * 1.5 * sizeof(float);
                    compositionBias = (float*)realloc(compositionBias, compBufferSize);
                    memset(compositionBias, 0, compBufferSize);
                }
                SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, aSeq.numSequence, aSeq.L, compositionBias);
            }

            char *data = resultAbReader->getData(i, thread_idx);
            while (*data != '\0') {
                Matcher::result_t resultAb = Matcher::parseAlignmentRecord(data, false);
                data = Util::skipLine(data);

                if (resultAb.backtrace.empty()) {
                    Debug(Debug::ERROR) << "Alignment must contain a backtrace\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int bResKey = resultAb.dbKey;
                size_t bResId = resultBcReader.getId(bResKey);
                if (isCa3m) {
                    unsigned int key;
                    CompressedA3M::extractMatcherResults(key, resultsBc, resultBcReader.getData(bResId, thread_idx),
                                                         resultBcReader.getEntryLen(bResId), *cReader, false);
                } else {
                    Matcher::readAlignmentResults(resultsBc, resultBcReader.getData(bResId, thread_idx), false);
                }

                std::stable_sort(resultsAc.begin(), resultsAc.end(), compareHitsByKeyScore);

                size_t lastCKey = SIZE_MAX;
                currBestAc.score = INT_MIN;
                for (size_t k = 0; k < resultsBc.size(); ++k) {
                    Matcher::result_t &resultBc = resultsBc[k];
                    if (resultBc.backtrace.size() == 0) {
                        Debug(Debug::ERROR) << "Alignment must contain a backtrace\n";
                        EXIT(EXIT_FAILURE);
                    }
                    translator.translateResult(resultAb, resultBc, resultAc);
                    if (resultAc.backtrace.size() == 0) {
                        continue;
                    }

                    if (Util::canBeCovered(par.covThr, par.covMode, resultAc.qLen, resultAc.dbLen) == false) {
                        continue;
                    }

                    unsigned int cSeqKey = resultBc.dbKey;
                    if (lastCKey != cSeqKey) {
                        if (currBestAc.score != INT_MIN) {
                            resultsAc.emplace_back(currBestAc);
                        }
                        currBestAc.score = INT_MIN;

                        size_t cSeqId = cReader->getId(cSeqKey);
                        cSeq.mapSequence(cSeqId, cSeqKey, cReader->getData(cSeqId, thread_idx), cReader->getSeqLen(cSeqId));
                        lastCKey = cSeqKey;
                    }
                    rescoreResultByBacktrace(resultAc, aSeq, cSeq, subMat, compositionBias, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
                    if (resultAc.score > par.minDiagScoreThr && resultAc.score > currBestAc.score) {
                        currBestAc = resultAc;
                    }
                }
                if (currBestAc.score != INT_MIN) {
                    resultsAc.emplace_back(currBestAc);
                }
                resultsBc.clear();
            }

            SORT_SERIAL(resultsAc.begin(), resultsAc.end(), Matcher::compareHits);
            writer.writeStart(thread_idx);
            for (size_t j = 0; j < resultsAc.size(); ++j) {
                size_t len = Matcher::resultToBuffer(buffer, resultsAc[j], true, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(resAbKey, thread_idx);
            resultsAc.clear();
        }
        free(compositionBias);
    }
    writer.close();
    resultBcReader.close();
    resultAbReader->close();
    delete resultAbReader;
    if (cReader != NULL) {
        cReader->close();
        delete cReader;
    }
    aReader.close();

    return EXIT_SUCCESS;
}
