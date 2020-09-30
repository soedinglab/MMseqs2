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
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
#include "MsaFilter.h"
#include "PSSMCalculator.h"
#include "PSSMMasker.h"
#include "FastSort.h"
#include "IntervalArray.h"
#include <stack>
#include <map>
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
                score += tSeq.profile_for_alignment[qSeq.numSequence[qPos] * tSeq.L + tPos];
            } else if (isQueryProf) {
                score += qSeq.profile_for_alignment[tSeq.numSequence[tPos] * qSeq.L + qPos];
            } else {
                score += subMat.subMatrix[qSeq.numSequence[qPos]][tSeq.numSequence[tPos]] + static_cast<short>((compositionBias[qPos] < 0.0) ? (compositionBias[qPos] - 0.5) : (compositionBias[qPos] + 0.5));
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
    result.seqId = identities;
}

static bool compareHitsByKeyScore(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.score > second.score)
        return true;
    if (second.score > first.score)
        return false;
    return false;
}

int expandaln(int argc, const char **argv, const Command& command, bool returnAlnRes) {
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

    int dbType = Parameters::DBTYPE_ALIGNMENT_RES;
    if (returnAlnRes == false) {
        dbType = Parameters::DBTYPE_HMM_PROFILE;
    }
    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, dbType);
    writer.open();

    BacktraceTranslator translator;
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);\

    EvalueComputation *evaluer = NULL;
    ProbabilityMatrix *probMatrix = NULL;
    if (returnAlnRes == false) {
        probMatrix = new ProbabilityMatrix(subMat);
    } else {
        evaluer = new EvalueComputation(cReader->getAminoAcidDBSize(), &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    }
    Debug::Progress progress(resultAbReader->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence aSeq(par.maxSeqLen, aSeqDbType, &subMat, 0, false, par.compBiasCorrection);
        Sequence cSeq(par.maxSeqLen, cSeqDbType, &subMat, 0, false, false);

        MultipleAlignment *aligner = NULL;
        MsaFilter *filter = NULL;
        PSSMCalculator *calculator = NULL;
        PSSMMasker *masker = NULL;
        std::vector<std::vector<unsigned char>> seqSet;
        std::string result;

        if (returnAlnRes == false) {
            aligner = new MultipleAlignment(par.maxSeqLen, &subMat);
            if (par.filterMsa) {
                filter = new MsaFilter(par.maxSeqLen, 300, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
            }
            calculator = new PSSMCalculator(&subMat, par.maxSeqLen, 300, par.pca, par.pcb);
            masker = new PSSMMasker(par.maxSeqLen, *probMatrix, subMat);
            seqSet.reserve(300);
            result.reserve(par.maxSeqLen * Sequence::PROFILE_READIN_SIZE);
        }

        size_t compBufferSize = (par.maxSeqLen + 1) * sizeof(float);
        float *compositionBias = (float*)malloc(compBufferSize);
        memset(compositionBias, 0, compBufferSize);

        char buffer[1024 + 32000];

        std::vector<Matcher::result_t> resultsBc;
        resultsBc.reserve(300);

        Matcher::result_t resultAc;
        resultAc.backtrace.reserve(par.maxSeqLen + 1);
        std::map<unsigned int, IntervalArray *> interval;
        std::stack<IntervalArray *> intervalBuffer;
        std::vector<Matcher::result_t> resultsAc;
        resultsAc.reserve(1000);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultAbReader->getSize(); ++i) {
            progress.updateProgress();

            unsigned int queryKey = resultAbReader->getDbKey(i);

            size_t aSeqId = aReader.getId(queryKey);
            aSeq.mapSequence(aSeqId, queryKey, aReader.getData(aSeqId, thread_idx), aReader.getSeqLen(aSeqId));

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
                if(returnAlnRes == false && resultAb.eval > par.evalProfile){
                    continue;
                }
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

                std::stable_sort(resultsBc.begin(), resultsBc.end(), compareHitsByKeyScore);

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
                    // A single target sequence can cover a query just a single time
                    // If a target has the same domain several times, then we only consider one
                    if(interval.find(cSeqKey) != interval.end()) {
                        if(interval[cSeqKey]->doesOverlap(resultAc.qStartPos, resultAc.qEndPos)){
                            continue;
                        }
                    } else {
                        size_t cSeqId = cReader->getId(cSeqKey);
                        cSeq.mapSequence(cSeqId, cSeqKey, cReader->getData(cSeqId, thread_idx), cReader->getSeqLen(cSeqId));
                        rescoreResultByBacktrace(resultAc, aSeq, cSeq, subMat, compositionBias, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
                        if(resultAc.score < -6){ // alignment too bad (fitted on regression benchmark EXPAND)
                            continue;
                        }

                        if(par.expansionMode == Parameters::EXPAND_RESCORE_BACKTRACE){
                            resultAc.eval = evaluer->computeEvalue(resultAc.score, aSeq.L);
                            resultAc.score = static_cast<int>(evaluer->computeBitScore(resultAc.score)+0.5);
                            resultAc.seqId = Util::computeSeqId(par.seqIdMode, resultAc.seqId, aSeq.L, cSeq.L, resultAc.backtrace.size());
                        }else{
                            resultAc.eval = resultAb.eval;
                            resultAc.score = resultAb.score;
                            resultAc.seqId = resultAb.seqId;
                        }
                        float queryCov = SmithWaterman::computeCov(resultAc.qStartPos, resultAc.qEndPos, resultAc.qLen);
                        float targetCov = SmithWaterman::computeCov(resultAc.dbStartPos, resultAc.dbEndPos, resultAc.dbLen);
                        bool hasCov = Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov);
                        bool hasSeqId = resultAc.seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                        bool hasEvalue = (resultAc.eval <= par.evalThr);
                        bool hasAlnLen = (static_cast<int>(resultAc.alnLength) >= par.alnLenThr);
                        if(hasCov && hasSeqId && hasEvalue && hasAlnLen){
                            if (returnAlnRes == false) {
                                seqSet.emplace_back(cSeq.numSequence, cSeq.numSequence + cSeq.L);
                            }
                            resultsAc.emplace_back(resultAc);
                            if(intervalBuffer.size() == 0){
                                interval[cSeqKey] = new IntervalArray();
                            } else {
                                interval[cSeqKey] = intervalBuffer.top();
                                intervalBuffer.pop();
                            }
                            interval[cSeqKey]->insert(resultAc.qStartPos, resultAc.qEndPos);
                        }

                    }
                }
                resultsBc.clear();
            }
            for (std::map<unsigned int, IntervalArray *>::iterator it = interval.begin(); it != interval.end(); it++ )
            {
                it->second->reset();
                intervalBuffer.push(it->second);
            }
            interval.clear();

            if (returnAlnRes) {
                SORT_SERIAL(resultsAc.begin(), resultsAc.end(), Matcher::compareHits);
                writer.writeStart(thread_idx);
                for (size_t j = 0; j < resultsAc.size(); ++j) {
                    size_t len = Matcher::resultToBuffer(buffer, resultsAc[j], true, true);
                    writer.writeAdd(buffer, len, thread_idx);
                }
                writer.writeEnd(queryKey, thread_idx);
                resultsAc.clear();
            } else {
                MultipleAlignment::MSAResult res = aligner->computeMSA(&aSeq, seqSet, resultsAc, true);
                resultsAc.clear();
                size_t filteredSetSize = par.filterMsa == false ? res.setSize
                                                                : filter->filter(res.setSize, res.centerLength, (int)(par.covMSAThr * 100), (int)(par.qid * 100), par.qsc, (int)(par.filterMaxSeqId * 100), par.Ndiff, (const char **) res.msaSequence, true);
                PSSMCalculator::Profile pssmRes = calculator->computePSSMFromMSA(filteredSetSize, aSeq.L, (const char **) res.msaSequence, par.wg);
                if (par.maskProfile == true) {
                    masker->mask(aSeq, pssmRes);
                }
                pssmRes.toBuffer(aSeq, subMat, result);
                writer.writeData(result.c_str(), result.length(), queryKey, thread_idx);
                result.clear();
                MultipleAlignment::deleteMSA(&res);
                seqSet.clear();
            }
        }
        free(compositionBias);
        if (returnAlnRes == false) {
            delete aligner;
            if (filter != NULL) {
                delete filter;
            }
            delete calculator;
            delete masker;
        }
        while(intervalBuffer.size()){
            delete intervalBuffer.top();
            intervalBuffer.pop();
        }
    }

    writer.close(returnAlnRes == false);
    if (probMatrix != NULL) {
        delete probMatrix;
    }
    if (evaluer != NULL) {
        delete evaluer;
    }
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

int expandaln(int argc, const char **argv, const Command& command) {
    return expandaln(argc, argv, command, true);
}

int expand2profile(int argc, const char **argv, const Command& command) {
    return expandaln(argc, argv, command, false);
}
