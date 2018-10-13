#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "Matcher.h"
#include "EvalueComputation.h"
#include "CompressedA3M.h"
#include "Sequence.h"
#include "Alignment.h"
#include "SubstitutionMatrix.h"

#include <cassert>

#ifdef OPENMP
#include <omp.h>
#endif

class BacktraceTranslator {
public:
    BacktraceTranslator() {
        //  AB state, BC state -> AC state, increment in AB cigar, increment in BC cigar

        // Covis' Rules
//        transitions[M][M] = Transition('M',1,1);
//        transitions[I][M] = Transition('I',1,0);
//        transitions[D][M] = Transition('D',1,1);
//
//        transitions[M][D] = Transition('D',0,1);
//        transitions[I][D] = Transition('D',0,1);
//        transitions[D][D] = Transition('D',0,1);
//
//        transitions[M][I] = Transition('I',1,1);
//        transitions[I][I] = Transition('I',1,0);
//        transitions[D][I] = Transition('\0',1,1);

        // Eli's Rules
        transitions[M][M] = Transition('M',1,1);
        transitions[I][M] = Transition('I',1,1);
        transitions[D][M] = Transition('\0', 1,0);

        transitions[M][D] = Transition('D',1,1);
        transitions[I][D] = Transition('\0', 1,1);
        transitions[D][D] = Transition('D',1,0);

        transitions[M][I] = Transition('I',0,1);
        transitions[I][I] = Transition('I',0,1);
        transitions[D][I] = Transition('D',1,0);
    }

    void translateResult(const Matcher::result_t& resultAB, Matcher::result_t& resultBC) {
        std::string btAC;
        btAC.reserve(std::max(resultAB.backtrace.size(), resultBC.backtrace.size()));

        int minPos = std::min(resultAB.dbStartPos, resultBC.qStartPos);
        size_t startBA = resultAB.dbStartPos - minPos;
        size_t startBC = resultBC.qStartPos - minPos;

        ssize_t offset = startBC - startBA;

        size_t currentAB = 0;
        size_t currentBC = 0;
        if (offset > 0) {
//            btAC.append(offset, 'D');
            currentAB = offset;
            resultBC.qStartPos += offset;
        }

        if (offset < 0) {
//            btAC.append(-offset, 'I');
            currentBC = -offset;
            resultBC.dbStartPos -= offset;
        }

        unsigned int lastM = 0;
        unsigned int qAlnLength = 0;
        unsigned int i = 0;
        while (currentAB < resultAB.backtrace.size() && currentBC < resultBC.backtrace.size()) {
            Transition& t = transitions[mapState(resultAB.backtrace[currentAB])][mapState(resultBC.backtrace[currentBC])];
            switch (t.newState) {
                case '\0':
                    i--;
                    break;
                case 'M':
                    lastM = i;
                    // FALLTHROUGH
                case 'D':
                    qAlnLength++;
                    // FALLTHROUGH
                case 'I':
                    btAC.append(1, t.newState);;
                    break;
                default:
                    Debug(Debug::ERROR) << "Invalid alignment state.\n";
                    EXIT(EXIT_FAILURE);

            }
            currentAB += t.incrementAB;
            currentBC += t.incrementBC;
            i++;
        }
        btAC.resize(lastM - 1);
        resultBC.qEndPos = resultBC.qStartPos + qAlnLength - 1;
        resultBC.backtrace.swap(btAC);
    }


//    void translate(size_t startBA, size_t startBC, const char* btAB, size_t btABLength, const char* btBC, size_t btBCLength, std::string& btAC) {
//        btAC.clear();
//
//        size_t currentAB = 0;
//        size_t currentBC = 0;
//        ssize_t offset = startBC - startBA;
//        if (offset > 0) {
//            btAC.append(offset, 'D');
//            currentAB = offset;
//        }
//
//        if (offset < 0) {
//            btAC.append(-offset, 'I');
//            currentBC = -offset;
//        }
//
//        while (currentAB < btABLength && currentBC < btBCLength) {
//            Transition& t = transitions[mapState(btAB[currentAB])][mapState(btBC[currentBC])];
//            if (t.newState != '\0') {
//                assert(t.newState == 'M' || t.newState == 'I' || t.newState == 'D');
//                btAC.append(1, t.newState);
//            }
//            currentAB += t.incrementAB;
//            currentBC += t.incrementBC;
//        }
//    }

private:
    struct Transition {
        Transition() : newState('\0'), incrementAB(0), incrementBC(0) {};
        Transition(char newState, uint8_t incrementAB, uint8_t incrementBC) : newState(newState), incrementAB(incrementAB), incrementBC(incrementBC) {};
        char newState;
        uint8_t incrementAB;
        uint8_t incrementBC;
    };

    Transition transitions[3][3];

    enum State {
        M = 0,
        I,
        D
    };

    State mapState(char state) {
        if (state == 'M') {
            return M;
        } else if (state == 'I') {
            return I;
        } else if  (state == 'D') {
            return D;
        } else {
            Debug(Debug::ERROR) << "Invalid alignment state.\n";
            EXIT(EXIT_FAILURE);
        }
    }
};

void rescoreResultByBacktrace(Matcher::result_t &result, Sequence &qSeq, Sequence &tSeq, SubstitutionMatrix &subMat, float *compositionBias, EvalueComputation &evaluer, int gapOpen, int gapExtend, int seqIdMode) {
    size_t qPos = result.qStartPos;
    size_t tPos = result.dbStartPos;
    int score = 0;
    char lastState = '\0';
    int identities = 0;
    for (size_t i = 0; i < result.backtrace.size(); ++i) {
        char state = result.backtrace[i];
        if (state == 'M') {
            if (tSeq.getSeqType() == Sequence::HMM_PROFILE) {
                score += tSeq.profile_for_alignment[qSeq.int_sequence[qPos] * tSeq.L + tPos];
            } else if (qSeq.getSeqType() == Sequence::HMM_PROFILE) {
                score += qSeq.profile_for_alignment[tSeq.int_sequence[tPos] * qSeq.L + qPos];
            } else {
                score += subMat.subMatrix[qSeq.int_sequence[qPos]][tSeq.int_sequence[tPos]] + compositionBias[qPos];
            }
            identities += qSeq.int_sequence[qPos] == tSeq.int_sequence[tPos] ? 1 : 0;
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
    result.score = static_cast<short>(evaluer.computeBitScore(score)+0.5);
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
    par.parseParameters(argc, argv, command, 5);

    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str());
    queryReader.open(DBReader<unsigned int>::NOSORT);
    const int queryDbType = queryReader.getDbtype();
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        queryReader.readMmapedDataInMemory();
    }

    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str());
    targetReader.open(DBReader<unsigned int>::NOSORT);
    const int targetDbType = targetReader.getDbtype();
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        targetReader.readMmapedDataInMemory();
    }

    if (queryDbType == Sequence::HMM_PROFILE && targetDbType == Sequence::HMM_PROFILE) {
        Debug(Debug::ERROR) << "Profile-profile is currently not supported.\n";
        return EXIT_FAILURE;
    }

    DBReader<unsigned int> *resultReader = NULL;
    DBReader<unsigned int> *ca3mSequenceReader = NULL;
    if (FileUtil::fileExists((par.db3 + "_ca3m.ffdata").c_str())) {
        Debug(Debug::INFO) << "Result database: " << par.db3 << "_ca3m\n";
        resultReader = new DBReader<unsigned int>((par.db3 + "_ca3m.ffdata").c_str(), (par.db3 + "_ca3m.ffindex").c_str());
        resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        ca3mSequenceReader = new DBReader<unsigned int>((par.db3 + "_sequence.ffdata").c_str(), (par.db3 + "_sequence.ffindex").c_str());
        ca3mSequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    } else {
        Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
        resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
        resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    }

    Debug(Debug::INFO) << "Expansion result database: " << par.db4 << "\n";
    DBReader<unsigned int> expansionReader(par.db4.c_str(), par.db4Index.c_str());
    expansionReader.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        expansionReader.readMmapedDataInMemory();
    }

    Debug(Debug::INFO) << "Output database: " << par.db5 << "\n";
    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), par.threads);
    writer.open();

    BacktraceTranslator translator;
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, par.scoreBias);
    EvalueComputation evaluer(targetReader.getAminoAcidDBSize(), &subMat, par.gapOpen, par.gapExtend);

    Debug(Debug::INFO) << "Computing expanded alignment result...\n";
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence qSeq(par.maxSeqLen, queryDbType, &subMat, 0, false, par.compBiasCorrection);
        Sequence tSeq(par.maxSeqLen, targetDbType, &subMat, 0, false, par.compBiasCorrection);
        float *compositionBias = new float[par.maxSeqLen];
        memset(compositionBias, 0, sizeof(float) * par.maxSeqLen);

        std::vector<Matcher::result_t> expanded;
        expanded.reserve(300);

        std::vector<Matcher::result_t> results;
        results.reserve(1000);

        char buffer[1024];

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader->getSize(); ++i) {
            Debug::printProgress(i);
            unsigned int queryKey = resultReader->getDbKey(i);

            size_t querySeqId = queryReader.getId(queryKey);
            qSeq.mapSequence(querySeqId, queryKey, queryReader.getData(querySeqId));

            if(par.compBiasCorrection == true && queryDbType == Sequence::AMINO_ACIDS){
                SubstitutionMatrix::calcLocalAaBiasCorrection(&subMat, qSeq.int_sequence, qSeq.L, compositionBias);
            }

            char *data = resultReader->getData(i);
            while (*data != '\0') {
                Matcher::result_t resultAB = Matcher::parseAlignmentRecord(data, false);

                if (resultAB.backtrace.size() == 0) {
                    Debug(Debug::ERROR) << "Alignment must contain a backtrace.\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int targetKey = resultAB.dbKey;
                size_t targetId = expansionReader.getId(targetKey);

                size_t targetSeqId = targetReader.getId(targetKey);
                tSeq.mapSequence(targetSeqId, targetKey, targetReader.getData(targetSeqId));

                if (ca3mSequenceReader != NULL) {
                    unsigned int key;
                    CompressedA3M::extractMatcherResults(key, expanded, expansionReader.getData(targetId),
                                                         expansionReader.getSeqLens(targetId), *ca3mSequenceReader, false);
                } else {
                    Matcher::readAlignmentResults(expanded, expansionReader.getData(targetId), false);
                }
                for (size_t k = 0; k < expanded.size(); ++k) {
                    Matcher::result_t &resultBC = expanded[k];
                    if (resultBC.backtrace.size() == 0) {
                        Debug(Debug::ERROR) << "Alignment must contain a backtrace.\n";
                        EXIT(EXIT_FAILURE);
                    }

                    translator.translateResult(resultAB, resultBC);
                    if (resultBC.backtrace.size() == 0) {
                        continue;
                    }

                    rescoreResultByBacktrace(resultBC, qSeq, tSeq, subMat, compositionBias,
                                             evaluer, par.gapOpen, par.gapExtend, par.seqIdMode);

                    if (Alignment::checkCriteria(resultBC, false, par.evalThr, par.seqIdThr, par.covMode, par.covThr)) {
                        results.emplace_back(resultBC);
                    }
                }
                expanded.clear();
                data = Util::skipLine(data);
            }

            std::vector<Matcher::result_t> *finalResults = &results;
            if (par.expansionMode == 1) {
                // keep only the best hit to same target
                std::sort(results.begin(), results.end(), compareHitsByKeyEvalScore);
                ssize_t lastKey = -1;
                for (size_t j = 0; j < results.size(); ++j) {
                    Matcher::result_t& res = results[i];
                    if (res.dbKey != lastKey) {
                        expanded.emplace_back(res);
                    }
                    lastKey = res.dbKey;
                }
                finalResults = &expanded;
            }
            std::sort(finalResults->begin(), finalResults->end(), Matcher::compareHits);

            writer.writeStart(thread_idx);
            for (size_t j = 0; j < finalResults->size(); ++j) {
                size_t len = Matcher::resultToBuffer(buffer, (*finalResults)[j], true, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(queryKey, thread_idx);
            expanded.clear();
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
