#ifndef BACKTRACE_TRANSLATOR_H
#define BACKTRACE_TRANSLATOR_H

#include "Matcher.h"

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

    void translateResult(const Matcher::result_t& resultAB, const Matcher::result_t& resultBC, Matcher::result_t& resultAC) {
        int startAab = resultAB.qStartPos;
        int startBab = resultAB.dbStartPos;
        int startBbc = resultBC.qStartPos;
        int startCbc = resultBC.dbStartPos;

        int minB = std::min(startBab, startBbc);
        int maxB = std::max(startBab, startBbc);

        int offsetBab;
        int offsetBbc;
        int startAac;
        int startCac;
        if (startBab < startBbc) {
            offsetBab = maxB - minB;
            offsetBbc = 0;
            startAac = startAab + offsetBab;
            startCac = startCbc;
        } else if (startBab > startBbc) {
            offsetBab = 0;
            offsetBbc = maxB - minB;
            startAac = startAab;
            startCac = startCbc + offsetBbc;
        } else {
            offsetBab = 0;
            offsetBbc = 0;
            startAac = startAab;
            startCac = startCbc;
        }

        resultAC.backtrace.clear();

        unsigned int lastM = 0;
        unsigned int qAlnLength = 0;
        unsigned int dbAlnLength = 0;
        unsigned int i = 0;
        while (offsetBab < static_cast<int>(resultAB.backtrace.size()) && offsetBbc < static_cast<int>(resultBC.backtrace.size())) {
            i++;
            Transition& t = transitions[mapState(resultAB.backtrace[offsetBab])][mapState(resultBC.backtrace[offsetBbc])];
            switch (t.newState) {
                case '\0':
//                    i = i > 0 ? i - 1 : 0;
                    i--;
                    goto next;
                case 'M':
                    lastM = i;
                    qAlnLength++;
                    dbAlnLength++;
                    break;
                case 'D':
                    qAlnLength++;
                    break;
                case 'I':
                    dbAlnLength++;
                    break;
                default:
                    Debug(Debug::ERROR) << "Invalid backtrace translation state.\n";
                    EXIT(EXIT_FAILURE);

            }
            resultAC.backtrace.append(1, t.newState);
            next:
            offsetBab += t.incrementAB;
            offsetBbc += t.incrementBC;
        }

        resultAC.dbKey = resultBC.dbKey;
        resultAC.score = resultBC.score;
        resultAC.qcov = resultBC.qcov;
        resultAC.dbcov = resultBC.dbcov;
        resultAC.seqId = resultBC.seqId;
        resultAC.eval = resultBC.eval;
        resultAC.alnLength = resultBC.alnLength;
        resultAC.qStartPos = startAac;
        resultAC.qEndPos = startAac + qAlnLength - 1;
        resultAC.qLen = resultAB.qLen;
        resultAC.dbStartPos = startCac;
        resultAC.dbEndPos = startCac + dbAlnLength - 1;
        resultAC.dbLen = resultBC.dbLen;
        resultAC.backtrace.resize(lastM);
    }

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

#endif
