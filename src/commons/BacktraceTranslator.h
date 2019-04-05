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

        // Martins Clovis Eli's Rules
        transitions[M][M] = Transition('M',1,1);
        transitions[I][M] = Transition('I',1,1);
        transitions[D][M] = Transition('D', 1,1);

        transitions[M][D] = Transition('D',1,1);
        transitions[I][D] = Transition('\0', 1,1);
        transitions[D][D] = Transition('D',1,1);

        transitions[M][I] = Transition('I',1,1);
        transitions[I][I] = Transition('I',0,1);
        transitions[D][I] = Transition('\0',1,1);
    }

    /* translate results takes an pairwise alignment between sequence AB (resultAB), and BC (resultBC) to infer an alignment between AC (resultAC)
     E.g. alignment example
     ######## Alignment AB
     A: ATT-G--  \
                   MMMIM
     B: ATTTGCA  /
     ######## Alignment BC
     B: ATTTGCA  \
                   MIMMM
     C: --T-GCA  /
     ######### infer    AC
     A: ATTG--   \
                   MM
     C: --TG--   /
    */
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
        int distanceInB = maxB - minB;
        // here we need to align they start position in B
        if (startBab < startBbc) {
            int aOffset = 0;
            int bOffset = 0;
            int btOffset = 0;
            while(bOffset < distanceInB && btOffset < static_cast<int>(resultAB.backtrace.size())){
                bOffset += (resultAB.backtrace[btOffset] == 'M' || resultAB.backtrace[btOffset] == 'D');
                aOffset += (resultAB.backtrace[btOffset] == 'M' || resultAB.backtrace[btOffset] == 'I');
                btOffset++;
            }
            offsetBbc = 0;
            offsetBab = btOffset;
            startAac = startAab + aOffset;
            startCac = startCbc;
        } else if (startBab > startBbc) {
            int bOffset = 0;
            int cOffset = 0;
            int btOffset = 0;
            while(bOffset < distanceInB && btOffset < static_cast<int>(resultBC.backtrace.size())){
                bOffset += (resultBC.backtrace[btOffset] == 'M'  || resultBC.backtrace[btOffset] == 'I');
                cOffset += (resultBC.backtrace[btOffset] == 'M'  || resultBC.backtrace[btOffset] == 'D');
                btOffset++;
            }
            offsetBab = 0;
            offsetBbc = btOffset;
            startAac = startAab;
            startCac = startCbc + cOffset;
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
        int backtraceABSize = static_cast<int>(resultAB.backtrace.size());
        int backtraceBCSize = static_cast<int>(resultBC.backtrace.size());
        while (offsetBab < backtraceABSize && offsetBbc < backtraceBCSize) {
            i++;
            State ab = mapState(resultAB.backtrace[offsetBab]);
            State bc = mapState(resultBC.backtrace[offsetBbc]);
            Transition& t = transitions[ab][bc];
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
                    dbAlnLength++;
                    break;
                case 'I':
                    qAlnLength++;
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
