#ifndef PSSMMASKER_H
#define PSSMMASKER_H

#include <cstring>
#include "tantan.h"
#include "BaseMatrix.h"
#include "PSSMCalculator.h"

class PSSMMasker {
public:
    PSSMMasker(size_t maxSeqLen, ProbabilityMatrix& probMatrix, BaseMatrix& subMat) : maxSeqLen(maxSeqLen), probMatrix(probMatrix), subMat(subMat), xAmioAcid(subMat.aa2num[static_cast<int>('X')]) {
        charSequence = (char*)malloc(sizeof(char) * maxSeqLen);
    }

    ~PSSMMasker() {
        free(charSequence);
    }

    void mask(Sequence& centerSequence, PSSMCalculator::Profile& pssmRes) {
        if ((size_t)centerSequence.L > maxSeqLen) {
            maxSeqLen = sizeof(char) * centerSequence.L * 1.5;
            charSequence = (char*)realloc(charSequence, maxSeqLen);
        }
        memcpy(charSequence, centerSequence.numSequence, sizeof(unsigned char) * centerSequence.L);
        tantan::maskSequences(charSequence, charSequence + centerSequence.L,
                              50 /*options.maxCycleLength*/,
                              probMatrix.probMatrixPointers,
                              0.005 /*options.repeatProb*/,
                              0.05 /*options.repeatEndProb*/,
                              0.9 /*options.repeatOffsetProbDecay*/,
                              0, 0,
                              0.9 /*options.minMaskProb*/,
                              probMatrix.hardMaskTable);

        for (int pos = 0; pos < centerSequence.L; pos++) {
            if (charSequence[pos] == xAmioAcid) {
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    pssmRes.prob[pos * Sequence::PROFILE_AA_SIZE + aa] = subMat.pBack[aa] * 0.5;
                }
                pssmRes.consensus[pos] = 'X';
            }
        }
    }
private:
    char *charSequence;
    size_t maxSeqLen;
    ProbabilityMatrix& probMatrix;
    BaseMatrix& subMat;
    const int xAmioAcid;
};

#endif
