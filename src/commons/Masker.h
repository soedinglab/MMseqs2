#ifndef MMSEQS_MASKER_H
#define MMSEQS_MASKER_H

#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "PSSMCalculator.h"
#include <cctype>

class Masker {
public:
    Masker(BaseMatrix &subMat);

    ~Masker();

    int maskSequence(Sequence & seq, bool maskTantan, double maskProb,
                     bool maskLowerCaseLetter, int maskNrepeating);

    void maskPssm(Sequence& centerSequence, float maskProb, PSSMCalculator::Profile& pssmRes);

    void applySoftmasking(unsigned char *charSequence, const unsigned char * numSequence, unsigned int seqLen);

    char maskLetterNum;

private:
    int maskRepeats(unsigned char *numSequence, const unsigned int seqLen, int maskNrepeating, char maskChar);

    void finalizeMasking(unsigned char * numSequence, const unsigned int seqLen);

    BaseMatrix &subMat;
    ProbabilityMatrix probMatrix;

    unsigned char * charSequence;
    size_t maxSeqLen;
};
#endif
