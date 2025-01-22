#include "Masker.h"
#include <algorithm> // for std::toupper

Masker::Masker(BaseMatrix &s) : subMat(s), probMatrix(s)
{
    maxSeqLen = 1;
    charSequence = (unsigned char *)malloc(maxSeqLen * sizeof(char));
    maskLetterNum = subMat.aa2num[(int)'X'];
}

Masker::~Masker() {
    free(charSequence);
}

int Masker::maskSequence(Sequence & seq, bool maskTantan, double maskProb,
                         bool maskLowerCaseLetter, int maskNrepeats) {

    int maskedResidues = 0;

    if(maskTantan){
        // 1. Apply tantan masking without influencing by repeat mask
        maskedResidues += tantan::maskSequences(seq.numSequence,
                                            seq.numSequence + seq.L,
                                            50 /*maxCycleLength*/,
                                            probMatrix.probMatrixPointers,
                                            0.005 /*repeatProb*/,
                                            0.05 /*repeatEndProb*/,
                                            0.9 /*repeatOffsetProbDecay*/,
                                            0, 0,
                                            maskProb /*minMaskProb*/,
                                            probMatrix.hardMaskTable);
    }
    if( maskNrepeats > 0){
        // 2. Generate the mask for repeats
        maskedResidues += maskRepeats(seq.numSequence, seq.L, maskNrepeats, maskLetterNum);
    }
    // 3. Handle lowercase masking
    if(maskLowerCaseLetter){
        if ((Parameters::isEqualDbtype(seq.getSequenceType(), Parameters::DBTYPE_AMINO_ACIDS) ||
             Parameters::isEqualDbtype(seq.getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES))) {
            const char *charSeq = seq.getSeqData();
            for (int i = 0; i < seq.L; i++) {
                if (std::islower((unsigned char)charSeq[i])) {
                    seq.numSequence[i] = maskLetterNum; // Apply masking
                    maskedResidues++;
                }
            }
        }
    }
    // 4. Finalize masking
    if(maskTantan || maskNrepeats || maskLowerCaseLetter){
        finalizeMasking(seq.numSequence, seq.L);
    }
    return maskedResidues;
}

void Masker::maskPssm(Sequence& centerSequence, float maskProb, PSSMCalculator::Profile& pssmRes) {
    if ((size_t)centerSequence.L > maxSeqLen) {
        maxSeqLen = sizeof(char) * centerSequence.L * 1.5;
        charSequence = (unsigned char*)realloc(charSequence, maxSeqLen);
    }
    memcpy(charSequence, centerSequence.numSequence, sizeof(unsigned char) * centerSequence.L);
    tantan::maskSequences(charSequence, charSequence + centerSequence.L,
                          50 /*options.maxCycleLength*/,
                          probMatrix.probMatrixPointers,
                          0.005 /*options.repeatProb*/,
                          0.05 /*options.repeatEndProb*/,
                          0.9 /*options.repeatOffsetProbDecay*/,
                          0, 0,
                          maskProb /*options.minMaskProb*/,
                          probMatrix.hardMaskTable);

    for (int pos = 0; pos < centerSequence.L; pos++) {
        if (charSequence[pos] == maskLetterNum) {
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + aa] = -1;
            }
        }
    }
}


int Masker::maskRepeats(unsigned char * numSequence, const unsigned int seqLen, int maskNrepeating, char maskChar) {

    unsigned int repeatCount = 0;
    int startOfRepeat = -1;
    char previousChar = '\0';
    int maskedResidues = 0; // Counter for masked residues

    for (unsigned int pos = 0; pos < seqLen; ++pos) {
        char currentChar = numSequence[pos];

        if (currentChar == previousChar) {
            repeatCount++;
        } else {
            if (repeatCount > (unsigned int)maskNrepeating) {
                for (unsigned int i = startOfRepeat; i < pos; ++i) {
                    numSequence[i] = maskChar;
                    maskedResidues++;
                }
            }
            repeatCount = 1;
            startOfRepeat = pos;
            previousChar = currentChar;
        }
    }

    // Handle the last run
    if (repeatCount > (unsigned int)maskNrepeating) {
        for (unsigned int i = startOfRepeat; i < seqLen; ++i) {
            numSequence[i] = maskChar;
            maskedResidues++;
        }
    }

    return maskedResidues;
}

void Masker::finalizeMasking(unsigned char * numSequence, const unsigned int seqLen) {
    unsigned char maskChar = probMatrix.hardMaskTable[0];

    for (unsigned int i = 0; i < seqLen; i++) {
        unsigned char code = numSequence[i];
        numSequence[i] = (code == maskChar || code == maskLetterNum) ? maskLetterNum : numSequence[i];
    }
}

void Masker::applySoftmasking(unsigned char *charSequence, const unsigned char * num_sequence, unsigned int seqLen) {
    for (unsigned int pos = 0; pos < seqLen; pos++) {
        // If masked, lowercase (soft) or uppercase (hard) could be applied here if needed.
        // For simplicity, we treat maskChar as masked and others as uppercase:
        charSequence[pos] = (num_sequence[pos] == maskLetterNum)
                            ? (char)std::tolower(charSequence[pos])
                            : (char)std::toupper(charSequence[pos]);
    }
}
