//
// Created by mad on 3/15/15.
//

#ifndef _MMSEQS_MULTIPLEALIGNMENT_H_
#define _MMSEQS_MULTIPLEALIGNMENT_H_


#include <stddef.h>
#include <vector>
#include "Sequence.h"
#include "Matcher.h"
#include "SubstitutionMatrix.h"


class MultipleAlignment {
public:
    struct MSAResult{
        size_t msaSequenceLength;
        size_t centerLength;
        size_t setSize;
        const char ** msaSequence;

        MSAResult(size_t msaSequenceLength, size_t centerLength, size_t setSize, const char **msa)
                : msaSequenceLength(msaSequenceLength), centerLength(centerLength), setSize(setSize), msaSequence(msa) {}
    };


    MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat);

    ~MultipleAlignment();

    MultipleAlignment::MSAResult computeMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs, bool noDeletionMSA);
    static void print(MSAResult msaResult);
private:
    Matcher * aligner;
    BaseMatrix * subMat;
    char *  msaData;
    char ** msaSequence;
    size_t maxSeqLen;
    unsigned int * queryGaps;
};


#endif //_MMSEQS_MULTIPLEALIGNMENT_H_
