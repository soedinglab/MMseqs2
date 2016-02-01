//
// Created by mad on 3/15/15.
//

#ifndef MMSEQS_MULTIPLEALIGNMENT_H
#define MMSEQS_MULTIPLEALIGNMENT_H


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


    MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat,
                                             Matcher *aligner);

    ~MultipleAlignment();
    // Compute center star multiple alignment from sequence input
    MultipleAlignment::MSAResult computeMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs, bool noDeletionMSA);
    static void print(MSAResult msaResult);
private:
    Matcher * aligner;
    BaseMatrix * subMat;
    char *  msaData;
    char ** msaSequence;
    size_t maxSeqLen;
    size_t maxMsaSeqLen;
    unsigned int * queryGaps;

    std::vector<Matcher::result_t> computeBacktrace(Sequence *center, std::vector<Sequence *> sequenes,
                                                                       size_t dbSetSize);

    void computeQueryGaps(unsigned int *queryGaps, Sequence *center, std::vector<Sequence *> seqs,
                          std::vector<Matcher::result_t> alignmentResults);

    size_t updateGapsInCenterSequence(char **msaSequence, Sequence *centerSeq, bool noDeletionMSA);

    void updateGapsInSequenceSet(char **centerSeqSize, size_t seqs, std::vector<Sequence *> vector,
                                                    std::vector<Matcher::result_t> queryGaps, unsigned int *noDeletionMSA,
                                                    bool b);

};


#endif //MMSEQS_MULTIPLEALIGNMENT_H
