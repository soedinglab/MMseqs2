//
// Created by mad on 3/15/15.
//

#ifndef MMSEQS_MULTIPLEALIGNMENT_H
#define MMSEQS_MULTIPLEALIGNMENT_H


#include <cstddef>
#include <vector>

#include "Matcher.h"

class SubstitutionMatrix;
class Sequence;


class MultipleAlignment {
public:
    enum alignment_element {
        ANY=20,   //number representing an X (any amino acid) internally
        NAA=20,   //number of amino acids (0-19)
        GAP=21    //number representing a gap internally
    };

    struct MSAResult{
        size_t msaSequenceLength;
        size_t centerLength;
        size_t setSize;
        char ** msaSequence;
		std::vector<Matcher::result_t> alignmentResults;

        MSAResult(size_t msaSequenceLength, size_t centerLength, size_t setSize, char **msa)
                : msaSequenceLength(msaSequenceLength), centerLength(centerLength), setSize(setSize), msaSequence(msa) {}

        MSAResult(size_t msaSequenceLength, size_t centerLength, size_t setSize, char **msa,std::vector<Matcher::result_t> alignmentResults)
                : msaSequenceLength(msaSequenceLength), centerLength(centerLength), setSize(setSize), msaSequence(msa), alignmentResults(alignmentResults) {}
    };


    MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat,
                                             Matcher *aligner);

    ~MultipleAlignment();
    // Compute center star multiple alignment from sequence input
    MultipleAlignment::MSAResult computeMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs, bool noDeletionMSA);
    static void print(MSAResult msaResult, SubstitutionMatrix * subMat);

    // init aligned memory for the MSA
    static char *initX(int len);

    MSAResult computeMSA(Sequence *pSequence, std::vector<Sequence *> vector, std::vector<Matcher::result_t> vector1,
                         bool i);
    // clean memory for MSA
    static void deleteMSA(MultipleAlignment::MSAResult * res);
	
	
private:
    Matcher * aligner;
    BaseMatrix * subMat;

    size_t maxSeqLen;
    size_t maxSetSize;
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

    MSAResult singleSequenceMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs);
	
};


#endif //MMSEQS_MULTIPLEALIGNMENT_H
