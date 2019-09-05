//
// Created by mad on 3/15/15.
//

#ifndef MMSEQS_MULTIPLEALIGNMENT_H
#define MMSEQS_MULTIPLEALIGNMENT_H

#include <cstddef>
#include <vector>

#include "Matcher.h"

class MultipleAlignment {
public:
    enum alignment_element {
        ANY=20,   //number representing an X (any amino acid) internally
        NAA=20,   //number of amino acids (0-19)
        GAP=21,   //number representing a gap internally
        ENDGAP=22 //number representing a ignored gaps (for some calculations like gap percentage)
    };

    struct MSAResult{
        size_t msaSequenceLength;
        size_t centerLength;
        size_t setSize;
        char ** msaSequence;
        char * keep;
		std::vector<Matcher::result_t> alignmentResults;

        MSAResult(size_t msaSequenceLength, size_t centerLength, size_t setSize, char **msa)
                : msaSequenceLength(msaSequenceLength), centerLength(centerLength), setSize(setSize), msaSequence(msa), keep(NULL) {}

        MSAResult(size_t msaSequenceLength, size_t centerLength, size_t setSize, char **msa,std::vector<Matcher::result_t> alignmentResults)
                : msaSequenceLength(msaSequenceLength), centerLength(centerLength), setSize(setSize), msaSequence(msa), keep(NULL), alignmentResults(alignmentResults) {}
    };


    MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat, Matcher *aligner);

    ~MultipleAlignment();

    // Compute center star multiple alignment from sequence input
    MSAResult computeMSA(Sequence *centerSeq, const std::vector<Sequence *> &edgeSeqs, bool noDeletionMSA);

    MSAResult computeMSA(Sequence *centerSeq, const std::vector<Sequence *> &edgeSeqs, const std::vector<Matcher::result_t> &alignmentResults, bool noDeletionMSA);

    static void print(MSAResult msaResult, SubstitutionMatrix * subMat);

    // init aligned memory for the MSA
    static char *initX(int len);

    // clean memory for MSA
    static void deleteMSA(MultipleAlignment::MSAResult * res);
	
	
private:
    Matcher * aligner;
    BaseMatrix * subMat;

    size_t maxSeqLen;
    size_t maxSetSize;
    size_t maxMsaSeqLen;
    unsigned int * queryGaps;

    std::vector<Matcher::result_t> computeBacktrace(Sequence *centerSeq, const std::vector<Sequence *> &sequences);

    void computeQueryGaps(unsigned int *queryGaps, Sequence *centerSeq, const std::vector<Sequence *> &seqs, const std::vector<Matcher::result_t> &alignmentResults);

    size_t updateGapsInCenterSequence(char **msaSequence, Sequence *centerSeq, bool noDeletionMSA);

    void updateGapsInSequenceSet(char **msaSequence, size_t centerSeqSize, const std::vector<Sequence *> &seqs,
                                 const std::vector<Matcher::result_t> &alignmentResults, unsigned int *queryGaps,
                                 bool noDeletionMSA);

    MSAResult singleSequenceMSA(Sequence *centerSeq);
	
};


#endif //MMSEQS_MULTIPLEALIGNMENT_H
