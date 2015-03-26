//
// Created by mad on 3/24/15.
//

#ifndef _MMSEQS_PSSM_H_
#define _MMSEQS_PSSM_H_


#include <SubstitutionMatrix.h>
#include <stddef.h>
#include "MultipleAlignment.h"

class PSSM {
public:

    PSSM(SubstitutionMatrix *subMat, size_t maxSeqLength);

    ~PSSM();

    void computePSSMFromMSA(MultipleAlignment::MSAResult msaResult);
    void print(size_t queryLength);
private:
    SubstitutionMatrix * subMat;
    size_t maxSeqLength;
    unsigned char * pssm;
    float * profile;
    bool * gapPosition;
    // compute position frequency matrix
    // M_{aa,pos}=\frac{1}{N} \sum_{i=1}^N I(X_{i,pos}=aa)
    void computeFrequencyMatrix(size_t setSize, char const **msaSeqs, size_t msaSeqLength);

    // compute position-specific scoring matrix PSSM score
    // 1.) convert PFM to PPM (position probability matrix)
    //     Both PPMs assume statistical independence between positions in the pattern
    // 2.) PSSM Log odds score
    //     M_{aa,pos}={log(M_{aa,pos} / b_{aa}).
    void computeLogPSSM(float *profile, size_t queryLength);

    float * seqWeight;
    float * frequency;
    float * frequency_with_pseudocounts;
    float * Neff_M;

    // normalize a fector to 1.0
    float NormalizeTo1(float *array, int length, double const *def_array);


    void PreparePseudocounts(float *frequency, float *frequency_with_pseudocounts, size_t queryLength, const float **R);

    // compute the Neff_M per column -p log(p)
    void computeNeff_M(float *frequency, float *seqWeight, float *Neff_M, size_t queryLength, size_t setSize, char const **msaSeqs);

    // Compute weight for sequence based on "Position-based Sequence Weights' (1994)
    void computeSequenceWeights(float *seqWeight, size_t queryLength, size_t setSize, const char **msaSeqs);

    // compute pseudocounts from Neff_M -p log(p) per column
    void computePseudoCounts(float *profile, float *frequency, float *frequency_with_pseudocounts, size_t length);

};


#endif //_MMSEQS_PSSM_H_
