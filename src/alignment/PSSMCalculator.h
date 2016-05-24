#ifndef MMSEQS_PSSM_H
#define MMSEQS_PSSM_H

#include <cstddef>
#include <string>

class SubstitutionMatrix;

class PSSMCalculator {
public:

    PSSMCalculator(SubstitutionMatrix *subMat, size_t maxSeqLength, float pca, float pcb);

    ~PSSMCalculator();

    std::pair<const char *, std::string> computePSSMFromMSA(size_t setSize, size_t queryLength, const char **msaSeqs,
                                    bool wg);

    void printProfile(size_t queryLength);
    void printPSSM(size_t queryLength);

private:
    SubstitutionMatrix * subMat;

    // contains sequence weights (global)
    float * seqWeight;

    // contains MSA AA matchWeight
    float * matchWeight;

    // contains MSA AA pseudocount weight
    float * pseudocountsWeight;

    // Entropy of MSA
    float * Neff_M;

    // Profile of MSA
    float * profile;

    // PSSM contains log odds PSSM values
    char * pssm;

    // pseudocount matrix (mem aligned)
    float * R;

    // number of sequences in subalignment i (only for DEBUGGING)
    int *nseqs;

    // weight contribution value for each sequence
    float **w_contrib;

    // weight of sequence k in column i, calculated from subalignment i
    float *wi;

    // number of different amino acids
    int *naa;

    size_t maxSeqLength;

    // compute position-specific scoring matrix PSSM score
    // 1.) convert PFM to PPM (position probability matrix)
    //     Both PPMs assume statistical independence between positions in the pattern
    // 2.) PSSM Log odds score
    //     M_{aa,pos}={log(M_{aa,pos} / b_{aa}).
    void computeLogPSSM(char *pssm, float *profile, size_t queryLength, float scoreBias);


    // normalize a fector to 1.0
    float NormalizeTo1(float *array, int length, double const *def_array);

    // prepare pseudocounts
    void preparePseudoCounts(float *frequency, float *frequency_with_pseudocounts, size_t queryLength, const float *R);

    // compute the Neff_M per column -p log(p)
    void computeNeff_M(float *frequency, float *seqWeight, float *Neff_M, size_t queryLength, size_t setSize, char const **msaSeqs);

    // Compute weight for sequence based on "Position-based Sequence Weights' (1994)
    void computeSequenceWeights(float *seqWeight, size_t queryLength, size_t setSize, const char **msaSeqs);

    // compute pseudocounts from Neff_M -p log(p) per column
    void computePseudoCounts(float *profile, float *frequency, float *frequency_with_pseudocounts, size_t length,float pca, float pcb);

    void computeMatchWeights(float * matchWeight, float * seqWeight, size_t setSize, size_t queryLength, const char **msaSeqs);

    void computeContextSpecificWeights(float * matchWeight, float *seqWeight, float * Neff_M, size_t queryLength, size_t setSize, const char **msaSeqs);

    float pca;
    float pcb;

    std::string computeConsensusSequence(float *pDouble, size_t queryLength, double *back, char *int2aa);
};


#endif //MMSEQS_PSSM_H
