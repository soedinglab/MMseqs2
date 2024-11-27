#ifndef MMSEQS_PSSM_H
#define MMSEQS_PSSM_H

#include <cstddef>
#include <string>
#include <CSProfile.h>
#include <vector>

#include "Matcher.h"


class BaseMatrix;
class Sequence;

class PSSMCalculator {
public:

    struct Profile{
        char * pssm;
        float * prob;
        const float * neffM;
#ifdef GAP_POS_SCORING
        const uint8_t *gDel;
        const uint8_t *gIns;
#endif
//        std::string consensus;
        unsigned char * consensus;

#ifdef GAP_POS_SCORING
        Profile(char *pssm, float *prob, float *neffM, const uint8_t *gDel, const uint8_t *gIns, unsigned char * consensus)
                : pssm(pssm), prob(prob), neffM(neffM), gDel(gDel), gIns(gIns), consensus(consensus) {}
#else
        Profile(char *pssm, float *prob, float *neffM, unsigned char * consensus)
                : pssm(pssm), prob(prob), neffM(neffM), consensus(consensus) {}
#endif
        void toBuffer(const unsigned char* centerSequence, size_t centerSeqLen, BaseMatrix& subMat, std::string& result);
        void toBuffer(Sequence& centerSequence, BaseMatrix& subMat, std::string& result);
    };

    PSSMCalculator(SubstitutionMatrix *subMat, size_t maxSeqLength, size_t maxSetSize, int pcmode,
                   MultiParam<PseudoCounts> pca, MultiParam<PseudoCounts> pcb
#ifdef GAP_POS_SCORING
                   , int gapOpen
                   , int gapPseudoCount
#endif
    );

    ~PSSMCalculator();

    Profile computePSSMFromMSA(size_t setSize, size_t queryLength, const char **msaSeqs, bool wg, float scoreBias);
#ifdef GAP_POS_SCORING
    Profile computePSSMFromMSA(size_t setSize, size_t queryLength, const char **msaSeqs, const std::vector<Matcher::result_t> &alnResults, bool wg, float scoreBias);
#endif

    void printProfile(size_t queryLength);
    void printPSSM(size_t queryLength);
    void profileToString(std::string& result, size_t queryLength);

    // prepare pseudocounts
    static void preparePseudoCounts(float *frequency, float *frequency_with_pseudocounts, size_t entrySize, size_t queryLength, const float **R);

    // compute pseudocounts from Neff_M -p log(p) per column
    static void computePseudoCounts(float *profile, float *frequency, float *frequency_with_pseudocounts, size_t entrySize, float *Neff_M, size_t length,float pca, float pcb);

    // Compute weight for sequence based on "Position-based Sequence Weights' (1994)
    static void computeSequenceWeights(float *seqWeight, size_t queryLength, size_t setSize, const char **msaSeqs);

    // compute position-specific scoring matrix PSSM score
    // 1.) convert PFM to PPM (position probability matrix)
    //     Both PPMs assume statistical independence between positions in the pattern
    // 2.) PSSM Log odds score
    //     M_{aa,pos}={log(M_{aa,pos} / b_{aa}).
    static void computeLogPSSM(BaseMatrix *subMat, char *pssm, const float *profile, float bitFactor, size_t queryLength, float scoreBias);

private:
    BaseMatrix* subMat;

    // cs profiles
    CSProfile * ps;

    // contains sequence weights (global)
    float * seqWeight;

    // sum of sequence weights
    float seqWeightTotal;

    // contains MSA AA matchWeight
    float * matchWeight;

    // counts for cs pseudo counts
    float * counts;

    // contains MSA AA pseudocount weight
    float * pseudocountsWeight;

    // Entropy of MSA
    float * Neff_M;

    // Profile of MSA
    float * profile;

    // PSSM contains log odds PSSM values
    char * pssm;

    // Consensus sequence
    unsigned char * consensusSequence;

#ifdef GAP_POS_SCORING
    // position-specific deletion penalties
    uint8_t *gDel;

    // position-specific gap open penalties for insertions
    uint8_t *gIns;

    // preallocated memory for computing of gap penalties
    std::vector<float> gapWeightsIns;

    // default gap opening penalty
    int gapOpen;

    // pseudo count for calculation of gap opening penalties
    int gapPseudoCount;
#endif

    // number of sequences in subalignment i (only for DEBUGGING)
    int *nseqs;

    // weight contribution value for each sequence
    float **w_contrib;
    // backing aligned memory
    unsigned char *w_contrib_backing;

    // weight of sequence k in column i, calculated from subalignment i
    float *wi;

    // number of different amino acids
    int *naa;

    float **f;

    int **n;
    // backing aligned memory
    unsigned char *n_backing;

    size_t maxSeqLength;
    size_t maxSetSize;

    // compute the Neff_M per column -p log(p)
    void computeNeff_M(float *frequency, float *seqWeight, float *Neff_M, size_t queryLength, size_t setSize, char const **msaSeqs);

    void computeMatchWeights(float * matchWeight, float * seqWeight, size_t setSize, size_t queryLength, const char **msaSeqs);

    void computeContextSpecificWeights(float * matchWeight, float *seqWeight, float * Neff_M, size_t queryLength, size_t setSize, const char **msaSeqs);

    int pcmode;
    MultiParam<PseudoCounts> pca;
    MultiParam<PseudoCounts> pcb;

    void computeConsensusSequence(unsigned char * consensusSeq, float *frequency, size_t queryLength, double *back, char *num2aa);
//    std::string computeConsensusSequence(float *pDouble, size_t queryLength, double *back, char *num2aa);
    void increaseSetSize(size_t newSetSize);

#ifdef GAP_POS_SCORING
    // compute position-specific gap penalties for both deletions and insertions
    void computeGapPenalties(size_t queryLength, size_t setSize, const char **msaSeqs, const std::vector<Matcher::result_t> &alnResults);
#endif

    void fillCounteProfile(float *counts, float *matchWeight, float *Neff_M, size_t queryLength);
};


#endif //MMSEQS_PSSM_H
