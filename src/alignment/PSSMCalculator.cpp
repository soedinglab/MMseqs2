//
// Created by mad on 3/24/15.
//
#include "PSSMCalculator.h"
#include <stddef.h>
#include "simd.h"

PSSMCalculator::PSSMCalculator(SubstitutionMatrix *subMat, size_t maxSeqLength) :subMat(subMat) {
    profile            = new float[Sequence::PROFILE_AA_SIZE * maxSeqLength];
    Neff_M             = new float[maxSeqLength];
    seqWeight          = new float[maxSeqLength];
    pssm = new char[Sequence::PROFILE_AA_SIZE * maxSeqLength];
    matchWeight        = (float *) malloc_simd_float(Sequence::PROFILE_AA_SIZE * maxSeqLength * sizeof(float));
    pseudocountsWeight = (float *) malloc_simd_float(Sequence::PROFILE_AA_SIZE * maxSeqLength * sizeof(float));
    R                  = (float *) malloc_simd_float(Sequence::PROFILE_AA_SIZE * Sequence::PROFILE_AA_SIZE * sizeof(float));
    for (size_t aa_i = 0; aa_i < Sequence::PROFILE_AA_SIZE; ++aa_i){
        for (size_t aa_j = 0; aa_j < Sequence::PROFILE_AA_SIZE; ++aa_j){
            R[aa_i * Sequence::PROFILE_AA_SIZE + aa_j] = subMat->subMatrixPseudoCounts[aa_i][aa_j];
        }
    }
}

PSSMCalculator::~PSSMCalculator() {
    delete [] profile;
    delete [] Neff_M;
    delete [] seqWeight;
    delete [] pssm;
    free(matchWeight);
    free(pseudocountsWeight);
    free(R);
}

char const * PSSMCalculator::computePSSMFromMSA(size_t setSize, size_t queryLength, const char **msaSeqs) {

    for(size_t pos = 0; pos < queryLength; pos++){
        if(msaSeqs[0][pos] == '-'){
            Debug(Debug::ERROR) << "Error in computePSSMFromMSA. First sequence of MSA is not allowed ot contain gaps.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // Quick and dirty calculation of the weight per sequence wg[k]
    computeSequenceWeights(seqWeight, queryLength, setSize, msaSeqs);
    // compute matchWeight based on sequence weight
    computeMatchWeights(setSize, queryLength, msaSeqs);
    // compute NEFF_M
    computeNeff_M(matchWeight, seqWeight, Neff_M, queryLength, setSize, msaSeqs);
    // add pseudocounts (compute the scalar product between matchWeight and substitution matrix with pseudo counts)
    preparePseudoCounts(matchWeight, pseudocountsWeight, queryLength, (const float *) R);
    computePseudoCounts(profile, matchWeight, pseudocountsWeight, queryLength);
    // create final Matrix
    computeLogPSSM(pssm, profile, queryLength, -0.2);
    return pssm;
}

void PSSMCalculator::printProfile(size_t queryLength){
    printf("Pos ");
    for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf("%2c    ", subMat->int2aa[aa]);
    }
    printf("\n");
    for(size_t i = 0; i < queryLength; i++){
        printf("%3zu ", i);
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            printf("%03.2f ", profile[i * Sequence::PROFILE_AA_SIZE + aa] );
        }
        printf("\n");
    }
}

void PSSMCalculator::printPSSM(size_t queryLength){
    printf("Pos ");
    for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->int2aa[aa]);
    }
    printf("\n");
    for(size_t i = 0; i <  queryLength; i++) {
        printf("%3zu ", i);
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            printf("%3d ", pssm[i * Sequence::PROFILE_AA_SIZE + aa] );
        }
        printf("\n");
    }
}

void PSSMCalculator::computeLogPSSM(char *pssm, float *profile,
                                    size_t queryLength, float scoreBias) {
    // calculate the substitution matrix
    for(size_t pos = 0; pos < queryLength; pos++) {
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            const float aaProb = profile[pos * Sequence::PROFILE_AA_SIZE + aa];
            const unsigned int idx = pos * Sequence::PROFILE_AA_SIZE + aa;
            profile[idx] = subMat->getBitFactor() * Util::flog2(aaProb / subMat->pBack[aa]) + scoreBias;
            pssm[idx] = (char) floor (profile[pos * Sequence::PROFILE_AA_SIZE + aa] + 0.5);
        }
    }
}

void PSSMCalculator::preparePseudoCounts(float *frequency, float *frequency_with_pseudocounts,
                                         size_t queryLength, float const *R) {
    for (size_t pos = 0; pos < queryLength; pos++) {
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            frequency_with_pseudocounts[pos * Sequence::PROFILE_AA_SIZE + aa] = ScalarProd20(&R[aa * Sequence::PROFILE_AA_SIZE],
                                                                                             &frequency[pos * Sequence::PROFILE_AA_SIZE]);
        }
    }
}

// Normalize a float array such that it sums to one
// If it sums to 0 then assign def_array elements to array (optional)
inline float PSSMCalculator::NormalizeTo1(float* array, int length, const double* def_array = NULL) {
    float sum = 0.0f;
    for (size_t k = 0; k < length; k++){
        sum += array[k];
    }
    if (sum != 0.0f) {
        float fac = 1.0 / sum;
        for (size_t i = 0; i < length; i++) {
            array[i] *= fac;
        }
    } else if (def_array) {
        for (size_t i = 0; i < length; i++) {
            array[i] = def_array[i];
        }
    }
    return sum;
}

void PSSMCalculator::computeNeff_M(float *frequency, float *seqWeight, float *Neff_M,
                                   size_t queryLength, size_t setSize, char const **msaSeqs) {
    float Neff_HMM = 0.0f;
    for (size_t pos = 0; pos < queryLength; pos++) {
        float sum = 0.0f;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa){
            float freq_pos_aa = frequency[pos * Sequence::PROFILE_AA_SIZE + aa];
            if (freq_pos_aa > 1E-10) {
                sum -= freq_pos_aa * Util::flog2(freq_pos_aa);
            }
        }
        Neff_HMM += Util::fpow2(sum);
    }
    Neff_HMM /= queryLength;
    float Nlim = fmax(10.0, Neff_HMM + 1.0);    // limiting Neff
    float scale = Util::flog2((Nlim - Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
    for (size_t pos = 0; pos < queryLength; pos++) {
        float w_M = -1.0 / setSize;
        for (size_t k = 0; k < setSize; ++k){
            if ( msaSeqs[k][pos] != '-') {
                w_M += seqWeight[k];
            }
        }
        Neff_M[pos] = (w_M < 0) ? 1.0 : Nlim - (Nlim - 1.0) * Util::fpow2(scale * w_M);
        //fprintf(stderr,"M  i=%3i  ncol=---  Neff_M=%5.2f  Nlim=%5.2f  w_M=%5.3f  Neff_M=%5.2f\n",pos,Neff_HMM,Nlim,w_M,Neff_M[pos]);
    }
}

void PSSMCalculator::computeSequenceWeights(float *seqWeight, size_t queryLength, size_t setSize, const char **msaSeqs) {
    unsigned int *number_res = new unsigned int[setSize];
    // initialized wg[k] with tiny pseudo counts
    std::fill(seqWeight, seqWeight + setSize,  1e-6);
    // count number of residues per sequence
    for (size_t k = 0; k < setSize; ++k)
    {
        unsigned int nr = 0;
        for (size_t pos = 0; pos < queryLength; pos++) {
            if (msaSeqs[k][pos] != '-') {
                nr++;
            }
        }
        number_res[k] = nr;
    }
    for (size_t pos = 0; pos < queryLength; pos++) {
        int nl[ Sequence::PROFILE_AA_SIZE ];  //nl[a] = number of seq's with amino acid a at position l
        //number of different amino acids (ignore X)
        std::fill(nl, nl + Sequence::PROFILE_AA_SIZE,  0);
        for (size_t k = 0; k < setSize; ++k) {
            if (msaSeqs[k][pos] != '-') {
                const unsigned int aa_pos = subMat->aa2int[(int)msaSeqs[k][pos]];
                if(aa_pos < Sequence::PROFILE_AA_SIZE){
                    nl[aa_pos]++;
                }
            }
        }
        //count distinct amino acids (ignore X)
        int distinct_aa_count = 0;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            if (nl[aa]) {
                ++distinct_aa_count;
            }
        }
//        if(distinct_aa_count == 0){
//            Debug(Debug::ERROR) << "Error in computeSequenceWeights. Distinct amino acid count is 0.\n";
//            EXIT(EXIT_FAILURE);
//        }
        // Compute sequence Weight
        // "Position-based Sequence Weights", Henikoff (1994)
        for (size_t k = 0; k < setSize; ++k) {
            if (msaSeqs[k][pos] != '-') {
                if(distinct_aa_count == 0){
                    seqWeight[k] += 0.0;
                } else {
                    const unsigned int aa_pos = subMat->aa2int[(int)msaSeqs[k][pos]];
                    if(aa_pos < Sequence::PROFILE_AA_SIZE){ // Treat score of X with other amino acid as 0.0
                        seqWeight[k] += 1.0f / (float(nl[aa_pos]) * float(distinct_aa_count) * (float(number_res[k]) + 30.0f));
                    }
                }
                // ensure that each residue of a short sequence contributes as much as a residue of a long sequence:
                // contribution is proportional to one over sequence length nres[k] plus 30.
            }
        }
    }
    NormalizeTo1(seqWeight, setSize);
    delete [] number_res;
}

void PSSMCalculator::computePseudoCounts(float *profile, float *frequency, float *frequency_with_pseudocounts, size_t queryLength) {
    float pca = 1.0f; //TODO
    float pcb = 1.5f; //TODO
    for (size_t pos = 0; pos < queryLength; pos++) {
        float tau = fmin(1.0, pca / (1.0 + Neff_M[pos] / pcb));
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            // compute proportion of pseudo counts and signal
            float pseudoCounts    = tau * frequency_with_pseudocounts[pos * Sequence::PROFILE_AA_SIZE + aa];
            float frequencySignal = (1.0 - tau) * frequency[pos * Sequence::PROFILE_AA_SIZE + aa];
            profile[pos * Sequence::PROFILE_AA_SIZE + aa] = frequencySignal + pseudoCounts;
        }
    }
}

void PSSMCalculator::computeMatchWeights(size_t setSize, size_t queryLength, const char **msaSeqs) {
    for (size_t pos = 0; pos < queryLength; pos++) {
        memset(matchWeight + pos * Sequence::PROFILE_AA_SIZE, 0,
               Sequence::PROFILE_AA_SIZE * sizeof(float));
        for (size_t k = 0; k < setSize; ++k){
            if(msaSeqs[k][pos] != '-'){
                unsigned int aa_pos = subMat->aa2int[(int)msaSeqs[k][pos]];
                if(aa_pos < Sequence::PROFILE_AA_SIZE) { // Treat score of X with other amino acid as 0.0
                    matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa_pos] += seqWeight[k];
                }
            }
        }
        NormalizeTo1(&matchWeight[pos * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE, subMat->pBack);
    }
}