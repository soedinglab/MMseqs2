//
// Created by mad on 3/24/15.
//
#include "PSSMCalculator.h"
#include "simd.h"
#include "MathUtil.h"
#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "Util.h"
#include "Debug.h"
#include "MultipleAlignment.h"

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
    // compute context specific counts
    //computeContextSpecificWeights(seqWeight, queryLength, setSize, msaSeqs);
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
            profile[idx] = subMat->getBitFactor() * MathUtil::flog2(aaProb / subMat->pBack[aa]) + scoreBias;
            const float pssmVal = profile[pos * Sequence::PROFILE_AA_SIZE + aa];
            pssm[idx] = static_cast<char>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
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
    for (int k = 0; k < length; k++){
        sum += array[k];
    }
    if (sum != 0.0f) {
        float fac = 1.0 / sum;
        for (int i = 0; i < length; i++) {
            array[i] *= fac;
        }
    } else if (def_array) {
        for (int i = 0; i < length; i++) {
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
                sum -= freq_pos_aa * MathUtil::flog2(freq_pos_aa);
            }
        }
        Neff_HMM += MathUtil::fpow2(sum);
    }
    Neff_HMM /= queryLength;
    float Nlim = fmax(10.0, Neff_HMM + 1.0);    // limiting Neff
    float scale = MathUtil::flog2((Nlim - Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
    for (size_t pos = 0; pos < queryLength; pos++) {
        float w_M = -1.0 / setSize;
        for (size_t k = 0; k < setSize; ++k){
            if ( msaSeqs[k][pos] != MultipleAlignment::GAP) {
                w_M += seqWeight[k];
            }
        }
        Neff_M[pos] = (w_M < 0) ? 1.0 : Nlim - (Nlim - 1.0) * MathUtil::fpow2(scale * w_M);
//        fprintf(stderr,"M  i=%3i  ncol=---  Neff_M=%5.2f  Nlim=%5.2f  w_M=%5.3f  Neff_M=%5.2f\n",pos,Neff_HMM,Nlim,w_M,Neff_M[pos]);
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
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
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
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
                const unsigned int aa_pos = msaSeqs[k][pos];
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
            if (msaSeqs[k][pos] != MultipleAlignment::GAP) {
                if(distinct_aa_count == 0){
                    seqWeight[k] += 0.0;
                } else {
                    const unsigned int aa_pos = msaSeqs[k][pos];
//                    std::cout << "k="<< k << "\t";
                    if(aa_pos < Sequence::PROFILE_AA_SIZE){ // Treat score of X with other amino acid as 0.0
                        seqWeight[k] += 1.0f / (float(nl[aa_pos]) * float(distinct_aa_count) * (float(number_res[k]) + 30.0f));
//                        std::cout << number_res[k] << "\t" << distinct_aa_count << "\t" <<  nl[aa_pos] << "\t";
                    }
                    printf("%0.6f\n", seqWeight[k]);
                }
                // ensure that each residue of a short sequence contributes as much as a residue of a long sequence:
                // contribution is proportional to one over sequence length nres[k] plus 30.
            }
        }
    }
    std::cout << setSize << std::endl;
    NormalizeTo1(seqWeight, setSize);
//    std::cout << " Seq. Weight: " << std::endl;
//    for (size_t k = 0; k < setSize; ++k) {
//        std::cout << " k="<< k << "\t" << seqWeight[k] << std::endl;
//    }
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
            if(msaSeqs[k][pos] != MultipleAlignment::GAP){
                unsigned int aa_pos = msaSeqs[k][pos];
                if(aa_pos < Sequence::PROFILE_AA_SIZE) { // Treat score of X with other amino acid as 0.0
                    matchWeight[pos * Sequence::PROFILE_AA_SIZE + aa_pos] += seqWeight[k];
                }
            }
        }
        NormalizeTo1(&matchWeight[pos * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE, subMat->pBack);
    }
}

void PSSMCalculator::computeContextSpecificWeights(float *wg, size_t queryLength, size_t setSize,
                                                   const char **msaSeqs) {
//    const char ** X = msaSeqs;
//    const int NAA = 20;
//    const int ANY = 20;
//    unsigned int NAA_VECSIZE = ((NAA+ 3 + VECSIZE_INT - 1) / VECSIZE_INT) * VECSIZE_INT; // round NAA+3 up to next multiple of VECSIZE_INT
//    int ** n = new int*[queryLength + 2];
//    for (int j = 0; j < queryLength; j++)
//        n[j] = (int *) malloc_simd_int(NAA_VECSIZE * sizeof(int));
//    for (int j = 0; j <= queryLength; j++)
//        for (int a = 0; a < NAA + 3; ++a)
//            n[j][a] = 0;
//    float ** w_contrib = new float*[queryLength + 2];
//    float * wi = new float[queryLength];
//    int * naa = new int[queryLength];
//    for (int j = 0; j < queryLength; j++){
//        w_contrib[j] = (float *) malloc_simd_int(NAA_VECSIZE * sizeof(float));
//        memset(w_contrib[j], 0,NAA_VECSIZE * sizeof(int));
//    }
//
//    float Neff_HMM = 0.0f;
//    float * Neff = new float[queryLength];
//    //////////////////////////////////////////////////////////////////////////////////////////////
//    // Main loop through alignment columns
//    for (int i = 0; i < queryLength; i++)  // Calculate wi[k] at position i as well as Neff[i]
//    {
//        bool change = 0;
//        // Check all sequences k and update n[j][a] and ri[j] if necessary
//        for (int k = 0; k < setSize; ++k) {
//            // Update amino acid and GAP / ENDGAP counts for sequences with AA in i-1 and GAP/ENDGAP in i or vice versa
//            if (X[k][i - 1] >= ANY && X[k][i] < ANY) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
//                change = 1;
//                nseqi++;
//                for (int j = 1; j <= queryLength; ++j)
//                    n[j][(int) X[k][j]]++;
//            } else if (X[k][i - 1] < ANY && X[k][i] >= ANY) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
//                change = 1;
//                nseqi--;
//                for (int j = 1; j <= queryLength; ++j)
//                    n[j][(int) X[k][j]]--;
//            }
//        }  //end for (k)
//        nseqs[i] = nseqi;
//
//        // Only if subalignment changed we need to update weights wi[k] and Neff[i]
//        if (change) {
//            // We gained a factor ~8.0 for the following computation of weights
//            // and profile by exchanging the inner two loops (j, k => k, j)
//            // and precomputing the weight contributions w_contrib[j][a].
//            // M. Steinegger and J. Soeding (29 July 2014)
//
//            // Initialize weights and numbers of residues for subalignment i
//            int ncol = 0;
//            for (int k = 0; k < setSize; ++k)
//                wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;
//
//            // Find min and max borders between which > fraction MAXENDGAPFRAC of sequences in subalignment contain an aa
//            int jmin;
//            int jmax;
//            for (jmin = 1; jmin <= queryLength && n[jmin][ENDGAP] > MAXENDGAPFRAC * nseqi;
//                 ++jmin) {
//            };
//            for (jmax = queryLength; jmax >= 1 && n[jmax][ENDGAP] > MAXENDGAPFRAC * nseqi;
//                 --jmax) {
//            };
//            ncol = jmax - jmin + 1;
//
//            // Check whether number of columns in subalignment is sufficient
//            if (ncol < NCOLMIN) {
//                // Take global weights
//                for (int k = 0; k < setSize; ++k)
//                    wi[k] = (X[k][i] < ANY)? wg[k] : 0.0f;
//            } else {
//                // Count number of different amino acids in column j
//                for (int j = jmin; j <= jmax; ++j){
//                    naa[j] = 0;
//                    for (int a = 0; a < ANY; ++a){
//                        naa[j] += (n[j][a] ? 1 : 0);
//                    }
//                }
//                // Compute the contribution of amino acid a to the weight
//                //for (a = 0; a < ANY; ++a)
//                //      w_contrib[j][a] = (n[j][a] > 0) ? 1.0/ float(naa[j]*n[j][a]): 0.0f;
//                for (int j = jmin; j <= jmax; ++j) {
//                    simd_float naa_j = simdi32_i2f(simdi32_set(naa[j]));
//                    const simd_int *nj = (const simd_int *) n[j];
//                    const int aa_size = (ANY + VECSIZE_INT - 1) / VECSIZE_INT;
//                    for (int a = 0; a < aa_size; ++a) {
//                        simd_float nja = simdi32_i2f(simdi_load(nj + a));
//                        simd_float res = simdf32_mul(nja, naa_j);
//                        simdf32_store(w_contrib[j] + (a * VECSIZE_INT), simdf32_rcp(res));
//                    }
//                    for (int a = ANY; a < NAA + 3; ++a)
//                        w_contrib[j][a] = 0.0f;  // set non-amino acid values to 0 to avoid checking in next loop for X[k][j]<ANY
//                }
//
//                // Compute pos-specific weights wi[k]
//                for (int k = 0; k < setSize; ++k) {
//                    if (X[k][i] >= ANY)
//                        continue;
//                    for (int j = jmin; j <= jmax; ++j)  // innermost, time-critical loop; O(L*setSize*L)
//                        wi[k] += w_contrib[j][(int) X[k][j]];
//                }
//            }
//
//            // Calculate Neff[i]
//            Neff[i] = 0.0;
//
//            // Allocate and reset amino acid frequencies
//            for (int j = jmin; j <= jmax; ++j)
//                memset(f[j], 0, ANY * sizeof(float));
//
//            // Update f[j][a]
//            for (int k = 0; k < setSize; ++k) {
//                if (X[k][i] >= ANY)
//                    continue;
//                for (int j = jmin; j <= jmax; ++j)  // innermost loop; O(L*setSize*L)
//                    f[j][(int) X[k][j]] += wi[k];
//            }
//
//            // Add contributions to Neff[i]
//            for (int j = jmin; j <= jmax; ++j) {
//                NormalizeTo1(f[j], NAA);
//                for (a = 0; a < 20; ++a)
//                    if (f[j][a] > 1E-10)
//                        Neff[i] -= f[j][a] * fast_log2(f[j][a]);
//            }
//
//            if (ncol > 0)
//                Neff[i] = fpow2(Neff[i] / ncol);
//            else
//                Neff[i] = 1.0;
//
//        }
//
//        else  //no update was necessary; copy values for i-1
//        {
//            Neff[i] = Neff[i - 1];
//        }
//    }
//
//    // Calculate amino acid frequencies q->f[i][a] from weights wi[k]
//    for (int a = 0; a < 20; ++a)
//        q->f[i][a] = 0;
//    for (int k = 0; k < setSize; ++k)
//            q->f[i][(int) X[k][i]] += wi[k];
//    NormalizeTo1(q->f[i], NAA, pb);

}
