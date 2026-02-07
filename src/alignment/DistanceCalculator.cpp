#include "DistanceCalculator.h"

#include "simd.h"

template<typename T>
unsigned int DistanceCalculator::computeInverseHammingDistance(const T *seq1, const T *seq2, unsigned int length){
    unsigned int diff = 0;
    unsigned int simdBlock = length/(VECSIZE_INT*4);
    simd_int * simdSeq1 = (simd_int *) seq1;
    simd_int * simdSeq2 = (simd_int *) seq2;
    for (unsigned int pos = 0; pos < simdBlock; pos++ ) {
        simd_int seq1vec = simdi_loadu(simdSeq1+pos);
        simd_int seq2vec = simdi_loadu(simdSeq2+pos);
        // int _mm_movemask_epi8(__m128i a) creates 16-bit mask from most significant bits of
        // the 16 signed or unsigned 8-bit integers in a and zero-extends the upper bits.
        simd_int seqComparision = simdi8_eq(seq1vec, seq2vec);
        int res = simdi8_movemask(seqComparision);
        diff += __builtin_popcount(res);  // subtract positions that should not contribute to coverage
    }
    // compute missing rest
    for (unsigned int pos = simdBlock*(VECSIZE_INT*4); pos < length; pos++ ) {
        diff += (seq1[pos] == seq2[pos]);
    }
    return diff;
}

template unsigned int DistanceCalculator::computeInverseHammingDistance<char>(const char *seq1, const char *seq2, unsigned int length);
