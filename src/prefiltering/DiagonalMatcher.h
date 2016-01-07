//
// Created by mad on 12/15/15.
//

#ifndef MMSEQS_DIAGONALMATCHER_H
#define MMSEQS_DIAGONALMATCHER_H

#include <SubstitutionMatrix.h>
#include <simd.h>
#include <CountInt32Array.h>
#include "QueryScore.h"
#include "SequenceLookup.h"
class DiagonalMatcher {

public:

    DiagonalMatcher(const unsigned int maxSeqLen, BaseMatrix *substitutionMatrix,
                    SequenceLookup *sequenceLookup);

    ~DiagonalMatcher();

    // This function computes the diagonal score for each CounterResult object
    // it assigns the diagonal score to the CounterResult object
    void processQuery(Sequence *seq, float *compositionBias, CounterResult *results,
                      size_t resultSize, unsigned int thr);

private:
    const static unsigned int DIAGONALCOUNT = 0xFFFF + 1;
    const static unsigned int PROFILESIZE = 32;

    unsigned int *score_arr;
    unsigned char *vectorSequence;
    char *queryProfile;
    CounterResult *** diagonalMatches;
    unsigned char * diagonalCounter;
    BaseMatrix *subMatrix;
    SequenceLookup *sequenceLookup;

    // this function bins the hit_t by diagonals by distributing each hit in an array of 256 * 16(sse)/32(avx2)
    // the function scoreDiagonalAndUpdateHits is called for each bin that reaches its maximum (16 or 32)
    void computeScores(const char *queryProfile,
                       const unsigned int queryLen,
                       CounterResult * results,
                       const size_t resultSize,
                       const short bias,
                       unsigned int thr);
    // scores a single diagonal
    const int scalarDiagonalScoring(const char *profile,
                                    const int bias,
                                    const unsigned int seqLen,
                                    const unsigned char *dbSeq);

    // scores the diagonal of  16/32 db sequences in parallel
    const simd_int vectorDiagonalScoring(const char *profile,
                                         const char bias, const unsigned int seqLen, const unsigned char *dbSeq);

    std::pair<unsigned char *, unsigned int> mapSequences(std::pair<unsigned char *, unsigned int> * seqs, unsigned int seqCount);

    // calles vectorDiagonalScoring or scalarDiagonalScoring depending on the hitSize
    // and updates diagonalScore of the hit_t objects
    void scoreDiagonalAndUpdateHits(const char *queryProfile, const unsigned int queryLen,
                                    const unsigned short diagonal, CounterResult **hits, const unsigned int hitSize,
                                    const short bias);

#ifdef AVX2
    const __m256i Shuffle(const __m256i &value, const __m256i &shuffle);
#endif

    unsigned short distanceFromDiagonal(const unsigned short diagonal);

    unsigned int computeSplit(unsigned int second, const unsigned int i);

    void extractScores(unsigned int *score_arr, simd_int score);

    unsigned char normalizeScore(const unsigned char score, const unsigned int len);

    short createProfile(Sequence *seq, float *biasCorrection, short **subMat, int alphabetSize);
};


#endif //MMSEQS_DIAGONALMATCHER_H
