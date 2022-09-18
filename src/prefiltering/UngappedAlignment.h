//
// Created by mad on 12/15/15.
//

#ifndef MMSEQS_DIAGONALMATCHER_H
#define MMSEQS_DIAGONALMATCHER_H

#include "SubstitutionMatrix.h"
#include "simd.h"
#include "CacheFriendlyOperations.h"
#include "SequenceLookup.h"
class UngappedAlignment {

public:

    UngappedAlignment(const unsigned int maxSeqLen, BaseMatrix *substitutionMatrix,
                      SequenceLookup *sequenceLookup);

    ~UngappedAlignment();

    // This function computes the diagonal score for each CounterResult object
    // it assigns the diagonal score to the CounterResult object
    void processQuery(Sequence *seq, float *compositionBias, CounterResult *results,
                      size_t resultSize);

    int scoreSingelSequenceByCounterResult(CounterResult &result);

    int scoreSingleSequence(std::pair<const unsigned char *, const unsigned int> dbSeq,
                            unsigned short diagonal,
                            unsigned short minDistToDiagonal);

    inline short getQueryBias() {
        return 0;
    }

#ifdef AVX2
    static __m256i Shuffle(const __m256i & value, const __m256i & shuffle)
    {
        const __m256i K0 = _mm256_setr_epi8(
                (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70,
                (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0);
        const __m256i K1 = _mm256_setr_epi8(
                (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0, (char)0xF0,
                (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70, (char)0x70);
        return _mm256_or_si256(_mm256_shuffle_epi8(value, _mm256_add_epi8(shuffle, K0)),
                               _mm256_shuffle_epi8(_mm256_permute4x64_epi64(value, 0x4E), _mm256_add_epi8(shuffle, K1)));
    }
#endif

private:
    const static unsigned int DIAGONALCOUNT = 0xFFFF + 1;
#ifdef AVX2
    const static unsigned int DIAGONALBINSIZE = 8;
#else
    const static unsigned int DIAGONALBINSIZE = 4;
#endif
    unsigned int *score_arr;
    char *queryProfile;
    unsigned int queryLen;
    CounterResult ** diagonalMatches;
    unsigned char * diagonalCounter;
    char * aaCorrectionScore;
    BaseMatrix *subMatrix;
    SequenceLookup *sequenceLookup;

    // this function bins the hit_t by diagonals by distributing each hit in an array of 256 * 16(sse)/32(avx2)
    // the function scoreDiagonalAndUpdateHits is called for each bin that reaches its maximum (16 or 32)
    void computeScores(const char *queryProfile,
                       const unsigned int queryLen,
                       CounterResult * results,
                       const size_t resultSize);

    // scores a single diagonal
    int scalarDiagonalScoring(const char *profile,
                              const unsigned int seqLen,
                              const unsigned char *dbSeq);

    template <unsigned int T>
    void unrolledDiagonalScoring(const char * profile,
                                 const unsigned int * seqLen,
                                 const unsigned char ** dbSeq,
                                 unsigned int * max);

    // calles vectorDiagonalScoring or scalarDiagonalScoring depending on the hitSize
    // and updates diagonalScore of the hit_t objects
    void scoreDiagonalAndUpdateHits(const char *queryProfile, const unsigned int queryLen,
                                    const short diagonal, CounterResult **hits, const unsigned int hitSize);

    unsigned short distanceFromDiagonal(const unsigned short diagonal);

    void extractScores(unsigned int *score_arr, simd_int score);

    void createProfile(Sequence *seq, float *biasCorrection, short **subMat);

    int computeSingelSequenceScores(const char *queryProfile, const unsigned int queryLen,
                                    std::pair<const unsigned char *, const unsigned int> &dbSeq,
                                    int diagonal, unsigned int minDistToDiagonal);

    int computeLongScore(const char * queryProfile, unsigned int queryLen,
                         std::pair<const unsigned char *, const unsigned int> &dbSeq,
                         unsigned short diagonal);


};


#endif //MMSEQS_DIAGONALMATCHER_H

