//
// Created by mad on 12/15/15.
//

#ifndef MMSEQS_DIAGONALMATCHER_H
#define MMSEQS_DIAGONALMATCHER_H

#include <simd.h>
#include "QueryScore.h"


class DiagonalMatcher {

public:

    DiagonalMatcher(const unsigned int maxSeqLen, SubstitutionMatrix *substitutionMatrix,
                        SequenceLookup *sequenceLookup);

    ~DiagonalMatcher();

    void processQuery(Sequence *query, hit_t *results, size_t resultSize);

private:
    const static unsigned int DIAGONALCOUNT = 256;
    const static unsigned int  PROFILESIZE = 32;

    unsigned int *score_arr;
    unsigned char *vectorSequence;
    char *queryProfile;
    hit_t * diagonalMatches[DIAGONALCOUNT][VECSIZE_INT*4];
    unsigned char diagonalCounter[DIAGONALCOUNT];
    SubstitutionMatrix *subMatrix;
    SequenceLookup *sequenceLookup;

    void computeScores(const char *queryProfile,
                       const unsigned int queryLen,
                       hit_t *results,
                       const size_t resultSize,
                       const short bias);

    const int scalarDiagonalScoring(const char *profile,
                                             const unsigned int profileSize,
                                             const int bias,
                                             const unsigned int seqLen,
                                             const unsigned char *dbSeq);

    std::pair<unsigned char *, unsigned int> mapSequences(hit_t *seqIds, unsigned int seqCount);

    void scoreDiagonalAndUpdateHits(const char *queryProfile, const unsigned int queryLen,
                                    const unsigned char diagonal, hit_t *hits, const unsigned int hitSize,
                                    const short bias);

    const simd_int vectorDiagonalScoring(const char *profile, const unsigned int profileSize,
                                         const char bias, const unsigned int seqLen, const unsigned char *dbSeq);

    const __m256i Shuffle(const __m256i &value, const __m256i &shuffle);
};


#endif //MMSEQS_DIAGONALMATCHER_H
