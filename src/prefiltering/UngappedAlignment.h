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
        return bias;
    }

private:
    const static unsigned int DIAGONALCOUNT = 0xFFFF + 1;
    const static unsigned int PROFILESIZE = 32;

    unsigned int *score_arr;
    unsigned char *vectorSequence;
    char *queryProfile;
    unsigned int queryLen;
    short bias;
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
                       const size_t resultSize,
                       const short bias);
    // scores a single diagonal
    int scalarDiagonalScoring(const char *profile,
                                    const int bias,
                                    const unsigned int seqLen,
                                    const unsigned char *dbSeq);

    // scores the diagonal of  16/32 db sequences in parallel
    simd_int vectorDiagonalScoring(const char *profile,
                                         const char bias, const unsigned int seqLen, const unsigned char *dbSeq);

    std::pair<unsigned char *, unsigned int> mapSequences(std::pair<unsigned char *, unsigned int> * seqs, unsigned int seqCount);

    // calles vectorDiagonalScoring or scalarDiagonalScoring depending on the hitSize
    // and updates diagonalScore of the hit_t objects
    void scoreDiagonalAndUpdateHits(const char *queryProfile, const unsigned int queryLen,
                                    const short diagonal, CounterResult **hits, const unsigned int hitSize,
                                    const short bias);

    unsigned short distanceFromDiagonal(const unsigned short diagonal);

    void extractScores(unsigned int *score_arr, simd_int score);

    short createProfile(Sequence *seq, float *biasCorrection, short **subMat, int alphabetSize);

    unsigned int diagonalLength(const short diagonal, const unsigned int len, const unsigned int second);

    int computeSingelSequenceScores(const char *queryProfile, const unsigned int queryLen,
                                    std::pair<const unsigned char *, const unsigned int> &dbSeq,
                                    int diagonal, unsigned int minDistToDiagonal, short bias);

    int computeLongScore(const char * queryProfile, unsigned int queryLen,
                         std::pair<const unsigned char *, const unsigned int> &dbSeq,
                         unsigned short diagonal, short bias);


};


#endif //MMSEQS_DIAGONALMATCHER_H
