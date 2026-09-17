//
// Created by mad on 12/15/15.
//

#ifndef MMSEQS_DIAGONALMATCHER_H
#define MMSEQS_DIAGONALMATCHER_H

#include "SubstitutionMatrix.h"
#include "CacheFriendlyOperations.h"
#include "SequenceLookup.h"
class UngappedAlignment {

public:

    UngappedAlignment(const unsigned int maxSeqLen, BaseMatrix *substitutionMatrix,
                      SequenceLookup *sequenceLookup,
                      const unsigned char *dbRemap = NULL);

    ~UngappedAlignment();

    void createProfile(Sequence *seq, float *biasCorrection,
                       const unsigned char *numSeqOverride = NULL);

    // This function computes the diagonal score for each CounterResult object
    // it assigns the diagonal score to the CounterResult object
    void align(CounterResult *results,
               size_t resultSize);

    int scoreSingelSequenceByCounterResult(CounterResult &result);

    int scoreSingleSequence(std::pair<const unsigned char *, const unsigned int> dbSeq,
                            unsigned short diagonal,
                            unsigned short minDistToDiagonal);

    inline short getQueryBias() {
        return 0;
    }

private:
    const static unsigned int DIAGONALCOUNT = 0xFFFF + 1;
    unsigned int *score_arr;
    char *queryProfile;
    unsigned int queryLen;
    CounterResult ** diagonalMatches;
    unsigned char * diagonalCounter;
    char * aaCorrectionScore;
    BaseMatrix *subMatrix;
    SequenceLookup *sequenceLookup;
    const unsigned char *dbRemap;       // remap table for packed bytes from lookup (NULL = no remap)
    unsigned char *remapBuffer;         // pre-allocated buffer for remapped sequences
    unsigned int remapBufferSize;

    // this function bins the hit_t by diagonals by distributing each hit in an array of 256 * 16(sse)/32(avx2)
    // the function scoreDiagonalAndUpdateHits is called for each bin that reaches its maximum (16 or 32)
    template <bool HasRemap>
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
    template <bool HasRemap>
    void scoreDiagonalAndUpdateHits(const char *queryProfile, const unsigned int queryLen,
                                    const short diagonal, CounterResult **hits, const unsigned int hitSize);

    // Fetch db sequence, applying remap when HasRemap=true. Zero overhead when HasRemap=false.
    template <bool HasRemap>
    inline const unsigned char* getDbSeq(DBLocalId seqId, unsigned int &outLen, unsigned int bufferOffset = 0);

    unsigned short distanceFromDiagonal(const unsigned short diagonal);

    int computeSingelSequenceScores(const char *queryProfile, const unsigned int queryLen,
                                    std::pair<const unsigned char *, const unsigned int> &dbSeq,
                                    int diagonal, unsigned int minDistToDiagonal);

    int computeLongScore(const char * queryProfile, unsigned int queryLen,
                         std::pair<const unsigned char *, const unsigned int> &dbSeq,
                         unsigned short diagonal);


};


#endif //MMSEQS_DIAGONALMATCHER_H
