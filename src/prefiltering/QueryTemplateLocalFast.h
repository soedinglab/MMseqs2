//
// Created by mad on 5/26/15.
//

#ifndef MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
#define MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H

#include "CacheFriendlyOperations.h"
#include "QueryTemplateMatcher.h"
#include "UngappedAlignment.h"


class QueryTemplateLocalFast : public virtual QueryTemplateMatcher {
public:
    QueryTemplateLocalFast(BaseMatrix *m, IndexTable *indexTable,
                           unsigned int *seqLens, short kmerThr,
                           double kmerMatchProb, int kmerSize, size_t dbSize,
                           unsigned int maxSeqLen, unsigned int effectiveKmerSize,
                           size_t maxHitsPerQuery, bool aaBiasCorrection, bool diagonalScoring, unsigned int minDiagScoreThr);
    ~QueryTemplateLocalFast();

    // returns result for the sequence
    // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
    std::pair<hit_t *, size_t>  matchQuery(Sequence * seq, unsigned int identityId);

    // find duplicates in the diagonal bins
    size_t evaluateBins(IndexEntryLocal **hitsByIndex, CounterResult *output,
                        size_t outputSize, unsigned short indexFrom, unsigned short indexTo, bool computeTotalScore);

    void updateScoreBins(CounterResult *result, size_t elementCount);
protected:

    unsigned int maxDbMatches;
    unsigned int dbSize;

    // result hit buffer
    //CacheFriendlyOperations * diagonalMatcher;
    static const unsigned int L2_CACH_SIZE = 262144;
    unsigned int activeCounter;
#define CacheFriendlyOperations(x)  CacheFriendlyOperations<x> * cachedOperation##x
    CacheFriendlyOperations(2);
    CacheFriendlyOperations(4);
    CacheFriendlyOperations(8);
    CacheFriendlyOperations(16);
    CacheFriendlyOperations(32);
    CacheFriendlyOperations(64);
    CacheFriendlyOperations(128);
    CacheFriendlyOperations(256);
    CacheFriendlyOperations(512);
    CacheFriendlyOperations(1024);
    CacheFriendlyOperations(2048);
#undef CacheFriendlyOperations

    // matcher for diagonal
    UngappedAlignment *ungappedAlignment;

    // score distribution of current query
    unsigned int *scoreSizes;

    // result hit buffer
    hit_t *resList;

    // i position to hits pointer
    IndexEntryLocal **indexPointer;

    // keeps data in inner loop
    IndexEntryLocal * __restrict databaseHits;

    // evaluated bins
    CounterResult * foundDiagonals;

    // last data pointer (for overflow check)
    IndexEntryLocal * lastSequenceHit;

    // the following variables are needed to calculate the Z-score computation
    double mu;

    //log match prob (mu) of poisson distribution
    double logMatchProb;

    //pre computed score factorials
    // S_fact = score!
    double *logScoreFactorial;

    // max seq. per query
    size_t maxHitsPerQuery;

    //pointer to seqLens
    float *seqLens;

    // match sequence against the IndexTable
    size_t match(Sequence *seq, float *pDouble);

    // extract result from databaseHits
    std::pair<hit_t *, size_t> getResult(CounterResult * results,
                                         size_t resultSize,
                                         const int l, const unsigned int id,
                                         const unsigned short thr,
                                         const bool diagonalScoring);
    // compute double hits
    size_t getDoubleDiagonalMatches();

    float *compositionBias;

    // diagonal scoring active
    bool diagonalScoring;
    unsigned int minDiagScoreThr;
    // size of max diagonalMatcher result objects
    size_t counterResultSize;

    void initDiagonalMatcher(size_t dbsize, unsigned int maxDbMatches);

    void deleteDiagonalMatcher(unsigned int activeCounter);

    size_t mergeElements(bool diagonalScoring, CounterResult *foundDiagonals, size_t hitCounter);

    size_t keepMaxScoreElementOnly(CounterResult *foundDiagonals, size_t resultSize);
};

#endif //MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
