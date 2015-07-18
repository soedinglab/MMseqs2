//
// Created by mad on 5/26/15.
//

#ifndef MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
#define MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H

#include <sys/cdefs.h>
#include "CountInt32Array.h"
#include "QueryTemplateMatcher.h"


class QueryTemplateMatcherExactMatch : public virtual QueryTemplateMatcher {
public:
    QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                   unsigned int *seqLens, short kmerThr,
                                   double kmerMatchProb, int kmerSize, size_t dbSize,
                                   unsigned int maxSeqLen, size_t maxHitsPerQuery);
    ~QueryTemplateMatcherExactMatch();

    // returns result for the sequence
    // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
    std::pair<hit_t *, size_t>  matchQuery(Sequence * seq, unsigned int identityId);

    // find duplicates in the diagonal bins
    size_t evaluateBins(CounterResult *inputOutput, size_t i);

    void updateScoreBins(CounterResult *result, size_t elementCount);
protected:

    const static unsigned int MAX_DB_MATCHES = 16777216;

    // result hit buffer
    CountInt32Array * counter;

    // score distribution of current query
    unsigned int *scoreSizes;

    // result hit buffer
    hit_t *resList;

    // keeps data in inner loop
    CounterResult * __restrict databaseHits;

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
    unsigned int *seqLens;

    // match sequence against the IndexTable
    void match(Sequence* seq);

    // extract result from databaseHits
    std::pair<hit_t *, size_t> getResult(const int l, const unsigned int id,
                                         const unsigned short thr);
};

#endif //MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
