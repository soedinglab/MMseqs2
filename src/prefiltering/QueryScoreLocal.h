//
//  QueryScoreLocal.h
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef QUERYSCORELOCAL_H
#define QUERYSCORELOCAL_H
#include "QueryScore.h"
#include <vector>
#include <map>
#include <sys/cdefs.h>
#include <zconf.h>


class QueryScoreLocal : public QueryScore {

public:
    QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr, double kmerMatchProb);

    ~QueryScoreLocal();

    // get the list of the sequences with the score > scoreThr
    std::pair<hit_t *, size_t> getResult(const unsigned int querySeqLen, unsigned int scoreThr, const unsigned int identityId);

    void reset();

    void updateScoreBins();

    unsigned int computeScoreThreshold(size_t maxHitsPerQuery);

private:

    //log match prob (mu) of poisson distribution
    double logMatchProb;

    //pre computed score factorials
    // S_fact = score!
    double *logScoreFactorial;

    // compute -log(p)
    double computeLogProbabiliy(const unsigned short rawScore, const unsigned int dbSeqLen);
};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
