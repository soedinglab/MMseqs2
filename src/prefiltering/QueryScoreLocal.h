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

    // compute look up table based on stirling approximation
    static void computeFactorial(double *output, const size_t range) {
        output[0] = log(1.0);
        for(size_t score = 1; score < range; score++){
            const double scoreDbl = static_cast<double>(score);

            const double S_fact = std::min(std::numeric_limits<double>::max(),
                                           sqrt(2 * M_PI * scoreDbl) * pow(scoreDbl / exp(1), scoreDbl) * exp(1 / (12  * scoreDbl)));
            output[score] = log(S_fact);
        }
    }

    // compute -log(p)
    static inline double computeLogProbability(const unsigned short rawScore, const unsigned int dbSeqLen,
                                               const double kmerMatchProb, const double logMatchProb,
                                               const double logScoreFactorial) {
        const double score = static_cast<double>(rawScore);
        const double dbSeqLenDbl = static_cast<double>(dbSeqLen);
        const double mu = kmerMatchProb * dbSeqLenDbl;
        const double mid_term = score * (logMatchProb + log(dbSeqLenDbl));
        const double first_term = -(mu * score /(score + 1));
        return first_term + mid_term - logScoreFactorial;
    }
private:

    //log match prob (mu) of poisson distribution
    double logMatchProb;

    //pre computed score factorials
    // S_fact = score!
    double *logScoreFactorial;

};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
