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


class QueryScoreLocal : public QueryScore {
    
public:
    QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr, double kmerMatchProb, size_t maxHitsPerQuery);
    
    ~QueryScoreLocal();
    
    std::pair<hit_t *, size_t> getResult (int querySeqLen, unsigned int identityId);
    
    void reset();
    
    bool checkForOverflowAndResizeArray();
    
    // NOT needed for Local scoring
    void setPrefilteringThresholds();

    void updateScoreBins();

private:

    size_t maxHitsPerQuery;
    double logMatchProb;
};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
