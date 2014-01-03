//
//  QueryScoreLocal.h
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef QUERYSCORESEMILOCAL_H
#define QUERYSCORESEMILOCAL_H
#include "QueryScore.h"

class QueryScoreSemiLocal : public QueryScore {
    
public:
    QueryScoreSemiLocal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
    : QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
    {
        this->lastScores = new LastScore[dbSize];
        memset (this->lastScores, 0, sizeof(LastScore) * dbSize);
    };

    ~QueryScoreSemiLocal();

     struct LastScore{ 
         unsigned short lastScore;
         short lastMatchPos;
     };
    
    void addScores (int* seqList, int seqListSize, unsigned short score);

    void reset();

private:
    // last scores for local matching
    LastScore * lastScores;
};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
