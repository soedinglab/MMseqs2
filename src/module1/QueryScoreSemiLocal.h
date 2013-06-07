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
    QueryScoreSemiLocal(int dbSize, unsigned short * seqLens, float prefThreshold, int k)
    : QueryScore(dbSize, seqLens, prefThreshold, k)    // Call the QueryScore constructor
    {
        this->lastScores = new LastScore[dbSize];
        memset (this->lastScores, 0, sizeof(LastScore) * dbSize);
    };

    ~QueryScoreSemiLocal();

     struct LastScore{ 
         unsigned short lastScore;
         unsigned short lastMatchPos;
     };
    
    void addScores (int* seqList, int seqListSize, unsigned short score);
    void reset();

private:
    // last scores for local matching
    LastScore * lastScores;
};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
