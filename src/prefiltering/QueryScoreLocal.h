//
//  QueryScoreLocal.h
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef QUERYSCORELOCAL_H
#define QUERYSCORELOCAL_H
#include "QueryScore.h"
#include <vector>
#include <map>

struct LocalResult{
    const unsigned short i;
    const unsigned short j;
    const unsigned short score;
    LocalResult( unsigned short i,
                unsigned short j,
                unsigned short score) :
    i(i), j(j), score(score) {};
};


class QueryScoreLocal : public QueryScore {
    
public:
    QueryScoreLocal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
    : QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
    {    };

    ~QueryScoreLocal();


    std::map<unsigned int,std::vector<LocalResult *> > localScoreResults;


    
    // add k-mer match score for all DB sequences from the list
    inline void addScores (unsigned int* __restrict seqList, int seqListSize,
                                unsigned short i, unsigned short j,
                                unsigned short score){
        for (int i = 0; i < seqListSize; i++){
            const unsigned int seqId = seqList[i];
            const int diagonal = i - j + 255; // contains now the last diagonal
            if(diagonal == scores[seqId]){ // two elements found (Not true for 0)
                // add score to
                localScoreResults[seqId].push_back(new LocalResult(i, j, score));
            }
            scores[seqId] = diagonal;
        }
        
        scoresSum += score * seqListSize;
        numMatches += seqListSize;
    };

    void reset();

};
#endif /* defined(QUERYSCORESEMILOCAL_H) */
