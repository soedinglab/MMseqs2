//
//  QueryScoreLocal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreLocal.h"

QueryScoreLocal::QueryScoreLocal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
{    };

QueryScoreLocal::~QueryScoreLocal(){

}

void QueryScoreLocal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    localResult.clear();
}

void QueryScoreLocal::setPrefilteringThresholds() { }

std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    size_t elementCounter = 0;
    for(std::unordered_map<unsigned int, unsigned short>::iterator it = localResult.begin();
        it != localResult.end(); it++){
        hit_t * result = (resList + elementCounter);
        result->seqId = it->first;
        result->zScore = it->second;
        result->prefScore = it->second;
        elementCounter++;
        if(elementCounter >= MAX_RES_LIST_LEN)
            break;
        
    }
    std::sort(resList, resList + elementCounter, compareHits);
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
