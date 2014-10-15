//
//  QueryScoreLocal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreLocal.h"
#include <algorithm>

QueryScoreLocal::QueryScoreLocal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
{
    localResultSize = 0;
    minKmerMatch = 69/4;

    gapOpenPenalty = minKmerMatch;
    gapExtendPenalty = 2;
    matchArray = new LocalMatch[dbSize];
};

QueryScoreLocal::~QueryScoreLocal(){
    delete [] matchArray;
    matchArray = NULL;
}

void QueryScoreLocal::reset() {
    memset (matchArray, 0, dbSize * sizeof(LocalMatch));
    scoresSum = 0;
    numMatches = 0;
    localResultSize = 0;
}

void QueryScoreLocal::setPrefilteringThresholds() { }

//bool QueryScoreLocal::compareLocalResult(LocalResult first, LocalResult second){
//    return (first.seqId < second.seqId) ||
//    (first.seqId == second.seqId && first.i < second.i) ||
//    (first.seqId == second.seqId && first.i == second.i && first.j < second.j);
////    (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i == second.i && first.j <second.j);
////           (first.seqId == second.seqId && first.diagonal < second.diagonal) ||
////           (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i < second.i) ||
////           (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i == second.i && first.j <second.j);
//}


bool QueryScoreLocal::compareLocalMatch(LocalMatch first, LocalMatch second){
    return (first.prefScore > second.prefScore);

}

std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    const size_t maxHit = 2000;

    std::partial_sort(matchArray,
                      matchArray + maxHit,
                      matchArray + dbSize,
                      compareLocalMatch
                      );
    size_t elementCounter = 0;
    for (size_t i = 0; i < maxHit; i++) {
        hit_t * currHit = resList + i;
        if(matchArray[i].prefScore > 0 ){
            currHit->prefScore = matchArray[i].prefScore;
            currHit->seqId = matchArray[i].seqId;
            currHit->zScore = matchArray[i].prefScore;
            elementCounter++;
        }
    }
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
