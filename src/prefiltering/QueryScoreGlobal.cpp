//
//  QueryScoreGlobal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreGlobal.h"

void QueryScoreGlobal::addScores (int* seqList, int seqListSize, short score){
    for (int i = 0; i < seqListSize; i++){
        scores[seqList[i]] = sadd16_signed(scores[seqList[i]], score);
    }
    scoresSum += score * seqListSize;
    numMatches += seqListSize;
}

void QueryScoreGlobal::addScoresRevSeq(int* seqList, int seqListSize, short score){
    for (int i = 0; i < seqListSize; i++){
        scores[seqList[i]] = ssub16_signed(scores[seqList[i]], score);
    }
    scoresSumRevSeq +=  score * seqListSize;
    numMatchesRevSeq += seqListSize;
}

void QueryScoreGlobal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (thresholds_128, 0, scores_128_size * 2);
    resList->clear();
    scoresSum = 0;
    scoresSumRevSeq = 0;
    numMatches = 0;
    numMatchesRevSeq = 0;
}

