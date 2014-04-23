//
//  QueryScoreLocal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreSemiLocal.h"

QueryScoreSemiLocal::~QueryScoreSemiLocal(){
    free(scores_128);
    delete resList;
    delete[] lastScores;
}

void QueryScoreSemiLocal::addScores(int *seqList, int seqListSize, unsigned short score){
    const int currentMatchPos=0;
    const int f = 1;
    for (int i = 0; i < seqListSize; i++){
        const int seqId=seqList[i];
        LastScore lastScore = this->lastScores[seqId];
        short scoreDrop = f*(currentMatchPos - lastScore.lastMatchPos);
        lastScore.lastScore          = std::max( 0,lastScore.lastScore - scoreDrop ); // saturated subtract
        lastScore.lastScore          = sadd16_signed(lastScore.lastScore , score);
        scores[seqId]                = std::max(scores[seqId],lastScore.lastScore);
        lastScore.lastMatchPos       = currentMatchPos;
    }
}

void QueryScoreSemiLocal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (this->lastScores, 0, sizeof(LastScore) * dbSize);    
}
