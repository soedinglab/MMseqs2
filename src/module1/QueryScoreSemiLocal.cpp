//
//  QueryScoreLocal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreSemiLocal.h"

inline unsigned short sadd16(unsigned short a, unsigned short b)
{
	unsigned int s = (unsigned int)a+b;
	return -(s>>16) | (unsigned short)s;
}

void QueryScoreSemiLocal::addScores(int *seqList, int seqListSize, unsigned short score){
    const int currentMatchPos=0;
    const int f = 1;
    for (int i = 0; i < seqListSize; i++){
        const int seqId=seqList[i];
        LastScore lastScore = this->lastScores[seqId];
        unsigned short scoreDrop = f*(currentMatchPos - lastScore.lastMatchPos);
        lastScore.lastScore          = std::max( 0,lastScore.lastScore - scoreDrop ); // saturated subtract
        lastScore.lastScore          = sadd16(lastScore.lastScore , score);
        scores[seqId]                = std::max(scores[seqId],lastScore.lastScore);
        lastScore.lastMatchPos       = currentMatchPos;
    }
}


void QueryScoreSemiLocal::reset() {
    memset (scores_128, 0, (dbSize/8 + 1) * 16);
    memset (this->lastScores, 0, sizeof(LastScore) * dbSize);
    
    resList->clear();
    hitList->clear();
}