//
//  QueryScoreGlobal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreGlobal.h"

inline unsigned short sadd16(unsigned short a, unsigned short b)
{
	unsigned int s = (unsigned int)a+b;
	return -(s>>16) | (unsigned short)s;
}
    
void QueryScoreGlobal::addScores (int* seqList, int seqListSize, unsigned short score){
    for (int i = 0; i < seqListSize; i++){
        const int seqId=seqList[i];
        scores[seqId]=sadd16(scores[seqId],score);
    }
}


void QueryScoreGlobal::reset() {
    memset (scores_128, 0, (dbSize/8 + 1) * 16);
    
    resList->clear();
    hitList->clear();
}