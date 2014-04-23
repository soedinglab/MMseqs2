#include "QueryScoreGlobal.h"

void QueryScoreGlobal::addScores (int* seqList, int seqListSize, unsigned short score){
    for (int i = 0; i < seqListSize; i++){
        scores[seqList[i]] = sadd16(scores[seqList[i]], score);
    }
    scoresSum += score * seqListSize;
    numMatches += seqListSize;
}

void QueryScoreGlobal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (thresholds_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
}

