#include "QueryScoreGlobal.h"

void QueryScoreGlobal::addScores (int* __restrict seqList, int seqListSize, unsigned short score){
    for (int i = 0; i < seqListSize; i++){
        const int seqId = seqList[i];
        scores[seqId] = sadd16(scores[seqId], score);
    }
    scoresSum += score * seqListSize;
    numMatches += seqListSize;
}

void QueryScoreGlobal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (thresholds_128, 0, scores_128_size * 2);
    resList->clear();
    scoresSum = 0;
    numMatches = 0;
}

