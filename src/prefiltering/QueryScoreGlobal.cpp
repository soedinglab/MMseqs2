#include "QueryScoreGlobal.h"



void QueryScoreGlobal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (thresholds_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
}

