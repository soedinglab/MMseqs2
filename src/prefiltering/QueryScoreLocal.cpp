//
//  QueryScoreLocal.cpp
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#include "QueryScoreLocal.h"

QueryScoreLocal::~QueryScoreLocal(){
    free(scores_128);
    delete resList;
}


void QueryScoreLocal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
//    memset (thresholds_128, 0, scores_128_size * 2);
    localScoreResults.clear();
    scoresSum = 0;
    numMatches = 0;
}
