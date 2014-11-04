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
{
    localResultSize = 0;
    minKmerMatch = 69/4;
    
    gapOpenPenalty = minKmerMatch;
    gapExtendPenalty = 2;
    localResults = new LocalResult[MAX_LOCAL_RESULT_SIZE];
};

QueryScoreLocal::~QueryScoreLocal(){
    delete [] localResults;
}

void QueryScoreLocal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    localResultSize = 0;
}

void QueryScoreLocal::setPrefilteringThresholds() { }

bool QueryScoreLocal::compareLocalResult(LocalResult first, LocalResult second){
    return (first.seqId < second.seqId);
}


std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    size_t elementCounter = 0;
    std::sort(localResults, localResults + localResultSize, compareLocalResult);
    
    unsigned short scoreMax  = 0;  // current best score
    LocalResult * currentLocalResult = localResults + 0;
    seqIdCurr = currentLocalResult->seqId;
    seqIdPrev = currentLocalResult->seqId;
    for (unsigned int i = 1; i < localResultSize; i++) {
        scoreMax += currentLocalResult->score;   // add score of current kmer match
        currentLocalResult = localResults + i;
        seqIdCurr = currentLocalResult->seqId;
        // if new sequence occures or end of data write the result back
        if( seqIdCurr != seqIdPrev ||  i == (localResultSize - 1) ){
            // write result to list
            hit_t * result = (resList + elementCounter);
            result->seqId = seqIdPrev;
            result->zScore = scoreMax;
            result->prefScore = scoreMax;
            elementCounter++;
            if(elementCounter >= MAX_RES_LIST_LEN)
                break;
            seqIdPrev = seqIdCurr;
            // reset values
            scoreMax  = 0;  // current best score
        }
        
    }
    
    
    std::sort(resList, resList + elementCounter, compareHits);
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
