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
    return (first.seqId < second.seqId) ||
    (first.seqId == second.seqId && first.i < second.i) ||
    (first.seqId == second.seqId && first.i == second.i && first.j < second.j);
    //    (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i == second.i && first.j <second.j);
    //           (first.seqId == second.seqId && first.diagonal < second.diagonal) ||
    //           (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i < second.i) ||
    //           (first.seqId == second.seqId && first.diagonal == second.diagonal && first.i == second.i && first.j <second.j);
}


std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    size_t elementCounter = 0;
    std::sort(localResults, localResults + localResultSize, compareLocalResult);
    
    
    unsigned short scoreMax  = 0;  // current best score
    unsigned short scoreCurr = 0;  // current score
    unsigned int   seqIdCurr  = 0;  // current sequence id
    unsigned int   seqIdPrev  = INT_MAX;  // current prev sequence id
    int   diagPrev = INT_MAX; // diagonal of previous kmer match
    unsigned short jPrev = 0;  // position of previous kmer match in db sequence (“target”)
    // load first element
    LocalResult * currentLocalResult = localResults + 0;
    seqIdCurr = currentLocalResult->seqId;
    seqIdPrev = currentLocalResult->seqId;
    //std::cout << "Extract Result (" << currentLocalResultPos << "): " << std::endl;
    for (unsigned int i = 1; i < localResultSize; i++) {
        int diagCurr = currentLocalResult->diagonal;  // diagonal i - j of current kmer match
        unsigned short jCurr = currentLocalResult->j;  // position of current kmer match in db sequence (“target”)
        //        std::cout << "id "<< seqIdCurr  << " i: " << currentLocalResult->i << " jCurr: " << jCurr << " jPrev: " << jPrev
        //        << " diag: " << currentLocalResult->diagonal << " diagPrev: " << diagPrev
        //        << " Score: "<< currentLocalResult->score << " scoreCurr: " << scoreCurr << " scoreMax: " << scoreMax;
        //        if(diagPrev != INT_MAX && std::abs(diagCurr - diagPrev) > 1000)
        //            std::cout <<
        
        
        if (diagPrev != diagCurr) {
            // If last match of previous diagonal is not compatible with first kmer match of new diagonal
            // make sure the scores of these two diagonals cannot be accumulated together
            scoreCurr = (jPrev >= jCurr? 0.0 : scoreCurr); // effect on performance 0.4721
            // Add score of first kmer match of diagonal
            //int scoreDiagonalSwitch = (diagPrev == INT_MAX) ? 0.0 : minKmerMatch; // not a good a feature 0.4563526602096326
            // both combined = 0.46486002388218123
            //if(diagPrev != INT_MAX && std::abs(diagCurr - diagPrev) > 1000)
            //    std::cout << "Wrong babe " << jCurr << " " << jPrev << " "
            //              << diagCurr << " " << diagPrev << " " << std::abs(diagCurr - diagPrev) << std::endl;
            // Subtract gap costs for gap
            //scoreDiagonalSwitch -= (diagPrev == INT_MAX) ?  0.0 :
            //                        this->gapOpenPenalty + this->gapExtendPenalty * std::abs(diagCurr - diagPrev);
            //scoreCurr = (scoreCurr > scoreCurr + scoreDiagonalSwitch) ? scoreCurr : scoreCurr + scoreDiagonalSwitch;
            
        }
        
        scoreCurr += currentLocalResult->score;   // add score of current kmer match
        scoreMax = std::max(scoreCurr, scoreMax);
        diagPrev = diagCurr;
        jPrev = jCurr;
        // load next result
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
            scoreCurr = 0;  // current score
            jPrev = 0;  // position of previous kmer match in db sequence (“target”)
            diagPrev = INT_MAX; // diagonal of previous kmer match
        }
        
    }
    
    
    std::sort(resList, resList + elementCounter, compareHits);
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
