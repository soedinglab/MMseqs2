//
//  QueryScoreLocal.cpp
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//
#include "QueryScoreLocal.h"
#include "QueryScore.h"

QueryScoreLocal::QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr,
                                 double kmerMatchProb, float zscoreThr, size_t binSize)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
{
    localResultSize = 0;
    localResults = new unsigned int[MAX_LOCAL_RESULT_SIZE];
    this->seqsLens = seqLens;
}

QueryScoreLocal::~QueryScoreLocal(){
    delete [] localResults;
}

void QueryScoreLocal::reset() {
//    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    localResultSize = 0;
}

void QueryScoreLocal::setPrefilteringThresholds() { }

std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    simd_int thr = simdi16_set(0);
    simd_int zero =  simdi_setzero();
    const size_t lenght = scores_128_size / (SIMD_SHORT_SIZE); // *2for short
    simd_int* __restrict s   = scores_128;
    const float log_qL = log(querySeqLen);
    unsigned int elementCounter = 0;
    for (size_t pos = 0; pos < lenght; pos++ ){
        // look for entries above the threshold
        simd_int currElements = simdi_load(s + pos);
        // shift out the diagonal information and move the score to the first byte
        currElements = simdi16_srli(currElements, 8);
        simd_int cmp = simdi16_gt(currElements, thr);
        const unsigned int cmp_set_bits = simdi8_movemask(cmp);
        // set zero so no memset is needed
        simdi_store(s + pos, zero);
        if (cmp_set_bits != 0){
            for(unsigned int i = 0; i < SIMD_SHORT_SIZE; i++){
                if( CHECK_BIT(cmp_set_bits,i*2)) {
                    unsigned short element = simdi16_extract (currElements,  i);
                    hit_t * result = (resList + elementCounter);
                    result->seqId = pos * SIMD_SHORT_SIZE + i;
                    // log(100) * log(100) /(log(qL)*log(dbL))
                    result->zScore = (element) * 21.20/(log_qL*log(seqsLens[result->seqId]));
                    //result->zScore = (element);
                    result->prefScore = element;
                    localResultSize += element;
                    elementCounter++;
                    if(elementCounter >= MAX_RES_LIST_LEN){
                        // because the memset will not be finished
                        memset (scores_128, 0, scores_128_size * 2);
                        goto OuterLoop;
                    }
                }
            }
        }
    }
    OuterLoop:
    // sort hits by score
    //TODO maybe sort of hits not needed if SW is included
    std::sort(resList, resList + elementCounter, compareHits);
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}
