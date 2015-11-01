//
//  QueryScoreLocal.cpp
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2015 Martin Steinegger. All rights reserved.
//
#include "QueryScoreLocal.h"
#include "QueryScore.h"
#include <random>

QueryScoreLocal::QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr, double kmerMatchProb)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb)    // Call the QueryScore constructor
{
    this->scoreSizes = new unsigned int[SCORE_RANGE];
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
    this->logMatchProb = log(kmerMatchProb);
    this->logScoreFactorial = new double[SCORE_RANGE];
    // pre compute the factorial of all possible scores
    computeFactorial(logScoreFactorial, SCORE_RANGE);
}

QueryScoreLocal::~QueryScoreLocal(){
    delete [] scoreSizes;
    delete [] logScoreFactorial;
}

void QueryScoreLocal::reset() {
//    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
}

std::pair<hit_t *, size_t> QueryScoreLocal::getResult(const unsigned int querySeqLen, unsigned int scoreThr, const unsigned int identityId) {
    // scoreThr get all great than scoreThr
    // adjust score threshold because we use > for check later
    scoreThr = (scoreThr == 0 ) ? 0 : scoreThr - 1;
    if(scoreThr > 127)
        scoreThr = 0;
    simd_int thr = simdi8_set( 127 - scoreThr );
    simd_int zero =  simdi_setzero();
    const size_t lenght = scores_128_size / (SIMD_SHORT_SIZE); // *2for short
    simd_int* __restrict s   = scores_128;
    unsigned int elementCounter = 0;

    // check if there is the identity of the query sequence in the database
    // the identity should be included in the results
    if (identityId != UINT_MAX){
        hit_t * result = (resList + 0);
        const unsigned char * data =  (unsigned char *) scores;
        const unsigned int seqIndex = identityId * 2;
        unsigned char rawScore  = data[seqIndex + 1];
        rawScore = (rawScore == 0) ? 1 : rawScore;
        result->seqId = identityId;
        result->prefScore = rawScore;
        result->pScore = -computeLogProbability(rawScore, seqLens[identityId],
                                                kmerMatchProb, logMatchProb, logScoreFactorial[rawScore]);
        elementCounter++;
    }

    for (size_t pos = 0; pos < lenght; pos++ ){
        // look for entries above the threshold
        simd_int currElements = simdi_load(s + pos);
        // make a saturated add (significant bit has to set so that move mask sets a bit)
        simd_int cmp = simdui8_adds(currElements, thr);
        const unsigned int cmp_set_bits = simdi8_movemask(cmp);
        // set zero so no memset is needed
        simdi_store(s + pos, zero);
#ifdef AVX2
        if (UNLIKELY((cmp_set_bits & 0xAAAAAAAA) != 0)){
#else
        if (UNLIKELY((cmp_set_bits & 0xAAAA) != 0)){
#endif
            for(size_t i = 0; i < SIMD_SHORT_SIZE; i++){
                if( CHECK_BIT(cmp_set_bits, i * 2 + 1) && pos * SIMD_SHORT_SIZE + i != identityId) {
                    unsigned short rawScore = simdi16_extract (currElements,  i) >> 8;
                    hit_t * result = (resList + elementCounter);
                    result->seqId = pos * SIMD_SHORT_SIZE + i;
                    //result->pScore = (((float)rawScore) - mu )/ sqrt(mu);
                    result->prefScore = rawScore;
                    result->pScore = -computeLogProbability(rawScore, seqLens[result->seqId],
                                                            kmerMatchProb, logMatchProb, logScoreFactorial[rawScore]);
                    //std::cout << result->pScore << std::endl;

                    //result->pScore = (rawScore);
                    //scoreSizes += rawScore;
                    elementCounter++;
                    if(elementCounter >= MAX_RES_LIST_LEN){
                        // because the memset will not be finished
                        for (size_t rest_pos = pos; rest_pos < lenght; rest_pos++ ) {
                            simdi_store(s + rest_pos, zero);
                        }
                        goto OuterLoop;
                    }
                }
            }
        }
    }
    OuterLoop:

    // sort hits by score
    // include the identity in results if its there
    if(elementCounter > 1){
        if (identityId != UINT_MAX){
            std::sort(resList + 1, resList + elementCounter, compareHits);
        }
        else{
            std::sort(resList, resList + elementCounter, compareHits);
        }
    }
    return std::make_pair(this->resList, elementCounter);
}

void QueryScoreLocal::updateScoreBins() {
    unsigned char * data =  (unsigned char *) scores;
    for(size_t i = 0; i < dbSize; i++){
        scoreSizes[data[(i * 2) + 1]]++;
    }
}

