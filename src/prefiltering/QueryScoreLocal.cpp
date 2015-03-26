//
//  QueryScoreLocal.cpp
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//
#include <random>
#include "QueryScoreLocal.h"
#include "QueryScore.h"

QueryScoreLocal::QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr, double kmerMatchProb, size_t maxHitsPerQuery)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb)    // Call the QueryScore constructor
{
    this->scoreSizes = new unsigned int[SCORE_RANGE];
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
    this->maxHitsPerQuery = maxHitsPerQuery;
    this->logMatchProb = log(kmerMatchProb);
}

QueryScoreLocal::~QueryScoreLocal(){
    delete [] scoreSizes;
}

void QueryScoreLocal::reset() {
//    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    memset(scoreSizes, 0, SCORE_RANGE * sizeof(unsigned int));
}

void QueryScoreLocal::setPrefilteringThresholds() { }

std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    size_t foundHits = 0;
    size_t scoreThr = 0;
    for(scoreThr = SCORE_RANGE - 1; scoreThr > 0 ; scoreThr--){
        foundHits += scoreSizes[scoreThr];
        if(foundHits >= maxHitsPerQuery)
            break;
    }
    scoreThr = (scoreThr == 0 ) ? 0 : scoreThr - 1;
    simd_int thr = simdi16_set( scoreThr );
    simd_int zero =  simdi_setzero();
    const size_t lenght = scores_128_size / (SIMD_SHORT_SIZE); // *2for short
    simd_int* __restrict s   = scores_128;
    unsigned int elementCounter = 0;
    for (size_t pos = 0; pos < lenght; pos++ ){
        // look for entries above the threshold
        simd_int currElements = simdi_load(s + pos);
        // shift out the diagonal information and move the score to the first byte
        currElements = simdi16_srli(currElements, 8);
        simd_int cmp = simdi16_gt(currElements, thr); // currElement > thr
        const unsigned int cmp_set_bits = simdi8_movemask(cmp);
        // set zero so no memset is needed
        simdi_store(s + pos, zero);
        if (cmp_set_bits != 0){
            for(unsigned int i = 0; i < SIMD_SHORT_SIZE; i++){
                if( CHECK_BIT(cmp_set_bits,i*2)) {
                    unsigned short rawScore = simdi16_extract (currElements,  i);
                    hit_t * result = (resList + elementCounter);
                    result->seqId = pos * SIMD_SHORT_SIZE + i;
                    // log(100) * log(100) /(log(qL)*log(dbL))
                    double dbSeqLen = seqLens[result->seqId];
                    float mu = kmerMatchProb * dbSeqLen;
                    //result->zScore = (((float)rawScore) - mu )/ sqrt(mu);

                    double score = (double) rawScore;
                    // compute -log(p)
                    // S_fact = score!
                    double S_fact = std::min(std::numeric_limits<double>::max(),
                            sqrt(2 * M_PI * score) * pow(score / exp(1), score) * exp(1 / (12  * score)));
                    double mid_term = score * (logMatchProb + log(dbSeqLen));
                    double first_term = -(mu * score /(score + 1));
                    double prob = first_term + mid_term - log(S_fact);
                    result->zScore = -prob;
                    //std::cout << result->zScore << std::endl;

                    //result->zScore = (rawScore);
                    result->prefScore = rawScore;
                    //scoreSizes += rawScore;
                    elementCounter++;
                    if(elementCounter >= MAX_RES_LIST_LEN){
                        // because the memset will not be finished
                        memset (s + pos, 0, scores_128_size * 2);
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

void QueryScoreLocal::updateScoreBins() {
    unsigned char * data =  (unsigned char *) scores;
    for(size_t i = 0; i < dbSize; i++){
        scoreSizes[data[(i * 2) + 1]]++;
    }
}
