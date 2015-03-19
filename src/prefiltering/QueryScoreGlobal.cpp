#include <random>
#include "QueryScoreGlobal.h"
#include "simd.h"

QueryScoreGlobal::QueryScoreGlobal(size_t dbSize, unsigned int *dbSeqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
: QueryScore(dbSize, dbSeqLens, k, kmerThr, kmerMatchProb)    // Call the QueryScore constructor
{
    
    thresholds_128 = (simd_int*) mem_align(ALIGN_INT, scores_128_size * sizeof(unsigned short));
    thresholds = (unsigned short * ) thresholds_128;
    memset (thresholds_128, 0, scores_128_size * sizeof(unsigned short));
    this->zscore_thr = zscoreThr;

    // initialize the points where a score threshold should be recalculated
    std::list<int> steps_list;
    float seqLen = this->seqLens[0];
    steps_list.push_back(0);
    // check the sequence length and decide if it changed enough to recalculate the score threshold here
    // minimum step length is 8 (one __m128 register)
    for (size_t i = 0; i < scores_128_size; i += 8){
        if (this->seqLens[i]/seqLen < 0.9){
            steps_list.push_back(i);
            seqLen = this->seqLens[i];
        }
    }
    steps_list.push_back(scores_128_size);
    
    nsteps = steps_list.size();
    steps = new int[nsteps];
    for (size_t i = 0; i < nsteps; i++){
        steps[i] = steps_list.front();
        steps_list.pop_front();
    }
    
    counter = 0;
    s_per_match = 0.0f;
    s_per_pos = 0.0f;
}

QueryScoreGlobal::~QueryScoreGlobal(){
    free(thresholds_128);
    delete[] seqLens;
    delete[] steps;
}


void QueryScoreGlobal::reset() {
    memset (scores_128, 0, scores_128_size * 2);
    memset (thresholds_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
}

void QueryScoreGlobal::setPrefilteringThresholds(){
    /* adding 0.000001 to some values should prevent nan values in case of untypical/less meaningful input parameters */
    // pseudo-count sum of sequence lengths
    // 100 000 * 350
    float seqLenSum_pc = 35000000.0;
    // pseudo-number of k-mer matches
    float numMatches_pc = seqLenSum_pc *  kmerMatchProb + 0.000001f;
    
    // pseudo-sum of k-mer scores
    float matchScoresSum_pc = numMatches_pc * (float) (kmerThr + 8);
    this->s_per_match = ((float)scoresSum + matchScoresSum_pc)/((float)numMatches + numMatches_pc) + 0.000001f;
    
    float scoresSum_pc = seqLenSum_pc * kmerMatchProb * s_per_match;
    this->s_per_pos = ((float)scoresSum + scoresSum_pc)/(seqLenSum + seqLenSum_pc) + 0.000001f;
    
    float seqLen;
    float mean;
    float stddev;
    float threshold;
    unsigned short threshold_ushrt;
    unsigned int ushrt_max = USHRT_MAX;
    
    simd_int* __restrict thr = thresholds_128;
    
    for (int i = 0; i < nsteps - 1; i++){
        seqLen = seqLens[steps[i]];
        mean = s_per_pos * seqLen;
        stddev = sqrt(seqLen * s_per_pos * s_per_match);
        threshold = zscore_thr * stddev + mean;
        
        // saturated conversion of float threshold into short threshold
        threshold_ushrt = (unsigned short) std::min(ushrt_max, (unsigned int) threshold);
        
        // store short treshold in a vector
        simd_int v = simdi16_set(threshold_ushrt);
        
        // set all thresolds to the current value until calculation of the next value is necessary
        for (int j = steps[i]; j < steps[i+1]; j += SIMD_SHORT_SIZE){
            *thr = v;
            thr++;
        }
    }
}

std::pair<hit_t *, size_t> QueryScoreGlobal::getResult (int querySeqLen, unsigned int identityId){
    
    size_t elementCounter = 0;
    const simd_int zero = simdi32_set(0);
    
    simd_int* __restrict p   = scores_128;
    simd_int* __restrict thr = thresholds_128;
    
    
    // check if there is the identity of the query sequence in the database
    // the identity should be included in the results
    if (identityId != UINT_MAX){
        elementCounter++;
    }
    
    // go through each vector
    const size_t lenght = scores_128_size/ (SIMD_SHORT_SIZE); // *2for short
    for (size_t pos = 0; pos < lenght; pos++ ){
        
        // look for entries above the threshold
        simd_int cmp = simdui16_subs(*p, *thr);
        cmp = simdi16_eq(cmp, zero);
        const unsigned int cmp_set_bits = simdi8_movemask(cmp);
        
        // here are some sequences above the prefiltering threshold
#ifdef AVX2
        if (cmp_set_bits != 0xffffffff){
#else
        if (cmp_set_bits != 0xffff){
#endif
            // and search for highest
            for(unsigned int i = 0; i < SIMD_SHORT_SIZE; i++){
                
                if(!CHECK_BIT(cmp_set_bits,i*2) && (pos * SIMD_SHORT_SIZE + i) != identityId){
                    const float zscore = getZscore(pos * SIMD_SHORT_SIZE + i);
                    hit_t * result = (resList + elementCounter);
                    result->seqId = pos * SIMD_SHORT_SIZE + i;
                    result->zScore = zscore;
                    result->prefScore = scores[pos * SIMD_SHORT_SIZE + i];
                    elementCounter++;
                    if(elementCounter >= MAX_RES_LIST_LEN)
                        goto OuterLoop;
                }
            }
        }
        p++;
        thr++;
    }
OuterLoop:
    // include the identity in results if its there
    if (identityId != UINT_MAX){
        const float zscore = getZscore(identityId);
        hit_t * result = (resList + 0);
        result->seqId = identityId;
        result->zScore = zscore;
        result->prefScore = scores[identityId];
        std::sort(resList + 1, resList + elementCounter, compareHits);
    }
    else
        std::sort(resList, resList + elementCounter, compareHits);
    
    return std::make_pair(resList, elementCounter);
}

float QueryScoreGlobal::getZscore(int seqId){
    return ( (float)scores[seqId] - s_per_pos * seqLens[seqId] ) / sqrt(s_per_pos * seqLens[seqId] * s_per_match);
}


