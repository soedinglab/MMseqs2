#include "QueryScoreGlobal.h"


QueryScoreGlobal::QueryScoreGlobal(int dbSize, unsigned short * dbSeqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr)
: QueryScore(dbSize, dbSeqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
{
    
    thresholds_128 = (__m128i*) Util::mem_align(16, scores_128_size * 2);
    thresholds = (unsigned short * ) thresholds_128;
    
    memset (thresholds_128, 0, scores_128_size * 2);
    
    // initialize sequence lenghts with each seqLens[i] = L_i - k + 1
    this->seqLens = new float[scores_128_size];
    memset (seqLens, 0, scores_128_size * 4);
    
    for (int i = 0; i < dbSize; i++){
        if (dbSeqLens[i] > (k - 1))
            this->seqLens[i] = (float) (dbSeqLens[i] - k + 1);
        else
            this->seqLens[i] = 1.0f;
    }
    
    this->seqLenSum = 0.0f;
    for (int i = 0; i < dbSize; i++)
        this->seqLenSum += this->seqLens[i];
    
    // initialize the points where a score threshold should be recalculated
    std::list<int> steps_list;
    float seqLen = this->seqLens[0];
    steps_list.push_back(0);
    // check the sequence length and decide if it changed enough to recalculate the score threshold here
    // minimum step length is 8 (one __m128 register)
    for (int i = 0; i < scores_128_size; i += 8){
        if (this->seqLens[i]/seqLen < 0.9){
            steps_list.push_back(i);
            seqLen = this->seqLens[i];
        }
    }
    steps_list.push_back(scores_128_size);
    
    nsteps = steps_list.size();
    steps = new int[nsteps];
    for (int i = 0; i < nsteps; i++){
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
    
    __m128i* __restrict thr = thresholds_128;
    __m128i v = _mm_setzero_si128();
    
    for (int i = 0; i < nsteps - 1; i++){
        seqLen = seqLens[steps[i]];
        mean = s_per_pos * seqLen;
        stddev = sqrt(seqLen * s_per_pos * s_per_match);
        threshold = zscore_thr * stddev + mean;
        
        // saturated conversion of float threshold into short threshold
        threshold_ushrt = (unsigned short) std::min(ushrt_max, (unsigned int) threshold);
        
        // store short treshold in a vector
        v =  _mm_set1_epi16(threshold_ushrt);
        
        // set all thresolds to the current value until calculation of the next value is necessary
        for (int j = steps[i]; j < steps[i+1]; j += 8){
            *thr = v;
            thr++;
        }
    }
}

std::pair<hit_t *, size_t> QueryScoreGlobal::getResult (int querySeqLen, unsigned int identityId){
    
    size_t elementCounter = 0;
    const __m128i zero = _mm_setzero_si128();
    
    __m128i* __restrict p   = scores_128;
    __m128i* __restrict thr = thresholds_128;
    
    
    // check if there is the identity of the query sequence in the database
    // the identity should be included in the results
    if (identityId != UINT_MAX){
        elementCounter++;
    }
    
    // go through each vector
    const size_t lenght = scores_128_size/8;
    for (size_t pos = 0; pos < lenght; pos++ ){
        
        // look for entries above the threshold
        __m128i cmp = _mm_subs_epu16(*p, *thr);
        cmp = _mm_cmpeq_epi16(cmp, zero);
        const unsigned int cmp_set_bits = _mm_movemask_epi8(cmp);
        
        // here are some sequences above the prefiltering threshold
        if (cmp_set_bits != 0xffff){
            // and search for highest
            for(int i = 0; i < 8; i++){
                
                if(!CHECK_BIT(cmp_set_bits,i*2) && (pos * 8 + i) != identityId){
                    const float zscore = getZscore(pos*8+i);
                    hit_t * result = (resList + elementCounter);
                    result->seqId = pos * 8 + i;
                    result->zScore = zscore;
                    result->prefScore = scores[pos*8+i];
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


