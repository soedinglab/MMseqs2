#include "QueryScore.h"
#include "../commons/Util.h"
#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define _mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#define _mm_extract_epi64(x, imm) _mm_cvtsi128_si64(_mm_srli_si128((x), 8 * (imm)))

QueryScore::QueryScore (int dbSize, unsigned short * dbSeqLens, int k, short kmerThr, float kmerMatchProb, float zscoreThr){

    this->dbSize = dbSize;
    this->kmerMatchProb = kmerMatchProb;
    this->kmerThr = kmerThr;
    this->zscore_thr = zscoreThr;

    this->scores_128_size = (dbSize + 7)/8 * 8;
    // 8 DB short int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (__m128i*) Util::mem_align(16, scores_128_size * 2);
    scores = (unsigned short * ) scores_128;

    // set scores to zero
    memset (scores_128, 0, scores_128_size * 2);

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

    this->resList = (hit_t *) Util::mem_align(16, MAX_RES_LIST_LEN * sizeof(hit_t) );

    scoresSum = 0;

    numMatches = 0;

    counter = 0;

    s_per_match = 0.0f;

    s_per_pos = 0.0f;
}

QueryScore::~QueryScore (){
    free(scores_128);
    free(thresholds_128);
    delete[] seqLens;
    delete[] steps;
    free(resList);
}

bool QueryScore::compareHits(hit_t first, hit_t second){
    return (first.zScore > second.zScore) ? true : false;
}

void QueryScore::setPrefilteringThresholds(){

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

    __m128i* thr = thresholds_128;
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

float QueryScore::getZscore(int seqId){
    return ( (float)scores[seqId] - s_per_pos * seqLens[seqId] ) / sqrt(s_per_pos * seqLens[seqId] * s_per_match);
}

std::pair<hit_t *, size_t> QueryScore::getResult (int querySeqLen, unsigned int identityId){
    size_t elementCounter = 0;
    const __m128i zero = _mm_setzero_si128(); 

    __m128i* p = scores_128;
    __m128i* thr = thresholds_128;

    __m128i cmp;
    
    // check if there is the identity of the query sequence in the database
    // the identity should be included in the results
    if (identityId != UINT_MAX){
        elementCounter++;
    }

    // go through each vector
    for (int pos = 0; pos < scores_128_size/8; pos++ ){

        // look for entries above the threshold
        cmp = _mm_subs_epu16(*p, *thr);
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

    return std::make_pair<hit_t *, size_t>(resList, elementCounter);
}

short QueryScore::sse2_extract_epi16(__m128i v, int pos) {
    switch(pos){
        case 0: return _mm_extract_epi16(v, 0);
        case 1: return _mm_extract_epi16(v, 1);
        case 2: return _mm_extract_epi16(v, 2);
        case 3: return _mm_extract_epi16(v, 3);
        case 4: return _mm_extract_epi16(v, 4);
        case 5: return _mm_extract_epi16(v, 5);
        case 6: return _mm_extract_epi16(v, 6);
        case 7: return _mm_extract_epi16(v, 7);
    }
    Debug(Debug::ERROR) << "Fatal error in QueryScore: position in the vector is not in the legal range (pos = " << pos << ")\n";
    exit(1);
    // never executed
    return 0;
}

void QueryScore::printVector(__m128i v){
    for (int i = 0; i < 8; i++)
        std::cout << (unsigned short) sse2_extract_epi16(v, i) << " ";
    std::cout << "\n";
}

void QueryScore::printScores(){
    std::cout << "Scores:\n";
    for (int i = 0; i < dbSize; i++)
        std::cout << scores[i] << "\n";
}
