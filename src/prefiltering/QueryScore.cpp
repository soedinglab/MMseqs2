#include "QueryScore.h"

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define _mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#define _mm_extract_epi64(x, imm) _mm_cvtsi128_si64(_mm_srli_si128((x), 8 * (imm)))


void *memalign(size_t boundary, size_t size)
{
    void *pointer;
    if (posix_memalign(&pointer,boundary,size) != 0)
    {
        std::cerr<<"Error in memalign: Could not allocate memory by memalign. Please report this bug to developers\n";
        exit(3);
    }
    return pointer;
}

QueryScore::QueryScore (int dbSize, unsigned short * dbSeqLens, int k, short kmerThr, float kmerMatchProb, float zscoreThr){

    this->dbSize = dbSize;
    this->kmerMatchProb = kmerMatchProb;
    this->kmerThr = kmerThr;
    this->zscore_thr = zscoreThr;

    this->scores_128_size = (dbSize + 7)/8 * 8;
    // 8 DB short int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (__m128i*) memalign(16, scores_128_size * 2);
    scores = (short * ) scores_128;

    // set scores to zero
    memset (scores_128, 0, scores_128_size * 2);

    thresholds_128 = (__m128i*) memalign(16, scores_128_size * 2);
    thresholds = (short * ) thresholds_128;

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

    nsteps = steps_list.size();
    steps = new int[nsteps];
    for (int i = 0; i < nsteps; i++){
        steps[i] = steps_list.front();
        steps_list.pop_front();
    }

    this->resList = new std::list<hit_t>();

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
    delete resList;
}

bool QueryScore::compareHits(hit_t first, hit_t second){
   if (first.prefScore > second.prefScore)
       return true;
   return false;
}

void QueryScore::setPrefilteringThresholds(){

    // pseudo-count sum of sequence lengths
    // 100 000 * 350
    float seqLenSum_pc = 35000000.0;
    // pseudo-number of k-mer matches
    float numMatches_pc = seqLenSum_pc *  kmerMatchProb;
   
    // pseudo-sum of k-mer scores
    float matchScoresSum_pc = numMatches_pc * (float) (kmerThr + 8);
    this->s_per_match = ((float)scoresSum + matchScoresSum_pc)/((float)numMatches + numMatches_pc);

    float scoresSum_pc = seqLenSum_pc * kmerMatchProb * s_per_match;
    this->s_per_pos = ((float)scoresSum + scoresSum_pc)/(seqLenSum + seqLenSum_pc);
    
    float seqLen;
    float mean;
    float stddev;
    float threshold;
    short short_threshold = SHRT_MAX;

    __m128i* thr = thresholds_128;

    for (int i = 0; i < nsteps - 1; i++){
        seqLen = seqLens[steps[i]];
        mean = s_per_pos * seqLen;
        stddev = sqrt(seqLen * s_per_pos * s_per_match);
        threshold = zscore_thr * stddev + mean;

        // saturated conversion of float threshold into short threshold
        // write short threshold into the __m128i vector
        int t = (int) threshold;
        // write 2x int threshold into 64 bit integer
        int64_t thr_int64 =  ((int64_t)t) << 32 | t;
        // write 64 bit integer threshold into __m64 data type
        __m64 i_m64 = _mm_cvtsi64_m64 (thr_int64);
        // saturated conversion of 2x 32 bit into 4x 16 bit 
        __m64 i_short_v = _mm_packs_pi32(i_m64, i_m64);
        // store short treshold in a vector
        __m128i v =  _mm_set1_epi64(i_short_v);

        // set all thresolds to the current value until calculation of the next value is necessary
        for (int j = steps[i]; j < steps[i+1]; j += 8){
            *thr = v;
            thr++;
        }
    }
    // set the rest of the thresholds from the last step position
    for (int j = steps[nsteps-1]; j < scores_128_size; j += 8){
        *thr = _mm_set1_epi16(short_threshold);
        thr++;
    }

}


void QueryScore::setPrefilteringThresholdsRevSeq(){
    this->s_per_match = (float)scoresSumRevSeq / (float)numMatchesRevSeq;
    this->matches_per_pos = (float)numMatchesRevSeq / seqLenSum;
    float factor = zscore_thr * s_per_match * sqrt( 2.0 * matches_per_pos);
    float threshold;
    short short_threshold = SHRT_MAX;

    __m128i* thr = thresholds_128;
    for (int i = 0; i < nsteps - 1; i++){
        threshold = factor * sqrt(seqLens[steps[i]]);

        // saturated conversion of float threshold into short threshold
        // write short threshold into the __m128i vector
        int t = (int) threshold;
        // write 2x int threshold into 64 bit integer
        int64_t thr_int64 =  ((int64_t)t) << 32 | t;
        // write 64 bit integer threshold into __m64 data type
        __m64 i_m64 = _mm_cvtsi64_m64 (thr_int64);
        // saturated conversion of 2x 32 bit into 4x 16 bit 
        __m64 i_short_v = _mm_packs_pi32(i_m64, i_m64);
        // store short treshold in a vector
        __m128i v =  _mm_set1_epi64(i_short_v);

        // set all thresolds to the current value until calculation of the next value is necessary
        for (int j = steps[i]; j < steps[i+1]; j += 8){
            *thr = v;
            thr++;
        }
    }
    // set the rest of the thresholds from the last step position
    for (int j = steps[nsteps-1]; j < scores_128_size; j += 8){
        *thr = _mm_set1_epi16(short_threshold);
        thr++;
    }
}

/*void QueryScore::setPrefilteringThresholds(){

    float s_per_pos = (float) ((float)scoresSum/(float)seqLenSum);
    float s_per_match = (float) ((float)scoresSum/(float)numMatches);
    
    for (int i = 0; i < dbSize; i++){
        float seqLen = (float) seqLens[i];
        float mean = s_per_pos * seqLen;
        float stddev = sqrt(seqLen * s_per_pos * s_per_match);
        int threshold = zscore_thr * stddev + mean;
        unsigned short ushort_threshold;
        if (threshold >= USHRT_MAX)
            ushort_threshold = USHRT_MAX;
        else
            ushort_threshold = (unsigned short) threshold;
        thresholds[i] = ushort_threshold;
    }
}*/

float QueryScore::getZscore(int seqId){
    return ( (float)scores[seqId] - s_per_pos * seqLens[seqId] ) / sqrt(s_per_pos * seqLens[seqId] * s_per_match);
}

float QueryScore::getZscoreRevSeq(int seqId){
    return  ( (float)scores[seqId] / (this->s_per_match * sqrt( 2.0 * this->matches_per_pos * seqLens[seqId] )) );
}

std::list<hit_t>* QueryScore::getResult (int querySeqLen,  float (QueryScore::*calcZscore)(int)){

//    const __m128i zero = _mm_setzero_si128(); 

    __m128i* p = scores_128;
    __m128i* thr = thresholds_128;

    __m128i cmp;
//    __m128i tmp;

    int set_cnt = 0;
    int set_bits_cnt = 0;
    // go through each vector
    for (int pos = 0; pos < scores_128_size/8; pos++ ){

        // look for entries above the threshold
/*        tmp = _mm_subs_epi16(*p, *thr);

        cmp = _mm_cmpeq_epi16(tmp, zero);*/
        cmp = _mm_cmpgt_epi16(*thr, *p);
        const unsigned int cmp_set_bits = _mm_movemask_epi8(cmp);
        
        // here are some sequences above the prefiltering threshold
        if (cmp_set_bits != 0xffff){
            set_cnt++;
            // and search for highest
            for(int i = 0; i < 8; i++){
                if(!CHECK_BIT(cmp_set_bits,i*2)){
                    set_bits_cnt++;
                    //hit_t hit = {pos * 8 + i, ((float)sse2_extract_epi16(*p,i))/(float)querySeqLen, 0.0};
                    float zscore = (this->*calcZscore)(pos*8+i); 
                    hit_t hit = {pos * 8 + i, zscore, 0.0};
                    resList->push_back(hit);
                }
            }
        }
        p++;
        thr++;
    }
    resList->sort(compareHits);
    return resList;
}

/*std::list<hit_t>* QueryScore::getResult (int querySeqLen, float (QueryScore::*calcZscore)(int)){

    setPrefilteringThresholds();
    float s_per_pos = (float)scoresSum/(float)seqLenSum;
    float s_per_match = (float)scoresSum/(float)numMatches;

    for (int i = 0; i < dbSize; i++){
        if (scores[i] > thresholds[i]){
            float zscore = (this->*calcZscore)(i);
            hit_t hit = {i, zscore, 0};
            resList->push_back(hit);
        }
    }
    resList->sort(compareHits);
    return resList;
}*/


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
    std::cerr << "Fatal error in QueryScore: position in the vector is not in the legal range (pos = " << pos << ")\n";
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
