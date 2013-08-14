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

QueryScore::QueryScore (int dbSize, unsigned short * seqLens, int k){

    this->dbSize = dbSize;

    // 8 DB short int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    scores = (unsigned short * ) scores_128;

    // set scores to zero
    memset (scores_128, 0, (dbSize/8 + 1) * 16);

    thresholds_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    thresholds = (unsigned short * ) thresholds_128;

    float_thresholds = new float[dbSize];

    memset (thresholds_128, 0, (dbSize/8 + 1) * 16);

    // initialize sequence lenghts with each seqLens[i] = L_i - k + 1
    this->seqLens_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    this->seqLens = (unsigned short *) this->seqLens_128;

    std::ofstream seq_lens_file;
    seq_lens_file.open("/net/cluster/user/maria/test/seq_lens.dat");
    for (int i = 0; i < dbSize; i++){
        if (seqLens[i] > (k - 1))
            this->seqLens[i] = seqLens[i] - k + 1;
        else
            this->seqLens[i] = 0;
        seq_lens_file << seqLens[i] << "\n";
    }
    seq_lens_file.close();

    this->seqLenSum = 0;
    for (int i = 0; i < dbSize; i++)
        this->seqLenSum += this->seqLens[i];

    this->resList = new std::list<hit_t>();

    scoresSum = 0;

    numMatches = 0;

    counter = 0;

    zscores = new float [dbSize];
    memset (zscores, 0.0, dbSize);
}

QueryScore::~QueryScore (){
    free(scores_128);
    delete resList;
}

    bool QueryScore::compareHitList(hit_t first, hit_t second){
        if (first.eval < second.eval)
            return true;
        return false;
    }

void QueryScore::setPrefilteringThresholds(){

/*    std::ofstream s_per_pos_file;
    s_per_pos_file.open("/net/cluster/user/maria/test/s_per_pos.dat", std::fstream::app);

    std::ofstream second_term_file;
    second_term_file.open("/net/cluster/user/maria/test/second_term.dat", std::fstream::app);

    std::ofstream norm_score1_file;
    norm_score1_file.open("/net/cluster/user/maria/test/norm_score1.dat", std::fstream::app);

    std::ofstream zscore_file;
    zscore_file.open("/net/cluster/user/maria/test/zscore.dat", std::fstream::app);

    std::ofstream std_dev_file;
    std_dev_file.open("/net/cluster/user/maria/test/std_dev.dat", std::fstream::app);
*/
    float s_per_pos = (float)scoresSum/(float)seqLenSum;
//    s_per_pos_file << s_per_pos << "\n";

//    s_per_pos_file.close();

    float s_per_match = (float)scoresSum/(float)numMatches;

    for (int i = 0; i < dbSize; i++){
        /* ######## ATTENTION pay attention to a possible unsigned short overflow! ######## */
        thresholds[i] = (unsigned short) (10.0 * sqrt(s_per_pos * (float)seqLens[i] * s_per_match)  + s_per_pos * (float)seqLens[i]);
//        float_thresholds[i] = 4.0 * sqrt(s_per_pos * (float)seqLens[i] * s_per_match)  + s_per_pos * (float)seqLens[i];
//        zscores[i] = ( (float)scores[i] - s_per_pos * (float)seqLens[i] ) / sqrt(s_per_pos * (float)seqLens[i] * s_per_match);
/*        if (scores[i] >= thresholds[i] && zscores[i] < 4.0){
            std::cout << (4.0 * sqrt(s_per_pos * (float)seqLens[i] * s_per_match)  + s_per_pos * (float)seqLens[i]) << " " << thresholds[i] << " " << scores[i] << "\n";
            std::cout << seqLens[i] << " " << s_per_pos << " " << s_per_match << " " << zscores[i] <<  "\n";
            std::cout << "\n";
        }*/
/*        float zscore = ( (float)scores[i] - s_per_pos * (float)seqLens[i] ) / sqrt(s_per_pos * (float)seqLens[i] * s_per_match);

        if (counter % 1000 == 0){
            second_term_file << ( s_per_pos * (float)seqLens[i] ) << "\n";
            norm_score1_file << ( (float)scores[i] - s_per_pos * (float)seqLens[i] ) << "\n";
            std_dev_file << (sqrt(s_per_pos * (float)seqLens[i] * s_per_match)) << "\n";
            zscore_file << zscore << "\n";
        }
        counter++;*/
    }
/*    second_term_file.close();
    norm_score1_file.close();
    std_dev_file.close();
    zscore_file.close();*/


}


std::list<hit_t>* QueryScore::getResult (int querySeqLen){
    setPrefilteringThresholds();

    float s_per_pos = (float)scoresSum/(float)seqLenSum;
    float s_per_match = (float)scoresSum/(float)numMatches;

    // minimum score for this sequence that satisfies the score per colum threshold
/*    int minScoreInt = (int) (prefThreshold * (float)querySeqLen);
    unsigned short minScore;
    if (minScoreInt < SHRT_MAX)
        minScore = (unsigned short) (prefThreshold * (float)querySeqLen);
    else
        minScore = SHRT_MAX;*/
    // set all elements of thr to the threshold score
//        const __m128i thr = _mm_set1_epi16(minScore);
    const __m128i zero = _mm_setzero_si128(); 

    __m128i* p = scores_128;
    __m128i* thr = thresholds_128;

    __m128i cmp;
    __m128i tmp;

    for (int pos = 0; pos < dbSize/8 + 1; pos++ ){

        // look for entries above the threshold
        tmp = _mm_subs_epu16(*p, *thr);

//                tmp = _mm_subs_epu16(*p, thr);
        cmp = _mm_cmpeq_epi16(tmp, zero);
        const int cmp_set_bits = _mm_movemask_epi8(cmp);
        // here are some sequences above the prefiltering threshold
        if (cmp_set_bits != 0){
            // and search for highest
            for(int i = 0; i < 8; i++){
                if(!CHECK_BIT(cmp_set_bits,i*2)){
                    //hit_t hit = {pos * 8 + i, ((float)sse2_extract_epi16(*p,i))/(float)querySeqLen, 0.0};
                    float zscore = ( (float)scores[pos * 8 + i] - s_per_pos * (float)seqLens[pos * 8 + i] ) / sqrt(s_per_pos * (float)seqLens[pos * 8 + i] * s_per_match);
                    hit_t hit = {pos * 8 + i, zscore, 0.0};
                    resList->push_back(hit);
                }
            }
        }
        p++;
    }

    resList->sort(compareHitList);
    return resList;
}

unsigned short QueryScore::sse2_extract_epi16(__m128i v, int pos) {
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


