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

QueryScore::QueryScore (int dbSize, unsigned short * seqLens, float prefThreshold, int k){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;

    // 8 DB short int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    scores = (unsigned short * ) scores_128;

    // set scores to zero
    memset (scores_128, 0, (dbSize/8 + 1) * 16);

    // initialize sequence lenghts with each seqLens[i] = L_i - k + 1
    this->seqLens_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    this->seqLens = (unsigned short *) this->seqLens_128;

    for (int i = 0; i < dbSize; i++)
        this->seqLens[i] = seqLens[i] - k + 1;
    
    this->seqLenSum = 0;
    for (int i = 0; i < dbSize; i++)
        this->seqLenSum += this->seqLens[i];

    this->resList = new std::list<hit_t>();

    dbFractCnt = 0.0;
    qSeqCnt = 0;
    counter = 0;

}

QueryScore::~QueryScore (){
    free(scores_128);
    delete[] seqLens;
    delete resList;
}

inline unsigned short sadd16(unsigned short a, unsigned short b)
{
	unsigned int s = (unsigned int)(a+b);
	return -(s>>16) | (unsigned short)s;
}

void QueryScore::addScores (int* seqList, int seqListSize, unsigned short score){
    for (int i = 0; i < seqListSize; i++){
        scores[seqList[i]] = sadd16(scores[seqList[i]], score);
    }
}

float QueryScore::getPrefilteringThreshold (){
/*
    __m128i sum_v = _mm_setzero_si128();
    __m128i sqSum_v = _mm_setzero_si128();
    __m128i uintScores = _mm_setzero_si128();
    __m128i sqScores = _mm_setzero_si128();
    const __m128i zero = _mm_setzero_si128();
    
    __m128i* p = scores_128;

    int maxSetSize = 100000/8;
    if (this->dbSize < 100000)
        maxSetSize = this->dbSize/8;

    for (int pos = 0; pos < maxSetSize + 1; pos++){
        // add 4 lower short scores converted to int to the sum_v
        uintScores = _mm_unpackhi_epi16(*p, zero);
        sum_v = _mm_add_epi32(sum_v, uintScores);
        // add the squared scores
        sqScores = _mm_mul_epu32(uintScores, uintScores);
        sqSum_v = _mm_add_epi64(sqSum_v, sqScores);

        uintScores = _mm_srli_si128(uintScores, 4);
        sqScores = _mm_mul_epu32(uintScores, uintScores);
        sqSum_v = _mm_add_epi64(sqSum_v, sqScores);
        
        // add 4 higher short scores converted to int to the sum_v
        uintScores = _mm_unpacklo_epi16(*p, zero);
        sum_v = _mm_add_epi32(sum_v, uintScores);
        // add the squared scores
        sqScores = _mm_mul_epu32(uintScores, uintScores);
        sqSum_v = _mm_add_epi64(sqSum_v, sqScores);

        uintScores = _mm_srli_si128(uintScores, 4);
        sqScores = _mm_mul_epu32(uintScores, uintScores);
        sqSum_v = _mm_add_epi64(sqSum_v, sqScores);

        p++;
    }

    unsigned int sum = 0;
    unsigned long int sqSum = 0;

    sum += _mm_extract_epi32(sum_v, 0); 
    sum += _mm_extract_epi32(sum_v, 1);
    sum += _mm_extract_epi32(sum_v, 2);
    sum += _mm_extract_epi32(sum_v, 3);

    sqSum += _mm_extract_epi32(sqSum_v, 0);
    sqSum += _mm_extract_epi32(sqSum_v, 1);
    sqSum += _mm_extract_epi32(sqSum_v, 2);
    sqSum += _mm_extract_epi32(sqSum_v, 3);

    double mu_S = (double)sum/(double)(maxSetSize*8);
    double sigma_S = sqrt((double)sqSum/(double)(maxSetSize*8) - (mu_S * mu_S));

//    std::cout << "\nmu_S: " << mu_S << "\n";
//    std::cout << "sigma_S: " << sigma_S << "\n";

    double sThr = mu_S + sigma_S*4.0;
//    std::cout << "sThr: " << sThr << "\n";

    */

    std::cout << "----- " << counter << " ------\n";
    double sum = 0.0;
    int cnt = 0;
    double s = 0.0;
    for (int i = 0; i < dbSize; i++){
        s = ((double)scores[i])/((double)seqLens[i]);
        if (s > 0.3){
            sum += s - 0.3;
            cnt++;
        }
    }

    double mu_S = sum/(double)cnt;

    std::cout << "mu_S: " << mu_S << "\n";
    std::cout << "cnt: " << cnt << "\n";

    std::cout << "f(x) = " << cnt << " * 0.01 * 1.0/" << mu_S << "*exp(-(x-0.3)/" << mu_S << ")\n";

    std::ofstream outfile;
    std::ostringstream fname;
    fname << "/cluster/user/maria/test/distr/" << counter << ".out";
    outfile.open (fname.str().c_str());
    for (int i = 0; i < dbSize; i++){
        s = ((double)scores[i])/((double)seqLens[i]);
        if (s > 0.3){
            outfile << i << "\t" << scores[i] << "\t" << s << "\n";
        }
    }
    outfile.close();

    double eval = 0.0;
    for (int i = 0; i < dbSize; i++){
        if (seqLens[i] == 0)
            continue;
        s = ((double)scores[i])/((double)seqLens[i]);
        eval = ((double)cnt) * exp ((0.3 - s)/mu_S);
        if (eval < 1.0){
            hit_t hit = {i, s, eval};
            resList->push_back(hit);
        }
    }

    resList->sort(compareHitList);

    counter++;
    return 0.0;
}

bool QueryScore::compareHitList(hit_t first, hit_t second){
    if (first.eval < second.eval)
        return true;
    else
        return false;
}

std::list<hit_t>* QueryScore::getResult (int querySeqLen){
//    getPrefilteringThreshold();
//    return resList;    
    // minimum score for this sequence that satisfies the score per colum threshold
    const unsigned short minScore = (unsigned short) (prefThreshold * (float)querySeqLen);
    // set all elements of thr to the threshold score
    const __m128i thr = _mm_set1_epi16(minScore);
    
    __m128i* p = scores_128;

    for (int pos = 0; pos < dbSize/8 + 1; pos++ ){
        // look for entries above the threshold
        const __m128i cmp = _mm_cmpgt_epi16(*p, thr);
        const int cmp_set_bits=_mm_movemask_epi8(cmp);
        // here are some sequences above the prefiltering threshold
        if (cmp_set_bits != 0){
            // and search for highest
            for(int i = 0; i < 8; i++){
                    if(CHECK_BIT(cmp_set_bits,i*2)){
                        hit_t hit = {pos * 8 + i, ((float)sse2_extract_epi16(*p,i))/(float)querySeqLen, 0.0};
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

void QueryScore::reset(){
    memset (scores_128, 0, (dbSize/8 + 1) * 16);
    resList->clear();
}

void QueryScore::printStats(){
    std::cout << "Average occupancy of the DB scores array: " << dbFractCnt/(double)qSeqCnt << "\n";
}

void QueryScore::printVector(__m128i v){
    for (int i = 0; i < 8; i++)
        std::cout << sse2_extract_epi16(v, i) << " ";
    std::cout << "\n";
}

