#include "QueryScore.h"

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))


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

QueryScore::QueryScore (int dbSize, float prefThreshold){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;

    // 8 DB int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    scores = (unsigned short * ) scores_128;
    
    // set scores to zero
    memset (scores_128, 0, (dbSize/8 + 1) * 16);

    //this->lastMatchPos = new short[dbSize];
    //memset (lastMatchPos, 0, sizeof(short) * dbSize);
    this->lastScores = new LastScore[dbSize];
    memset (this->lastScores, 0, sizeof(LastScore) * dbSize);

    this->hitList = new DynamicArray();
    this->resList = new std::list<hit_t>();

    dbFractCnt = 0.0;
    qSeqCnt = 0;

}

QueryScore::~QueryScore (){
    delete lastScores;
    free(scores_128);
    delete hitList;
    delete resList;
}




    
std::list<hit_t>* QueryScore::getResult (int querySeqLen){
    // minimum score for this sequence that satisfies the score per colum threshold
    const unsigned short minScore = (unsigned short) (prefThreshold * (float)querySeqLen);
    // set all elements of thr to the threshold score
    const __m128i thr = _mm_set1_epi16(minScore);
    // temporary help vectors
    
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
                        hit_t hit = {pos * 8 + i, ((float)sse2_extract_epi16(*p,i))/(float)querySeqLen};
                        resList->push_back(hit);
                    }
                }
        }
        p++;
    }

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



void QueryScore::printStats(){
    std::cout << "Average occupancy of the DB scores array: " << dbFractCnt/(double)qSeqCnt << "\n";
}

void QueryScore::printVector(__m128i v){
    for (int i = 0; i < 8; i++)
        std::cout << sse2_extract_epi16(v, i) << " ";
    std::cout << "\n";
}



