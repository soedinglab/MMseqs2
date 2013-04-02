#include "QueryScoreSSE.h"

QueryScore::QueryScore (int dbSize, float prefThreshold){

    this->dbSize = dbSize;
    this->prefThreshold = prefThreshold;

    // 8 DB int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores = (__m128i*) memalign(16, (dbSize/8 + 1) * 16);
    // set scores to zero
    memset (scores, 0, (dbSize/8 + 1) * 16);

    //this->lastMatchPos = new short[dbSize];
    //memset (lastMatchPos, 0, sizeof(short) * dbSize);

    this->hitList = new DynamicArray();
    this->resList = new std::list<hit_t>();

    dbFractCnt = 0.0;
    qSeqCnt = 0;

}

QueryScore::~QueryScore (){
    free(scores);
    delete hitList;
    delete resList;
}


void QueryScore::addScores (int* seqList, int seqListSize, unsigned short score){
    __m128i tmp = _mm_setzero_si128();
    int seqId = 0;
    for (int i = 0; i < seqListSize; i++){
        seqId = seqList[i];
        // set the score to add at the right position
        tmp = sse2_insert_epi16(tmp, score, seqId%8);
        // saturated add of the score to the scores array
        *(scores + seqId/8) = _mm_adds_epu16(*(scores + seqId/8), tmp);
        tmp = _mm_setzero_si128();
    }
}

std::list<hit_t>* QueryScore::getResult (int querySeqLen){
    // minimum score for this sequence that satisfies the score per colum threshold
    unsigned short minScore = (unsigned short) (prefThreshold * (float)querySeqLen);
    // set all elements of thr to the threshold score
    __m128i thr = _mm_set1_epi16(minScore);
    // temporary help vectors
    __m128i zero = _mm_setzero_si128();
    __m128i tmp;

    unsigned short cmp;
    unsigned short score;
    
    __m128i* p = scores;
    // go through all vectors
    for (int pos = 0; pos < dbSize/8 + 1; pos++){
        
        // look for entries above the threshold
        tmp = _mm_subs_epu16(*p, thr);
        tmp = _mm_cmpeq_epi16(tmp, zero);
        cmp = _mm_movemask_epi8(tmp);

        // here are some sequences above the prefiltering threshold
        if (cmp != 65535){
            // go through the vector and search for 
            for (int i = 0; i < 8; i++){
                score = sse2_extract_epi16(*p, i);
                if (score >= minScore){
                    // record the sequence ID with the score above the threshold
                    hit_t hit = {pos * 8 + i, (float)score/(float)querySeqLen};
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

__m128i QueryScore::sse2_insert_epi16(__m128i v, unsigned short val, int pos) {
    switch(pos){
        case 0: return _mm_insert_epi16(v, val, 0);
        case 1: return _mm_insert_epi16(v, val, 1);
        case 2: return _mm_insert_epi16(v, val, 2);
        case 3: return _mm_insert_epi16(v, val, 3);
        case 4: return _mm_insert_epi16(v, val, 4);
        case 5: return _mm_insert_epi16(v, val, 5);
        case 6: return _mm_insert_epi16(v, val, 6);
        case 7: return _mm_insert_epi16(v, val, 7);
    }
    std::cerr << "Fatal error in QueryScore: position in the vector is not in the legal range (pos = " << pos << ")\n";
    exit(1);
    // never executed
    return _mm_insert_epi16(v, val, 0);
}

void QueryScore::reset(){
    memset (scores, 0, (dbSize/8 + 1) * 16);
    resList->clear();
    hitList->clear();
}

void QueryScore::printStats(){
    std::cout << "Average occupancy of the DB scores array: " << dbFractCnt/(double)qSeqCnt << "\n";
} 

void QueryScore::printVector(__m128i v){
    for (int i = 0; i < 8; i++)
        std::cout << sse2_extract_epi16(v, i) << " ";
    std::cout << "\n";
}
