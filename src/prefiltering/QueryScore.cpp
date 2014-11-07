#include "QueryScore.h"
#include "simd.h"


QueryScore::QueryScore (int dbSize, unsigned short * dbSeqLens, int k, short kmerThr,
                        float kmerMatchProb, float zscoreThr){
    
    this->dbSize = dbSize;
    this->kmerMatchProb = kmerMatchProb;
    this->kmerThr = kmerThr;
    this->zscore_thr = zscoreThr;
    this->scores_128_size = (dbSize + SIMD_SHORT_SIZE -1)/ SIMD_SHORT_SIZE * SIMD_SHORT_SIZE;
    // 8 DB short int entries are stored in one __m128i vector
    // one __m128i vector needs 16 byte
    scores_128 = (simd_int*) mem_align(ALIGN_INT, scores_128_size * 2);
    scores = (unsigned short * ) scores_128;
    // set scores to zero
    memset (scores_128, 0, scores_128_size * 2);
    
    this->resList = (hit_t *) mem_align(ALIGN_INT, MAX_RES_LIST_LEN * sizeof(hit_t) );
    
    scoresSum = 0;
    
    numMatches = 0;
    
}

QueryScore::~QueryScore (){
    free(scores_128);
    free(resList);
}

bool QueryScore::compareHits(hit_t first, hit_t second){
    return (first.zScore > second.zScore) ? true : false;
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
    EXIT(1);
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


