//
// Created by mad on 12/15/15.

#include "simd.h"
#include "UngappedAlignment.h"

// Bin size of the diagonal scoring. It follows the vector width, so it cannot live in the
// header any more (that would pull simd.h into every includer of UngappedAlignment.h).
#ifdef AVX512
static const unsigned int DIAGONALBINSIZE = 16;
#elif defined(AVX2)
static const unsigned int DIAGONALBINSIZE = 8;
#else
static const unsigned int DIAGONALBINSIZE = 4;
#endif


static void extractScores(unsigned int *score_arr, simd_int score) {
#ifdef AVX512
    // AVX512 has no 32-bit extract, so take the 256-bit half first
    #define EXTRACT_AVX(i) score_arr[i] = _mm256_extract_epi32(_mm512_extracti64x4_epi64(score, i <= 7 ? 0 : 1), i <= 7 ? i : i - 8)
    EXTRACT_AVX(0);  EXTRACT_AVX(1);  EXTRACT_AVX(2);  EXTRACT_AVX(3);
    EXTRACT_AVX(4);  EXTRACT_AVX(5);  EXTRACT_AVX(6);  EXTRACT_AVX(7);
    EXTRACT_AVX(8);  EXTRACT_AVX(9);  EXTRACT_AVX(10);  EXTRACT_AVX(11);
    EXTRACT_AVX(12);  EXTRACT_AVX(13);  EXTRACT_AVX(14);  EXTRACT_AVX(15);
#undef EXTRACT_AVX
#elif defined(AVX2)
    #define EXTRACT_AVX(i) score_arr[i] = _mm256_extract_epi32(score, i)
    EXTRACT_AVX(0);  EXTRACT_AVX(1);  EXTRACT_AVX(2);  EXTRACT_AVX(3);
    EXTRACT_AVX(4);  EXTRACT_AVX(5);  EXTRACT_AVX(6);  EXTRACT_AVX(7);
#undef EXTRACT_AVX
#else
#define EXTRACT_SSE(i) score_arr[i] = _mm_extract_epi32(score, i)
    EXTRACT_SSE(0);  EXTRACT_SSE(1);   EXTRACT_SSE(2);  EXTRACT_SSE(3);
#undef EXTRACT_SSE
#endif
}

UngappedAlignment::UngappedAlignment(const unsigned int maxSeqLen,
                                     BaseMatrix *substitutionMatrix, SequenceLookup *sequenceLookup,
                                     const unsigned char *dbRemap)
        : subMatrix(substitutionMatrix), sequenceLookup(sequenceLookup), dbRemap(dbRemap) {
    score_arr = new unsigned int[DIAGONALBINSIZE];
    diagonalCounter = new unsigned char[DIAGONALCOUNT];
    queryProfile   = (char *) malloc_simd_int((Sequence::PROFILE_AA_SIZE + 1) * maxSeqLen);
    memset(queryProfile, 0, (Sequence::PROFILE_AA_SIZE + 1) * maxSeqLen);
    aaCorrectionScore = (char *) malloc_simd_int(maxSeqLen);
    diagonalMatches = new CounterResult*[DIAGONALCOUNT * DIAGONALBINSIZE];
    if (dbRemap != NULL) {
        remapBufferSize = maxSeqLen * DIAGONALBINSIZE;
        remapBuffer = (unsigned char*)malloc(remapBufferSize);
    } else {
        remapBuffer = NULL;
        remapBufferSize = 0;
    }
}

UngappedAlignment::~UngappedAlignment() {
    delete [] diagonalMatches;
    free(aaCorrectionScore);
    free(queryProfile);
    delete [] diagonalCounter;
    delete [] score_arr;
    if (remapBuffer != NULL) {
        free(remapBuffer);
    }
}

void UngappedAlignment::align(CounterResult *results, size_t resultSize) {
    if (dbRemap != NULL) {
        computeScores<true>(queryProfile, queryLen, results, resultSize);
    } else {
        computeScores<false>(queryProfile, queryLen, results, resultSize);
    }
}


int UngappedAlignment::scalarDiagonalScoring(const char * profile,
                                             const unsigned int seqLen,
                                             const unsigned char * dbSeq) {
    int max = 0;
    int score = 0;
    for(unsigned int pos = 0; pos < seqLen; pos++){
        int curr = *((profile + pos * (Sequence::PROFILE_AA_SIZE + 1)) + dbSeq[pos]);
        score = curr + score;
        score = (score < 0) ? 0 : score;
        max = (score > max)? score : max;
    }
    return max;
}

template <unsigned int T>
void UngappedAlignment::unrolledDiagonalScoring(const char * profile,
                                                const unsigned int * seqLen,
                                                const unsigned char ** dbSeq,
                                                unsigned int * max) {
    unsigned int maxScores[DIAGONALBINSIZE];
    simd_int zero = simdi32_set(0);
    simd_int maxVec = simdi32_set(0);
    simd_int score = simdi32_set(0);

    for(unsigned int pos = 0; pos < seqLen[0]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore0 =  profileColumn[dbSeq[0][pos]];
        int subScore1 =  profileColumn[dbSeq[1][pos]];
        int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];
#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, subScore1, subScore0);
#elif defined(AVX512)
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, subScore1, subScore0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, subScore1, subScore0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[0]; pos < seqLen[1]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore1 =  profileColumn[dbSeq[1][pos]];
        int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];
#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, subScore1, 0);
#elif defined(AVX512)
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, subScore1, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, subScore1, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[1]; pos < seqLen[2]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];
#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, 0, 0);
#elif defined(AVX512)
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, 0, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, 0, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[2]; pos < seqLen[3]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore3 =  profileColumn[dbSeq[3][pos]];
#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, 0, 0, 0);
#elif defined(AVX512)
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, subScore4, subScore3, 0, 0, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, 0, 0, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
#if defined(AVX2)
    for(unsigned int pos = seqLen[3]; pos < seqLen[4]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[4]; pos < seqLen[5]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[5]; pos < seqLen[6]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[6]; pos < seqLen[7]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
#elif defined(AVX512)
    for(unsigned int pos = seqLen[3]; pos < seqLen[4]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, subScore4, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[4]; pos < seqLen[5]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, subScore5, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[5]; pos < seqLen[6]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, subScore6, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[6]; pos < seqLen[7]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            subScore7, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[7]; pos < seqLen[8]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore8 =  profileColumn[dbSeq[8][pos]];
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, subScore8,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[8]; pos < seqLen[9]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore9 =  profileColumn[dbSeq[9][pos]];
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, subScore9, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[9]; pos < seqLen[10]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore10 =  profileColumn[dbSeq[10][pos]];
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, subScore10, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[10]; pos < seqLen[11]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore11 =  profileColumn[dbSeq[11][pos]];
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, subScore11, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[11]; pos < seqLen[12]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore12 =  profileColumn[dbSeq[12][pos]];
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, subScore12, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[12]; pos < seqLen[13]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore13 =  profileColumn[dbSeq[13][pos]];
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, subScore13, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[13]; pos < seqLen[14]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore14 =  profileColumn[dbSeq[14][pos]];
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, subScore14, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[14]; pos < seqLen[15]; pos++){
        const char * profileColumn = (profile + pos * T);
        int subScore15 =  profileColumn[dbSeq[15][pos]];
        simd_int subScores = _mm512_set_epi32(
            subScore15, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        maxVec = simdui8_max(maxVec, score);
    }
#endif

    extractScores(maxScores, maxVec);

    for(size_t i = 0; i < DIAGONALBINSIZE; i++){
        max[i] = std::max(maxScores[i], max[i]);
    }
}

template <bool HasRemap>
void UngappedAlignment::scoreDiagonalAndUpdateHits(const char * queryProfile,
                                                   const unsigned int queryLen,
                                                   const short diagonal,
                                                   CounterResult ** hits,
                                                   const unsigned int hitSize) {
    //    unsigned char minDistToDiagonal = distanceFromDiagonal(diagonal);
    //    unsigned char maxDistToDiagonal = (minDistToDiagonal == 0) ? 0 : (DIAGONALCOUNT - minDistToDiagonal);
    //    unsigned int i_splits = computeSplit(queryLen, minDistToDiagonal);
    unsigned short minDistToDiagonal = distanceFromDiagonal(diagonal);

    if(queryLen >= 32768){
        for (size_t hitIdx = 0; hitIdx < hitSize; hitIdx++) {
            const DBLocalId seqId = hits[hitIdx]->id;
            unsigned int dbLen;
            const unsigned char *dbPtr = getDbSeq<HasRemap>(seqId, dbLen);
            std::pair<const unsigned char *, const unsigned int> dbSeq = std::make_pair(dbPtr, dbLen);
            int max = computeLongScore(queryProfile, queryLen, dbSeq, diagonal);
            hits[hitIdx]->count = static_cast<unsigned char>(std::min(255, max));
        }
        return;
    }
    memset(score_arr, 0, sizeof(unsigned int) * DIAGONALBINSIZE);
    if (hitSize == DIAGONALBINSIZE) {
        struct DiagonalSeq{
            unsigned char * seq;
            unsigned int seqLen;
            unsigned int id;
            static bool compareDiagonalSeqByLen(const DiagonalSeq &first, const DiagonalSeq &second) {
                return first.seqLen < second.seqLen;
            }
        };
        DiagonalSeq seqs[DIAGONALBINSIZE];
        unsigned int bufOffset = 0;
        for (unsigned int seqIdx = 0; seqIdx < hitSize; seqIdx++) {
            unsigned int tmpLen;
            const unsigned char *tmpPtr = getDbSeq<HasRemap>(hits[seqIdx]->id, tmpLen, bufOffset);
            if(tmpLen >= 32768){
                // hack to avoid too long sequences
                // this sequences will be processed by computeLongScore later
                seqs[seqIdx].seq = (unsigned char *) tmpPtr;
                seqs[seqIdx].seqLen = 0;
                seqs[seqIdx].id = seqIdx;
            }else{
                seqs[seqIdx].seq = (unsigned char *) tmpPtr;
                seqs[seqIdx].seqLen = tmpLen;
                seqs[seqIdx].id = seqIdx;
                if (HasRemap) {
                    bufOffset += tmpLen;
                }
            }
        }
        std::sort(seqs, seqs+DIAGONALBINSIZE, DiagonalSeq::compareDiagonalSeqByLen);
        unsigned int targetMaxLen = seqs[DIAGONALBINSIZE-1].seqLen;
        if (diagonal >= 0 && minDistToDiagonal < queryLen) {
            const unsigned char * tmpSeqs[DIAGONALBINSIZE];
            unsigned int seqLength[DIAGONALBINSIZE];
            unsigned int minSeqLen = std::min(targetMaxLen, queryLen - minDistToDiagonal);
            for(size_t i = 0; i < DIAGONALBINSIZE; i++) {
                tmpSeqs[i] = seqs[i].seq;
                seqLength[i] = std::min(seqs[i].seqLen, minSeqLen);
            }
            unrolledDiagonalScoring<Sequence::PROFILE_AA_SIZE + 1>(queryProfile + (minDistToDiagonal * (Sequence::PROFILE_AA_SIZE + 1)),
                                                                   seqLength, tmpSeqs, score_arr);

        } else if (diagonal < 0 && minDistToDiagonal < targetMaxLen) {
            const unsigned char * tmpSeqs[DIAGONALBINSIZE];
            unsigned int seqLength[DIAGONALBINSIZE];
            unsigned int minSeqLen = std::min(targetMaxLen - minDistToDiagonal, queryLen);
            for(size_t i = 0; i < DIAGONALBINSIZE; i++) {
                tmpSeqs[i] = seqs[i].seq + minDistToDiagonal;
                seqLength[i] = (seqs[i].seqLen > minDistToDiagonal) ? std::min(seqs[i].seqLen - minDistToDiagonal, minSeqLen) : 0;
            }
            unrolledDiagonalScoring<Sequence::PROFILE_AA_SIZE + 1>(queryProfile, seqLength,
                                                                   tmpSeqs, score_arr);
        }

        // update score
        for(size_t hitIdx = 0; hitIdx < hitSize; hitIdx++){
            CounterResult *const hit = hits[seqs[hitIdx].id];
            hit->count = static_cast<unsigned char>(std::min(static_cast<unsigned int>(255),
                                                           score_arr[hitIdx]));
            if(seqs[hitIdx].seqLen == 0){
                unsigned int dbLen2;
                const unsigned char *dbPtr2 = getDbSeq<HasRemap>(hit->id, dbLen2);
                if(dbLen2 >= 32768){
                    std::pair<const unsigned char *, const unsigned int> dbSeq2 = std::make_pair(dbPtr2, dbLen2);
                    int max = computeLongScore(queryProfile, queryLen, dbSeq2, diagonal);
                    hit->count = static_cast<unsigned char>(std::min(255, max));
                }
            }
        }
    }else {
        for (size_t hitIdx = 0; hitIdx < hitSize; hitIdx++) {
            const DBLocalId seqId = hits[hitIdx]->id;
            unsigned int dbLen;
            const unsigned char *dbPtr = getDbSeq<HasRemap>(seqId, dbLen);
            std::pair<const unsigned char *, const unsigned int> dbSeq = std::make_pair(dbPtr, dbLen);
            int max;
            if(dbSeq.second >= 32768){
                max = computeLongScore(queryProfile, queryLen, dbSeq, diagonal);
            }else{
                max = computeSingelSequenceScores(queryProfile, queryLen, dbSeq, diagonal, minDistToDiagonal);
            }
            hits[hitIdx]->count = static_cast<unsigned char>(std::min(255, max));
        }
    }
}

int UngappedAlignment::computeLongScore(const char * queryProfile, unsigned int queryLen,
                                        std::pair<const unsigned char *, const unsigned int> &dbSeq,
                                        unsigned short diagonal){
    int totalMax=0;
    for(unsigned int devisions = 1; devisions <= 1+ dbSeq.second /32768; devisions++ ){
        int realDiagonal = (-devisions * 65536  + diagonal);
        int minDistToDiagonal = abs(realDiagonal);
        int max = computeSingelSequenceScores(queryProfile, queryLen, dbSeq, realDiagonal, minDistToDiagonal);
        totalMax = std::max(totalMax, max);
    }
    for(unsigned int devisions = 0; devisions <= queryLen/65536; devisions++ ) {
        int realDiagonal = (devisions*65536+diagonal);
        int minDistToDiagonal = abs(realDiagonal);
        int max = computeSingelSequenceScores(queryProfile, queryLen, dbSeq, realDiagonal, minDistToDiagonal);
        totalMax = std::max(totalMax, max);
    }
    return totalMax;
}

template <bool HasRemap>
void UngappedAlignment::computeScores(const char *queryProfile,
                                      const unsigned int queryLen,
                                      CounterResult * results,
                                      const size_t resultSize) {
    memset(diagonalCounter, 0, DIAGONALCOUNT * sizeof(unsigned char));
    for(size_t i = 0; i < resultSize; i++){
//        // skip all that count not find enough diagonals
//        if(results[i].count < thr){
//            continue;
//        }
        const unsigned short currDiag = results[i].diagonal;
        // skip results that already have a diagonal score
        if(results[i].count != 0){
            continue;
        }
        diagonalMatches[currDiag * DIAGONALBINSIZE + diagonalCounter[currDiag]] = &results[i];
        diagonalCounter[currDiag]++;
        if(diagonalCounter[currDiag] == DIAGONALBINSIZE) {
            scoreDiagonalAndUpdateHits<HasRemap>(queryProfile, queryLen, static_cast<short>(currDiag),
                                       &diagonalMatches[currDiag * DIAGONALBINSIZE], diagonalCounter[currDiag]);
            diagonalCounter[currDiag] = 0;
        }
    }
    // process rest
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        if(diagonalCounter[i] > 0){
            scoreDiagonalAndUpdateHits<HasRemap>(queryProfile, queryLen, static_cast<short>(i),
                                       &diagonalMatches[i * DIAGONALBINSIZE], diagonalCounter[i]);
        }
        diagonalCounter[i] = 0;
    }
}

unsigned short UngappedAlignment::distanceFromDiagonal(const unsigned short diagonal) {
    const unsigned short zero = 0;
    const unsigned short dist1 =  zero - diagonal;
    const unsigned short dist2 =  diagonal - zero;
    return std::min(dist1 , dist2);
}



template <bool HasRemap>
inline const unsigned char* UngappedAlignment::getDbSeq(DBLocalId seqId, unsigned int &outLen, unsigned int bufferOffset) {
    std::pair<const unsigned char *, const unsigned int> raw = sequenceLookup->getSequence(seqId);
    outLen = raw.second;
    if (HasRemap) {
        unsigned int needed = bufferOffset + raw.second;
        if (needed > remapBufferSize) {
            remapBufferSize = needed * 2;
            remapBuffer = (unsigned char*)realloc(remapBuffer, remapBufferSize);
        }
        unsigned char *dst = remapBuffer + bufferOffset;
        for (unsigned int i = 0; i < raw.second; i++) {
            dst[i] = dbRemap[raw.first[i]];
        }
        return dst;
    }
    return raw.first;
}

void UngappedAlignment::createProfile(Sequence *seq,
                                      float * biasCorrection,
                                      const unsigned char *numSeqOverride) {
    queryLen = seq->L;
    memset(queryProfile, 0, (Sequence::PROFILE_AA_SIZE + 1) * seq->L);
    if(Parameters::isEqualDbtype(seq->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)) {
        // profile path: aaCorrectionScore not used
    } else if (biasCorrection != NULL) {
        for (int pos = 0; pos < seq->L; pos++) {
            float aaCorrBias = biasCorrection[pos];
            aaCorrBias = (aaCorrBias < 0.0) ? aaCorrBias/4 - 0.5 : aaCorrBias/4 + 0.5;
            aaCorrectionScore[pos] = static_cast<char>(aaCorrBias);
        }
    } else {
        memset(aaCorrectionScore, 0, seq->L);
    }
    // create profile
    if(Parameters::isEqualDbtype(seq->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)) {
        const int8_t * profile_aln = seq->getAlignmentProfile();
        for (int pos = 0; pos < seq->L; pos++) {
            for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                queryProfile[pos * (Sequence::PROFILE_AA_SIZE + 1) + aa_num] = (profile_aln[aa_num * seq->L + pos]);
            }
        }
    }else{
        const unsigned char *numSeq = (numSeqOverride != NULL) ? numSeqOverride : seq->numSequence;
        for (int pos = 0; pos < seq->L; pos++) {
            unsigned int aaIdx = numSeq[pos];
            for (int i = 0; i < subMatrix->alphabetSize; i++) {
                queryProfile[pos * (Sequence::PROFILE_AA_SIZE + 1) + i] = (subMatrix->subMatrix[aaIdx][i] + aaCorrectionScore[pos]);
            }
        }
    }
}

int UngappedAlignment::computeSingelSequenceScores(const char *queryProfile, const unsigned int queryLen,
                                                   std::pair<const unsigned char *, const unsigned int> &dbSeq,
                                                   int diagonal, unsigned int minDistToDiagonal) {
    int max = 0;
    if(diagonal >= 0 && minDistToDiagonal < queryLen){
        unsigned int minSeqLen = std::min(dbSeq.second, queryLen - minDistToDiagonal);
        int scores = scalarDiagonalScoring(queryProfile + (minDistToDiagonal * (Sequence::PROFILE_AA_SIZE+1)), minSeqLen, dbSeq.first);
        max = std::max(scores, max);
    }else if(diagonal < 0 && minDistToDiagonal < dbSeq.second){
        unsigned int minSeqLen = std::min(dbSeq.second - minDistToDiagonal, queryLen);
        int scores = scalarDiagonalScoring(queryProfile, minSeqLen, dbSeq.first + minDistToDiagonal);
        max = std::max(scores, max);
    }
    return max;
}


int UngappedAlignment::scoreSingelSequenceByCounterResult(CounterResult &result) {
    unsigned int dbLen;
    const unsigned char *dbPtr;
    if (dbRemap != NULL) {
        dbPtr = getDbSeq<true>(result.id, dbLen);
    } else {
        dbPtr = getDbSeq<false>(result.id, dbLen);
    }
    unsigned short minDistToDiagonal = distanceFromDiagonal(result.diagonal);
    return scoreSingleSequence(std::make_pair(dbPtr, dbLen), result.diagonal, minDistToDiagonal);
}

int UngappedAlignment::scoreSingleSequence(std::pair<const unsigned char *, const unsigned int> dbSeq,
                                           unsigned short diagonal,
                                           unsigned short minDistToDiagonal) {
    if(queryLen >= 32768 || dbSeq.second >= 32768) {
        return computeLongScore(queryProfile, queryLen, dbSeq, diagonal);
    } else {
        return computeSingelSequenceScores(queryProfile,queryLen ,dbSeq, static_cast<short>(diagonal), minDistToDiagonal);
    }
}
