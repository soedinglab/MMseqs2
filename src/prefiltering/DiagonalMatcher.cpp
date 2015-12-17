//
// Created by mad on 12/15/15.
//

#include <SubstitutionMatrix.h>
#include "DiagonalMatcher.h"
#include "QueryScore.h"


DiagonalMatcher::DiagonalMatcher(const unsigned int maxSeqLen,
                                 SubstitutionMatrix * substitutionMatrix,
                                 SequenceLookup * sequenceLookup) {
    score_arr = new unsigned int[VECSIZE_INT*4];
    vectorSequence = (unsigned char *) malloc_simd_int(VECSIZE_INT * 4 * maxSeqLen);
    queryProfile   = (char *) malloc_simd_int(PROFILESIZE * maxSeqLen);
    memset(queryProfile, 0, PROFILESIZE * maxSeqLen);
    this->subMatrix = substitutionMatrix;
    this->sequenceLookup = sequenceLookup;
}

DiagonalMatcher::~DiagonalMatcher() {
    delete [] score_arr;
    free(vectorSequence);
}

void DiagonalMatcher::processQuery(Sequence * query, hit_t * results, size_t resultSize) {
    short **subMat = subMatrix->subMatrix;

    short bias = 0;
    for (size_t i = 0; i < subMatrix->alphabetSize; i++) {
        for (size_t j = 0; j < subMatrix->alphabetSize; j++) {
            if (subMat[i][j] < bias) {
                bias = subMat[i][j];
            }
        }
    }
    bias = abs(bias);
    memset(queryProfile, bias, PROFILESIZE * query->L);
    // create profile
    for (size_t pos = 0; pos < query->L; pos++) {
        unsigned int aaIdx = query->int_sequence[pos];
        for (size_t i = 0; i < subMatrix->alphabetSize; i++) {
            queryProfile[pos * PROFILESIZE + i] = subMat[aaIdx][i]+bias;
        }
    }
    computeScores(queryProfile, query->L, results, resultSize, bias);

}

const int DiagonalMatcher::scalarDiagonalScoring(const char * profile,
                                                          const unsigned int profileSize,
                                                          const int bias,
                                                          const unsigned int seqLen,
                                                          const unsigned char * dbSeq) {
    int max = 0;
    int sum = 0;
    for(unsigned int pos = 0; pos < seqLen; pos++){
//            int curr = *((profile + (pos * 20)) + number[i][pos]);
        int curr = *((profile + pos * profileSize) +  dbSeq[pos]);
        sum = (curr - bias) + max;
        std::cout << (int) dbSeq[pos] << "\t" << curr << "\t" << bias << "\t" << (curr - bias) << std::endl;
        sum = (sum < 0) ? 0 : sum;
        max = (sum > max)? sum : max;
    }
    return std::min(max, 255);
}



inline const __m256i DiagonalMatcher::Shuffle(const __m256i & value, const __m256i & shuffle)
{
    const __m256i K0 = _mm256_setr_epi8(
            0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70,
            0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0);

    const __m256i K1 = _mm256_setr_epi8(
            0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0, 0xF0,
            0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70, 0x70);

    return _mm256_or_si256(_mm256_shuffle_epi8(value, _mm256_add_epi8(shuffle, K0)),
                           _mm256_shuffle_epi8(_mm256_permute4x64_epi64(value, 0x4E), _mm256_add_epi8(shuffle, K1)));
}

const simd_int  DiagonalMatcher::vectorDiagonalScoring(const char *profile,
                                                           const unsigned int profileSize,
                                                           const char bias,
                                                           const unsigned int seqLen,
                                                           const unsigned char *dbSeq) {
    simd_int vscore        = simdi_setzero();
    simd_int vMaxScore     = simdi_setzero();
    const simd_int vBias   = simdi8_set(bias);
    const simd_int sixten  = simdi8_set(16);
    const simd_int fiveten = simdi8_set(15);

    std::cout << std::endl;
    for(unsigned int pos = 0; pos < seqLen; pos++){
        simd_int template01 = simdi_load((simd_int *)&dbSeq[pos*VECSIZE_INT*4]);
#ifdef AVX2
        __m256i score_matrix_vec01 = _mm256_load_si256((simd_int *)&profile[pos * PROFILESIZE]);
        //TODO check

        __m256i score_vec_8bit = Shuffle(score_matrix_vec01, template01);
//        __m256i score_vec_8bit = _mm256_shuffle_epi8(score_matrix_vec01, template01);
//        __m256i lookup_mask01  = _mm256_cmpgt_epi8(sixten, template01); // 16 > t
//        score_vec_8bit = _mm256_and_si256(score_vec_8bit, lookup_mask01);
#elif SSE
        // each position has 32 byte
        // 20 scores and 12 zeros
        // load score 0 - 15
        __m128i score_matrix_vec01 = _mm_load_si128((__m128i *)&profile[pos * 32]);
        // load score 16 - 32
        __m128i score_matrix_vec16 = _mm_load_si128((__m128i *)&profile[pos * 32 + 16]);
        // parallel score lookup
        // _mm_shuffle_epi8
        // for i ... 16
        //   score01[i] = score_matrix_vec01[template01[i]%16]
        __m128i score01 =_mm_shuffle_epi8(score_matrix_vec01,template01);
        __m128i score16 =_mm_shuffle_epi8(score_matrix_vec16,template01);
        // t[i] < 16 => 0 - 15
        // example: template01: 02 15 12 18 < 16 16 16 16 => FF FF FF 00
        __m128i lookup_mask01 = _mm_cmplt_epi8(template01, sixten);
        // 15 < t[i] => 16 - xx
        // example: template01: 16 16 16 16 < 02 15 12 18 => 00 00 00 FF
        __m128i lookup_mask16 = _mm_cmplt_epi8(fiveten, template01);
        // score01 & lookup_mask01 => Score   Score   Score   NoScore
        score01 = _mm_and_si128(lookup_mask01,score01);
        // score16 & lookup_mask16 => NoScore NoScore NoScore Score
        score16 = _mm_and_si128(lookup_mask16,score16);
        //     Score   Score   Score NoScore
        // + NoScore NoScore NoScore   Score
        // =   Score   Score   Score   Score
        __m128i score_vec_8bit = _mm_add_epi8(score01,score16);
#endif

        std::cout << (int)((char *)&template01)[0] << "\t" <<  (int)((char *)&score_vec_8bit)[0] << "\t" << (int)((char *)&vBias)[0] << std::endl;
        vscore    = simdui8_adds(vscore, score_vec_8bit);
        vscore    = simdui8_subs(vscore, vBias);
        vMaxScore = simdui8_max(vMaxScore, vscore);
    }

    return vMaxScore;
}

std::pair<unsigned char *, unsigned int> DiagonalMatcher::mapSequences(hit_t * seqIds, unsigned int seqCount) {
    std::pair<unsigned char *, unsigned int> seqs[VECSIZE_INT*4];
    unsigned int maxLen = 0;
    for(unsigned int seqIdx = 0; seqIdx < seqCount;  seqIdx++) {
        std::pair<const unsigned char *, const unsigned int> tmp = sequenceLookup->getSequence(seqIds[seqIdx].seqId);
        seqs[seqIdx] = std::make_pair((unsigned  char *)tmp.first, (unsigned int) tmp.second);
        maxLen = (maxLen < seqs[seqIdx].second) ? seqs[seqIdx].second  : maxLen;
    }
    for(unsigned int pos = 0; pos < maxLen;  pos++){
        for(unsigned int seqIdx = 0; seqIdx < VECSIZE_INT*4;  seqIdx++){
            vectorSequence[pos * VECSIZE_INT * 4 + seqIdx] = (seqIdx < seqCount) ? seqs[seqIdx].first[pos] : 0;
        }
    }
    return std::make_pair(vectorSequence, maxLen);
}

void DiagonalMatcher::scoreDiagonalAndUpdateHits(const char * queryProfile,
                                                 const unsigned int queryLen,
                                                 const unsigned char diagonal,
                                                 hit_t * hits,
                                                 const unsigned int hitSize,
                                                 const short bias) {



    unsigned int i_splits = std::max((unsigned int)1, queryLen / DIAGONALCOUNT);
    if(hitSize > (VECSIZE_INT * 4) / 4){
        std::pair<unsigned char *, unsigned int> seq = mapSequences(hits, hitSize);
        simd_int vMaxScore = simdi_setzero();
        unsigned int j_splits = std::max((unsigned int)1, seq.second / DIAGONALCOUNT);
        for(unsigned int i = 0; i < i_splits; i++) {
            for(unsigned int j = 0; j < j_splits; j++) {
                simd_int ret = vectorDiagonalScoring(queryProfile + (i * DIAGONALCOUNT * PROFILESIZE + diagonal * PROFILESIZE ), PROFILESIZE, bias, queryLen,
                                                                   seq.first + (j * DIAGONALCOUNT * VECSIZE_INT*4));
                vMaxScore = simdui8_max(ret, vMaxScore);
            }
        }
#ifdef AVX2
#define EXTRACT_AVX(i) score_arr[i] = _mm256_extract_epi8(vMaxScore, i)
        EXTRACT_AVX(0);  EXTRACT_AVX(1);  EXTRACT_AVX(2);  EXTRACT_AVX(3);
        EXTRACT_AVX(4);  EXTRACT_AVX(5);  EXTRACT_AVX(6);  EXTRACT_AVX(7);
        EXTRACT_AVX(8);  EXTRACT_AVX(9);  EXTRACT_AVX(10);  EXTRACT_AVX(11);
        EXTRACT_AVX(12);  EXTRACT_AVX(13);  EXTRACT_AVX(14);  EXTRACT_AVX(15);
        EXTRACT_AVX(16);  EXTRACT_AVX(17);  EXTRACT_AVX(18);  EXTRACT_AVX(19);
        EXTRACT_AVX(20);  EXTRACT_AVX(21);  EXTRACT_AVX(22);  EXTRACT_AVX(23);
        EXTRACT_AVX(24);  EXTRACT_AVX(25);  EXTRACT_AVX(26);  EXTRACT_AVX(27);
        EXTRACT_AVX(28);  EXTRACT_AVX(29);  EXTRACT_AVX(30);  EXTRACT_AVX(31);
#undef EXTRACT_AVX
#elif SSE
        #define EXTRACT_SSE(i) score_arr[i] = _mm_extract_epi8(vMaxScore, i)
    EXTRACT_SSE(0);  EXTRACT_SSE(1);  EXTRACT_SSE(2);  EXTRACT_SSE(3);
    EXTRACT_SSE(4);  EXTRACT_SSE(5);  EXTRACT_SSE(6);  EXTRACT_SSE(7);
    EXTRACT_SSE(8);  EXTRACT_AVX(9);  EXTRACT_SSE(10);  EXTRACT_SSE(11);
    EXTRACT_SSE(12);  EXTRACT_SSE(13);  EXTRACT_SSE(14);  EXTRACT_SSE(15);
#undef EXTRACT_SSE
#endif

        // update score
        for(size_t hitIdx = 0; hitIdx < hitSize; hitIdx++){
            hits[hitIdx].diagonalScore = score_arr[hitIdx];
        }
    }else {
        for (size_t hitIdx = 0; hitIdx < hitSize; hitIdx++) {
            const unsigned int seqId = hits[hitIdx].seqId;
            std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(seqId);
            unsigned int j_splits = std::max((unsigned int)1, dbSeq.second / DIAGONALCOUNT);
            const unsigned int seqLen = std::min(queryLen, dbSeq.second);
            int max = 0;
            for(unsigned int i = 0; i < i_splits; i++) {
                for (unsigned int j = 0; j < j_splits; j++) {
                    int scores = scalarDiagonalScoring(queryProfile +  (i * 256 * PROFILESIZE + diagonal * PROFILESIZE), PROFILESIZE,
                                                                bias, seqLen, dbSeq.first );
                    max = std::max(scores, max);
                }
            }
            hits[hitIdx].diagonalScore = max;
        }
    }
}
void DiagonalMatcher::computeScores(const char *queryProfile, const unsigned int queryLen, hit_t *results,
                                    const size_t resultSize, const short bias) {
    memset(diagonalCounter, 0, DIAGONALCOUNT * sizeof(unsigned char));
    for(size_t i = 0; i < resultSize; i++){
        const unsigned char currDiag = results[i].diagonal;
        const unsigned char writePos = diagonalCounter[currDiag];
        diagonalMatches[currDiag][writePos] = &results[i];
        diagonalCounter[currDiag]++;
        if(diagonalCounter[currDiag] >= (VECSIZE_INT * 4) ) {
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, currDiag,
                                       *diagonalMatches[currDiag], diagonalCounter[currDiag], bias);
            diagonalCounter[currDiag] = 0;
        }
    }
    // process rest
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        if(diagonalCounter[i] > 0){
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, i,
                                       *diagonalMatches[i], diagonalCounter[i], bias);
        }
        diagonalCounter[i] = 0;
    }
}