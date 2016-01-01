//
// Created by mad on 12/15/15.

#include "DiagonalMatcher.h"

DiagonalMatcher::DiagonalMatcher(const unsigned int maxSeqLen,
                                 BaseMatrix * substitutionMatrix,
                                 SequenceLookup * sequenceLookup) {
    score_arr = new unsigned int[VECSIZE_INT*4];
    vectorSequence = (unsigned char *) malloc_simd_int(VECSIZE_INT * 4 * maxSeqLen);
    queryProfile   = (char *) malloc_simd_int(PROFILESIZE * maxSeqLen);
    memset(queryProfile, 0, PROFILESIZE * maxSeqLen);

    diagonalMatches = new hit_t**[DIAGONALCOUNT];
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        diagonalMatches[i] = new hit_t *[VECSIZE_INT*4];
    }
    diagonalCounter = new unsigned char[DIAGONALCOUNT];
    this->subMatrix = substitutionMatrix;
    this->sequenceLookup = sequenceLookup;
}

DiagonalMatcher::~DiagonalMatcher() {
    delete [] score_arr;
    delete [] diagonalCounter;
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        delete [] diagonalMatches[i];
    }
    delete [] diagonalMatches;
    free(vectorSequence);
    free(queryProfile);
}

void DiagonalMatcher::processQuery(Sequence *query, float *biasCorrection,
                                   std::pair<hit_t *, size_t> results) {
    short **subMat = subMatrix->subMatrix2Bit;
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
        float aaBias = biasCorrection[pos];
        char aaBiasCorrection = (char) (aaBias < 0.0) ? aaBias - 0.5: aaBias + 0.5;

        for (size_t i = 0; i < subMatrix->alphabetSize; i++) {
            queryProfile[pos * PROFILESIZE + i] = (subMat[aaIdx][i] + aaBiasCorrection + bias);
//            std::cout << aaIdx << "\t" << (int) queryProfile[pos * PROFILESIZE + i] << "\t" << (int) subMat[aaIdx][i] << "\t" << (int) aaBiasCorrection << "\t" << (int) bias << std::endl;
        }
    }
    computeScores(queryProfile, query->L, results, bias);
}

const int DiagonalMatcher::scalarDiagonalScoring(const char * profile,
                                                 const int bias,
                                                 const unsigned int seqLen,
                                                 const unsigned char * dbSeq) {
    int max = 0;
    int score = 0;
    for(unsigned int pos = 0; pos < seqLen; pos++){
        int curr = *((profile + pos * PROFILESIZE) + dbSeq[pos]);
        score = (curr - bias) + score;
        score = (score < 0) ? 0 : score;
        std::cout << (int) dbSeq[pos] << "\t" << curr << "\t" << max << "\t" << score <<  "\t" << (curr - bias) << std::endl;
        max = (score > max)? score : max;
    }
    return std::min(max, 255);
}

#ifdef AVX2
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
#endif

const simd_int  DiagonalMatcher::vectorDiagonalScoring(const char *profile,
                                                       const char bias,
                                                       const unsigned int seqLen,
                                                       const unsigned char *dbSeq) {
    simd_int vscore        = simdi_setzero();
    simd_int vMaxScore     = simdi_setzero();
    const simd_int vBias   = simdi8_set(bias);
#ifdef SSE
    const simd_int sixten  = simdi8_set(16);
    const simd_int fiveten = simdi8_set(15);
#endif
    for(unsigned int pos = 0; pos < seqLen; pos++){
        simd_int template01 = simdi_load((simd_int *)&dbSeq[pos*VECSIZE_INT*4]);
#ifdef AVX2
        __m256i score_matrix_vec01 = _mm256_load_si256((simd_int *)&profile[pos * PROFILESIZE]);
        //TODO check

        __m256i score_vec_8bit = Shuffle(score_matrix_vec01, template01);
        //        __m256i score_vec_8bit = _mm256_shuffle_epi8(score_matrix_vec01, template01);
        //        __m256i lookup_mask01  = _mm256_cmpgt_epi8(sixten, template01); // 16 > t
        //        score_vec_8bit = _mm256_and_si256(score_vec_8bit, lookup_mask01);
#elif defined(SSE)
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

        vscore    = simdui8_adds(vscore, score_vec_8bit);
        vscore    = simdui8_subs(vscore, vBias);
        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        vMaxScore = simdui8_max(vMaxScore, vscore);

    }
    return vMaxScore;
}

std::pair<unsigned char *, unsigned int> DiagonalMatcher::mapSequences(std::pair<unsigned char *, unsigned int> * seqs,
                                                                       unsigned int seqCount) {
    unsigned int maxLen = 0;
    for(unsigned int seqIdx = 0; seqIdx < seqCount;  seqIdx++) {
        maxLen = std::max(seqs[seqIdx].second, maxLen);
    }
    memset(vectorSequence, 21, maxLen * VECSIZE_INT * 4 * sizeof(unsigned char));
    for(unsigned int seqIdx = 0; seqIdx < VECSIZE_INT*4;  seqIdx++){
        const unsigned char * seq  = seqs[seqIdx].first;
        const unsigned int seqSize = seqs[seqIdx].second;
        for(unsigned int pos = 0; pos < seqSize;  pos++){
            vectorSequence[pos * VECSIZE_INT * 4 + seqIdx] = seq[pos];
        }
    }
    return std::make_pair(vectorSequence, maxLen);
}

void DiagonalMatcher::scoreDiagonalAndUpdateHits(const char * queryProfile,
                                                 const unsigned int queryLen,
                                                 const unsigned short diagonal,
                                                 hit_t ** hits,
                                                 const unsigned int hitSize,
                                                 const short bias) {
    //    unsigned char minDistToDiagonal = distanceFromDiagonal(diagonal);
    //    unsigned char maxDistToDiagonal = (minDistToDiagonal == 0) ? 0 : (DIAGONALCOUNT - minDistToDiagonal);
    //    unsigned int i_splits = computeSplit(queryLen, minDistToDiagonal);
    unsigned short minDistToDiagonal = distanceFromDiagonal(diagonal);
    if(hitSize > (VECSIZE_INT * 4) / 16){
        std::pair<unsigned char *, unsigned int> seqs[VECSIZE_INT*4];
        for(unsigned int seqIdx = 0; seqIdx < hitSize;  seqIdx++) {
            std::pair<const unsigned char *, const unsigned int> tmp = sequenceLookup->getSequence(
                                                                                                   hits[seqIdx]->seqId);
            seqs[seqIdx] = std::make_pair((unsigned char *) tmp.first, (unsigned int) tmp.second);
        }
        std::pair<unsigned char *, unsigned int> seq = mapSequences(seqs, hitSize);

        simd_int vMaxScore = simdi_setzero();

        if(minDistToDiagonal < queryLen){
            unsigned int minSeqLen = std::min(seq.second, queryLen - minDistToDiagonal);
            simd_int ret = vectorDiagonalScoring(queryProfile + (minDistToDiagonal * PROFILESIZE ), bias, minSeqLen, seq.first);
            vMaxScore = simdui8_max(ret, vMaxScore);
        }
        if(minDistToDiagonal < seq.second){
            unsigned int minSeqLen = std::min(seq.second - diagonal, queryLen);
            simd_int ret = vectorDiagonalScoring(queryProfile, bias, minSeqLen, seq.first + minDistToDiagonal * VECSIZE_INT * 4 );
            vMaxScore = simdui8_max(ret, vMaxScore);
        }
        extractScores(score_arr, vMaxScore);
        // update score
        for(size_t hitIdx = 0; hitIdx < hitSize; hitIdx++){
            hits[hitIdx]->diagonalScore = normalizeScore(score_arr[hitIdx], seqs[hitIdx].second);
        }

    }else {
        for (size_t hitIdx = 0; hitIdx < hitSize; hitIdx++) {
            const unsigned int seqId = hits[hitIdx]->seqId;
            std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(seqId);
            int max = 0;

            if(minDistToDiagonal < queryLen){
                unsigned int minSeqLen = std::min(dbSeq.second, queryLen - minDistToDiagonal);
                int scores = scalarDiagonalScoring(queryProfile + (minDistToDiagonal * PROFILESIZE), bias, minSeqLen, dbSeq.first);
                max = std::max(scores, max);
            }
            if(minDistToDiagonal < dbSeq.second){
                unsigned int minSeqLen = std::min(dbSeq.second - minDistToDiagonal, queryLen);
                int scores = scalarDiagonalScoring(queryProfile, bias, minSeqLen, dbSeq.first + minDistToDiagonal);
                max = std::max(scores, max);
            }
            hits[hitIdx]->diagonalScore = normalizeScore(static_cast<unsigned char>(std::min(255, max)), dbSeq.second);
        }
    }
}
void DiagonalMatcher::computeScores(const char *queryProfile,
                                    const unsigned int queryLen,
                                    std::pair<hit_t *, size_t > results,
                                    const short bias) {
    memset(diagonalCounter, 0, DIAGONALCOUNT * sizeof(unsigned char));
    for(size_t i = 0; i < results.second; i++){
        const unsigned short currDiag = results.first[i].diagonal;
        diagonalMatches[currDiag][diagonalCounter[currDiag]] = &results.first[i];
        diagonalCounter[currDiag]++;
        if(diagonalCounter[currDiag] >= (VECSIZE_INT * 4) ) {
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, currDiag,
                                       diagonalMatches[currDiag], diagonalCounter[currDiag], bias);
            diagonalCounter[currDiag] = 0;
        }
    }
    // process rest
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        if(diagonalCounter[i] > 0){
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, i,
                                       diagonalMatches[i], diagonalCounter[i], bias);
        }
        diagonalCounter[i] = 0;
    }
}

unsigned short DiagonalMatcher::distanceFromDiagonal(const unsigned short diagonal) {
    const unsigned short zero = 0;
    const unsigned short dist1 =  zero - diagonal;
    const unsigned short dist2 =  diagonal - zero;
    return std::min(dist1 , dist2);
}

unsigned int DiagonalMatcher::computeSplit(unsigned int size, const unsigned int i) {
    unsigned int result = size / DIAGONALCOUNT + 1;
    // adjust for wrong memory access
    return ((result - 1)*DIAGONALCOUNT + i < size) ? result : result - 1;
}

void DiagonalMatcher::extractScores(unsigned int *score_arr, simd_int score) {
#ifdef AVX2
#define EXTRACT_AVX(i) score_arr[i] = _mm256_extract_epi8(score, i)
    EXTRACT_AVX(0);  EXTRACT_AVX(1);  EXTRACT_AVX(2);  EXTRACT_AVX(3);
    EXTRACT_AVX(4);  EXTRACT_AVX(5);  EXTRACT_AVX(6);  EXTRACT_AVX(7);
    EXTRACT_AVX(8);  EXTRACT_AVX(9);  EXTRACT_AVX(10);  EXTRACT_AVX(11);
    EXTRACT_AVX(12);  EXTRACT_AVX(13);  EXTRACT_AVX(14);  EXTRACT_AVX(15);
    EXTRACT_AVX(16);  EXTRACT_AVX(17);  EXTRACT_AVX(18);  EXTRACT_AVX(19);
    EXTRACT_AVX(20);  EXTRACT_AVX(21);  EXTRACT_AVX(22);  EXTRACT_AVX(23);
    EXTRACT_AVX(24);  EXTRACT_AVX(25);  EXTRACT_AVX(26);  EXTRACT_AVX(27);
    EXTRACT_AVX(28);  EXTRACT_AVX(29);  EXTRACT_AVX(30);  EXTRACT_AVX(31);
#undef EXTRACT_AVX
#elif defined(SSE)
#define EXTRACT_SSE(i) score_arr[i] = _mm_extract_epi8(score, i)
    EXTRACT_SSE(0);  EXTRACT_SSE(1);   EXTRACT_SSE(2);  EXTRACT_SSE(3);
    EXTRACT_SSE(4);  EXTRACT_SSE(5);   EXTRACT_SSE(6);  EXTRACT_SSE(7);
    EXTRACT_SSE(8);  EXTRACT_SSE(9);   EXTRACT_SSE(10); EXTRACT_SSE(11);
    EXTRACT_SSE(12); EXTRACT_SSE(13);  EXTRACT_SSE(14); EXTRACT_SSE(15);
#undef EXTRACT_SSE
#endif
}

unsigned char DiagonalMatcher::normalizeScore(const unsigned char score, const unsigned int len) {
    float log2Len = Util::flog2(static_cast<float>(len));
    float floatScore = static_cast<float>(score);
    return static_cast<unsigned char>((log2Len > floatScore) ? 0.0 : floatScore - log2Len );
}