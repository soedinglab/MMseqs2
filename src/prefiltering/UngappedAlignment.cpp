//
// Created by mad on 12/15/15.

#include "UngappedAlignment.h"

UngappedAlignment::UngappedAlignment(const unsigned int maxSeqLen,
                                     BaseMatrix *substitutionMatrix, SequenceLookup *sequenceLookup)
        : subMatrix(substitutionMatrix), sequenceLookup(sequenceLookup) {
    score_arr = new unsigned int[DIAGONALBINSIZE];
    diagonalCounter = new unsigned char[DIAGONALCOUNT];
    queryProfile   = (char *) malloc_simd_int((Sequence::PROFILE_AA_SIZE + 1) * maxSeqLen);
    memset(queryProfile, 0, (Sequence::PROFILE_AA_SIZE + 1) * maxSeqLen);
    aaCorrectionScore = (char *) malloc_simd_int(maxSeqLen);
    diagonalMatches = new CounterResult*[DIAGONALCOUNT * DIAGONALBINSIZE];
}

UngappedAlignment::~UngappedAlignment() {
    delete [] diagonalMatches;
    free(aaCorrectionScore);
    free(queryProfile);
    delete [] diagonalCounter;
    delete [] score_arr;
}

void UngappedAlignment::processQuery(Sequence *seq,
                                     float *biasCorrection,
                                     CounterResult *results,
                                     size_t resultSize) {
    createProfile(seq, biasCorrection, subMatrix->subMatrix);
    queryLen = seq->L;
    computeScores(queryProfile, seq->L, results, resultSize);
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
//        std::cout << (int) dbSeq[pos] << "\t" << curr << "\t" << max << "\t" << score <<  "\t" << (curr - bias) << std::endl;
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
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, subScore1, subScore0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }

    for(unsigned int pos = seqLen[0]; pos < seqLen[1]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        int subScore1 =  profileColumn[dbSeq[1][pos]];
        int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];

#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, subScore1, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, subScore1, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }

    for(unsigned int pos = seqLen[1]; pos < seqLen[2]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];

#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, subScore2, 0, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, subScore2, 0, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }

    for(unsigned int pos = seqLen[2]; pos < seqLen[3]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        //int subScore2 =  profileColumn[dbSeq[2][pos]];
        int subScore3 =  profileColumn[dbSeq[3][pos]];

#ifdef AVX2
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, subScore3, 0, 0, 0);
#else
        simd_int subScores = _mm_set_epi32(subScore3, 0, 0, 0);
#endif
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }
#ifdef AVX2
    for(unsigned int pos = seqLen[3]; pos < seqLen[4]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        //int subScore2 =  profileColumn[dbSeq[2][pos]];
        //int subScore3 =  profileColumn[dbSeq[3][pos]];
        int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, subScore4, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[4]; pos < seqLen[5]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        //int subScore2 =  profileColumn[dbSeq[2][pos]];
        //int subScore3 =  profileColumn[dbSeq[3][pos]];
        //int subScore4 =  profileColumn[dbSeq[4][pos]];
        int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, subScore5, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[5]; pos < seqLen[6]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        //int subScore2 =  profileColumn[dbSeq[2][pos]];
        //int subScore3 =  profileColumn[dbSeq[3][pos]];
        //int subScore4 =  profileColumn[dbSeq[4][pos]];
        //int subScore5 =  profileColumn[dbSeq[5][pos]];
        int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, subScore6, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }
    for(unsigned int pos = seqLen[6]; pos < seqLen[7]; pos++){
        const char * profileColumn = (profile + pos * T);
        //int subScore0 =  profileColumn[dbSeq[0][pos]];
        //int subScore1 =  profileColumn[dbSeq[1][pos]];
        //int subScore2 =  profileColumn[dbSeq[2][pos]];
        //int subScore3 =  profileColumn[dbSeq[3][pos]];
        //int subScore4 =  profileColumn[dbSeq[4][pos]];
        //int subScore5 =  profileColumn[dbSeq[5][pos]];
        //int subScore6 =  profileColumn[dbSeq[6][pos]];
        int subScore7 =  profileColumn[dbSeq[7][pos]];
        simd_int subScores = _mm256_set_epi32(subScore7, 0, 0, 0, 0, 0, 0, 0);
        score = simdi32_add(score, subScores);
        score = simdi32_max(score, zero);
        //        std::cout << (int)((char *)&template01)[0] << "\t" <<  SSTR(((char *)&score_vec_8bit)[0]) << "\t" << SSTR(((char *)&vMaxScore)[0]) << "\t" << SSTR(((char *)&vscore)[0]) << std::endl;
        maxVec = simdui8_max(maxVec, score);
    }
#endif

    extractScores(maxScores, maxVec);

    for(size_t i = 0; i < DIAGONALBINSIZE; i++){
        max[i] = std::max(maxScores[i], max[i]);
    }
}

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
            const unsigned int seqId = hits[hitIdx]->id;
            std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(seqId);
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
        for (unsigned int seqIdx = 0; seqIdx < hitSize; seqIdx++) {
            std::pair<const unsigned char *, const unsigned int> tmp = sequenceLookup->getSequence(
                    hits[seqIdx]->id);
            if(tmp.second >= 32768){
                // hack to avoid too long sequences
                // this sequences will be processed by computeLongScore later
                seqs[seqIdx].seq = (unsigned char *) tmp.first;
                seqs[seqIdx].seqLen = 1;
                seqs[seqIdx].id = seqIdx;
            }else{
                seqs[seqIdx].seq = (unsigned char *) tmp.first;
                seqs[seqIdx].seqLen = (unsigned int) tmp.second;
                seqs[seqIdx].id = seqIdx;
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
                seqLength[i] = std::min(seqs[i].seqLen - minDistToDiagonal, minSeqLen);
            }
            unrolledDiagonalScoring<Sequence::PROFILE_AA_SIZE + 1>(queryProfile, seqLength,
                                                                   tmpSeqs, score_arr);
        }

        // update score
        for(size_t hitIdx = 0; hitIdx < hitSize; hitIdx++){
            hits[seqs[hitIdx].id]->count = static_cast<unsigned char>(std::min(static_cast<unsigned int>(255),
                                                                      score_arr[hitIdx]));
            if(seqs[hitIdx].seqLen == 1){
                std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(hits[hitIdx]->id);
                if(dbSeq.second >= 32768){
                    int max = computeLongScore(queryProfile, queryLen, dbSeq, diagonal);
                    hits[seqs[hitIdx].id]->count = static_cast<unsigned char>(std::min(255, max));
                }
            }
        }
    }else {
        for (size_t hitIdx = 0; hitIdx < hitSize; hitIdx++) {
            const unsigned int seqId = hits[hitIdx]->id;
            std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(seqId);
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
        diagonalMatches[currDiag * DIAGONALBINSIZE + diagonalCounter[currDiag]] = &results[i];
        diagonalCounter[currDiag]++;
        if(diagonalCounter[currDiag] == DIAGONALBINSIZE) {
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, static_cast<short>(currDiag),
                                       &diagonalMatches[currDiag * DIAGONALBINSIZE], diagonalCounter[currDiag]);
            diagonalCounter[currDiag] = 0;
        }
    }
    // process rest
    for(size_t i = 0; i < DIAGONALCOUNT; i++){
        if(diagonalCounter[i] > 0){
            scoreDiagonalAndUpdateHits(queryProfile, queryLen, static_cast<short>(i),
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

void UngappedAlignment::extractScores(unsigned int *score_arr, simd_int score) {
#ifdef AVX2
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


void UngappedAlignment::createProfile(Sequence *seq,
                                      float * biasCorrection,
                                      short **subMat) {

    if(Parameters::isEqualDbtype(seq->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)) {
        memset(queryProfile, 0, (Sequence::PROFILE_AA_SIZE + 1) * seq->L);
    }else{
        memset(queryProfile, 0, (Sequence::PROFILE_AA_SIZE + 1) * seq->L);
        for (int pos = 0; pos < seq->L; pos++) {
            float aaCorrBias = biasCorrection[pos];
            aaCorrBias = (aaCorrBias < 0.0) ? aaCorrBias/4 - 0.5 : aaCorrBias/4 + 0.5;
            aaCorrectionScore[pos] = static_cast<char>(aaCorrBias);
        }
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
        for (int pos = 0; pos < seq->L; pos++) {
            unsigned int aaIdx = seq->numSequence[pos];
            for (int i = 0; i < subMatrix->alphabetSize; i++) {
                queryProfile[pos * (Sequence::PROFILE_AA_SIZE + 1) + i] = (subMat[aaIdx][i] + aaCorrectionScore[pos]);
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
    std::pair<const unsigned char *, const unsigned int> dbSeq =  sequenceLookup->getSequence(result.id);
    unsigned short minDistToDiagonal = distanceFromDiagonal(result.diagonal);
    return scoreSingleSequence(dbSeq, result.diagonal, minDistToDiagonal);
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




