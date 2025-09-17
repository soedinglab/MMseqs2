#include "Fwbw.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "Alignment.h"
#include "Matcher.h"
#include "Util.h"
#include "Parameters.h"
#include "simd.h"
#include "Sequence.h"
#include "Timer.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <numeric>
#include <cmath>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

struct States {
    const static uint8_t STOP=0;
    const static uint8_t M=1;
    const static uint8_t I=2;
    const static uint8_t D=3;
};

struct FWBWState {
    const static bool FORWARD = true;
    const static bool BACKWARD = false;
};

inline void calculate_max4(float& max, float& term1, float& term2, float& term3, float& term4, uint8_t& state) {
    if (term1 > term2) { max = term1; state = States::STOP; }
    else { max = term2; state = States::M; }
    if (term3 > max) { max = term3; state = States::I; }
    if (term4 > max) { max = term4; state = States::D; }
}

inline simd_float simdf32_prefixsum(simd_float a) {
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 4)));
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 8)));
#ifdef AVX2
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 16)));
#endif
    return a;
// Fallback scalar implementation
//     float buf[8];
//     simdf32_storeu(buf, a);

//     buf[1] += buf[0];
//     buf[2] += buf[1];
//     buf[3] += buf[2];
// #ifdef AVX2
//     buf[4] += buf[3];
//     buf[5] += buf[4];
//     buf[6] += buf[5];
//     buf[7] += buf[6];
// #endif

//     return simdf32_loadu(buf);
}

// FwBwAligner Constructor for general case: use profile scoring matrix
FwBwAligner::FwBwAligner(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
                : temperature(temperature), length(length), gapOpen(gapOpen), gapExtend(gapExtend), mact(mact), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {    
    blockCapacity = colsCapacity / length;
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    // Score Matrix (targetProfile 21xcolsCapacity)
    scoreForwardProfile = malloc_matrix<float>(21, colsCapacity);
    scoreForwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    scoreBackwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));
    // Blosum matrix
    blosum = malloc_matrix<float>(21, 21);
    for (int i = 0; i < subMat.alphabetSize; ++i) {
        for (int j = 0; j < subMat.alphabetSize; ++j) {
            blosum[i][j] = static_cast<float>(subMat.subMatrix[i][j]);
        }
    }
    if (backtrace) {
        btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
        S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
        S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    }
}

// FwBwAligner Constructor for user-defined scoring matrix
FwBwAligner::FwBwAligner(float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
                    : temperature(temperature), length(length), gapOpen(gapOpen), gapExtend(gapExtend), mact(mact), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {
    
    // scoreForward
    scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));

    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));

    if (backtrace != 0) {
        blosum = nullptr;
        btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
        S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
        S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    }
}

FwBwAligner::~FwBwAligner(){
    // matrices used both in lolalign and general cases
    free(zm);
    free(zmBlockPrev);
    free(zmBlockCurr);
    free(zeBlock);
    free(zfBlock);
    free(zmFirst);
    free(zeFirst);
    free(zfFirst);
    free(vj);
    free(wj);
    free(zInit);
    free(exp_ge_arr);

    // matrices used in only one case
    if (scoreForward != nullptr) {
       free(scoreForward);
    }
    if (scoreForwardProfile != nullptr) {
        free(scoreForwardProfile);
    }
    if (scoreForwardProfile_exp != nullptr) {
        free(scoreForwardProfile_exp);
    }
    if (scoreBackwardProfile_exp != nullptr) {
        free(scoreBackwardProfile_exp);
    }
    if (btMatrix != nullptr) {
        free(btMatrix);
    }
    if (blosum != nullptr) {
        free(blosum);
    }
    if (S_prev != nullptr) {
        free(S_prev);
    }
    if (S_curr != nullptr) {
        free(S_curr);
    }
    
}

//Reallocatation or Resizing
void FwBwAligner::reallocateProfile(size_t newColsCapacity) { // reallocate profile when colSeqLen(queryLen in general) exceeds colsCapacity
    free(scoreForwardProfile); scoreForwardProfile = malloc_matrix<float>(21, newColsCapacity);
    free(scoreForwardProfile_exp); scoreForwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
    free(scoreBackwardProfile_exp); scoreBackwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
}

template<>
void FwBwAligner::resizeMatrix<true,true>(size_t newRowLen, size_t newColLen) {
    //profile = true(scoreForwardProfile, scoreForwardProfile_exp, scoreBackwardProfile_exp)
    //backtrace = true(btMatrix, S_prev, S_curr)

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //backtrace
    free(S_prev); S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(S_curr); S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
}

template<>
void FwBwAligner::resizeMatrix<true,false>(size_t newRowLen, size_t newColLen) {
    //profile = true(scoreForwardProfile, scoreForwardProfile_exp, scoreBackwardProfile_exp)
    
    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
}

template<>
void FwBwAligner::resizeMatrix<false,true>(size_t newRowLen, size_t newColLen) {
    //profile = false(scoreForward)
    //backtrace = true(btMatrix, S_prev, S_curr)

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //profile
    free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
    //backtrace
    free(S_prev); S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(S_curr); S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
}

template<>
void FwBwAligner::resizeMatrix<false,false>(size_t newRowLen, size_t newColLen) {
    //profile = false(scoreForward)
    //backtrace = false

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //profile
    free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
}

void FwBwAligner::resetParams(float newGapOpen, float newGapExtend, float newTemperature) {
    gapOpen = newGapOpen;
    gapExtend = newGapExtend;
    temperature = newTemperature;
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
}


void FwBwAligner::initAlignment(unsigned char* targetAANum, size_t targetLen, size_t queryLen) {
    rowSeqAANum = targetAANum; rowSeqLen = targetLen;
    resizeMatrix<true, true>(targetLen, queryLen);
}

void FwBwAligner::initScoreMatrix(float** inputScoreMatrix, int* gaps) {
    //gaps: [rowStart, rowEnd, colStart, colEnd]
    rowSeqLen = gaps[1]-gaps[0]; colSeqLen = gaps[3]-gaps[2]; //row=tlen, col=qlen
    colSeqLen_padding = ((colSeqLen + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //colpadding
    blocks = (colSeqLen / length) + (colSeqLen % length != 0);
    simd_float vTemp = simdf32_set(temperature);
    if (colSeqLen > colsCapacity) {
        size_t newColsCapacity = ((colSeqLen + length-1)/length)* length;
        free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, newColsCapacity);
    }
    size_t colLoopCount = colSeqLen/VECSIZE_FLOAT;
    size_t colEndPos = colLoopCount*VECSIZE_FLOAT;
    for (size_t i = 0; i < rowSeqLen; ++i){
        for (size_t j = 0; j < colEndPos; j+=VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_loadu(&inputScoreMatrix[i+gaps[0]][j+gaps[2]]);
            vScoreForward = simdf32_div(vScoreForward, vTemp);
            simdf32_store(&scoreForward[i][j], vScoreForward);
        }
        for(size_t j = colEndPos; j < colSeqLen; ++j){
            scoreForward[i][j] = inputScoreMatrix[i+gaps[0]][j+gaps[2]] / temperature;

        }
        for (size_t j = colSeqLen; j < colSeqLen_padding; ++j){
            scoreForward[i][j] = FLT_MIN_EXP;
        } 
    }
}

void FwBwAligner::initProfile(unsigned char* colAANum, size_t colAALen) {
    colSeqAANum = colAANum; colSeqLen = colAALen;
    colSeqLen_padding = ((colSeqLen + VECSIZE_FLOAT -1) / VECSIZE_FLOAT) * VECSIZE_FLOAT;
    blocks = (colSeqLen / length) + (colSeqLen % length != 0);
    if (colSeqLen > colsCapacity) {
        size_t newColsCapacity = ((colSeqLen + length-1)/length)* length;
        reallocateProfile(newColsCapacity);
    }
    
    //scoreForward : 21 * qlen
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < colAALen; ++j) {
            float score = blosum[i][colSeqAANum[j]]/temperature;
            scoreForwardProfile[i][j] = score;
        }
        std::fill(&scoreForwardProfile[i][colSeqLen], &scoreForwardProfile[i][colSeqLen_padding], FLT_MIN_EXP);
    }
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < colSeqLen_padding; j += VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_load(&scoreForwardProfile[i][j]);
            vScoreForward = simdf32_exp(vScoreForward);
            simdf32_store(&scoreForwardProfile_exp[i][j], vScoreForward);
        }
        for (size_t j=0; j < colSeqLen; ++j){
            size_t reverse_j = colSeqLen - 1 - j;
            scoreBackwardProfile_exp[i][reverse_j] = scoreForwardProfile_exp[i][j];
        }
        //remainder 
        for (size_t j=colSeqLen; j < colSeqLen_padding; ++j){
            scoreBackwardProfile_exp[i][j] = 0;
        }
    }
}

template <bool profile>
void FwBwAligner::forward() {
    //Init zInit
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    max_zm = -std::numeric_limits<float>::max(); 
    simd_float vMax_zm = simdf32_set(max_zm);
    P = nullptr; // reset p. do we need this?
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, colSeqLen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));
            
        memcpy(zmFirst + 1, zInit[0], rowSeqLen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], rowSeqLen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], rowSeqLen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= rowSeqLen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (profile) {
                    vScoreMatrix = simdf32_exp(simdf32_load(&scoreForwardProfile[rowSeqAANum[i-1]][start + j - 1]));
                } else {
                    vScoreMatrix = simdf32_exp(simdf32_load(&scoreForward[i-1][start + j - 1]));
                } 
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            zmMaxRowBlock = simdf32_hmax(vZmax_tmp);
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);

            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, simdf32_set(exp_go)),
                                        simdf32_mul(vZf, simdf32_set(exp_ge))
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                // simd_float vZeUpdate = simdf32_fmadd(vZeI0, vExp_ge_arr, simdf32_div(vCumsumZm, vWj));
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            simd_float vCurrMax = simdf32_set(current_max);
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT){
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), vCurrMax);
                vMax_zm = simdf32_max(vMax_zm, vZmCurr);
                simdf32_store(&zm[i - 1][start + j-1], vZmCurr);
            }     

            #if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #endif


            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, vCurrMax);

            zInit[0][i-1] = vNextZinit[0]; zInit[1][i-1] = vNextZinit[1]; zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            
            if (i < rowSeqLen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max    
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
    }

    sum_exp= 0.0; 
    simd_float vSum_exp = simdf32_setzero();

    //Calculate max_zm
    for (size_t i = 0; i < VECSIZE_FLOAT; ++i) {
        max_zm = std::max(max_zm, vMax_zm[i]);
    }
    vMax_zm = simdf32_set(max_zm);

    //Calculate sum_exp
    for (size_t i = 0; i < rowSeqLen; ++i) {
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < colSeqLen_padding; j+= VECSIZE_FLOAT) {
            simd_float vZmForward = simdf32_load(&zm[i][j]);
            vZmForward = simdf32_exp(simdf32_sub(vZmForward, vMax_zm));
            vSum_exp = simdf32_add(vSum_exp, vZmForward);
        }
    }
    sum_exp += simdf32_hadd(vSum_exp);
}

template <bool profile>
void FwBwAligner::backward()  {
    //Init zInit
    size_t vecsize_float = static_cast<size_t>(VECSIZE_FLOAT);
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, colSeqLen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));

        memcpy(zmFirst + 1, zInit[0], rowSeqLen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], rowSeqLen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], rowSeqLen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= rowSeqLen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (profile) {
                    vScoreMatrix = simdf32_load(&scoreBackwardProfile_exp[rowSeqAANum[rowSeqLen - i]][start + j - 1]);
                } else {
                    size_t reverse_i = rowSeqLen - i;
                    size_t reverse_j = colSeqLen - start - j + 1;
                    simd_float vScoreBackward;
                    if (reverse_j >= vecsize_float) {
                        vScoreBackward = simdf32_loadu(&scoreForward[reverse_i][reverse_j - vecsize_float]);
                    } else {
                        size_t elements_to_fill = reverse_j;
                        vScoreBackward = simdf32_set(FLT_MIN_EXP);
                        // fill from the back
                        for (size_t k = 0; k < elements_to_fill; ++k) {
                            vScoreBackward[VECSIZE_FLOAT - elements_to_fill + k] = scoreForward[reverse_i][k];
                        }
                    }
                    vScoreBackward = simdf32_reverse(vScoreBackward);
                    vScoreMatrix = simdf32_exp(vScoreBackward);
                }
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            zmMaxRowBlock = simdf32_hmax(vZmax_tmp);
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);

            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, simdf32_set(exp_go)),
                                        simdf32_mul(vZf, simdf32_set(exp_ge))
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            size_t adjusted_memcpycols = memcpy_cols - memcpy_cols % VECSIZE_FLOAT;
            size_t forwardBlockStart = colSeqLen - start;
            simd_float vCurrMax = simdf32_set(current_max);
            for (size_t j = 1; j <= adjusted_memcpycols; j += VECSIZE_FLOAT) {
                forwardBlockStart -= vecsize_float;
                size_t simd_index = forwardBlockStart;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), vCurrMax);

                simd_float vZmForward = simdf32_loadu(&zm[rowSeqLen - i][simd_index]);

                simd_float vZmCurr_reverse = simdf32_reverse(vZmCurr);
                simd_float vZmForward_Backward = simdf32_add(vZmForward, vZmCurr_reverse);
                simdf32_storeu(&zm[rowSeqLen - i][simd_index], vZmForward_Backward);
            }

            // Handle remainder
            if (memcpy_cols != length) {
                size_t remainder = memcpy_cols % VECSIZE_FLOAT;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[adjusted_memcpycols+1]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[adjusted_memcpycols+1], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), simdf32_set(current_max));
                for (size_t k = 0; k < remainder; ++k) {
                    size_t j_index = remainder - k - 1;
                    zm[rowSeqLen - i][j_index] += vZmCurr[k];
                }
            }

#if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#endif

            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, vCurrMax);

            zInit[0][i-1] = vNextZinit[0]; zInit[1][i-1] = vNextZinit[1]; zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            
            if (i < rowSeqLen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max     
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
    }
}

template<>
// profile: true for profile, false for query*target matrix
// backtrace: 0 for no backtrace, 1 for local alignment, 2 for semi-global alignment, 3 for global alignment. Default: <true, 1>
void FwBwAligner::runFwBw<true, 0>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
}

template<>
void FwBwAligner::runFwBw<false, 0>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
}

template<>
void FwBwAligner::runFwBw<true, 1>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<1>();
}

template<>
void FwBwAligner::runFwBw<false, 1>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<1>();
}

template<>
void FwBwAligner::runFwBw<true, 2>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<2>();
}

template<>
void FwBwAligner::runFwBw<false, 2>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<2>();
}

template<>
void FwBwAligner::runFwBw<true, 3>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<3>();
}

template<>
void FwBwAligner::runFwBw<false, 3>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<3>();
}

template<bool profile>
void FwBwAligner::computeProbabilityMatrix() {
    float logsumexp_zm = max_zm + log(sum_exp);

    simd_float vLogsumexp_zm = simdf32_set(logsumexp_zm);
    size_t colLoopCount = colSeqLen / VECSIZE_FLOAT; 
    size_t colLoopEndPos = colLoopCount * VECSIZE_FLOAT;

    P = zm; //reuse zm matrix to store the probability matrix
    maxP = 0.0;
    simd_float vMaxP = simdf32_setzero();
    for (size_t i = 0; i < rowSeqLen; ++i) {
        // Fill the probability matrix
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < colLoopEndPos; j += VECSIZE_FLOAT) {
            simd_float vZmForward_Backward = simdf32_load(&zm[i][j]);
            simd_float scoreForwardVal;
            if (profile) {
                scoreForwardVal = simdf32_load(&scoreForwardProfile[rowSeqAANum[i]][j]);
            } else {
                scoreForwardVal = simdf32_load(&scoreForward[i][j]);
            }
            // simd_float scoreForwardVal = simdf32_load(&scoreForward[targetNum[i]][j]);
            simd_float P_val = simdf32_exp(simdf32_sub(vZmForward_Backward, simdf32_add(scoreForwardVal, vLogsumexp_zm)));
            simdf32_store(&P[i][j], P_val);
            vMaxP = simdf32_max(vMaxP, P_val);
        }    
        for (size_t j = colLoopEndPos; j < colSeqLen; ++j) {
            if (profile) {
                P[i][j] = exp(zm[i][j] - scoreForwardProfile[rowSeqAANum[i]][j] - logsumexp_zm);
            } else {
                P[i][j] = exp(zm[i][j] - scoreForward[i][j] - logsumexp_zm);
            }
            // P[i][j] = exp(zm[i][j] - scoreForward[targetNum[i]][j] - logsumexp_zm);
            maxP = std::max(maxP, P[i][j]);
        }
    }

    // Calculate the maximum probability
    for (size_t k = 0; k < VECSIZE_FLOAT; ++k) {
        maxP = std::max(maxP, vMaxP[k]);
    }

}

template<int backtrace>
void FwBwAligner::computeBacktrace() {
    // MAC algorithm from HH-suite
    uint8_t val;
    size_t max_i = 0;
    size_t max_j = 0;
    float term1, term2, term3, term4 = 0.0f;
    float score_MAC = -std::numeric_limits<float>::max();

    memset(S_curr, 0, (colSeqLen + 1) * sizeof(float));

    switch (backtrace) {
        case 1: // local
            memset(S_prev, 0, (colSeqLen + 1) * sizeof(float));
            break;
        case 2: // semiglobal
            memset(S_prev, 0, (colSeqLen + 1) * sizeof(float));
            break;
        case 3: // global
            std::fill(S_prev, S_prev + colSeqLen + 1, -std::numeric_limits<float>::max());
            S_prev[0] = 0.0;
            S_curr[0] = -std::numeric_limits<float>::max();
            break;
    }


    for (size_t i = 0; i <= rowSeqLen; ++i) {
        btMatrix[i][0] = States::STOP;
    }
    for (size_t j = 0; j <= colSeqLen; ++j) {
        btMatrix[0][j] = States::STOP;
    }

    for (size_t i = 1; i <= rowSeqLen; ++i) {
        for (size_t j = 1; j <= colSeqLen; ++j) {
            term1 = P[i - 1][j - 1] - mact; // STOP
            term2 = S_prev[j - 1] + P[i - 1][j - 1] - mact; // M
            term4 = S_prev[j] - 0.5 * mact; // D
            term3 = S_curr[j - 1] - 0.5 * mact; // I
            calculate_max4(S_curr[j], term1, term2, term3, term4, val);
            btMatrix[i][j] = val;
            
            switch (backtrace) {
                case 1: // local
                    if (S_curr[j] > score_MAC) {
                        max_i = i;
                        max_j = j;
                        score_MAC = S_curr[j];
                    }
                    break;
                case 2: // semiglobal
                    if ((i == rowSeqLen || j == colSeqLen) && S_curr[j] > score_MAC) { // only calculate for last column if j is last column
                        max_i = i;
                        max_j = j;
                        score_MAC = S_curr[j];
                    }
                    break;
                case 3: // global
                    break;
            }
        }
        std::swap(S_prev, S_curr);
        if (backtrace == 3 && i==1) {
            S_prev[0] = -std::numeric_limits<float>::max();
        }
    }
    // Set max_i and max_j for global alignment
    if (backtrace == 3) {
        max_i = rowSeqLen;
        max_j = colSeqLen;
        score_MAC = S_curr[colSeqLen];
    }
    // traceback 
    alignResult = {};
    alignResult.cigar = "";
    alignResult.cigar.reserve(colSeqLen + rowSeqLen);
    alignResult.score1 = maxP;
    alignResult.score2 = score_MAC;

    alignResult.qEndPos1 = max_j - 1;
    alignResult.dbEndPos1 = max_i - 1;
    uint32_t aaIds = 0;
    bool exitLoop = false;
    while (max_i > 0 && max_j > 0 && !exitLoop) {
        uint8_t state = btMatrix[max_i][max_j];
        switch (state) {
            case States::M:
                --max_i;
                --max_j;
                alignResult.qStartPos1 = max_j;
                alignResult.dbStartPos1 = max_i;
                alignResult.cigar.push_back('M');
                aaIds += (rowSeqAANum[max_i] == colSeqAANum[max_j]);
                break;

            case States::I:
                --max_j;
                alignResult.cigar.push_back('I');
                break;

            case States::D:
                --max_i;
                alignResult.cigar.push_back('D');
                break;

            default:
                exitLoop = true; 
                break;
        }
    }
    while (!alignResult.cigar.empty() && alignResult.cigar.back() != 'M') {
        alignResult.cigar.pop_back();
    }
    alignResult.cigarLen = alignResult.cigar.length();
    std::reverse(alignResult.cigar.begin(), alignResult.cigar.end());
    alignResult.identicalAACnt = aaIds;
    alignResult.qCov = SmithWaterman::computeCov(alignResult.qStartPos1, alignResult.qEndPos1, colSeqLen);
    alignResult.dbCov = SmithWaterman::computeCov(alignResult.dbStartPos1, alignResult.dbEndPos1, rowSeqLen);
}

FwBwAligner::s_align FwBwAligner::getFwbwAlnResult() {
    return alignResult;
}

int fwbw(int argc, const char **argv, const Command &command) {
    //Prepare the parameters & DB
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> alnRes (par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    alnRes.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter fwbwAlnWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    fwbwAlnWriter.open();
    const int querySeqType = qdbr.getDbtype();
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        Debug(Debug::ERROR) << "Invalid datatype. Nucleotide.\n";
        EXIT(EXIT_FAILURE);
    }
    SubstitutionMatrix subMat = SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias); // Check : par.scoreBias = 0.0
    
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(alnRes.getSize()) / static_cast<double>(flushSize)));
    Debug(Debug::INFO) << "Processing " << iterations << " iterations\n";
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(alnRes.getSize() - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);

#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            size_t length = par.blocklen;
            Sequence qSeq(par.maxSeqLen, qdbr.getDbtype(), &subMat, 0, false, false);
            Sequence dbSeq(par.maxSeqLen, tdbr.getDbtype(), &subMat, 0, false, false);
            
            const size_t assignSeqLen = VECSIZE_FLOAT * sizeof(float) * 20; 
            FwBwAligner fwbwaligner(subMat, -par.fwbwGapopen, -par.fwbwGapextend, par.temperature, par.mact, assignSeqLen, assignSeqLen, length, true);
            char entrybuffer[1024 + 32768*4];
            std::string alnResultsOutString;
            char buffer[1024 + 32768*4];
            std::vector<Matcher::result_t> localFwbwResults;
            localFwbwResults.reserve(300);

#pragma omp for schedule(dynamic,1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                unsigned int key = alnRes.getDbKey(id);
                const size_t queryId = qdbr.getId(key);
                char *alnData = alnRes.getData(id, thread_idx);
                localFwbwResults.clear();

                const char* querySeq = qdbr.getData(queryId, thread_idx);
                size_t queryLen = qdbr.getSeqLen(queryId);

                qSeq.mapSequence(queryId, key, querySeq, queryLen);
                fwbwaligner.initProfile(qSeq.numSequence, queryLen);
                fwbwAlnWriter.writeStart(thread_idx);

                while (*alnData != '\0'){
                    Util::parseKey(alnData, entrybuffer);
                    unsigned int targetKey = (unsigned int) strtoul(entrybuffer, NULL, 10);
                    const size_t targetId = tdbr.getId(targetKey);
                    const char* targetSeq = tdbr.getData(targetId, thread_idx);
                    size_t targetLen = tdbr.getSeqLen(targetId);

                    dbSeq.mapSequence(targetId, targetKey, targetSeq, targetLen);
                    //Init target & Resizing memory
                    fwbwaligner.initAlignment(dbSeq.numSequence, targetLen, queryLen); 
                    switch(par.fwbwBacktraceMode) {
                        case 0: fwbwaligner.runFwBw<true,0>(); break;
                        case 1: fwbwaligner.runFwBw<true,1>(); break;
                        // case 2: fwbwaligner.runFwBw<true,2>(); break; //hidden
                        // case 3: fwbwaligner.runFwBw<true,3>(); break; //hidden
                    }

                    // Map s_align values to result_t 
                    FwBwAligner::s_align fwbwAlignment = fwbwaligner.getFwbwAlnResult();
                    
                    float qcov = fwbwAlignment.qCov;
                    float dbcov = fwbwAlignment.dbCov;
                    float evalue = 0;
                    const int score = fwbwAlignment.score2;
                    const unsigned int qStartPos = fwbwAlignment.qStartPos1;
                    const unsigned int dbStartPos = fwbwAlignment.dbStartPos1;
                    const unsigned int qEndPos = fwbwAlignment.qEndPos1;
                    const unsigned int dbEndPos = fwbwAlignment.dbEndPos1;
                    std::string backtrace = fwbwAlignment.cigar;
                    unsigned int alnLength = backtrace.size();
                    float seqId = Util::computeSeqId(par.seqIdMode, fwbwAlignment.identicalAACnt, queryLen, targetLen, alnLength);
                    Matcher::result_t res = Matcher::result_t(targetKey, score, qcov, dbcov, seqId, evalue, alnLength, qStartPos, qEndPos, queryLen, dbStartPos, dbEndPos, targetLen, backtrace);
                    if (Alignment::checkCriteria(res, 0, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        localFwbwResults.emplace_back(res);
                    }
                    alnData = Util::skipLine(alnData);
                }

                // sort local results. They will currently be sorted by first fwbwscore, then targetlen, then by targetkey.
                SORT_SERIAL(localFwbwResults.begin(), localFwbwResults.end(), Matcher::compareHits);
                for (size_t result = 0; result < localFwbwResults.size(); result++) {
                    size_t len = Matcher::resultToBuffer(buffer, localFwbwResults[result], true, true);
                    alnResultsOutString.append(buffer, len);
                }
                fwbwAlnWriter.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), alnRes.getDbKey(id), thread_idx);
                alnResultsOutString.clear();
                localFwbwResults.clear();            
            }
        }
        alnRes.remapData();
        
    }
    fwbwAlnWriter.close();
    alnRes.close();
    qdbr.close();
    tdbr.close();

    return EXIT_SUCCESS;
}