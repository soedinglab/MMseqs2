#ifndef FWBW
#define FWBW

#include "SubstitutionMatrix.h"
#include "IndexReader.h"
#include "DBReader.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>

class FwBwAligner {
public:
    typedef struct {
        float score1;
        float score2;
        int32_t dbStartPos1;
        int32_t dbEndPos1;
        int32_t	qStartPos1;
        int32_t qEndPos1;
        int32_t ref_end2;
        float qCov;
        float dbCov;
        std::string cigar;
        double evalue;
        int identicalAACnt;
        int32_t cigarLen;
        int word;
    } s_align;


    FwBwAligner(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity = 320, size_t colsCapacity = 320, size_t length = 16, int backtrace = 0);
    FwBwAligner(float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity = 320, size_t colsCapacity = 320, size_t length = 16,  int backtrace = 0);

    ~FwBwAligner();

    s_align align();
    size_t getRowsCapacity() const { return rowsCapacity; }
    size_t getColsCapacity() const { return colsCapacity; }
    size_t getBlockLength() const { return length; }
    
    //Reallocation & Resizing
    void reallocateProfile(size_t newColsCapacity);
    template<bool profile, int backtrace>
    void resizeMatrix(size_t newRowsCapacity, size_t newColsCapacity);
    void resetParams(float newGapOpen, float newGapExtend, float newTemperature);
    
    //Initilization
    void initProfile(unsigned char* colAANum, size_t colAALen);
    void initAlignment(unsigned char* targetNum, size_t targetLen, size_t queryLen);
    void initScoreMatrix(float** inputScoreMatrix, int* gaps);

    float getProbability(int i, int j) const {
        return P[i][j];
    }

    float** getProbabiltyMatrix() {
        return P;
    }

    unsigned char* colSeqAANum = nullptr;
    unsigned char* rowSeqAANum = nullptr;

    template<bool profile, int backtrace>
    void runFwBw();

    s_align getFwbwAlnResult();

    float maxP;
    float temperature;

private:

    size_t length;
    float** zm;
    float** P;
    float* zmFirst;
    float* zeFirst;
    float* zfFirst;

    float** zInit;
    // Score Matrix
    float** scoreForwardProfile = nullptr; // profile true
    float** scoreForwardProfile_exp = nullptr; // profile true
    float** scoreBackwardProfile_exp = nullptr; // profile true

    float** scoreForward = nullptr; // profile false
    //backtrace
    uint8_t** btMatrix = nullptr; // backtrace true

    float* zmBlockPrev;
    float* zmBlockCurr;
    float* zeBlock;
    float* zfBlock;

    float* vj;
    float* wj;
    float* exp_ge_arr;

    // const SubstitutionMatrix & subMat;
    float gapOpen;
    float gapExtend;
    float mact; // backtrace true

    // float temperature;
    size_t rowsCapacity;
    size_t colsCapacity;

    size_t blockCapacity;
    size_t blocks;
    size_t colSeqLen;
    size_t rowSeqLen;

    float** blosum= nullptr; // Profile true
    float* S_prev = nullptr; // backtrace true
    float* S_curr = nullptr; // backtrace true
    float exp_go;
    float exp_ge;
    float max_zm;
    float sum_exp;
    // float maxP;
    
    
    //simd_float vSum_exp;
    size_t colSeqLen_padding;
    
    s_align alignResult;

    template<bool profile>
    void forward();

    template<bool profile>
    void backward();

    template<bool profile>
    void computeProbabilityMatrix();

    template<int backtrace>
    void initBacktraceMatrix();

    template<int backtrace>
    void computeBacktrace();
};


#endif //FWBW_H