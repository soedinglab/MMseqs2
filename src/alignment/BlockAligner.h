#ifndef BLOCKALIGNER_H
#define BLOCKALIGNER_H

#include "block_aligner.h"
#include "StripedSmithWaterman.h"

class EvalueComputation;

class BlockAligner {
public:
    BlockAligner(
        int dbtype,
        size_t maxSequenceLength,
        BaseMatrix *m, SubstitutionMatrix::FastMatrix* fastMatrix, EvalueComputation * evaluer,
        bool compBiasCorrection, float compBiasCorrectionScale,
        int8_t gapOpen, int8_t gapExtend
    );

    ~BlockAligner();

    void initQuery(Sequence* query);

    typedef struct UngappedAln_res {
        int bitScore;
        float qcov;
        float tcov;
        double eval;
        int alnLen;
        int qStart;
        int tStart;
        int qEnd;
        int tEnd;
        int diagonalLen;
        int score;
        int diagonal;
        UngappedAln_res() : bitScore(0), qcov(0), tcov(0), eval(0), alnLen(0), qStart(0), tStart(0), qEnd(0), tEnd(0), diagonalLen(0), score(0), diagonal(0) {}
        UngappedAln_res(int bitScore, float qcov, float tcov, double eval, int alnLen, size_t qStart, size_t tStart, size_t qEnd, size_t tEnd, int diagonalLen, int score, int diagonal) :
            bitScore(bitScore), qcov(qcov), tcov(tcov), eval(eval), alnLen(alnLen), qStart(qStart), tStart(tStart), qEnd(qEnd), tEnd(tEnd), diagonalLen(diagonalLen), score(score), diagonal(diagonal) {}
    } UngappedAln_res;

    BlockAligner::UngappedAln_res ungappedAlign(
        Sequence* target,
        const unsigned short diagonal
    );
    BlockAligner::UngappedAln_res hammingDistance(
        Sequence* target, const unsigned short diagonal);



    s_align align(
        Sequence* currentTarget,
        size_t qStart,
        size_t tStart,
        std::string& backtrace,
        int xdrop,
        float covThr,
        int covMode
    );

    s_align gappedLocalAlign(
        Sequence* currentTarget,
        int q_idx, int t_idx,
        Cigar* cigar, int32_t x_drop,
        float covThr,
        int covMode
    );

    s_align bandedalignForward(
        Sequence* currentTarget,
        size_t qIdx,
        size_t tIdx,
        std::string& backtrace,
        int xdrop
    );

    s_align bandedalign(
        Sequence* currentTarget,
        size_t qIdx,
        size_t tIdx,
        std::string& backtrace,
        int xdrop,
        float covThr,
        int covMode
    );

    s_align bandedalignBackward(
        Sequence* currentTarget,
        size_t qIdx,
        size_t tIdx,
        std::string& backtrace,
        int xdrop
    );

    s_align gappedLocalAlignBackward(
        Sequence* currentTarget,
        int qIdx, int tIdx,
        Cigar* cigar, int32_t x_drop
    );
    s_align gappedLocalAlignForward(
        Sequence* currentTarget,
        int qIdx, int tIdx,
        Cigar* cigar, int32_t x_drop
    );

   
private:
    size_t maxSequenceLength;    
// holds values of the current active query
    Sequence * currentQuery;

    PaddedBytes* query;
    PaddedBytes* target;
    PosBias* queryBias;
    PosBias* targetBias;
    AAMatrix* matrix;
    BlockHandle blockTrace;
    BlockHandle blockNoTrace;
    Cigar* cigar;
    Gaps gaps;
    SizeRange range;

    bool compBiasCorrection;
    // weight for the correlation score, if set to 0.0 it is turned off
    float compBiasCorrectionScale;

    int dbtype;
    // costs to open and extend a gap

    // substitution matrix
    BaseMatrix* subMat=nullptr;
    SubstitutionMatrix::FastMatrix* fastMatrix ;
    // evalue
    EvalueComputation * evaluer;
    
    const char* querySeq;
    const unsigned char * queryNumSeq;
    int queryLength;
    int querySeqType;
    int16_t* queryCompBias;
    int16_t* targetCompBias;
    int16_t* queryCompBiasRevArr;
    float *tmpCompBias;
    int8_t* queryRevNumSeq;
    int16_t* queryCompBiasRev;
    
};

#endif
