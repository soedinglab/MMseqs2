#ifndef BLOCKALIGNER_H
#define BLOCKALIGNER_H

#include "block_aligner.h"
#include "StripedSmithWaterman.h"

class EvalueComputation;

class BlockAligner {
public:
    BlockAligner(
        size_t maxSequenceLength,
        uintptr_t min, uintptr_t max,
        int8_t gapOpen, int8_t gapExtend,
        BaseMatrix& subMat,
        int dbtype = Parameters::DBTYPE_AMINO_ACIDS
    );

    ~BlockAligner();

    void initQuery(
        const char* querySeq,
        unsigned int queryLength,
        int querySeqType
    );

    typedef struct LocalAln {
        size_t a_start;
        size_t b_start;
        size_t a_end;
        size_t b_end;
        int32_t score;
        LocalAln() : a_start(0), b_start(0), a_end(0), b_end(0), score(0) {}
        LocalAln(size_t a_start, size_t b_start, size_t a_end, size_t b_end, int32_t score) :
            a_start(a_start), b_start(b_start), a_end(a_end), b_end(b_end), score(score) {}
    } LocalAln;

    LocalAln ungappedAlign(
        const char* targetData, size_t targetLen,
        int diagonal
    );

    s_align align(
        const char* targetSeq,
        unsigned int targetLength,
        int diagonal,
        std::string& backtrace,
        EvalueComputation *evaluer,
        int xdrop
    );

private:
    PaddedBytes* a;
    PaddedBytes* b;
    AAProfile* bProfile;
    AAMatrix* matrix;
    BlockHandle blockTrace;
    BlockHandle blockNoTrace;
    Cigar* cigar;
    uint32_t* sAlnCigar;

    SizeRange range;
    Gaps gaps;
    BaseMatrix& subMat;
    int dbtype;
    SubstitutionMatrix::FastMatrix fastMatrix;

    const char* querySeq;
    unsigned int queryLength;
    int querySeqType;
    std::string queryConsensus;
};

#endif
