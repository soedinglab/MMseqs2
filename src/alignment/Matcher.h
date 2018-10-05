#ifndef MATCHER_H
#define MATCHER_H

//
// Written by Martin Steinegger & Maria Hauser, mhauser@genzentrum.lmu.de
//
// Calls SSE2 parallelized calculation of Smith-Waterman alignment and non-parallelized traceback afterwards.
//

#include <cfloat>
#include <algorithm>
#include <vector>

#include "Sequence.h"
#include "BaseMatrix.h"
#include "StripedSmithWaterman.h"
#include "EvalueComputation.h"
#include "BandedNucleotideAligner.h"

class Matcher{

public:
    static const unsigned int SCORE_ONLY = 0;
    static const unsigned int SCORE_COV = 1;
    static const unsigned int SCORE_COV_SEQID = 2;
    const static int ALN_RES_WITH_OUT_BT_COL_CNT = 10;

    const static int ALN_RES_WITH_BT_COL_CNT = 11;

    struct result_t {
        unsigned int dbKey;
        int score;
        float qcov;
        float dbcov;
        float seqId;
        double eval;
        unsigned int alnLength;
        int qStartPos;
        int qEndPos;
        unsigned int qLen;
        int dbStartPos;
        int dbEndPos;
        unsigned int dbLen;
        std::string backtrace;
        result_t(unsigned int dbkey,int score,
                 float qcov, float dbcov,
                 float seqId, double eval,
                 unsigned int alnLength,
                 int qStartPos,
                 int qEndPos,
                 unsigned int qLen,
                 int dbStartPos,
                 int dbEndPos,
                 unsigned int dbLen,
                 std::string backtrace) : dbKey(dbkey), score(score), qcov(qcov),
                                          dbcov(dbcov), seqId(seqId), eval(eval), alnLength(alnLength),
                                          qStartPos(qStartPos), qEndPos(qEndPos), qLen(qLen),
                                          dbStartPos(dbStartPos), dbEndPos(dbEndPos), dbLen(dbLen),
                                          backtrace(backtrace) {};
        result_t(){};
    };

    Matcher(int querySeqType, int maxSeqLen, BaseMatrix *m,
            EvalueComputation * evaluer, bool aaBiasCorrection,
            int gapOpen, int gapExtend);

    ~Matcher();

    // run SSE2 parallelized Smith-Waterman alignment calculation and traceback
    result_t getSWResult(Sequence* dbSeq, const int diagonal, const int covMode, const float covThr, const double evalThr,
                         unsigned int alignmentMode, unsigned int seqIdMode, bool isIdentical);

    // need for sorting the results
    static bool compareHits (const result_t &first, const result_t &second){
        //return (first.eval < second.eval);
        if(first.eval < second.eval )
            return true;
        if(second.eval < first.eval )
            return false;
        if(first.score > second.score )
            return true;
        if(second.score > first.score )
            return false;
        return false;
    }

    // map new query into memory (create queryProfile, ...)
    void initQuery(Sequence* query);

    static result_t parseAlignmentRecord(char *data, bool readCompressed=false);

    static void readAlignmentResults(std::vector<result_t> &result, char *data, bool readCompressed = false);

    static float estimateSeqIdByScorePerCol(uint16_t score, unsigned int qLen, unsigned int tLen);

    static std::string compressAlignment(const std::string &bt);

    static std::string uncompressAlignment(const std::string &cbt);


    static size_t resultToBuffer(char * buffer, const result_t &result, bool addBacktrace, bool compress  = true);

    static size_t computeAlnLength(size_t anEnd, size_t start, size_t dbEnd, size_t dbStart);


private:

    // costs to open a gap
    int gapOpen;
    // costs to extend a gap
    int gapExtend;

    // calculate the query queryProfile for SIMD registers processing 8 elements
    int maxSeqLen;

    // holds values of the current active query
    Sequence * currentQuery;

    // aligner Class
    SmithWaterman * aligner;
    // aligner for nucl
    BandedNucleotideAligner * nuclaligner;
    // substitution matrix
    BaseMatrix* m;
    // evalue
    EvalueComputation * evaluer;
    // byte version of substitution matrix
    int8_t * tinySubMat;
    // set substituion matrix
    void setSubstitutionMatrix(BaseMatrix *m);

};

#endif
