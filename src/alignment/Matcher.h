#ifndef _H
#define _H

//
// Written by Maria Hauser, mhauser@genzentrum.lmu.de
//
// Calls SSE2 parallelized calculation of Smith-Waterman alignment and non-parallelized traceback afterwards.
//

#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <float.h>
#include <algorithm>

#include "Sequence.h"
#include "BaseMatrix.h"
#include "smith_waterman_sse2.h"

class Matcher{

    public:



    static const unsigned int SCORE_ONLY = 0;
    static const unsigned int SCORE_COV = 1;
    static const unsigned int SCORE_COV_SEQID = 2;

        struct result_t {
            std::string dbKey;
            int score;
            float qcov;
            float dbcov;
            float seqId;
            double eval;
            result_t(std::string dbkey,int score,
                 float qcov, float dbcov,
                 float seqId, double eval) : dbKey(dbkey), score(score), qcov(qcov), dbcov(dbcov), seqId(seqId), eval(eval) {};
        };

        Matcher(BaseMatrix* m, int maxSeqLen);

        ~Matcher();

        // run SSE2 parallelized Smith-Waterman alignment calculation and traceback
        result_t getSWResult(Sequence* dbSeq,const size_t seqDbSize,const double evalThr, const unsigned int mode);

        // need for sorting the results
        static bool compareHits (result_t first, result_t second){ if (first.score > second.score) return true; return false; }
    
        // map new query into memory (create profile, ...)
        void initQuery(Sequence* query);
    private:

        // calculate the query profile for SIMD registers processing 8 elements
        int maxSeqLen;
    
        // holds values of the current active query
        Sequence * currentQuery;
    
        // aligner Class
        SmithWaterman * aligner;
        // parameter for alignment
        const unsigned short GAP_OPEN = 10;
        const unsigned short GAP_EXTEND = 1;
        // substitution matrix
        BaseMatrix* m;
        // byte version of substitution matrix
        int8_t * tinySubMat;



};

#endif
