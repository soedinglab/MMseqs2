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
        typedef struct {
            std::string dbKey;
            int score;
            float qcov;
            float dbcov;
            float seqId;
            double eval;
        } result_t;

        Matcher(BaseMatrix* m, int maxSeqLen);

        ~Matcher();

        // run SSE2 parallelized Smith-Waterman alignment calculation and traceback
        result_t getSWResult(Sequence* query, Sequence* dbSeq, int seqDbSize);

        static bool compareHits (result_t first, result_t second){ if (first.score > second.score) return true; return false; }


    private:

        int maxAllocatedLen;

        void* H_workspace;

        void* E_workspace;

        void* F_workspace;

        // calculate the query profile for SIMD registers processing 8 elements
        void calcQueryProfileWord(Sequence* query);

        int maxSeqLen;

 
    // query profile for SIMD registers
    unsigned short* queryProfileWord;

    BaseMatrix* m;

    void* workspace;

    void* workspace_memory;

};

#endif
