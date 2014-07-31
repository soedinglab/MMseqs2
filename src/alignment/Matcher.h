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

#include "../commons/Sequence.h"
#include "../commons/BaseMatrix.h"
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

        /////////////////////////////////////////////////////////////////////////////////////
        // fast 2^x
        // ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
        // Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
        // Speed: 2.1E-8s (2.3E-8s) per call! (exp(): 8.5E-8, pow(): 1.7E-7)
        // Internal representation of float number according to IEEE 754: 
        //   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
//                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000 
//   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
/////////////////////////////////////////////////////////////////////////////////////
    inline float fpow2(float x)
    {
        if (x>FLT_MAX_EXP) return FLT_MAX;
        if (x<FLT_MIN_EXP) return 0.0f;
        int *px = (int*)(&x);                 // store address of float as pointer to long int
        float tx = (x-0.5f) + (3<<22);        // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127), 
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        int lx = *((int*)&tx) - 0x4b400000;   // integer value of x 
        float dx = x-(float)(lx);             // float remainder of x
        x = 1.0f + dx*(0.693019f             // polynomial approximation of 2^x for x in the range [0, 1]
                + dx*(0.241404f             // Gives relative deviation < 4.6E-6 
                    + dx*(0.0520749f            // Speed: 2.1E-8s
                        + dx* 0.0134929f )));
        *px += (lx<<23);                      // add integer power of 2 to exponent
        return x;
    }

    // query profile for SIMD registers
    unsigned short* queryProfileWord;

    BaseMatrix* m;

    void* workspace;

    void* workspace_memory;

};

#endif
