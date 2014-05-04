#ifndef QUERY_SCORE_H
#define QUERY_SCORE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Calculates the overall prefiltering score for the template database sequences and returns all sequences 
// with the prefiltering score >= prefiltering threshold.
//


#include <emmintrin.h>
#include <mmintrin.h>

#include <stdlib.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <bitset>

#include "../commons/Debug.h"

typedef struct {
    size_t seqId;
    float zScore;
    unsigned short prefScore;
} hit_t;


class QueryScore {
    public:

        QueryScore (int dbSize, unsigned short * seqLens, int k, short kmerThr, float kmerMatchProb, float zscoreThr);

        virtual ~QueryScore ();

        // add k-mer match score for all DB sequences from the list
        virtual void addScores (int* seqList, int seqListSize, unsigned short score) = 0;

        void setPrefilteringThresholds();

        void setPrefilteringThresholdsRevSeq();

        float getZscore(int seqPos);

       // get the list of the sequences with the score > z-score threshold 
        std::pair<hit_t *, size_t> getResult (int querySeqLen, unsigned int identityId);

        // reset the prefiltering score counter for the next query sequence
        virtual void reset () = 0;

        void printScores();

        // maximal resultList
        static const size_t MAX_RES_LIST_LEN = 150000;

    private:
        static bool compareHits(hit_t first, hit_t second);

        short sse2_extract_epi16(__m128i v, int pos);

        void printVector(__m128i v);

        int counter;

        short kmerThr;

        double kmerMatchProb;

        float s_per_match;

        float s_per_pos;

        float zscore_thr;


    protected:
        inline unsigned short sadd16(unsigned short a, unsigned short b)
        {
            unsigned int s = (unsigned int)(a+b);
            return -(s>>16) | (unsigned short)s;
        }

        inline short sadd16_signed(short x, short y)
        {   
            unsigned short ux = x;
            unsigned short uy = y;
            unsigned short res = ux + uy;

            /* Calculate overflowed result. (Don't change the sign bit of ux) */
            ux = (ux >> 15) + SHRT_MAX;

            /* Force compiler to use cmovns instruction */
            if ((short) ((ux ^ uy) | ~(uy ^ res)) >= 0)
            {   
                res = ux;
            }

            return res;
        }

        inline short ssub16_signed (short x, short y)
        {
            unsigned short ux = x;
            unsigned short uy = y;
            unsigned short res = ux - uy;

            ux = (ux >> 15) + SHRT_MAX;

            /* Force compiler to use cmovns instruction */
            if ((short)((ux ^ uy) & (ux ^ res)) < 0)
            {
                res = ux;
            }

            return res;
        }

        // size of the database in scores_128 vector (the rest of the last _m128i vector is filled with zeros)
        int scores_128_size;
        // position in the array: sequence id
        // entry in the array: prefiltering score
        __m128i* scores_128;
        unsigned short  * scores;

        __m128i* thresholds_128;
        unsigned short  * thresholds;

        // float because it is needed for statistical calculations
        float * seqLens;
        float seqLenSum;

        int* steps;
        int nsteps;

        size_t scoresSum;

        int numMatches;

        float matches_per_pos;

        // number of sequences in the target DB
        int dbSize;

        // list of all DB sequences with the prefiltering score > z-score threshold with the corresponding scores
        hit_t * resList;
};

#endif
