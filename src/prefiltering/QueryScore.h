#ifndef QUERY_SCORE_H
#define QUERY_SCORE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Calculates the overall prefiltering score for the template database sequences and returns all sequences 
// with the prefiltering score >= prefiltering threshold.
//


#include <emmintrin.h>

#include <stdlib.h>
#include <list>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <math.h>
#include <fstream>
#include <sstream>

typedef struct {
    int seqId;
    float prefScore;
    float eval;
} hit_t;


class QueryScore {
    public:

        QueryScore (int dbSize, unsigned short * seqLens, int k);

        ~QueryScore ();

        // add k-mer match score for all DB sequences from the list
        virtual void addScores (int* seqList, int seqListSize, unsigned short score) = 0;

        // get the list of the sequences with the score > z-score threshold and the corresponding 
        std::list<hit_t>* getResult (int querySeqLen);

        // reset the prefiltering score counter for the next query sequence
        virtual void reset () = 0;

        int numMatches;


    private:
        void setPrefilteringThresholds();

        static bool compareHitList(hit_t first, hit_t second);

        void addElementToResults(int seqId);

        unsigned short sse2_extract_epi16(__m128i v, int pos);

        void printVector(__m128i v);

        int counter;

        std::ofstream s_per_pos_file;
        std::ofstream second_term_file;
        std::ofstream norm_score1_file;
        std::ofstream zscore_file;


    protected:
        inline unsigned short sadd16(unsigned short a, unsigned short b)
        {
            unsigned int s = (unsigned int)(a+b);
            return -(s>>16) | (unsigned short)s;
        }

        // position in the array: sequence id
        // entry in the array: prefiltering score
        __m128i* scores_128;
        unsigned short  * scores;

        __m128i* thresholds_128;
        unsigned short  * thresholds;

        __m128i* seqLens_128;
        unsigned short * seqLens;

        int* steps;
        int nsteps;

        size_t seqLenSum;

        size_t scoresSum;

        // number of sequences in the target DB
        int dbSize;

        // list of all DB sequences with the prefiltering score > z-score threshold with the corresponding scores
        std::list<hit_t>* resList;

};

#endif
