#ifndef QUERY_SCORE_H
#define QUERY_SCORE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Calculates the overall prefiltering score for the template database sequences and returns all sequences 
// with the prefiltering score >= prefiltering threshold.
//

#include <list>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include "DynamicArray.h"

typedef struct {
    int seqId;
    float prefScore;
} hit_t;


class QueryScore {
    public:

        QueryScore (int dbSize, float prefThreshold);

        ~QueryScore ();

        // add k-mer match score for all DB sequences from the list
        void addScores (int* seqList, int seqListSize, short score);

        // increment the query position 
        //        void moveToNextQueryPos();

        // get the list of the sequences with the score > prefThreshold and the corresponding 
        std::list<hit_t>* getResult (int querySeqLen);

        // reset the prefiltering score counter for the next query sequence
        void reset ();

        void printStats();

    private:

        // number of sequences in the target DB
        int dbSize;

        // prefiltering threshold
        float prefThreshold;

        // position in the array: sequence id
        // entry in the array: prefiltering score
        int* scores;

        // sorted list of all DB sequences with the prefiltering score >= prefThreshold
        DynamicArray* hitList;

        // list of all DB sequences with the prefiltering score >= prefThreshold with the corresponding scores
        std::list<hit_t>* resList;

        double dbFractCnt;

        int qSeqCnt;

        void addElementToResults(int seqId);

};

#endif
