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

typedef struct {
    short seqId;
    short prefScore;
} hit_t;


class QueryScore {
    public:

        QueryScore (int dbSize, short prefThreshold);

        ~QueryScore ();

        // add k-mer match score for all DB sequences from the list
        void addScores (int* hitList, int hitListSize, short score);

        // increment the query position (only one similar k-mer match per db sequence and query position is added to the prefiltering score)
        void moveToNextQueryPos();

        // get the list of the sequences with the score > prefThreshold and the corresponding 
        std::list<hit_t>* getResult ();

        // reset the prefiltering score counter for the next query sequence
        void reset ();

    private:

        // number of sequences in the target DB
        int dbSize;

        // prefiltering threshold
        short prefThreshold;

        // position in the array: sequence id
        // entry in the array: prefiltering score
        short* scores;

        // mapping db id -> position of the last match in the query sequence
        short* lastMatchPos;

        // current position in the query sequence (prevents multiple matches for the same query position, different similar k-mers and the same db sequence)
        short currQueryPos;

        // sorted list of all DB sequences with the prefiltering score >= prefThreshold
        std::list<int>* hitList;

        // list of all DB sequences with the prefiltering score >= prefThreshold with the corresponding scores
        std::list<hit_t>* resList;

        void addElementToResults(int seqId);

};

#endif
