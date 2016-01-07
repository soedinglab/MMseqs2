#ifndef QUERY_SCORE_H
#define QUERY_SCORE_H

// Written by Martin Steinegger martin.steinegger@campus.lmu.de & Maria Hauser mhauser@genzentrum.lmu.de
//
// Calculates the overall prefiltering score for the template database sequences and returns all sequences
// with the prefiltering score >= prefiltering threshold.
//
#include <stdlib.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/cdefs.h>
#include <zconf.h>
#include <stdlib.h>     /* abs */

#include "Debug.h"
#include "Util.h"
#include "simd.h"
#include "IndexTable.h"


typedef struct {
    unsigned int seqId;
    float pScore;
    unsigned short diagonal;
    unsigned char prefScore;
} hit_t;


class QueryScore {
public:

    QueryScore(size_t dbSize, unsigned int *seqLens, unsigned int seedLength, short kmerThr, float kmerMatchProb);

    virtual ~QueryScore ();

    inline void addScoresLocal(IndexEntryLocal *seqList, const unsigned short i, const int unsigned seqListSize) {
        unsigned char * data =  (unsigned char *) scores;
        for (unsigned int seqIdx = 0; LIKELY(seqIdx < seqListSize); seqIdx++){
            const IndexEntryLocal entry = seqList[seqIdx];
            const unsigned int seqIndex = entry.seqId * 2;
            const unsigned short currDiagonal = i - entry.position_j;
            const unsigned char dbDiagonal = data[seqIndex];
            const unsigned char oldScore   = data[seqIndex + 1];
            const unsigned char scoreToAdd = (UNLIKELY(currDiagonal == dbDiagonal ) && LIKELY(oldScore < 255)) ? 1 : 0;
            if(UNLIKELY(currDiagonal == dbDiagonal )){
                std::cout << "seq=" << seqIndex << "\ti=" << i  << "\tj=" << entry.position_j << "\tdiag=" << (unsigned short) currDiagonal << std::endl;
            }
            const unsigned char newScore = oldScore + scoreToAdd;
            data[seqIndex]     = currDiagonal;
            data[seqIndex + 1] = newScore;
        }
        numMatches += seqListSize;
    }

//    inline void addScoresLocal (IndexEntryLocal * __restrict seqList, const unsigned short i,
//                                const int seqListSize, unsigned short score){
//        for (unsigned int seqIdx = 0; LIKELY(seqIdx < seqListSize); seqIdx++){
//            IndexEntryLocal entry = seqList[seqIdx];
//            const unsigned short j = entry.position_j;
//            const unsigned int seqId = entry.seqId;
//            const unsigned char diagonal = i - j + 256;
//            if (UNLIKELY(diagonal == scores[seqId])){
//                if(UNLIKELY(scoreSizes >= MAX_LOCAL_RESULT_SIZE)){
//                    break;
//                }
//                LocalResult * currLocalResult = localResults + scoreSizes++;
//                currLocalResult->seqId = seqId;
//                currLocalResult->score = score;
//            } else {
//                scores[seqId] = diagonal;
//            }
//        }
//        numMatches += seqListSize;
//
//    }

    // add k-mer match score for all DB sequences from the list
    inline void addScores (unsigned int* __restrict seqList, int seqListSize, unsigned short score){
        for (int i = 0; i < seqListSize; i++){
            const int seqId = seqList[i];
            scores[seqId] = Util::sadd16(scores[seqId], score);
        }
        scoresSum += score * seqListSize;
        numMatches += seqListSize;
    }

    // reset the prefiltering score counter for the next query sequence
    virtual void reset () = 0;

    // prints the score results
    void printScores();

    // returns the number of Matches for this Query
    unsigned int getNumMatches(){
        return numMatches;
    };

    // returns the current local Result size
    size_t getLocalResultSize(){
        size_t retValue = 0;
        for(size_t i = 1; i < SCORE_RANGE; i++){
            retValue += scoreSizes[i] * i;
        }
        return retValue;
    }

    // maximal resultList
    static const size_t MAX_RES_LIST_LEN = 150000;

    short sse2_extract_epi16(__m128i v, int pos);

    void printVector(__m128i v);

    // sorting functions for the hitlist
    static bool compareHitsByPValue(hit_t first, hit_t second);
    static bool compareHitsByDiagonalScore(hit_t first, hit_t second);


    const static size_t SCORE_RANGE = 256;

    // score distribution of current query
    unsigned int *scoreSizes;


protected:

    const unsigned int SIMD_SHORT_SIZE = VECSIZE_INT * 2;  // *2 for short

    // size of the database in scores_128 vector (the rest of the last _m128i vector is filled with zeros)
    int scores_128_size;

    short kmerThr;

    double kmerMatchProb;

    // position in the array: sequence id
    // entry in the array: prefiltering score
    simd_int* __restrict scores_128;

    unsigned short  * __restrict scores;

    size_t numMatches;
    size_t distanceCnt;
    size_t normalDistCnt;
    // number of sequences in the target DB
    size_t dbSize;

    // list of all DB sequences with the prefiltering score > z-score threshold with the corresponding scores
    hit_t * resList;

    // sum up the scores over the query
    size_t scoresSum;

    // float because it is needed for statistical calculations
    float * seqLens;

    // total sum of amino acid minus k
    float seqLenSum;

};

#endif
