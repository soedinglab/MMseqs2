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
#include <math.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <stdlib.h>     /* abs */

#include "Debug.h"
#include "Util.h"
#include "IndexTable.h"


typedef struct {
    unsigned int seqId;
    float zScore;
    unsigned short prefScore;
} hit_t;

typedef struct {
    unsigned short diagonal;
    unsigned short prefScore;
    unsigned int seqId;
} LocalMatch;

//struct LocalResult{
//    unsigned int seqId;
//    unsigned short i;
//    unsigned short j;
//    unsigned short diagonal;
//    unsigned short score;
//
//    LocalResult(unsigned int seqId,
//                unsigned short i,
//                unsigned short j,
//                unsigned short diagonal,
//                unsigned short score) :
//    seqId(seqId), i(i), j(j), diagonal(diagonal), score(score) {};
//    LocalResult() : seqId(0), i(0), j(0), diagonal(0), score(0) {};
//
//}; // 2 + 2 + 2  + 2 + 4 = 12 byte => with padding 16 byte




class QueryScore {
public:
    
    QueryScore (int dbSize, unsigned short * seqLens, int k, short kmerThr,
                float kmerMatchProb, float zscoreThr);
    
    virtual ~QueryScore ();
    
    
    inline void addScoresLocal (IndexEntryLocal * __restrict seqList, const unsigned short i,
                                const int seqListSize, unsigned short score){
        
        for (int seqIdx = 0; seqIdx < seqListSize; seqIdx++){
            IndexEntryLocal entry = seqList[seqIdx];
            const unsigned short j = entry.position_j;
            const unsigned int seqId = entry.seqId;
            const unsigned short currDiagonal = i - j + 32768;
            //std::cout <<  i << " " << j << " " << entry.seqId << " " <<  diagonal << std::endl;
            LocalMatch * currMatch = matchArray + seqId;
            const unsigned short matchScore   = Util::sadd16(score, currMatch->prefScore);
            const unsigned short prevDiagonal = currMatch->diagonal;
            currMatch->diagonal   = currDiagonal;
            currMatch->prefScore  = (prevDiagonal == currDiagonal) ? matchScore : currMatch->prefScore;
            currMatch->seqId = seqId;
        }
        numMatches += seqListSize;

    }
    
    // add k-mer match score for all DB sequences from the list
    inline void addScores (unsigned int* __restrict seqList, int seqListSize, unsigned short score){
        for (int i = 0; i < seqListSize; i++){
            const int seqId = seqList[i];
            scores[seqId] = Util::sadd16(scores[seqId], score);
        }
        scoresSum += score * seqListSize;
        numMatches += seqListSize;
    }
    
    
    virtual void setPrefilteringThresholds() = 0;
    
    // get the list of the sequences with the score > z-score threshold
    virtual std::pair<hit_t *, size_t> getResult (int querySeqLen, unsigned int identityId) = 0;
    
    // reset the prefiltering score counter for the next query sequence
    virtual void reset () = 0;
    
    // prints the score results
    void printScores();
    
    // returns the number of Matches for this Query
    unsigned int getNumMatches(){
        return numMatches;
    };
    
    // maximal resultList
    static const size_t MAX_RES_LIST_LEN = 150000;
    
    short sse2_extract_epi16(__m128i v, int pos);

    void printVector(__m128i v);

    static bool compareHits(hit_t first, hit_t second);
    
protected:
//    std::unordered_map<unsigned int , unsigned short > localResult;
//    LocalResult * localResults;
//    max LocalResult size
//    const unsigned int MAX_LOCAL_RESULT_SIZE = 100000;
    
    // keeps track of diagnoals while Matching
    LocalMatch * matchArray;
    
    // current position in Localresults while adding Score
    unsigned int localResultSize;

    // size of the database in scores_128 vector (the rest of the last _m128i vector is filled with zeros)
    int scores_128_size;
    
    short kmerThr;
    
    double kmerMatchProb;
    
    // position in the array: sequence id
    // entry in the array: prefiltering score
    __m128i* __restrict scores_128;
    
    unsigned short  * __restrict scores;
    
    unsigned int numMatches;
        
    // number of sequences in the target DB
    int dbSize;
    
    // list of all DB sequences with the prefiltering score > z-score threshold with the corresponding scores
    hit_t * resList;
    
    // sum up the scores over the query
    size_t scoresSum;
    
    //
    float zscore_thr;


};

#endif
