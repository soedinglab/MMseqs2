#ifndef QUERY_SCORE_H
#define QUERY_SCORE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
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
#include <unordered_map>
#include <vector>
#include <stdlib.h>     /* abs */

#include "Debug.h"
#include "Util.h"
#include "simd.h"
#include "IndexTable.h"


typedef struct {
    size_t seqId;
    float zScore;
    unsigned short prefScore;
} hit_t;

struct LocalResult{
    unsigned int seqId;
    unsigned short score;
    LocalResult(unsigned int seqId,
                unsigned short score) :
    seqId(seqId), score(score) {};
    LocalResult() : seqId(0), score(0) {};
}; // 4 + 2 = 6 byte => with padding 8 byte

class QueryScore {
public:
    
    QueryScore (int dbSize, unsigned short * seqLens, int k, short kmerThr, float kmerMatchProb, float zscoreThr);
    
    virtual ~QueryScore ();
    
    
    inline void addScoresLocal (IndexEntryLocal * __restrict seqList, const unsigned short i,
                                const int seqListSize, unsigned short score){
        
        //const unsigned short checkIfMatchedBefore = (i % 2) ? 0x8000 : 0x7FFF; // 1000000000000 0111111111111111
        for (int seqIdx = 0; LIKELY(seqIdx < seqListSize); seqIdx++){
            IndexEntryLocal entry = seqList[seqIdx];
            const unsigned short j = entry.position_j;
            const unsigned int seqId = entry.seqId;
            const unsigned short diagonal = i - j + 32768;
            if (UNLIKELY(diagonal == scores[seqId])){
                if(UNLIKELY(localResultSize >= MAX_LOCAL_RESULT_SIZE)){
                    break;
                }
                LocalResult * currLocalResult = localResults + localResultSize++;
                currLocalResult->seqId = seqId;
                currLocalResult->score = score;
            } else {
                scores[seqId] = diagonal;
            }
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

    // returns the current local Result size
    unsigned int getLocalResultSize(){
        return localResultSize;
    }
    
    // maximal resultList
    static const size_t MAX_RES_LIST_LEN = 150000;
    
    short sse2_extract_epi16(__m128i v, int pos);
    
    void printVector(__m128i v);
    
    static bool compareHits(hit_t first, hit_t second);
    
protected:

    const unsigned int SIMD_SHORT_SIZE = VECSIZE_INT * 2;  // *2 for short

    LocalResult * localResults;
    
    // current position in Localresults while adding Score
    unsigned int localResultSize;
    
    // max LocalResult size
    const unsigned int MAX_LOCAL_RESULT_SIZE = 100000000;
    
    // size of the database in scores_128 vector (the rest of the last _m128i vector is filled with zeros)
    int scores_128_size;
    
    short kmerThr;
    
    double kmerMatchProb;
    
    // position in the array: sequence id
    // entry in the array: prefiltering score
    simd_int* __restrict scores_128;
    
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
