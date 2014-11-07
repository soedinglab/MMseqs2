//
//  QueryScoreGlobal.h
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef QUERYSCOREGLOBAL_H
#define QUERYSCOREGLOBAL_H

#include "QueryScore.h"

class QueryScoreGlobal : public QueryScore {
    
public:
    QueryScoreGlobal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr);

    ~QueryScoreGlobal();
        
    void reset();
    
    void setPrefilteringThresholds();
    
    float getZscore(int seqPos);
    
    std::pair<hit_t *, size_t> getResult (int querySeqLen, unsigned int identityId);
    
private:
    int counter;
    
    float s_per_match;
    
    float s_per_pos;
        
    simd_int* __restrict thresholds_128;
    
    unsigned short  * __restrict thresholds;
    
    // float because it is needed for statistical calculations
    float * seqLens;
    
    float seqLenSum;
    
    int* steps;
    
    int nsteps;
            
    float matches_per_pos;

    
};
#endif /* defined(QUERYSCOREGLOBAL_H) */
