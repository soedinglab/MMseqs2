//
//  QueryScoreLocal.cpp
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//
#include "QueryScoreLocal.h"
#include "QueryScore.h"

QueryScoreLocal::QueryScoreLocal(size_t dbSize, unsigned int *seqLens, int k, short kmerThr,
                                 double kmerMatchProb, float zscoreThr, size_t binSize)
: QueryScore(dbSize, seqLens, k, kmerThr, kmerMatchProb, zscoreThr)    // Call the QueryScore constructor
{
    localResultSize = 0;
    localResults = new unsigned int[MAX_LOCAL_RESULT_SIZE];
    this->binSize = binSize;
    binData = new unsigned int[BIN_COUNT * binSize];
    counter = new CountInt32Array(dbSize, this->binSize / 4);
}

QueryScoreLocal::~QueryScoreLocal(){
    delete [] localResults;
    delete [] binData;
    delete counter;
}

void QueryScoreLocal::reset() {
//    memset (scores_128, 0, scores_128_size * 2);
    scoresSum = 0;
    numMatches = 0;
    localResultSize = 0;
}

void QueryScoreLocal::setPrefilteringThresholds() { }

std::pair<hit_t *, size_t> QueryScoreLocal::getResult (int querySeqLen, unsigned int identityId){
    size_t elementCounter = 0;
    std::sort(localResults, localResults + localResultSize);
    unsigned short scoreMax  = 1;  // current best score
    unsigned int seqIdPrev = localResults[0];
    for (unsigned int i = 1; i < localResultSize; i++) {
        unsigned int seqIdCurr = localResults[i];
        // if new sequence occurs or last element
        if (seqIdCurr != seqIdPrev || (i + 1 == localResultSize)) {
            // write result to list
            hit_t *result = (resList + elementCounter);
            result->seqId = seqIdPrev;
            result->zScore = scoreMax;
            result->prefScore = scoreMax;
            elementCounter++;
            if (elementCounter >= MAX_RES_LIST_LEN)
                break;
            seqIdPrev = seqIdCurr;
            // reset values
            scoreMax = 0;  // current best score
        }
        scoreMax += 1;   // add score of current kmer match
    }
    // sort hits by score
    //TODO maybe sort of hits not needed if SW is included
    std::sort(resList, resList + elementCounter, compareHits);
    std::pair<hit_t *, size_t>  pair = std::make_pair(this->resList, elementCounter);
    return pair;
}

void QueryScoreLocal::evaluateBins(){
    for(size_t bin = 0; bin < BIN_COUNT; bin++){
        const unsigned int *binStartPos = (binData + bin * binSize);
        const size_t N = (diagonalBins[bin] - binStartPos);
        localResultSize += counter->countElements(binStartPos, N, localResults + localResultSize);
    }
}

void QueryScoreLocal::setupBinPointer() {
    // Example binCount = 3
    // bin start             |-----------------------|-----------------------| bin end
    //    segments[bin_step][0]
    //                            segments[bin_step][1]
    //                                                    segments[bin_step][2]
    size_t curr_pos = 0;
    for(size_t bin = 0; bin < BIN_COUNT; bin++){
        diagonalBins[bin] = binData + curr_pos;
        curr_pos += binSize;
    }
}


void QueryScoreLocal::reallocBinMemory(const unsigned int binCount, const size_t binSize) {
    delete [] binData;
    binData = new unsigned int[binCount * binSize];
}

bool QueryScoreLocal::checkForOverflowAndResizeArray() {
    const unsigned int * bin_ref_pointer = binData;
    unsigned int * lastPosition = (binData + BIN_COUNT * binSize) - 1;
    for (size_t bin = 0; bin < BIN_COUNT; bin++) {
        const unsigned int *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t n = (diagonalBins[bin] - binStartPos);
        // if one bin has more elements than binSize
        // or the current bin pointer is at the end of the binDataFrame
        // reallocate new memory
        if( n > binSize || (diagonalBins[bin] - lastPosition) == 0) {
            // overflow detected
            // find nearest upper power of 2^(x)
            std::cout << "Diagonal Found overlow" << std::endl;
            this->binSize = pow(2, ceil(log(binSize + 1)/log(2)));
            reallocBinMemory(BIN_COUNT, this->binSize);
            return true;
        }
    }
    return false;
}
