#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H
#include <cstddef>

// <match score, k-mer index>
struct ScoreMatrix {
    const short *        score;
    const unsigned int * index;
    size_t elementSize;
    size_t rowSize;
    ScoreMatrix(short * scoreMatrix,
                unsigned int *indexMatrix,
                size_t elementSize,
                size_t rowSize):
                              score(scoreMatrix),
                              index(indexMatrix),
                              elementSize(elementSize),
                              rowSize(rowSize){}
    
    
};

#endif