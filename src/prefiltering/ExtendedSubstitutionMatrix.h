#ifndef EXTENDEDSUBSTITUIONMATRIXH
#define EXTENDEDSUBSTITUIONMATRIXH
#include <string>
#include <vector>
#include <algorithm>

#include "ScoreMatrix.h" // ScoreMatrix
#include "BaseMatrix.h"

class ExtendedSubstitutionMatrix
{
public:
    static ScoreMatrix calcScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize);
    static void freeScoreMatrix(ScoreMatrix& matrix);

    static short calcScore(int * i_seq,int * j_seq,size_t seq_size,short **subMatrix);

private:
    static std::vector<std::vector<int> > buildInput(size_t dimension,size_t range);
    static void createCartesianProduct(
                                std::vector<std::vector<int> > & output,  // final result
                                std::vector<int>&  current_result,   // current result
                                std::vector<std::vector<int> >::const_iterator current_input, // current input
                                std::vector<std::vector<int> >::const_iterator end); // final input

};
#endif
