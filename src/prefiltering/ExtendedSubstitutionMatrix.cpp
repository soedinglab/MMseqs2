#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include <iostream>
#include <iterator>
#include <math.h>
#include <math.h>

#include <stdlib.h>

struct sort_by_score {
    bool operator()(const std::pair<short,unsigned int> left, const std::pair<short,unsigned int> right) {
        return left.first > right.first;
    }
};


ExtendedSubstitutionMatrix::ExtendedSubstitutionMatrix(short ** subMatrix, 
                                                       const size_t kmerSize,
                                                       const size_t alphabetSize){
    Indexer indexer( (int) alphabetSize, (int) kmerSize);
    this->size = pow(alphabetSize, kmerSize);
    // create permutation 
    std::vector<std::vector<int> > input(buildInput(kmerSize,alphabetSize));
    this->scoreMatrix = (std::pair<short,unsigned int> **) new std::pair<short,unsigned int> *[this->size];
    for(size_t i = 0; i < this->size;i++){
        this->scoreMatrix[i]=(std::pair<short,unsigned int> *) new std::pair<short,unsigned int> [this->size];
    }
    std::vector<std::vector<int> > permutation;
    std::vector<int> outputTemp;
    createCartesianProduct(permutation, outputTemp, input.begin(), input.end());
    
    // fill matrix  
    for(std::vector<int>::size_type i = 0; i != permutation.size(); i++) {
        unsigned int i_index=indexer.int2index(&permutation[i][0]);
        
        for(std::vector<int>::size_type j = 0; j != permutation.size(); j++) {
            unsigned int j_index=indexer.int2index(&permutation[j][0]);
            short score=calcScore(&permutation[i][0],&permutation[j][0],kmerSize,subMatrix);
            scoreMatrix[i_index][j].first=score;
            scoreMatrix[i_index][j].second=j_index;
        }
        std::sort (scoreMatrix[i_index], scoreMatrix[i_index]+this->size,sort_by_score());
    }
    
}


ExtendedSubstitutionMatrix::~ExtendedSubstitutionMatrix(){
    for(size_t i = 0; i < this->size; i++){
        delete[]  scoreMatrix[i];
    }
    delete[] scoreMatrix;
}

short ExtendedSubstitutionMatrix::calcScore(int * i_seq,int * j_seq,size_t seq_size,short **subMatrix){
    short score = 0;
    for(size_t i = 0; i < seq_size; i++){
        score+= subMatrix[i_seq[i]][j_seq[i]];
    }
    return score;
}

// Creates the input
std::vector<std::vector<int> > ExtendedSubstitutionMatrix::buildInput(size_t dimension,size_t range) {
    std::vector<std::vector<int> >  dimension_vector;
    
    for(size_t i = 0; i < dimension; i++) {
        std::vector<int> range_vector;
        for(size_t j = 0; j < range; j++) {
            range_vector.push_back(j);
        }
        dimension_vector.push_back(range_vector);
    }
    return dimension_vector;
}

// recursive algorithm to to produce cart. prod.
// At any given moment, "me" points to some Vi in the middle of the
// input data set.
//   for int i in *me:
//      add i to current result
//      recurse on next "me"
//
void ExtendedSubstitutionMatrix::createCartesianProduct(
                                              std::vector<std::vector<int> > & output,  // final result
                                              std::vector<int>&  current_result,   // current result
                                              std::vector<std::vector<int> >::const_iterator current_input, // current input
                                              std::vector<std::vector<int> >::const_iterator end) // final input
{
    if(current_input == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. Add the current result (rvi)
        // to the total set of results (rvvvi).
        output.push_back(current_result);
        return;
    }
    
    // need an easy name for my vector-of-ints
    const std::vector<int>& mevi = *current_input;
    for(std::vector<int>::const_iterator it = mevi.begin();it != mevi.end();it++) {
        current_result.push_back(*it);  // add ME
        createCartesianProduct(output, current_result, current_input+1, end);
        current_result.pop_back(); // clean current result off for next round
    }
}



