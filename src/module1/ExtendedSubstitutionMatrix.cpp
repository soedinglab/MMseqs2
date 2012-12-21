#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include <iostream>
#include <iterator>
#include <math.h>
#include <math.h>

#include <stdlib.h>




struct sort_by_score {
    bool operator()(const std::pair<short,size_t> &left, const std::pair<short,size_t> &right) {
        return left.first > right.first;
    }
};


ExtendedSubstitutionMatrix::ExtendedSubstitutionMatrix(short ** subMatrix, 
                                                       const size_t kmer_size,
                                                       const size_t alphabet_size){
    
    Indexer indexer(alphabet_size, kmer_size);
    
    this->size = pow(alphabet_size, kmer_size);
    
    // create permutation 
    std::vector<std::vector<int> > input(build_input(kmer_size,alphabet_size));    
    this->scoreMatrix = new std::vector<std::pair<short,size_t> >*[this->size];
    for(int i = 0; i < this->size;i++){
        this->scoreMatrix[i]=new std::vector<std::pair<short,size_t> >();  
    }
    std::vector<std::vector<int> > output;
    std::vector<int> outputTemp;
    cart_product(output, outputTemp, input.begin(), input.end());
    
    // fill matrix  
    for(std::vector<int>::size_type i = 0; i != output.size(); i++) {
        size_t i_index=indexer.int2index(&output[i][0]);
        std::vector<std::pair<short,size_t> > * i_vector= scoreMatrix[i_index];
        
        for(std::vector<int>::size_type j = 0; j != output.size(); j++) {
            size_t j_index=indexer.int2index(&output[j][0]);
            int score=calc_score(&output[i][0],&output[j][0],kmer_size,subMatrix);
            i_vector->push_back(std::make_pair(score,j_index));
        }
        std::sort (i_vector->begin(), i_vector->end(),sort_by_score()); 
        
    }
    
}


ExtendedSubstitutionMatrix::~ExtendedSubstitutionMatrix(){
    for(int i = 0; i < this->size; i++){
        delete  scoreMatrix[i];
    }
    delete[] scoreMatrix;
}

int ExtendedSubstitutionMatrix::calc_score(int * i_seq,int * j_seq,size_t seq_size,short **subMatrix){
    short score = 0;
    for(int i = 0; i < seq_size; i++){
        score+= subMatrix[i_seq[i]][j_seq[i]];
    }
    return score;
}

// Creates the input
std::vector<std::vector<int> > ExtendedSubstitutionMatrix::build_input(size_t dimension,size_t range) {
    std::vector<std::vector<int> >  vvi;
    
    for(int i = 0; i < dimension; i++) {
        std::vector<int> vi;
        for(int j = 0; j < range; j++) {
            vi.push_back(j);
        }
        vvi.push_back(vi);
    }
    return vvi;
}



// recursive algorithm to to produce cart. prod.
// At any given moment, "me" points to some Vi in the middle of the
// input data set.
//   for int i in *me:
//      add i to current result
//      recurse on next "me"
//
void ExtendedSubstitutionMatrix::cart_product(
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
        cart_product(output, current_result, current_input+1, end);
        current_result.pop_back(); // clean current result off for next round
    }
}



