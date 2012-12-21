#ifndef EXTENDEDSUBSTITUIONMATRIXH 
#define EXTENDEDSUBSTITUIONMATRIXH 
#include <string>
#include <vector>
#include <algorithm>

class ExtendedSubstitutionMatrix 
{
public:
    ExtendedSubstitutionMatrix(short ** subMatrix, 
                               const size_t kmer_size,
                               const size_t alpherbet_size);
    
    ~ExtendedSubstitutionMatrix();
    size_t size;
    // <k-mer index, match score>
    std::vector<std::pair<short,size_t> > ** scoreMatrix;
private: 
    std::vector<std::vector<int> > build_input(size_t dimension,size_t range);
    void cart_product(
                 std::vector<std::vector<int> > & output,  // final result
                 std::vector<int>&  current_result,   // current result
                 std::vector<std::vector<int> >::const_iterator current_input, // current input
                 std::vector<std::vector<int> >::const_iterator end); // final input
    int calc_score(int * i_seq,int * j_seq,size_t seq_size,short **subMatrix);
    
};
#endif
