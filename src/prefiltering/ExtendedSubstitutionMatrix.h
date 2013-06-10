#ifndef EXTENDEDSUBSTITUIONMATRIXH 
#define EXTENDEDSUBSTITUIONMATRIXH 
#include <string>
#include <vector>
#include <algorithm>

class ExtendedSubstitutionMatrix 
{
public:
    ExtendedSubstitutionMatrix(short ** subMatrix, 
                               const size_t kmerSize,
                               const size_t alphabetSize);
    
    ~ExtendedSubstitutionMatrix();
    size_t size;
    // <match score, k-mer index>
    std::pair<short,unsigned int> ** scoreMatrix;
private: 
    std::vector<std::vector<int> > buildInput(size_t dimension,size_t range);
    void createCartesianProduct(
                 std::vector<std::vector<int> > & output,  // final result
                 std::vector<int>&  current_result,   // current result
                 std::vector<std::vector<int> >::const_iterator current_input, // current input
                 std::vector<std::vector<int> >::const_iterator end); // final input
    short calcScore(int * i_seq,int * j_seq,size_t seq_size,short **subMatrix);
    
};
#endif
