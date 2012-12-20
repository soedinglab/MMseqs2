#ifndef KMERGENERATOR_H 
#define KMERGENERATOR_H 
#include <string>
#include <vector>
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"


typedef struct {
    size_t count;
    std::vector<std::pair<short, size_t> >* score_kmer_list;
} kmer_list;


class KmerGenerator 
{
    public: 
        KmerGenerator(size_t kmer_size,size_t alpherbet_size, short threshold, 
                  ExtendedSubstitutionMatrix * three,ExtendedSubstitutionMatrix * two );
        ~KmerGenerator();
        /*calculates the kmer list */
        kmer_list generateKmerList(const int * int_seq);
      

    private:
    
        /*creates the product between two vectors and write it to the output vector */
        int calcProduct(const std::vector<std::pair<short, size_t> > * vec1,
                           const std::vector<std::pair<short, size_t> > * vec2,
                           std::vector<std::pair<short, size_t> > * outputvec, 
                           const short cutoff1,const short possible_rest,const size_t pow);
    
        
        /* maximum return values */
        const static size_t VEC_LIMIT = 8000;
        /* min score  */
        short threshold;
        /* size of kmer  */
        size_t kmer_size;
        /* partition steps of the kmer size in (2,3)  */
        size_t divide_steps_count;
        /* divider of the steps (2,3) */
        size_t * divide_steps; 
        size_t * kmer_index;
        size_t * pow_per_step;
        short * max_score_per_vec;
        short * possible_rest;
        Indexer * indexer;
        ExtendedSubstitutionMatrix ** matrixLookup; 
        ExtendedSubstitutionMatrix * three; 
        ExtendedSubstitutionMatrix * two; 
        std::vector<std::pair<short, size_t> > ** outputvec;
    
        /* kmer splitting stragety (3,2)
           fill up the divide step and calls init_result_list */
        void calc_devide_stragety();
        /* init the output vectors for the kmer calculation*/
        void init_result_lists(size_t divide_steps);
    
};
#endif

