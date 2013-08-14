#ifndef KMERGENERATOR_H 
#define KMERGENERATOR_H 
#include <string>
#include <vector>
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"


typedef struct {
    size_t count;
    // score, k-mer index
    std::pair<short, unsigned int> * scoreKmerList;
} KmerGeneratorResult;



class KmerGenerator 
{
    public: 
        KmerGenerator(size_t kmerSize,size_t alphabetSize, short threshold,
                  ExtendedSubstitutionMatrix * three,ExtendedSubstitutionMatrix * two );
        ~KmerGenerator();
        /*calculates the kmer list */
        KmerGeneratorResult generateKmerList(const int * intSeq);


    private:
    
        /*creates the product between two arrays and write it to the output array */
        int calculateArrayProduct( const std::pair<short,unsigned int> * array1,
                              const size_t array1Size,
                              const std::pair<short,unsigned int> * array2,
                              const size_t array2Size,
                              std::pair<short,unsigned int> * outputvec,
                              const short cutoff1,const short possibleRest,
                              const unsigned int pow);
    
    
        /* maximum return values */
        const static size_t VEC_LIMIT = 20000;
        /* min score  */
        short threshold;
        /* size of kmer  */
        size_t kmerSize;
        /* partition steps of the kmer size in (2,3)  */
        size_t divideStepCount;
        /* divider of the steps (2,3) */
        unsigned int * divideStep;
        unsigned int * kmerIndex;
        unsigned int * stepMultiplicator;
        short * highestScorePerArray;
        short * possibleRest;
        Indexer * indexer;
        ExtendedSubstitutionMatrix ** matrixLookup; 
        ExtendedSubstitutionMatrix * three; 
        ExtendedSubstitutionMatrix * two; 
        std::pair<short, unsigned int>  ** outputArray;
    
        /* kmer splitting stragety (3,2)
           fill up the divide step and calls init_result_list */
        void calcDivideStrategy();
        /* init the output vectors for the kmer calculation*/
        void initResultList(size_t divideSteps);
    
};
#endif

