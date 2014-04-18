#ifndef KMERGENERATOR_H 
#define KMERGENERATOR_H 
#include <string>
#include <vector>
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "../commons/Util.h"
#include "../commons/Debug.h"



typedef struct {
    size_t count;
    // score, k-mer index
    short        * score;
    unsigned int * index;
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
        int calculateArrayProduct(const short        * scoreArray1,
                                  const unsigned int * indexArray1,
                                  const size_t array1Size,
                                  const short        * scoreArray2,
                                  const unsigned int * indexArray2,
                                  const size_t array2Size,
                                  short              * outputScoreArray,
                                  unsigned int       * outputIndexArray,
                                  const short cutoff1,
                                  const short possibleRest,
                                  const unsigned int pow);
    
    
        /* maximum return values */
        const static size_t MAX_KMER_RESULT_SIZE = 20000;
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
        short        ** outputScoreArray;
        unsigned int ** outputIndexArray;

        /* kmer splitting stragety (3,2)
           fill up the divide step and calls init_result_list */
        void calcDivideStrategy();
        /* init the output vectors for the kmer calculation*/
        void initResultList(size_t divideSteps);
    
};
#endif

