#ifndef KMERGENERATOR_H 
#define KMERGENERATOR_H 
#include <string>
#include <vector>
#include "Indexer.h"
#include "ScoreMatrix.h"
#include "Debug.h"


class KmerGenerator 
{
    public: 
        KmerGenerator(size_t kmerSize,size_t alphabetSize, short threshold);
        ~KmerGenerator();
        /*calculates the kmer list */
        std::pair<size_t *, size_t> generateKmerList(const unsigned char * intSeq, bool addIdentity = false);

        /* kmer splitting stragety (3,2)
         fill up the divide step and calls init_result_list */
        void setDivideStrategy(ScoreMatrix * three, ScoreMatrix * two );

        /* kmer splitting stragety (1)
         fill up the divide step and calls init_result_list */
        void setDivideStrategy(ScoreMatrix ** one);

	    void setThreshold(short threshold);
    private:
    
        /*creates the product between two arrays and write it to the output array */
        size_t calculateArrayProduct(const short        * __restrict scoreArray1,
                                  const size_t       * __restrict indexArray1,
                                  const size_t array1Size,
                                  const short        * __restrict scoreArray2,
                                  const unsigned int * __restrict indexArray2,
                                  const size_t array2Size,
                                  short              * __restrict outputScoreArray,
                                  size_t             * __restrict outputIndexArray,
                                  const short cutoff1,
                                  const short possibleRest,
                                  const size_t pow);
    
    
        /* maximum return values */
        /* 48   MB */
        const static size_t MAX_KMER_RESULT_SIZE = 262144*32;
        /* min score  */
        short threshold;
        /* size of kmer  */
        size_t kmerSize;
        /* partition steps of the kmer size in (2,3)  */
        size_t divideStepCount;
        /* divider of the steps (2,3) */
        unsigned int * divideStep;
        size_t * kmerIndex;
        size_t * stepMultiplicator;
        short * highestScorePerArray;
        short * possibleRest;
        Indexer * indexer;
        ScoreMatrix  ** matrixLookup;
        short        ** outputScoreArray;
        size_t       ** outputIndexArray;


        /* init the output vectors for the kmer calculation*/
        void initDataStructure();
    
};
#endif

