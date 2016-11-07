#ifndef COUNTIN32ARRAY_H
#define COUNTIN32ARRAY_H

#include <cstdlib>
#include <cmath>
#include <cstring>
#include "IndexTable.h"

#define IS_REPRESENTIBLE_IN_D_BITS(D, N)                \
  (((unsigned long) N >= (1UL << (D - 1)) && (unsigned long) N < (1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (31                                     \
                 + IS_REPRESENTIBLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(16, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(17, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(18, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(19, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(20, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(21, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(22, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(23, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(24, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(25, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(26, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(27, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(28, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(29, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(30, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(31, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(32, N)    \
                 )                                      \
   )

struct  __attribute__((__packed__))  CounterResult {
    unsigned int  id;
    unsigned short diagonal;
    unsigned char count;
};

template<unsigned int BINSIZE> class CacheFriendlyOperations{
public:
    // 00000000000000000000000111111111
    static const unsigned int MASK_0_5 = BINSIZE - 1;
    static const unsigned int MASK_0_5_BIT = BITS_TO_REPRESENT(MASK_0_5);

    CacheFriendlyOperations(size_t maxElement,
                    size_t initBinSize);

    ~CacheFriendlyOperations();

    size_t countElements(IndexEntryLocal **input, CounterResult *output,
                         size_t outputSize, unsigned short indexFrom, unsigned short indexTo, bool computeTotalScore);
    // merge elements in CounterResult
    // assumption is that each element (diagonalMatcher.id) exists maximal two times
    size_t mergeElementsByScore(CounterResult *inputOutputArray, const size_t N);

    // merge elements in CounterResult by diagonal
    // it combines elements with same ids that occurs after each other
    size_t mergeElementsByDiagonal(CounterResult *inputOutputArray, const size_t N);
    size_t keepMaxScoreElementOnly(CounterResult *inputOutputArray, const size_t N);
private:
    // this bit array should fit in L1/L2
    size_t duplicateBitArraySize;
    unsigned char * duplicateBitArray;
    // needed for lower bit hashing function
    const static unsigned int BINCOUNT = MASK_0_5 + 1;
    size_t binSize;
    // pointer for hashing
    CounterResult ** bins;
    // array to keep the bin elements
    CounterResult * binDataFrame;


    struct  __attribute__((__packed__))  TmpResult {
        unsigned int  id;
        unsigned short diagonal;
        unsigned char score;
    };
    // needed to temporary keep ids
    TmpResult *tmpElementBuffer;
    // detect if overflow occurs
    bool checkForOverflowAndResizeArray(CounterResult **bins,
                                        const unsigned int binCount,
                                        const size_t binSize);

    // extend the bin size
    void reallocBinMemory(unsigned int const binCount, size_t const binSize);

    // reset pointer to the bin start pos
    void setupBinPointer(CounterResult **bins, const unsigned int binCount,
                         CounterResult *binDataFrame, const size_t binSize);

    // hash input array based on MASK_0_5
    void hashElements(CounterResult *inputArray, size_t N, CounterResult **hashBins);

    // hash index entry and compute diagonal
    void hashIndexEntry(unsigned short position_i, IndexEntryLocal *inputArray,
                        size_t N, CounterResult **hashBins, CounterResult * lastPosition);

    // detect duplicates in diagonal
    size_t findDuplicates(CounterResult **bins, unsigned int binCount,
                          CounterResult * output, size_t outputSize, bool findDuplicates);

    // merge by id and combine score
    size_t mergeDuplicates(CounterResult **bins, unsigned int binCount, CounterResult *output);

    //
    size_t mergeDiagonalDuplicates(CounterResult **bins, unsigned int binCount, CounterResult *output);

    size_t keepMaxElement(CounterResult **pResult, const unsigned int bincount, CounterResult *pCounterResult);
};

#undef BITS_TO_REPRESENT
#undef IS_REPRESENTIBLE_IN_D_BITS
#endif
