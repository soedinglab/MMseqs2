#ifndef COUNTIN32ARRAY_H
#define COUNTIN32ARRAY_H

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/cdefs.h>
#include <IndexTable.h>

struct  __attribute__((__packed__))  CounterResult {
    unsigned int  id;
    unsigned short diagonal;
    unsigned char count;
};



class CountInt32Array {
public:
    // 00000000000000000000000111111111
    static const unsigned int MASK_0_5 = 0x0000007F;
    static const unsigned int MASK_0_5_BIT = 7;

    CountInt32Array(size_t maxElement,
                    size_t initBinSize);

    ~CountInt32Array();

    size_t countElements(IndexEntryLocal **input, CounterResult *output, unsigned int outputSize,
                         unsigned short indexFrom, unsigned short indexTo);
    // merge elements in CounterResult
    // assumption is that each element (counter.id) exists maximal two times
    size_t mergeElements(CounterResult *inputOutputArray, const size_t N);

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
                          CounterResult * output, unsigned int outputSize);


    size_t mergeDuplicates(CounterResult **bins, unsigned int binCount, CounterResult *output);
};

#endif
