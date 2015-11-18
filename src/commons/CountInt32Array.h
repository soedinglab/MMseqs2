#ifndef COUNTIN32ARRAY_H
#define COUNTIN32ARRAY_H

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/cdefs.h>
#include <zconf.h>

struct  __attribute__((__packed__))  CounterResult {
    unsigned int  id;
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

    size_t countElements(CounterResult *inputOutputArray, const size_t N);
    // merge elements in CounterResult
    // assumption is that each element (counter.id) exists maximal two times
    size_t mergeElements(CounterResult *inputOutputArray, const size_t N);


    // detect duplicates in diagonal
    size_t findDuplicates(const CounterResult * binData, const size_t binSize, CounterResult **bins, unsigned int binCount,
                          CounterResult * output);
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
    // needed to temporary keep ids
    unsigned int *tmpElementBuffer;
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

    size_t mergeDuplicates(CounterResult **bins, unsigned int binCount, CounterResult *output);
};

#endif
