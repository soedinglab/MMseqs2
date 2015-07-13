#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/cdefs.h>
#include <zconf.h>

#ifndef COUNTIN32ARRAY_H
#define COUNTIN32ARRAY_H

struct CounterResult {
    unsigned int  id;
    unsigned char count;
};
class CountInt32Array {
public:
    // 00000000000000000000000111111111
    static const unsigned int MASK_0_5 = 0x0000003F;
    static const unsigned int MASK_0_5_BIT = 6;
    static const unsigned int MASK_6_11 = 0x00000FC0;
    static const unsigned int MASK_6_11_BIT = 6;

    CountInt32Array(size_t maxElement,
                    size_t initBinSize);

    ~CountInt32Array();

    size_t countElements(CounterResult *inputArray, const size_t N, CounterResult * output);

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

    unsigned int *tmpElementBuffer;

    bool checkForOverflowAndResizeArray(CounterResult **bins,
                                        const unsigned int binCount,
                                        const size_t binSize);


    void reallocBinMemory(unsigned int const binCount, size_t const binSize);


    void setupBinPointer(CounterResult **bins, const unsigned int binCount,
                         CounterResult *binDataFrame, const size_t binSize);

    void hashElements(CounterResult *inputArray, size_t N, CounterResult **hashBins);

    size_t findDuplicates(CounterResult **bins, unsigned int binCount,
                          CounterResult * output);

    unsigned int highest_bit_set(size_t num);

};
#endif