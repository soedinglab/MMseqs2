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

    size_t countElements(unsigned int *inputArray, const size_t N, CounterResult * output);

private:

    // this bit array should fit in L1/L2
    size_t duplicateBitArraySize;
    unsigned char * duplicateBitArray;
    // needed for lower bit hashing function
    unsigned int binCount;
    size_t binSize;
    unsigned int **bins;
    // array to keep the bin elements
    unsigned int * binDataFrame;

    // sub bin data structure
    unsigned int **subBins;
    unsigned int subBinCount;
    size_t subBinSize;
    unsigned int * subBinDataFrame;
    unsigned char * lookup;
    unsigned int *tmpElementBuffer;

    bool checkForOverflowAndResizeArray(unsigned int **bins,
                                        const unsigned int binCount,
                                        const size_t binSize);


    void reallocBinMemory(unsigned int const binCount, size_t const binSize);


    void setupBinPointer(unsigned int **bins, const unsigned int binCount,
                         unsigned int *binDataFrame, const size_t binSize);

    void hashElements(unsigned int const *inputArray, size_t N,
                      unsigned int **hashBins);

    size_t findDuplicates(unsigned int **bins, unsigned int binCount,
                          unsigned int **subBins, unsigned int subBinCount,
                          CounterResult * output);

    unsigned int highest_bit_set(size_t num);


    size_t countDuplicates(unsigned int **subBins, unsigned int subBinCount,
                           unsigned int *inputArray, const unsigned int elementCount,
                           CounterResult * output);

};
#endif