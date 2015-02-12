#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef COUNTIN32ARRAY_H
#define COUNTIN32ARRAY_H
class CountInt32Array {
public:
    // 00000000000000000000000111111111
    static const unsigned int MASK_9 = 0x000001FF;
    static const unsigned int MASK_9_BIT = 9;

    CountInt32Array(unsigned int maxElement,
            size_t initBinSize);

    ~CountInt32Array();

    size_t countElements(unsigned int const *inputArray, size_t N,
            unsigned int *outputArray);

private:

    // this bit array should fit in L1/L2
    size_t duplicateBitArraySize;
    unsigned char * duplicateBitArray;
    // needed for lower bit hashing function
    unsigned int binCount;
    size_t binSize;
    unsigned int **bins;
    unsigned int * binDataFrame;

    bool checkForOverflowAndResizeArray(unsigned int **bins,
            const unsigned int binCount,
            const size_t binSize);


    void reallocBinMemory(unsigned int const binCount, size_t const binSize);


    void setupBinPointer(unsigned int **bins, const unsigned int binCount,
            unsigned int *binDataFrame, const size_t binSize);

    void hashElements(unsigned int const *inputArray, size_t N,
            unsigned int **hashBins, const unsigned int MASK);

    size_t findDuplicates(unsigned int **bins,
                               unsigned int binCount,
                               unsigned int *outputArray);
};
#endif