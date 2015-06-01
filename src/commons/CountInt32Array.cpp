#include "CountInt32Array.h"

#include <iostream>

CountInt32Array::CountInt32Array(unsigned int maxElement, size_t initBinSize) {
    // find nearest upper power of 2^(x)
    size_t size = pow(2, ceil(log(maxElement)/log(2)));
    size = size  >> (MASK_9_BIT + 3); // space needed in bit array
    size = std::max(size, (size_t) 1); // minimal size of 1;
    duplicateBitArraySize = size;
    duplicateBitArray = new unsigned char[size];
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));
    // init data for bins
    binCount = MASK_9 + 1;
    // find nearest upper power of 2^(x)
    initBinSize = pow(2, ceil(log(initBinSize)/log(2)));
    binSize = initBinSize;
    // pointer to start points of the bins
    bins = new unsigned int *[binCount];
    binDataFrame = new unsigned int[binCount * binSize];
}

CountInt32Array::~CountInt32Array(){
    delete [] duplicateBitArray;
    delete [] bins;
    delete [] binDataFrame;
}

size_t CountInt32Array::countElements(unsigned int *inputArray, const size_t N, unsigned int *outputArray,
                                      const unsigned int * lastOutputArrayPtr) {
    newStart:
    setupBinPointer(bins, binCount, binDataFrame, binSize);
    hashElements(inputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, binCount, binSize) == true) // overflowed occurred
        goto newStart;
    return findDuplicates(this->bins, binCount, outputArray, lastOutputArrayPtr);
}

size_t CountInt32Array::findDuplicates(unsigned int **bins, unsigned int binCount,
                                       unsigned int *outputArray,
                                       const unsigned int * lastOutputArrayPtr) {
    size_t pos = 0;
    const unsigned int * bin_ref_pointer = binDataFrame;
    for (size_t bin = 0; bin < binCount; bin++) {
        const unsigned int *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t N = (bins[bin] - binStartPos);
        size_t startPos = pos;
        // find duplicates
        for (size_t n = 0; n < N; n++) {
            const unsigned int element = binStartPos[n];
            const unsigned int hashBinElement = element >> (MASK_9_BIT);
            const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
            const unsigned char bitPosMask = 1 << (hashBinElement & 7);  // 7 = 00000111
            if (duplicateBitArray[byteArrayPos] & bitPosMask){
                if(outputArray + pos == lastOutputArrayPtr)
                    goto outer;
                outputArray[pos] = element;
                pos++;

            }
            duplicateBitArray[byteArrayPos] |= bitPosMask;
        }
        // append first position by iterating checking the results
//        for(size_t p = startPos; p < endPos; p++){
//            const unsigned int element = outputArray[p];
//
//            const unsigned int hashBinElement = element >> (MASK_9_BIT);
//            const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
//            const unsigned char bitPosMask = 1 << (hashBinElement & 7); // 7 = 00000111
//            if (duplicateBitArray[byteArrayPos] & bitPosMask){
//                outputArray[pos] = element;
//                pos++;
//            }
//            // unset position
//            duplicateBitArray[byteArrayPos] &= ~(bitPosMask);
//        }
        // clean memory the fast way if size is smaller duplicateBitArraySize
        outer:
        if(duplicateBitArraySize > N){
            for (size_t n = 0; n < N; n++) {
                const unsigned int element = binStartPos[n];
                const unsigned int byteArrayPos = element >> (MASK_9_BIT + 3) ;
                duplicateBitArray[byteArrayPos] = 0;
            }
        }else{
            memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));
        }

    }

    return pos;
}

bool CountInt32Array::checkForOverflowAndResizeArray(unsigned int **bins,
        const unsigned int binCount,
        const size_t binSize) {
    const unsigned int * bin_ref_pointer = binDataFrame;
    unsigned int * lastPosition = (binDataFrame + binCount * binSize) - 1;
    for (size_t bin = 0; bin < binCount; bin++) {
        const unsigned int *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t n = (bins[bin] - binStartPos);
        // if one bin has more elements than BIN_SIZE
        // or the current bin pointer is at the end of the binDataFrame
        // reallocate new memory
        if( n > binSize || (bins[bin] - lastPosition) == 0) {
            // overflow detected
            // find nearest upper power of 2^(x)
            std::cout << "Found overlow" << std::endl;
            this->binSize = pow(2, ceil(log(binSize + 1)/log(2)));
            reallocBinMemory(binCount, this->binSize);
            return true;
        }
    }
    return false;
}


void CountInt32Array::reallocBinMemory(const unsigned int binCount, const size_t binSize) {
    delete [] binDataFrame;
    binDataFrame = new unsigned int[binCount * binSize];
}

void CountInt32Array::setupBinPointer(unsigned int **bins, const unsigned int binCount,
        unsigned int *binDataFrame, const size_t binSize)
{
    // Example binCount = 3
    // bin start             |-----------------------|-----------------------| bin end
    //    segments[bin_step][0]
    //                            segments[bin_step][1]
    //                                                    segments[bin_step][2]
    size_t curr_pos = 0;
    for(size_t bin = 0; bin < binCount; bin++){
        bins[bin] = binDataFrame + curr_pos;
        curr_pos += binSize;
    }
}


void CountInt32Array::hashElements(unsigned int const *inputArray, size_t N,
        unsigned int **hashBins)
{
    unsigned int * lastPosition = (binDataFrame + binCount * binSize) - 1;
    for(size_t n = 0; n < N; n++) {
        const unsigned int element = inputArray[n];
        const unsigned int bin_id = (element & MASK_9);
        *hashBins[bin_id] = element;
        // do not write over boundary of the data frame
        hashBins[bin_id] += ((hashBins[bin_id] - lastPosition) != 0) ? 1 : 0;
    }
}