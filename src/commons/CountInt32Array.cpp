#include "CountInt32Array.h"

#include <iostream>

CountInt32Array::CountInt32Array(unsigned int maxElement, size_t initBinSize) {
    // find nearest upper power of 2^(x)
    size_t size = pow(2, ceil(log(maxElement)/log(2)));
    size = size  >> (MASK_9_BIT + 3);; // space needed in bit array
    duplicateBitArraySize = size;
    duplicateBitArray = new unsigned char[size];
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));

    // init data for bins
    binCount = MASK_9 + 1;
    // find nearest upper power of 2^(x)
    initBinSize = pow(2, ceil(log(initBinSize)/log(2)));
    binSize = initBinSize;
    // +2 is needed so keep the 0 element as reference point and the last pointer for the overflow check
    bins = new unsigned int *[binCount + 2];
    binDataFrame = new unsigned int[(binCount + 2) * binSize];
}

CountInt32Array::~CountInt32Array(){
    delete [] duplicateBitArray;
    delete [] bins;
    delete [] binDataFrame;
}

size_t CountInt32Array::countElements(unsigned int const *inputArray, size_t N,
        unsigned int *outputArray) {
    newStart:
    setupBinPointer(bins, binCount, binDataFrame, binSize);
    hashElements(inputArray, N, this->bins, MASK_9); //TODO what happens in case of overflow
    if(checkForOverflowAndResizeArray(bins, binCount, binSize) == true) // overflowed occurred
        goto newStart;
    return findDuplicates(this->bins, binCount, outputArray);
}

size_t CountInt32Array::findDuplicates(unsigned int **bins,
        unsigned int binCount,
        unsigned int *outputArray) {
    size_t pos = 0;
    const unsigned int * bin_ref_pointer = bins[0];
    for (size_t bin = 0; bin < binCount; bin++) {
        const unsigned int *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t N = (bins[bin + 1] - binStartPos);
        size_t startPos = pos;
        // find duplicates
        for (size_t n = 0; n < N; n++) {
            const unsigned int element = binStartPos[n];
            const unsigned int hashBinElement = element >> (MASK_9_BIT);
            const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
            const unsigned char bitPosMask = 1 << (hashBinElement & 7);  // 7 = 00000111
            if (duplicateBitArray[byteArrayPos] & bitPosMask){
                outputArray[pos] = element;
                pos++;
            }
            duplicateBitArray[byteArrayPos] |= bitPosMask;
        }
        // append first position by iterating checking the results
        size_t endPos = pos;
        for(size_t p = startPos; p < endPos; p++){
            const unsigned int element = outputArray[p];
            
            const unsigned int hashBinElement = element >> (MASK_9_BIT);
            const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
            const unsigned char bitPosMask = 1 << (hashBinElement & 7); // 7 = 00000111
            if (duplicateBitArray[byteArrayPos] & bitPosMask){
                outputArray[pos] = element;
                pos++;
            }
            // unset position
            duplicateBitArray[byteArrayPos] &= ~(bitPosMask);
        }
        // clean memory the fast way if size is smaller duplicateBitArraySize
        if(N < duplicateBitArraySize){
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
    const unsigned int * bin_ref_pointer = bins[0];
    for (size_t bin = 0; bin < binCount; bin++) {
        const unsigned int *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t N = (bins[bin + 1] - binStartPos);
        if(N > binSize ){
            // overflow detected
            // find nearest upper power of 2^(x)
            std::cout << "Found overlow" << std::endl;
            this->binSize = pow(2, ceil(log(N)/log(2)));
            reallocBinMemory(binCount, this->binSize);
            return true;
        }
    }
    return false;
}


void CountInt32Array::reallocBinMemory(unsigned int const binCount, size_t const binSize) {
    delete [] binDataFrame;
    binDataFrame = new unsigned int[(binCount + 1) * binSize];
}

void CountInt32Array::setupBinPointer(unsigned int **bins, const unsigned int binCount,
        unsigned int *binDataFrame, const size_t binSize)
{
    // Example binCount = 3
    // bin start             |-----------------------|-----------------------| bin end
    //    segments[bin_step][0]
    //    segments[bin_step][1]
    //                            segments[bin_step][2]
    //                                                    segments[bin_step][3]
    //                                                    segments[bin_step][4] last bin for overflow check
    size_t curr_pos = 0;
    bins[0] = binDataFrame;
    for(size_t bin = 0; bin <= binCount; bin++){
        bins[bin+1] = binDataFrame + curr_pos;
        curr_pos += binSize;
    }
    bins[binCount + 1] = binDataFrame + curr_pos; // do not write over this position
}


void CountInt32Array::hashElements(unsigned int const *inputArray, size_t N,
        unsigned int **hashBins, const unsigned int MASK)
{
    for(size_t n = 0; n < N; n++) {
        const unsigned int element = inputArray[n];
        const unsigned int bin_id = (element & MASK) + 1; // dont write in 0
        if((hashBins[bin_id] - hashBins[bin_id+1]) == 0) // check that we dont write outside the memory
            break;
        *hashBins[bin_id] = element;
        hashBins[bin_id]++;
    }
}