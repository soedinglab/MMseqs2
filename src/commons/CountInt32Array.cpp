#include "CountInt32Array.h"
#include <iostream>
#include "Util.h"

CountInt32Array::CountInt32Array(size_t maxElement, size_t initBinSize) {
    // find nearest upper power of 2^(x)
    size_t size = pow(2, ceil(log(maxElement)/log(2)));
    size = std::max(size  >> MASK_0_5_BIT, (size_t) 1); // space needed in bit array
    duplicateBitArraySize = size;
    duplicateBitArray = new unsigned short[size];
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned short));
    // find nearest upper power of 2^(x)
    initBinSize = pow(2, ceil(log(initBinSize)/log(2)));
    binSize = initBinSize;
    tmpElementBuffer = new CounterResult[binSize];

    bins = new CounterResult*[BINCOUNT];
    binDataFrame = new CounterResult[BINCOUNT * binSize];
}

CountInt32Array::~CountInt32Array(){
    delete [] duplicateBitArray;
    delete [] binDataFrame;
    delete [] tmpElementBuffer;
    delete [] bins;
}

size_t CountInt32Array::countElements(CounterResult *inputOutputArray, const size_t N) {
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    hashElements(inputOutputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return findDuplicates(this->bins, this->BINCOUNT, inputOutputArray);
}

size_t CountInt32Array::mergeElements(CounterResult *inputOutputArray, const size_t N) {
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    hashElements(inputOutputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return mergeDuplicates(this->bins, this->BINCOUNT, inputOutputArray);
}

size_t CountInt32Array::mergeDuplicates(CounterResult **bins, unsigned int binCount,
                                        CounterResult * output) {
    size_t doubleElementCount = 0;
    const CounterResult *bin_ref_pointer = binDataFrame;
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));

    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t currBinSize = (bins[bin] - binStartPos);
        // merge double hits
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            const unsigned char currScore = element.count;
            const unsigned char dbScore = duplicateBitArray[hashBinElement];
            const unsigned char newScore = (currScore > 0xFF - dbScore) ? 0xFF : dbScore + currScore;
            duplicateBitArray[hashBinElement] = newScore;
        }
        // extract final scores and set dubplicateBitArray to 0
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            output[doubleElementCount].id    = element.id;
            output[doubleElementCount].count = duplicateBitArray[hashBinElement];
            // memory overflow can not happen since input array = output array
            doubleElementCount += (UNLIKELY(duplicateBitArray[hashBinElement] != 0  ) ) ? 1 : 0;;
            duplicateBitArray[hashBinElement] = 0;
        }
    }
    return doubleElementCount;
}

size_t CountInt32Array::findDuplicates(CounterResult **bins, unsigned int binCount,
                                       CounterResult * output) {
    size_t doubleElementCount = 0;
    const CounterResult * bin_ref_pointer = binDataFrame;
    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t currBinSize = (bins[bin] - binStartPos);
        unsigned char * duplicateBitArrayChar = (unsigned char *) duplicateBitArray;
        size_t elementCount = 0;
        // find duplicates
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            //const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
            //const unsigned char bitPosMask = 1 << (hashBinElement & 7);  // 7 = 00000111
            // check if duplicate element was found before
            const unsigned char currDiagonal = element.count;
            //currDiagonal = (currDiagonal == 0) ? 200 : currDiagonal;
            const unsigned int idx = hashBinElement*2;
            const unsigned char dbDiagonal = duplicateBitArrayChar[idx];
            const unsigned char prevScore = duplicateBitArrayChar[idx + 1];
            tmpElementBuffer[elementCount].id = element.id;
            tmpElementBuffer[elementCount].score = element.score + prevScore;
            const unsigned int step = (UNLIKELY(currDiagonal == dbDiagonal  ) ) ? 1 : 0;
            elementCount += step;
            // if same diagonal score = 0 else score = score;
            duplicateBitArrayChar[idx + 1] =  (1-step) *  element.score;

            // set element corresponding bit in byte
            duplicateBitArrayChar[idx] = currDiagonal;
        }

        // set memory to zero
        for (size_t n = 0; n < elementCount; n++) {
            const unsigned int element = tmpElementBuffer[n].id >> (MASK_0_5_BIT);
            duplicateBitArray[element] = 0;
        }

        // sum up score
        for (size_t n = 0; n < elementCount; n++) {
            const unsigned int element = tmpElementBuffer[n].id >> (MASK_0_5_BIT);
            duplicateBitArray[element]  += (duplicateBitArray[element] > 0xFF - tmpElementBuffer[n].score) ? 0xFF : tmpElementBuffer[n].score;
        }

        // extract results
        for (size_t n = 0; n < elementCount; n++) {
            const unsigned int element = tmpElementBuffer[n].id;
            const unsigned int hashBinElement = element >> (MASK_0_5_BIT);
            output[doubleElementCount].id    = element;
            output[doubleElementCount].count = duplicateBitArray[hashBinElement];
            // memory overflow can not happen since input array = output array
            doubleElementCount += (duplicateBitArray[hashBinElement] != 0) ? 1 : 0;
            duplicateBitArray[hashBinElement] = 0;
        }
        // clean memory faster if current bin size is smaller duplicateBitArraySize
        if(currBinSize < duplicateBitArraySize/16){
            for (size_t n = 0; n < currBinSize; n++) {
                const unsigned int byteArrayPos = binStartPos[n].id >> (MASK_0_5_BIT);
                duplicateBitArray[byteArrayPos] = 0;
            }
        }else{
            memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned short));
        }
    }
    return doubleElementCount;
}

bool CountInt32Array::checkForOverflowAndResizeArray(CounterResult **bins,
                                                     const unsigned int binCount,
                                                     const size_t binSize) {
    const CounterResult * bin_ref_pointer = binDataFrame;
    CounterResult * lastPosition = (binDataFrame + binCount * binSize) - 1;
    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t n = (bins[bin] - binStartPos);
        // if one bin has more elements than BIN_SIZE
        // or the current bin pointer is at the end of the binDataFrame
        // reallocate new memory
        if( n > binSize || bins[bin] >= lastPosition) {
            // overflow detected
            // find nearest upper power of 2^(x)
//            std::cout << "Found overlow " << n << std::endl;
            this->binSize = pow(2, ceil(log(binSize + 1)/log(2)));
            reallocBinMemory(binCount, this->binSize);
            return true;
        }
    }
    return false;
}

void CountInt32Array::reallocBinMemory(const unsigned int binCount, const size_t binSize) {
    delete [] binDataFrame;
    delete [] tmpElementBuffer;
    binDataFrame     = new CounterResult[binCount * binSize];
    tmpElementBuffer = new CounterResult[binSize];
}

void CountInt32Array::setupBinPointer(CounterResult **bins, const unsigned int binCount,
                                      CounterResult *binDataFrame, const size_t binSize)
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

void CountInt32Array::hashElements(CounterResult *inputArray, size_t N, CounterResult **hashBins)
{
    CounterResult * lastPosition = (binDataFrame + BINCOUNT * binSize) - 1;
    for(size_t n = 0; n < N; n++) {
        const CounterResult element = inputArray[n];
        const unsigned int bin_id = (element.id & MASK_0_5);
        hashBins[bin_id]->id    = element.id;
        hashBins[bin_id]->score = element.score;
        hashBins[bin_id]->count = element.count;
        // do not write over boundary of the data frame
        hashBins[bin_id] += (hashBins[bin_id] >= lastPosition) ? 0 : 1;
    }
}