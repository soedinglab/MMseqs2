#include "CacheFriendlyOperations.h"
#include <new>
#include <iostream>
#include "IndexTable.h"
#include "Util.h"

template<unsigned int BINSIZE> CacheFriendlyOperations<BINSIZE>::CacheFriendlyOperations(size_t maxElement, size_t initBinSize) {
    // find nearest upper power of 2^(x)
    size_t size = pow(2, ceil(log(maxElement)/log(2)));
    size = std::max(size  >> MASK_0_5_BIT, (size_t) 1); // space needed in bit array
    duplicateBitArraySize = size;
    duplicateBitArray = new(std::nothrow) unsigned char[size];
    Util::checkAllocation(duplicateBitArray, "Can not allocate duplicateBitArray memory in CacheFriendlyOperations");
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));
    // find nearest upper power of 2^(x)
    initBinSize = pow(2, ceil(log(initBinSize)/log(2)));
    binSize = initBinSize;
    tmpElementBuffer = new(std::nothrow) TmpResult[binSize];
    Util::checkAllocation(tmpElementBuffer, "Can not allocate tmpElementBuffer memory in CacheFriendlyOperations");

    bins = new CounterResult*[BINCOUNT];
    binDataFrame = new(std::nothrow) CounterResult[BINCOUNT * binSize];
    Util::checkAllocation(binDataFrame, "Can not allocate binDataFrame memory in CacheFriendlyOperations");

}

template<unsigned int BINSIZE> CacheFriendlyOperations<BINSIZE>::~CacheFriendlyOperations<BINSIZE>(){
    delete [] duplicateBitArray;
    delete [] binDataFrame;
    delete [] tmpElementBuffer;
    delete [] bins;
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::countElements(IndexEntryLocal **input, CounterResult *output,
                                                                              size_t outputSize, unsigned short indexFrom, unsigned short indexTo,
                                                                              bool computeTotalScore)
{
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    CounterResult * lastPosition = (binDataFrame + BINCOUNT * binSize) - 1;

    for(unsigned int i = indexFrom; i < indexTo; i++){
        const size_t N = input[i + 1] - input[i];
        hashIndexEntry(i, input[i], N, this->bins, lastPosition);
    }
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return findDuplicates(this->bins, this->BINCOUNT, output, outputSize, computeTotalScore);
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::mergeElementsByScore(CounterResult *inputOutputArray, const size_t N) {
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    hashElements(inputOutputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return mergeDuplicates(this->bins, this->BINCOUNT, inputOutputArray);
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::mergeElementsByDiagonal(CounterResult *inputOutputArray, const size_t N) {
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    hashElements(inputOutputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return mergeDiagonalDuplicates(this->bins, this->BINCOUNT, inputOutputArray);
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::keepMaxScoreElementOnly(CounterResult *inputOutputArray, const size_t N) {
    newStart:
    setupBinPointer(bins, BINCOUNT, binDataFrame, binSize);
    hashElements(inputOutputArray, N, this->bins);
    if(checkForOverflowAndResizeArray(bins, BINCOUNT, binSize) == true) // overflowed occurred
        goto newStart;
    return keepMaxElement(this->bins, this->BINCOUNT, inputOutputArray);
}


template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::mergeDiagonalDuplicates(CounterResult **bins, unsigned int binCount,
                                                                                        CounterResult * output) {
    size_t doubleElementCount = 0;
    const CounterResult *bin_ref_pointer = binDataFrame;

    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t currBinSize = (bins[bin] - binStartPos);
        size_t n = currBinSize - 1;
        // write diagonals + 1 in reverse order in the byte array
        while ( n != static_cast<size_t>(-1) )
        {
            const unsigned int element = binStartPos[n].id >> (MASK_0_5_BIT);
            duplicateBitArray[element] = static_cast<unsigned char>(tmpElementBuffer[n].diagonal) + 1;
            --n;
        }
        // combine diagonals
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            output[doubleElementCount].id    = element.id;
            output[doubleElementCount].count = element.count;
            output[doubleElementCount].diagonal = element.diagonal;
//            std::cout << output[doubleElementCount].id << " " << (int)output[doubleElementCount].count << " " << (int)static_cast<unsigned char>(output[doubleElementCount].diagonal) << std::endl;
            // memory overflow can not happen since input array = output array
            doubleElementCount += (duplicateBitArray[hashBinElement] != static_cast<unsigned char>(tmpElementBuffer[n].diagonal)) ? 1 : 0;

            duplicateBitArray[hashBinElement] = static_cast<unsigned char>(element.diagonal);
        }
    }
    return doubleElementCount;
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::mergeDuplicates(CounterResult **bins, unsigned int binCount,
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
            output[doubleElementCount].diagonal = element.diagonal;
            // memory overflow can not happen since input array = output array
            doubleElementCount += (UNLIKELY(duplicateBitArray[hashBinElement] != 0  ) ) ? 1 : 0;
            duplicateBitArray[hashBinElement] = static_cast<unsigned char>(tmpElementBuffer[n].diagonal);
        }
    }
    return doubleElementCount;
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::findDuplicates(CounterResult **bins,
                                                                               unsigned int binCount,
                                                                               CounterResult * output,
                                                                               size_t outputSize,
                                                                               bool computeTotalScore) {
    size_t doubleElementCount = 0;
    const CounterResult * bin_ref_pointer = binDataFrame;
    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t currBinSize = (bins[bin] - binStartPos);
        size_t elementCount = 0;
        // find duplicates
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            //const unsigned int byteArrayPos = hashBinElement >> 3; // equal to  hashBinElement / 8
            //const unsigned char bitPosMask = 1 << (hashBinElement & 7);  // 7 = 00000111
            // check if duplicate element was found before
            const unsigned char currDiagonal = element.diagonal;
            //currDiagonal = (currDiagonal == 0) ? 200 : currDiagonal;
            const unsigned char prevDiagonal = duplicateBitArray[hashBinElement];
            tmpElementBuffer[elementCount].id = element.id;
            tmpElementBuffer[elementCount].diagonal = element.diagonal;
            elementCount += (UNLIKELY(currDiagonal == prevDiagonal)) ? 1 : 0;
            // set element corresponding bit in byte
            duplicateBitArray[hashBinElement] = currDiagonal;
        }
        // check for overflow
        if(doubleElementCount + elementCount >= outputSize){
            return doubleElementCount;
        }
//        // set memory to zero
        if(computeTotalScore){
            for (size_t n = 0; n < elementCount; n++) {
                const unsigned int element = tmpElementBuffer[n].id >> (MASK_0_5_BIT);
                duplicateBitArray[element] = 0;
            }
            // sum up score
            for (size_t n = 0; n < elementCount; n++) {
                const unsigned int element = tmpElementBuffer[n].id >> (MASK_0_5_BIT);
                duplicateBitArray[element] += (duplicateBitArray[element] < 255) ? 1 : 0;
            }
            // extract results
            for (size_t n = 0; n < elementCount; n++) {
                const unsigned int element = tmpElementBuffer[n].id;
                const unsigned int hashBinElement = element >> (MASK_0_5_BIT);
                output[doubleElementCount].id    = element;
                output[doubleElementCount].count = duplicateBitArray[hashBinElement];
                output[doubleElementCount].diagonal = tmpElementBuffer[n].diagonal;

                // memory overflow can not happen since input array = output array
                doubleElementCount += (duplicateBitArray[hashBinElement] != 0) ? 1 : 0;
                duplicateBitArray[hashBinElement] = 0;
            }
        }else{
            // set duplicate bit array to first diagonal + 1
            // so (duplicateBitArray[hashBinElement] != tmpElementBuffer[n].diagonal) is true
            size_t n = elementCount - 1;
            while ( n != static_cast<size_t>(-1) )
            {
                const unsigned int element = tmpElementBuffer[n].id >> (MASK_0_5_BIT);
                duplicateBitArray[element] = static_cast<unsigned char>(tmpElementBuffer[n].diagonal) + 1;
                --n;
            }

            // extract results
            for (size_t n = 0; n < elementCount; n++) {
                const unsigned int element = tmpElementBuffer[n].id;
                const unsigned int hashBinElement = element >> (MASK_0_5_BIT);
                output[doubleElementCount].id    = element;
                output[doubleElementCount].count = tmpElementBuffer[n].score;
                output[doubleElementCount].diagonal = tmpElementBuffer[n].diagonal;
    //            const unsigned char diagonal = static_cast<unsigned char>(tmpElementBuffer[n].diagonal);
                // memory overflow can not happen since input array = output array
    //            if(duplicateBitArray[hashBinElement] != tmpElementBuffer[n].diagonal){
    //                std::cout << "seq="<< output[doubleElementCount].id << "\tDiag=" << (int) output[doubleElementCount].diagonal
    //                <<  " dup.Array=" << (int)duplicateBitArray[hashBinElement] << " tmp.Arr="<< (int)tmpElementBuffer[n].diagonal << std::endl;
    //            }
                doubleElementCount += (duplicateBitArray[hashBinElement] != static_cast<unsigned char>(tmpElementBuffer[n].diagonal)) ? 1 : 0;

                duplicateBitArray[hashBinElement] = static_cast<unsigned char>(tmpElementBuffer[n].diagonal);
            }
        }
        // clean memory faster if current bin size is smaller duplicateBitArraySize
        if(currBinSize < duplicateBitArraySize/16){
            for (size_t n = 0; n < currBinSize; n++) {
                const unsigned int byteArrayPos = binStartPos[n].id >> (MASK_0_5_BIT);
                duplicateBitArray[byteArrayPos] = 0;
            }
        }else{
            memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));
        }
    }
    return doubleElementCount;
}

template<unsigned int BINSIZE> bool CacheFriendlyOperations<BINSIZE>::checkForOverflowAndResizeArray(CounterResult **bins,
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

template<unsigned int BINSIZE> void CacheFriendlyOperations<BINSIZE>::reallocBinMemory(const unsigned int binCount, const size_t binSize) {
    delete [] binDataFrame;
    delete [] tmpElementBuffer;
    binDataFrame     = new(std::nothrow) CounterResult[binCount * binSize];
    memset(binDataFrame, 0, sizeof(CounterResult) * binSize * binCount);
    Util::checkAllocation(binDataFrame, "Can not allocate reallocBinMemory memory in CacheFriendlyOperations::reallocBinMemory");
    tmpElementBuffer = new(std::nothrow) TmpResult[binSize];
    memset(tmpElementBuffer, 0, sizeof(TmpResult) * binSize);
    Util::checkAllocation(tmpElementBuffer, "Can not allocate tmpElementBuffer memory in CacheFriendlyOperations::reallocBinMemory");
}

template<unsigned int BINSIZE> void CacheFriendlyOperations<BINSIZE>::setupBinPointer(CounterResult **bins, const unsigned int binCount,
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

template<unsigned int BINSIZE> void CacheFriendlyOperations<BINSIZE>::hashElements(CounterResult *inputArray, size_t N, CounterResult **hashBins)
{
    CounterResult * lastPosition = (binDataFrame + BINCOUNT * binSize) - 1;
    for(size_t n = 0; n < N; n++) {
        const CounterResult element = inputArray[n];
        const unsigned int bin_id  = (element.id & MASK_0_5);
        hashBins[bin_id]->id       = element.id;
        hashBins[bin_id]->diagonal = element.diagonal;
        hashBins[bin_id]->count    = element.count;

        // do not write over boundary of the data frame
        hashBins[bin_id] += (hashBins[bin_id] >= lastPosition) ? 0 : 1;
    }
}

template<unsigned int BINSIZE> void CacheFriendlyOperations<BINSIZE>::hashIndexEntry(unsigned short position_i, IndexEntryLocal *inputArray,
                                                                             size_t N, CounterResult **hashBins, CounterResult * lastPosition)
{
    for(size_t n = 0; n < N; n++) {
        const IndexEntryLocal element = inputArray[n];
        const unsigned int bin_id = (element.seqId & MASK_0_5);
        hashBins[bin_id]->id    = element.seqId;
        hashBins[bin_id]->diagonal = position_i - element.position_j;
        // do not write over boundary of the data frame
//        std::cout << hashBins[bin_id]->id  << " " << position_i << " "
//        << element.position_j << " " << hashBins[bin_id]->diagonal << " " << position_i - element.position_j   << std::endl;
        hashBins[bin_id] += (hashBins[bin_id] >= lastPosition) ? 0 : 1;
    }
}

template<unsigned int BINSIZE> size_t CacheFriendlyOperations<BINSIZE>::keepMaxElement(CounterResult **bins,
                                                                               unsigned int binCount,
                                                                               CounterResult * output) {
    size_t doubleElementCount = 0;
    const CounterResult *bin_ref_pointer = binDataFrame;
    memset(duplicateBitArray, 0, duplicateBitArraySize * sizeof(unsigned char));
    for (size_t bin = 0; bin < binCount; bin++) {
        const CounterResult *binStartPos = (bin_ref_pointer + bin * binSize);
        const size_t currBinSize = (bins[bin] - binStartPos);
        // found max element and store it in duplicateBitArray
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            const unsigned char currScore = element.count;
            const unsigned char dbScore = duplicateBitArray[hashBinElement];
            const unsigned char maxScore = (currScore > dbScore) ? currScore : dbScore;
            duplicateBitArray[hashBinElement] = maxScore;
        }
        // extract final scores and set dubplicateBitArray to 0
        for (size_t n = 0; n < currBinSize; n++) {
            const CounterResult element = binStartPos[n];
            const unsigned int hashBinElement = element.id >> (MASK_0_5_BIT);
            output[doubleElementCount].id = element.id;
            output[doubleElementCount].count = element.count;
            output[doubleElementCount].diagonal = element.diagonal;
            // memory overflow can not happen since input array = output array
            bool found = (UNLIKELY(duplicateBitArray[hashBinElement] == element.count)) ? 1 : 0;
            doubleElementCount += found;
            duplicateBitArray[hashBinElement] = duplicateBitArray[hashBinElement] * (1 - found);
        }
    }
    return doubleElementCount;
}

template class CacheFriendlyOperations<2048>;
template class CacheFriendlyOperations<1024>;
template class CacheFriendlyOperations<512>;
template class CacheFriendlyOperations<256>;
template class CacheFriendlyOperations<128>;
template class CacheFriendlyOperations<64>;
template class CacheFriendlyOperations<32>;
template class CacheFriendlyOperations<16>;
template class CacheFriendlyOperations<8>;
template class CacheFriendlyOperations<4>;
template class CacheFriendlyOperations<2>;