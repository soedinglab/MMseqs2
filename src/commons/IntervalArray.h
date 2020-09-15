//
// Created by Martin Steinegger on 8/4/20.
//

#ifndef MMSEQS_INTERVALARRAY_H
#define MMSEQS_INTERVALARRAY_H

#include "MathUtil.h"
#include <utility>
#include <list>
#include <iostream>
#include <string.h>

class IntervalArray {

public:
    struct Range{
        unsigned int index;
        unsigned int start;
        unsigned int end;
        Range(){};
        Range(unsigned int index, unsigned int start, unsigned int end)
                : index(index), start(start), end(end)
        {}

        static bool compareStart(const Range &x, const Range &y){
            return (x.start < y.start);
        }

        static bool compareEndStart(const Range &x, const Range &y){
            return (x.start < y.end);
        }

    };
    IntervalArray(){
        const int size = 64;
        array=(unsigned char *)calloc(size, sizeof(unsigned char));
        arraySizeInBytes = size;
        maxSizeInByte = size*8;
    }


    ~IntervalArray(){
        free(array);
        array = NULL;
    }

    void reset(){
        int ceilMax = MathUtil::ceilIntDivision(std::max(1,maxSizeInByte), 8);
        memset(array, 0, ceilMax * sizeof(unsigned char));
        ranges.clear();
        maxSizeInByte = 0;
    }

    unsigned char getLowMask(unsigned int rest) {
        const unsigned char mask[8] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0,0x80};
        return mask[rest];
    }


    unsigned char getHighRest(unsigned int rest) {
        const unsigned char mask[8] = { 0xFF, 0x7F, 0x3F, 0x1F, 0x0F,0x07, 0x03, 0x01};
        return mask[rest];
    }


    void insert(int low, int high){
        if(low > high){
            std::swap(low, high);
        }
        //insert(array, low, high);
        maxSizeInByte = std::max(high+1, maxSizeInByte + 1);
        int ceilMax = MathUtil::ceilIntDivision(maxSizeInByte, 8);
        if(ceilMax >= arraySizeInBytes){
            int prevSize = arraySizeInBytes;
            arraySizeInBytes = std::max(arraySizeInBytes * 2, ceilMax+1);
            array = (unsigned char *) realloc(array,  arraySizeInBytes);
            memset(array+prevSize, 0,  arraySizeInBytes-prevSize);
        }
        bool lowFound = isSet(low);
        bool highFound =  isSet(high);;
        if((lowFound == true && highFound == true)){
            return;
        }
        unsigned int startPos=low/8;
        unsigned int startRest=low%8;
        unsigned int endPos=high/8;
        unsigned int endRest=high%8;
        for(size_t pos = startPos+1; pos < endPos; pos++ ){
            array[pos] = 0xFF;
        }
        unsigned char lowMask = getLowMask(startRest);
        unsigned char highMask = getHighRest(7-endRest);
        if(startPos == endPos){
            array[startPos] |= lowMask & highMask;
        }else{
            array[startPos] |= lowMask;
            array[endPos] |= highMask;
        }
    }

    bool checkOverlap(int i1Low, int i1High, int i2Low, int i2High)
    {
        return (i1Low <= i2High && i2Low<= i1High) ? true : false;
    }

    size_t getRangesSize(){
        return ranges.size();
    }

    int findIndex(int low, int high){

        if(ranges.size() == 0){
            return -1;
        }
        // check for overlaps between the intervals
        if(low > high){
            std::swap(low, high);
        }

        Range val;
        val.start = low;
        std::vector<Range>::iterator it;
        it = std::lower_bound(ranges.begin(), ranges.end(), val, Range::compareStart);
        if(it == ranges.end()){
            it--;
        }
        bool overlap = checkOverlap(it->start, it->end, low, high);
        if(overlap == false && ranges.begin()!=it){
            --it;
            overlap = checkOverlap(it->start, it->end, low, high);
        }

        if(overlap == false){
            return -1;
        }else{
            return it->index;
        }
    }

    bool doesOverlap(int low, int high)
    {
        // check for overlaps between the intervals
        if(low > high){
            std::swap(low, high);
        }
        bool lowFound;
        if(low >= maxSizeInByte){
            return false;
        }else{
            lowFound = isSet(low);
        }
        bool highFound;
        if(high >= maxSizeInByte){
            highFound = false;
        }else{
            highFound = isSet(high);;
        }
        if(lowFound || highFound){
            return true;
        }
        // check if interval is contained in low, high
        high = std::min(high, maxSizeInByte);
        unsigned int startPos=low/8;
        unsigned int startRest=low%8;
        unsigned int endPos=high/8;
        unsigned int endRest=high%8;
        int foundOverlap = 0;
        for(size_t pos = startPos+1; pos < endPos && foundOverlap == 0; pos++ ){
            foundOverlap += (array[pos]>0);
        }
        unsigned char lowMask = getLowMask(startRest);
        unsigned char highMask = getHighRest(7-endRest);
        if(startPos == endPos){
            foundOverlap += (array[startPos] & lowMask & highMask);
        }else{
            foundOverlap += array[startPos] & lowMask;
            foundOverlap += array[endPos] & highMask;
        }
        return (foundOverlap==0)? false : true;
    }

    bool isSet(int pos){
        unsigned int posIdx = pos/8;
        unsigned int posRest= pos%8;
        unsigned char value = array[posIdx];
        unsigned char check = (1U << posRest);
        return check & value;
    }


    void print(){
        bool started = false;
        for(int pos = 0; pos < maxSizeInByte; pos++){
            if(isSet(pos)  && started == false){
                started = true;
                std::cout << "[" << pos << ", ";
            }
            if(isSet(pos)==false && started== true){
                started = false;
                std::cout << pos - 1 << "]" << std::endl;
            }
        }
        if(started == true){
            std::cout << maxSizeInByte - 1 << "]" << std::endl;
        }
    }


    Range getRange(int index) {
        return ranges[index];
    }

private:
    std::vector<Range> ranges;
    unsigned char * array;
    int arraySizeInBytes;
    int maxSizeInByte;

};


#endif //MMSEQS_INTERVALARRAY_H
