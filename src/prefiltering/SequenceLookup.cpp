//
// Created by mad on 12/14/15.
//
#include <new>
#include <cstring>
#include <sys/mman.h>
#include "Debug.h"
#include "Util.h"
#include "SequenceLookup.h"

SequenceLookup::SequenceLookup(size_t dbSize,
                               size_t entrySize) {
    sequenceCount = dbSize;
    sequence = new(std::nothrow) char*[sequenceCount + 1];
    Util::checkAllocation(sequence, "Could not allocate sequence memory in SequenceLookup");
    dataSize = entrySize;
    data = new(std::nothrow)     char[dataSize + 1];
    Util::checkAllocation(data, "Could not allocate data memory in SequenceLookup");

    currWritePos = data;
    sequence[0] = data;
    sequence[dbSize] = &data[entrySize];
    externalData = false;
}

SequenceLookup::SequenceLookup(size_t dbSize) {
    sequenceCount = dbSize;
    data = NULL;
    sequence = new(std::nothrow) char*[sequenceCount + 1];
    Util::checkAllocation(sequence, "Could not allocate sequence memory in SequenceLookup");
    externalData = true;
}

SequenceLookup::~SequenceLookup() {
    delete [] sequence;
    if(externalData == false){
        delete [] data;
    }else{
        munlock(data, dataSize);
    }
}

void SequenceLookup::addSequence(Sequence *  seq, size_t offset){
    sequence[seq->getId()] = data + offset;
    for(int pos = 0; pos < seq->L; pos++){
        unsigned char aa = seq->int_sequence[pos];
        sequence[seq->getId()][pos] = aa;
    }
}

void SequenceLookup::addSequence(Sequence * seq) {
    addSequence(seq, currWritePos - data);
    currWritePos = currWritePos + seq->L;
}

std::pair<const unsigned char *, const unsigned int> SequenceLookup::getSequence(size_t id) {
    const unsigned int N = (sequence[id + 1] - sequence[id])/ sizeof(char);
    return std::make_pair(( const unsigned char *)sequence[id], N);
}

const char *SequenceLookup::getData() {
    return data;
}

int64_t SequenceLookup::getDataSize() {
    return dataSize;
}

size_t SequenceLookup::getSequenceCount() {
    return sequenceCount;
}

void SequenceLookup::initLookupByExternalData(char * seqData,
                                              unsigned int * seqSizes) {
    // copy data to data element
    data = seqData;
    char * it = data;
    magicByte = 0;
    // set the pointers
    for (size_t i = 0; i < sequenceCount; i++){
        const unsigned int entriesCount = (unsigned int) seqSizes[i];
        sequence[i] = (char *) it;
        it += entriesCount;
        magicByte += sequence[i][0]; // this will read 4kb
    }
    dataSize = it - data;
    mlock(data, dataSize);
    sequence[sequenceCount] = it;
}