//
// Created by mad on 12/14/15.
//

#include <string.h>
#include "SequenceLookup.h"

SequenceLookup::SequenceLookup(size_t dbSize, size_t entrySize) {
    sequence = new char*[dbSize + 1];
    data = new char[entrySize + 1];
    currWritePos = data;
    sequence[0] = data;
    sequence[dbSize] = &data[entrySize];
}

SequenceLookup::~SequenceLookup() {
    delete [] sequence;
    delete [] data;
}

void SequenceLookup::addSequence(Sequence * seq) {
    sequence[seq->getId()] = currWritePos;
    for(size_t pos = 0; pos < seq->L; pos++){
        unsigned char aa = seq->int_sequence[pos];
        sequence[seq->getId()][pos] = aa;
    }
    currWritePos = currWritePos + seq->L;
}

std::pair<const unsigned char *, const unsigned int> SequenceLookup::getSequence(size_t id) {
    const unsigned int N = (sequence[id + 1] - sequence[id])/ sizeof(char);
    return std::make_pair(( const unsigned char *)sequence[id], N);
}