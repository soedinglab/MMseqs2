//
// Created by mad on 12/14/15.
//

#include <string.h>
#include "SequenceIndex.h"

SequenceIndex::SequenceIndex(size_t dbSize, size_t entrySize) {
    sequence = new char*[dbSize + 1];
    data = new char[entrySize + 1];
    currWritePos = data;
    sequence[0] = data;
    sequence[dbSize] = &data[entrySize];
}

SequenceIndex::~SequenceIndex() {
    delete [] sequence;
    delete [] data;
}

void SequenceIndex::addSequence(Sequence * seq) {
    sequence[seq->getId()] = currWritePos;
    for(size_t pos = 0; pos < seq->L; pos++){
        sequence[seq->getId()][pos] = seq->int_sequence[pos];
    }
    currWritePos = currWritePos + seq->L;
}

std::pair<char *, unsigned int> SequenceIndex::getSequence(size_t id) {
    unsigned int N = (sequence[id + 1] - sequence[id])/ sizeof(char);
    return std::make_pair(sequence[id], N);
}