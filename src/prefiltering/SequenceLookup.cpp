//
// Created by mad on 12/14/15.
//

#include <cstring>
#include "Debug.h"
#include "Util.h"
#include "SequenceLookup.h"

SequenceLookup::SequenceLookup(size_t dbSize, size_t entrySize) {
    sequenceCount = dbSize;
    sequence = new char*[sequenceCount + 1];
    dataSize = entrySize;
    data = new char[dataSize + 1];
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
    for(int pos = 0; pos < seq->L; pos++){
        unsigned char aa = seq->int_sequence[pos];
        sequence[seq->getId()][pos] = aa;
    }
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

void SequenceLookup::initLookupByExternalData(FILE * datafile,
                                              size_t seqSizesOffset,
                                              size_t seqDataOffset) {
    // copy data to data element

    fseek (datafile, seqDataOffset, SEEK_SET);
    size_t errCode = fread(data, 1,  dataSize * sizeof(char),datafile);
    if (errCode != dataSize * sizeof(char)) {
        Debug(Debug::ERROR) << "IndexTable error while reading entries.\n";
        EXIT (EXIT_FAILURE);
    }
    fseek (datafile , seqSizesOffset, SEEK_SET);
    unsigned int *seqSizes = new unsigned int[sequenceCount];

    errCode = fread (seqSizes, 1, sequenceCount * sizeof(unsigned int), datafile);
    if (errCode != sequenceCount * sizeof(unsigned int)) {
        Debug(Debug::ERROR) << "IndexTable error while reading entries.\n";
        EXIT (EXIT_FAILURE);
    }

    char * it = data;
    // set the pointers
    for (size_t i = 0; i < sequenceCount; i++){
        const unsigned int entriesCount = (unsigned int) seqSizes[i];
        sequence[i] = (char *) it;
        it += entriesCount;
    }
    delete [] seqSizes;
}
