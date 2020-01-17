//
// Created by mad on 12/14/15.
//

#ifndef MMSEQS_SEQUENCEINDEX_H
#define MMSEQS_SEQUENCEINDEX_H


#include <cstddef>
#include "Sequence.h"

class SequenceLookup {
public:
    SequenceLookup(size_t dbSize, size_t entrySize);
    SequenceLookup(size_t dbSize);
    ~SequenceLookup();

    // add sequence at offset
    void addSequence(unsigned char *seq, int L, size_t index, size_t offset);

    // add sequence to index
    void addSequence(Sequence * seq);

    // get sequence data
    std::pair<const unsigned char *, const unsigned int> getSequence(size_t id);

    const char *getData();

    int64_t getDataSize();

    size_t getSequenceCount();

    size_t *getOffsets();

    void initLookupByExternalData(char *seqData, size_t dataSize, size_t *seqOffsets);
    void initLookupByExternalDataCopy(char *seqData, size_t *seqOffsets);

private:
    size_t sequenceCount;

    // data contains sequence data
    char *data;
    size_t dataSize;

    size_t *offsets;

    // write position
    size_t currentIndex;
    size_t currentOffset;

    // if data are read from mmap
    bool externalData;
};


#endif //MMSEQS_SEQUENCEINDEX_H
