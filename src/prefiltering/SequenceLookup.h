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
    void addSequence(Sequence * seq, size_t offset);

    // add sequence to index
    void addSequence(Sequence * seq);

    // get sequence data
    std::pair<const unsigned char *, const unsigned int> getSequence(size_t id);

    const char *getData();

    int64_t getDataSize();

    size_t getSequenceCount();

    void initLookupByExternalData(char * seqData,
                                  unsigned int * seqSizes);

private:

    // data contains sequence data
    char * data;

    // pointer to the start of the sequence
    char ** sequence;

    // write position
    char * currWritePos;
    size_t sequenceCount;
    size_t dataSize;

    // if data are read from mmap
    bool externalData;

    // magic byte to cache database
    size_t magicByte;

};


#endif //MMSEQS_SEQUENCEINDEX_H
