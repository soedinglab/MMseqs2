//
// Created by mad on 12/14/15.
//

#ifndef MMSEQS_SEQUENCEINDEX_H
#define MMSEQS_SEQUENCEINDEX_H


#include <stddef.h>
#include <Sequence.h>

class SequenceLookup {

public:
    SequenceLookup(size_t dbSize, size_t entrySize);

    ~SequenceLookup();

    // add sequence to index
    void addSequence(Sequence * seq);

    // get sequence data
    std::pair<const unsigned char *, const unsigned int> getSequence(size_t id);

    const char *getData();

    int64_t getDataSize();

    size_t getSequenceCount();

    void initLookupByExternalData(FILE *datafile, size_t seqSizesOffset, size_t seqDataOffset);

private:

    // data contains sequence data
    char * data;

    // pointer to the start of the sequence
    char ** sequence;

    // write position
    char * currWritePos;
    size_t sequenceCount;
    size_t dataSize;

};


#endif //MMSEQS_SEQUENCEINDEX_H
