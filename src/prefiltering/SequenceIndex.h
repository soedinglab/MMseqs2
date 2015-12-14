//
// Created by mad on 12/14/15.
//

#ifndef MMSEQS_SEQUENCEINDEX_H
#define MMSEQS_SEQUENCEINDEX_H


#include <stddef.h>
#include <Sequence.h>

class SequenceIndex {

public:
    SequenceIndex(size_t dbSize, size_t entrySize);

    ~SequenceIndex();

    // add sequence to index
    void addSequence(Sequence * seq);

    // get sequence data
    std::pair<char *, unsigned int> getSequence(size_t id);

private:

    // data contains sequence data
    char * data;

    // pointer to the start of the sequence
    char ** sequence;

    // write position
    char * currWritePos;
};


#endif //MMSEQS_SEQUENCEINDEX_H
