#ifndef SEQUENCE_H
#define SEQUENCE_H

// Wrote by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger Martin.Steinegger@campus.lmu.de
// 
// Represents a database sequence object, holds its representation in the int array form.
//


#include <string>
#include <cstring>

class Sequence
{
public:
    Sequence(size_t maxLen,int* aa2int,char* int2aa);  
    
    ~Sequence();

    // Map char -> int
    void mapSequence(const char *seq);

    // checks if there is still a k-mer left 
    bool hasNextKmer(size_t kmerSize);

    // returns next k-mer
    const int* nextKmer(size_t kmerSize);

    // resets the sequence position pointer to the start of the sequence
    void reset() { currItPos = -1; }

    // set the id of the sequence object
    void setId(size_t id);

    size_t id;
    size_t L;              // length of sequence
    int * int_sequence;    // int sequence 
    
    // current iterator position
    size_t currItPos;
    size_t maxLen;
private:

    void print(); // for debugging 
    int  * aa2int; // ref to mapping from aa -> int
    char * int2aa; // ref mapping from int -> aa
};
#endif
