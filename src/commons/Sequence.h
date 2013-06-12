#ifndef SEQUENCE_H
#define SEQUENCE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger Martin.Steinegger@campus.lmu.de
// 
// Represents a database sequence object, holds its representation in the int array form.
//


#include <string>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <cstdlib>

typedef struct {
        float kmersPerPos;
        int dbMatches;
} statistics_t;

class Sequence
{
    public:
        Sequence(size_t maxLen,int* aa2int,char* int2aa);  

        ~Sequence();

        // Map char -> int
        void mapSequence(const char *seq);

        // checks if there is still a k-mer left 
        bool hasNextKmer(int kmerSize);

        // returns next k-mer
        const int* nextKmer(int kmerSize);

        // resets the sequence position pointer to the start of the sequence
        void resetCurrPos() { currItPos = -1; }

        void print(); // for debugging 

        int id;
        char* dbKey;

        int L;              // length of sequence
        int * int_sequence;    // int sequence 

        int  * aa2int; // ref to mapping from aa -> int
        char * int2aa; // ref mapping from int -> aa

        statistics_t* stats;
    private:
        // current iterator position
        int currItPos;
        size_t maxLen;
};
#endif
