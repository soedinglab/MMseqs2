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
        size_t dbMatches;
} statistics_t;

class Sequence
{
    public:
        Sequence(size_t maxLen,int* aa2int,char* int2aa, int seqType);  

        ~Sequence();

        // Map char -> int
        void mapSequence(int id, char* dbKey, const char *seq);

        // checks if there is still a k-mer left 
        bool hasNextKmer(int kmerSize);

        // returns next k-mer
        const int* nextKmer(int kmerSize);

        // resets the sequence position pointer to the start of the sequence
        void resetCurrPos() { currItPos = -1; }

        void print(); // for debugging 

        int getId() { return id; }

        char* getDbKey() { return dbKey; }

        // reverse the sequence for the match statistics calculation
        void reverse();

        static const int AMINO_ACIDS = 0;
        static const int NUCLEOTIDES = 1;

        // length of sequence
        int L;
        // each amino acid coded as integer
        int * int_sequence;  

        int  * aa2int; // ref to mapping from aa -> int
        char * int2aa; // ref mapping from int -> aa

        statistics_t* stats;

    private:
        void mapProteinSequence(const char *seq);
        void mapNucleotideSequence(const char *seq);
        
        int id;
        char* dbKey;
        // current iterator position
        int currItPos;
        // AMINO_ACIDS or NUCLEOTIDES
        int seqType;
        // maximum possible length of sequence
        size_t maxLen;
};
#endif
