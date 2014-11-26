#ifndef SEQUENCE_H
#define SEQUENCE_H
#define READ_BUFFER_SIZE 20971520

// Written by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger Martin.Steinegger@campus.lmu.de
// 
// Represents a database sequence object, holds its representation in the int array form.
//


#include <string>
#include <cstring>
#include <cstdio>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../commons/BaseMatrix.h"      // for pBack
#include "../commons/ScoreMatrix.h"      // for ScoreMatrix
#include "KmerIterator.h"      // for ScoreMatrix

const int8_t seed_4[]        = {1, 1, 1, 1};
const int8_t seed_4_spaced[] = {1, 1, 1, 0, 1};
const int8_t seed_5[]        = {1, 1, 1, 1, 1};
//const char seed_5_spaced[] = {1, 1, 1, 0, 1, 1};
//const char seed_5_spaced[] = {1, 1, 0, 1, 0, 1, 1};// better than 111011
const int8_t seed_5_spaced[] = {1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1}; // just 0.001 %worse ROC5 but way faster
const int8_t seed_6[]        = {1, 1, 1, 1, 1, 1};
//const char seed_6_spaced[] = {1, 1, 1, 0, 1, 1, 0, 1};
const int8_t seed_6_spaced[] = {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}; // better than 11101101
const int8_t seed_7[]        = {1, 1, 1, 1, 1, 1, 1};
//const char seed_7_spaced[] = {1, 1, 1, 1, 0, 1, 0, 1, 0, 1};
const int8_t seed_7_spaced[] = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1};

class Sequence
{
    public:
        Sequence(size_t maxLen, int* aa2int, char* int2aa,
                 int seqType,
                 const unsigned int kmerSize = 0,
                 const unsigned int stepSize = 1,
                 const bool spaced = false,
                 BaseMatrix * subMat = NULL );
        ~Sequence();

        // Map char -> int
        void mapSequence(int id, char* dbKey, const char *seq);
    
        // Map Profile HMM
        void mapProfile(const char *data);

        // checks if there is still a k-mer left
        bool hasNextKmer();

        // returns next k-mer
        const int* nextKmer();

        // resets the sequence position pointer to the start of the sequence
        void resetCurrPos() { kmerIterator->reset(); }

        void print(); // for debugging

        int getId() { return id; }
    
        int getCurrentPosition() { return kmerIterator->getPosition(); }


        char* getDbKey() { return dbKey; }
    
        int getSeqType() { return seqType; }


        // reverse the sequence for the match statistics calculation
        void reverse();

        static const int AMINO_ACIDS = 0;
        static const int NUCLEOTIDES = 1;
        static const int HMM_PROFILE = 2;

    
        // Kmer Iterator
        KmerIterator * kmerIterator;
    
        // length of sequence
        int L;
        // each amino acid coded as integer
        int * int_sequence;  

        // Contains profile information
        short           * profile_score;
        unsigned int    * profile_index;
        size_t profile_row_size;
        static const size_t PROFILE_AA_SIZE = 20;
        ScoreMatrix ** profile_matrix;
    
        int  * aa2int; // ref to mapping from aa -> int
        char * int2aa; // ref mapping from int -> aa

        std::pair<const int8_t *, unsigned int> getSpacedPattern(bool spaced, unsigned int kmerSize);

    
    private:
        void mapProteinSequence(const char *seq);
        void mapNucleotideSequence(const char *seq);
        int id;
        char* dbKey;
        // AMINO_ACIDS or NUCLEOTIDES
        int seqType;
        // Matrix for Profile calculation
        BaseMatrix * subMat;
        // maximum possible length of sequence
        size_t maxLen;
        // read next kmer profile in profile_matrix
        void nextProfileKmer(const unsigned int * kmerPos);
        // kmer Size
        unsigned int kmerSize;
        // sequence window will be filled by newxtKmer (needed for spaced patterns)
        int * kmerWindow;
};
#endif
