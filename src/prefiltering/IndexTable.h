#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>

#include "../commons/Sequence.h"
#include "Indexer.h"

class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize, int skip);

        ~IndexTable();

        // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
        void addKmerCount (Sequence* s);

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        void removeDuplicateEntries();

        // init the arrays for the sequence lists 
        void init();

        // get list of DB sequences containing this k-mer
        int* getDBSeqList (int kmer, int* matchedListSize);

        void print();

        // alphabetSize**kmerSize
        int tableSize;

    private:
        int ipow (int base, int exponent);

        // Index table: contains pointers to the point in the entries array where starts the list of sequence ids for a certain k-mer
        int** table;

        // Index table entries: ids of sequences containing a certain k-mer, stored sequentially in the memory
        int* entries;

        int* sizes;

        // only for init: current position in the DB id array of the index table where the next sequence id can be written
        int* currPos;

        Indexer* idxer;
    
        int alphabetSize;

        int kmerSize;

        Sequence* s;

        // number of skipped k-mers
        int skip;
    
        // amount of sequences in Index
        size_t size;
    

        // number of entries in all sequence lists
        int64_t tableEntriesNum;
};

#endif
