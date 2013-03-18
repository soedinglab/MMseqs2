#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <iostream>
#include <algorithm>

#include "Sequence.h"
#include "Indexer.h"

class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize);

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

//        void checkSizeAndCapacity();

 //       void reduceMemoryUsage();

    private:
        
        int ipow (int base, int exponent);

        // Index table: contains pointers to the arrays (stored in DynamicArray structure) of DB sequences containing a certain k-mer
//        DynamicArray** table

        int** table;

        int* sizes;

        // only for init: current position in the DB id array of the index table where the next sequence id can be written
        int* currPos;

        Indexer* idxer;

        int alphabetSize;

        int kmerSize;

        Sequence* s;

};

#endif
