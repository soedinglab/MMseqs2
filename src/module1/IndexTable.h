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
#include "DynamicArray.h"

class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize);

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        // convert the sequence lists for each k-mer to array
        void init();

        // get list of DB sequences containing this k-mer
        int* getDBSeqList (int kmer, int* matchedListSize);

        // alphabetSize**kmerSize
        int tableSize;

        void checkSizeAndCapacity();

        void reduceMemoryUsage();

    private:
        
        int ipow (int base, int exponent);

        // Index table: contains pointers to the arrays (stored in DynamicArray structure) of DB sequences containing a certain k-mer
        DynamicArray** table;

        Indexer* idxer;

        int alphabetSize;

        int kmerSize;

        Sequence* s;

};

#endif
