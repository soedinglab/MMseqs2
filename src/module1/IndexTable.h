#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <vector>
#include <iostream>

#include "Sequence.h"
#include "Indexer.h"

class IndexTable {

    public:

        IndexTable (int alphabetSize, int kmerSize);

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        // convert the sequence lists for each k-mer to array
        void init();

        // get list of DB sequences containing this k-mer
        int* getDBSeqList (int kmer, int* matchedListSize);

    private:
        
        int ipow (int base, int exponent);

        // Index table: contains pointers to the arrays of DB sequences containing a certain k-mer
        int** table;

        // Index table before init: DB sequences for a certain k-mer are stored in lists.
        // In init(), these lists are copied into fixed size arrays and deleted.
        std::vector<int>** tableDummy;

        Indexer* idxer;

        int alphabetSize;

        int kmerSize;

        // ALPHABET_SIZE**KMER_SIZE
        int tableSize;

        // for each possible k-mer, the size of the corresponding sequences list
        int* listSizes;

};

#endif
