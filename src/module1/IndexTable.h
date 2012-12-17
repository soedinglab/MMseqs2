#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <list>
#include <iostream>

#include "Sequence.h"
#include "Indexer.h"

class IndexTable {

    public:

        IndexTable (size_t alphabetSize, size_t kmerSize);

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        // convert the sequence lists for each k-mer to array
        void init();

        // get list of DB sequences containing this k-mer
        size_t* getDBSeqList (size_t kmer, size_t* matchedListSize);

    private:
        
        size_t ipow (size_t base, size_t exponent);

        // Index table: contains posize_ters to the arrays of DB sequences containing a certain k-mer
        size_t** table;

        // Index table before init: DB sequences for a certain k-mer are stored in lists.
        // In init(), these lists are copied size_to fixed size arrays and deleted.
        std::list<size_t>** tableDummy;

        Indexer* idxer;

        size_t alphabetSize;

        size_t kmerSize;

        // ALPHABET_SIZE**KMER_SIZE
        size_t tableSize;

        // for each possible k-mer, the size of the corresponding sequences list
        size_t* listSizes;

};

#endif
