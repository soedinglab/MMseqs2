#ifndef INDEX_TABLE_LOCAL_H
#define INDEX_TABLE_LOCAL_H

// Written by Martin Steinegger martin.steinegger@campus.lmu.de
// 
// Index table for local sequence search stores the list of DB sequences and position of certain k-mer

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>


#include "IndexTable.h"
#include "Sequence.h"
#include "Indexer.h"


class IndexTableLocal : public virtual IndexTable {

    public:

        IndexTableLocal (int alphabetSize, int kmerSize, int skip);

        ~IndexTableLocal();

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        // remove dublicate entries from list
        void removeDuplicateEntries();

        // print the index table (debugging)
        void print(char * int2aa);
    
};

#endif
