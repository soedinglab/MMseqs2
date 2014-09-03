#ifndef INDEX_TABLE_GLOBAL_H
#define INDEX_TABLE_GLOBAL_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Index table stores the list of DB sequences containing a certain k-mer, for each k-mer. 
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>


#include "IndexTable.h"
#include "Sequence.h"
#include "Indexer.h"

class IndexTableGlobal : public virtual IndexTable {

    public:

        IndexTableGlobal (int alphabetSize, int kmerSize, int skip);

        ~IndexTableGlobal();

        // add k-mers of the sequence to the index table
        void addSequence (Sequence* s);

        // remove dublicate entries from list
        void removeDuplicateEntries();

        // print the index table (debugging)
        void print(char * int2aa);
    
    
};

#endif
