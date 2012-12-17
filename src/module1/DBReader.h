#ifndef DBREADER_H
#define DBREADER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Manages ffindex DB read access.
//

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <cstdlib>
#include <iostream>


class DBReader {

    public:

        DBReader(char* dataFileName, char* indexFileName);

        void open();

        void close();

        char* getData(int id);

        char* getDataByDBKey(char* key);

        size_t getSize();

        char* getDbKey(int id);

    private:
        char* dataFileName;

        char* indexFileName;

        FILE* dataFile;

        FILE* indexFile;
        
        char* data;
        
        ffindex_index_t* index; 
        // number of entries in the index
        size_t size;
        // size of all data stored in ffindex
        size_t dataSize;
        // mapping id -> ffindex_entry
        ffindex_entry_t** id2entry;

};


#endif
