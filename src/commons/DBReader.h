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
#include <fstream>
#include <algorithm>


class DBReader {

    public:

        DBReader(const char* dataFileName, const char* indexFileName);

        void open(int sort = 0);

        void close();

        char* getData(size_t id);

        char* getDataByDBKey(char* key);

        size_t getSize();

        char* getDbKey(size_t id);

        unsigned short* getSeqLens();

    private:

        void sort(size_t* ids, size_t* workspace);

        void merge(size_t* ids, size_t iLeft, size_t iRight, size_t iEnd, size_t* workspace);

        void calcLocalIdMapping();

        size_t* id2local;

        size_t* local2id;

        const char* dataFileName;

        const char* indexFileName;

        unsigned short* seqLens;

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
