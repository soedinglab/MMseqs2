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
#include <climits>
#include <map>
#include <cstring>
#include <sys/mman.h>

struct StrCompare : public std::binary_function<const char*, const char*, bool> {
    public:
        bool operator() (const char* str1, const char* str2) const
        { return std::strcmp(str1, str2) < 0; }
};

class DBReader {

    public:

        DBReader(const char* dataFileName, const char* indexFileName, int mode = DATA_AND_INDEX);

        ~DBReader();

        void open(int sort);

        void close();

        char* getDataFileName() { return dataFileName; }

        char* getIndexFileName() { return indexFileName; }
    
        unsigned int getAminoAcidDBSize(){ return aaDbSize; }

        char* getData(size_t id);
    
        char* getDataByDBKey(char* key);
    
        size_t getSize();

        char* getDbKey(size_t id);
    
        FILE* getDataAsFile(char* key);

        unsigned short* getSeqLens();

        // does a binary search in the ffindex and returns index of the entry with dbKey
        // returns UINT_MAX if the key is not contained in index
        size_t getId (const char* dbKey);


        static const int NOSORT = 0;
        static const int SORT = 1;
    
        static const int INDEXONLY = 1;
        static const int DATA_AND_INDEX = 0;

    private:

        void sort(size_t* ids, size_t* workspace);

        void merge(size_t* ids, size_t iLeft, size_t iRight, size_t iEnd, size_t* workspace);

        void calcLocalIdMapping();

        void checkClosed();

        size_t* id2local;

        size_t* local2id;

        char* dataFileName;

        char* indexFileName;

        unsigned short* seqLens;

        FILE* dataFile;
    
        int dataMode;

        FILE* indexFile;
        
        char* data;
        
        ffindex_index_t* index; 
        // number of entries in the index
        size_t size;
        // size of all data stored in ffindex
        size_t dataSize;
        // mapping id -> ffindex_entry
        ffindex_entry_t** id2entry;
        // amino acid size
        unsigned int aaDbSize;
        // flag to check if db was closed
        int closed;

};


#endif
