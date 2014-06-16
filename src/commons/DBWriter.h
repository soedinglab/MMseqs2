#ifndef DBWRITER_H
#define DBWRITER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Manages ffindex DB write access. 
// For parallel write access, one ffindex DB per thread is generated. 
// After the parallel calculation is done, all ffindexes are merged into one.
//

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <cstdlib>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <sys/mman.h>

class DBWriter {
    public:

        DBWriter(const char* dataFileName, const char* indexFileName, int maxThreadNum = 1);

        ~DBWriter();

        void open();

        int close();
    
        char* getDataFileName() { return dataFileName; }
    
        char* getIndexFileName() { return indexFileName; }

        void write(char* data, int dataSize, char* key, int threadIdx = 0);

    private:

        void initFFIndexWrite(const char* dataFileName, const char* indexFileName, FILE** dataFile, FILE** indexFile);

        void checkClosed();

        char* dataFileName;

        char* indexFileName;

        FILE** dataFiles;

        FILE** indexFiles;

        char** dataFileNames;

        char** indexFileNames;

        size_t* offsets;

        int maxThreadNum;

        int closed;
};

#endif
