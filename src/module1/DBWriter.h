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

class DBWriter {
    public:

        DBWriter(const char* dataFileName, const char* indexFileName, int maxThreadNum);

        void open();

        void close();

        void write(char* data, int dataSize, char* key, int threadIdx);

    private:

        void initFFIndexWrite(char* dataFileName, char* indexFileName, FILE** dataFile, FILE** indexFile);

        const char* dataFileName;

        const char* indexFileName;

        FILE** dataFiles;

        FILE** indexFiles;

        size_t* offsets;

        int maxThreadNum;
};

#endif
