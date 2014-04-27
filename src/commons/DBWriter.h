#ifndef DBWRITER_H
#define DBWRITER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
// 
// Manages ffindex DB write access. 
// For parallel write access, one ffindex DB per thread is generated. 
// After the parallel calculation is done, all ffindexes are merged into one.
//
#include "DBReader.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <cstdlib>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <stdio.h>

class DBWriter {
    public:

        DBWriter(const char* dataFileName, const char* indexFileName, int maxThreadNum = 1);

        ~DBWriter();

        void open();

        int close();
    
        char* getDataFileName() { return dataFileName; }
    
        char* getIndexFileName() { return indexFileName; }

        void write(char* data, int dataSize, char* key, int threadIdx = 0);
    
        void mergeFiles(DBReader * qdbr,
                        std::vector<std::pair<std::string, std::string> > files,
                        size_t maxLineLength);
    
        static void errorIfFileExist(const char * file);

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
