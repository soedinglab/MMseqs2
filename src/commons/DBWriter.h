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
#include <sys/mman.h>

class DBWriter {
    public:
    
        static const size_t ASCII_MODE = 0;
        static const size_t BINARY_MODE = 1;

        DBWriter(const char* dataFileName, const char* indexFileName, int maxThreadNum = 1, size_t mode = ASCII_MODE);

        ~DBWriter();

        void open();

        int close();
    
        char* getDataFileName() { return dataFileName; }
    
        char* getIndexFileName() { return indexFileName; }

        void write(char* data, int64_t dataSize, char* key, int threadIdx = 0);
    
        void mergeFiles(DBReader * qdbr,
                        std::vector<std::pair<std::string, std::string> > files,
                        size_t maxLineLength);

        void swapResults(std::string inputDb, size_t splitSize);

        static void errorIfFileExist(const char * file);

        void writeFile(FILE * file, char* key, int thrIdx);

        void sortDatafileByIdOrder(DBReader *qdbr);

        static void initFFIndexWrite(const char* dataFileName, const char* indexFileName, const char * dataFileMode, FILE** dataFile, FILE** indexFile);

        static void mergeFFindexFile(const char * outFileName, const char * outFileNameIndex, const char * datafileMode,
                                     const char **dataFileNames, const char **indexFileNames, int fileCount );

        void mergeFilePair(const char *inData1, const char *inIndex1, const char *inData2, const char *inIndex2);
private:

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
    
        std::string datafileMode;

};

#endif
