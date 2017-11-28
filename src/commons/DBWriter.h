#ifndef DBWRITER_H
#define DBWRITER_H

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de
// 
// Manages ffindex DB write access.
// For parallel write access, one ffindex DB per thread is generated.
// After the parallel calculation is done, all ffindexes are merged into one.
//

#include <string>
#include <vector>
#include "DBReader.h"

template <typename T> class DBReader;

class DBWriter {
    public:
        static const size_t ASCII_MODE = 0;
        static const size_t BINARY_MODE = 1;
        static const size_t LEXICOGRAPHIC_MODE = 2;


        DBWriter(const char* dataFileName, const char* indexFileName, unsigned int threads = 1, size_t mode = ASCII_MODE);

        ~DBWriter();

        void open(size_t bufferSize = 64 * 1024 * 1024);

        void close(int dbType = -1);
    
        char* getDataFileName() { return dataFileName; }
    
        char* getIndexFileName() { return indexFileName; }

        void writeStart(unsigned int thrIdx = 0);
        void writeAdd(const char* data, size_t dataSize, unsigned int thrIdx = 0);
        void writeEnd(unsigned int key, unsigned int thrIdx = 0, bool addNullByte = true);

        void writeData(const char *data, size_t dataSize, unsigned int key, unsigned int threadIdx = 0, bool addNullByte = true);

        size_t indexToBuffer(char *buff1, unsigned int key, size_t offsetStart, size_t len);

        void alignToPageSize();

        void mergeFiles(DBReader<unsigned int>& qdbr,
                        const std::vector<std::pair<std::string, std::string> >& files,
                        const std::vector<std::string>& prefixes);

        void sortDatafileByIdOrder(DBReader<unsigned int>& qdbr);

        static void mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                                 const std::vector<std::pair<std::string, std::string>> &files,
                                 bool lexicographicOrder = false);

        static void mergeResults(const char *outFileName, const char *outFileNameIndex,
                                 const char **dataFileNames, const char **indexFileNames,
                                 unsigned long fileCount, bool lexicographicOrder = false);

        void mergeFilePair(const char *inData1, const char *inIndex1, const char *inData2, const char *inIndex2);


private:
    template <typename T>
    static void writeIndex(FILE *outFile, size_t indexSize, T *index, unsigned int *seqLen);

    void checkClosed();

    char* dataFileName;
    char* indexFileName;

    FILE** dataFiles;
    char** dataFilesBuffer;
    size_t bufferSize;
    FILE** indexFiles;

    char** dataFileNames;
    char** indexFileNames;

    size_t* starts;
    size_t* offsets;

    const unsigned int threads;
    const size_t mode;

    bool closed;

    std::string datafileMode;


};

#endif
