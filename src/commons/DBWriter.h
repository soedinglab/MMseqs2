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
#include <zstd/lib/zstd.h>
#include "DBReader.h"

template <typename T> class DBReader;

class DBWriter {
public:


    DBWriter(const char* dataFileName, const char* indexFileName, unsigned int threads, size_t mode, int dbtype);

    ~DBWriter();

    void open(size_t bufferSize = 64 * 1024 * 1024);

    void close(bool merge = false);

    char* getDataFileName() { return dataFileName; }

    char* getIndexFileName() { return indexFileName; }

    void writeStart(unsigned int thrIdx = 0);
    size_t writeAdd(const char* data, size_t dataSize, unsigned int thrIdx = 0);
    void writeEnd(unsigned int key, unsigned int thrIdx = 0, bool addNullByte = true, bool addIndexEntry = true);

    void writeData(const char *data, size_t dataSize, unsigned int key, unsigned int threadIdx = 0, bool addNullByte = true, bool addIndexEntry = true);

    static size_t indexToBuffer(char *buff1, unsigned int key, size_t offsetStart, size_t len);

    void alignToPageSize();

    void mergeFiles(DBReader<unsigned int>& qdbr,
                    const std::vector<std::pair<std::string, std::string> >& files,
                    const std::vector<std::string>& prefixes);

    void sortDatafileByIdOrder(DBReader<unsigned int>& qdbr);

    static void mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                             const std::vector<std::pair<std::string, std::string>> &files,
                             bool lexicographicOrder = false);


    static void mergeResultsIndexOnly(const char *outFileName, const char *outFileNameIndex,
                                      const char **dataFileNames, const char **indexFileNames,
                                      const unsigned long fileCount, const bool lexicographicOrder = false);

    static void mergeIndex(const char **indexFileNames, std::vector<size_t> threadDataFileSizes,
                           const unsigned long fileCount);

    static void sortIndex(const char *inFileNameIndex, const char *outFileNameIndex, const bool lexicographicOrder);


    static void mergeResultsNormal(const char *outFileName, const char *outFileNameIndex,
                                   const char **dataFileNames, const char **indexFileNames,
                                   unsigned long fileCount, bool lexicographicOrder = false);

    void mergeFilePair(const std::vector<std::pair<std::string, std::string>> fileNames);

    void writeIndexEntry(unsigned int key, size_t offset, size_t length, unsigned int thrIdx);

    static void writeDbtypeFile(const char* path, int dbtype, bool isCompressed);

    size_t getStart(unsigned int threadIdx){
        return starts[threadIdx];
    }

    size_t getOffset(unsigned int threadIdx){
        return offsets[threadIdx];
    }

    template <typename T>
    static void writeIndex(FILE *outFile, size_t indexSize, T *index, unsigned int *seqLen);

    template <typename T>
    static void writeIndexEntryToFile(FILE *outFile, char *buff1, T &index,  unsigned int seqLen);

private:
    size_t addToThreadBuffer(const void *data, size_t itmesize, size_t nitems, int threadIdx);
    void writeThreadBuffer(unsigned int idx, size_t dataSize);

    void checkClosed();

    char* dataFileName;
    char* indexFileName;

    FILE** dataFiles;
    char** dataFilesBuffer;
    size_t bufferSize;
    FILE** indexFiles;

    char** dataFileNames;
    char** indexFileNames;
    char** compressedBuffers;
    size_t * compressedBufferSizes;
    char** threadBuffer;
    size_t * threadBufferSize;
    size_t * threadBufferOffset;

    size_t* starts;
    size_t* offsets;
    int* state;
    static const int INIT_STATE=0;
    static const int NOTCOMPRESSED=1;
    static const int COMPRESSED=2;

    ZSTD_CStream** cstream;

    const unsigned int threads;
    const size_t mode;
    int dbtype;

    bool closed;

    std::string datafileMode;


};

#endif
