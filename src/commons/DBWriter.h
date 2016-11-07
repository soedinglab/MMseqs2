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

        DBWriter(const char* dataFileName, const char* indexFileName, unsigned int threads = 1, size_t mode = ASCII_MODE);

        ~DBWriter();

        void open(size_t bufferSize = 64 * 1024 * 1024);

        void close();
    
        char* getDataFileName() { return dataFileName; }
    
        char* getIndexFileName() { return indexFileName; }

        void writeData(const char *data, size_t dataSize, const char *key, unsigned int threadIdx = 0);

        void mergeFiles(DBReader<unsigned int>& qdbr,
                        const std::vector<std::pair<std::string, std::string> >& files,
                        const std::vector<std::string>& prefixes);

        void sortDatafileByIdOrder(DBReader<unsigned int>& qdbr);

        static void mergeResults(const char *outFileName, const char *outFileNameIndex, const char **dataFileNames,
                                        const char **indexFileNames, unsigned int fileCount);

        void mergeFilePair(const char *inData1, const char *inIndex1, const char *inData2, const char *inIndex2);


private:

    struct IndexType {
        typedef DBReader<unsigned int>::Index int_type;
        typedef DBReader<std::string>::Index string_type;
    };
    static void writeIndex(FILE *outFile, IndexType::string_type *index, size_t indexSize, unsigned int *seqLen);
    static void writeIndex(FILE *outFile, size_t indexSize, IndexType::int_type *index, unsigned int *seqLen);


    void checkClosed();

    char* dataFileName;
    char* indexFileName;

    FILE** dataFiles;
    char** dataFilesBuffer;
    size_t bufferSize;
    FILE** indexFiles;

    char** dataFileNames;
    char** indexFileNames;

    size_t* offsets;

    unsigned int threads;

    bool closed;

    std::string datafileMode;


};

#endif
