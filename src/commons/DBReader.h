#ifndef DBREADER_H
#define DBREADER_H

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de
//
// Manages DB read access.
//

#include <cstdlib>
#include <utility>
#include <string>

class DBReader {

public:
    struct Index {
        unsigned int id;
        char * data;
        static bool compareById(Index x, Index y){
            return (x.id <= y.id);
        }
    };
    DBReader(const char* dataFileName, const char* indexFileName, int mode = DATA_AND_INDEX);

    ~DBReader();

    void open(int sort);

    void close();

    char* getDataFileName() { return dataFileName; }

    char* getIndexFileName() { return indexFileName; }

    size_t getAminoAcidDBSize(){ return aaDbSize; }

    char* getData(size_t id);

    char* getDataByDBKey(const char* key);

    size_t getSize();

    std::string getDbKey(size_t id);

    unsigned int * getSeqLens();

    size_t getSeqLens(size_t id);

    void remapData();

    size_t bsearch(const Index * index, size_t size, unsigned int value);

    // does a binary search in the ffindex and returns index of the entry with dbKey
    // returns UINT_MAX if the key is not contained in index
    size_t getId (const char* dbKey);


    static const int NOSORT = 0;
    static const int SORT = 1;

    static const int INDEXONLY = 1;
    static const int DATA_AND_INDEX = 0;

    const char * getData(){
        return data;
    }

    size_t getDataSize(){
        return dataSize;
    }

    static char *mmapData(FILE *file, size_t *dataSize);

    void readIndex(char *indexFileName, Index *index, char *data, unsigned int *entryLength);

    static size_t countLine(const char *name);

private:


    struct compareIndexLengthPairById {
        bool operator() (const std::pair<Index, unsigned  int>& lhs, const std::pair<Index, unsigned  int>& rhs) const{
            return (lhs.first.id < rhs.first.id);
        }
    };

    void checkClosed();

    char* dataFileName;

    char* indexFileName;

    FILE* dataFile;

    int dataMode;

    char* data;

    unsigned int * ids;
    // number of entries in the index
    size_t size;
    // size of all data stored in ffindex
    size_t dataSize;
    // amino acid size
    size_t aaDbSize;
    // flag to check if db was closed
    int closed;

    Index * index;
    unsigned int *seqLens;
    size_t maxId;
};


#endif
