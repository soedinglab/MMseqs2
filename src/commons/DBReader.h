#ifndef DBREADER_H
#define DBREADER_H

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de
//
// Manages DB read access.
//

#include <cstddef>
#include <utility>
#include <string>

template <typename T>
class DBReader {

public:
    struct Index {
        T id;
        size_t offset;
        static bool compareById(Index x, Index y){
            return (x.id <= y.id);
        }
    };
    DBReader(const char* dataFileName, const char* indexFileName, int mode = USE_DATA|USE_INDEX);

    virtual ~DBReader();

    void open(int sort);

    void close();

    const char* getDataFileName() { return dataFileName; }

    const char* getIndexFileName() { return indexFileName; }

    size_t getAminoAcidDBSize(){ return aaDbSize; }

    char* getData(size_t id);

    char* getDataByDBKey(T key);

    size_t getSize();

    T getDbKey(size_t id);

    unsigned int * getSeqLens();

    size_t getSeqLens(size_t id);

    void remapData();

    size_t bsearch(const Index * index, size_t size, T value);

    // does a binary search in the ffindex and returns index of the entry with dbKey
    // returns UINT_MAX if the key is not contained in index
    size_t getId (T dbKey);

    static const int NOSORT = 0;
    static const int SORT_BY_LENGTH = 1;
    static const int LINEAR_ACCCESS = 2;
    static const int SORT_BY_ID = 3;
    static const int SORT_BY_LINE = 4; // the local IDs correspond to the line number in the original index file

    static const int USE_INDEX    = 0;
    static const int USE_DATA     = 1;
    static const int USE_WRITABLE = 2;

    const char * getData(){
        return data;
    }

    size_t getDataSize(){
        return dataSize;
    }

    char *mmapData(FILE *file, size_t *dataSize);

    void readIndex(char *indexFileName, Index *index, char *data, unsigned int *entryLength);

    void readIndexId(T* id, char * line, char** cols);

    void readMmapedDataInMemory();

    void sortIndex();

    void unmapData();

    FILE* getDatafile(){
        return dataFile;
    }

    size_t getDataOffset(T i);


    Index* getIndex() {
        return index;
    }

    void printMagicNumber();
    
    T getLastKey();

private:

    struct compareIndexLengthPairById {
        bool operator() (const std::pair<Index, unsigned  int>& lhs, const std::pair<Index, unsigned  int>& rhs) const{
            return (lhs.first.id < rhs.first.id);
        }
    };
	
    struct compareIndexLengthPairByIdKeepTrack {
        bool operator() (const std::pair<Index, std::pair<size_t, unsigned int> >& lhs, const std::pair<Index, std::pair<size_t, unsigned int> >& rhs) const{
            return (lhs.first.id < rhs.first.id);
        }
    };

    struct comparePairBySeqLength {
        bool operator() (const std::pair<unsigned int, unsigned  int>& lhs, const std::pair<unsigned int, unsigned  int>& rhs) const{
            return (lhs.second > rhs.second);
        }
    };

    struct comparePairByOffset{
        bool operator() (const std::pair<unsigned int, size_t >& lhs, const std::pair<unsigned int, size_t >& rhs) const{
            return (lhs.second < rhs.second);
        }
    };

    void checkClosed();

    char* dataFileName;

    char* indexFileName;

    FILE* dataFile;

    int dataMode;

    char* data;

    // number of entries in the index
    size_t size;
    // size of all data stored in ffindex
    size_t dataSize;
    // amino acid size
    size_t aaDbSize;
    // flag to check if db was closed
    int closed;

    Index * index;

    unsigned int * seqLens;
    unsigned int * id2local;
    unsigned int * local2id;

    T lastKey;
    
    bool dataMapped;
    int accessType;

    // needed to avoid compiler to optimize away the loop
    size_t magicBytes;
};

#endif
