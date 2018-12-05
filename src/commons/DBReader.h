#ifndef DBREADER_H
#define DBREADER_H

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de
//
// Manages DB read access.
//

#include <cstddef>
#include <utility>
#include <string>
#include "Sequence.h"
#include "Parameters.h"

#define ZSTD_STATIC_LINKING_ONLY // ZSTD_findDecompressedSize
#include <zstd/lib/zstd.h>

template <typename T>
class DBReader {

public:
    struct Index {
        T id;
        size_t offset;
        static bool compareById(const Index& x, const Index& y){
            return (x.id <= y.id);
        }
    };
    // = USE_DATA|USE_INDEX
    DBReader(const char* dataFileName, const char* indexFileName, int threads, int mode);

    DBReader(Index* index, unsigned int *seqLens, size_t size, size_t aaDbSize, T lastKey, int dbType, unsigned int maxSeqLen);

    void setDataFile(const char* dataFileName);

    virtual ~DBReader();

    bool open(int sort);

    void close();

    const char* getDataFileName() { return dataFileName; }

    const char* getIndexFileName() { return indexFileName; }

    size_t getAminoAcidDBSize(){ return aaDbSize; }

    char* getData(size_t id, int thrIdx);

    char* getDataCompressed(size_t id, int thrIdx);

    char* getDataUncompressed(size_t id);

    void touchData(size_t id);

    char* getDataByDBKey(T key, int thrIdx);

    size_t getSize();

    T getDbKey(size_t id);

    unsigned int * getSeqLens();

    size_t getSeqLens(size_t id);

    size_t maxCount(char c);

    void remapData();

    size_t bsearch(const Index * index, size_t size, T value);

    // does a binary search in the ffindex and returns index of the entry with dbKey
    // returns UINT_MAX if the key is not contained in index
    size_t getId (T dbKey);

    static const int NOSORT = 0;
    static const int SORT_BY_LENGTH = 1;
    static const int LINEAR_ACCCESS = 2;
    static const int SORT_BY_ID     = 3;
    static const int SORT_BY_LINE   = 4; // the local IDs correspond to the line number in the original index file
    static const int SHUFFLE        = 5;
    static const int HARDNOSORT = 6; // do not even sort by ids.
    static const int SORT_BY_ID_OFFSET = 7;


    static const int USE_INDEX    = 0;
    static const int USE_DATA     = 1;
    static const int USE_WRITABLE = 2;
    static const int USE_FREAD    = 4;

    // compressed
    static const int UNCOMPRESSED    = 0;
    static const int COMPRESSED     = 1;


    const char * getData(){
        return data;
    }

    size_t getDataSize(){
        return dataSize;
    }

    char *mmapData(FILE *file, size_t *dataSize);

    bool readIndex(char *indexFileName, Index *index, unsigned int *entryLength);

    void readIndexId(T* id, char * line, char** cols);

    void readMmapedDataInMemory();

    void mlock();

    void sortIndex(bool isSortedById);

    void unmapData();

    size_t getDataOffset(T i);


    Index* getIndex() {
        return index;
    }

    Index* getIndex(size_t id) {
        return index + local2id[id];
    }
    

    void printMagicNumber();
    
    T getLastKey();

    static size_t indexMemorySize(const DBReader<unsigned int> &idx);

    static char* serialize(const DBReader<unsigned int> &idx);

    static DBReader<unsigned int> *unserialize(const char* data);

    static int parseDbType(const char *name);

    int getDbtype(){
        return dbtype;
    }

    const char* getDbTypeName() {
        return getDbTypeName(dbtype);
    }

    static const char* getDbTypeName(int dbtype) {
        switch(dbtype) {
            case Parameters::DBTYPE_AMINO_ACIDS: return "Aminoacid";
            case Parameters::DBTYPE_NUCLEOTIDES: return "Nucleotide";
            case Parameters::DBTYPE_HMM_PROFILE: return "Profile";
            case Parameters::DBTYPE_PROFILE_STATE_SEQ: return "Profile state";
            case Parameters::DBTYPE_PROFILE_STATE_PROFILE: return "Profile profile";
            default: return "Unknown";
        }
    }

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
            if(lhs.second > rhs.second)
                return true;
            if(rhs.second > lhs.second)
                return false;
            if(lhs.first < rhs.first  )
                return true;
            if(rhs.first < lhs.first )
                return false;
            return false;
        }
    };

    struct comparePairByIdAndOffset {
        bool operator() (const std::pair<unsigned int, Index>& lhs, const std::pair<unsigned int, Index>& rhs) const{
            if(lhs.second.id < rhs.second.id)
                return true;
            if(rhs.second.id < lhs.second.id)
                return false;
            if(lhs.second.offset < rhs.second.offset  )
                return true;
            if(rhs.second.offset < lhs.second.offset )
                return false;
            return false;
        }
    };


    struct comparePairByOffset{
        bool operator() (const std::pair<unsigned int, size_t >& lhs, const std::pair<unsigned int, size_t >& rhs) const{
            return (lhs.second < rhs.second);
        }
    };

    void setData(char *data, size_t dataSize);

    void setMode(const int mode);

    size_t getOffset(size_t id);

    size_t findNextOffsetid(size_t id);

    int isCompressed(){
        return isCompressed(dbtype);
    }


    static int isCompressed(int dbtype);


private:

    void checkClosed();

    char* data;

    int threads;

    int dataMode;

    char* dataFileName;

    char* indexFileName;

    // number of entries in the index
    size_t size;
    // size of all data stored in ffindex
    size_t dataSize;
    // amino acid size
    size_t aaDbSize;
    // Last Key in Index
    T lastKey;
    // max seqLen
    unsigned int maxSeqLen;
    // flag to check if db was closed
    int closed;
    // stores the dbtype (if dbtype file exists)
    int dbtype;
    int compression;
    char ** compressedBuffers;
    size_t * compressedBufferSizes;
    ZSTD_DStream ** dstream;

    Index * index;

    unsigned int * seqLens;
    unsigned int * id2local;
    unsigned int * local2id;

    bool dataMapped;
    int accessType;

    bool externalData;

    bool didMlock;

    // needed to prevent the compiler from optimizing away the loop
    char magicBytes;

};

#endif
