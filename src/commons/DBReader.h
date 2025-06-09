#ifndef DBREADER_H
#define DBREADER_H

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de
//
// Manages DB read access.
//
#include "MemoryTracker.h"
#include <cstddef>
#include <utility>
#include <vector>
#include <string>
#include "Sequence.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"

#define ZSTD_STATIC_LINKING_ONLY // ZSTD_findDecompressedSize
#include <zstd.h>

struct DBFiles {
    enum Files {
        DATA              = (1ull << 0),
        DATA_INDEX        = (1ull << 1),
        DATA_DBTYPE       = (1ull << 2),
        HEADER            = (1ull << 3),
        HEADER_INDEX      = (1ull << 4),
        HEADER_DBTYPE     = (1ull << 5),
        LOOKUP            = (1ull << 6),
        SOURCE            = (1ull << 7),
        TAX_MAPPING       = (1ull << 8),
        TAX_NAMES         = (1ull << 9),
        TAX_NODES         = (1ull << 10),
        TAX_MERGED        = (1ull << 11),
        CA3M_DATA         = (1ull << 12),
        CA3M_INDEX        = (1ull << 13),
        CA3M_SEQ          = (1ull << 14),
        CA3M_SEQ_IDX      = (1ull << 15),
        CA3M_HDR          = (1ull << 16),
        CA3M_HDR_IDX      = (1ull << 17),
        TAX_BINARY        = (1ull << 18),


        GENERIC           = DATA | DATA_INDEX | DATA_DBTYPE,
        HEADERS           = HEADER | HEADER_INDEX | HEADER_DBTYPE,
        TAXONOMY          = TAX_MAPPING | TAX_NAMES | TAX_NODES | TAX_MERGED | TAX_BINARY,
        SEQUENCE_DB       = GENERIC | HEADERS | TAXONOMY | LOOKUP | SOURCE,
        SEQUENCE_ANCILLARY= SEQUENCE_DB & (~GENERIC),
        SEQUENCE_NO_DATA_INDEX = SEQUENCE_DB & (~DATA_INDEX),

        ALL               = (size_t) -1,
    };
};

template <typename T>
class DBReader : public MemoryTracker {
public:
    struct Index {
        T id;
        size_t offset;
        unsigned int length;

        // we need a non-strict-weak ordering function here
        // so our upper_bound call works correctly
        static bool compareByIdOnly(const Index &x, const Index &y) {
            return x.id <= y.id;
        }

        static bool compareById(const Index &x, const Index &y) {
            if (x.id < y.id)
                return true;
            if (y.id < x.id)
                return false;
            if (x.offset < y.offset)
                return true;
            if (y.offset < x.offset)
                return false;
            if (x.length < y.length)
                return true;
            if (y.length < x.length)
                return false;
            return false;
        }

        static bool compareByOffset(const Index &x, const Index &y) {
            if (x.offset < y.offset)
                return true;
            if (y.offset < x.offset)
                return false;
            if (x.id < y.id)
                return true;
            if (y.id < x.id)
                return false;
            if (x.length < y.length)
                return true;
            if (y.length < x.length)
                return false;
            return false;
        }

        // strict-weak ordering by length, then offset, then id
        static bool compareByLength(const Index &x, const Index &y) {
            if (x.length < y.length)
                return true;
            if (y.length < x.length)
                return false;
            if (x.offset < y.offset)
                return true;
            if (y.offset < x.offset)
                return false;
            return x.id < y.id;
        }
    };

    struct LookupEntry {
        T id;
        std::string entryName;
        unsigned int fileNumber;

        // we need a non-strict-weak ordering function here
        // so our upper_bound call works correctly
        static bool compareByIdOnly(const LookupEntry& x, const LookupEntry& y) {
            return x.id <= y.id;
        }

        static bool compareById(const LookupEntry& x, const LookupEntry& y) {
            if (x.id < y.id)
                return true;
            if (y.id < x.id)
                return false;
            if (x.entryName < y.entryName)
                return true;
            if (y.entryName < x.entryName)
                return false;
            if (x.fileNumber < y.fileNumber)
                return true;
            if (y.fileNumber < x.fileNumber)
                return false;
            return false;
        }

        static bool compareByAccessionOnly(const LookupEntry& x, const LookupEntry& y){
            return x.entryName.compare(y.entryName) <= 0;
        }

        static bool compareByAccession(const LookupEntry& x, const LookupEntry& y) {
            if (x.entryName < y.entryName)
                return true;
            if (y.entryName < x.entryName)
                return false;
            if (x.id < y.id)
                return true;
            if (y.id < x.id)
                return false;
            if (x.fileNumber < y.fileNumber)
                return true;
            if (y.fileNumber < x.fileNumber)
                return false;
            return false;
        }
    };

    // = USE_DATA|USE_INDEX
    DBReader(const char* dataFileName, const char* indexFileName, int threads, int mode);

    DBReader(Index* index, size_t size, size_t aaDbSize, T lastKey,
             int dbType, unsigned int maxSeqLen, int threads);

    void setDataFile(const char* dataFileName);

    virtual ~DBReader();

    bool open(int sort);

    void close();

    const char* getDataFileName() { return dataFileName; }

    const char* getIndexFileName() { return indexFileName; }

    size_t getAminoAcidDBSize();

    size_t getDataSize() { return dataSize; }
    void setDataSize(size_t newSize) { dataSize = newSize; }

    char* getData(size_t id, int thrIdx);

    char* getDataCompressed(size_t id, int thrIdx);

    char* getUnpadded(size_t id, int thrIdx);

    char* getDataUncompressed(size_t id);

    void touchData(size_t id);

    char* getDataByDBKey(T key, int thrIdx);

    char * getDataByOffset(size_t offset);

    size_t getSize() const;

    unsigned int getMaxSeqLen(){ 
            return (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE ) ) ?
                    (std::max(maxSeqLen, 1u)) / Sequence::PROFILE_READIN_SIZE :
                    (std::max(maxSeqLen, 2u));
    }

    T getDbKey(size_t id);


    size_t getSeqLen(size_t id){
        if (id >= size){
            Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
            Debug(Debug::ERROR) << "getSeqLen: local id (" << id << ") >= db size (" << size << ")\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int length;
        if (local2id != NULL) {
            length=index[local2id[id]].length;
        }else{
            length=index[id].length;
        }

        if(Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE ) ){
            // -1 null byte
            return (std::max(length, 1u) - 1u) / Sequence::PROFILE_READIN_SIZE;
        }else{
            // -2 newline and null byte
            return (std::max(length, 2u) - 2u);
        }
    }

    size_t getEntryLen(size_t id){
        if (id >= size){
            Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
            Debug(Debug::ERROR) << "getEntryLen: local id (" << id << ") >= db size (" << size << ")\n";
            EXIT(EXIT_FAILURE);
        }
        if (local2id != NULL) {
            return index[local2id[id]].length;
        }else{
            return index[id].length;
        }
    }

    size_t maxCount(char c);

    void remapData();

    size_t bsearch(const Index * index, size_t size, T value);

    // does a binary search in the index and returns index of the entry with dbKey
    // returns UINT_MAX if the key is not contained in index
    size_t getId (T dbKey);

    // does a binary search in the lookup and returns index of the entry
    size_t getLookupSize() const;
    size_t getLookupIdByKey(T dbKey);
    size_t getLookupIdByAccession(const std::string& accession);
    T getLookupKey(size_t id);
    std::string getLookupEntryName(size_t id);
    unsigned int getLookupFileNumber(size_t id);
    LookupEntry* getLookup() { return lookup; };

    static const int NOSORT = 0;
    static const int SORT_BY_LENGTH = 1;
    static const int LINEAR_ACCCESS = 2;
    static const int SORT_BY_ID     = 3;
    static const int SORT_BY_LINE   = 4; // the local IDs correspond to the line number in the original index file
    static const int SHUFFLE        = 5;
    static const int HARDNOSORT = 6; // do not even sort by ids.
    static const int SORT_BY_ID_OFFSET = 7;
    static const int SORT_BY_OFFSET = 8; // only offset sorting saves memory and does not support random access
    static const int SORT_BY_WEIGHTS= 9;


    static const unsigned int USE_INDEX      = 0;
    static const unsigned int USE_DATA       = 1;
    static const unsigned int USE_WRITABLE   = 2;
    static const unsigned int USE_FREAD      = 4;
    static const unsigned int USE_LOOKUP     = 8;
    static const unsigned int USE_LOOKUP_REV = 16;


    // compressed
    static const int UNCOMPRESSED    = 0;
    static const int COMPRESSED     = 1;

    char * getDataForFile(size_t fileIdx){
        return dataFiles[fileIdx];
    }

    size_t getDataFileCnt(){
        return dataFileCnt;
    }

    size_t getDataSizeForFile(size_t fileIdx){
        return dataSizeOffset[fileIdx+1]-dataSizeOffset[fileIdx];
    }

    std::vector<std::string> getDataFileNames(){
        return dataFileNames;
    }

    size_t getTotalDataSize(){
        return totalDataSize;
    }

    static void moveDatafiles(const std::vector<std::string>& files, const std::string& destination);

    static void moveDb(const std::string &srcDbName, const std::string &dstDbName);

    static void removeDb(const std::string &databaseName);

    static void lookupEntryToBuffer(std::string& buffer, const LookupEntry& entry);

    static void aliasDb(const std::string &databaseName, const std::string &alias, DBFiles::Files dbFilesFlags = DBFiles::ALL);
    static void softlinkDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags = DBFiles::ALL);
    static void copyDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags = DBFiles::ALL);

    char *mmapData(FILE *file, size_t *dataSize);

    bool readIndex(char *data, size_t indexDataSize, Index *index, size_t & dataSize);

    void readLookup(char *data, size_t dataSize, LookupEntry *lookup);

    void readIndexId(T* id, char * line, const char** cols);

    unsigned int indexIdToNum(T* id);

    void readMmapedDataInMemory();

    void mlock();

    void sortIndex(bool isSortedById);
    void sortIndex(float *weights);
    bool isSortedByOffset();

    void unmapData();

    size_t getDataOffset(T i);


    Index* getIndex() {
        return index;
    }

    Index* getIndex(size_t id) {
        if (local2id != NULL) {
            return index + local2id[id];
        }
        return index + id;
    }


    void printMagicNumber();

    T getLastKey();

    static size_t indexMemorySize(const DBReader<unsigned int> &idx);

    static char* serialize(const DBReader<unsigned int> &idx);

    static DBReader<unsigned int> *unserialize(const char* data, int threads);

    int getDbtype() const {
        return dbtype;
    }

    static inline uint16_t getExtendedDbtype(int dbtype) {
        // remove first (compressed) and last bit (compatbility for compressed)
        return (uint16_t)((uint32_t)dbtype >> 16) & 0x7FFE;
    }

    static inline int setExtendedDbtype(int dbtype, uint16_t extended) {
        return dbtype | ((extended & 0x7FFE) << 16);
    }

    const char* getDbTypeName() const {
        return Parameters::getDbTypeName(dbtype);
    }


    struct sortIndecesById {
        sortIndecesById(const Index * ind) : _ind(ind) {}
        bool operator() (unsigned int i, unsigned int j) const { 
            return (_ind[i].id < _ind[j].id); 
        }
        const Index * _ind;
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

    struct comparePairByWeight {
        bool operator() (const std::pair<unsigned int, float>& lhs, const std::pair<unsigned int, float>& rhs) const{
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

    void setSequentialAdvice();

    void decomposeDomainByAminoAcid(size_t worldRank, size_t worldSize, size_t *startEntry, size_t *numEntries);

private:
    void checkClosed() const;

    int threads;

    int dataMode;

    char* dataFileName;
    char* indexFileName;

    // number of entries in the index
    size_t size;

    // offset for all datafiles
    char** dataFiles;
    size_t * dataSizeOffset;
    size_t dataFileCnt;
    size_t totalDataSize;
    std::vector<std::string> dataFileNames;


    // summed up size of all entries
    size_t dataSize;
    // Last Key in Index
    T lastKey;
    // max seqLen
    unsigned int maxSeqLen;
    // flag to check if db was closed
    int closed;
    // stores the dbtype (if dbtype file exists)
    int dbtype;
    int compression;
    int padded;
    char ** compressedBuffers;
    size_t * compressedBufferSizes;
    ZSTD_DStream ** dstream;

    Index * index;
    size_t lookupSize;
    LookupEntry * lookup;
    bool sortedByOffset;

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
