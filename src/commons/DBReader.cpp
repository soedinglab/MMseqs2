#include "DBReader.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include <cstddef>
#include <random>

#include <sys/mman.h>
#include <sys/stat.h>
#include <omptl/omptl_algorithm>


#include "MemoryMapped.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

template <typename T>
DBReader<T>::DBReader(const char* dataFileName_, const char* indexFileName_, int threads, int dataMode) :
        data(NULL), threads(threads), dataMode(dataMode), dataFileName(strdup(dataFileName_)),
        indexFileName(strdup(indexFileName_)), size(0), dataSize(0), aaDbSize(0), lastKey(T()), closed(1), dbtype(-1),
        compressedBuffers(NULL), compressedBufferSizes(NULL), index(NULL), seqLens(NULL), id2local(NULL), local2id(NULL),
        dataMapped(false), accessType(0), externalData(false), didMlock(false)
{}

template <typename T>
DBReader<T>::DBReader(DBReader<T>::Index *index, unsigned int *seqLens, size_t size, size_t aaDbSize, T lastKey,
        int dbType, unsigned int maxSeqLen, int threads) :
        data(NULL), threads(threads), dataMode(USE_INDEX), dataFileName(NULL), indexFileName(NULL),
        size(size), dataSize(0), aaDbSize(aaDbSize), lastKey(lastKey), maxSeqLen(maxSeqLen), closed(1), dbtype(dbType),
        compressedBuffers(NULL), compressedBufferSizes(NULL), index(index), seqLens(seqLens), id2local(NULL), local2id(NULL),
        dataMapped(false), accessType(NOSORT), externalData(true), didMlock(false)
{}

template <typename T>
void DBReader<T>::setDataFile(const char* dataFileName_)  {
    if (dataFileName != NULL) {
        unmapData();
        free(dataFileName);
    }

    dataMode |= USE_DATA;
    dataFileName = strdup(dataFileName_);
}


template <typename T>
void DBReader<T>::readMmapedDataInMemory(){
    if ((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
        Debug(Debug::INFO) << "Touch data file " << dataFileName << " ... ";
        magicBytes = Util::touchMemory(data, dataSize);
        Debug(Debug::INFO) << "Done.\n";
    }
}

template <typename T>
void DBReader<T>::mlock(){
    if (dataMode & USE_DATA) {
        if (didMlock == false) {
            ::mlock(data, dataSize);
            didMlock = true;
        }
    }
}


template <typename T>
void DBReader<T>::printMagicNumber(){
    Debug(Debug::INFO) << magicBytes << "\n";
}

template <typename T> DBReader<T>::~DBReader(){
    if(dataFileName != NULL) {
        free(dataFileName);
    }

    if(indexFileName != NULL) {
        free(indexFileName);
    }
}



template <typename T> bool DBReader<T>::open(int accessType){
    // count the number of entries
    this->accessType = accessType;
    if (dataFileName != NULL) {
        dbtype = parseDbType(dataFileName);
    }

    if (dataMode & USE_DATA) {
        FILE* dataFile = fopen(dataFileName, "r");
        if (dataFile == NULL) {
            Debug(Debug::ERROR) << "Could not open data file " << dataFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        data = mmapData(dataFile, &dataSize);
        fclose(dataFile);
        dataMapped = true;
    }
    bool isSortedById = false;
    if (externalData == false) {
        if(FileUtil::fileExists(indexFileName)==false){
            Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        size = FileUtil::countLines(indexFileName);
        index = new Index[this->size];
        seqLens = new unsigned int[size];

        isSortedById = readIndex(indexFileName, index, seqLens);

        if (accessType != HARDNOSORT) {
            sortIndex(isSortedById);
        }

        // init seq lens array and dbKey mapping
        aaDbSize = 0;
        for (size_t i = 0; i < size; i++){
            unsigned int size = seqLens[i];
            aaDbSize += size;
        }
    }

    compression = isCompressed(dbtype);
    if(compression == COMPRESSED){
        compressedBufferSizes = new size_t[threads];
        compressedBuffers = new char*[threads];
        dstream = new ZSTD_DStream*[threads];
        for(int i = 0; i < threads; i++){
            // allocated buffer
            compressedBufferSizes[i] = std::max(maxSeqLen+1, 1024u);
            compressedBuffers[i] = (char*) malloc(compressedBufferSizes[i]);
            if(compressedBuffers[i]==NULL){
                Debug(Debug::ERROR) << "Could not allocate compressedBuffer!\n";
                EXIT(EXIT_FAILURE);
            }
            dstream[i] = ZSTD_createDStream();
            if (dstream==NULL) {
                Debug(Debug::ERROR) << "ZSTD_createDStream() error \n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    closed = 0;
    return isSortedById;
}

template<typename T>
void DBReader<T>::sortIndex(bool) {
}

template<>
void DBReader<std::string>::sortIndex(bool isSortedById) {
    if (accessType == SORT_BY_ID){
        if (isSortedById) {
            return;
        }

        std::pair<Index, unsigned int> *sortArray = new std::pair<Index, unsigned int>[size];
        for (size_t i = 0; i < size; i++) {
            sortArray[i] = std::make_pair(index[i], seqLens[i]);
        }
        omptl::sort(sortArray, sortArray + size, compareIndexLengthPairById());
        for (size_t i = 0; i < size; ++i) {
            index[i].id = sortArray[i].first.id;
            index[i].offset = sortArray[i].first.offset;
            seqLens[i] = sortArray[i].second;
        }
        delete[] sortArray;
    }else{
        if(accessType != NOSORT && accessType != HARDNOSORT){
            Debug(Debug::ERROR) << "DBReader<std::string> can not be opened in sort mode\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

template<>
void DBReader<unsigned int>::sortIndex(bool isSortedById) {
    // First, we sort the index by IDs and we keep track of the original
    // ordering in mappingToOriginalIndex array
    size_t* mappingToOriginalIndex=NULL;
    if(accessType==SORT_BY_LINE){
        mappingToOriginalIndex = new size_t[size];
    }
    if(isSortedById == false){
        std::pair<Index, std::pair<size_t,unsigned int> > *sortArray = new std::pair<Index, std::pair<size_t,unsigned int> >[size];
        for (size_t i = 0; i < size; i++) {
            sortArray[i] = std::make_pair(index[i], std::make_pair(i,seqLens[i]));
        }
        omptl::sort(sortArray, sortArray + size, compareIndexLengthPairByIdKeepTrack());
        for (size_t i = 0; i < size; ++i) {
            index[i].id = sortArray[i].first.id;
            index[i].offset = sortArray[i].first.offset;
            seqLens[i] = (sortArray[i].second).second;
        }
        if(accessType==SORT_BY_LINE){
            for (size_t i = 0; i < size; ++i) {
                mappingToOriginalIndex[i] = (sortArray[i].second).first;
            }
        }
        delete[] sortArray;
    }else {
        if(accessType== SORT_BY_LINE){
            for (size_t i = 0; i < size; ++i) {
                mappingToOriginalIndex[i] = i;
            }
        }
    }
    if (accessType == SORT_BY_LENGTH) {
        // sort the entries by the length of the sequences
        std::pair<unsigned int, unsigned int> *sortForMapping = new std::pair<unsigned int, unsigned int>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, seqLens[i]);
        }
        //this sort has to be stable to assure same clustering results
        omptl::sort(sortForMapping, sortForMapping + size, comparePairBySeqLength());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
            seqLens[i] = sortForMapping[i].second;
        }
        delete[] sortForMapping;
    } else if (accessType == SHUFFLE) {
        size_t *tmpIndex = new size_t[size];
        for (size_t i = 0; i < size; i++) {
            tmpIndex[i] = i;
        }

        std::mt19937 rnd(0);
        std::shuffle(tmpIndex, tmpIndex + size, rnd);

        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[tmpIndex[i]] = i;
            local2id[i] = tmpIndex[i];
        }
        delete[] tmpIndex;

        unsigned int *tmpSize = new unsigned int[size];
        memcpy(tmpSize, seqLens, size * sizeof(unsigned int));
        for (size_t i = 0; i < size; i++) {
            seqLens[i] = tmpSize[local2id[i]];
        }
        delete[] tmpSize;
    } else if (accessType == LINEAR_ACCCESS) {
        // sort the entries by the offset of the sequences
        std::pair<unsigned int, size_t> *sortForMapping = new std::pair<unsigned int, size_t>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, index[i].offset);
        }
        omptl::sort(sortForMapping, sortForMapping + size, comparePairByOffset());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
        }
        delete[] sortForMapping;
        unsigned int *tmpSizeArray = new unsigned int[size];
        memcpy(tmpSizeArray, seqLens, size * sizeof(unsigned int));
        for (size_t i = 0; i < size; i++) {
            seqLens[i] = tmpSizeArray[local2id[i]];
        }
        delete[] tmpSizeArray;
    } else if (accessType == SORT_BY_ID_OFFSET) {
        // sort the entries by the offset of the sequences
        std::pair<unsigned int, Index> *sortForMapping = new std::pair<unsigned int, Index>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, index[i]);
        }
        omptl::sort(sortForMapping, sortForMapping + size, comparePairByIdAndOffset());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
        }
        delete[] sortForMapping;
        unsigned int *tmpSizeArray = new unsigned int[size];
        memcpy(tmpSizeArray, seqLens, size * sizeof(unsigned int));
        for (size_t i = 0; i < size; i++) {
            seqLens[i] = tmpSizeArray[local2id[i]];
        }
        delete[] tmpSizeArray;
    } else if (accessType == SORT_BY_LINE) {
        // sort the entries by the original line number in the index file
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = mappingToOriginalIndex[i];
            local2id[mappingToOriginalIndex[i]] = i;
        }
        unsigned int *tmpSizeArray = new unsigned int[size];
        memcpy(tmpSizeArray, seqLens, size * sizeof(unsigned int));
        for (size_t i = 0; i < size; i++) {
            seqLens[i] = tmpSizeArray[local2id[i]];
        }
        delete[] tmpSizeArray;
    }
    if(mappingToOriginalIndex){
        delete [] mappingToOriginalIndex;
    }
}

template <typename T> char* DBReader<T>::mmapData(FILE * file, size_t *dataSize) {
    struct stat sb;
    if (fstat(fileno(file), &sb) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << dataFileName << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }

    *dataSize = sb.st_size;
    int fd =  fileno(file);

    char *ret;
    if ((dataMode & USE_FREAD) == 0) {
        int mode;
        if (dataMode & USE_WRITABLE) {
            mode = PROT_READ | PROT_WRITE;
        } else {
            mode = PROT_READ;
        }
        ret = static_cast<char*>(mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0));
        if (ret == MAP_FAILED){
            int errsv = errno;
            Debug(Debug::ERROR) << "Failed to mmap memory dataSize=" << *dataSize <<" File=" << dataFileName << ". Error " << errsv << ".\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
        ret = static_cast<char*>(malloc(*dataSize));
        Util::checkAllocation(ret, "Not enough system memory to read in the whole data file.");
        size_t result = fread(ret, 1, *dataSize, file);
        if (result != *dataSize) {
            Debug(Debug::ERROR) << "Failed to read in datafile (" << dataFileName << "). Error " << errno << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return ret;
}

template <typename T> void DBReader<T>::remapData(){
    if ((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
        unmapData();
        FILE* dataFile = fopen(dataFileName, "r");
        data = mmapData(dataFile, &dataSize);
        fclose(dataFile);
        dataMapped = true;
    }
}

template <typename T> void DBReader<T>::close(){
    if(dataMode & USE_DATA){
        unmapData();
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE || accessType == SHUFFLE){
        delete [] id2local;
        delete [] local2id;
    }

    if(compressedBuffers){
        for(int i = 0; i < threads; i++){
            ZSTD_freeDStream(dstream[i]);
            delete [] compressedBuffers[i];
        }
        delete [] compressedBuffers;
        delete [] compressedBufferSizes;
        delete [] dstream;
    }

    if(externalData == false) {
        delete[] index;
        delete[] seqLens;
    }
    closed = 1;
}

template <typename T> size_t DBReader<T>::bsearch(const Index * index, size_t N, T value)
{
    Index val;
    val.id = value;
    return std::upper_bound(index, index + N, val, Index::compareById) - index;
}

template <typename T> char* DBReader<T>::getDataCompressed(size_t id, int thrIdx) {
    char *data = getDataUncompressed(id);

    unsigned int cSize = *(reinterpret_cast<unsigned int *>(data));

    size_t totalSize = 0;
    const void *cBuff = static_cast<void *>(data + sizeof(unsigned int));
    const char *dataStart = data + sizeof(unsigned int);
    bool isCompressed = (dataStart[cSize] == 0) ? true : false;
    if(isCompressed){
        ZSTD_inBuffer input = {cBuff, cSize, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {compressedBuffers[thrIdx], compressedBufferSizes[thrIdx], 0};
            // size of next compressed block
            size_t toRead = ZSTD_decompressStream(dstream[thrIdx], &output, &input);
            if (ZSTD_isError(toRead)) {
                Debug(Debug::ERROR) << "ERROR: " << id << " ZSTD_decompressStream " << ZSTD_getErrorName(toRead) << "\n";
                EXIT(EXIT_FAILURE);
            }
            totalSize += output.pos;
        }
        compressedBuffers[thrIdx][totalSize] = '\0';
    }else{
        memcpy(compressedBuffers[thrIdx], cBuff, cSize);
        compressedBuffers[thrIdx][cSize] = '\0';
    }
    return compressedBuffers[thrIdx];
}

template <typename T> char* DBReader<T>::getData(size_t id, int thrIdx){
    if(compression == COMPRESSED){
        return getDataCompressed(id, thrIdx);
    }else{
        return getDataUncompressed(id);
    }
}

template <typename T> char* DBReader<T>::getDataUncompressed(size_t id){
    checkClosed();
    if(!(dataMode & USE_DATA)) {
        Debug(Debug::ERROR) << "DBReader is just open in INDEXONLY mode. Call of getData is not allowed" << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }

    if ((size_t) (index[id].offset) >= dataSize){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: global id (" << id << ")\n";
        Debug(Debug::ERROR) << "Size of data: " << dataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << index[id].offset << "\n";
        EXIT(EXIT_FAILURE);
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE || accessType == SHUFFLE){
        return data + index[local2id[id]].offset;
    }else{
        return data + index[id].offset;
    }
}

template <typename T>
void DBReader<T>::touchData(size_t id) {
    if((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
        char *data = getDataUncompressed(id);
        size_t size = getSeqLens(id);
        magicBytes = Util::touchMemory(data, size);
    }
}

template <typename T> char* DBReader<T>::getDataByDBKey(T dbKey, int thrIdx) {
    size_t id = getId(dbKey);
    if(compression == COMPRESSED ){
        return (id != UINT_MAX) ? getDataCompressed(id, thrIdx) : NULL;
    }else{
        return (id != UINT_MAX) ? data + index[id].offset : NULL;
    }
}

template <typename T> size_t DBReader<T>::getSize (){
    checkClosed();
    return size;
}

template <typename T> T DBReader<T>::getDbKey (size_t id){
    checkClosed();
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE || accessType == SHUFFLE){
        id = local2id[id];
    }
    return index[id].id;
}

template <typename T> size_t DBReader<T>::getId (T dbKey){
    size_t id = bsearch(index, size, dbKey);
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE || accessType == SHUFFLE){
        return  (id < size && index[id].id == dbKey) ? id2local[id] : UINT_MAX;
    }
    return (id < size && index[id].id == dbKey ) ? id : UINT_MAX;
}

template <typename T> unsigned int* DBReader<T>::getSeqLens(){
    return seqLens;
}

template <typename T> size_t DBReader<T>::getSeqLens(size_t id){
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getSeqLens: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return seqLens[id];
}

template <typename T> size_t DBReader<T>::maxCount(char c) {
    checkClosed();

    size_t max = 0;
    if (compression == COMPRESSED) {
        size_t entries = getSize();
#ifdef OPENMP
        size_t localThreads = std::min(entries, static_cast<size_t>(threads));
#endif
#pragma omp parallel num_threads(localThreads)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10) reduction(max:max)
            for (size_t id = 0; id < entries; id++) {
                char *data = getData(id, thread_idx);
                size_t count = 0;
                for (size_t i = 0; i < seqLens[id]; ++i) {
                    if (data[i] == c) {
                        count++;
                    }
                }
                max = std::max(max, count);
            }
        }
        return max;
    }

    size_t count = 0;
    for (size_t i = 0; i < dataSize; ++i) {
        if (data[i] == c) {
            count++;
        }

        if (data[i] == '\0') {
            max = std::max(max, count);
            count = 0;
        }
    }

    return max;
}

template <typename T> void DBReader<T>::checkClosed(){
    if (closed == 1){
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

template<typename T>
bool DBReader<T>::readIndex(char *indexFileName, Index *index, unsigned int *entryLength) {
    MemoryMapped indexData(indexFileName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!indexData.isValid()){
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t i = 0;

    size_t currPos = 0;
    char* indexDataChar = (char *) indexData.getData();
    char * cols[3];
    T prevId=T(); // makes 0 or empty string
    size_t isSorted = true;
    maxSeqLen=0;
    while (currPos < indexData.size()){
        if (i >= this->size) {
            Debug(Debug::ERROR) << "Corrupt memory, too many entries!\n";
            EXIT(EXIT_FAILURE);
        }
        Util::getWordsOfLine(indexDataChar, cols, 3 );
        readIndexId(&index[i].id, indexDataChar, cols);
        isSorted *= (index[i].id >= prevId);
        size_t offset = Util::fast_atoi<size_t>(cols[1]);
        size_t length = Util::fast_atoi<size_t>(cols[2]);
        index[i].offset = offset;
        entryLength[i] = length;
        maxSeqLen = std::max(static_cast<unsigned int>(length), maxSeqLen);
        indexDataChar = Util::skipLine(indexDataChar);
        currPos = indexDataChar - (char *) indexData.getData();
        lastKey = std::max(index[i].id, lastKey);
        prevId = index[i].id;
        i++;
    }
    indexData.close();
    return isSorted;
}

template<typename T> T DBReader<T>::getLastKey() {
    return lastKey;
}

template<>
void DBReader<std::string>::readIndexId(std::string* id, char* line, char** cols){
    ptrdiff_t keySize =  ((cols[1] - 1) - line) ;
    id->assign(line, keySize);
}
template<>
void DBReader<unsigned int>::readIndexId(unsigned int* id, char*, char** cols) {
    *id = Util::fast_atoi<unsigned int>(cols[0]);
}

template <typename T> void DBReader<T>::unmapData() {
    if(dataMapped == true){
        if (didMlock == true) {
            munlock(data, dataSize);
            didMlock = false;
        }

        if ((dataMode & USE_FREAD) == 0) {
            if(munmap(data, dataSize) < 0){
                Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSize <<" File=" << dataFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            free(data);
        }
        dataMapped = false;
    }
}

template <typename T>  size_t DBReader<T>::getDataOffset(T i) {
    size_t id = bsearch(index, size, i);
    return index[id].offset;
}

template <>
size_t DBReader<unsigned int>::indexMemorySize(const DBReader<unsigned int> &idx) {
    size_t memSize = // size + aaDbSize
            2 * sizeof(size_t)
            // maxSeqLen + lastKey + dbtype
            + 3 * sizeof(unsigned int)
            // index
            + idx.size * sizeof(DBReader<unsigned int>::Index)
            // seqLens
            + idx.size * sizeof(unsigned int);

    return memSize;
}

template <>
char* DBReader<unsigned int>::serialize(const DBReader<unsigned int> &idx) {
    char* data = (char*) malloc(indexMemorySize(idx));
    char* p = data;
    memcpy(p, &idx.size, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.aaDbSize, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.lastKey, sizeof(unsigned int));
    p += sizeof(unsigned int);
    memcpy(p, &idx.dbtype, sizeof(int));
    p += sizeof(unsigned int);
    memcpy(p, &idx.maxSeqLen, sizeof(unsigned int));
    p += sizeof(unsigned int);
    memcpy(p, idx.index, idx.size * sizeof(DBReader<unsigned int>::Index));
    p += idx.size * sizeof(DBReader<unsigned int>::Index);
    memcpy(p, idx.seqLens, idx.size * sizeof(unsigned int));

    return data;
}

template <>
DBReader<unsigned int> *DBReader<unsigned int>::unserialize(const char* data, int threads) {
    const char* p = data;
    size_t size = *((size_t*)p);
    p += sizeof(size_t);
    size_t aaDbSize = *((size_t*)p);
    p += sizeof(size_t);
    unsigned int lastKey = *((unsigned int*)p);
    p += sizeof(unsigned int);
    int dbType = *((int*)p);
    p += sizeof(int);
    unsigned int maxSeqLen = *((unsigned int*)p);
    p += sizeof(unsigned int);
    DBReader<unsigned int>::Index *idx = (DBReader<unsigned int>::Index *)p;
    p += size * sizeof(DBReader<unsigned int>::Index);
    unsigned int *seqLens = (unsigned int *)p;

    return new DBReader<unsigned int>(idx, seqLens, size, aaDbSize, lastKey, dbType, maxSeqLen, threads);
}

template <typename T>
int DBReader<T>::parseDbType(const char *name) {
    std::string dbTypeFile = std::string(name) + ".dbtype";
    int dbtype = -1;
    if (FileUtil::fileExists(dbTypeFile.c_str()) == true) {
        size_t fileSize = FileUtil::getFileSize(dbTypeFile);
        if (fileSize != sizeof(int)) {
            Debug(Debug::ERROR) << "File size of " << dbTypeFile << " seems to be wrong!\n";
            Debug(Debug::ERROR) << "It should have 4 bytes but it has " <<  fileSize << " bytes.";
            EXIT(EXIT_FAILURE);
        }
        FILE *dbtypeDataFile = fopen(dbTypeFile.c_str(), "r");
        if (dbtypeDataFile == NULL) {
            Debug(Debug::ERROR) << "Could not open data file " << dbTypeFile << "!\n";
            EXIT(EXIT_FAILURE);
        }
        size_t result = fread(&dbtype, 1, fileSize, dbtypeDataFile);
        if (result != fileSize) {
            Debug(Debug::ERROR) << "Could not read " << dbTypeFile << "!\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(dbtypeDataFile);
    }
    return dbtype;
}

template<typename T>
void DBReader<T>::setData(char *data, size_t dataSize) {
    this->data = data;
    this->dataSize = dataSize;
}

template<typename T>
void DBReader<T>::setMode(const int mode) {
    this->dataMode = mode;
}

template<typename T>
size_t DBReader<T>::getOffset(size_t id) {
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE || accessType == SHUFFLE){
        id = local2id[id];
    }
    return index[id].offset;
}

template<typename T>
size_t DBReader<T>::findNextOffsetid(size_t id) {
    size_t idOffset = getOffset(id);
    size_t nextOffset = SIZE_MAX;
    for(size_t i = 0; i < size; i++){
        if(index[i].offset > idOffset && index[i].offset < nextOffset){
            nextOffset=index[i].offset;
        }
    }
    return nextOffset;
}

template<typename T>
int DBReader<T>::isCompressed(int dbtype) {
    return (dbtype & (1 << 31))? COMPRESSED : UNCOMPRESSED;
}

template class DBReader<unsigned int>;
template class DBReader<std::string>;
