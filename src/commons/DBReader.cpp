#include "DBReader.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include <cstddef>

#include <sys/mman.h>
#include <sys/stat.h>
#include <omptl/omptl_algorithm>

#include "MemoryMapped.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

template <typename T> DBReader<T>::DBReader(const char* dataFileName_, const char* indexFileName_, int dataMode)
{
    dataSize = 0;
    this->dataMode = dataMode;
    this->dataFileName = strdup(dataFileName_);
    this->indexFileName = strdup(indexFileName_);
    closed = 1;
    accessType = 0;
}

template <typename T> void DBReader<T>::readMmapedDataInMemory(){
    size_t bytes = 0;
    for(size_t i = 0; i < dataSize; i++){
        bytes += data[i];
    }
    this->magicBytes = bytes;
}

template <typename T>
void DBReader<T>::printMagicNumber(){
    Debug(Debug::INFO) << magicBytes << "\n";
}

template <typename T> DBReader<T>::~DBReader(){
    free(dataFileName);
    free(indexFileName);
}

template <typename T> void DBReader<T>::open(int accessType){
    // count the number of entries
    this->size = FileUtil::countLines(indexFileName);
    this->accessType = accessType;

    if (dataMode & USE_DATA) {
        dataFile = fopen(dataFileName, "r");
        if (dataFile == NULL) {
            Debug(Debug::ERROR) << "Could not open data file " << dataFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        data = mmapData(dataFile, &dataSize);
        dataMapped = true;
    }

    index = new Index[this->size];
    seqLens = new unsigned int[size];

    readIndex(indexFileName, index, data, seqLens);

    sortIndex();

    // init seq lens array and dbKey mapping
    aaDbSize = 0;
    for (size_t i = 0; i < size; i++){
        unsigned int size = seqLens[i];
        aaDbSize += size;
    }

    closed = 0;
}

template<typename T>
void DBReader<T>::sortIndex() {
}

template<>
void DBReader<std::string>::sortIndex() {
    if (accessType == SORT_BY_ID){
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
        if(accessType != NOSORT){
            Debug(Debug::ERROR) << "DBReader<std::string> can not be opend in sort mode\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

template<>
void DBReader<unsigned int>::sortIndex() {
	
	// First, we sort the index by IDs and we keep track of the original
	// ordering in mappingToOriginalIndex array
	size_t* mappingToOriginalIndex = new size_t[size];
	
    std::pair<Index, std::pair<size_t,unsigned int> > *sortArray = new std::pair<Index, std::pair<size_t,unsigned int> >[size];
    for (size_t i = 0; i < size; i++) {
        sortArray[i] = std::make_pair(index[i], std::make_pair(i,seqLens[i]));
    }
    omptl::sort(sortArray, sortArray + size, compareIndexLengthPairByIdKeepTrack());
    for (size_t i = 0; i < size; ++i) {
        index[i].id = sortArray[i].first.id;
        index[i].offset = sortArray[i].first.offset;
        seqLens[i] = (sortArray[i].second).second;
		mappingToOriginalIndex[i] = (sortArray[i].second).first;
    }
	
    delete[] sortArray;
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
        omptl::sort(sortForMapping, sortForMapping + size, comparePairBySeqLength());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
            seqLens[i] = sortForMapping[i].second;
        }
        delete[] sortForMapping;
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
    delete [] mappingToOriginalIndex;
}

template <typename T> char* DBReader<T>::mmapData(FILE * file, size_t *dataSize){
    struct stat sb;
    fstat(fileno(file), &sb);
    *dataSize = sb.st_size;
    int fd =  fileno(file);
    int mode;
    if(dataMode & USE_WRITABLE) {
        mode = PROT_READ | PROT_WRITE;
    } else {
        mode = PROT_READ;
    }
    void * ret = mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0);
    if(ret == MAP_FAILED){
        Debug(Debug::ERROR) << "Failed to mmap memory dataSize=" << *dataSize <<" File=" << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    return static_cast<char*>(ret);
}

template <typename T> void DBReader<T>::remapData(){
    if(dataMode & USE_DATA){
        if(munmap(data, dataSize) < 0){
            Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSize <<" File=" << dataFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        data = mmapData(dataFile, &dataSize);
    }
}

template <typename T> void DBReader<T>::close(){
    if(dataMode & USE_DATA){
        fclose(dataFile);
        unmapData();
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE){
        delete [] id2local;
        delete [] local2id;
    }
    delete [] index;
    delete [] seqLens;
    closed = 1;
}

template <typename T> size_t DBReader<T>::bsearch(const Index * index, size_t N, T value)
{
    Index val;
    val.id = value;
    return std::upper_bound(index, index + N, val, Index::compareById) - index;
}

template <typename T> char* DBReader<T>::getData(size_t id){
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
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE){
        return data + index[local2id[id]].offset;
    }else{
        return data + index[id].offset;
    }
}

template <typename T> char* DBReader<T>::getDataByDBKey(T dbKey) {
    size_t id = getId(dbKey);
    return (id != UINT_MAX) ? data + index[id].offset : NULL;
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
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS || accessType == SORT_BY_LINE){
        id = local2id[id];
    }
    return index[id].id;
}

template <typename T> size_t DBReader<T>::getId (T dbKey){
    size_t id = bsearch(index, size, dbKey);
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS  || accessType == SORT_BY_LINE){
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
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return seqLens[id];
}

template <typename T> void DBReader<T>::checkClosed(){
    if (closed == 1){
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

template<typename T>
void DBReader<T>::readIndex(char *indexFileName, Index *index, char *data, unsigned int *entryLength) {
    MemoryMapped indexData(indexFileName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!indexData.isValid()){
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t i = 0;

    size_t currPos = 0;
    char* indexDataChar = (char *) indexData.getData();
    char * cols[3];

    while (currPos < indexData.size()){
        if (i >= this->size) {
            Debug(Debug::ERROR) << "Corrupt memory, too many entries!\n";
            EXIT(EXIT_FAILURE);
        }
        Util::getWordsOfLine(indexDataChar, cols, 3 );
        readIndexId(&index[i].id, indexDataChar, cols);
        size_t offset = strtoull(cols[1], NULL, 10);
        size_t length = strtoull(cols[2], NULL, 10);
        index[i].offset = offset;
        entryLength[i] = length;
        indexDataChar = Util::skipLine(indexDataChar);
        currPos = indexDataChar - (char *) indexData.getData();
        lastKey = std::max(index[i].id,lastKey); 
        i++;

    }
    indexData.close();
}

template<typename T> T DBReader<T>::getLastKey()
{
    return lastKey;
}

template<>
void DBReader<std::string>::readIndexId(std::string* id, char * line, char** cols){
    ptrdiff_t keySize =  ((cols[1] - 1) - line) ;
    id->assign(line, keySize);
}
template<>
void DBReader<unsigned int>::readIndexId(unsigned int* id, char * line, char** cols) {
    *id = (unsigned int) strtoul(cols[0], NULL, 10);
}

template <typename T> void DBReader<T>::unmapData() {
    if(dataMapped == true){
        if(munmap(data, dataSize) < 0){
            Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSize <<" File=" << dataFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        dataMapped = false;
    }
}

template <typename T>  size_t DBReader<T>::getDataOffset(T i) {
    size_t id = bsearch(index, size, i);
    return index[id].offset;
}

template class DBReader<unsigned int>;
template class DBReader<std::string>;
