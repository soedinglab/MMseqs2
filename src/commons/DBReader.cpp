#include "DBReader.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>

#include <sys/mman.h>
#include <sys/stat.h>

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
    Debug(Debug::WARNING) << "Magic number " << bytes << "\n";
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
        std::sort(sortArray, sortArray + size, compareIndexLengthPairById());
        for (size_t i = 0; i < size; ++i) {
            index[i].id = sortArray[i].first.id;
            index[i].data = sortArray[i].first.data;
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
    std::pair<Index, unsigned int> *sortArray = new std::pair<Index, unsigned int>[size];
    for (size_t i = 0; i < size; i++) {
        sortArray[i] = std::make_pair(index[i], seqLens[i]);
    }

    std::sort(sortArray, sortArray + size, compareIndexLengthPairById());
    for (size_t i = 0; i < size; ++i) {
        index[i].id = sortArray[i].first.id;
        index[i].data = sortArray[i].first.data;
        seqLens[i] = sortArray[i].second;
    }
    delete[] sortArray;
    if (accessType == SORT_BY_LENGTH) {
        // sort the enties by the length of the sequences
        std::pair<unsigned int, unsigned int> *sortForMapping = new std::pair<unsigned int, unsigned int>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, seqLens[i]);
        }
        std::sort(sortForMapping, sortForMapping + size, comparePairBySeqLength());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
            seqLens[i] = sortForMapping[i].second;
        }
        delete[] sortForMapping;
    } else if (accessType == LINEAR_ACCCESS) {
        // sort the enties by the length of the sequences
        std::pair<unsigned int, size_t> *sortForMapping = new std::pair<unsigned int, size_t>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, (size_t) index[i].data);
        }
        std::sort(sortForMapping, sortForMapping + size, comparePairByOffset());
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
    }
}

template <typename T> char* DBReader<T>::mmapData(FILE * file, size_t *dataSize){
    struct stat sb;
    fstat(fileno(file), &sb);
    *dataSize = sb.st_size;
    int fd =  fileno(file);
    int mode = PROT_READ;
    if(dataMode & USE_WRITABLE) {
        mode |= PROT_WRITE;
    }
    return (char*)mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0);
}

template <typename T> void DBReader<T>::remapData(){
    if(dataMode & USE_DATA){
        munmap(data, dataSize);
        data = mmapData(dataFile, &dataSize);
    }
}

template <typename T> void DBReader<T>::close(){
    if(dataMode & USE_DATA){
        fclose(dataFile);
        unmapData();
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS){
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

    if ((size_t) (index[id].data - data) >= dataSize){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: global id (" << id << ")\n";
        Debug(Debug::ERROR) << "Size of data: " << dataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << (size_t)  (index[id].data - data) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS){
        return index[local2id[id]].data;
    }else{
        return index[id].data;
    }
}


template <typename T> char* DBReader<T>::getDataByDBKey(T dbKey){
    checkClosed();
    if(!(dataMode & USE_DATA)) {
        Debug(Debug::ERROR) << "DBReader is just open in INDEX_ONLY mode. Call of getData is not allowed" << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t id = bsearch(index, size, dbKey);
    return (index[id].id == dbKey) ? index[id].data : NULL;
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
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS){
        id = local2id[id];
    }
    return index[id].id;
}

template <typename T> size_t DBReader<T>::getId (T dbKey){
    size_t id = bsearch(index, size, dbKey);
    if(accessType == SORT_BY_LENGTH || accessType == LINEAR_ACCCESS){
        return  (index[id].id == dbKey) ? id2local[id] : UINT_MAX;
    }
    return (index[id].id == dbKey) ? id : UINT_MAX;
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
    std::ifstream indexFile(indexFileName);

    if (indexFile.fail()) {
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }

    char *save;
    size_t i = 0;
    std::string line;
    while (std::getline(indexFile, line)) {
        char *l = (char *) line.c_str();
        readIndexId(&index[i].id, l, &save);
        size_t offset = strtoull(strtok_r(NULL, "\t", &save), NULL, 10);
        size_t length = strtoull(strtok_r(NULL, "\t", &save), NULL, 10);

        if (i >= this->size) {
            Debug(Debug::ERROR) << "Corrupt memory, too many entries!\n";
            EXIT(EXIT_FAILURE);
        }

        if (dataMode & USE_DATA) {
            index[i].data = data + offset;
        } else {
            index[i].data = (char *) offset;
        }

        entryLength[i] = length;

        i++;
    }

    indexFile.close();
}

template<>
void DBReader<std::string>::readIndexId(std::string* id, char* line, char** save) {
    *id = strtok_r(line, "\t", save);
}
template<>
void DBReader<unsigned int>::readIndexId(unsigned int* id, char* line, char** save) {
    *id = (unsigned int) strtoul(strtok_r(line, "\t", save), NULL, 10);
}

template <typename T> void DBReader<T>::unmapData() {
    if(dataMapped == true){
        munmap(data, dataSize);
        dataMapped = false;
    }
}

template <typename T>  size_t DBReader<T>::getDataOffset(T i) {
    size_t id = bsearch(index, size, i);
    return index[id].data - data;
}

template class DBReader<unsigned int>;
template class DBReader<std::string>;
