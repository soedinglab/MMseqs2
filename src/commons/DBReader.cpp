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

template <typename T> DBReader<T>::DBReader(const char* dataFileName_, const char* indexFileName_, int dataMode)
{
    dataSize = 0;
    this->dataMode = dataMode;
    this->dataFileName = new char [strlen(dataFileName_) + 1];
    memcpy(dataFileName, dataFileName_, sizeof(char) * (strlen(dataFileName_) + 1));
    this->indexFileName = new char [strlen(indexFileName_) + 1];
    memcpy(indexFileName, indexFileName_, sizeof(char) * (strlen(indexFileName_) + 1));
    closed = 1;
}

template <typename T> DBReader<T>::~DBReader(){
    delete[] dataFileName;
    delete[] indexFileName;
}

template <typename T> void DBReader<T>::open(int sort){
    // count the number of entries
    this->size = countLine(indexFileName);

    // open ffindex
    if(dataMode == DATA_AND_INDEX){
        dataFile = fopen(dataFileName, "r");
        if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", dataFileName);  EXIT(EXIT_FAILURE); }
        data = mmapData(dataFile, &dataSize);
        dataMapped = true;
    }
    index = new Index[this->size];
    seqLens = new unsigned int[size];
    readIndex(indexFileName, index, data, seqLens);
    std::pair<Index, unsigned  int> * sortArray = new std::pair<Index, unsigned  int>[size];
    for(size_t i = 0; i < size; i++){
        sortArray[i] = std::make_pair( index[i], seqLens[i] );
    }
    std::sort(sortArray, sortArray + size, compareIndexLengthPairById() );

    for(size_t i = 0; i < size; ++i )
    {
        index[i].id = sortArray[i].first.id;
        index[i].data = sortArray[i].first.data;
        seqLens[i] = sortArray[i].second;
    }
    delete [] sortArray;
    // init seq lens array and dbKey mapping
    aaDbSize = 0;
    for (size_t i = 0; i < size; i++){
        unsigned int size = seqLens[i];
        aaDbSize += size;
    }
    if (aaDbSize == 0){
        Debug(Debug::ERROR) << "Invalid database in data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }

    closed = 0;
}

template <typename T> char* DBReader<T>::mmapData(FILE * file, size_t *dataSize){
    struct stat sb;
    fstat(fileno(file), &sb);
    *dataSize = sb.st_size;
    int fd =  fileno(file);
    return (char*)mmap(NULL, *dataSize, PROT_READ, MAP_PRIVATE, fd, 0);
}

template <typename T> void DBReader<T>::remapData(){
    if(dataMode == DATA_AND_INDEX){
        munmap(data, dataSize);
        data = mmapData(dataFile, &dataSize);
    }
}

template <typename T> void DBReader<T>::close(){
    if(dataMode == DATA_AND_INDEX){
        fclose(dataFile);
        unmapData();
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

template <typename T> char* DBReader<T>::getData (size_t id){
    checkClosed();
    if(dataMode == INDEXONLY){
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
    return index[id].data;
}


template <typename T> char* DBReader<T>::getDataByDBKey (T dbKey){
    checkClosed();
    if(dataMode ==INDEXONLY){
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
    return index[id].id;
}

template <typename T> size_t DBReader<T>::getId (T dbKey){
    size_t id = bsearch(index, size, dbKey);
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

template <typename T> size_t DBReader<T>::countLine(const char *indexFileName) {
    size_t cnt = 0;
    std::string line;
    std::ifstream index_file(indexFileName);
    if (index_file.is_open()) {
        while (std::getline(index_file, line)) {
            cnt++;
        }
        index_file.close();
    }
    else {
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    return cnt;
}

template <typename T> void DBReader<T>::readIndex(char *indexFileName, Index * index, char *data, unsigned int * entryLength) {
    std::ifstream index_file(indexFileName);
    std::string line;
    if (index_file.is_open()) {
        size_t i = 0;

        while (std::getline(index_file, line)) {
            std::string id;
            size_t  offset, length;
            std::stringstream ss(line);
            ss >> id >> offset >> length;
            if(i >= this->size){
                Debug(Debug::ERROR) << "Corrupte memory\n";
                EXIT(EXIT_FAILURE);
            }
            index[i].id = strtol(id.c_str(),NULL, 0);
            if(dataMode  ==  DATA_AND_INDEX) {
                index[i].data = data + (offset);
            }else{
                index[i].data = NULL;
            }
            entryLength[i] = length;
            i++;
        }
        index_file.close();
    }else {
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

template <typename T> void DBReader<T>::unmapData() {
    if(dataMapped){
        munmap(data, dataSize);
        dataMapped = false;
    }
}


template class DBReader<unsigned int>;
template class DBReader<std::string>;
