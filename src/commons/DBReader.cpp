#include "DBReader.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include "Debug.h"
#include "Util.h"

DBReader::DBReader(const char* dataFileName_, const char* indexFileName_,int dataMode /*DATA_AND_INDEX*/)
{
    dataSize = 0;
    this->dataMode = dataMode;
    this->dataFileName = new char [strlen(dataFileName_) + 1];
    memcpy(dataFileName, dataFileName_, sizeof(char) * (strlen(dataFileName_) + 1));
    this->indexFileName = new char [strlen(indexFileName_) + 1];
    memcpy(indexFileName, indexFileName_, sizeof(char) * (strlen(indexFileName_) + 1));
    closed = 1;
}

DBReader::~DBReader(){
    delete[] dataFileName;
    delete[] indexFileName;
}

void DBReader::open(int sort){
    // count the number of entries
    this->size = countLine(indexFileName);

    // open ffindex
    if(dataMode == DATA_AND_INDEX){
        dataFile = fopen(dataFileName, "r");
        if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", dataFileName);  EXIT(EXIT_FAILURE); }
        data = mmapData(dataFile, &dataSize);
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

char * DBReader::mmapData(FILE * file, size_t *dataSize){
    struct stat sb;
    fstat(fileno(file), &sb);
    *dataSize = sb.st_size;
    int fd =  fileno(file);
    return (char*)mmap(NULL, *dataSize, PROT_READ, MAP_PRIVATE, fd, 0);
}

void DBReader::remapData(){
    if(dataMode == DATA_AND_INDEX){
        munmap(data, dataSize);
        data = mmapData(dataFile, &dataSize);
    }
}

void DBReader::close(){
    if(dataMode == DATA_AND_INDEX){
        fclose(dataFile);
        munmap(data, dataSize);
    }
    delete [] index;
    delete [] seqLens;
    closed = 1;
}



size_t DBReader::bsearch(const Index * index, size_t N, unsigned int value)
{
    Index val;
    val.id = value;
    return std::upper_bound(index, index + N, val, Index::compareById) - index;

//    while(key == min){
//
//    }

}



char* DBReader::getData (size_t id){
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

    if (index[id].data - data >= dataSize){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: global id (" << id << ")\n";
        Debug(Debug::ERROR) << "Size of data: " << dataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << (size_t)  (index[id].data - data) << "\n";
        EXIT(EXIT_FAILURE);
    }
    return index[id].data;
}


char* DBReader::getDataByDBKey (const char* dbKey){
    checkClosed();
    if(dataMode ==INDEXONLY){
        Debug(Debug::ERROR) << "DBReader is just open in INDEX_ONLY mode. Call of getData is not allowed" << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t key = strtol(dbKey,NULL, 0);
    size_t id = bsearch(index, size, key);

    return (index[id].id == key) ? index[id].data : NULL;
}

size_t DBReader::getSize (){
    checkClosed();
    return size;
}

std::string DBReader::getDbKey (size_t id){
    checkClosed();
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return std::to_string(index[id].id);
}

size_t DBReader::getId (const char* dbKey){
    size_t key = strtol(dbKey,NULL, 0);
    size_t id = bsearch(index, size, key);
    return (index[id].id == key) ? id : UINT_MAX;
}

unsigned int * DBReader::getSeqLens(){
    return seqLens;
}

size_t DBReader::getSeqLens(size_t id){
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return seqLens[id];
}

void DBReader::checkClosed(){
    if (closed == 1){
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

size_t DBReader::countLine(const char *indexFileName) {
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

void DBReader::readIndex(char *indexFileName, Index * index, char *data, unsigned int * entryLength) {
    std::ifstream index_file(indexFileName);
    std::string line;
    if (index_file.is_open()) {
        size_t i = 0;

        while (std::getline(index_file, line)) {
            size_t id, offset, length;
            std::stringstream ss(line);
            ss >> id >> offset >> length;
            index[i].id = id;
            index[i].data = data + (offset);
            entryLength[i] = length;
            i++;
        }
        index_file.close();
    }else {
        Debug(Debug::ERROR) << "Could not open index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}