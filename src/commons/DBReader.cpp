#include "DBReader.h"

DBReader::DBReader(const char* dataFileName_, const char* indexFileName_):
    dataFileName(dataFileName_),
    indexFileName(indexFileName_)
{
    dataSize = 0;
}

void DBReader::open(){
    // open ffindex
    dataFile = fopen(dataFileName, "r");
    indexFile = fopen(indexFileName, "r");

    if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", indexFileName);  exit(EXIT_FAILURE); }
    if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", dataFileName);  exit(EXIT_FAILURE); }

    data = ffindex_mmap_data(dataFile, &dataSize);

    index = ffindex_index_parse(indexFile, 0);

    if(index == NULL)
    {
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", indexFileName);
        exit(EXIT_FAILURE);
    }

    size = index->n_entries;

    seqLens = new unsigned short [size];

    for (size_t i = 0; i < size; i++){
        seqLens[i] = (unsigned short)(ffindex_get_entry_by_index(index, i)->length);
    }
}

void DBReader::close(){
    fclose(dataFile);
    fclose(indexFile);
}

char* DBReader::getData (int id){
    if (id >= size){
        std::cerr << "getData: id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }
    if (ffindex_get_entry_by_index(index, id)->offset >= dataSize){ 
        std::cerr << "Invalid database read for id=" << id << "\n";
        exit(1);
    }
    return data + (ffindex_get_entry_by_index(index, id)->offset);
}

char* DBReader::getDataByDBKey (char* key){
    return ffindex_get_data_by_name(data, index, key);
}

size_t DBReader::getSize (){
    return size;
}

char* DBReader::getDbKey (int id){
    if (id >= size){
        std::cerr << "getDbKey: id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }
    return &(ffindex_get_entry_by_index(index, id)->name[0]);
}
