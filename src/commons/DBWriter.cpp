#include "DBWriter.h"

DBWriter::DBWriter (const char* dataFileName_, const char* indexFileName_, int maxThreadNum_):
    dataFileName(dataFileName_),
    indexFileName(indexFileName_),
    maxThreadNum(maxThreadNum_)
{
    dataFiles = new FILE*[maxThreadNum];
    dataFileNames = new char*[maxThreadNum];
    indexFiles = new FILE*[maxThreadNum];
    indexFileNames = new char*[maxThreadNum];
    offsets = new size_t[maxThreadNum];
    for (int i = 0; i < maxThreadNum; i++)
        offsets[i] = 0;
}

DBWriter::~DBWriter(){
    delete[] dataFiles;
    delete[] indexFiles;
    delete[] offsets;
}


void DBWriter::open(){

    for (int i = 0; i < maxThreadNum; i++){
        std::stringstream dataFileNameSs;
        dataFileNameSs << dataFileName << "." << i;

        std::stringstream indexFileNameSs;
        indexFileNameSs << indexFileName << "." << i;

        std::string dataFileNameStr = dataFileNameSs.str();
        std::string indexFileNameStr = indexFileNameSs.str();

        char* dataFileNameCStr = (char*)dataFileNameStr.c_str();
        char* indexFileNameCStr = (char*)indexFileNameStr.c_str();

        dataFileNames[i] = new char[strlen(dataFileNameCStr)+1];
        indexFileNames[i] = new char[strlen(indexFileNameCStr)+1];
        strcpy(dataFileNames[i], dataFileNameCStr);
        strcpy(indexFileNames[i], indexFileNameCStr);

        initFFIndexWrite(dataFileNameCStr, indexFileNameCStr, &dataFiles[i], &indexFiles[i]);
    }

}

int DBWriter::close(){
    // merge ffindexes from each thread into one ffindex
    
    FILE* data_file;
    FILE* index_file;

    initFFIndexWrite(dataFileName, indexFileName, &data_file, &index_file);

    size_t offset = 0;
    for(int i = 0; i < maxThreadNum; i++)
    {
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);

        FILE* data_file_to_add  = fopen(dataFileNames[i], "r");  if(  data_file_to_add == NULL) { perror(dataFileNames[i]); return EXIT_FAILURE; }
        FILE* index_file_to_add = fopen(indexFileNames[i], "r"); if( index_file_to_add == NULL) { perror(indexFileNames[i]); return EXIT_FAILURE; }
        size_t data_size;
        char *data_to_add = ffindex_mmap_data(data_file_to_add, &data_size);
        if (data_size > 0){
            ffindex_index_t* index_to_add = ffindex_index_parse(index_file_to_add, 0);
            ffindex_insert_ffindex(data_file, index_file, &offset, data_to_add, index_to_add);
        }

        fclose(data_file_to_add);
        fclose(index_file_to_add);
    }
    fclose(data_file);
    fclose(index_file);

    return EXIT_SUCCESS;
}

void DBWriter::write(char* data, int dataSize, char* key, int thrIdx){
    if (thrIdx >= maxThreadNum){
        std::cerr << "ERROR: Thread index " << thrIdx << " > maximum thread number " << maxThreadNum << "\n";
        exit(1);
    }
    ffindex_insert_memory(dataFiles[thrIdx], indexFiles[thrIdx], &offsets[thrIdx], data, dataSize, key);
}

void DBWriter::initFFIndexWrite(const char* dataFileName, const char* indexFileName, FILE** dataFile, FILE** indexFile){
    
    struct stat st;
    if(stat(dataFileName, &st) == 0) { errno = EEXIST; perror(dataFileName); exit(EXIT_FAILURE); }
    if(stat(indexFileName, &st) == 0) { errno = EEXIST; perror(indexFileName); exit(EXIT_FAILURE); }

    *dataFile = fopen(dataFileName, "w");
    *indexFile = fopen(indexFileName, "w");

    if( *indexFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", indexFileName);  exit(EXIT_FAILURE); }
    if( *dataFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", dataFileName);  exit(EXIT_FAILURE); }
}
