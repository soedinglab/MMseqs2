#include "DBWriter.h"

DBWriter::DBWriter (char* dataFileName_, char* indexFileName_, int maxThreadNum_):
    dataFileName(dataFileName_),
    indexFileName(indexFileName_),
    maxThreadNum(maxThreadNum_)
{
    dataFiles = new FILE*[maxThreadNum];
    indexFiles = new FILE*[maxThreadNum];
    offsets = new size_t[maxThreadNum];
    for (int i = 0; i < maxThreadNum; i++)
        offsets[i] = 0;
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
        
        initFFIndexWrite(dataFileNameCStr, indexFileNameCStr, &dataFiles[i], &indexFiles[i]);
    }

}

void DBWriter::close(){
    // merge ffindexes from each thread into one ffindex
    char* merge_command  = (char*) malloc(FILENAME_MAX * 5);
    for(int i = 0; i < maxThreadNum; i++)
    {
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);
        snprintf(merge_command, FILENAME_MAX, "ffindex_build -as %s %s -d %s.%d -i %s.%d",
                dataFileName, indexFileName, dataFileName, i, indexFileName, i);
        system(merge_command);
    }
    delete[] dataFiles;
    delete[] indexFiles;
    delete[] offsets;
}

void DBWriter::write(char* data, int dataSize, char* key, int thrIdx){
    if (thrIdx >= maxThreadNum){
        std::cerr << "ERROR: Thread index " << thrIdx << " > maximum thread number " << maxThreadNum << "\n";
        exit(1);
    }
    ffindex_insert_memory(dataFiles[thrIdx], indexFiles[thrIdx], &offsets[thrIdx], data, dataSize, key);
}

void DBWriter::initFFIndexWrite(char* dataFileName, char* indexFileName, FILE** dataFile, FILE** indexFile){
    
    struct stat st;
    if(stat(dataFileName, &st) == 0) { errno = EEXIST; perror(dataFileName); exit(EXIT_FAILURE); }
    if(stat(indexFileName, &st) == 0) { errno = EEXIST; perror(indexFileName); exit(EXIT_FAILURE); }

    *dataFile = fopen(dataFileName, "w");
    *indexFile = fopen(indexFileName, "w");

    if( *indexFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", indexFileName);  exit(EXIT_FAILURE); }
    if( *dataFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", dataFileName);  exit(EXIT_FAILURE); }
}
