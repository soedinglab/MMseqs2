#include "DBWriter.h"

DBWriter::DBWriter (const char* dataFileName_, const char* indexFileName_, int maxThreadNum_)
{
    this->dataFileName = new char [strlen(dataFileName_) + 1];
    memcpy(dataFileName, dataFileName_, sizeof(char) * (strlen(dataFileName_) + 1));

    this->indexFileName = new char [strlen(indexFileName_) + 1];
    memcpy(indexFileName, indexFileName_, sizeof(char) * (strlen(indexFileName_) + 1));
    this->maxThreadNum = maxThreadNum_;
    dataFiles = new FILE*[maxThreadNum];
    dataFileNames = new char*[maxThreadNum];
    indexFiles = new FILE*[maxThreadNum];
    indexFileNames = new char*[maxThreadNum];
    offsets = new size_t[maxThreadNum];
    for (int i = 0; i < maxThreadNum; i++)
        offsets[i] = 0;
    closed = 1;
}

DBWriter::~DBWriter(){
    delete[] dataFileName;
    delete[] indexFileName;
    delete[] dataFiles;
    delete[] indexFiles;
    delete[] dataFileNames;
    delete[] indexFileNames;
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

    closed = 0;
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
            // count the number of entries
            char line[1000];
            int cnt = 0;
            std::ifstream index_file_cnt(indexFileNames[i]);
            if (index_file_cnt.is_open()) {
                while ( index_file_cnt.getline (line, 1000) ){
                    cnt++;
                }
                index_file_cnt.close();
            }
            else{
                std::cerr << "Could not open ffindex index file " << indexFileNames[i] << "\n";
                exit(EXIT_FAILURE);
            }
            // merge data and indexes
            ffindex_index_t* index_to_add = ffindex_index_parse(index_file_to_add, cnt);
            ffindex_insert_ffindex(data_file, index_file, &offset, data_to_add, index_to_add);
            free(index_to_add);
        }
        munmap(data_to_add, data_size);

        fclose(data_file_to_add);
        fclose(index_file_to_add);
        if (remove(dataFileNames[i]) != 0)
            std::cerr << "Error while removing file " << dataFileNames[i] << "\n";
        if (remove(indexFileNames[i]) != 0)
            std::cerr << "Error while removing file " << indexFileNames[i] << "\n";

    }
    fclose(data_file);
    fclose(index_file);

    // sort the index file
    char line[1000];
    int cnt = 0;
    std::ifstream index_file_cnt(indexFileName);
    if (index_file_cnt.is_open()) {
        while ( index_file_cnt.getline (line, 1000) ){
            cnt++;
        }
        index_file_cnt.close();
    }
    else{
        std::cerr << "Could not open ffindex index file " << indexFileName << "\n";
        exit(EXIT_FAILURE);
    }

    index_file = fopen(indexFileName, "r+");
    if(index_file == NULL) { perror(indexFileName); return EXIT_FAILURE; }

    ffindex_index_t* index = ffindex_index_parse(index_file, cnt);
    if(index == NULL) { perror("ffindex_index_parse failed"); return (EXIT_FAILURE); }

    fclose(index_file);

    ffindex_sort_index_file(index);
    index_file = fopen(indexFileName, "w");
    if(index_file == NULL) { perror(indexFileName); return EXIT_FAILURE; }
    ffindex_write(index, index_file);
    fclose(index_file);
    free(index);

    for (int i = 0; i < maxThreadNum; i++){
        delete[] dataFileNames[i];
        delete[] indexFileNames[i];
    }

    closed = 1;

    return EXIT_SUCCESS;
}

void DBWriter::write(char* data, int dataSize, char* key, int thrIdx){
    checkClosed();
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

    if( *dataFile == NULL)  { perror(dataFileName); exit(EXIT_FAILURE); } 
    if( *indexFile == NULL) { perror(indexFileName); exit(EXIT_FAILURE); }
}

void DBWriter::checkClosed(){
    if (closed == 1){
        std::cerr << "Trying to write to a closed database.\n";
        exit(EXIT_FAILURE);
    }
}
