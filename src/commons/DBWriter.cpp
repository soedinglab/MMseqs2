#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

DBWriter::DBWriter (const char* dataFileName_,
                    const char* indexFileName_,
                    int maxThreadNum_,
                    size_t mode /* default ASCII */)
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
    if(mode == ASCII_MODE)
        datafileMode = "w";
    else if(mode == BINARY_MODE)
        datafileMode = "wb";
    else {
        Debug(Debug::ERROR) <<  "No right mode for DBWriter " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    
    closed = 1;
}

DBWriter::~DBWriter(){
    delete[] dataFileName;
    delete[] indexFileName;
    delete[] dataFiles;
    delete[] indexFiles;
    for (int i = 0; i < maxThreadNum; i++){
        delete[] dataFileNames[i];
        delete[] indexFileNames[i];
    }
    delete[] dataFileNames;
    delete[] indexFileNames;
    delete[] offsets;
}




void DBWriter::sortDatafileByIdOrder(DBReader *dbr) {
    Debug(Debug::INFO) << "Sorting the results...  " << dataFileName << " .. ";

#pragma omp parallel for schedule(static)
    for (size_t id = 0; id < dbr->getSize(); id++) {
        int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
        char * data = dbr->getData(id);
        this->write(data, strlen(data), dbr->getDbKey(id), thread_idx);
    }
    Debug(Debug::INFO) << "Done\n";
}


void DBWriter::mergeFiles(DBReader * qdbr,
                          std::vector<std::pair<std::string, std::string> > files,
                          size_t maxLineLength)
{
    Debug(Debug::INFO) << "Merging the results... to " << dataFileName << " .. ";
    const size_t file_count = files.size();

    // open DBReader
    DBReader * filesToMerge [file_count];
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file] = new DBReader(files[file].first.c_str(),
                                          files[file].second.c_str());
        filesToMerge[file]->open(DBReader::NOSORT);
    }
    for (size_t id = 0; id < qdbr->getSize(); id++){
        std::string mergeResultsOutString;
        mergeResultsOutString.reserve(maxLineLength);
        // get all data for the id from all files
        for(size_t file = 0; file < file_count; file++){
            mergeResultsOutString.append(filesToMerge[file]->getData(id));
        }
        // create merged string
        if (mergeResultsOutString.length() >= maxLineLength ){
            Debug(Debug::ERROR) << "ERROR: Buffer overflow at id: " << qdbr->getDbKey(id) << " during the merging.\n";
            Debug(Debug::ERROR) << "Output buffer size < prefiltering result size! (" << maxLineLength << " < " << mergeResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            continue; // read next id
        }
        // write result
        char* mergeResultsOutData = (char *) mergeResultsOutString.c_str();
        this->write(mergeResultsOutData, mergeResultsOutString.length(), qdbr->getDbKey(id), 0);
    }
    // close all reader
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file]->close();
        delete filesToMerge[file];
    }
    Debug(Debug::INFO) << "Done";
}



void DBWriter::swapResults(std::string inputDb, size_t splitSize) {

    DBReader dbr(inputDb.c_str(), std::string(inputDb+".index").c_str());
    dbr.open(DBReader::NOSORT);

    char dbKey[255+1];
    Debug(Debug::WARNING) << "Start to swap results. Write to " << this->dataFileName << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string,std::string > > filenames_to_delete;
    std::map<std::string, std::string * > swapMap;
    typedef std::map<std::string, std::string * >::iterator SwapIt;
    for(size_t split = 0; split < splitSize; split++){
        // create splite file name
        std::string out_name       = std::string(this->dataFileName) + "_" + SSTR(split);
        std::string out_name_index = (out_name + ".index");
        std::cout << "Process split " << split  << " ... ";
        // create and open db write
        DBWriter splitWrite(out_name.c_str(), out_name_index.c_str(), 1);
        splitWrite.open();
        filenames_to_delete.push_back(std::pair<std::string, std::string>(out_name,out_name_index));

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decompose_domain(dbr.getSize(), split, splitSize, &startIndex, &domainSize);
        for(size_t i = startIndex; i < (startIndex + domainSize); i++){
            char * outerKey = dbr.getDbKey(i);
            char * data = dbr.getData(i);
            if(*data == '\0'){ // check if file contains entry
                Debug(Debug::ERROR) << "\nSequence " << outerKey
                        << " does not containe any sequence!\n";
                continue;
            }

            while (*data != '\0')
            {
                // extract key from results (ids must be always at the first position)
                Util::parseKey(data, dbKey);
                std::string * entry = NULL;
                SwapIt it = swapMap.find(dbKey);
                if(it == swapMap.end()|| it->second == NULL){
                    entry = new std::string();
                    entry->reserve(1620);
                    swapMap[dbKey] = entry;
                }else{
                    entry = swapMap[dbKey];
                }
                // write data to map
                entry->append(outerKey);
                // next db key
                char * endPosOfId    = data + Util::skipNoneWhitespace(data);
                data = Util::skipLine(data);
                entry->append(endPosOfId, data);
            }
        }
        // write results and delete memory
        for(SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            splitWrite.write((char *) iterator->second->c_str(),  iterator->second->size(), (char *) iterator->first.c_str() ,0);
            entries_num++;
            // remove just the value (string *) not the keys
            // the keys are needed for the merging step later
            delete iterator->second;
            iterator->second = NULL;
        }
        splitWrite.close();
        std::cout << "Done." << std::endl;
    }
    dbr.close();

    this->open();
    //write new index with all ids (A -> B) of site B
    std::string tmp_name       = std::string(this->dataFileName) + "_all_ids_index";
    std::string tmp_name_index = (tmp_name + ".index");
    FILE* all_index = fopen(tmp_name_index.c_str(), "w");
    for(SwapIt it = swapMap.begin(); it != swapMap.end(); it++) {
        fprintf(all_index, "%s\t%zd\t%zd\n", it->first.c_str(), 0, 0);
    }
    fclose(all_index);
    swapMap.clear();
    // make temp. DBReader with all ids
    DBReader all_ids(inputDb.c_str() /*can be everything that exists */, tmp_name_index.c_str());
    all_ids.open(DBReader::NOSORT);
    this->mergeFiles(&all_ids, filenames_to_delete, 1000000);
    all_ids.close();
    remove(tmp_name_index.c_str());
    for (size_t i = 0; i < filenames_to_delete.size(); i++) {
        remove(filenames_to_delete[i].first.c_str());
        remove(filenames_to_delete[i].second.c_str());
    }

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

        FILE* data_file_to_add  = fopen(dataFileNames[i], "r");
        if( data_file_to_add == NULL) { perror(dataFileNames[i]); return EXIT_FAILURE; }
        FILE* index_file_to_add = fopen(indexFileNames[i], "r");
        if( index_file_to_add == NULL) { perror(indexFileNames[i]); return EXIT_FAILURE; }
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
                EXIT(EXIT_FAILURE);
            }
            // merge data and indexes
            ffindex_index_t* index_to_add = ffindex_index_parse(index_file_to_add, cnt);
            ffindex_insert_ffindex(data_file, index_file, &offset, data_to_add, index_to_add);
            munmap(index_to_add->index_data, index_to_add->index_data_size);
            free(index_to_add);
            munmap(data_to_add,data_size);
        }

        fclose(data_file_to_add);
        fclose(index_file_to_add);
        if (remove(dataFileNames[i]) != 0)
            Debug(Debug::ERROR) << "Error while removing file " << dataFileNames[i] << "\n";
        if (remove(indexFileNames[i]) != 0)
            Debug(Debug::ERROR) << "Error while removing file " << indexFileNames[i] << "\n";

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
        Debug(Debug::ERROR) <<  "Could not open ffindex index file " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
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

    closed = 1;

    return EXIT_SUCCESS;
}

void DBWriter::writeFile(FILE * file, char* key, int thrIdx){
    ffindex_insert_filestream(dataFiles[thrIdx], indexFiles[thrIdx],
                                 &offsets[thrIdx], file, key);

}

void DBWriter::write(char* data, int64_t dataSize, char* key, int thrIdx){
    checkClosed();
    if (thrIdx >= maxThreadNum){
        Debug(Debug::ERROR) <<  "ERROR: Thread index " << thrIdx << " > maximum thread number " << maxThreadNum << "\n";
        EXIT(1);
    }
    ffindex_insert_memory(dataFiles[thrIdx], indexFiles[thrIdx], &offsets[thrIdx], data, dataSize, key);
}

void DBWriter::errorIfFileExist(const char * file){
    struct stat st;
    if(stat(file, &st) == 0) { errno = EEXIST; perror(file); EXIT(EXIT_FAILURE); }
}

void DBWriter::initFFIndexWrite(const char* dataFileName,
                                const char* indexFileName,
                                FILE** dataFile, FILE** indexFile){
    DBWriter::errorIfFileExist(dataFileName);
    DBWriter::errorIfFileExist(indexFileName);

    *dataFile = fopen(dataFileName, datafileMode.c_str());
    *indexFile = fopen(indexFileName, "w");

    if( *dataFile == NULL)  { perror(dataFileName); EXIT(EXIT_FAILURE); }
    if( *indexFile == NULL) { perror(indexFileName); EXIT(EXIT_FAILURE); }
}

void DBWriter::checkClosed(){
    if (closed == 1){
        Debug(Debug::ERROR) <<  "Trying to write to a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}
