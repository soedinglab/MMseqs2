#include <cstdlib>

#include <sstream>
#include <fstream>
#include <sys/mman.h>

#include "DBWriter.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

extern "C" {
#include "ffindex.h"
}

#ifdef OPENMP
#include <omp.h>
#endif

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

void DBWriter::sortDatafileByIdOrder(DBReader<unsigned int>*dbr) {
    Debug(Debug::INFO) << "Sorting the results...  " << dataFileName << " .. ";

#pragma omp parallel for schedule(static)
    for (size_t id = 0; id < dbr->getSize(); id++) {
        int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
        char * data = dbr->getData(id);
        this->write(data, strlen(data), SSTR(dbr->getDbKey(id)).c_str(), thread_idx);
    }
    Debug(Debug::INFO) << "Done\n";
}


void DBWriter::mergeFiles(DBReader<unsigned int>* qdbr,
                          std::vector<std::pair<std::string, std::string> > files,
                          size_t maxLineLength)
{
    Debug(Debug::INFO) << "Merging the results... to " << dataFileName << " .. ";
    const size_t file_count = files.size();
    // open DBReader
    DBReader<unsigned int>* filesToMerge [file_count];
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file] = new DBReader<unsigned int>(files[file].first.c_str(),
                                          files[file].second.c_str());
        filesToMerge[file]->open(DBReader<unsigned int>::NOSORT);
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
        this->write(mergeResultsOutData, mergeResultsOutString.length(), SSTR(qdbr->getDbKey(id)).c_str(), 0);
    }
    // close all reader
    for(size_t file = 0; file < file_count; file++){
        filesToMerge[file]->close();
        delete filesToMerge[file];
    }
    Debug(Debug::INFO) << "Done";
}



void DBWriter::swapResults(std::string inputDb, size_t splitSize) {
    DBReader<unsigned int> dbr(inputDb.c_str(), std::string(inputDb+".index").c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);

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
        Debug(Debug::INFO) << "Process split " << split  << " ... ";
        // create and open db write
        DBWriter splitWrite(out_name.c_str(), out_name_index.c_str(), 1);
        splitWrite.open();
        filenames_to_delete.push_back(std::pair<std::string, std::string>(out_name,out_name_index));

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decomposeDomain(dbr.getSize(), split, splitSize, &startIndex, &domainSize);
        for(size_t i = startIndex; i < (startIndex + domainSize); i++){
            std::string outerKey = SSTR(dbr.getDbKey(i));
            char * data = dbr.getData(i);
            if(*data == '\0'){ // check if file contains entry
                Debug(Debug::ERROR) << "\nSequence " << outerKey
                << " does not contain any sequence!\n";
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
        Debug(Debug::INFO) << "Done.\n";
    }
    dbr.close();

    this->open();
    //write new index with all ids (A -> B) of site B
    std::string tmp_name       = std::string(this->dataFileName) + "_all_ids_index";
    std::string tmp_name_index = (tmp_name + ".index");
    FILE* all_index = fopen(tmp_name_index.c_str(), "w");
    for(SwapIt it = swapMap.begin(); it != swapMap.end(); it++) {
        fprintf(all_index, "%s\t%d\t%d\n", it->first.c_str(), 0, 0);
    }
    fclose(all_index);
    swapMap.clear();
    // make temp. DBReader with all ids
    DBReader<unsigned int> all_ids(inputDb.c_str() /*can be everything that exists */, tmp_name_index.c_str());
    all_ids.open(DBReader<unsigned int>::NOSORT);
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

        initFFIndexWrite(dataFileNameCStr, indexFileNameCStr, datafileMode.c_str(), &dataFiles[i], &indexFiles[i]);
    }

    closed = 0;
}

int DBWriter::close(){
    // close all datafiles
    for(int i = 0; i < maxThreadNum; i++){
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);
    }
    mergeFFindexFile(dataFileName, indexFileName, datafileMode.c_str(), (const char **) dataFileNames, (const char **) indexFileNames, maxThreadNum);
    closed = 1;

    return EXIT_SUCCESS;
}

void DBWriter::write(const char* data, int64_t dataSize, const char* key, int thrIdx){
    checkClosed();
    if (thrIdx >= maxThreadNum){
        Debug(Debug::ERROR) <<  "ERROR: Thread index " << thrIdx << " > maximum thread number " << maxThreadNum << "\n";
        EXIT(EXIT_FAILURE);
    }

    // legacy ffindex uses char* instead of const char*, the data is not changed however so the const cast is safe
    ffindex_insert_memory(dataFiles[thrIdx], indexFiles[thrIdx], &offsets[thrIdx], const_cast<char*>(data), dataSize, const_cast<char*>(key));
}


void DBWriter::initFFIndexWrite(const char* dataFileName,
                                const char* indexFileName,
                                const char* datafileMode,
                                FILE** dataFile, FILE** indexFile){
    FileUtil::errorIfFileExist(dataFileName);
    FileUtil::errorIfFileExist(indexFileName);

    *dataFile = fopen(dataFileName, datafileMode);
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

void DBWriter::mergeFFindexFile(const char * outFileName, const char * outFileNameIndex, const char * datafileMode,
                                const char **dataFileNames, const char **indexFileNames, int fileCount ) {
    // merge ffindexes from each thread into one ffindex
    // merge each data file
    std::ofstream data_file_stream(dataFileNames[0], std::ios_base::binary | std::ios_base::app);
    for(int i = 1; i < fileCount; i++)
    {
        std::ifstream data_to_add_stream(dataFileNames[i], std::ios_base::binary);
        data_file_stream.seekp(0, std::ios_base::end);
        data_file_stream << data_to_add_stream.rdbuf();
        data_to_add_stream.close();
        //        fclose(index_file_to_add);
        if (remove(dataFileNames[i]) != 0)
            Debug(Debug::ERROR) << "Error while removing file " << dataFileNames[i] << "\n";
    }
    data_file_stream.close();
    // rename file to datafile
    std::rename(dataFileNames[0],  outFileName);

    // merge index
    FILE* index_file = fopen(outFileNameIndex, "w");
    if( index_file == NULL)  { perror(outFileNameIndex); EXIT(EXIT_FAILURE); }
    size_t globalOffset = 0;
    for(int fileIdx = 0; fileIdx < fileCount; fileIdx++){
        DBReader<std::string> reader(indexFileNames[fileIdx], indexFileNames[fileIdx], DBReader<std::string>::USE_INDEX);
        reader.open(DBReader<std::string>::NOSORT);
        size_t tmpOffset = 0;
        for(size_t i = 0; i < reader.getSize(); i++){
            size_t currOffset = reinterpret_cast<size_t>(reader.getIndex()[i].data);
            fprintf(index_file, "%s\t%zd\t%zd\n", reader.getIndex()[i].id.c_str(), globalOffset + currOffset, reader.getSeqLens(i));
            tmpOffset += reader.getSeqLens(i);
        }
        globalOffset += tmpOffset;
        reader.close();
        if (remove(indexFileNames[fileIdx]) != 0)
            Debug(Debug::ERROR) << "Error while removing file " << indexFileNames[fileIdx] << "\n";
    }
    fclose(index_file);

    // sort the index
    DBReader<std::string> indexReader(outFileNameIndex, outFileNameIndex, DBReader<std::string>::USE_INDEX);
    indexReader.open(DBReader<std::string>::SORT_BY_ID);
    index_file = fopen(outFileNameIndex, "w");
    for(size_t i = 0; i < indexReader.getSize(); i++){
        size_t currOffset = reinterpret_cast<size_t>(indexReader.getIndex()[i].data);
        fprintf(index_file, "%s\t%zd\t%zd\n", indexReader.getIndex()[i].id.c_str(), currOffset, indexReader.getSeqLens(i));
    }

    fclose(index_file);
    indexReader.close();
}

void DBWriter::mergeFilePair(const char *inData1, const char *inIndex1,
                             const char *inData2, const char *inIndex2) {

    DBReader<unsigned int> in1(inData1, inIndex1);
    in1.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> in2(inData2, inIndex2);
    in2.open(DBReader<unsigned int>::NOSORT);
    Debug(Debug::WARNING) << "Merge file " << inData1 << " and " << inData2 << "\n";
    size_t dbSize = in1.getSize();
    char ** buffer = new char*[maxThreadNum]; //6MB

#pragma omp parallel for schedule(static)
    for(int i = 0; i < maxThreadNum; i++){
        buffer[i] = new char[6400000]; //6MB
    }
#pragma omp parallel for schedule(static)
    for(size_t i = 0; i < dbSize; i++){
        int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

        unsigned int dbKey = in1.getDbKey(i);
        const char * data1 = in1.getData(i);
        const char * data2 = in2.getData(i);
        size_t entry1Size = in1.getSeqLens(i);
        size_t entry2Size = in2.getSeqLens(i);
        int64_t dataSize = entry1Size + entry2Size;
        if(dataSize > 6400000){
            Debug(Debug::ERROR) <<  "Entrie " << dbKey << " of " << inIndex2 << " and " << inIndex2 <<" is " << dataSize << " bytes long. "
                                    "The allowed max size is 102400000 byte. \n";
            EXIT(EXIT_FAILURE);
        }
        memcpy(buffer[thread_idx], data1, entry1Size - 1); // -1 for the nullbyte
        memcpy(buffer[thread_idx] + entry1Size -1, data2, entry2Size- 1);
        this->write(buffer[thread_idx], dataSize - 2, SSTR(dbKey).c_str(), thread_idx);
    }
    for(int i = 0; i < maxThreadNum; i++) {
        delete [] buffer[i];
    }
    delete [] buffer;
    in1.close();
    in2.close();
}
