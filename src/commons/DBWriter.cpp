#include <cstdlib>

#include <sstream>
#include <fstream>

#include "DBWriter.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

DBWriter::DBWriter(const char *dataFileName_,
                   const char *indexFileName_,
                   int threads,
                   size_t mode) {
    dataFileName = strdup(dataFileName_);
    indexFileName = strdup(indexFileName_);

    this->threads = threads;

    dataFiles = new FILE *[threads];
    dataFileNames = new char *[threads];

    indexFiles = new FILE *[threads];
    indexFileNames = new char *[threads];

    offsets = new size_t[threads];

    for (int i = 0; i < threads; i++) {
        offsets[i] = 0;
    }

    if (mode == ASCII_MODE) {
        datafileMode = "w";
    } else if (mode == BINARY_MODE) {
        datafileMode = "wb";
    } else {
        Debug(Debug::ERROR) << "No right mode for DBWriter " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }

    closed = 1;
}

DBWriter::~DBWriter() {
    free(dataFileName);
    free(indexFileName);

    delete[] dataFiles;
    delete[] indexFiles;

    for (int i = 0; i < threads; i++) {
        free(dataFileNames[i]);
        free(indexFileNames[i]);
    }
    delete[] dataFileNames;
    delete[] indexFileNames;

    delete[] offsets;
}

void DBWriter::sortDatafileByIdOrder(DBReader<unsigned int> &dbr) {
    Debug(Debug::INFO) << "Sorting the results...  " << dataFileName << " .. ";

#pragma omp parallel for schedule(static)
    for (size_t id = 0; id < dbr.getSize(); id++) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        char *data = dbr.getData(id);
        write(data, strlen(data), SSTR(dbr.getDbKey(id)).c_str(), thread_idx);
    }
    
    Debug(Debug::INFO) << "Done\n";
}

void DBWriter::mergeFiles(DBReader<unsigned int> &qdbr, std::vector<std::pair<std::string, std::string>> files) {
    Debug(Debug::INFO) << "Merging the results... to " << dataFileName << " .. ";

    // open DBReader
    const size_t fileCount = files.size();
    DBReader<unsigned int> *filesToMerge[fileCount];
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i] = new DBReader<unsigned int>(files[i].first.c_str(),
                                                     files[i].second.c_str());
        filesToMerge[i]->open(DBReader<unsigned int>::NOSORT);
    }

    for (size_t id = 0; id < qdbr.getSize(); id++) {
        std::ostringstream ss;
        // get all data for the id from all files
        for (size_t i = 0; i < fileCount; i++) {
            ss << filesToMerge[i]->getData(id);
        }

        // write result
        std::string result = ss.str();
        write(result.c_str(), result.length(), SSTR(qdbr.getDbKey(id)).c_str(), 0);
    }

    // close all reader
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i]->close();
        delete filesToMerge[i];
    }

    Debug(Debug::INFO) << "Done";
}


void DBWriter::swapResults(std::string inputDb, size_t splitSize) {
    std::pair<std::string, std::string> name = Util::databaseNames(inputDb);
    DBReader<unsigned int> dbr(name.first.c_str(), name.second.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);

    char dbKey[255 + 1];
    Debug(Debug::WARNING) << "Start to swap results. Write to " << dataFileName << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string, std::string> > filenames_to_delete;
    std::map<std::string, std::string *> swapMap;
    typedef std::map<std::string, std::string *>::iterator SwapIt;
    for (size_t split = 0; split < splitSize; split++) {
        Debug(Debug::INFO) << "Process split " << split << " ... ";
        // create and open db write
        // create splite file name
        std::string splitName(dataFileName);
        splitName.append("_");
        splitName.append(SSTR(split));
        std::pair<std::string, std::string> splitNames = Util::databaseNames(splitName);
        DBWriter splitWrite(splitNames.first.c_str(), splitNames.second.c_str(), 1);
        splitWrite.open();
        filenames_to_delete.push_back(std::pair<std::string, std::string>(splitNames.first, splitNames.second));

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decomposeDomain(dbr.getSize(), split, splitSize, &startIndex, &domainSize);
        for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
            std::string outerKey = SSTR(dbr.getDbKey(i));
            char *data = dbr.getData(i);
            if (*data == '\0') { // check if file contains entry
                Debug(Debug::ERROR) << "\nSequence " << outerKey
                << " does not contain any sequence!\n";
                continue;
            }

            while (*data != '\0') {
                // extract key from results (ids must be always at the first position)
                Util::parseKey(data, dbKey);
                std::string *entry = NULL;
                SwapIt it = swapMap.find(dbKey);
                if (it == swapMap.end() || it->second == NULL) {
                    entry = new std::string();
                    entry->reserve(1620);
                    swapMap[dbKey] = entry;
                } else {
                    entry = swapMap[dbKey];
                }
                // write data to map
                entry->append(outerKey);
                // next db key
                char *endPosOfId = data + Util::skipNoneWhitespace(data);
                data = Util::skipLine(data);
                entry->append(endPosOfId, data);
            }
        }
        // write results and delete memory
        for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            splitWrite.write((char *) iterator->second->c_str(), iterator->second->size(),
                             (char *) iterator->first.c_str(), 0);
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
    std::string tmp_name(dataFileName);
    tmp_name.append("_all_ids_index");
    std::string tmp_name_index(tmp_name);
    tmp_name_index.append(".index");
    FILE *all_index = fopen(tmp_name_index.c_str(), "w");
    for (SwapIt it = swapMap.begin(); it != swapMap.end(); it++) {
        fprintf(all_index, "%s\t%d\t%d\n", it->first.c_str(), 0, 0);
    }
    fclose(all_index);
    swapMap.clear();

    // make temp. DBReader with all ids
    DBReader<unsigned int> all_ids(inputDb.c_str() /*can be everything that exists */, tmp_name_index.c_str());
    all_ids.open(DBReader<unsigned int>::NOSORT);
    mergeFiles(all_ids, filenames_to_delete);
    all_ids.close();

    remove(tmp_name_index.c_str());
    for (size_t i = 0; i < filenames_to_delete.size(); i++) {
        remove(filenames_to_delete[i].first.c_str());
        remove(filenames_to_delete[i].second.c_str());
    }

}

// allocates heap memory, careful
char* makeResultFilename(const char* name, size_t split) {
    std::stringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void DBWriter::open() {
    for (size_t i = 0; i < threads; i++) {
        dataFileNames[i] = makeResultFilename(dataFileName, i);
        indexFileNames[i] = makeResultFilename(indexFileName, i);

        FileUtil::errorIfFileExist(dataFileNames[i]);
        FileUtil::errorIfFileExist(indexFileNames[i]);

        dataFiles[i] = fopen(dataFileNames[i], datafileMode.c_str());
        indexFiles[i] = fopen(indexFileNames[i], "w");

        if (dataFiles[i] == NULL) {
            perror(dataFileNames[i]);
            EXIT(EXIT_FAILURE);
        }

        if (indexFiles[i] == NULL) {
            perror(indexFileNames[i]);
            EXIT(EXIT_FAILURE);
        }
    }

    closed = 0;
}

void DBWriter::close() {
    // close all datafiles
    for (int i = 0; i < threads; i++) {
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);
    }

    mergeResults(dataFileName, indexFileName,
                 (const char **) dataFileNames, (const char **) indexFileNames, threads);
    closed = 1;
}

void DBWriter::write(const char *data, size_t dataSize, const char *key, int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "ERROR: Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t offsetStart = offsets[thrIdx];
    size_t written = fwrite(data, sizeof(char), dataSize, dataFiles[thrIdx]);
    if (written != dataSize) {
        Debug(Debug::ERROR) << "Could not write to data file " << dataFileName[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
    offsets[thrIdx] += written;

    // entries are always separated by a null byte
    char nullByte = '\0';
    written = fwrite(&nullByte, sizeof(char), 1, dataFiles[thrIdx]);
    if (written != 1) {
        Debug(Debug::ERROR) << "Could not write to data file " << dataFileName[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
    offsets[thrIdx] += 1;

    size_t length = offsets[thrIdx] - offsetStart;
    fprintf(indexFiles[thrIdx], "%s\t%zd\t%zd\n", key, offsetStart, length);
}

void DBWriter::checkClosed() {
    if (closed == 1) {
        Debug(Debug::ERROR) << "Trying to write to a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::mergeResults(const char *outFileName, const char *outFileNameIndex,
                            const char **dataFileNames, const char **indexFileNames, int fileCount) {
    // merge results from each thread into one result file
    // merge each data file
    for(int i = 1; i < fileCount; i++)
    {
        std::ofstream data_file_stream(dataFileNames[0], std::ios_base::binary | std::ios_base::app);
        std::ifstream data_to_add_stream(dataFileNames[i], std::ios_base::binary);
        data_file_stream.seekp(0, std::ios_base::end);
        data_file_stream << data_to_add_stream.rdbuf();
        data_to_add_stream.close();
        if (std::remove(dataFileNames[i]) != 0) {
            Debug(Debug::WARNING) << "Could not remove file " << dataFileNames[i] << "\n";
        }
        data_file_stream.close();
    }

    // rename file to datafile
    std::rename(dataFileNames[0], outFileName);

    FILE *index_file = fopen(outFileNameIndex, "w");
    if (index_file == NULL) {
        perror(outFileNameIndex);
        EXIT(EXIT_FAILURE);
    }

    // merge index
    size_t globalOffset = 0;
    for (int fileIdx = 0; fileIdx < fileCount; fileIdx++) {
        DBReader<std::string> reader(indexFileNames[fileIdx], indexFileNames[fileIdx],
                                     DBReader<std::string>::USE_INDEX);
        reader.open(DBReader<std::string>::NOSORT);
        if(reader.getSize() > 0){
            size_t tmpOffset = 0;
            for(size_t i = 0; i < reader.getSize(); i++){
                const char *id = reader.getIndex()[i].id.c_str();
                size_t currOffset = reinterpret_cast<size_t>(reader.getIndex()[i].data);
                size_t seqLens = reader.getSeqLens(i);
                fprintf(index_file, "%s\t%zd\t%zd\n", id, globalOffset + currOffset, seqLens);
                tmpOffset += reader.getSeqLens(i);
            }
            globalOffset += tmpOffset;
        }
        reader.close();

        if (std::remove(indexFileNames[fileIdx]) != 0) {
            Debug(Debug::WARNING) << "Could not remove file " << indexFileNames[fileIdx] << "\n";
        }
    }
    fclose(index_file);

    // sort the index
    DBReader<std::string> indexReader(outFileNameIndex, outFileNameIndex, DBReader<std::string>::USE_INDEX);
    indexReader.open(DBReader<std::string>::SORT_BY_ID);
    DBReader<std::string>::Index *index = indexReader.getIndex();
    index_file = fopen(outFileNameIndex, "w");
    for (size_t i = 0; i < indexReader.getSize(); i++) {
        const char *id = index[i].id.c_str();
        size_t currOffset = reinterpret_cast<size_t>(index[i].data);
        size_t seqLens = indexReader.getSeqLens(i);
        fprintf(index_file, "%s\t%zd\t%zd\n", id, currOffset, seqLens);
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
    char **buffer = new char *[threads]; //6MB

#pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++) {
        buffer[i] = new char[6400000]; //6MB
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dbSize; i++) {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        unsigned int dbKey = in1.getDbKey(i);
        const char *data1 = in1.getData(i);
        const char *data2 = in2.getData(i);
        size_t entry1Size = in1.getSeqLens(i);
        size_t entry2Size = in2.getSeqLens(i);
        size_t dataSize = entry1Size + entry2Size;
        if(dataSize > 6400000){
            Debug(Debug::ERROR) <<  "Entry " << dbKey << " of " << inIndex2
                                << " and " << inIndex2 << " is " << dataSize
                                << " bytes long. The allowed max size is 102400000 byte. \n";
            EXIT(EXIT_FAILURE);
        }
        memcpy(buffer[thread_idx], data1, entry1Size - 1); // -1 for the nullbyte
        memcpy(buffer[thread_idx] + entry1Size - 1, data2, entry2Size - 1);
        write(buffer[thread_idx], dataSize - 2, SSTR(dbKey).c_str(), thread_idx);
    }
    for (int i = 0; i < threads; i++) {
        delete[] buffer[i];
    }
    delete[] buffer;
    in1.close();
    in2.close();
}
