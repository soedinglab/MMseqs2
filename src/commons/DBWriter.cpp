#include <cstdlib>

#include <sstream>
#include <fstream>
#include <sys/time.h>

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
                   unsigned int threads,
                   size_t mode) {
    dataFileName = strdup(dataFileName_);
    indexFileName = strdup(indexFileName_);

    this->threads = threads;

    dataFiles = new FILE *[threads];
    dataFileNames = new char *[threads];

    indexFiles = new FILE *[threads];
    indexFileNames = new char *[threads];

    offsets = new size_t[threads];

    for (unsigned int i = 0; i < threads; i++) {
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

    closed = true;
}

DBWriter::~DBWriter() {
    free(dataFileName);
    free(indexFileName);

    delete[] dataFiles;
    delete[] indexFiles;

    for (unsigned int i = 0; i < threads; i++) {
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
    DBReader<unsigned int> **filesToMerge = new DBReader<unsigned int>*[fileCount];
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i] = new DBReader<unsigned int>(files[i].first.c_str(),
                                                     files[i].second.c_str());
        filesToMerge[i]->open(DBReader<unsigned int>::NOSORT);
    }

    for (size_t id = 0; id < qdbr.getSize(); id++) {
        std::ostringstream ss;
        // get all data for the id from all files
        for (size_t i = 0; i < fileCount; i++) {
			char *data = filesToMerge[i]->getDataByDBKey(qdbr.getDbKey(id));
			if (data != NULL)
            ss << filesToMerge[i]->getDataByDBKey(qdbr.getDbKey(id));
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
    delete [] filesToMerge;

    Debug(Debug::INFO) << "Done";
}

// allocates heap memory, careful
char* makeResultFilename(const char* name, size_t split) {
    std::stringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void DBWriter::open() {
    for (unsigned int i = 0; i < threads; i++) {
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

    closed = false;
}

void DBWriter::close() {
    // close all datafiles
    for (unsigned int i = 0; i < threads; i++) {
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);
    }

    mergeResults(dataFileName, indexFileName,
                 (const char **) dataFileNames, (const char **) indexFileNames, threads);
    closed = true;
}

void DBWriter::write(const char *data, size_t dataSize, const char *key, unsigned int thrIdx) {
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
    if (closed == true) {
        Debug(Debug::ERROR) << "Trying to write to a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::mergeResults(const char *outFileName, const char *outFileNameIndex,
                            const char **dataFileNames, const char **indexFileNames, unsigned int fileCount) {

    struct timeval start, end;
    gettimeofday(&start, NULL);
    // merge results from each thread into one result file
    // merge each data file
    for(unsigned int i = 1; i < fileCount; i++)
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
    for (unsigned int fileIdx = 0; fileIdx < fileCount; fileIdx++) {
        DBReader<std::string> reader(indexFileNames[fileIdx], indexFileNames[fileIdx],
                                     DBReader<std::string>::USE_INDEX);
        reader.open(DBReader<std::string>::NOSORT);
        if(reader.getSize() > 0){
            size_t tmpOffset = 0;
            for(size_t i = 0; i < reader.getSize(); i++){
                const char *id = reader.getIndex()[i].id.c_str();
                size_t currOffset = reader.getIndex()[i].offset;
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
        size_t currOffset = index[i].offset;
        size_t seqLens = indexReader.getSeqLens(i);
        fprintf(index_file, "%s\t%zd\t%zd\n", id, currOffset, seqLens);
    }
    fclose(index_file);
    indexReader.close();
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for merging files: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) <<" s\n";
}

void DBWriter::mergeFilePair(const char *inData1, const char *inIndex1,
                             const char *inData2, const char *inIndex2) {
    FILE *file1 = fopen(inData1, "r");
    FILE *file2 = fopen(inData2, "r");
    int c1, c2;
    while((c1=fgetc(file1)) != EOF) {
        char file1Char = (char) c1;
        if(file1Char == '\0'){
            while((c2=fgetc(file2)) != EOF && c2 != (int) '\0') {
                char file2Char = (char) c2;
                fwrite(&file2Char, sizeof(char), 1, dataFiles[0]);
            }
            char nullByte = '\0';
            fwrite(&nullByte, sizeof(char), 1, dataFiles[0]);
        }else{
            fwrite(&file1Char, sizeof(char), 1, dataFiles[0]);
        }
    }
    fclose(file1);
    fclose(file2);
    Debug(Debug::WARNING) << "Merge file " << inData1 << " and " << inData2 << "\n";
    DBReader<unsigned int> reader1(inIndex1, inIndex1,
                                 DBReader<std::string>::USE_INDEX);
    reader1.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> reader2(inIndex2, inIndex2,
                                  DBReader<std::string>::USE_INDEX);
    reader2.open(DBReader<unsigned int>::NOSORT);
    size_t currOffset = 0;
    for(size_t id = 0; id < reader1.getSize(); id++){
        // add lenght for file1 and file2 and substrace -1 for one null byte
        size_t seqLen = reader1.getSeqLens(id) + reader2.getSeqLens(id) - 1;
        unsigned int key = reader1.getIndex()[id].id;
        fprintf(indexFiles[0], "%s\t%zd\t%zd\n", SSTR(key).c_str(), currOffset, seqLen);
        currOffset += seqLen;
    }

}
