#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <unistd.h>
#include "DBWriter.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Concat.h"
#include "itoa.h"

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
    dataFilesBuffer = new char *[threads];

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
        delete [] dataFilesBuffer[i];
        free(dataFileNames[i]);
        free(indexFileNames[i]);
    }
    delete[] dataFilesBuffer;
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
        writeData(data, strlen(data), SSTR(dbr.getDbKey(id)).c_str(), thread_idx);
    }

    Debug(Debug::INFO) << "Done\n";
}

void DBWriter::mergeFiles(DBReader<unsigned int> &qdbr,
                          const std::vector<std::pair<std::string, std::string>>& files,
                          const std::vector<std::string>& prefixes) {
    Debug(Debug::INFO) << "Merging the results to " << dataFileName << "\n";

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
            if (data != NULL) {
                if(i < prefixes.size()) {
                    ss << prefixes[i];
                }
                ss << filesToMerge[i]->getDataByDBKey(qdbr.getDbKey(id));
            }
        }
        // write result
        std::string result = ss.str();
        writeData(result.c_str(), result.length(), SSTR(qdbr.getDbKey(id)).c_str(), 0);
    }

    // close all reader
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i]->close();
        delete filesToMerge[i];
    }
    delete [] filesToMerge;

    Debug(Debug::INFO) << "Done\n";
}

// allocates heap memory, careful
char* makeResultFilename(const char* name, size_t split) {
    std::stringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void DBWriter::open(size_t bufferSize) {
    for (unsigned int i = 0; i < threads; i++) {
        dataFileNames[i] = makeResultFilename(dataFileName, i);
        indexFileNames[i] = makeResultFilename(indexFileName, i);

        FileUtil::errorIfFileExist(dataFileNames[i]);
        FileUtil::errorIfFileExist(indexFileNames[i]);

        dataFiles[i] = fopen(dataFileNames[i], datafileMode.c_str());
        dataFilesBuffer[i] = new char[bufferSize];
        this->bufferSize = bufferSize;
        // set buffer to 64
        if(setvbuf (dataFiles[i], dataFilesBuffer[i] , _IOFBF , bufferSize ) != 0){
            Debug(Debug::ERROR) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }

        indexFiles[i] = fopen(indexFileNames[i], "w");
        if(setvbuf ( indexFiles[i]  , NULL , _IOFBF , bufferSize ) != 0){
            Debug(Debug::ERROR) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }
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

void DBWriter::writeData(const char *data, size_t dataSize, const char *key, unsigned int thrIdx) {
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

    if(fileCount > 1 ) {
        FILE *outFile = fopen(outFileName, "w");
        FILE **infiles = new FILE *[fileCount];
        for (unsigned int i = 0; i < fileCount; i++) {
            infiles[i] = fopen(dataFileNames[i], "r");
        }
        Concat::concatFiles(infiles, fileCount, outFile);
        for (unsigned int i = 0; i < fileCount; i++) {
            fclose(infiles[i]);
            if (std::remove(dataFileNames[i]) != 0) {
                Debug(Debug::WARNING) << "Could not remove file " << dataFileNames[i] << "\n";
            }
        }
        delete[] infiles;
        fclose(outFile);
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
            if (reader.getSize() > 0) {
                size_t tmpOffset = 0;
                DBReader<std::string>::Index * index = reader.getIndex();
                for (size_t i = 0; i < reader.getSize(); i++) {
                    size_t currOffset = reader.getIndex()[i].offset;
                    index[i].offset = globalOffset + currOffset;
                    tmpOffset += reader.getSeqLens(i);
                }
                writeIndex(index_file, index, reader.getSize(), reader.getSeqLens());
                globalOffset += tmpOffset;
            }
            reader.close();
            if (std::remove(indexFileNames[fileIdx]) != 0) {
                Debug(Debug::WARNING) << "Could not remove file " << indexFileNames[fileIdx] << "\n";
            }
        }
        fclose(index_file);
    }else {  // fileCount < 1
        std::rename(dataFileNames[0], outFileName);
        std::rename(indexFileNames[0], outFileNameIndex);
    }
    // sort the index
    DBReader<std::string> indexReader(outFileNameIndex, outFileNameIndex, DBReader<std::string>::USE_INDEX);
    indexReader.open(DBReader<std::string>::SORT_BY_ID);
    DBReader<std::string>::Index *index = indexReader.getIndex();
    FILE *index_file  = fopen(outFileNameIndex, "w");
    writeIndex(index_file, index, indexReader.getSize(), indexReader.getSeqLens());
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
#if HAVE_POSIX_FADVISE
    if (posix_fadvise (fileno(file1), 0, 0, POSIX_FADV_SEQUENTIAL) != 0){
       Debug(Debug::ERROR) << "posix_fadvise returned an error\n";
    }
    if (posix_fadvise (fileno(file2), 0, 0, POSIX_FADV_SEQUENTIAL) != 0){
       Debug(Debug::ERROR) << "posix_fadvise returned an error\n";
    }
#endif
    int c1, c2;
    char * buffer = dataFilesBuffer[0];
    size_t writePos = 0;
    int dataFilefd =  fileno(dataFiles[0]);
    while((c1=getc_unlocked(file1)) != EOF) {
        if(c1 == '\0'){
            while((c2=getc_unlocked(file2)) != EOF && c2 != '\0') {
                buffer[writePos] = (char) c2;
                writePos++;
                if(writePos == bufferSize){
                    write(dataFilefd, buffer, bufferSize);
                    writePos = 0;
                }
            }
            buffer[writePos] = '\0';
            writePos++;
            if(writePos == bufferSize){
                write(dataFilefd, buffer, bufferSize);
                writePos = 0;
            }
        }else{
            buffer[writePos] = (char) c1;;
            writePos++;
            if(writePos == bufferSize){
                write(dataFilefd, buffer, bufferSize);
                writePos = 0;
            }
        }
    }
    if(writePos != 0){ // if there are data in the buffer that are not yet written
        write(dataFilefd, (const void *) dataFilesBuffer[0], writePos);
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
    DBReader<unsigned int>::Index* index1 = reader1.getIndex();
    unsigned int * seqLen1 = reader1.getSeqLens();
    unsigned int * seqLen2 = reader2.getSeqLens();
    for(size_t id = 0; id < reader1.getSize(); id++){
        // add lenght for file1 and file2 and substrace -1 for one null byte
        size_t seqLen = seqLen1[id] + seqLen2[id] - 1;
        seqLen1[id] = seqLen;
        index1[id].offset = currOffset;
        currOffset += seqLen;
    }
    writeIndex(indexFiles[0], reader1.getSize(), index1, seqLen1);
    reader2.close();
    reader1.close();
}

void DBWriter::writeIndex(FILE *outFile, IndexType::string_type  *index, size_t indexSize,  unsigned int *seqLen){
    char buff1[1024];

    for(size_t id = 0; id < indexSize; id++){
        strncpy(buff1, index[id].id.c_str(), index[id].id.length());
        buff1[index[id].id.length()] = '\t';
        uint64_t currOffset = index[id].offset;
        char * tmpBuff = u64toa_sse2(currOffset, buff1 + (index[id].id.length() + 1) );
        *(tmpBuff-1) = '\t';
        uint32_t sLen = seqLen[id];
        tmpBuff = u32toa_sse2(sLen, tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        fwrite(buff1, sizeof(char), strlen(buff1), outFile);
    }
}

void DBWriter::writeIndex(FILE *outFile,   size_t indexSize, IndexType::int_type  *index,  unsigned int *seqLen){
    char buff1[1024];
    for(size_t id = 0; id < indexSize; id++){
        char * tmpBuff = u32toa_sse2((uint32_t)index[id].id,buff1);
        *(tmpBuff-1) = '\t';
        size_t currOffset = index[id].offset;
        tmpBuff = u64toa_sse2(currOffset, tmpBuff);
        *(tmpBuff-1) = '\t';
        uint32_t sLen = seqLen[id];
        tmpBuff = u32toa_sse2(sLen,tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        fwrite(buff1, sizeof(char), strlen(buff1), outFile);
    }
}
