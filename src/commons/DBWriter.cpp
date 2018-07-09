#include "DBWriter.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Concat.h"
#include "itoa.h"

#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <sys/time.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

DBWriter::DBWriter(const char *dataFileName_, const char *indexFileName_, unsigned int threads, size_t mode)
        : threads(threads), mode(mode) {
    dataFileName = strdup(dataFileName_);
    indexFileName = strdup(indexFileName_);

    dataFiles = new FILE *[threads];
    dataFilesBuffer = new char *[threads];

    dataFileNames = new char *[threads];

    indexFiles = new FILE *[threads];
    indexFileNames = new char *[threads];

    starts = new size_t[threads];
    std::fill(starts, starts + threads, 0);
    offsets = new size_t[threads];
    std::fill(offsets, offsets + threads, 0);

    if ((mode & BINARY_MODE) != 0) {
        datafileMode = "wb";
    } else {
        datafileMode = "w";
    }

    closed = true;
}

DBWriter::~DBWriter() {
    delete[] offsets;
    delete[] starts;
    delete[] indexFileNames;
    delete[] indexFiles;
    delete[] dataFileNames;
    delete[] dataFilesBuffer;
    delete[] dataFiles;
    free(indexFileName);
    free(dataFileName);
}

void DBWriter::sortDatafileByIdOrder(DBReader<unsigned int> &dbr) {
    Debug(Debug::INFO) << "Sorting the results...  " << dataFileName << " .. ";

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

#pragma omp for schedule(static)
        for (size_t id = 0; id < dbr.getSize(); id++) {
            char *data = dbr.getData(id);
            size_t length = dbr.getSeqLens(id);
            writeData(data, (length == 0 ? 0 : length - 1), dbr.getDbKey(id), thread_idx);
        }
    };
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
        unsigned int key = qdbr.getDbKey(id);
        std::ostringstream ss;
        // get all data for the id from all files
        for (size_t i = 0; i < fileCount; i++) {
            const char *data = filesToMerge[i]->getDataByDBKey(key);
            if (data != NULL) {
                if(i < prefixes.size()) {
                    ss << prefixes[i];
                }
                ss << data;
            }
        }
        // write result
        std::string result = ss.str();
        writeData(result.c_str(), result.length(), key, 0);
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
    std::ostringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void DBWriter::open(size_t bufferSize) {
    for (unsigned int i = 0; i < threads; i++) {
        dataFileNames[i] = makeResultFilename(dataFileName, i);
        indexFileNames[i] = makeResultFilename(indexFileName, i);

        dataFiles[i] = fopen(dataFileNames[i], datafileMode.c_str());
        if (dataFiles[i] == NULL) {
            Debug(Debug::ERROR) << "Could not open " << dataFileNames[i] << " for writing!\n";
            EXIT(EXIT_FAILURE);
        }

        int fd = fileno(dataFiles[i]);
        int flags;
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Could not set mode for " << dataFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        dataFilesBuffer[i] = new char[bufferSize];
        this->bufferSize = bufferSize;

        // set buffer to 64
        if (setvbuf(dataFiles[i], dataFilesBuffer[i], _IOFBF, bufferSize) != 0) {
            Debug(Debug::WARNING) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }

        indexFiles[i] = fopen(indexFileNames[i], "w");
        if (indexFiles[i] == NULL) {
            Debug(Debug::ERROR) << "Could not open " << indexFileNames[i] << " for writing!\n";
            EXIT(EXIT_FAILURE);
        }

        fd = fileno(indexFiles[i]);
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Could not set mode for " << indexFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        if (setvbuf(indexFiles[i], NULL, _IOFBF, bufferSize) != 0) {
            Debug(Debug::WARNING) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
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

void DBWriter::close(int dbType) {
    // close all datafiles
    for (unsigned int i = 0; i < threads; i++) {
        if (fsync(fileno(dataFiles[i])) < 0) {
            Debug(Debug::ERROR) << "Could not sync data file " << dataFileNames[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(dataFiles[i]);
        if (fsync(fileno(indexFiles[i])) < 0) {
            Debug(Debug::ERROR) << "Could not sync index file " << indexFileNames[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(indexFiles[i]);
    }

    if (dbType > -1){
        std::string dataFile = dataFileName;
        std::string dbTypeFile = (dataFile+".dbtype").c_str();
        FILE * dbtypeDataFile = fopen(dbTypeFile.c_str(), "wb");
        if (dbtypeDataFile == NULL) {
            Debug(Debug::ERROR) << "Could not open data file " << dbTypeFile << "!\n";
            EXIT(EXIT_FAILURE);
        }
        size_t written = fwrite(&dbType, sizeof(int), 1, dbtypeDataFile);
        if (written != 1) {
            Debug(Debug::ERROR) << "Could not write to data file " << dbTypeFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fsync(fileno(dbtypeDataFile)) < 0) {
            Debug(Debug::ERROR) << "Could not sync type file " << dbTypeFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(dbtypeDataFile);
    }

    mergeResults(dataFileName, indexFileName,
                 (const char **) dataFileNames, (const char **) indexFileNames, threads, ((mode & LEXICOGRAPHIC_MODE) != 0));

    for (unsigned int i = 0; i < threads; i++) {
        delete [] dataFilesBuffer[i];
        free(dataFileNames[i]);
        free(indexFileNames[i]);
    }
    closed = true;
}

void DBWriter::writeStart(unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "ERROR: Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }

    starts[thrIdx] = offsets[thrIdx];
}

void DBWriter::writeAdd(const char* data, size_t dataSize, unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "ERROR: Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }

    size_t written = fwrite(data, sizeof(char), dataSize, dataFiles[thrIdx]);
    if (written != dataSize) {
        Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
    offsets[thrIdx] += written;
}

void DBWriter::writeEnd(unsigned int key, unsigned int thrIdx, bool addNullByte) {
    size_t written;
    // entries are always separated by a null byte
    if(addNullByte == true){
        char nullByte = '\0';
        written = fwrite(&nullByte, sizeof(char), 1, dataFiles[thrIdx]);
        if (written != 1) {
            Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        offsets[thrIdx] += 1;
    }

    size_t length = offsets[thrIdx] - starts[thrIdx];

    char buffer[1024];
    size_t len = indexToBuffer(buffer, key, starts[thrIdx], length );
    written = fwrite(buffer, sizeof(char), len, indexFiles[thrIdx]);
    if (written != len) {
        Debug(Debug::ERROR) << "Could not write to data file " << indexFiles[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::writeData(const char *data, size_t dataSize, unsigned int key, unsigned int thrIdx, bool addNullByte) {
    writeStart(thrIdx);
    writeAdd(data, dataSize, thrIdx);
    writeEnd(key, thrIdx, addNullByte);
}

size_t DBWriter::indexToBuffer(char *buff1, unsigned int key, size_t offsetStart, size_t len){
    char * basePos = buff1;
    char * tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(key), buff1);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(offsetStart), tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(len), tmpBuff);
    *(tmpBuff-1) = '\n';
    *(tmpBuff) = '\0';
    return tmpBuff - basePos;
}

void DBWriter::alignToPageSize() {
    if (threads > 1) {
        Debug(Debug::ERROR) << "Data file can only be aligned in single threaded mode.\n";
        EXIT(EXIT_FAILURE);
    }

    size_t currentOffset = offsets[0];
    size_t pageSize = Util::getPageSize();
    size_t newOffset = ((pageSize - 1) & currentOffset) ? ((currentOffset + pageSize) & ~(pageSize - 1)) : currentOffset;
    char nullByte = '\0';
    for (size_t i = currentOffset; i < newOffset; ++i) {
        size_t written = fwrite(&nullByte, sizeof(char), 1, dataFiles[0]);
        if (written != 1) {
            Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[0] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    offsets[0] = newOffset;
}


void DBWriter::checkClosed() {
    if (closed == true) {
        Debug(Debug::ERROR) << "Trying to write to a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                            const std::vector<std::pair<std::string, std::string >> &files,
                            const  bool lexicographicOrder) {
    const char **datafilesNames = new const char *[files.size()];
    const char **indexFilesNames = new const char *[files.size()];
    for (size_t i = 0; i < files.size(); i++) {
        datafilesNames[i] = files[i].first.c_str();
        indexFilesNames[i] = files[i].second.c_str();
    }
    mergeResults(outFileName.c_str(), outFileNameIndex.c_str(), datafilesNames, indexFilesNames, files.size(), lexicographicOrder);
    delete[] datafilesNames;
    delete[] indexFilesNames;

}

template <>
void DBWriter::writeIndex(FILE *outFile, size_t indexSize, DBReader<unsigned int>::Index *index,  unsigned int *seqLen){
    char buff1[1024];
    for (size_t id = 0; id < indexSize; id++) {
        char * tmpBuff = Itoa::u32toa_sse2((uint32_t)index[id].id,buff1);
        *(tmpBuff-1) = '\t';
        size_t currOffset = index[id].offset;
        tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
        *(tmpBuff-1) = '\t';
        uint32_t sLen = seqLen[id];
        tmpBuff = Itoa::u32toa_sse2(sLen,tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
    }
}

template <>
void DBWriter::writeIndex(FILE *outFile, size_t indexSize, DBReader<std::string>::Index *index,  unsigned int *seqLen){
    char buff1[1024];
    for (size_t id = 0; id < indexSize; id++) {
        size_t keyLen = index[id].id.length();
        char * tmpBuff = (char*)memcpy((void*)buff1, (void*)index[id].id.c_str(), keyLen);
        tmpBuff+=keyLen;
        *(tmpBuff) = '\t';
        tmpBuff++;
        size_t currOffset = index[id].offset;
        tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
        *(tmpBuff-1) = '\t';
        uint32_t sLen = seqLen[id];
        tmpBuff = Itoa::u32toa_sse2(sLen,tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
    }
}


void DBWriter::mergeResults(const char *outFileName, const char *outFileNameIndex,
                            const char **dataFileNames, const char **indexFileNames,
                            const unsigned long fileCount, const bool lexicographicOrder) {
    if (fileCount <= 1) {
        if (std::rename(dataFileNames[0], outFileName) != 0) {
            Debug(Debug::ERROR) << "Could not move result " << dataFileNames[0] << " to final location " << outFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        if (std::rename(indexFileNames[0], outFileNameIndex) != 0) {
            Debug(Debug::ERROR) << "Could not move result index " << indexFileNames[0] << " to final location " << outFileNameIndex << "!\n";
            EXIT(EXIT_FAILURE);
        }
        return;
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);
    // merge result data files into the first result data file
    std::vector<size_t> dataSizes;
    FILE **files = new FILE*[fileCount];
    for (unsigned int i = 0; i < fileCount; ++i) {
        files[i] = fopen(dataFileNames[i], i == 0 ? "a" : "r");
        if (files[i] == NULL) {
            Debug(Debug::ERROR) << "Could not open result file " << dataFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        struct stat sb;
        if (fstat(fileno(files[i]), &sb) < 0) {
            int errsv = errno;
            Debug(Debug::ERROR) << "Failed to fstat file " << dataFileNames[i] << ". Error " << errsv << ".\n";
            EXIT(EXIT_FAILURE);
        }
        dataSizes.push_back(sb.st_size);
    }
    Concat::concatFiles(files + 1, fileCount - 1, files[0]);
    for (unsigned int i = 0; i < fileCount; ++i) {
        if (fsync(fileno(files[i])) < 0) {
            Debug(Debug::ERROR) << "Could not sync data file " << dataFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(files[i]);
        if (i == 0) {
            if (std::rename(dataFileNames[i], outFileName) != 0) {
                Debug(Debug::ERROR) << "Could not move result " << dataFileNames[i] << " to final location " << outFileName << "!\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            if (std::remove(dataFileNames[i]) != 0) {
                Debug(Debug::WARNING) << "Could not remove file " << dataFileNames[i] << "\n";
            }
        }
    }
    delete[] files;

    // merge index files
    FILE *indexFile = fopen(outFileNameIndex, "w");
    if (indexFile == NULL) {
        perror(outFileNameIndex);
        EXIT(EXIT_FAILURE);
    }
    size_t globalOffset = 0;
    for (unsigned int fileIdx = 0; fileIdx < fileCount; ++fileIdx) {
        DBReader<unsigned int> reader(indexFileNames[fileIdx], indexFileNames[fileIdx], DBReader<unsigned int>::USE_INDEX);
        reader.open(DBReader<unsigned int>::HARDNOSORT);
        if (reader.getSize() > 0) {
            DBReader<unsigned int>::Index * index = reader.getIndex();
            for (size_t i = 0; i < reader.getSize() - 1; ++i) {
                size_t currOffset = index[i].offset;
                index[i].offset = globalOffset + currOffset;
            }

            size_t currOffset = index[reader.getSize() - 1].offset;
            index[reader.getSize() - 1].offset = globalOffset + currOffset;
            writeIndex(indexFile, reader.getSize(), index, reader.getSeqLens());
            globalOffset += dataSizes[fileIdx];
        }
        reader.close();

        if (std::remove(indexFileNames[fileIdx]) != 0) {
            Debug(Debug::WARNING) << "Could not remove file " << indexFileNames[fileIdx] << "\n";
        }
    }
    if (fsync(fileno(indexFile)) < 0) {
        Debug(Debug::ERROR) << "Could not sync index file " << outFileNameIndex << "!\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(indexFile);


    if (lexicographicOrder == false) {
        // sort the index
        DBReader<unsigned int> indexReader(outFileNameIndex, outFileNameIndex, DBReader<unsigned int>::USE_INDEX);
        indexReader.open(DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int>::Index *index = indexReader.getIndex();
        FILE *indexFile  = fopen(outFileNameIndex, "w");
        if (indexFile == NULL) {
            perror(outFileNameIndex);
            EXIT(EXIT_FAILURE);
        }
        writeIndex(indexFile, indexReader.getSize(), index, indexReader.getSeqLens());
        if (fsync(fileno(indexFile)) < 0) {
            Debug(Debug::ERROR) << "Could not sync index file " << outFileNameIndex << "!\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(indexFile);
        indexReader.close();
    } else {
        DBReader<std::string> indexReader(outFileNameIndex, outFileNameIndex, DBReader<std::string>::USE_INDEX);
        indexReader.open(DBReader<std::string>::SORT_BY_ID);
        DBReader<std::string>::Index *index = indexReader.getIndex();
        FILE *indexFile  = fopen(outFileNameIndex, "w");
        if (indexFile == NULL) {
            perror(outFileNameIndex);
            EXIT(EXIT_FAILURE);
        }
        writeIndex(indexFile, indexReader.getSize(), index, indexReader.getSeqLens());
        if (fsync(fileno(indexFile)) < 0) {
            Debug(Debug::ERROR) << "Could not sync index file " << outFileNameIndex << "!\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(indexFile);
        indexReader.close();
    }
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for merging files: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) <<" s\n";
}


void DBWriter::mergeFilePair(const std::vector<std::pair<std::string, std::string>> fileNames) {
    FILE ** files = new FILE*[fileNames.size()];
    for (size_t i = 0; i < fileNames.size();i++) {
        files[i] = FileUtil::openFileOrDie(fileNames[i].first.c_str(), "r", true);
#if HAVE_POSIX_FADVISE
        int status;
        if ((status = posix_fadvise (fileno(files[i]), 0, 0, POSIX_FADV_SEQUENTIAL)) != 0){
           Debug(Debug::ERROR) << "posix_fadvise returned an error: " << strerror(status) << "\n";
        }
#endif
    }

    int c1;
    char * buffer = dataFilesBuffer[0];
    size_t writePos = 0;
    int dataFilefd = fileno(dataFiles[0]);
    do {
        for (size_t i = 0; i < fileNames.size(); ++i) {
            while ((c1 = getc_unlocked(files[i])) != EOF) {
                if (c1 == '\0') {
                    break;
                }
                buffer[writePos] = (char) c1;
                writePos++;
                if (writePos == bufferSize) {
                    size_t written = write(dataFilefd, buffer, bufferSize);
                    if (written != bufferSize) {
                        Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[0] << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    writePos = 0;
                }
            }
        }
        buffer[writePos] = '\0';
        writePos++;
        if (writePos == bufferSize) {
            size_t written = write(dataFilefd, buffer, bufferSize);
            if (written != bufferSize) {
                Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[0] << "\n";
                EXIT(EXIT_FAILURE);
            }
            writePos = 0;
        }
    } while (c1!=EOF);

    if (writePos != 0) {
        // if there is data in the buffer that is not yet written
        size_t written = write(dataFilefd, (const void *) dataFilesBuffer[0], writePos);
        if (written != writePos) {
            Debug(Debug::ERROR) << "Could not write to data file " << dataFileNames[0] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    for (size_t i = 0; i < fileNames.size(); ++i) {
        if (fsync(fileno(files[i])) < 0) {
            Debug(Debug::ERROR) << "Could not sync data file " << fileNames[i].first << "!\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(files[i]);
    }

    Debug(Debug::INFO) << "Merge file " << fileNames[0].first << " and " << fileNames[0].second << "\n";
    DBReader<unsigned int> reader1(fileNames[0].first.c_str(), fileNames[0].second.c_str(),
                                   DBReader<unsigned int>::USE_INDEX);
    reader1.open(DBReader<unsigned int>::NOSORT);
    unsigned int *seqLen1 = reader1.getSeqLens();
    DBReader<unsigned int>::Index *index1 = reader1.getIndex();

    for (size_t i = 1; i < fileNames.size(); i++) {
        DBReader<unsigned int> reader2(fileNames[i].first.c_str(), fileNames[i].second.c_str(),
                                       DBReader<unsigned int>::USE_INDEX);
        reader2.open(DBReader<unsigned int>::NOSORT);
        unsigned int *seqLen2 = reader2.getSeqLens();
        size_t currOffset = 0;

        for (size_t id = 0; id < reader1.getSize(); id++) {
            // add length for file1 and file2 and subtract -1 for one null byte
            size_t seqLen = seqLen1[id] + seqLen2[id] - 1;
            seqLen1[id] = seqLen;
            index1[id].offset = currOffset;
            currOffset += seqLen;
        }
        reader2.close();
    }

    writeIndex(indexFiles[0], reader1.getSize(), index1, seqLen1);
    reader1.close();
}
