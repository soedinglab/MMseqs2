#include "DBWriter.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Concat.h"
#include "itoa.h"
#include "Timer.h"
#include "Parameters.h"

#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

DBWriter::DBWriter(const char *dataFileName_, const char *indexFileName_, unsigned int threads, size_t mode, int dbtype)
        : threads(threads), mode(mode), dbtype(dbtype) {
    dataFileName = strdup(dataFileName_);
    indexFileName = strdup(indexFileName_);

    dataFiles = new FILE *[threads];
    dataFilesBuffer = new char *[threads];

    dataFileNames = new char *[threads];

    indexFiles = new FILE *[threads];
    indexFileNames = new char *[threads];
    compressedBuffers=NULL;
    compressedBufferSizes=NULL;
    if((mode & Parameters::WRITER_COMPRESSED_MODE) != 0){
        compressedBuffers = new char*[threads];
        compressedBufferSizes = new size_t[threads];
        cstream = new ZSTD_CStream*[threads];
        state = new int[threads];
        threadBuffer = new char*[threads];
        threadBufferSize = new size_t[threads];
        threadBufferOffset = new size_t[threads];
    }

    starts = new size_t[threads];
    std::fill(starts, starts + threads, 0);
    offsets = new size_t[threads];
    std::fill(offsets, offsets + threads, 0);
    if((mode & Parameters::WRITER_COMPRESSED_MODE) != 0 ){
        datafileMode = "wb+";
    } else {
        datafileMode = "wb";
    }

    closed = true;
}

size_t DBWriter::addToThreadBuffer(const void *data, size_t itmesize, size_t nitems, int threadIdx) {
    size_t bytesToWrite = (itmesize*nitems);
    size_t bytesLeftInBuffer = threadBufferSize[threadIdx] - threadBufferOffset[threadIdx];
    if( (itmesize*nitems)  >= bytesLeftInBuffer ){
        size_t newBufferSize = std::max(threadBufferSize[threadIdx] + bytesToWrite, threadBufferSize[threadIdx] * 2 );
        threadBufferSize[threadIdx] = newBufferSize;
        threadBuffer[threadIdx] = (char*) realloc(threadBuffer[threadIdx], newBufferSize);
        if(compressedBuffers[threadIdx] == NULL){
            Debug(Debug::ERROR) << "Realloc of buffer for " << threadIdx << " failed. Buffer size = " << threadBufferSize[threadIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    memcpy(threadBuffer[threadIdx] + threadBufferOffset[threadIdx], data, bytesToWrite);
    threadBufferOffset[threadIdx] += bytesToWrite;
    return bytesToWrite;
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
    if(compressedBuffers){
        delete [] threadBuffer;
        delete [] threadBufferSize;
        delete [] threadBufferOffset;
        delete [] compressedBuffers;
        delete [] compressedBufferSizes;
        delete [] cstream;
        delete [] state;
    }
}

void DBWriter::sortDatafileByIdOrder(DBReader<unsigned int> &dbr) {
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

#pragma omp for schedule(static)
        for (size_t id = 0; id < dbr.getSize(); id++) {
            char *data = dbr.getData(id, thread_idx);
            size_t length = dbr.getEntryLen(id);
            writeData(data, (length == 0 ? 0 : length - 1), dbr.getDbKey(id), thread_idx);
        }
    };

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
                                                     files[i].second.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        filesToMerge[i]->open(DBReader<unsigned int>::NOSORT);
    }
    std::string result;

    for (size_t id = 0; id < qdbr.getSize(); id++) {
        unsigned int key = qdbr.getDbKey(id);
        // get all data for the id from all files
        for (size_t i = 0; i < fileCount; i++) {
            const char *data = filesToMerge[i]->getDataByDBKey(key, 0);
            if (data != NULL) {
                if(i < prefixes.size()) {
                    result.append( prefixes[i]);
                }
                result.append(data);
            }
        }
        // write result
        writeData(result.c_str(), result.length(), key, 0);
        result.clear();
    }

    // close all reader
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i]->close();
        delete filesToMerge[i];
    }
    delete [] filesToMerge;


}

// allocates heap memory, careful
char* makeResultFilename(const char* name, size_t split) {
    std::ostringstream ss;
    ss << name << "." << split;
    std::string s = ss.str();
    return strdup(s.c_str());
}

void DBWriter::open(size_t bufferSize) {
    if (bufferSize == SIZE_MAX) {
        if (Util::getTotalSystemMemory() < (8ull * 1024 * 1024 * 1024)) {
            // reduce this buffer if our system does not have much memory
            // createdb runs into trouble since it creates 2x32 splits with 64MB each (=4GB)
            // 8MB should be enough
            bufferSize = 8ull * 1024 * 1024;
        } else {
            bufferSize = 32ull * 1024 * 1024;
        }
    }
    for (unsigned int i = 0; i < threads; i++) {
        dataFileNames[i] = makeResultFilename(dataFileName, i);
        indexFileNames[i] = makeResultFilename(indexFileName, i);

        dataFiles[i] = FileUtil::openAndDelete(dataFileNames[i], datafileMode.c_str());
        int fd = fileno(dataFiles[i]);
        int flags;
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Can not set mode for " << dataFileNames[i] << "!\n";
            EXIT(EXIT_FAILURE);
        }

        dataFilesBuffer[i] = new(std::nothrow) char[bufferSize];
        incrementMemory(bufferSize);
        Util::checkAllocation(dataFilesBuffer[i], "Cannot allocate buffer for DBWriter");
        this->bufferSize = bufferSize;

        // set buffer to 64
        if (setvbuf(dataFiles[i], dataFilesBuffer[i], _IOFBF, bufferSize) != 0) {
            Debug(Debug::WARNING) << "Write buffer could not be allocated (bufferSize=" << bufferSize << ")\n";
        }

        indexFiles[i] = FileUtil::openAndDelete(indexFileNames[i], "w");
        fd = fileno(indexFiles[i]);
        if ((flags = fcntl(fd, F_GETFL, 0)) < 0 || fcntl(fd, F_SETFD, flags | FD_CLOEXEC) == -1) {
            Debug(Debug::ERROR) << "Can not set mode for " << indexFileNames[i] << "!\n";
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

        if((mode & Parameters::WRITER_COMPRESSED_MODE) != 0){
            compressedBufferSizes[i] = 2097152;
            threadBufferSize[i] = 2097152;
            state[i] = false;
            compressedBuffers[i] = (char*) malloc(compressedBufferSizes[i]);
            incrementMemory(compressedBufferSizes[i]);
            threadBuffer[i] = (char*) malloc(threadBufferSize[i]);
            incrementMemory(threadBufferSize[i]);
            cstream[i] = ZSTD_createCStream();
        }
    }

    closed = false;
}

void DBWriter::writeDbtypeFile(const char* path, int dbtype, bool isCompressed) {
    if (dbtype == Parameters::DBTYPE_OMIT_FILE) {
        return;
    }

    std::string name = std::string(path) + ".dbtype";
    FILE* file = FileUtil::openAndDelete(name.c_str(), "wb");
    dbtype = isCompressed ? dbtype | (1 << 31) : dbtype & ~(1 << 31);
    size_t written = fwrite(&dbtype, sizeof(int), 1, file);
    if (written != 1) {
        Debug(Debug::ERROR) << "Can not write to data file " << name << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(file);
}


void DBWriter::close(bool merge) {
    // close all datafiles
    for (unsigned int i = 0; i < threads; i++) {
        fclose(dataFiles[i]);
        fclose(indexFiles[i]);
    }

    if(compressedBuffers){
        for (unsigned int i = 0; i < threads; i++) {
            free(compressedBuffers[i]);
            decrementMemory(compressedBufferSizes[i]);
            free(threadBuffer[i]);
            decrementMemory(threadBufferSize[i]);
            ZSTD_freeCStream(cstream[i]);
        }
    }

    mergeResults(dataFileName, indexFileName, (const char **) dataFileNames, (const char **) indexFileNames,
                 threads, merge, ((mode & Parameters::WRITER_LEXICOGRAPHIC_MODE) != 0));

    writeDbtypeFile(dataFileName, dbtype, (mode & Parameters::WRITER_COMPRESSED_MODE) != 0);

    for (unsigned int i = 0; i < threads; i++) {
        delete [] dataFilesBuffer[i];
        decrementMemory(bufferSize);
        free(dataFileNames[i]);
        free(indexFileNames[i]);
    }
    closed = true;
}

void DBWriter::writeStart(unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }
    starts[thrIdx] = offsets[thrIdx];
    if((mode & Parameters::WRITER_COMPRESSED_MODE) != 0){
        state[thrIdx] = INIT_STATE;
        threadBufferOffset[thrIdx]=0;
        int cLevel = 3;
        size_t const initResult = ZSTD_initCStream(cstream[thrIdx], cLevel);
        if (ZSTD_isError(initResult)) {
            Debug(Debug::ERROR) << "ZSTD_initCStream() error in thread " << thrIdx << ". Error "
                                << ZSTD_getErrorName(initResult) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

size_t DBWriter::writeAdd(const char* data, size_t dataSize, unsigned int thrIdx) {
    checkClosed();
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "Thread index " << thrIdx << " > maximum thread number " << threads << "\n";
        EXIT(EXIT_FAILURE);
    }
    bool isCompressedDB = (mode & Parameters::WRITER_COMPRESSED_MODE) != 0;
    if(isCompressedDB && state[thrIdx] == INIT_STATE && dataSize < 60){
        state[thrIdx] = NOTCOMPRESSED;
    }
    size_t totalWriten = 0;
    if(isCompressedDB && (state[thrIdx] == INIT_STATE || state[thrIdx] == COMPRESSED) ) {
        state[thrIdx] = COMPRESSED;
        // zstd seems to have a hard time with elements < 60
        ZSTD_inBuffer input = { data, dataSize, 0 };
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {compressedBuffers[thrIdx], compressedBufferSizes[thrIdx], 0};
            size_t toRead = ZSTD_compressStream( cstream[thrIdx], &output, &input);   /* toRead is guaranteed to be <= ZSTD_CStreamInSize() */
            if (ZSTD_isError(toRead)) {
                Debug(Debug::ERROR) << "ZSTD_compressStream() error in thread " << thrIdx << ". Error "
                                    << ZSTD_getErrorName(toRead) << "\n";
                EXIT(EXIT_FAILURE);
            }
            size_t written = addToThreadBuffer(compressedBuffers[thrIdx], sizeof(char), output.pos, thrIdx);
            if (written != output.pos) {
                Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            offsets[thrIdx] += written;
            totalWriten += written;
        }
    }else{
        size_t written;
        if(isCompressedDB){
            written = addToThreadBuffer(data, sizeof(char), dataSize,  thrIdx);
        }else{
            written = fwrite(data, sizeof(char), dataSize, dataFiles[thrIdx]);
        }
        if (written != dataSize) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        offsets[thrIdx] += written;
    }

    return totalWriten;
}

void DBWriter::writeEnd(unsigned int key, unsigned int thrIdx, bool addNullByte, bool addIndexEntry) {
    // close stream
    bool isCompressedDB = (mode & Parameters::WRITER_COMPRESSED_MODE) != 0;
    if(isCompressedDB) {
        size_t compressedLength = 0;
        if(state[thrIdx] == COMPRESSED) {
            ZSTD_outBuffer output = {compressedBuffers[thrIdx], compressedBufferSizes[thrIdx], 0};
            size_t remainingToFlush = ZSTD_endStream(cstream[thrIdx], &output); /* close frame */

            //        std::cout << compressedLength << std::endl;
            if (ZSTD_isError(remainingToFlush)) {
                Debug(Debug::ERROR) << "ZSTD_endStream() error in thread " << thrIdx << ". Error "
                                    << ZSTD_getErrorName(remainingToFlush) << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (remainingToFlush) {
                Debug(Debug::ERROR) << "Stream not flushed\n";
                EXIT(EXIT_FAILURE);
            }
            size_t written = addToThreadBuffer(compressedBuffers[thrIdx], sizeof(char), output.pos, thrIdx);
            compressedLength = threadBufferOffset[thrIdx];
            offsets[thrIdx] += written;
            if (written != output.pos) {
                Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }else {
            compressedLength = offsets[thrIdx] - starts[thrIdx];
        }
        unsigned int compressedLengthInt = static_cast<unsigned int>(compressedLength);
        size_t written2 = fwrite(&compressedLengthInt, sizeof(unsigned int), 1, dataFiles[thrIdx]);
        if (written2 != 1) {
            Debug(Debug::ERROR) << "Can not write entry length to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        offsets[thrIdx] +=  sizeof(unsigned int);
        writeThreadBuffer(thrIdx, compressedLength);
    }


    size_t totalWritten = 0;
// entries are always separated by a null byte
    if (addNullByte == true) {
        char nullByte = '\0';
        if(isCompressedDB && state[thrIdx]==NOTCOMPRESSED){
            nullByte = static_cast<char>(0xFF);
        }
        const size_t written = fwrite(&nullByte, sizeof(char), 1, dataFiles[thrIdx]);
        if (written != 1) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        totalWritten += written;
        offsets[thrIdx] += 1;
    }

    if (addIndexEntry == true) {
        size_t length = offsets[thrIdx] - starts[thrIdx];
// keep original size in index
        if (isCompressedDB && state[thrIdx]==COMPRESSED) {
            ZSTD_frameProgression progression = ZSTD_getFrameProgression(cstream[thrIdx]);
            length = progression.consumed + totalWritten;
        }
        if (isCompressedDB && state[thrIdx]==NOTCOMPRESSED) {
            length -= sizeof(unsigned int);
        }
        writeIndexEntry(key, starts[thrIdx], length, thrIdx);
    }
}

void DBWriter::writeIndexEntry(unsigned int key, size_t offset, size_t length, unsigned int thrIdx){
    char buffer[1024];
    size_t len = indexToBuffer(buffer, key, offset, length );
    size_t written = fwrite(buffer, sizeof(char), len, indexFiles[thrIdx]);
    if (written != len) {
        Debug(Debug::ERROR) << "Can not write to data file " << dataFileName[thrIdx] << "\n";
        EXIT(EXIT_FAILURE);
    }
}


void DBWriter::writeData(const char *data, size_t dataSize, unsigned int key, unsigned int thrIdx, bool addNullByte, bool addIndexEntry) {
    writeStart(thrIdx);
    writeAdd(data, dataSize, thrIdx);
    writeEnd(key, thrIdx, addNullByte, addIndexEntry);
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

void DBWriter::alignToPageSize(int thrIdx) {
    size_t currentOffset = offsets[thrIdx];
    size_t pageSize = Util::getPageSize();
    size_t newOffset = ((pageSize - 1) & currentOffset) ? ((currentOffset + pageSize) & ~(pageSize - 1)) : currentOffset;
    char nullByte = '\0';
    for (size_t i = currentOffset; i < newOffset; ++i) {
        size_t written = fwrite(&nullByte, sizeof(char), 1, dataFiles[thrIdx]);
        if (written != 1) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileNames[thrIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    offsets[thrIdx] = newOffset;
}


void DBWriter::checkClosed() {
    if (closed == true) {
        Debug(Debug::ERROR) << "Trying to read a closed database. Datafile=" << dataFileName  << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::mergeResults(const std::string &outFileName, const std::string &outFileNameIndex,
                            const std::vector<std::pair<std::string, std::string >> &files,
                            const bool lexicographicOrder) {
    const char **datafilesNames = new const char *[files.size()];
    const char **indexFilesNames = new const char *[files.size()];
    for (size_t i = 0; i < files.size(); i++) {
        datafilesNames[i] = files[i].first.c_str();
        indexFilesNames[i] = files[i].second.c_str();
    }
    mergeResults(outFileName.c_str(), outFileNameIndex.c_str(), datafilesNames, indexFilesNames, files.size(), true, lexicographicOrder);
    delete[] datafilesNames;
    delete[] indexFilesNames;

    // leave only one dbtype file behind
    if (files.size() > 0) {
        std::string typeSrc = files[0].first + ".dbtype";
        std::string typeDest = outFileName + ".dbtype";
        if (FileUtil::fileExists(typeSrc.c_str())) {
            std::rename(typeSrc.c_str(), typeDest.c_str());
        }
        for (size_t i = 1; i < files.size(); i++) {
            std::string typeFile = files[i].first + ".dbtype";
            if (FileUtil::fileExists(typeFile.c_str())) {
                FileUtil::remove(typeFile.c_str());
            }
        }
    }
}

template <>
void DBWriter::writeIndexEntryToFile(FILE *outFile, char *buff1, DBReader<unsigned int>::Index &index){
    char * tmpBuff = Itoa::u32toa_sse2((uint32_t)index.id,buff1);
    *(tmpBuff-1) = '\t';
    size_t currOffset = index.offset;
    tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
    *(tmpBuff-1) = '\t';
    uint32_t sLen = index.length;
    tmpBuff = Itoa::u32toa_sse2(sLen,tmpBuff);
    *(tmpBuff-1) = '\n';
    *(tmpBuff) = '\0';
    fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
}

template <>
void DBWriter::writeIndexEntryToFile(FILE *outFile, char *buff1, DBReader<std::string>::Index &index)
{
    size_t keyLen = index.id.length();
    char * tmpBuff = (char*)memcpy((void*)buff1, (void*)index.id.c_str(), keyLen);
    tmpBuff+=keyLen;
    *(tmpBuff) = '\t';
    tmpBuff++;
    size_t currOffset = index.offset;
    tmpBuff = Itoa::u64toa_sse2(currOffset, tmpBuff);
    *(tmpBuff-1) = '\t';
    uint32_t sLen = index.length;
    tmpBuff = Itoa::u32toa_sse2(sLen,tmpBuff);
    *(tmpBuff-1) = '\n';
    *(tmpBuff) = '\0';
    fwrite(buff1, sizeof(char), (tmpBuff - buff1), outFile);
}

template <>
void DBWriter::writeIndex(FILE *outFile, size_t indexSize, DBReader<unsigned int>::Index *index) {
    char buff1[1024];
    for (size_t id = 0; id < indexSize; id++) {
        writeIndexEntryToFile(outFile, buff1, index[id]);
    }
}

template <>
void DBWriter::writeIndex(FILE *outFile, size_t indexSize, DBReader<std::string>::Index *index){
    char buff1[1024];
    for (size_t id = 0; id < indexSize; id++) {
        writeIndexEntryToFile(outFile, buff1, index[id]);
    }
}


void DBWriter::mergeResults(const char *outFileName, const char *outFileNameIndex,
                            const char **dataFileNames, const char **indexFileNames,
                            unsigned long fileCount, bool mergeDatafiles, bool lexicographicOrder) {
    Timer timer;
    std::vector<std::vector<std::string>> dataFilenames;
    for (unsigned int i = 0; i < fileCount; ++i) {
        dataFilenames.emplace_back(FileUtil::findDatafiles(dataFileNames[i]));
    }

    // merge results into one result file
    if (dataFilenames.size() > 1) {
        std::vector<FILE*> datafiles;
        std::vector<size_t> mergedSizes;
        for (unsigned int i = 0; i < dataFilenames.size(); i++) {
            std::vector<std::string>& filenames = dataFilenames[i];
            size_t cumulativeSize = 0;
            for (size_t j = 0; j < filenames.size(); ++j) {
                FILE* fh = fopen(filenames[j].c_str(), "r");
                if (fh == NULL) {
                    Debug(Debug::ERROR) << "Can not open result file " << filenames[j] << "!\n";
                    EXIT(EXIT_FAILURE);
                }
                struct stat sb;
                if (fstat(fileno(fh), &sb) < 0) {
                    int errsv = errno;
                    Debug(Debug::ERROR) << "Failed to fstat file " << filenames[j] << ". Error " << errsv << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                datafiles.emplace_back(fh);
                cumulativeSize += sb.st_size;
            }
            mergedSizes.push_back(cumulativeSize);
        }

        if (mergeDatafiles) {
            FILE *outFh = FileUtil::openAndDelete(outFileName, "w");
            Concat::concatFiles(datafiles, outFh);
            fclose(outFh);
        }

        for (unsigned int i = 0; i < datafiles.size(); ++i) {
            fclose(datafiles[i]);
        }

        if (mergeDatafiles) {
            for (unsigned int i = 0; i < dataFilenames.size(); i++) {
                std::vector<std::string>& filenames = dataFilenames[i];
                for (size_t j = 0; j < filenames.size(); ++j) {
                    FileUtil::remove(filenames[j].c_str());
                }
            }
        }

        // merge index
        mergeIndex(indexFileNames, dataFilenames.size(), mergedSizes);
    } else {
        std::vector<std::string>& filenames = dataFilenames[0];
        if (filenames.size() == 1) {
            // In single thread dbreader mode it will create a .0
            // that should be moved to the final destination dest instead of dest.0
            FileUtil::move(filenames[0].c_str(), outFileName);
        } else {
            DBReader<unsigned int>::moveDatafiles(filenames, outFileName);
        }
    }

    DBWriter::sortIndex(indexFileNames[0], outFileNameIndex, lexicographicOrder);
    FileUtil::remove(indexFileNames[0]);
    Debug(Debug::INFO) << "Time for merging to " << FileUtil::baseName(outFileName) << ": " << timer.lap() << "\n";
}

void DBWriter::mergeIndex(const char** indexFilenames, unsigned int fileCount, const std::vector<size_t> &dataSizes) {
    FILE *index_file = fopen(indexFilenames[0], "a");
    if (index_file == NULL) {
        perror(indexFilenames[0]);
        EXIT(EXIT_FAILURE);
    }
    size_t globalOffset = dataSizes[0];
    for (unsigned int fileIdx = 1; fileIdx < fileCount; fileIdx++) {
        DBReader<unsigned int> reader(indexFilenames[fileIdx], indexFilenames[fileIdx], 1, DBReader<unsigned int>::USE_INDEX);
        reader.open(DBReader<unsigned int>::HARDNOSORT);
        if (reader.getSize() > 0) {
            DBReader<unsigned int>::Index * index = reader.getIndex();
            for (size_t i = 0; i < reader.getSize(); i++) {
                size_t currOffset = index[i].offset;
                index[i].offset = globalOffset + currOffset;
            }
            writeIndex(index_file, reader.getSize(), index);
        }
        reader.close();
        FileUtil::remove(indexFilenames[fileIdx]);

        globalOffset += dataSizes[fileIdx];
    }
    fclose(index_file);
}

void DBWriter::sortIndex(const char *inFileNameIndex, const char *outFileNameIndex, const bool lexicographicOrder){
    if (lexicographicOrder == false) {
        // sort the index
        DBReader<unsigned int> indexReader(inFileNameIndex, inFileNameIndex, 1, DBReader<unsigned int>::USE_INDEX);
        indexReader.open(DBReader<unsigned int>::NOSORT);
        DBReader<unsigned int>::Index *index = indexReader.getIndex();
        FILE *index_file  = FileUtil::openAndDelete(outFileNameIndex, "w");
        writeIndex(index_file, indexReader.getSize(), index);
        fclose(index_file);
        indexReader.close();

    } else {
        DBReader<std::string> indexReader(inFileNameIndex, inFileNameIndex, 1, DBReader<std::string>::USE_INDEX);
        indexReader.open(DBReader<std::string>::SORT_BY_ID);
        DBReader<std::string>::Index *index = indexReader.getIndex();
        FILE *index_file  = FileUtil::openAndDelete(outFileNameIndex, "w");
        writeIndex(index_file, indexReader.getSize(), index);
        fclose(index_file);
        indexReader.close();
    }
}

void DBWriter::writeThreadBuffer(unsigned int idx, size_t dataSize) {
    size_t written = fwrite(threadBuffer[idx], 1, dataSize, dataFiles[idx]);
    if (written != dataSize) {
        Debug(Debug::ERROR) << "writeThreadBuffer: Could not write to data file " << dataFileNames[idx] << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void DBWriter::createRenumberedDB(const std::string& dataFile, const std::string& indexFile, const std::string& origData, const std::string& origIndex, int sortMode) {
    DBReader<unsigned int>* lookupReader = NULL;
    FILE *sLookup = NULL;
    if (origData.empty() == false && origIndex.empty() == false) {
        lookupReader = new DBReader<unsigned int>(origData.c_str(), origIndex.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
        lookupReader->open(DBReader<unsigned int>::NOSORT);
        sLookup = FileUtil::openAndDelete((dataFile + ".lookup").c_str(), "w");
    }

    DBReader<unsigned int> reader(dataFile.c_str(), indexFile.c_str(), 1, DBReader<unsigned int>::USE_INDEX);
    reader.open(sortMode);
    std::string indexTmp = indexFile + "_tmp";
    FILE *sIndex = FileUtil::openAndDelete(indexTmp.c_str(), "w");

    char buffer[1024];
    std::string strBuffer;
    strBuffer.reserve(1024);
    DBReader<unsigned int>::LookupEntry* lookup = NULL;
    if (lookupReader != NULL) {
        lookup = lookupReader->getLookup();
    }
    for (size_t i = 0; i < reader.getSize(); i++) {
        DBReader<unsigned int>::Index *idx = (reader.getIndex(i));
        size_t len = DBWriter::indexToBuffer(buffer, i, idx->offset, idx->length);
        int written = fwrite(buffer, sizeof(char), len, sIndex);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Can not write to data file " << indexFile << "_tmp\n";
            EXIT(EXIT_FAILURE);
        }
        if (lookupReader != NULL) {
            size_t lookupId = lookupReader->getLookupIdByKey(idx->id);
            DBReader<unsigned int>::LookupEntry copy = lookup[lookupId];
            copy.id = i;
            copy.entryName = SSTR(idx->id);
            lookupReader->lookupEntryToBuffer(strBuffer, copy);
            written = fwrite(strBuffer.c_str(), sizeof(char), strBuffer.size(), sLookup);
            if (written != (int) strBuffer.size()) {
                Debug(Debug::ERROR) << "Could not write to lookup file " << indexFile << "_tmp\n";
                EXIT(EXIT_FAILURE);
            }
            strBuffer.clear();
        }
    }
    fclose(sIndex);
    reader.close();
    std::rename(indexTmp.c_str(), indexFile.c_str());

    if (lookupReader != NULL) {
        fclose(sLookup);
        lookupReader->close();
    }
}
