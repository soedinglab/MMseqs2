#include "DBReader.h"
#include "FastSort.h"
#include <algorithm>
#include <climits>
#include <cstring>
#include <cstddef>
#include <random>

#include <sys/mman.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <unistd.h>

#include "MemoryMapped.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "itoa.h"
#include "BatchReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

template <typename T>
DBReader<T>::DBReader(const char* dataFileName_, const char* indexFileName_, int threads, int dataMode) :
threads(threads), dataMode(dataMode), dataFileName(strdup(dataFileName_)),
        indexFileName(strdup(indexFileName_)), size(0), dataFiles(NULL), dataSizeOffset(NULL), dataFileCnt(0),
        totalDataSize(0), dataSize(0), lastKey(T()), closed(1), dbtype(Parameters::DBTYPE_GENERIC_DB),
        compressedBuffers(NULL), compressedBufferSizes(NULL), indexSize(0), index(NULL), id2local(NULL), local2id(NULL),
        index2dbfile(NULL), clusterIndex(false), batchReader(NULL), loadedBuffer(), prefetchBuffer(), dataMapped(false), accessType(0), externalData(false), didMlock(false)
{}

template <typename T>
DBReader<T>::DBReader(DBReader<T>::Index *index, size_t size, size_t dataSize, T lastKey,
        int dbType, size_t maxSeqLen, int threads) :
        threads(threads), dataMode(USE_INDEX), dataFileName(NULL), indexFileName(NULL),
        size(size), dataFiles(NULL), dataSizeOffset(NULL), dataFileCnt(0), totalDataSize(0), dataSize(dataSize), lastKey(lastKey),
        maxSeqLen(maxSeqLen), closed(1), dbtype(dbType), compressedBuffers(NULL), compressedBufferSizes(NULL), indexSize(size), index(index), sortedByOffset(true),
        id2local(NULL), local2id(NULL), index2dbfile(NULL), clusterIndex(false), batchReader(NULL), loadedBuffer(), prefetchBuffer(), dataMapped(false), accessType(NOSORT), externalData(true), didMlock(false)
{}

template <typename T>
void DBReader<T>::setDataFile(const char* dataFileName_)  {
    if (dataFileName != NULL) {
        unmapData();
        free(dataFileName);
    }

    dataMode |= USE_DATA;
    dataFileName = strdup(dataFileName_);
}

template <typename T>
void DBReader<T>::readMmapedDataInMemory(){
    if ((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
        //Debug(Debug::INFO) << "Touch data file " << dataFileName << "\n";
        for(size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++){
            size_t dataSize = dataSizeOffset[fileIdx+1]-dataSizeOffset[fileIdx];
            magicBytes += Util::touchMemory(dataFiles[fileIdx], dataSize);
        }

    }
}

template <typename T>
void DBReader<T>::mlock(){
    if (dataMode & USE_DATA) {
        if (didMlock == false) {
            for(size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++) {
                size_t dataSize = dataSizeOffset[fileIdx+1]-dataSizeOffset[fileIdx];
                ::mlock(dataFiles[fileIdx], dataSize);
            }
        }
        didMlock = true;
    }
}

template <typename T>
void DBReader<T>::printMagicNumber(){
    Debug(Debug::INFO) << magicBytes << "\n";
}

template <typename T> DBReader<T>::~DBReader(){
    if(dataFileName != NULL) {
        free(dataFileName);
    }

    if(indexFileName != NULL) {
        free(indexFileName);
    }

    if(dataSizeOffset != NULL){
        delete [] dataSizeOffset;
    }

    if(dataFiles != NULL){
        delete [] dataFiles;
    }
}

template <typename T> bool DBReader<T>::open(int accessType){
    // count the number of entries
    this->accessType = accessType;
    bool isLinfinity = false;
    if (dataFileName != NULL) {
        // the data file <clu>.<node> of a --createdb-mode 3 cluster database has the dbtype of <clu>
        std::string dbtypeName(dataFileName);
        const size_t nodePos = dbtypeName.find_last_of('.');
        if ((dataMode & USE_LINFINITY) && nodePos != std::string::npos && Util::isNumber(dbtypeName.substr(nodePos + 1))
            && FileUtil::fileExists((dbtypeName + ".dbtype").c_str()) == false) {
            dbtypeName.resize(nodePos);
        }
        dbtype = FileUtil::parseDbType(dbtypeName.c_str());
        isLinfinity = getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_LINFINITY;
        clusterIndex = isLinfinity && Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_CLUSTER_RES);
        if (isLinfinity && (dataMode & USE_LINFINITY) == 0) {
            Debug(Debug::ERROR) << "Database " << dataFileName << " was created with --createdb-mode 3, which this command does not support\n";
            EXIT(EXIT_FAILURE);
        }
        // a mode 3 database is only read in key order
        if (isLinfinity && accessType != NOSORT && accessType != HARDNOSORT && accessType != LINEAR_ACCCESS) {
            accessType = NOSORT;
            this->accessType = NOSORT;
        }
    }
    if (dataMode & USE_DATA) {
        dataFileNames = FileUtil::findDatafiles(dataFileName);
        if (dataFileNames.empty()) {
            Debug(Debug::ERROR) << "No datafile could be found for " << dataFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        totalDataSize = 0;
        dataFileCnt = dataFileNames.size();
        dataSizeOffset = new size_t[dataFileNames.size() + 1];
        dataFiles = new char*[dataFileNames.size()];
        for(size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++){
            if (isLinfinity) {
                dataFiles[fileIdx] = NULL;
                dataSizeOffset[fileIdx] = 0;
                continue;
            }
            FILE* dataFile = fopen(dataFileNames[fileIdx].c_str(), "r");
            if (dataFile == NULL) {
                Debug(Debug::ERROR) << "Cannot open data file " << dataFileName << "!\n";
                EXIT(EXIT_FAILURE);
            }
            size_t dataSize;
            dataFiles[fileIdx] = mmapData(dataFile, &dataSize);
            dataSizeOffset[fileIdx]=totalDataSize;
            totalDataSize += dataSize;
            if (fclose(dataFile) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << dataFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        dataSizeOffset[dataFileNames.size()]=totalDataSize;
        dataMapped = true;
        if (accessType == LINEAR_ACCCESS || accessType == SORT_BY_OFFSET) {
            setSequentialAdvice();
        }
    }
    if (dataMode & USE_LOOKUP || dataMode & USE_LOOKUP_REV) {
        std::string lookupFilename = (std::string(dataFileName) + ".lookup");
        MemoryMapped lookupData(lookupFilename, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (lookupData.isValid() == false) {
            Debug(Debug::ERROR) << "Cannot open lookup file " << lookupFilename << "!\n";
            EXIT(EXIT_FAILURE);
        }
        char* lookupDataChar = (char *) lookupData.getData();
        size_t lookupDataSize = lookupData.size();
        lookupSize = Util::ompCountLines(lookupDataChar, lookupDataSize, threads);
        lookup = new(std::nothrow) LookupEntry[this->lookupSize];
        incrementMemory(sizeof(LookupEntry) * this->lookupSize);
        readLookup(lookupDataChar, lookupDataSize, lookup);
        if (dataMode & USE_LOOKUP) {
            SORT_PARALLEL(lookup, lookup + lookupSize, LookupEntry::compareById);
        } else {
            SORT_PARALLEL(lookup, lookup + lookupSize, LookupEntry::compareByAccession);
        }
        lookupData.close();
    }

    if (dataMode & USE_SOURCE || dataMode & USE_SOURCE_REV) {
        std::string sourceFilename = (std::string(dataFileName) + ".source");
        MemoryMapped sourceData(sourceFilename, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (sourceData.isValid() == false) {
            Debug(Debug::ERROR) << "Cannot open source file " << sourceFilename << "!\n";
            EXIT(EXIT_FAILURE);
        }
        char* sourceDataChar = (char *) sourceData.getData();
        size_t sourceDataSize = sourceData.size();
        sourceSize = Util::ompCountLines(sourceDataChar, sourceDataSize, threads);
        source = new(std::nothrow) SourceEntry[this->sourceSize];
        incrementMemory(sizeof(SourceEntry) * this->sourceSize);
        readSource(sourceDataChar, sourceDataSize, source);
        if (dataMode & USE_SOURCE) {
            SORT_SERIAL(source, source + sourceSize, SourceEntry::compareById);
        } else {
            SORT_SERIAL(source, source + sourceSize, SourceEntry::compareByFileName);
        }
        sourceData.close();
    }

    bool isSortedById = false;
    std::vector<size_t> fileBytes;
    if (externalData == false && clusterIndex) {
        // one range of all entries; only the entry count and the data size are read here
        FILE *indexFile = FileUtil::openFileOrDie(indexFileName, "r", true);
        uint64_t lastOffset = 0;
        if (fseek(indexFile, -static_cast<long>(sizeof(uint64_t)), SEEK_END) != 0 || fread(&lastOffset, sizeof(uint64_t), 1, indexFile) != 1) {
            Debug(Debug::ERROR) << "Cannot read index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t indexBytes = ftell(indexFile);
        if (fclose(indexFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        size = indexBytes / (2 * sizeof(uint64_t)) - 1;
        dataSize = lastOffset;
        maxSeqLen = 0;
        lastKey = T();
        sortedByOffset = true;
        indexSize = 1;
        index = new(std::nothrow) Index[indexSize];
        Util::checkAllocation(index, "Cannot allocate index memory in DBReader");
        incrementMemory(sizeof(Index) * indexSize);
        index[0].id = T();
        index[0].offset = 0;
        index[0].length = 0;
        index2dbfile = new(std::nothrow) unsigned int[indexSize];
        Util::checkAllocation(index2dbfile, "Cannot allocate index memory in DBReader");
        incrementMemory(indexSize * sizeof(unsigned int));
        index2dbfile[0] = 0;
        headerOffsetsFileNames.push_back(indexFileName);
        fileBytes.push_back(dataSize);
        fileBytes.push_back(indexBytes);
    } else if (externalData == false) {
        MemoryMapped indexData(indexFileName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!indexData.isValid()){
            Debug(Debug::ERROR) << "Cannot open index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        char* indexDataChar = (char *) indexData.getData();
        size_t indexDataSize = indexData.size();
        size = Util::ompCountLines(indexDataChar, indexDataSize, threads);

        indexSize = size;
        index = new(std::nothrow) Index[indexSize];
        Util::checkAllocation(index, "Cannot allocate index memory in DBReader");
        incrementMemory(sizeof(Index) * indexSize);

        if (isLinfinity) {
            // the headers are opened with <db>.nodeinfo as index, one line "firstKey count" per file
            const bool isHeader = Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_GENERIC_DB);
            if (isCompressed(dbtype) == COMPRESSED && isHeader == false) {
                Debug(Debug::ERROR) << "Database " << dataFileName << " was created with --createdb-mode 3 and cannot be compressed\n";
                EXIT(EXIT_FAILURE);
            }
            index2dbfile = new(std::nothrow) unsigned int[indexSize];
            Util::checkAllocation(index2dbfile, "Cannot allocate index memory in DBReader");
            incrementMemory(indexSize * sizeof(unsigned int));

            const char *cols[5];
            fileBytes.assign(dataFileNames.size() + (isHeader ? indexSize : 0), 0);
            size = 0;
            dataSize = 0;
            maxSeqLen = 0;
            sortedByOffset = false;
            for (size_t i = 0; i < indexSize; i++) {
                Util::getWordsOfLine(indexDataChar, cols, 5);
                readIndexId(&index[i].id, indexDataChar, cols);
                index[i].offset = isHeader ? 0 : Util::fast_atoi<size_t>(cols[1]);
                index[i].length = isHeader ? 0 : Util::fast_atoi<size_t>(cols[2]);
                const size_t count = Util::fast_atoi<size_t>(cols[isHeader ? 1 : 3]);
                index2dbfile[i] = isHeader ? static_cast<unsigned int>(i) : Util::fast_atoi<unsigned int>(cols[4]);
                if (isHeader) {
                    headerOffsetsFileNames.push_back(std::string(dataFileName) + ".index." + SSTR(i));
                    fileBytes[dataFileNames.size() + i] = (count + 1) * sizeof(uint64_t);
                }
                if ((count == 0 && isHeader == false) || indexIdToNum(&index[i].id) != size) {
                    Debug(Debug::ERROR) << "Keys in index file " << indexFileName << " are not consecutive at line " << i << "\n";
                    EXIT(EXIT_FAILURE);
                }
                if ((dataMode & USE_DATA) && index2dbfile[i] >= dataFileNames.size()) {
                    Debug(Debug::ERROR) << "Index file " << indexFileName << " refers to missing data file " << dataFileName << "." << index2dbfile[i] << "\n";
                    EXIT(EXIT_FAILURE);
                }
                if (dataMode & USE_DATA) {
                    fileBytes[index2dbfile[i]] += count * index[i].length;
                }
                size += count;
                dataSize += count * index[i].length;
                maxSeqLen = std::max(maxSeqLen, index[i].length);
                indexDataChar = Util::skipLine(indexDataChar);
            }
            lastKey = (size == 0) ? 0 : size - 1;
        } else {
            bool isSortedById = readIndex(indexDataChar, indexDataSize, index, dataSize);

            // sortIndex also handles access modes that don't require sorting
            sortIndex(isSortedById);

            size_t prevOffset = 0; // makes 0 or empty string
            sortedByOffset = true;
            for (size_t i = 0; i < size; i++) {
                sortedByOffset = sortedByOffset && index[i].offset >= prevOffset;
                prevOffset = index[i].offset;
            }
        }
        indexData.close();
    }
    if (externalData == false && isLinfinity && (dataMode & USE_DATA)) {
        const size_t pageCacheLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.1);
        for (size_t fileIdx = 0; fileIdx < dataFileNames.size() + headerOffsetsFileNames.size(); fileIdx++) {
            const std::string &fileName = (fileIdx < dataFileNames.size()) ? dataFileNames[fileIdx] : headerOffsetsFileNames[fileIdx - dataFileNames.size()];
            int fd = -1;
#ifdef O_DIRECT
            if ((dataMode & USE_DIRECT_IO) || ((dataMode & USE_BUFFERED_IO) == 0 && fileBytes[fileIdx] > pageCacheLimit)) {
                fd = ::open(fileName.c_str(), O_RDONLY | O_CLOEXEC | O_DIRECT);
            }
#endif
            if (fd < 0) {
                fd = ::open(fileName.c_str(), O_RDONLY | O_CLOEXEC);
            }
            if (fd < 0) {
                Debug(Debug::ERROR) << "Cannot open data file " << fileName << "!\n";
                EXIT(EXIT_FAILURE);
            }
            dataFds.push_back(fd);
        }
        batchReader = new BatchReader(threads, dataFds[0], dataMode & USE_IO_URING);
        if (isCompressed(dbtype) == COMPRESSED) {
            frameOffsets.resize(dataFileNames.size());
            frameSizes.resize(dataFileNames.size());
        }
    }

    compression = isLinfinity ? UNCOMPRESSED : isCompressed(dbtype);
    padded = (getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU);
    // A DB may pack two alphabets into a single byte (primary * auxAlphabetSize + aux, so
    // values run past 20); such a DB is flagged with Parameters::DBTYPE_EXTENDED_AUX_SEQ and
    // stores those raw bytes on disk whether it is padded or not. getUnpadded()'s numeric->ASCII
    // re-encoding assumes single-alphabet codes in 0..20 and must be skipped for that layout.
    packedAlphabet = (getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_AUX_SEQ) != 0;

    if(compression == COMPRESSED || padded || frameOffsets.empty() == false){
        compressedBufferSizes = new size_t[threads];
        compressedBuffers = new char*[threads];
        dstream = new ZSTD_DStream*[threads];
        for(int i = 0; i < threads; i++){
            // allocated buffer
            compressedBufferSizes[i] = std::max<size_t>(maxSeqLen+2, 1024);
            compressedBuffers[i] = (char*) malloc(compressedBufferSizes[i]);
            incrementMemory(compressedBufferSizes[i]);
            if(compressedBuffers[i]==NULL){
                Debug(Debug::ERROR) << "Cannot allocate compressedBuffer!\n";
                EXIT(EXIT_FAILURE);
            }
            dstream[i] = ZSTD_createDStream();
            if (dstream==NULL) {
                Debug(Debug::ERROR) << "ZSTD_createDStream() error \n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    closed = 0;
    return isSortedById;
}

template<typename T>
bool DBReader<T>::isSortedByOffset(){
    return sortedByOffset;
}

template <typename T> char* DBReader<T>::mmapData(FILE * file, size_t *dataSize) {
    struct stat sb;
    if (fstat(fileno(file), &sb) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << dataFileName << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }

    *dataSize = sb.st_size;
    int fd =  fileno(file);

    char *ret;
    if(*dataSize > 0){
        if ((dataMode & USE_FREAD) == 0) {
            int mode;
            if (dataMode & USE_WRITABLE) {
                mode = PROT_READ | PROT_WRITE;
            } else {
                mode = PROT_READ;
            }
            ret = static_cast<char*>(mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0));
            if (ret == MAP_FAILED){
                int errsv = errno;
                Debug(Debug::ERROR) << "Failed to mmap memory dataSize=" << *dataSize <<" File=" << dataFileName << ". Error " << errsv << ".\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            ret = static_cast<char*>(malloc(*dataSize));
            Util::checkAllocation(ret, "Not enough system memory to read in the whole data file.");
            incrementMemory(*dataSize);

            size_t result = fread(ret, 1, *dataSize, file);
            if (result != *dataSize) {
                Debug(Debug::ERROR) << "Failed to read in datafile (" << dataFileName << "). Error " << errno << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        return ret;
    }else{
        return NULL;
    }
}

template <typename T> void DBReader<T>::remapData(){
    if ((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0 && index2dbfile == NULL) {
        unmapData();
        for(size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++){
            FILE* dataFile = fopen(dataFileNames[fileIdx].c_str(), "r");
            if (dataFile == NULL) {
                Debug(Debug::ERROR) << "Cannot open data file " << dataFileNames[fileIdx] << "!\n";
                EXIT(EXIT_FAILURE);
            }
            size_t dataSize = 0;
            dataFiles[fileIdx] = mmapData(dataFile, &dataSize);
            if (fclose(dataFile) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << dataFileNames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        dataMapped = true;
    }
}

template <typename T> void DBReader<T>::close(){
    if (dataMode & USE_LOOKUP || dataMode & USE_LOOKUP_REV) {
        delete[] lookup;
    }

    if(dataMode & USE_DATA){
        unmapData();
    }

    if (id2local != NULL) {
        delete[] id2local;
        decrementMemory(size*sizeof(unsigned int));
    }
    if (local2id != NULL) {
        delete[] local2id;
        decrementMemory(size*sizeof(unsigned int));
    }

    if(compressedBuffers){
        for(int i = 0; i < threads; i++){
            ZSTD_freeDStream(dstream[i]);
            free(compressedBuffers[i]);
            decrementMemory(compressedBufferSizes[i]);
        }
        delete [] compressedBuffers;
        delete [] compressedBufferSizes;
        delete [] dstream;
    }

    if(externalData == false) {
        delete[] index;
        decrementMemory(indexSize * sizeof(Index));
    }

    if (index2dbfile != NULL) {
        delete[] index2dbfile;
        index2dbfile = NULL;
        decrementMemory(indexSize * sizeof(unsigned int));
    }
    if (batchReader != NULL) {
        releaseDataPages();
    }
    delete batchReader;
    batchReader = NULL;
    for (size_t fileIdx = 0; fileIdx < dataFds.size(); fileIdx++) {
        ::close(dataFds[fileIdx]);
    }
    dataFds.clear();
    headerOffsetsFileNames.clear();
    for (size_t fileIdx = 0; fileIdx < frameOffsets.size(); fileIdx++) {
        decrementMemory(frameOffsets[fileIdx].size() * sizeof(size_t));
    }
    frameOffsets.clear();
    frameSizes.clear();
    closed = 1;
}

template <typename T> size_t DBReader<T>::bsearch(const Index * index, size_t N, T value)
{
    Index val;
    val.id = value;
    return std::upper_bound(index, index + N, val, Index::compareByIdOnly) - index;
}


template <typename T> char* DBReader<T>::getUnpadded(size_t id, int thrIdx) {
    char *data = getDataUncompressed(id);
    size_t seqLen = getSeqLen(id);

    if (packedAlphabet) {
        // Packed two-alphabet bytes (up to 251) are the on-disk representation of the unpadded DB
        // as well, and are decoded downstream by Sequence's remap tables. Hand them back verbatim;
        // running them through CODE_TO_CHAR below would index a 21-entry table out of bounds.
        // Only the trailing SIMD padding is dropped, which getSeqLen() already excludes.
        memcpy(compressedBuffers[thrIdx], data, seqLen);
        compressedBuffers[thrIdx][seqLen + 0] = '\n';
        compressedBuffers[thrIdx][seqLen + 1] = '\0';
        return compressedBuffers[thrIdx];
    }

    static const char CODE_TO_CHAR[21] = {
            'A', /*  0 */ 'C', /*  1 */ 'D', /*  2 */
            'E', /*  3 */ 'F', /*  4 */ 'G', /*  5 */
            'H', /*  6 */ 'I', /*  7 */ 'K', /*  8 */
            'L', /*  9 */ 'M', /* 10 */ 'N', /* 11 */
            'P', /* 12 */ 'Q', /* 13 */ 'R', /* 14 */
            'S', /* 15 */ 'T', /* 16 */ 'V', /* 17 */
            'W', /* 18 */ 'Y', /* 19 */ 'X'  /* 20 */
    };

    for(size_t i = 0; i < seqLen; i++){
        unsigned char code = static_cast<unsigned char>(data[i]);
        unsigned char baseCode = (code >= 32) ? code - 32 : code;
        // restore masked characters as lowercase with bit twiddling
        compressedBuffers[thrIdx][i] = (code >= 32) ? (CODE_TO_CHAR[baseCode] | ' ') : CODE_TO_CHAR[baseCode];
    }
    compressedBuffers[thrIdx][seqLen + 0] = '\n';
    compressedBuffers[thrIdx][seqLen + 1] = '\0';
    return compressedBuffers[thrIdx];
}

template <typename T> char* DBReader<T>::getDataCompressed(size_t id, int thrIdx) {
    char *data = getDataUncompressed(id);

    unsigned int cSize = *(reinterpret_cast<unsigned int *>(data));

    size_t totalSize = 0;
    const void *cBuff = static_cast<void *>(data + sizeof(unsigned int));
    const char *dataStart = data + sizeof(unsigned int);
    bool isCompressed = (dataStart[cSize] == 0) ? true : false;
    if(isCompressed){
        ZSTD_inBuffer input = {cBuff, cSize, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {compressedBuffers[thrIdx], compressedBufferSizes[thrIdx], 0};
            // size of next compressed block
            size_t toRead = ZSTD_decompressStream(dstream[thrIdx], &output, &input);
            if (ZSTD_isError(toRead)) {
                Debug(Debug::ERROR) << id << " ZSTD_decompressStream " << ZSTD_getErrorName(toRead) << "\n";
                EXIT(EXIT_FAILURE);
            }
            totalSize += output.pos;
        }
        compressedBuffers[thrIdx][totalSize] = '\0';
    }else{
        memcpy(compressedBuffers[thrIdx], cBuff, cSize);
        compressedBuffers[thrIdx][cSize] = '\0';
    }
    return compressedBuffers[thrIdx];
}

template <typename T> size_t DBReader<T>::getAminoAcidDBSize() {
    checkClosed();
    if (headerOffsetsFileNames.empty() == false) {
        Debug(Debug::ERROR) << "Call of getAminoAcidDBSize is not allowed for the headers of a --createdb-mode 3 database " << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)){
        // Get the actual profile column without the null byte per entry
        return (dataSize / Sequence::PROFILE_READIN_SIZE) - size;
    } else {
        // Get the actual number of residues witout \n\0 per entry
        return dataSize - (2 * size);
    }
}

template <typename T> char* DBReader<T>::getData(size_t id, int thrIdx){
    if(compression == COMPRESSED){
        return getDataCompressed(id, thrIdx);
    }else if (padded) {
        return getUnpadded(id, thrIdx);
    } else {
        return getDataUncompressed(id);
    }
}

template <typename T> char* DBReader<T>::getDataUncompressed(size_t id){
    checkClosed();
    if(!(dataMode & USE_DATA)) {
        Debug(Debug::ERROR) << "DBReader is just open in INDEXONLY mode. Call of getData is not allowed" << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getData: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }


    if (local2id != NULL) {
        return getDataByOffset(index[local2id[id]].offset);
    }else if (index2dbfile != NULL) {
        return getDataFromLoadedBuffer(id);
    }else{
        return getDataByOffset(index[id].offset);
    }
}

template <typename T> char* DBReader<T>::getDataByOffset(size_t offset) {
    if (offset >= totalDataSize){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "Size of data: " << totalDataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << offset << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t cnt = 0;
    while ((offset >= dataSizeOffset[cnt] && offset < dataSizeOffset[cnt+1]) == false ) {
        cnt++;
    }
    size_t fileOffset = offset - dataSizeOffset[cnt];
    return dataFiles[cnt]+fileOffset;
}

template <typename T>
void DBReader<T>::touchData(size_t id) {
    if((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
        char *data = getDataUncompressed(id);
        size_t currDataOffset = getOffset(id);
        size_t nextDataOffset = findNextOffsetid(id);
        size_t dataSize = nextDataOffset-currDataOffset;
        magicBytes = Util::touchMemory(data, dataSize);
    }
}

template <typename T> char* DBReader<T>::getDataByDBKey(T dbKey, int thrIdx) {
    size_t id = getId(dbKey);
    if(compression == COMPRESSED ){
        return (id != DB_ENTRY_NOT_FOUND) ? getDataCompressed(id, thrIdx) : NULL;
    } if(padded) {
        return (id != DB_ENTRY_NOT_FOUND) ? getUnpadded(id, thrIdx) : NULL;
    } else{
        return (id != DB_ENTRY_NOT_FOUND) ? getDataByOffset(index[id].offset) : NULL;
    }
}

template <typename T> size_t DBReader<T>::getLookupSize() const {
    checkClosed();
    return lookupSize;
}

template <typename T> size_t DBReader<T>::getSourceSize() const {
    checkClosed();
    return sourceSize;
}

template <typename T> size_t DBReader<T>::getSize() const {
    checkClosed();
    return size;
}

template <typename T> T DBReader<T>::getDbKey (size_t id){
    checkClosed();
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    if (local2id != NULL) {
        id = local2id[id];
    } else if (index2dbfile != NULL) {
        T key = T();
        key += clusterIndex ? loadedBuffer.keys[findLoadedHeaderOffsetsPos(id)] : id;
        return key;
    }
    return index[id].id;
}

template <typename T> size_t DBReader<T>::getLookupIdByKey(T dbKey) {
    if ((dataMode & USE_LOOKUP) == 0) {
        Debug(Debug::ERROR) << "DBReader for datafile=" << dataFileName << ".lookup was not opened with lookup mode\n";
        EXIT(EXIT_FAILURE);
    }
    LookupEntry val;
    val.id = dbKey;
    size_t id = std::upper_bound(lookup, lookup + lookupSize, val, LookupEntry::compareByIdOnly) - lookup;

    return (id < lookupSize && lookup[id].id == dbKey) ? id : SIZE_MAX;
}

template <typename T> size_t DBReader<T>::getLookupIdByAccession(const std::string& accession) {
    if ((dataMode & USE_LOOKUP_REV) == 0) {
        Debug(Debug::ERROR) << "DBReader for datafile=" << dataFileName << ".lookup was not opened with lookup mode\n";
        EXIT(EXIT_FAILURE);
    }
    LookupEntry val;
    val.entryName = accession;
    size_t id = std::upper_bound(lookup, lookup + lookupSize, val, LookupEntry::compareByAccessionOnly) - lookup;

    return (id < lookupSize && lookup[id].entryName.compare(accession) == 0) ? id : SIZE_MAX;
}

template <typename T> T DBReader<T>::getLookupKey(size_t id){
    if (id >= lookupSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".lookup\n";
        Debug(Debug::ERROR) << "getLookupKey: local id (" << id << ") >= db size (" << lookupSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return lookup[id].id;
}

template <typename T> std::string DBReader<T>::getLookupEntryName (size_t id){
    if (id >= lookupSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".lookup\n";
        Debug(Debug::ERROR) << "getLookupEntryName: local id (" << id << ") >= db size (" << lookupSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return lookup[id].entryName;
}

template <typename T> DBKeyType DBReader<T>::getLookupFileNumber(size_t id){
    if (id >= lookupSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".lookup\n";
        Debug(Debug::ERROR) << "getLookupFileNumber: local id (" << id << ") >= db size (" << lookupSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return lookup[id].fileNumber;
}

template<>
void DBReader<DBKeyType>::lookupEntryToBuffer(std::string& buffer, const LookupEntry& entry) {
    buffer.append(SSTR(entry.id));
    buffer.append(1, '\t');
    buffer.append(entry.entryName);
    buffer.append(1, '\t');
    buffer.append(SSTR(entry.fileNumber));
    buffer.append(1, '\n');
}

template<>
void DBReader<std::string>::lookupEntryToBuffer(std::string& buffer, const LookupEntry& entry) {
    buffer.append(entry.id);
    buffer.append(1, '\t');
    buffer.append(entry.entryName);
    buffer.append(1, '\t');
    buffer.append(SSTR(entry.fileNumber));
    buffer.append(1, '\n');
}

template <typename T> T DBReader<T>::getSourceKey(size_t id){
    if (id >= sourceSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".source\n";
        Debug(Debug::ERROR) << "getSource id: local id (" << id << ") >= db size (" << sourceSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return source[id].id;
}

template <typename T> std::string DBReader<T>::getSourceFileName (size_t id){
    if (id >= sourceSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".source\n";
        Debug(Debug::ERROR) << "getSourceFileName: local id (" << id << ") >= db size (" << sourceSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return source[id].fileName;
}

template <typename T> size_t DBReader<T>::getSourceIdByFileName(const std::string& fName) {
    if ((dataMode & USE_SOURCE_REV) == 0) {
        Debug(Debug::ERROR) << "DBReader for datafile=" << dataFileName << ".source was not opened with source mode\n";
        EXIT(EXIT_FAILURE);
    }
    SourceEntry val;
    val.fileName = fName;
    size_t id = std::upper_bound(source, source + sourceSize, val, SourceEntry::compareByFileNameOnly) - source;
    if (id >= sourceSize) {
        Debug(Debug::ERROR) << "Source file " << fName << " exceed source size\n";
        EXIT(EXIT_FAILURE);
    }
    if (source[id].fileName.compare(fName) != 0) {
        Debug(Debug::ERROR) << "Source file " << fName << " not found\n";
        EXIT(EXIT_FAILURE);
    }
    return (id < sourceSize && source[id].fileName.compare(fName) == 0) ? id : SIZE_MAX;
}

template <typename T> void DBReader<T>::sortSourceById(){
    if (source != NULL) {
        SORT_SERIAL(source, source + sourceSize, SourceEntry::compareById);
    }
}

template <typename T> void DBReader<T>::sortSourceByFileName(){
    if (source != NULL) {
        SORT_SERIAL(source, source + sourceSize, SourceEntry::compareByFileName);
    }
}

template <typename T> size_t DBReader<T>::getId (T dbKey){
    if (headerOffsetsFileNames.empty() == false) {
        Debug(Debug::ERROR) << "Call of getId is not allowed for the headers of a --createdb-mode 3 database " << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (index2dbfile != NULL) {
        const size_t id = indexIdToNum(&dbKey);
        return (id < size) ? id : DB_ENTRY_NOT_FOUND;
    }
    size_t id = bsearch(index, size, dbKey);
    if (id2local != NULL) {
        return (id < size && index[id].id == dbKey) ? id2local[id] : DB_ENTRY_NOT_FOUND;
    }
    return (id < size && index[id].id == dbKey ) ? id : DB_ENTRY_NOT_FOUND;
}

template <typename T> size_t DBReader<T>::maxCount(char c) {
    checkClosed();

    size_t max = 0;
    if (compression == COMPRESSED) {
        size_t entries = getSize();
#ifdef OPENMP
        size_t localThreads = std::max(std::min(entries, static_cast<size_t>(threads)), (size_t)1);
#endif
#pragma omp parallel num_threads(localThreads)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10) reduction(max:max)
            for (size_t id = 0; id < entries; id++) {
                char *data = getData(id, thread_idx);
                size_t count = 0;
                for (size_t i = 0; i < getEntryLen(id); ++i) {
                    if (data[i] == c) {
                        count++;
                    }
                }
                max = std::max(max, count);
            }
        }
        return max;
    }

    size_t count = 0;
    for(size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++) {
        size_t dataSize = dataSizeOffset[fileIdx+1] - dataSizeOffset[fileIdx];
        char * data = dataFiles[fileIdx];
        for (size_t i = 0; i < dataSize; ++i) {
            if (data[i] == c) {
                count++;
            }

            if (data[i] == '\0') {
                max = std::max(max, count);
                count = 0;
            }
        }
    }

    return max;
}

template <typename T> void DBReader<T>::checkClosed() const {
    if (closed == 1){
        Debug(Debug::ERROR) << "Trying to read a closed database.\n";
        EXIT(EXIT_FAILURE);
    }
}

template <typename T> size_t DBReader<T>::findIdRange(size_t id) {
    T key = T();
    key += id;
    const size_t range = bsearch(index, indexSize, key);
    return (range < indexSize && index[range].id == key) ? range : range - 1;
}

template <typename T> size_t DBReader<T>::findLoadedHeaderOffsetsPos(size_t id) {
    if (id < loadedBuffer.keyBegin || id >= loadedBuffer.keyEnd || loadedBuffer.headerOffsets.empty()) {
        Debug(Debug::ERROR) << "getData: local id (" << id << ") was not read by loadData from " << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    // each loaded file keeps one offset more than its keys, the end of the last name
    return id - loadedBuffer.keyBegin + findIdRange(id) - findIdRange(loadedBuffer.keyBegin);
}

template <typename T> size_t DBReader<T>::planReadRequests(size_t keyBegin, size_t keyEnd, std::vector<BatchReader::ReadRequest> &readRequests, std::vector<size_t> &headerOffsets, std::vector<size_t> &keys, std::vector<BatchReader::ReadRequest> &compressedReadRequests) {
    // headers: without read requests plan reading the offsets, with the read offsets plan reading the names
    std::vector<BatchReader::ReadRequest> offsetsRequests;
    offsetsRequests.swap(readRequests);
    headerOffsets.clear();
    keys.clear();
    compressedReadRequests.clear();
    const size_t offsetBytes = (clusterIndex ? 2 : 1) * sizeof(uint64_t);
    if (batchReader == NULL || keyBegin >= keyEnd) {
        return 0;
    }
    const size_t pageSize = Util::getPageSize();
    const size_t firstRange = findIdRange(keyBegin);
    for (size_t range = firstRange; range < indexSize && indexIdToNum(&index[range].id) < keyEnd; range++) {
        const size_t firstId = indexIdToNum(&index[range].id);
        const size_t endId = (range + 1 < indexSize) ? indexIdToNum(&index[range + 1].id) : size;
        size_t begin = index[range].offset + (std::max(keyBegin, firstId) - firstId) * index[range].length;
        size_t end = index[range].offset + (std::min(keyEnd, endId) - firstId) * index[range].length;
        size_t fileIdx = index2dbfile[range];
        if (headerOffsetsFileNames.empty() == false && offsetsRequests.empty()) {
            begin = (std::max(keyBegin, firstId) - firstId) * offsetBytes;
            end = (std::min(keyEnd, endId) - firstId + 1) * offsetBytes;
            fileIdx += dataFileNames.size();
        } else if (clusterIndex) {
            const BatchReader::ReadRequest &offsetsRequest = offsetsRequests[range - firstRange];
            const uint64_t *cells = reinterpret_cast<const uint64_t *>(offsetsRequest.destination + keyBegin * offsetBytes - offsetsRequest.fileOffset);
            for (size_t i = 0; i <= keyEnd - keyBegin; i++) {
                keys.push_back(cells[2 * i]);
                headerOffsets.push_back(cells[2 * i + 1]);
            }
            begin = headerOffsets.front();
            end = headerOffsets.back();
        } else if (headerOffsetsFileNames.empty() == false) {
            const BatchReader::ReadRequest &offsetsRequest = offsetsRequests[range - firstRange];
            const uint64_t *offsets = reinterpret_cast<const uint64_t *>(offsetsRequest.destination + (std::max(keyBegin, firstId) - firstId) * sizeof(uint64_t) - offsetsRequest.fileOffset);
            const size_t keyCount = std::min(keyEnd, endId) - std::max(keyBegin, firstId);
            headerOffsets.insert(headerOffsets.end(), offsets, offsets + keyCount + 1);
            begin = offsets[0];
            end = offsets[keyCount];
        }
        const int fd = dataFds[fileIdx];
        const bool isCompressedNames = frameOffsets.empty() == false && offsetsRequests.empty() == false;
        if (isCompressedNames && frameOffsets[fileIdx].empty()) {
            FILE *file = FileUtil::openFileOrDie(dataFileNames[fileIdx].c_str(), "r", true);
            const size_t footerSize = 9;
            char footer[footerSize];
            bool isValid = fseek(file, -static_cast<long>(footerSize), SEEK_END) == 0 && fread(footer, sizeof(char), footerSize, file) == footerSize;
            const size_t fileSize = isValid ? ftell(file) : 0;
            uint32_t frameCount = 0;
            uint32_t magic = 0;
            if (isValid) {
                memcpy(&frameCount, footer, sizeof(uint32_t));
                memcpy(&magic, footer + 5, sizeof(uint32_t));
            }
            const size_t seekTableSize = 8 + static_cast<size_t>(frameCount) * 8 + footerSize;
            isValid = isValid && magic == 0x8F92EAB1 && (footer[4] & 0xFC) == 0 && seekTableSize <= fileSize;
            std::vector<uint32_t> seekTable(isValid ? 2 + 2 * static_cast<size_t>(frameCount) : 0);
            isValid = isValid && fseek(file, fileSize - seekTableSize, SEEK_SET) == 0 && fread(seekTable.data(), sizeof(uint32_t), seekTable.size(), file) == seekTable.size()
                      && seekTable[0] == (ZSTD_MAGIC_SKIPPABLE_START | 0xE) && seekTable[1] == seekTableSize - 8;
            if (fclose(file) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << dataFileNames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            frameSizes[fileIdx] = (isValid && frameCount > 0) ? seekTable[3] : 1;
            frameOffsets[fileIdx].assign(1, 0);
            for (size_t frame = 0; isValid && frame < frameCount; frame++) {
                const size_t frameSize = seekTable[3 + 2 * frame];
                isValid = frameSize > 0 && (frameSize == frameSizes[fileIdx] || (frame + 1 == frameCount && frameSize < frameSizes[fileIdx]));
                frameOffsets[fileIdx].push_back(frameOffsets[fileIdx].back() + seekTable[2 + 2 * frame]);
            }
            if (isValid == false || frameOffsets[fileIdx].back() + seekTableSize != fileSize) {
                Debug(Debug::ERROR) << "Invalid zstd seek table in " << dataFileNames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            incrementMemory(frameOffsets[fileIdx].size() * sizeof(size_t));
        }
        const size_t alignment = isCompressedNames ? frameSizes[fileIdx] : pageSize;
        if (readRequests.empty() || readRequests.back().fileToRead != fd) {
            const char *fileName = (fileIdx < dataFileNames.size()) ? dataFileNames[fileIdx].c_str() : headerOffsetsFileNames[fileIdx - dataFileNames.size()].c_str();
            BatchReader::ReadRequest request = { fd, fileName, NULL, begin / alignment * alignment, 0, 0 };
            readRequests.push_back(request);
        }
        BatchReader::ReadRequest &request = readRequests.back();
        request.minBytesToRead = end - request.fileOffset;
        request.bytesToRead = ((end + alignment - 1) / alignment * alignment - request.fileOffset + pageSize - 1) / pageSize * pageSize;
        if (isCompressedNames) {
            const size_t endFrame = (end + alignment - 1) / alignment;
            if (endFrame >= frameOffsets[fileIdx].size()) {
                Debug(Debug::ERROR) << "Header offsets exceed the zstd seek table of " << dataFileNames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            const size_t compressedEnd = frameOffsets[fileIdx][endFrame];
            BatchReader::ReadRequest compressedRequest = { fd, request.fileName, NULL, frameOffsets[fileIdx][request.fileOffset / alignment] / pageSize * pageSize, 0, 0 };
            compressedRequest.minBytesToRead = compressedEnd - compressedRequest.fileOffset;
            compressedRequest.bytesToRead = (compressedEnd + pageSize - 1) / pageSize * pageSize - compressedRequest.fileOffset;
            compressedReadRequests.push_back(compressedRequest);
        }
    }
    size_t bytes = 0;
    for (size_t i = 0; i < readRequests.size(); i++) {
        bytes += readRequests[i].bytesToRead;
    }
    for (size_t i = 0; i < compressedReadRequests.size(); i++) {
        bytes += compressedReadRequests[i].bytesToRead;
    }
    return bytes;
}

template <typename T> size_t DBReader<T>::getBufferSizeToLoad(size_t keyBegin, size_t keyEnd) {
    if (headerOffsetsFileNames.empty() == false) {
        Debug(Debug::ERROR) << "Call of getBufferSizeToLoad is not allowed for the headers of a --createdb-mode 3 database " << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<BatchReader::ReadRequest> readRequests;
    std::vector<size_t> headerOffsets;
    std::vector<size_t> keys;
    std::vector<BatchReader::ReadRequest> compressedReadRequests;
    return planReadRequests(keyBegin, keyEnd, readRequests, headerOffsets, keys, compressedReadRequests);
}

template <typename T> bool DBReader<T>::prepareBufferToLoad(BufferInfo &bufferToLoadInto, size_t keyBegin, size_t keyEnd) {
    std::vector<BatchReader::ReadRequest> readRequests;
    std::vector<size_t> headerOffsets;
    std::vector<size_t> keys;
    std::vector<BatchReader::ReadRequest> compressedReadRequests;
    // headers: the first pass reads the offsets of the names into the buffer
    for (bool readOffsets = headerOffsetsFileNames.empty() == false; ; readOffsets = false) {
        const size_t bytes = planReadRequests(keyBegin, keyEnd, readRequests, headerOffsets, keys, compressedReadRequests);
        if (readRequests.empty()) {
            return false;
        }
        if (bytes > bufferToLoadInto.allocatedBytes) {
            batchReader->waitAsyncRead();
            free(bufferToLoadInto.buffer);
            decrementMemory(bufferToLoadInto.allocatedBytes);
            // O_DIRECT reads need a page aligned buffer
            if (posix_memalign(reinterpret_cast<void **>(&bufferToLoadInto.buffer), Util::getPageSize(), bytes) != 0) {
                Debug(Debug::ERROR) << "Cannot allocate " << bytes << " bytes to read " << dataFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
            incrementMemory(bytes);
            bufferToLoadInto.allocatedBytes = bytes;
        }
        char *destination = bufferToLoadInto.buffer;
        for (size_t i = 0; i < readRequests.size(); i++) {
            readRequests[i].destination = destination;
            destination += readRequests[i].bytesToRead;
        }
        for (size_t i = 0; i < compressedReadRequests.size(); i++) {
            compressedReadRequests[i].destination = destination;
            destination += compressedReadRequests[i].bytesToRead;
        }
        if (readOffsets == false) {
            break;
        }
        batchReader->readSync(readRequests);
    }
    bufferToLoadInto.readRequests.swap(readRequests);
    bufferToLoadInto.headerOffsets.swap(headerOffsets);
    bufferToLoadInto.keys.swap(keys);
    bufferToLoadInto.compressedReadRequests.swap(compressedReadRequests);
    bufferToLoadInto.keyBegin = keyBegin;
    bufferToLoadInto.keyEnd = keyEnd;
    return true;
}

template <typename T> void DBReader<T>::loadData(size_t keyBegin, size_t keyEnd) {
    if (prefetchBuffer.readRequests.empty() == false && prefetchBuffer.keyBegin == keyBegin && prefetchBuffer.keyEnd == keyEnd) {
        batchReader->waitAsyncRead();
        std::swap(loadedBuffer, prefetchBuffer);
    } else if (prepareBufferToLoad(loadedBuffer, keyBegin, keyEnd)) {
        batchReader->readSync(frameOffsets.empty() ? loadedBuffer.readRequests : loadedBuffer.compressedReadRequests);
    }
    decompressLoadedBuffer();
}

template <typename T> void DBReader<T>::prefetchData(size_t keyBegin, size_t keyEnd) {
    if (prepareBufferToLoad(prefetchBuffer, keyBegin, keyEnd)) {
        batchReader->readAsync(frameOffsets.empty() ? prefetchBuffer.readRequests : prefetchBuffer.compressedReadRequests);
    }
}

template <typename T> void DBReader<T>::decompressLoadedBuffer() {
    if (frameOffsets.empty() || loadedBuffer.readRequests.empty()) {
        return;
    }
    // header read requests are one per <db>.nodeinfo row
    const size_t firstRange = findIdRange(loadedBuffer.keyBegin);
    std::vector<std::pair<size_t, size_t> > requestFrames;
    for (size_t i = 0; i < loadedBuffer.readRequests.size(); i++) {
        const BatchReader::ReadRequest &request = loadedBuffer.readRequests[i];
        const size_t frameSize = frameSizes[index2dbfile[firstRange + i]];
        for (size_t frame = request.fileOffset / frameSize; frame * frameSize < request.fileOffset + request.minBytesToRead; frame++) {
            requestFrames.push_back(std::make_pair(i, frame));
        }
    }
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < requestFrames.size(); i++) {
            const size_t requestIdx = requestFrames[i].first;
            const size_t frame = requestFrames[i].second;
            const BatchReader::ReadRequest &request = loadedBuffer.readRequests[requestIdx];
            const BatchReader::ReadRequest &compressedRequest = loadedBuffer.compressedReadRequests[requestIdx];
            const size_t fileIdx = index2dbfile[firstRange + requestIdx];
            const size_t frameSize = frameSizes[fileIdx];
            const size_t compressedOffset = frameOffsets[fileIdx][frame];
            const size_t size = ZSTD_decompressDCtx(dstream[thread_idx], request.destination + (frame * frameSize - request.fileOffset), frameSize,
                                                    compressedRequest.destination + (compressedOffset - compressedRequest.fileOffset), frameOffsets[fileIdx][frame + 1] - compressedOffset);
            if (ZSTD_isError(size) || size < std::min(frameSize, request.fileOffset + request.minBytesToRead - frame * frameSize)) {
                Debug(Debug::ERROR) << request.fileName << " ZSTD_decompressDCtx " << (ZSTD_isError(size) ? ZSTD_getErrorName(size) : "frame is shorter than its seek table entry") << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
}

template <typename T> char* DBReader<T>::getDataFromLoadedBuffer(size_t key) {
    if (batchReader != NULL && key < size) {
        const size_t range = findIdRange(key);
        size_t fileOffset = index[range].offset + (key - indexIdToNum(&index[range].id)) * index[range].length;
        size_t length = index[range].length;
        if (headerOffsetsFileNames.empty() == false) {
            const size_t headerOffsetsPos = findLoadedHeaderOffsetsPos(key);
            fileOffset = loadedBuffer.headerOffsets[headerOffsetsPos];
            length = loadedBuffer.headerOffsets[headerOffsetsPos + 1] - fileOffset;
        }
        const int fd = dataFds[index2dbfile[range]];
        for (size_t i = 0; i < loadedBuffer.readRequests.size(); i++) {
            const BatchReader::ReadRequest &request = loadedBuffer.readRequests[i];
            if (request.fileToRead == fd && fileOffset >= request.fileOffset && fileOffset + length <= request.fileOffset + request.minBytesToRead) {
                return request.destination + (fileOffset - request.fileOffset);
            }
        }
    }
    Debug(Debug::ERROR) << "getData: local id (" << key << ") was not read by loadData from " << dataFileName << "\n";
    EXIT(EXIT_FAILURE);
}

template<typename T>
bool DBReader<T>::readIndex(char *data, size_t indexDataSize, Index *index, size_t & dataSize) {
#ifdef OPENMP
    int threadCnt = 1;
    const int totalThreadCnt = threads;
    if (totalThreadCnt >= 4) { 
	threadCnt = 4;
    }
#endif


    size_t isSortedById = true;
    size_t globalIdOffset = 0;
    unsigned int localMaxSeqLen = 0;
    size_t localDataSize = 0;

    DBKeyType localLastKey = 0;
    const size_t BATCH_SIZE = 1048576;
#pragma omp parallel num_threads(threadCnt) reduction(max: localMaxSeqLen, localLastKey) reduction(+: localDataSize) reduction(min:isSortedById)
    {
        size_t currPos = 0;
        char* indexDataChar = (char *) data;
        const char * cols[3];
        size_t lineStartId = __sync_fetch_and_add(&(globalIdOffset), BATCH_SIZE);
        T prevId=T(); // makes 0 or empty string
        size_t currLine = 0;

        while (currPos < indexDataSize){
            if (currLine >= this->size) {
                Debug(Debug::ERROR) << "Corrupt memory, too many entries: " << currLine << " >= " << this->size << "\n";
                EXIT(EXIT_FAILURE);
            }
            if(currLine == lineStartId){
                for(size_t startIndex = lineStartId; startIndex < lineStartId + BATCH_SIZE && currPos < indexDataSize; startIndex++){
                    Util::getWordsOfLine(indexDataChar, cols, 3);
                    readIndexId(&index[startIndex].id, indexDataChar, cols);
                    isSortedById *= (index[startIndex].id >= prevId);
                    size_t offset = Util::fast_atoi<size_t>(cols[1]);
                    size_t length = Util::fast_atoi<size_t>(cols[2]);
                    localDataSize += length;
                    index[startIndex].offset = offset;
                    index[startIndex].length = length;
                    localMaxSeqLen = std::max(static_cast<unsigned int>(length), localMaxSeqLen);
                    indexDataChar = Util::skipLine(indexDataChar);
                    currPos = indexDataChar - (char *) data;
                    localLastKey = std::max(localLastKey, indexIdToNum(&index[startIndex].id));
                    prevId = index[startIndex].id;
                    currLine++;
                }
                lineStartId = __sync_fetch_and_add(&(globalIdOffset), BATCH_SIZE);
            }else{
                indexDataChar = Util::skipLine(indexDataChar);
                currPos = indexDataChar - (char *) data;
                currLine++;
            }

        }
    }
    dataSize = localDataSize;
    maxSeqLen = localMaxSeqLen;
    lastKey = localLastKey;
    return isSortedById;
}

template<typename T> T DBReader<T>::getLastKey() {
    return lastKey;
}

template<>
void DBReader<std::string>::readIndexId(std::string* id, char* line, const char** cols){
    ptrdiff_t keySize =  ((cols[1] - 1) - line) ;
    id->assign(line, keySize);
}
template<>
void DBReader<DBKeyType>::readIndexId(DBKeyType* id, char*, const char** cols) {
    *id = Util::fast_atoi<DBKeyType>(cols[0]);
}

template<>
DBKeyType DBReader<std::string>::indexIdToNum(std::string * id){
    return id->size();
}
template<>
DBKeyType DBReader<DBKeyType>::indexIdToNum(DBKeyType * id) {
    return *id;
}

template <typename T> void DBReader<T>::unmapData() {
    if (dataMapped == true) {
        for(size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
            size_t fileSize = dataSizeOffset[fileIdx+1] -dataSizeOffset[fileIdx];
            if(fileSize > 0) {
                if (didMlock == true) {
                    munlock(dataFiles[fileIdx], fileSize);
                }
                if ((dataMode & USE_FREAD) == 0) {
                    if (munmap(dataFiles[fileIdx], fileSize) < 0) {
                        Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << fileSize << " File=" << dataFileName
                                            << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                } else {
                    free(dataFiles[fileIdx]);
                    decrementMemory(dataSize);
                }
            }
        }
    }

    didMlock = false;
    dataMapped = false;
}

template <typename T>  size_t DBReader<T>::getDataOffset(T i) {
    size_t id = bsearch(index, size, i);
    return index[id].offset;
}

template <>
size_t DBReader<DBKeyType>::indexMemorySize(const DBReader<DBKeyType> &idx) {
    size_t memSize = // size + dataSize
            2 * sizeof(size_t)
            // maxSeqLen + lastKey + dbtype
            + sizeof(DBKeyType) + sizeof(int) + sizeof(unsigned int)
            // index
            + idx.size * sizeof(DBReader<DBKeyType>::Index);

    return memSize;
}

template <>
char* DBReader<DBKeyType>::serialize(const DBReader<DBKeyType> &idx) {
    char* data = (char*) malloc(indexMemorySize(idx));
    char* p = data;
    memcpy(p, &idx.size, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.dataSize, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.lastKey, sizeof(DBKeyType));
    p += sizeof(DBKeyType);
    memcpy(p, &idx.dbtype, sizeof(int));
    p += sizeof(int);
    memcpy(p, &idx.maxSeqLen, sizeof(unsigned int));
    p += sizeof(unsigned int);
    memcpy(p, idx.index, idx.size * sizeof(DBReader<DBKeyType>::Index));
    p += idx.size * sizeof(DBReader<DBKeyType>::Index);
    return data;
}

template <>
DBReader<DBKeyType> *DBReader<DBKeyType>::unserialize(const char* data, int threads) {
    const char* p = data;
    size_t size = *((size_t*)p);
    p += sizeof(size_t);
    size_t dataSize = *((size_t*)p);
    p += sizeof(size_t);
    DBKeyType lastKey = *((DBKeyType*)p);
    p += sizeof(DBKeyType);
    int dbType = *((int*)p);
    p += sizeof(int);
    unsigned int maxSeqLen = *((unsigned int*)p);
    p += sizeof(unsigned int);
    DBReader<DBKeyType>::Index *idx = (DBReader<DBKeyType>::Index *)p;
    p += size * sizeof(DBReader<DBKeyType>::Index);

    return new DBReader<DBKeyType>(idx, size, dataSize, lastKey, dbType, maxSeqLen, threads);
}

template<typename T>
void DBReader<T>::setData(char *data, size_t dataSize) {
    if(dataFiles == NULL){
        dataFiles = new char*[1];
        dataSizeOffset = new size_t[2];
        dataSizeOffset[0] = 0;
        dataSizeOffset[1] = dataSize;
        totalDataSize = dataSize;
        dataFileCnt = 1;
        dataFiles[0] = data;
    }else{
        Debug(Debug::ERROR) << "DataFiles is already set." << "\n";
        EXIT(EXIT_FAILURE);
    }
}

template<typename T>
void DBReader<T>::setMode(const int mode) {
    this->dataMode = mode;
}

template<typename T>
size_t DBReader<T>::getOffset(size_t id) {
    if (id >= size){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "getOffset: local id (" << id << ") >= db size (" << size << ")\n";
        EXIT(EXIT_FAILURE);
    }
    if (local2id != NULL) {
        id = local2id[id];
    }
    return index[id].offset;
}

template<typename T>
size_t DBReader<T>::findNextOffsetid(size_t id) {
    size_t idOffset = getOffset(id);
    size_t nextOffset = SIZE_MAX;
    for(size_t i = 0; i < size; i++){
        if(index[i].offset > idOffset && index[i].offset < nextOffset){
            nextOffset=index[i].offset;
        }
    }
    // if the offset is the last element in the index
    if(nextOffset == SIZE_MAX){
        nextOffset = dataSizeOffset[dataFileCnt];
    }
    return nextOffset;
}

template<typename T>
int DBReader<T>::isCompressed(int dbtype) {
    return (dbtype & (1 << 31)) ? COMPRESSED : UNCOMPRESSED;
}


template <typename T> void DBReader<T>::releaseDataPages() {
    if (batchReader != NULL) {
        batchReader->waitAsyncRead();
        free(loadedBuffer.buffer);
        decrementMemory(loadedBuffer.allocatedBytes);
        loadedBuffer = BufferInfo();
        free(prefetchBuffer.buffer);
        decrementMemory(prefetchBuffer.allocatedBytes);
        prefetchBuffer = BufferInfo();
        return;
    }
#ifdef MADV_DONTNEED
    // Util::madviseLogged goes through posix_madvise, which glibc turns into a no-op for
    // POSIX_MADV_DONTNEED, so this has to call madvise itself
    if (dataMapped == false || (dataMode & USE_FREAD) || (dataMode & USE_WRITABLE)) {
        return;
    }
    for (size_t i = 0; i < dataFileCnt; i++) {
        madvise(dataFiles[i], dataSizeOffset[i + 1] - dataSizeOffset[i], MADV_DONTNEED);
    }
#endif
}
template<typename T>
void DBReader<T>::setSequentialAdvice() {
    for(size_t i = 0; i < dataFileCnt; i++){
        size_t dataSize = dataSizeOffset[i+1] - dataSizeOffset[i];
        Util::madviseLogged(dataFiles[i], dataSize, POSIX_MADV_SEQUENTIAL, dataFileName);
    }
}

template<typename T>
void DBReader<T>::readLookup(char *data, size_t dataSize, DBReader::LookupEntry *lookup) {
    size_t i = 0;
    size_t currPos = 0;
    char* lookupData = (char *) data;
    const char * cols[3];
    while (currPos < dataSize){
        if (i >= this->lookupSize) {
            Debug(Debug::ERROR) << "Corrupt memory, too many entries!\n";
            EXIT(EXIT_FAILURE);
        }
        Util::getWordsOfLine(lookupData, cols, 3);
        lookup[i].id = Util::fast_atoi<size_t>(cols[0]);
        lookup[i].entryName = std::string(cols[1], (cols[2] - cols[1]) - 1);
        lookup[i].fileNumber = Util::fast_atoi<DBKeyType>(cols[2]);
        lookupData = Util::skipLine(lookupData);

        currPos = lookupData - (char *) data;

        i++;
    }
}

template <typename T>
void DBReader<T>::readSource(char *data, size_t dataSize, DBReader::SourceEntry *source) {
    size_t i=0;
    size_t currPos = 0;
    char* sourceData = (char *) data;
    const char * cols[3];
    while (currPos < dataSize){
        if (i >= this->sourceSize) {
            Debug(Debug::ERROR) << "Corrupt memory, too many entries!\n";
            EXIT(EXIT_FAILURE);
        }
        Util::getFieldsOfLine(sourceData, cols, 3);
        source[i].id = Util::fast_atoi<size_t>(cols[0]);
        std::string fileName = std::string(cols[1], (cols[2] - cols[1]));
        size_t lastDotPosition = fileName.rfind('.');

        if (lastDotPosition != std::string::npos) {
            fileName = fileName.substr(0, lastDotPosition);
        }
        source[i].fileName = fileName;
        sourceData = Util::skipLine(sourceData);

        currPos = sourceData - (char *) data;

        i++;
    }
}

// TODO: Move to DbUtils?

template<typename T>
void DBReader<T>::moveDatafiles(const std::vector<std::string>& files, const std::string& destination) {
    for (size_t i = 0; i < files.size(); i++) {
        std::string extention = files[i].substr(files[i].find_last_of(".") + 1);
        if (Util::isNumber(extention)) {
            std::string dst = (destination + "." + extention);
            FileUtil::move(files[i].c_str(), dst.c_str());
        } else {
            if (files.size() > 1) {
                Debug(Debug::ERROR) << "Both merged and unmerged database exist at the same path\n";
                EXIT(EXIT_FAILURE);
            }
            
            FileUtil::move(files[i].c_str(), destination.c_str());
        }
    }
}

template<typename T>
void DBReader<T>::moveDb(const std::string &srcDbName, const std::string &dstDbName) {
    std::vector<std::string> files = FileUtil::findDatafiles(srcDbName.c_str());
    moveDatafiles(files, dstDbName);

    // --createdb-mode 3 keeps one header index per data file and the clusthashfast representative bits per node
    moveDatafiles(FileUtil::findDatafiles((srcDbName + ".index").c_str()), dstDbName + ".index");
    if (FileUtil::fileExists((srcDbName + ".index").c_str())) {
        FileUtil::move((srcDbName + ".index").c_str(), (dstDbName + ".index").c_str());
    }
    moveDatafiles(FileUtil::findDatafiles((srcDbName + ".bits").c_str()), dstDbName + ".bits");
    if (FileUtil::fileExists((srcDbName + ".dbtype").c_str())) {
        FileUtil::move((srcDbName + ".dbtype").c_str(), (dstDbName + ".dbtype").c_str());
    }
    if (FileUtil::fileExists((srcDbName + ".lookup").c_str())) {
        FileUtil::move((srcDbName + ".lookup").c_str(), (dstDbName + ".lookup").c_str());
    }
    if (FileUtil::fileExists((srcDbName + ".nodeinfo").c_str())) {
        FileUtil::move((srcDbName + ".nodeinfo").c_str(), (dstDbName + ".nodeinfo").c_str());
    }
}

template<typename T>
void DBReader<T>::removeDb(const std::string &databaseName){
    std::vector<std::string> files = FileUtil::findDatafiles(databaseName.c_str());
    for (size_t i = 0; i < files.size(); ++i) {
        FileUtil::remove(files[i].c_str());
    }
    std::vector<std::string> indexFiles = FileUtil::findDatafiles((databaseName + ".index").c_str());
    for (size_t i = 0; i < indexFiles.size(); ++i) {
        FileUtil::remove(indexFiles[i].c_str());
    }
    if (FileUtil::fileExists((databaseName + ".index").c_str())) {
        FileUtil::remove((databaseName + ".index").c_str());
    }
    std::vector<std::string> bitsFiles = FileUtil::findDatafiles((databaseName + ".bits").c_str());
    for (size_t i = 0; i < bitsFiles.size(); ++i) {
        FileUtil::remove(bitsFiles[i].c_str());
    }
    std::string dbTypeFile = databaseName + ".dbtype";
    if (FileUtil::fileExists(dbTypeFile.c_str())) {
        FileUtil::remove(dbTypeFile.c_str());
    }
    std::string sourceFile = databaseName + ".source";
    if (FileUtil::fileExists(sourceFile.c_str())) {
        FileUtil::remove(sourceFile.c_str());
    }
    std::string lookupFile = databaseName + ".lookup";
    if (FileUtil::fileExists(lookupFile.c_str())) {
        FileUtil::remove(lookupFile.c_str());
    }
    std::string nodeInfoFile = databaseName + ".nodeinfo";
    if (FileUtil::fileExists(nodeInfoFile.c_str())) {
        FileUtil::remove(nodeInfoFile.c_str());
    }
}

typedef void (*DbAction)(const std::string &, const std::string &);
void copyLinkDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags, DbAction action) {
    if (dbFilesFlags & DBFiles::DATA) {
        std::vector<std::string> names = FileUtil::findDatafiles(databaseName.c_str());
        if (names.size() == 1) {
            action(names[0], outDb);
        } else {
            for (size_t i = 0; i < names.size(); i++) {
                std::string::size_type idx = names[i].rfind('.');
                std::string ext;
                if (idx != std::string::npos) {
                    ext = names[i].substr(idx);
                } else {
                    Debug(Debug::ERROR) << "File extension was not found but it is expected to be there!\n"
                                        << "Filename: " << names[i] << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                action(names[i], outDb + ext);
            }
        }
    }

    struct DBSuffix {
        DBFiles::Files flag;
        const char* suffix;
    };

    const DBSuffix suffices[] = {
        { DBFiles::DATA_INDEX,    ".index"            },
        { DBFiles::DATA_INDEX,    ".bits"             },
        { DBFiles::DATA_DBTYPE,   ".dbtype"           },
        { DBFiles::HEADER,        "_h"                },
        { DBFiles::HEADER_INDEX,  "_h.index"          },
        { DBFiles::HEADER_INDEX,  ".nodeinfo"         },
        { DBFiles::HEADER_DBTYPE, "_h.dbtype"         },
        { DBFiles::LOOKUP,        ".lookup"           },
        { DBFiles::SOURCE,        ".source"           },
        { DBFiles::TAX_MAPPING,   "_mapping"          },
        { DBFiles::TAX_NAMES,     "_names.dmp"        },
        { DBFiles::TAX_NODES,     "_nodes.dmp"        },
        { DBFiles::TAX_MERGED,    "_merged.dmp"       },
        { DBFiles::TAX_MERGED,    "_taxonomy"         },
        { DBFiles::CA3M_DATA,     "_ca3m.ffdata"      },
        { DBFiles::CA3M_INDEX,    "_ca3m.ffindex"     },
        { DBFiles::CA3M_SEQ,      "_sequence.ffdata"  },
        { DBFiles::CA3M_SEQ_IDX,  "_sequence.ffindex" },
        { DBFiles::CA3M_HDR,      "_header.ffdata"    },
        { DBFiles::CA3M_HDR_IDX,  "_header.ffindex"   },
    };

    // a --createdb-mode 3 database keeps the file numbers of its headers and header index files
    for (size_t i = 0; i < ARRAY_SIZE(suffices); ++i) {
        if (dbFilesFlags & suffices[i].flag) {
            std::vector<std::string> names = FileUtil::findDatafiles((databaseName + suffices[i].suffix).c_str());
            for (size_t j = 0; j < names.size(); j++) {
                action(names[j], outDb + names[j].substr(databaseName.size()));
            }
        }
    }
}

template<typename T>
void DBReader<T>::aliasDb(const std::string &databaseName, const std::string &alias, DBFiles::Files dbFilesFlags) {
    copyLinkDb(databaseName, alias, dbFilesFlags, FileUtil::symlinkAlias);
}

template<typename T>
void DBReader<T>::softlinkDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags) {
    copyLinkDb(databaseName, outDb, dbFilesFlags, FileUtil::symlinkAbs);
}

template<typename T>
void DBReader<T>::copyDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags) {
    copyLinkDb(databaseName, outDb, dbFilesFlags, FileUtil::copyFile);
}

template<typename T>
std::vector<size_t> DBReader<T>::getFileKeyOffsets(size_t fileCount) {
    std::vector<size_t> offsets(fileCount + 1, 0);
    for (size_t row = 0; row < indexSize; row++) {
        offsets[index2dbfile[row] + 1] = (row + 1 < indexSize) ? indexIdToNum(&index[row + 1].id) : size;
    }
    for (size_t file = 0; file < fileCount; file++) {
        offsets[file + 1] = std::max(offsets[file + 1], offsets[file]);
    }
    return offsets;
}

template<typename T>
void DBReader<T>::decomposeDomainByAminoAcid(size_t worldRank, size_t worldSize, size_t *startEntry, size_t *numEntries){
    if (headerOffsetsFileNames.empty() == false) {
        Debug(Debug::ERROR) << "Call of decomposeDomainByAminoAcid is not allowed for the headers of a --createdb-mode 3 database " << dataFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t dataSize = getDataSize();
    const size_t dbEntries = getSize();
    if (worldSize > dataSize) {
        // Assume the domain numEntries is greater than the world numEntries.
        Debug(Debug::ERROR) << "World Size: " << worldSize << " dbSize: " << dataSize << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (worldSize == 1) {
        *startEntry = 0;
        *numEntries = dbEntries;
        return;
    }
    
    if (index2dbfile != NULL) {
        // createdb-mode 3(linfinity)
        const size_t chunkSize = ceil(static_cast<double>(dataSize) / static_cast<double>(worldSize));
        size_t currentRank = 0;
        size_t sumCharsAssigned = 0;
        *startEntry = 0;
        *numEntries = 0;
        for (size_t row = 0; row < indexSize; row++) {
            const size_t rowEntries = ((row + 1 < indexSize) ? indexIdToNum(&index[row + 1].id) : dbEntries) - indexIdToNum(&index[row].id);
            if (row > 0 && index[row].length != index[row - 1].length && sumCharsAssigned >= (currentRank + 1) * chunkSize) {
                currentRank++;
            }
            sumCharsAssigned += rowEntries * index[row].length;
            if (currentRank < worldRank) {
                *startEntry += rowEntries;
            } else if (currentRank == worldRank) {
                *numEntries += rowEntries;
            }
        }
        return;
    }

    if (dbEntries <= worldSize) {
        *startEntry = worldRank < dbEntries ? worldRank : 0;
        *numEntries = worldRank < dbEntries ? 1 : 0;
        return;
    }

    size_t chunkSize = ceil(static_cast<double>(dataSize) / static_cast<double>(worldSize));

    size_t *entriesPerWorker = (size_t*)calloc(worldSize, sizeof(size_t));

    size_t currentRank = 0;
    size_t sumCharsAssignedToCurrRank = 0;
    for (size_t i = 0; i < dbEntries; ++i) {
        if (sumCharsAssignedToCurrRank >= chunkSize) {
            sumCharsAssignedToCurrRank = 0;
            currentRank++;
        }
        sumCharsAssignedToCurrRank += index[i].length;
        entriesPerWorker[currentRank] += 1;
    }

    *startEntry = 0;
    *numEntries = entriesPerWorker[worldRank];
    for (size_t j = 0; j < worldRank; ++j) {
        *startEntry += entriesPerWorker[j];
    }
    free(entriesPerWorker);
}

template class DBReader<DBKeyType>;
template class DBReader<std::string>;
