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

#if defined(__linux__)
#include <sys/syscall.h>
#if defined(HAVE_LINUX_IO_URING) && defined(__NR_io_uring_setup) && defined(__NR_io_uring_enter)
#include <linux/io_uring.h>
#define HAVE_IO_URING 1
#endif
#endif

// a per-thread bounce buffer starts here and grows on demand, so one huge entry cannot preallocate threads*maxSeqLen
static const size_t BOUNCE_BUFFER_PREALLOC = 64 * 1024;
// smallest room left for one ZSTD_decompressStream call, so a nearly full buffer cannot creep forward
static const size_t DECOMPRESS_MIN_ROOM = 4096;

#include "MemoryMapped.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

template <typename T>
DBReader<T>::DBReader(const char* dataFileName_, const char* indexFileName_, int threads, int dataMode) :
threads(threads), dataMode(dataMode), dataFileName(strdup(dataFileName_)),
        indexFileName(strdup(indexFileName_)), size(0), dataFiles(NULL), dataFds(NULL), directBuffers(NULL), ioBatch(NULL),
        dataSizeOffset(NULL), dataFileCnt(0),
        totalDataSize(0), dataSize(0), lastKey(T()), closed(1), dbtype(Parameters::DBTYPE_GENERIC_DB),
        compressedBuffers(NULL), compressedBufferSizes(NULL), index(NULL), id2local(NULL), local2id(NULL),
        dataMapped(false), accessType(0), externalData(false), didMlock(false), ioBufferedBatch(false), ioCacheAdvice(false), directIoAlign(0)
{}

template <typename T>
DBReader<T>::DBReader(DBReader<T>::Index *index, size_t size, size_t dataSize, T lastKey,
        int dbType, unsigned int maxSeqLen, int threads) :
        threads(threads), dataMode(USE_INDEX), dataFileName(NULL), indexFileName(NULL),
        size(size), dataFiles(NULL), dataFds(NULL), directBuffers(NULL), ioBatch(NULL),
        dataSizeOffset(NULL), dataFileCnt(0), totalDataSize(0), dataSize(dataSize), lastKey(lastKey),
        maxSeqLen(maxSeqLen), closed(1), dbtype(dbType), compressedBuffers(NULL), compressedBufferSizes(NULL), index(index), sortedByOffset(true),
        id2local(NULL), local2id(NULL), dataMapped(false), accessType(NOSORT), externalData(true), didMlock(false),
        ioBufferedBatch(false), ioCacheAdvice(false), directIoAlign(0)
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
    if ((dataMode & USE_DATA) && (dataMode & (USE_FREAD | USE_DIRECT_IO)) == 0) {
        //Debug(Debug::INFO) << "Touch data file " << dataFileName << "\n";
        // the budget has to cover every file, or a multi file db clears it once per file and touches all of them
        if (Util::canTouchMemory(totalDataSize) == false) {
            Debug(Debug::WARNING) << "Can not touch " << totalDataSize << " into main memory\n";
            return;
        }
        for(size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++){
            size_t dataSize = dataSizeOffset[fileIdx+1]-dataSizeOffset[fileIdx];
            magicBytes += Util::touchMemory(dataFiles[fileIdx], dataSize);
        }

    }
}

template <typename T>
void DBReader<T>::mlock(){
    if ((dataMode & USE_DATA) && (dataMode & USE_DIRECT_IO) == 0) {
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
    if (dataFileName != NULL) {
        dbtype = FileUtil::parseDbType(dataFileName);
    }
    if ((dataMode & USE_DIRECT_IO) && (dataMode & USE_WRITABLE)) {
        // writes would land in a bounce buffer and be dropped on the next read
        Debug(Debug::ERROR) << "USE_WRITABLE cannot be combined with USE_DIRECT_IO\n";
        EXIT(EXIT_FAILURE);
    }
    if (dataMode & USE_DATA) {
        dataFileNames = FileUtil::findDatafiles(dataFileName);
        if (dataFileNames.empty()) {
            Debug(Debug::ERROR) << "No datafile could be found for " << dataFileName << "!\n";
            EXIT(EXIT_FAILURE);
        }
        resolveIoPolicy(accessType);
        dataFileCnt = dataFileNames.size();
        dataSizeOffset = new size_t[dataFileNames.size() + 1];
        dataFiles = new char*[dataFileNames.size()];
        if (dataMode & USE_DIRECT_IO) {
            openDataFds();
        } else {
            mapDataFiles();
        }
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
    if (externalData == false) {
        MemoryMapped indexData(indexFileName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!indexData.isValid()){
            Debug(Debug::ERROR) << "Cannot open index file " << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        char* indexDataChar = (char *) indexData.getData();
        size_t indexDataSize = indexData.size();
        size = Util::ompCountLines(indexDataChar, indexDataSize, threads);

        index = new(std::nothrow) Index[size];
        Util::checkAllocation(index, "Cannot allocate index memory in DBReader");
        incrementMemory(sizeof(Index) * size);

        bool isSortedById = readIndex(indexDataChar, indexDataSize, index, dataSize);
        indexData.close();

        // sortIndex also handles access modes that don't require sorting
        sortIndex(isSortedById);

        // adjacent pairs carry no loop dependency, and index[0] alone is trivially sorted
        sortedByOffset = true;
#pragma omp parallel for schedule(static) reduction(&&: sortedByOffset) num_threads(threads)
        for (size_t i = 1; i < size; i++) {
            sortedByOffset = sortedByOffset && (index[i - 1].offset <= index[i].offset);
        }
    }

    compression = isCompressed(dbtype);
    padded = (getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU);

    if (dataMode & USE_DATA) {
        // loadBatch runs inside a parallel region, so the worker vector cannot be grown lazily there
        freeIoBatch();
        ioBatch = new IoBatch();
        ioBatch->workers.resize(threads);
    }

    if ((dataMode & USE_DATA) && (dataMode & USE_DIRECT_IO)) {
        if (compression == COMPRESSED) {
            // a compressed index stores the decompressed length, so it cannot bound the read
            Debug(Debug::ERROR) << "USE_DIRECT_IO cannot read the compressed database " << dataFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        // readIndex already tracked max index[i].length in maxSeqLen, so do not walk the index again
        allocateDirectBuffers();
    }

    if(compression == COMPRESSED || padded){
        compressedBufferSizes = new size_t[threads];
        compressedBuffers = new char*[threads];
        dstream = new ZSTD_DStream*[threads];
        for(int i = 0; i < threads; i++){
            // allocated buffer, grown on demand from here by getDataCompressed and getUnpadded
            compressedBufferSizes[i] = 0;
            compressedBuffers[i] = NULL;
            const size_t wanted = std::max(static_cast<size_t>(maxSeqLen) + 2, (size_t) 1024);
            growBuffer(&compressedBuffers[i], &compressedBufferSizes[i],
                       std::min(wanted, BOUNCE_BUFFER_PREALLOC), 1, 0);
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

template<typename T>
bool DBReader<T>::isSortedByEntryLengthDescending(int threads) {
    if (size < 2) {
        return true;
    }
    bool sorted = true;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(&&:sorted)
    for (size_t i = 1; i < size; ++i) {
        sorted = sorted && (getEntryLen(i - 1) >= getEntryLen(i));
    }
    return sorted;
}

template<>
void DBReader<std::string>::sortIndex(bool isSortedById) {
    if (accessType == SORT_BY_ID){
        if (isSortedById) {
            return;
        }
        SORT_PARALLEL(index, index + size, Index::compareById);
    } else {
        if(accessType != NOSORT && accessType != HARDNOSORT){
            Debug(Debug::ERROR) << "DBReader<std::string> cannot be opened in sort mode\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

template<>
void DBReader<DBKeyType>::sortIndex(float *weights) {

    this->accessType=DBReader::SORT_BY_WEIGHTS;
    std::pair<size_t, float> *sortForMapping = new std::pair<size_t, float>[size];
    id2local = new DBLocalId[size];
    local2id = new DBLocalId[size];
    incrementMemory(sizeof(DBLocalId) * 2 * size);
    for (size_t i = 0; i < size; i++) {
        id2local[i] = i;
        local2id[i] = i;
        sortForMapping[i] = std::make_pair(i, weights[i]);
    }
    //this sort has to be stable to assure same clustering results
    SORT_PARALLEL(sortForMapping, sortForMapping + size, comparePairByWeight());
    for (size_t i = 0; i < size; i++) {
        id2local[sortForMapping[i].first] = i;
        local2id[i] = sortForMapping[i].first;
    }
    delete[] sortForMapping;
}

template<>
void DBReader<DBKeyType>::sortIndex(bool isSortedById) {

    // First, we sort the index by IDs and we keep track of the original
    // ordering in mappingToOriginalIndex array
    size_t* mappingToOriginalIndex=NULL;
    if (accessType == SORT_BY_LINE) {
        mappingToOriginalIndex = new size_t[size];
    }
    
    if ((isSortedById == false) && (accessType != HARDNOSORT) && (accessType != SORT_BY_OFFSET)) {
        // create an array of the joint original indeces --> this will be sorted:
        // permutation of 0..size-1; DBLocalId keeps this 4 bytes/entry in the default build.
        DBLocalId *sortedIndices = new DBLocalId[size];
        for (size_t i = 0; i < size; ++i) {
            sortedIndices[i] = i;
        }
        // sort sortedIndices based on index.id:
        SORT_PARALLEL(sortedIndices, sortedIndices + size, sortIndecesById(index));

        // re-order will destroy sortedIndices so copy it, if needed:
        if (accessType == SORT_BY_LINE) {
            for (size_t i = 0; i < size; ++i) {
                mappingToOriginalIndex[i] = sortedIndices[i];
            }
        }

        // re-order in-place according to sortedIndices (ruined in the process)
        // based on: https://stackoverflow.com/questions/7365814/in-place-array-reordering
        Index indexAndOffsetBuff;

        for (size_t i = 0; i < size; i++) {
            // fill buffers with what will be overwritten:
            indexAndOffsetBuff.id = index[i].id;
            indexAndOffsetBuff.offset = index[i].offset;
            indexAndOffsetBuff.length = index[i].length;

            size_t j = i;
            while (1) {
                // The inner loop won't re-process already processed elements
                size_t k = sortedIndices[j];
                sortedIndices[j] = j; // mutating sortedIndices in the process
                if (k == i) {
                    break;
                }
                // overwite at destination place:
                index[j].id = index[k].id;
                index[j].offset = index[k].offset;
                index[j].length = index[k].length;
                // re-write what was overwritten at its destination: 
                j = k;
                index[j].id = indexAndOffsetBuff.id;
                index[j].offset = indexAndOffsetBuff.offset;
                index[j].length = indexAndOffsetBuff.length;
            }
        }
        delete[] sortedIndices;
    } else if (accessType == SORT_BY_LINE) {
        for (size_t i = 0; i < size; ++i) {
            mappingToOriginalIndex[i] = i;
        }
    }
    if (accessType == SORT_BY_LENGTH) {
        
        // sort the entries by the length of the sequences
        std::pair<size_t, unsigned int> *sortForMapping = new std::pair<size_t, unsigned int>[size];
        id2local = new DBLocalId[size];
        local2id = new DBLocalId[size];
        incrementMemory(sizeof(DBLocalId) * 2 * size);
        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, index[i].length);
        }
        //this sort has to be stable to assure same clustering results
        SORT_PARALLEL(sortForMapping, sortForMapping + size, comparePairBySeqLength());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
        }
        delete[] sortForMapping;
    } else if (accessType == SHUFFLE) {
        size_t *tmpIndex = new size_t[size];
        for (size_t i = 0; i < size; i++) {
            tmpIndex[i] = i;
        }

        std::mt19937 rnd(0);
        std::shuffle(tmpIndex, tmpIndex + size, rnd);

        id2local = new DBLocalId[size];
        local2id = new DBLocalId[size];
        incrementMemory(sizeof(DBLocalId) * 2 * size);

        for (size_t i = 0; i < size; i++) {
            id2local[tmpIndex[i]] = i;
            local2id[i] = tmpIndex[i];
        }
        delete[] tmpIndex;

    } else if (accessType == LINEAR_ACCCESS) {
        // do not sort if its already in correct order
        bool isSortedByOffset = true;
        size_t prevOffset = index[0].offset;
        for (size_t i = 0; i < size; i++) {
            isSortedByOffset &= (prevOffset <= index[i].offset);
            prevOffset = index[i].offset;
        }
        if(isSortedByOffset == true && isSortedById == true){
            accessType = NOSORT;
            return;
        }

        // sort the entries by the offset of the sequences
        std::pair<size_t, size_t> *sortForMapping = new std::pair<size_t, size_t>[size];
        id2local = new DBLocalId[size];
        local2id = new DBLocalId[size];
        incrementMemory(sizeof(DBLocalId) * 2 * size);

        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, index[i].offset);
        }
        SORT_PARALLEL(sortForMapping, sortForMapping + size, comparePairByOffset());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
        }
        delete[] sortForMapping;
    } else if (accessType == SORT_BY_ID_OFFSET) {
        // sort the entries by the offset of the sequences
        std::pair<size_t, Index> *sortForMapping = new std::pair<size_t, Index>[size];
        id2local = new DBLocalId[size];
        local2id = new DBLocalId[size];
        incrementMemory(sizeof(DBLocalId) * 2 * size);

        for (size_t i = 0; i < size; i++) {
            id2local[i] = i;
            local2id[i] = i;
            sortForMapping[i] = std::make_pair(i, index[i]);
        }
        SORT_PARALLEL(sortForMapping, sortForMapping + size, comparePairByIdAndOffset());
        for (size_t i = 0; i < size; i++) {
            id2local[sortForMapping[i].first] = i;
            local2id[i] = sortForMapping[i].first;
        }
        delete[] sortForMapping;
    } else if (accessType == SORT_BY_LINE) {
        // sort the entries by the original line number in the index file
        id2local = new DBLocalId[size];
        local2id = new DBLocalId[size];
        incrementMemory(sizeof(DBLocalId) * 2 * size);

        for (size_t i = 0; i < size; i++) {
            id2local[i] = mappingToOriginalIndex[i];
            local2id[mappingToOriginalIndex[i]] = i;
        }
    } else if (accessType == SORT_BY_OFFSET) {
        // sort index based on index.offset (no id sorting):
        SORT_PARALLEL(index, index + size, Index::compareByOffset);
    }
    if (mappingToOriginalIndex) {
        delete [] mappingToOriginalIndex;
    }
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

template <typename T>
void DBReader<T>::growBuffer(char** buffer, size_t* capacity, size_t needed, size_t alignment, size_t keepBytes) {
    if (needed <= *capacity) {
        return;
    }
    // doubling bounds the number of allocations, posix_memalign is slow under contention
    size_t newCapacity = std::max(needed, *capacity * 2);
    if (alignment > 1) {
        newCapacity = ((newCapacity + alignment - 1) / alignment) * alignment;
    }
    char* mem;
    if (alignment > 1) {
        void* aligned = NULL;
        if (posix_memalign(&aligned, alignment, newCapacity) != 0) {
            aligned = NULL;
        }
        mem = static_cast<char*>(aligned);
    } else {
        mem = static_cast<char*>(malloc(newCapacity));
    }
    if (mem == NULL) {
        Debug(Debug::ERROR) << "Cannot allocate " << newCapacity << " byte buffer in DBReader!\n";
        EXIT(EXIT_FAILURE);
    }
    if (keepBytes > 0) {
        memcpy(mem, *buffer, keepBytes);
    }
    size_t freed = 0;
    if (*buffer != NULL) {
        free(*buffer);
        freed = *capacity;
    }
    // threads grow their own buffers concurrently, so the shared counter needs an atomic update here
    __sync_fetch_and_add(&totalMemorySizeInst, newCapacity - freed);
    *buffer = mem;
    *capacity = newCapacity;
}

// one descriptor per data file; totalDataSize is recomputed so a mid-run switch keeps the index offsets
template <typename T> void DBReader<T>::openDataFds() {
    totalDataSize = 0;
    dataFds = new int[dataFileNames.size()];
    for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
        // buffered descriptors have no alignment constraint, so read exactly the entry bytes
        directIoAlign = ioBufferedBatch
            ? 1
            : std::max(directIoAlign, FileUtil::getDirectIoAlignment(dataFileNames[fileIdx]));
        size_t fileSize;
        dataFds[fileIdx] = openDirect(dataFileNames[fileIdx].c_str(), &fileSize);
        dataFiles[fileIdx] = NULL;
        dataSizeOffset[fileIdx] = totalDataSize;
        totalDataSize += fileSize;
    }
    dataSizeOffset[dataFileNames.size()] = totalDataSize;
}

// keeps a dup'd descriptor when cache advice was asked for: fadvise needs one that outlives the mapping
template <typename T> void DBReader<T>::mapDataFiles() {
    const bool keepCacheAdviceFd = ioCacheAdvice
                                   && isCompressed(dbtype) != COMPRESSED
                                   && (dataMode & USE_FREAD) == 0;
    if (keepCacheAdviceFd && dataFds == NULL) {
        dataFds = new int[dataFileNames.size()];
        for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
            dataFds[fileIdx] = -1;
        }
    }
    totalDataSize = 0;
    for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
        FILE *dataFile = fopen(dataFileNames[fileIdx].c_str(), "r");
        if (dataFile == NULL) {
            const int errsv = errno;
            Debug(Debug::ERROR) << "Cannot open data file " << dataFileNames[fileIdx]
                                << ". Error " << errsv << " (" << strerror(errsv) << ").\n";
            EXIT(EXIT_FAILURE);
        }
        size_t fileSize;
        dataFiles[fileIdx] = mmapData(dataFile, &fileSize);
        if (keepCacheAdviceFd) {
            dataFds[fileIdx] = ::dup(fileno(dataFile));
            if (dataFds[fileIdx] < 0) {
                const int errsv = errno;
                Debug(Debug::ERROR) << "Cannot duplicate data fd for " << dataFileNames[fileIdx]
                                    << ". Error " << errsv << ".\n";
                EXIT(EXIT_FAILURE);
            }
        }
        dataSizeOffset[fileIdx] = totalDataSize;
        totalDataSize += fileSize;
        if (fclose(dataFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << dataFileNames[fileIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    dataSizeOffset[dataFileNames.size()] = totalDataSize;
}

// one bounce buffer per thread, wide enough for the widest aligned read a single entry can need
template <typename T> void DBReader<T>::allocateDirectBuffers() {
    if (directBuffers != NULL) {
        return;
    }
    const size_t maxRead = ((static_cast<size_t>(maxSeqLen) + 2 * directIoAlign - 1) / directIoAlign) * directIoAlign;
    const size_t initialSize = std::min(maxRead, BOUNCE_BUFFER_PREALLOC);
    directBuffers = new DirectBuffer[threads];
    for (int i = 0; i < threads; i++) {
        directBuffers[i].buffer = NULL;
        directBuffers[i].size = 0;
        directBuffers[i].file = 0;
        directBuffers[i].offset = 0;
        directBuffers[i].length = 0;
        growBuffer(&directBuffers[i].buffer, &directBuffers[i].size, initialSize, directIoAlign, 0);
    }
}

template <typename T> void DBReader<T>::freeDirectBuffers() {
    if (directBuffers == NULL) {
        return;
    }
    for (int i = 0; i < threads; i++) {
        free(directBuffers[i].buffer);
        decrementMemory(directBuffers[i].size);
    }
    delete[] directBuffers;
    directBuffers = NULL;
}

template <typename T> bool DBReader<T>::useDescriptorIo(bool keepPageCache) {
    // set before the switch: openDirect reads it to decide whether the descriptors carry O_DIRECT
    ioBufferedBatch = keepPageCache;
    const bool direct = true;
    if ((dataMode & USE_DATA) == 0 || closed == 1 || dataMapped == false) {
        return false;
    }
    if (direct == ((dataMode & USE_DIRECT_IO) != 0)) {
        return true;
    }
    // these hand back reader-owned buffers, which the descriptor path has no way to serve
    if (direct && (compression == COMPRESSED || padded || (dataMode & (USE_WRITABLE | USE_FREAD)))) {
        return false;
    }
    unmapData();
    if (direct) {
        // unmapData leaves cache-advice descriptors alone, and openDataFds allocates its own array
        if (dataFds != NULL) {
            for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
                if (dataFds[fileIdx] >= 0) {
                    ::close(dataFds[fileIdx]);
                }
            }
            delete[] dataFds;
            dataFds = NULL;
        }
        dataMode |= USE_DIRECT_IO;
        openDataFds();
        allocateDirectBuffers();
    } else {
        dataMode &= ~USE_DIRECT_IO;
        freeDirectBuffers();
        mapDataFiles();
    }
    // the batch arenas were aligned for the path that is going away
    freeIoBatch();
    ioBatch = new IoBatch();
    ioBatch->workers.resize(threads);
    dataMapped = true;
    Debug(Debug::INFO) << "IO path for " << dataFileName << " switched to "
                       << (direct ? (ioBufferedBatch ? "buffered descriptors" : "O_DIRECT") : "mmap")
                       << "\n";
    return true;
}

// mmap is the default; MMSEQS_IO_POLICY pins every eligible reader.
template <typename T> void DBReader<T>::resolveIoPolicy(int) {
    if (dataMode & USE_DIRECT_IO) {
        return;
    }
    // compressed, padded, writable and fread readers hand back reader-owned buffers the direct path cannot serve
    if (isCompressed(dbtype) == COMPRESSED
        || (getExtendedDbtype(dbtype) & Parameters::DBTYPE_EXTENDED_GPU)
        || (dataMode & (USE_WRITABLE | USE_FREAD))) {
        return;
    }
}

template <typename T> int DBReader<T>::openDirect(const char *fileName, size_t *dataSize) {
#if defined(O_DIRECT)
    int fd = ioBufferedBatch ? ::open(fileName, O_RDONLY) : ::open(fileName, O_RDONLY | O_DIRECT);
    if (fd < 0 && (errno == EINVAL || errno == ENOTSUP || errno == EOPNOTSUPP)) {
        // tmpfs and some network filesystems reject O_DIRECT, the reads stay correct without it
        fd = ::open(fileName, O_RDONLY);
    }
#else
    int fd = ::open(fileName, O_RDONLY);
#endif
    if (fd < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Cannot open data file " << fileName << " for fd IO. Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
#ifdef HAVE_POSIX_FADVISE
    if (ioBufferedBatch) {
        // the batched path is deliberately random, so do not spend cache on speculative neighbours
        posix_fadvise(fd, 0, 0, POSIX_FADV_RANDOM);
    }
#endif
#if !defined(O_DIRECT) && defined(F_NOCACHE)
    if (ioBufferedBatch == false) {
        fcntl(fd, F_NOCACHE, 1);
    }
#endif
    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << fileName << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    *dataSize = sb.st_size;
    return fd;
}

// raw syscalls, not liburing, so the build gains no dependency; a refused ring falls back to pread
#if defined(__linux__) && defined(HAVE_IO_URING)
struct DBReaderRing {
    int fd;
    unsigned *sqTail, *sqMask, *sqArray;
    unsigned *cqHead, *cqTail, *cqMask;
    struct io_uring_sqe *sqes;
    struct io_uring_cqe *cqes;
    void *sqPtr; size_t sqSize;
    void *cqPtr; size_t cqSize;
    void *sqePtr; size_t sqeSize;
    unsigned entries;
};

static bool dbReaderRingInit(DBReaderRing &r, unsigned entries) {
    struct io_uring_params p;
    memset(&p, 0, sizeof(p));
    const int fd = syscall(__NR_io_uring_setup, entries, &p);
    if (fd < 0) {
        return false;
    }
    // IORING_OP_READ is 5.6+, so a 5.1 ring that accepts setup still has to use pread
    const size_t probeSize = sizeof(struct io_uring_probe)
                           + (IORING_OP_READ + 1) * sizeof(struct io_uring_probe_op);
    struct io_uring_probe *probe = (struct io_uring_probe *) calloc(1, probeSize);
    if (probe == NULL) { ::close(fd); return false; }
    const bool readSupported =
        syscall(__NR_io_uring_register, fd, IORING_REGISTER_PROBE, probe, IORING_OP_READ + 1) >= 0
        && probe->ops_len > IORING_OP_READ
        && (probe->ops[IORING_OP_READ].flags & IO_URING_OP_SUPPORTED) != 0;
    free(probe);
    if (readSupported == false) { ::close(fd); return false; }
    size_t sqSize = p.sq_off.array + p.sq_entries * sizeof(unsigned);
    size_t cqSize = p.cq_off.cqes + p.cq_entries * sizeof(struct io_uring_cqe);
    if (p.features & IORING_FEAT_SINGLE_MMAP) {
        sqSize = (cqSize > sqSize) ? cqSize : sqSize;
        cqSize = sqSize;
    }
    void *sq = mmap(NULL, sqSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd, IORING_OFF_SQ_RING);
    if (sq == MAP_FAILED) { ::close(fd); return false; }
    void *cq = sq;
    if ((p.features & IORING_FEAT_SINGLE_MMAP) == 0) {
        cq = mmap(NULL, cqSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd, IORING_OFF_CQ_RING);
        if (cq == MAP_FAILED) { munmap(sq, sqSize); ::close(fd); return false; }
    }
    const size_t sqeSize = p.sq_entries * sizeof(struct io_uring_sqe);
    void *sqes = mmap(NULL, sqeSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_POPULATE, fd, IORING_OFF_SQES);
    if (sqes == MAP_FAILED) {
        if (cq != sq) { munmap(cq, cqSize); }
        munmap(sq, sqSize);
        ::close(fd);
        return false;
    }
    r.fd = fd;
    r.sqPtr = sq;   r.sqSize = sqSize;
    r.cqPtr = cq;   r.cqSize = cqSize;
    r.sqePtr = sqes; r.sqeSize = sqeSize;
    r.entries = p.sq_entries;
    r.sqTail  = (unsigned *) ((char *) sq + p.sq_off.tail);
    r.sqMask  = (unsigned *) ((char *) sq + p.sq_off.ring_mask);
    r.sqArray = (unsigned *) ((char *) sq + p.sq_off.array);
    r.cqHead  = (unsigned *) ((char *) cq + p.cq_off.head);
    r.cqTail  = (unsigned *) ((char *) cq + p.cq_off.tail);
    r.cqMask  = (unsigned *) ((char *) cq + p.cq_off.ring_mask);
    r.cqes    = (struct io_uring_cqe *) ((char *) cq + p.cq_off.cqes);
    r.sqes    = (struct io_uring_sqe *) sqes;
    return true;
}

static void dbReaderRingFree(DBReaderRing &r) {
    munmap(r.sqePtr, r.sqeSize);
    if (r.cqPtr != r.sqPtr) { munmap(r.cqPtr, r.cqSize); }
    munmap(r.sqPtr, r.sqSize);
    ::close(r.fd);
}
#endif

// a lane is submitted whole, so its reads fly while the caller works on the other lane
static const unsigned DBREADER_BATCH_DEPTH = 1024;
// one lane's arena, so a batch of a few hundred short entries never has to be split
static const size_t DBREADER_BATCH_ARENA = 1 * 1024 * 1024;

template <typename T>
struct DBReader<T>::IoBatch {
    struct Slot {
        size_t arenaOffset;
        size_t delta;
        size_t avail;
        size_t id;
    };
    // without coalescing, 20 short entries in one 4 KiB block become 20 O_DIRECT reads of that block
    struct Request {
        size_t arenaOffset;
        size_t file;
        size_t offset;
        size_t length;
        size_t required;
    };
    struct Lane {
        char *arena;
        size_t capacity;
        std::vector<Slot> slots;
        std::vector<Request> requests;
#if defined(__linux__) && defined(HAVE_IO_URING)
        DBReaderRing ring;
#endif
        bool ringReady;
        bool ringProbed;
        // requests handed to the kernel, completions reaped, and how many are between the two
        size_t queued;
        size_t done;
        unsigned inflight;
        Lane() : arena(NULL), capacity(0), ringReady(false), ringProbed(false), queued(0), done(0), inflight(0) {}
    };
    struct Worker {
        Lane lane[BATCH_LANES];
    };
    std::vector<Worker> workers;
};

template <typename T> void DBReader<T>::freeIoBatch() {
    if (ioBatch == NULL) {
        return;
    }
    for (size_t i = 0; i < ioBatch->workers.size(); i++) {
        for (unsigned int l = 0; l < BATCH_LANES; l++) {
            // a read still in flight would land in the arena after it is freed
            awaitBatch(i, l);
            typename IoBatch::Lane &lane = ioBatch->workers[i].lane[l];
            if (lane.arena != NULL) {
                free(lane.arena);
                decrementMemory(lane.capacity);
            }
#if defined(__linux__) && defined(HAVE_IO_URING)
            if (lane.ringReady) {
                dbReaderRingFree(lane.ring);
            }
#endif
        }
    }
    delete ioBatch;
    ioBatch = NULL;
}

template <typename T>
size_t DBReader<T>::startBatch(const size_t *ids, size_t n, unsigned int thrIdx, unsigned int laneIdx) {
    if (n == 0) {
        return 0;
    }
    if (static_cast<int>(thrIdx) >= threads || laneIdx >= BATCH_LANES) {
        Debug(Debug::ERROR) << "startBatch: thread " << thrIdx << " lane " << laneIdx
                            << " is outside " << threads << " threads with " << BATCH_LANES << " lanes\n";
        EXIT(EXIT_FAILURE);
    }
    // the arena is about to be rewritten, so nothing of the last batch may still be landing in it
    awaitBatch(thrIdx, laneIdx);
    // open() builds this, because startBatch is called from inside a parallel region
    typename IoBatch::Lane &lane = ioBatch->workers[thrIdx].lane[laneIdx];
    lane.slots.clear();
    lane.requests.clear();
    lane.queued = 0;
    lane.done = 0;
    // on mmap the bytes are already addressable, so the batch is only a list of ids
    if ((dataMode & USE_DIRECT_IO) == 0) {
        for (size_t i = 0; i < n; i++) {
            typename IoBatch::Slot slot;
            slot.arenaOffset = 0;
            slot.delta = 0;
            slot.avail = 0;
            slot.id = ids[i];
            lane.slots.push_back(slot);
        }
        return n;
    }
    return startBatchDirect(ids, n, thrIdx, laneIdx);
}

template <typename T>
void DBReader<T>::awaitBatch(unsigned int thrIdx, unsigned int laneIdx) {
    typename IoBatch::Lane &lane = ioBatch->workers[thrIdx].lane[laneIdx];
    if (lane.done < lane.requests.size()) {
        pumpBatch(thrIdx, laneIdx, true);
    }
}

template <typename T>
size_t DBReader<T>::startBatchDirect(const size_t *ids, size_t n, unsigned int thrIdx, unsigned int laneIdx) {
    typename IoBatch::Lane &lane = ioBatch->workers[thrIdx].lane[laneIdx];
    // buffered descriptors need no alignment, but posix_memalign still demands a pointer-sized one
    const size_t arenaAlign = std::max(directIoAlign, sizeof(void *));
    if (lane.arena == NULL) {
        void *mem = NULL;
        if (posix_memalign(&mem, arenaAlign, DBREADER_BATCH_ARENA) != 0 || mem == NULL) {
            Debug(Debug::ERROR) << "Cannot allocate " << DBREADER_BATCH_ARENA << " byte batch arena\n";
            EXIT(EXIT_FAILURE);
        }
        lane.arena = static_cast<char *>(mem);
        lane.capacity = DBREADER_BATCH_ARENA;
        incrementMemory(lane.capacity);
    }

    // one Slot per input id even after merging spans, so batchAt keeps the caller's order
    size_t used = 0;
    size_t taken = 0;
    for (; taken < n; taken++) {
        const size_t id = ids[taken];
        if (id >= size) {
            Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << indexFileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t localId = (local2id != NULL) ? local2id[id] : id;
        const size_t offset = index[localId].offset;
        const size_t length = index[localId].length;
        const size_t file = fileIdxByOffset(offset);
        const size_t fileOffset = offset - dataSizeOffset[file];
        const size_t alignedOffset = fileOffset & ~(directIoAlign - 1);
        const size_t delta = fileOffset - alignedOffset;
        const size_t readLen = (delta + length + directIoAlign - 1) & ~(directIoAlign - 1);
        const size_t readEnd = alignedOffset + readLen;

        bool merged = false;
        if (lane.requests.empty() == false) {
            typename IoBatch::Request &request = lane.requests.back();
            const size_t requestEnd = request.offset + request.length;
            // Only merge forward. Arbitrary callers still work; they simply get one request per id.
            if (request.file == file && alignedOffset >= request.offset && alignedOffset <= requestEnd) {
                const size_t mergedEnd = std::max(requestEnd, readEnd);
                const size_t extra = mergedEnd - requestEnd;
                if (used + extra <= lane.capacity) {
                    request.length += extra;
                    request.required = std::max(request.required, fileOffset + length - request.offset);
                    used += extra;

                    typename IoBatch::Slot slot;
                    slot.arenaOffset = request.arenaOffset;
                    slot.delta = fileOffset - request.offset;
                    slot.avail = length;
                    slot.id = id;
                    lane.slots.push_back(slot);
                    merged = true;
                }
            }
        }
        if (merged) {
            continue;
        }

        if (used + readLen > lane.capacity) {
            if (taken > 0) {
                break;
            }
            // one entry wider than the arena still has to be read, so grow to exactly it
            free(lane.arena);
            decrementMemory(lane.capacity);
            void *mem = NULL;
            if (posix_memalign(&mem, arenaAlign, readLen) != 0 || mem == NULL) {
                Debug(Debug::ERROR) << "Cannot allocate " << readLen << " byte batch arena\n";
                EXIT(EXIT_FAILURE);
            }
            lane.arena = static_cast<char *>(mem);
            lane.capacity = readLen;
            incrementMemory(lane.capacity);
        }

        typename IoBatch::Request request;
        request.arenaOffset = used;
        request.file = file;
        request.offset = alignedOffset;
        request.length = readLen;
        request.required = delta + length;
        lane.requests.push_back(request);

        typename IoBatch::Slot slot;
        slot.arenaOffset = used;
        slot.delta = delta;
        slot.avail = length;
        slot.id = id;
        lane.slots.push_back(slot);
        used += readLen;
    }

    const size_t requestCount = lane.requests.size();
#if defined(__linux__) && defined(HAVE_IO_URING)
    if (lane.ringProbed == false) {
        lane.ringProbed = true;
        lane.ringReady = dbReaderRingInit(lane.ring, DBREADER_BATCH_DEPTH);
    }
    if (lane.ringReady) {
        pumpBatch(thrIdx, laneIdx, false);
        return taken;
    }
#endif
    // near EOF a short read is fine if it covers every slot in the span
    for (size_t i = 0; i < requestCount; i++) {
        const typename IoBatch::Request &request = lane.requests[i];
        ssize_t got;
        do {
            got = pread(dataFds[request.file], lane.arena + request.arenaOffset,
                        request.length, request.offset);
        } while (got < 0 && errno == EINTR);
        if (got < 0 || static_cast<size_t>(got) < request.required) {
            Debug(Debug::ERROR) << "Failed batch read of " << request.required << " byte from "
                                << dataFileName << " at offset " << request.offset
                                << ". Error " << errno << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    lane.done = requestCount;
    return taken;
}

template <typename T>
void DBReader<T>::pumpBatch(unsigned int thrIdx, unsigned int laneIdx, bool untilDone) {
#if defined(__linux__) && defined(HAVE_IO_URING)
    typename IoBatch::Lane &lane = ioBatch->workers[thrIdx].lane[laneIdx];
    DBReaderRing &r = lane.ring;
    const unsigned sqMask = *r.sqMask;
    const unsigned cqMask = *r.cqMask;
    const size_t requestCount = lane.requests.size();
    do {
        unsigned queued = 0;
        while (lane.inflight + queued < r.entries && lane.queued < requestCount) {
            const typename IoBatch::Request &request = lane.requests[lane.queued];
            const unsigned sqeIdx = static_cast<unsigned>(lane.queued % r.entries);
            struct io_uring_sqe *sqe = &r.sqes[sqeIdx];
            memset(sqe, 0, sizeof(*sqe));
            sqe->opcode = IORING_OP_READ;
            sqe->fd = dataFds[request.file];
            sqe->off = request.offset;
            sqe->addr = (uint64_t) (uintptr_t) (lane.arena + request.arenaOffset);
            sqe->len = static_cast<unsigned>(request.length);
            sqe->user_data = lane.queued;
            const unsigned tail = *r.sqTail;
            r.sqArray[tail & sqMask] = sqeIdx;
            __atomic_store_n(r.sqTail, tail + 1, __ATOMIC_RELEASE);
            queued++;
            lane.queued++;
        }
        unsigned waitFor = 0;
        if (untilDone) {
            waitFor = (lane.queued >= requestCount) ? (lane.inflight + queued) : 1u;
        }
        long ret;
        do {
            ret = syscall(__NR_io_uring_enter, r.fd, queued, waitFor,
                          waitFor > 0 ? IORING_ENTER_GETEVENTS : 0u, NULL, 0);
        } while (ret < 0 && errno == EINTR);
        if (ret < 0) {
            Debug(Debug::ERROR) << "io_uring_enter failed for " << dataFileName << ". Error " << errno << ".\n";
            EXIT(EXIT_FAILURE);
        }
        lane.inflight += queued;
        unsigned head = *r.cqHead;
        const unsigned cqTail = __atomic_load_n(r.cqTail, __ATOMIC_ACQUIRE);
        while (head != cqTail) {
            struct io_uring_cqe *cqe = &r.cqes[head & cqMask];
            const size_t requestIdx = (cqe->res < 0) ? 0 : static_cast<size_t>(cqe->user_data);
            const typename IoBatch::Request &request = lane.requests[requestIdx];
            if (cqe->res < 0) {
                Debug(Debug::ERROR) << "Failed to read from " << dataFileName << ". Error " << -cqe->res << ".\n";
                EXIT(EXIT_FAILURE);
            }
            if (static_cast<size_t>(cqe->res) < request.required) {
                Debug(Debug::ERROR) << "Short batch read of " << cqe->res << " byte from "
                                    << dataFileName << " at offset " << request.offset << "\n";
                EXIT(EXIT_FAILURE);
            }
            head++;
            lane.inflight--;
            lane.done++;
        }
        __atomic_store_n(r.cqHead, head, __ATOMIC_RELEASE);
    } while (untilDone && lane.done < requestCount);
#else
    (void) thrIdx;
    (void) laneIdx;
    (void) untilDone;
#endif
}

template <typename T>
const char *DBReader<T>::batchAt(unsigned int thrIdx, unsigned int laneIdx, size_t k) {
    const typename IoBatch::Lane &lane = ioBatch->workers[thrIdx].lane[laneIdx];
    const typename IoBatch::Slot &slot = lane.slots[k];
    if ((dataMode & USE_DIRECT_IO) == 0) {
        return getData(slot.id, thrIdx);
    }
    return lane.arena + slot.arenaOffset + slot.delta;
}



template <typename T>
void DBReader<T>::dropCacheRange(size_t beginOffset, size_t endOffset) {
#ifdef HAVE_POSIX_FADVISE
    if ((dataMode & USE_DATA) == 0 || (dataMode & USE_FREAD) || compression == COMPRESSED) {
        return;
    }
    if (dataMode & USE_WRITABLE) {
        return;
    }
    if ((dataMode & USE_DIRECT_IO) && ioBufferedBatch == false) {
        return;
    }
    const size_t pageSize = Util::getPageSize();
    for (size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++) {
        const size_t from = std::max(beginOffset, dataSizeOffset[fileIdx]);
        const size_t to = std::min(endOffset, dataSizeOffset[fileIdx + 1]);
        if (from >= to) {
            continue;
        }
        const size_t localFrom = ((from - dataSizeOffset[fileIdx]) + pageSize - 1) & ~(pageSize - 1);
        const size_t localTo = (to - dataSizeOffset[fileIdx]) & ~(pageSize - 1);
        if (localTo <= localFrom) {
            continue;
        }
        // fadvise skips pages that are still mapped, so the mapping has to go first
        if (dataFiles != NULL && dataFiles[fileIdx] != NULL) {
            ::madvise(dataFiles[fileIdx] + localFrom, localTo - localFrom, MADV_DONTNEED);
        }
        const int fd = (dataFds != NULL) ? dataFds[fileIdx] : -1;
        if (fd >= 0) {
            posix_fadvise(fd, localFrom, localTo - localFrom, POSIX_FADV_DONTNEED);
        }
    }
#else
    (void) beginOffset;
    (void) endOffset;
#endif
}

template <typename T>
void DBReader<T>::dropCacheAll() {
#ifdef HAVE_POSIX_FADVISE
    if ((dataMode & USE_DATA) == 0 || (dataMode & USE_FREAD) || compression == COMPRESSED) {
        return;
    }
    // MADV_DONTNEED on a writable MAP_PRIVATE mapping silently discards the caller's edits
    if (dataMode & USE_WRITABLE) {
        return;
    }
    // O_DIRECT never populated the cache, so there is nothing to drop
    if ((dataMode & USE_DIRECT_IO) && ioBufferedBatch == false) {
        return;
    }
    for (size_t fileIdx = 0; fileIdx < dataFileCnt; fileIdx++) {
        // fadvise skips pages that are still mapped, so the mapping has to go first
        if (dataFiles != NULL && dataFiles[fileIdx] != NULL) {
            const size_t fileSize = dataSizeOffset[fileIdx + 1] - dataSizeOffset[fileIdx];
            ::madvise(dataFiles[fileIdx], fileSize, MADV_DONTNEED);
        }
        const int fd = (dataFds != NULL) ? dataFds[fileIdx] : -1;
        if (fd >= 0) {
            posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
        }
    }
#endif
}

template <typename T> char* DBReader<T>::readDirect(size_t offset, size_t length, int thrIdx) {
    if (thrIdx < 0) {
#ifdef OPENMP
        thrIdx = omp_get_thread_num();
#else
        thrIdx = 0;
#endif
    }
    if (thrIdx >= threads) {
        Debug(Debug::ERROR) << "readDirect: thread index (" << thrIdx << ") >= threads (" << threads << ")\n";
        EXIT(EXIT_FAILURE);
    }
    if (offset >= totalDataSize) {
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "Size of data: " << totalDataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << offset << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t cnt = fileIdxByOffset(offset);
    size_t fileOffset = offset - dataSizeOffset[cnt];
    size_t alignedOffset = fileOffset & ~(directIoAlign - 1);
    size_t delta = fileOffset - alignedOffset;
    size_t readLen = (delta + length + directIoAlign - 1) & ~(directIoAlign - 1);
    DirectBuffer &buf = directBuffers[thrIdx];
    // an aligned read holds many short entries, so offset-ordered access mostly hits here
    if (buf.length > 0 && buf.file == cnt && alignedOffset >= buf.offset
        && (alignedOffset - buf.offset) + delta + length <= buf.length) {
        return buf.buffer + (alignedOffset - buf.offset) + delta;
    }
    if (readLen > buf.size) {
        // the cached block is gone with the old allocation
        buf.length = 0;
        growBuffer(&buf.buffer, &buf.size, readLen, directIoAlign, 0);
    }
    // reading past the end of the file returns a short count, the entry itself is still fully covered
    ssize_t read;
    do {
        read = pread(dataFds[cnt], buf.buffer, readLen, alignedOffset);
    } while (read < 0 && errno == EINTR);
    if (read < 0 || static_cast<size_t>(read) < delta + length) {
        int errsv = errno;
        buf.length = 0;
        Debug(Debug::ERROR) << "Failed to read " << length << " bytes at offset " << offset
                            << " from " << dataFileNames[cnt] << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    buf.file = cnt;
    buf.offset = alignedOffset;
    buf.length = static_cast<size_t>(read);
    return buf.buffer + delta;
}

template <typename T> void DBReader<T>::remapData(){
    if ((dataMode & USE_DATA) && (dataMode & (USE_FREAD | USE_DIRECT_IO)) == 0) {
        unmapData();
        for(size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++){
            FILE* dataFile = fopen(dataFileNames[fileIdx].c_str(), "r");
            if (dataFile == NULL) {
                // EMFILE and ENOENT mean very different things when a long run reopens the db often
                const int errsv = errno;
                Debug(Debug::ERROR) << "Cannot reopen data file " << dataFileNames[fileIdx]
                                    << " during remap. errno=" << errsv << " ("
                                    << strerror(errsv) << ")\n";
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
        // the cache-advice descriptors outlive unmapData on purpose, so this is where they go
        if (dataFds != NULL) {
            for (size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
                if (dataFds[fileIdx] >= 0) {
                    ::close(dataFds[fileIdx]);
                }
            }
            delete[] dataFds;
            dataFds = NULL;
        }
    }

    if (id2local != NULL) {
        delete[] id2local;
        decrementMemory(size*sizeof(DBLocalId));
    }
    if (local2id != NULL) {
        delete[] local2id;
        decrementMemory(size*sizeof(DBLocalId));
    }

    freeDirectBuffers();

    freeIoBatch();

    if(compressedBuffers){
        for(int i = 0; i < threads; i++){
            ZSTD_freeDStream(dstream[i]);
            free(compressedBuffers[i]);
            decrementMemory(compressedBufferSizes[i]);
        }
        delete [] compressedBuffers;
        delete [] compressedBufferSizes;
        delete [] dstream;
        compressedBuffers = NULL;
        compressedBufferSizes = NULL;
        dstream = NULL;
    }

    if(externalData == false) {
        delete[] index;
        decrementMemory(size*sizeof(Index));
    }
    closed = 1;
}

template <typename T> size_t DBReader<T>::bsearch(const Index * index, size_t N, T value)
{
    Index val;
    val.id = value;
    return std::upper_bound(index, index + N, val, Index::compareByIdOnly) - index;
}


template <typename T> char* DBReader<T>::getUnpadded(size_t id, int thrIdx) {
    char *data = getDataUncompressed(id, thrIdx);
    size_t seqLen = getSeqLen(id);

    static const char CODE_TO_CHAR[21] = {
            'A', /*  0 */ 'C', /*  1 */ 'D', /*  2 */
            'E', /*  3 */ 'F', /*  4 */ 'G', /*  5 */
            'H', /*  6 */ 'I', /*  7 */ 'K', /*  8 */
            'L', /*  9 */ 'M', /* 10 */ 'N', /* 11 */
            'P', /* 12 */ 'Q', /* 13 */ 'R', /* 14 */
            'S', /* 15 */ 'T', /* 16 */ 'V', /* 17 */
            'W', /* 18 */ 'Y', /* 19 */ 'X'  /* 20 */
    };

    growBuffer(&compressedBuffers[thrIdx], &compressedBufferSizes[thrIdx], seqLen + 2, 1, 0);
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
    char *data = getDataUncompressed(id, thrIdx);

    unsigned int cSize = *(reinterpret_cast<unsigned int *>(data));

    size_t totalSize = 0;
    const void *cBuff = static_cast<void *>(data + sizeof(unsigned int));
    const char *dataStart = data + sizeof(unsigned int);
    bool isCompressed = (dataStart[cSize] == 0) ? true : false;
    if(isCompressed){
        ZSTD_inBuffer input = {cBuff, cSize, 0};
        while (input.pos < input.size) {
            // the decompressed size is only known once the frame is consumed, so append and grow
            growBuffer(&compressedBuffers[thrIdx], &compressedBufferSizes[thrIdx],
                       totalSize + 1 + DECOMPRESS_MIN_ROOM, 1, totalSize);
            ZSTD_outBuffer output = {compressedBuffers[thrIdx] + totalSize,
                                     compressedBufferSizes[thrIdx] - totalSize - 1, 0};
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
        growBuffer(&compressedBuffers[thrIdx], &compressedBufferSizes[thrIdx], cSize + 1, 1, 0);
        memcpy(compressedBuffers[thrIdx], cBuff, cSize);
        compressedBuffers[thrIdx][cSize] = '\0';
    }
    return compressedBuffers[thrIdx];
}

template <typename T> size_t DBReader<T>::getAminoAcidDBSize() {
    checkClosed();
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
        return getDataUncompressed(id, thrIdx);
    }
}

template <typename T> char* DBReader<T>::getDataUncompressed(size_t id, int thrIdx){
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


    size_t localId = (local2id != NULL) ? local2id[id] : id;
    if (dataMode & USE_DIRECT_IO) {
        return readDirect(index[localId].offset, index[localId].length, thrIdx);
    }
    return getDataByOffset(index[localId].offset);
}

template <typename T> char* DBReader<T>::getDataByOffset(size_t offset) {
    if (dataMode & USE_DIRECT_IO) {
        Debug(Debug::ERROR) << "getDataByOffset is not supported in USE_DIRECT_IO mode, the entry length is unknown\n";
        EXIT(EXIT_FAILURE);
    }
    if (offset >= totalDataSize){
        Debug(Debug::ERROR) << "Invalid database read for database data file=" << dataFileName << ", database index=" << indexFileName << "\n";
        Debug(Debug::ERROR) << "Size of data: " << totalDataSize << "\n";
        Debug(Debug::ERROR) << "Requested offset: " << offset << "\n";
        EXIT(EXIT_FAILURE);
    }
    size_t cnt = fileIdxByOffset(offset);
    size_t fileOffset = offset - dataSizeOffset[cnt];
    return dataFiles[cnt]+fileOffset;
}

template <typename T>
void DBReader<T>::touchData(size_t id) {
    if((dataMode & USE_DATA) && (dataMode & (USE_FREAD | USE_DIRECT_IO)) == 0) {
        char *data = getDataUncompressed(id);
        size_t currDataOffset = getOffset(id);
        size_t nextDataOffset = findNextOffsetid(id);
        size_t dataSize = nextDataOffset-currDataOffset;
        magicBytes = Util::touchMemory(data, dataSize);
    }
}

template <typename T> char* DBReader<T>::getDataByDBKey(T dbKey, int thrIdx) {
    size_t id = getId(dbKey);
    if (id == DB_ENTRY_NOT_FOUND) {
        return NULL;
    }
    // getId returns a local id, so the offset lookup has to go through getData, not through index[id]
    return getData(id, thrIdx);
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
    size_t id = bsearch(index, size, dbKey);
    if (id2local != NULL) {
        return (id < size && index[id].id == dbKey) ? id2local[id] : DB_ENTRY_NOT_FOUND;
    }
    return (id < size && index[id].id == dbKey ) ? id : DB_ENTRY_NOT_FOUND;
}

template <typename T> size_t DBReader<T>::maxCount(char c) {
    checkClosed();

    size_t max = 0;
    // the direct io mode has no mapping to scan, so it counts per entry like the compressed mode does
    if (compression == COMPRESSED || (dataMode & USE_DIRECT_IO)) {
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
    if (dataMapped == true && (dataMode & USE_DIRECT_IO)) {
        for(size_t fileIdx = 0; fileIdx < dataFileNames.size(); fileIdx++) {
            if (::close(dataFds[fileIdx]) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << dataFileNames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        delete[] dataFds;
        dataFds = NULL;
        // the cached blocks belong to the descriptors that were just closed
        if (directBuffers != NULL) {
            for (int i = 0; i < threads; i++) {
                directBuffers[i].length = 0;
            }
        }
    } else if (dataMapped == true) {
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
            // DONTNEED on a reused address zero-fills it, so a stale pointer would corrupt rather than fault
            dataFiles[fileIdx] = NULL;
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


template<typename T>
void DBReader<T>::setSequentialAdvice() {
    if (dataMode & (USE_FREAD | USE_DIRECT_IO)) {
        return;
    }
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

    if (FileUtil::fileExists((srcDbName + ".index").c_str())) {
        FileUtil::move((srcDbName + ".index").c_str(), (dstDbName + ".index").c_str());
    }
    if (FileUtil::fileExists((srcDbName + ".dbtype").c_str())) {
        FileUtil::move((srcDbName + ".dbtype").c_str(), (dstDbName + ".dbtype").c_str());
    }
    if (FileUtil::fileExists((srcDbName + ".lookup").c_str())) {
        FileUtil::move((srcDbName + ".lookup").c_str(), (dstDbName + ".lookup").c_str());
    }
}

template<typename T>
void DBReader<T>::removeDb(const std::string &databaseName){
    std::vector<std::string> files = FileUtil::findDatafiles(databaseName.c_str());
    for (size_t i = 0; i < files.size(); ++i) {
        FileUtil::remove(files[i].c_str());
    }
    std::string index = databaseName + ".index";
    if (FileUtil::fileExists(index.c_str())) {
        FileUtil::remove(index.c_str());
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
        { DBFiles::DATA_DBTYPE,   ".dbtype"           },
        { DBFiles::HEADER,        "_h"                },
        { DBFiles::HEADER_INDEX,  "_h.index"          },
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

    for (size_t i = 0; i < ARRAY_SIZE(suffices); ++i) {
        std::string file = databaseName + suffices[i].suffix;
        if (dbFilesFlags & suffices[i].flag && FileUtil::fileExists(file.c_str())) {
            action(file, outDb + suffices[i].suffix);
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
void DBReader<T>::decomposeDomainByAminoAcid(size_t worldRank, size_t worldSize, size_t *startEntry, size_t *numEntries){
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
