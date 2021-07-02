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
        indexFileName(strdup(indexFileName_)), size(0), dataFiles(NULL), dataSizeOffset(NULL), dataFileCnt(0),
        totalDataSize(0), dataSize(0), lastKey(T()), closed(1), dbtype(Parameters::DBTYPE_GENERIC_DB),
        compressedBuffers(NULL), compressedBufferSizes(NULL), index(NULL), id2local(NULL), local2id(NULL),
        dataMapped(false), accessType(0), externalData(false), didMlock(false)
{
    if (threads > 1) {
        FileUtil::fixRlimitNoFile();
    }
}

template <typename T>
DBReader<T>::DBReader(DBReader<T>::Index *index, size_t size, size_t dataSize, T lastKey,
        int dbType, unsigned int maxSeqLen, int threads) :
        threads(threads), dataMode(USE_INDEX), dataFileName(NULL), indexFileName(NULL),
        size(size), dataFiles(NULL), dataSizeOffset(NULL), dataFileCnt(0), totalDataSize(0), dataSize(dataSize), lastKey(lastKey),
        maxSeqLen(maxSeqLen), closed(1), dbtype(dbType), compressedBuffers(NULL), compressedBufferSizes(NULL), index(index), sortedByOffset(true),
        id2local(NULL), local2id(NULL), dataMapped(false), accessType(NOSORT), externalData(true), didMlock(false)
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
    if (dataFileName != NULL) {
        dbtype = FileUtil::parseDbType(dataFileName);
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

        size_t prevOffset = 0; // makes 0 or empty string
        sortedByOffset = true;
        for (size_t i = 0; i < size; i++) {
            sortedByOffset = sortedByOffset && index[i].offset >= prevOffset;
            prevOffset = index[i].offset;
        }
    }

    compression = isCompressed(dbtype);
    if(compression == COMPRESSED){
        compressedBufferSizes = new size_t[threads];
        compressedBuffers = new char*[threads];
        dstream = new ZSTD_DStream*[threads];
        for(int i = 0; i < threads; i++){
            // allocated buffer
            compressedBufferSizes[i] = std::max(maxSeqLen+1, 1024u);
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
void DBReader<T>::sortIndex(bool) {
}

template<typename T>
bool DBReader<T>::isSortedByOffset(){
    return sortedByOffset;
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
void DBReader<unsigned int>::sortIndex(bool isSortedById) {
    // First, we sort the index by IDs and we keep track of the original
    // ordering in mappingToOriginalIndex array
    size_t* mappingToOriginalIndex=NULL;
    if (accessType == SORT_BY_LINE) {
        mappingToOriginalIndex = new size_t[size];
    }
    
    if ((isSortedById == false) && (accessType != HARDNOSORT) && (accessType != SORT_BY_OFFSET)) {
        // create an array of the joint original indeces --> this will be sorted:
        unsigned int *sortedIndices = new unsigned int[size];
        for (unsigned int i = 0; i < size; ++i) {
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

        for (unsigned int i = 0; i < size; i++) {
            // fill buffers with what will be overwritten:
            indexAndOffsetBuff.id = index[i].id;
            indexAndOffsetBuff.offset = index[i].offset;
            indexAndOffsetBuff.length = index[i].length;

            unsigned int j = i;
            while (1) {
                // The inner loop won't re-process already processed elements
                unsigned int k = sortedIndices[j];
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
        std::pair<unsigned int, unsigned int> *sortForMapping = new std::pair<unsigned int, unsigned int>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        incrementMemory(sizeof(unsigned int) * 2 * size);
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

        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        incrementMemory(sizeof(unsigned int) * 2 * size);

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
        std::pair<unsigned int, size_t> *sortForMapping = new std::pair<unsigned int, size_t>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        incrementMemory(sizeof(unsigned int) * 2 * size);

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
        std::pair<unsigned int, Index> *sortForMapping = new std::pair<unsigned int, Index>[size];
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        incrementMemory(sizeof(unsigned int) * 2 * size);

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
        id2local = new unsigned int[size];
        local2id = new unsigned int[size];
        incrementMemory(sizeof(unsigned int) * 2 * size);

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

template <typename T> void DBReader<T>::remapData(){
    if ((dataMode & USE_DATA) && (dataMode & USE_FREAD) == 0) {
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
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
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
    }else{
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
        return (id != UINT_MAX) ? getDataCompressed(id, thrIdx) : NULL;
    }else{
        return (id != UINT_MAX) ? getDataByOffset(index[id].offset) : NULL;
    }
}

template <typename T> size_t DBReader<T>::getLookupSize() const {
    checkClosed();
    return lookupSize;
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

template <typename T> unsigned int DBReader<T>::getLookupFileNumber(size_t id){
    if (id >= lookupSize){
        Debug(Debug::ERROR) << "Invalid database read for id=" << id << ", database index=" << dataFileName << ".lookup\n";
        Debug(Debug::ERROR) << "getLookupFileNumber: local id (" << id << ") >= db size (" << lookupSize << ")\n";
        EXIT(EXIT_FAILURE);
    }
    return lookup[id].fileNumber;
}

template<>
void DBReader<unsigned int>::lookupEntryToBuffer(std::string& buffer, const LookupEntry& entry) {
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

template <typename T> size_t DBReader<T>::getId (T dbKey){
    size_t id = bsearch(index, size, dbKey);
    if (id2local != NULL) {
        return (id < size && index[id].id == dbKey) ? id2local[id] : UINT_MAX;
    }
    return (id < size && index[id].id == dbKey ) ? id : UINT_MAX;
}

template <typename T> size_t DBReader<T>::maxCount(char c) {
    checkClosed();

    size_t max = 0;
    if (compression == COMPRESSED) {
        size_t entries = getSize();
#ifdef OPENMP
        size_t localThreads = std::min(entries, static_cast<size_t>(threads));
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

    unsigned int localLastKey = 0;
    const unsigned int BATCH_SIZE = 1048576;
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
void DBReader<unsigned int>::readIndexId(unsigned int* id, char*, const char** cols) {
    *id = Util::fast_atoi<unsigned int>(cols[0]);
}

template<>
unsigned int DBReader<std::string>::indexIdToNum(std::string * id){
    return id->size();
}
template<>
unsigned int DBReader<unsigned int>::indexIdToNum(unsigned int * id) {
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
size_t DBReader<unsigned int>::indexMemorySize(const DBReader<unsigned int> &idx) {
    size_t memSize = // size + dataSize
            2 * sizeof(size_t)
            // maxSeqLen + lastKey + dbtype
            + 3 * sizeof(unsigned int)
            // index
            + idx.size * sizeof(DBReader<unsigned int>::Index)
            // seqLens
            + idx.size * sizeof(unsigned int);

    return memSize;
}

template <>
char* DBReader<unsigned int>::serialize(const DBReader<unsigned int> &idx) {
    char* data = (char*) malloc(indexMemorySize(idx));
    char* p = data;
    memcpy(p, &idx.size, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.dataSize, sizeof(size_t));
    p += sizeof(size_t);
    memcpy(p, &idx.lastKey, sizeof(unsigned int));
    p += sizeof(unsigned int);
    memcpy(p, &idx.dbtype, sizeof(int));
    p += sizeof(unsigned int);
    memcpy(p, &idx.maxSeqLen, sizeof(unsigned int));
    p += sizeof(unsigned int);
    memcpy(p, idx.index, idx.size * sizeof(DBReader<unsigned int>::Index));
    p += idx.size * sizeof(DBReader<unsigned int>::Index);
    return data;
}

template <>
DBReader<unsigned int> *DBReader<unsigned int>::unserialize(const char* data, int threads) {
    const char* p = data;
    size_t size = *((size_t*)p);
    p += sizeof(size_t);
    size_t dataSize = *((size_t*)p);
    p += sizeof(size_t);
    unsigned int lastKey = *((unsigned int*)p);
    p += sizeof(unsigned int);
    int dbType = *((int*)p);
    p += sizeof(int);
    unsigned int maxSeqLen = *((unsigned int*)p);
    p += sizeof(unsigned int);
    DBReader<unsigned int>::Index *idx = (DBReader<unsigned int>::Index *)p;
    p += size * sizeof(DBReader<unsigned int>::Index);

    return new DBReader<unsigned int>(idx, size, dataSize, lastKey, dbType, maxSeqLen, threads);
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
#ifdef HAVE_POSIX_MADVISE
    for(size_t i = 0; i < dataFileCnt; i++){
        size_t dataSize = dataSizeOffset[i+1] - dataSizeOffset[i];
        if (dataSize > 0 && posix_madvise (dataFiles[i], dataSize, POSIX_MADV_SEQUENTIAL) != 0){
            Debug(Debug::ERROR) << "posix_madvise returned an error " << dataFileName << "\n";
        }
    }
#endif
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
        lookup[i].fileNumber = Util::fast_atoi<size_t>(cols[2]);
        lookupData = Util::skipLine(lookupData);

        currPos = lookupData - (char *) data;

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

void copyLinkDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags, bool link) {
    if (dbFilesFlags & DBFiles::DATA) {
        std::vector<std::string> names = FileUtil::findDatafiles(databaseName.c_str());
        if (names.size() == 1) {
            if (link) {
                FileUtil::symlinkAbs(names[0].c_str(), outDb.c_str());
            } else {
                FileUtil::copyFile(names[0].c_str(), outDb.c_str());
            }
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
                if (link) {
                    FileUtil::symlinkAbs(names[i], outDb + ext);
                } else {
                    FileUtil::copyFile(names[i].c_str(), (outDb + ext).c_str());
                }
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
            if (link) {
                FileUtil::symlinkAbs(file, outDb + suffices[i].suffix);
            } else {
                FileUtil::copyFile(file.c_str(), (outDb + suffices[i].suffix).c_str());
            }
        }
    }
}


template<typename T>
void DBReader<T>::softlinkDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags) {
    copyLinkDb(databaseName, outDb, dbFilesFlags, true);
}

template<typename T>
void DBReader<T>::copyDb(const std::string &databaseName, const std::string &outDb, DBFiles::Files dbFilesFlags) {
    copyLinkDb(databaseName, outDb, dbFilesFlags, false);
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

template class DBReader<unsigned int>;
template class DBReader<std::string>;
