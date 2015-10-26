//
// Created by lars on 10.06.15.
//

#include "AlignmentSymmetry.h"
#include <Debug.h>
#include <Log.h>
#include <sys/time.h>
#include <Parameters.h>
#include "AffinityClustering.h"

AlignmentSymmetry::AlignmentSymmetry(DBReader *seqDbr, DBReader *alnDbr, DBWriter *alnWr, int threads) {
    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->threads = threads;
    this->dbSize=alnDbr->getSize();
    this->alnWr=alnWr;
}

void AlignmentSymmetry::execute() {
    //time
    struct timeval start, end;
    const char * data = alnDbr->getData();
    size_t dataSize = alnDbr->getDataSize();
    size_t elementCount = Util::count_lines(data, dataSize);
    unsigned int * elements = new unsigned int[elementCount];
    unsigned int ** elementLookupTable = new unsigned int*[dbSize];
    size_t *elementOffsets = new size_t[dbSize + 1];
    //time
     gettimeofday(&start, NULL);
    //time
#pragma omp for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        const char *clusterId = seqDbr->getDbKey(i).c_str();
        const size_t alnId = alnDbr->getId(clusterId);
        char *data = alnDbr->getData(alnId);
        size_t dataSize = alnDbr->getSeqLens()[alnId];
        size_t elementCount = Util::count_lines(data, dataSize);
        elementOffsets[i] = elementCount;
    }
    // make offset table
    computeOffsetTable(elementOffsets, dbSize);
    // set element edge pointers by using the offset table
    setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);
    // fill elements
    readInData(alnDbr, seqDbr, elementLookupTable);
    // set element edge pointers by using the offset table
    setupElementLookupPointer(elements, elementLookupTable, elementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Parallel read in: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nSort entries.\n";
    // sort each element vector for bsearch
#pragma omp parallel for schedule(dynamic, 1000)
    for(size_t i = 0; i < dbSize; i++) {
        Log::printProgress(i);
        std::sort(elementLookupTable[i], elementLookupTable[i] + (elementOffsets[i+1] - elementOffsets[i]));
    }
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for Sorting entries: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";
    size_t * newElementOffsets = new size_t[dbSize + 1];
    memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));
    size_t newElementCount = findMissingLinks(elementLookupTable, newElementOffsets, dbSize,threads);
    delete [] elements;
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for finding missing connections: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nFind missing connections.\n";

    // resize elements
    elements = new unsigned int[newElementCount];
    memset(elements, 0, sizeof( unsigned int) * (newElementCount ));
    Debug(Debug::WARNING) << "\nFound "<< newElementCount - elementCount << " new connections.\n";
    setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nReconstruct initial order.\n";
    readInData(alnDbr, seqDbr, elementLookupTable);
    // set element edge pointers by using the offset table
    setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
//time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nAdd missing connections.\n";
    addMissingLinks(elementLookupTable, elementOffsets, dbSize);
    setupElementLookupPointer(elements, elementLookupTable, newElementOffsets, dbSize);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    Debug(Debug::WARNING) << "\nReconstruct set.\n";
    reconstructSet(alnDbr, seqDbr, alnWr, elementOffsets, newElementOffsets, elementLookupTable);
    //time
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);
    //time
    delete [] elementOffsets;
    delete [] newElementOffsets;
    delete [] elementLookupTable;
    delete [] elements;
}

void AlignmentSymmetry::readInData(DBReader *alnDbr, DBReader *seqDbr, unsigned int **elementLookupTable) {

    size_t dbSize = seqDbr->getSize();
#pragma omp parallel
    {
        char * dbKey = new char[255+1];
#pragma omp for schedule(dynamic, 100)
        for(size_t i = 0; i < dbSize; i++) {
            Log::printProgress(i);
            // seqDbr is descending sorted by length
            // the assumption is that clustering is B -> B (not A -> B)
            const char *clusterId = seqDbr->getDbKey(i).c_str();
            char *data = alnDbr->getDataByDBKey(clusterId);

            if (*data == '\0') { // check if file contains entry
                Debug(Debug::ERROR) << "ERROR: Sequence " << i
                << " does not contain any sequence for key " << clusterId
                << "!\n";

                continue;
            }
            while (*data != '\0' ) {
                Util::parseKey(data, dbKey);

                size_t curr_element = seqDbr->getId(dbKey);
                if (curr_element == UINT_MAX || curr_element > seqDbr->getSize()) {
                    Debug(Debug::ERROR) << "ERROR: Element " << dbKey
                    << " contained in some alignment list, but not contained in the sequence database!\n";
                    EXIT(EXIT_FAILURE);
                }
                *elementLookupTable[i] = curr_element;
                elementLookupTable[i]++;
                data = Util::skipLine(data);

            }
        }
        delete [] dbKey;
    }
}

void AlignmentSymmetry::readInData(DBReader *alnDbr, DBReader *seqDbr, unsigned int **elementLookupTable, unsigned short **elementScoreTable, int scoretype) {

    size_t dbSize = seqDbr->getSize();
#pragma omp parallel
    {
        char * dbKey = new char[255+1];
        char *similarity = new char[255+1];
#pragma omp for schedule(dynamic, 100)
        for(size_t i = 0; i < dbSize; i++) {
            Log::printProgress(i);
            // seqDbr is descending sorted by length
            // the assumption is that clustering is B -> B (not A -> B)
            const char *clusterId = seqDbr->getDbKey(i).c_str();
            char *data = alnDbr->getDataByDBKey(clusterId);

            if (*data == '\0') { // check if file contains entry
                Debug(Debug::ERROR) << "ERROR: Sequence " << i
                << " does not contain any sequence for key " << clusterId
                << "!\n";

                continue;
            }
            while (*data != '\0' ) {
                Util::parseKey(data, dbKey);

                size_t curr_element = seqDbr->getId(dbKey);
                if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                    Util::parseByColumnNumber(data, similarity, 1); //column 1 = alignmentscore
                    *elementScoreTable[i] = (short)(atof(std::string(similarity).c_str()));
                }else {
                    Util::parseByColumnNumber(data, similarity, 2); //column 2 = sequence identity
                    *elementScoreTable[i] = (short)(atof(std::string(similarity).c_str())*1000);
                }

                elementScoreTable[i]++;
                if (curr_element == UINT_MAX || curr_element > seqDbr->getSize()) {
                    Debug(Debug::ERROR) << "ERROR: Element " << dbKey
                    << " contained in some alignment list, but not contained in the sequence database!\n";
                    EXIT(EXIT_FAILURE);
                }
                *elementLookupTable[i] = curr_element;
                elementLookupTable[i]++;
                data = Util::skipLine(data);

            }
        }
        delete [] dbKey;
    }
}

size_t AlignmentSymmetry::findMissingLinks(unsigned int ** elementLookupTable, size_t * offsetTable, size_t dbSize, int threads) {
    // init memory for parallel merge
    unsigned short * tmpSize = new unsigned short[threads * dbSize];
    memset(tmpSize, 0, threads * dbSize * sizeof(unsigned short));
#pragma omp parallel for schedule(dynamic, 1000)
    for(size_t setId = 0; setId < dbSize; setId++) {
        Log::printProgress(setId);
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        size_t elementSize = (offsetTable[setId +1] - offsetTable[setId]);
        for(size_t elementId = 0; elementId < elementSize; elementId++) {
            const unsigned int currElm = elementLookupTable[setId][elementId];
            const unsigned int currElementSize = (offsetTable[currElm +1] - offsetTable[currElm]);
            const bool found = std::binary_search(elementLookupTable[currElm], elementLookupTable[currElm] + currElementSize,
                                                  setId);
            // this is a new connection
            if(found == false){
                tmpSize[currElm * threads + thread_idx] += 1; //TODO false sharing?
            }
        }
    }
    // merge size arrays
    size_t symmetricElementCount = 0;
    for(size_t setId = 0; setId < dbSize; setId++) {
        size_t elementSize = (offsetTable[setId + 1] - offsetTable[setId]);
        offsetTable[setId] = elementSize;
        for (size_t thread_idx = 0; thread_idx < threads; thread_idx++) {
            offsetTable[setId] += tmpSize[setId * threads + thread_idx];
        }
        symmetricElementCount += offsetTable[setId];
    }
    computeOffsetTable(offsetTable, dbSize);

    // clear memory
    delete [] tmpSize;

    return symmetricElementCount;
}

void AlignmentSymmetry::computeOffsetTable(size_t *elementSizes, size_t dbSize) {
    size_t elementLenght = elementSizes[0];
    elementSizes[0] = 0;
    for(size_t i = 0; i < dbSize; i++) {
        size_t tmp = elementSizes[i+1];
        elementSizes[i+1] = elementSizes[i] + elementLenght;
        elementLenght = tmp;
    }
}

void AlignmentSymmetry::setupElementLookupPointer(unsigned int * elements, unsigned int ** elementLookupTable, size_t * elementOffset, size_t dbSize) {
    for(size_t i = 0; i < dbSize; i++) {
        elementLookupTable[i] = elements + elementOffset[i];
    }
}

void AlignmentSymmetry::setupElementLookupPointerShort(unsigned short * elements, unsigned short ** elementLookupTable, size_t * elementOffset, size_t dbSize) {
    for(size_t i = 0; i < dbSize; i++) {
        elementLookupTable[i] = elements + elementOffset[i];
    }
}

void AlignmentSymmetry::addMissingLinks(unsigned int **elementLookupTable,
                                        size_t * offsetTable, size_t dbSize) {

    // iterate over all connections and check if it exists in the corresponding set
    // if not add it
    for(size_t setId = 0; setId < dbSize; setId++) {
        Log::printProgress(setId);
        size_t elementSize = (offsetTable[setId +1] - offsetTable[setId]);
        for(size_t elementId = 0; elementId < elementSize; elementId++) {
            const unsigned int currElm = elementLookupTable[setId][elementId];
            const unsigned int currElementSize = (offsetTable[currElm +1] - offsetTable[currElm]);
            bool found = false;
            for(size_t pos = 0; pos < currElementSize && found == false; pos++){
                found = (elementLookupTable[currElm][pos] == setId);
            }
            // this is a new connection
            if(found == false){ // add connection if it could not be found
                // find pos to write
                size_t pos;
                for(pos = currElementSize; elementLookupTable[currElm][pos] != 0; pos++ );
                elementLookupTable[currElm][pos] = setId;
            }
        }
    }
}

void AlignmentSymmetry::addMissingLinks(unsigned int **elementLookupTable,
                                        size_t * offsetTable, size_t dbSize, unsigned short **elementScoreTable) {

    // iterate over all connections and check if it exists in the corresponding set
    // if not add it
    for(size_t setId = 0; setId < dbSize; setId++) {
        Log::printProgress(setId);
        size_t elementSize = (offsetTable[setId +1] - offsetTable[setId]);
        for(size_t elementId = 0; elementId < elementSize; elementId++) {
            const unsigned int currElm = elementLookupTable[setId][elementId];
            const unsigned int currElementSize = (offsetTable[currElm +1] - offsetTable[currElm]);
            bool found = false;
            for(size_t pos = 0; pos < currElementSize && found == false; pos++){
                found = (elementLookupTable[currElm][pos] == setId);
            }
            // this is a new connection
            if(found == false){ // add connection if it could not be found
                // find pos to write
                size_t pos;
                for(pos = currElementSize; elementLookupTable[currElm][pos] != 0; pos++ );
                elementLookupTable[currElm][pos] = setId;
                elementScoreTable[currElm][pos]=elementScoreTable[setId][elementId];
            }
        }
    }
}

void AlignmentSymmetry::reconstructSet(DBReader *alnDbr, DBReader *seqDbr, DBWriter *alnWr,
                                       const size_t *oldElementOffset,
                                       const size_t *newElementOffset, unsigned int **elementLookupTable) {
    size_t dbSize = seqDbr->getSize();
#pragma omp parallel
    {
        char *buffer = new char[6400000]; //6MB
//        char * dbKey = new char[255+1];

#pragma omp for schedule(dynamic, 100)
        for (size_t set_i = 0; set_i < dbSize; set_i++) {
            Log::printProgress(set_i);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // seqDbr is descending sorted by length
            // the assumption is that clustering is B -> B (not A -> B)
            const char *clusterId = seqDbr->getDbKey(set_i).c_str();
            size_t setiIdlength = strlen(clusterId);

            const size_t alnId = alnDbr->getId(clusterId);
            const char *data = alnDbr->getData(alnId);
            const size_t dataSize = alnDbr->getSeqLens()[alnId];
            memcpy(buffer, data, dataSize - 1); // -1 for the nullbyte
            const unsigned int newElementSize = (newElementOffset[set_i +1] - newElementOffset[set_i]);
            const unsigned int oldElementSize = (oldElementOffset[set_i +1] - oldElementOffset[set_i]);
            size_t dataPos = dataSize - 1;
            if(newElementSize > oldElementSize){
                for(size_t j = oldElementSize; j < newElementSize; j++){
                    const unsigned int set_j = elementLookupTable[set_i][j];
                    const char *clusterId = seqDbr->getDbKey(set_j).c_str();
                    const size_t alnId = alnDbr->getId(clusterId);
                    char * setJData = alnDbr->getData(alnId);

                    const unsigned int setjElementSize = (oldElementOffset[set_j +1] - oldElementOffset[set_j]);
                    size_t pos;
                    for(pos = 0; pos < setjElementSize; pos++){
                        if(elementLookupTable[set_j][pos] == set_i){
                            break;
                        }
                    }
                    if(pos == setjElementSize){
                        Debug(Debug::ERROR) << "ERROR: Could not find set_i=" << set_i << " in set_j=" << set_j << "\n";
//                        for(pos = 0; pos < setjElementSize; pos++) {
//                            std::cout << elementLookupTable[set_j][pos] << std::endl;
//                        }
                        EXIT(EXIT_FAILURE);
                    }
                    for(size_t skip = 0; skip < pos; skip++){
                        setJData = Util::skipLine(setJData);
                    }
                    const char * beforePosSetJData = setJData;
//                    Util::parseKey((char*)beforePosSetJData, dbKey);
//                    size_t curr_element = seqDbr->getId(dbKey);
//                    if(set_i != curr_element){
//                        Debug(Debug::ERROR) << "ERROR: id "<< curr_element << " does not match " << set_i << "\n";
////                        for(pos = 0; pos < setjElementSize; pos++) {
////                            std::cout << elementLookupTable[set_j][pos] << std::endl;
////                        }
//                    }
                    const char * afterPosSetJData = Util::skipLine(setJData);
                    const size_t entrySize = afterPosSetJData - beforePosSetJData;
                    if(dataPos > 6400000){
                        Debug(Debug::ERROR) << "ERROR: dataPos in reconstructSet (set=" << set_i << ") is " << dataPos
                        << " is bigger than  " << 6400000 << "!\n";
                        std::cout << buffer << std::endl;
                        EXIT(EXIT_FAILURE);
                    }
                    size_t clusterIdlength = strlen(clusterId);
                    memcpy(buffer + dataPos, clusterId, clusterIdlength );
                    dataPos += clusterIdlength;
                    memcpy(buffer + dataPos, beforePosSetJData + setiIdlength, entrySize -setiIdlength );
                    dataPos += entrySize - setiIdlength;
//                    memcpy(buffer + dataPos, beforePosSetJData, entrySize );
//                    dataPos += entrySize;
                }
            }else if (newElementSize < oldElementSize){
                Debug(Debug::ERROR) << "ERROR: newElementSize " << newElementSize
                << " is smaller than oldElementSize " << oldElementSize << "!\n";
                EXIT(EXIT_FAILURE);
            }

            alnWr->write(buffer, dataPos, (char *) clusterId, thread_idx);
        }
        delete[] buffer;
//        delete[] dbKey;
    }
}
