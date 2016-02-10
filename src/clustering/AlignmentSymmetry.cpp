//
// Created by lars on 10.06.15.
//

#include "AlignmentSymmetry.h"
#include <climits>
#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "Log.h"
#include "AffinityClustering.h"
#include <algorithm>


#ifdef OPENMP
#include <omp.h>
#endif

AlignmentSymmetry::AlignmentSymmetry() { }


void AlignmentSymmetry::readInData(DBReader<unsigned int>*alnDbr, DBReader<unsigned int>*seqDbr, unsigned int **elementLookupTable) {

    size_t dbSize = seqDbr->getSize();
#pragma omp parallel for schedule(static)
    for(size_t i = 0; i < dbSize; i++) {
        Log::printProgress(i);
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        unsigned int clusterId = seqDbr->getDbKey(i);
        char *data = alnDbr->getDataByDBKey(clusterId);

        if (*data == '\0') { // check if file contains entry
            Debug(Debug::ERROR) << "ERROR: Sequence " << i
            << " does not contain any sequence for key " << clusterId
            << "!\n";

            continue;
        }
        while (*data != '\0' ) {
            char dbKey[255 + 1];
            Util::parseKey(data, dbKey);
            unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            unsigned int curr_element = seqDbr->getId(key);
            if (curr_element == UINT_MAX || curr_element > seqDbr->getSize()) {
                Debug(Debug::ERROR) << "ERROR: Element " << key
                << " contained in some alignment list, but not contained in the sequence database!\n";
                EXIT(EXIT_FAILURE);
            }
            *elementLookupTable[i] = curr_element;
            elementLookupTable[i]++;
            data = Util::skipLine(data);

        }
    }
}

void AlignmentSymmetry::readInData(DBReader<unsigned int>*alnDbr, DBReader<unsigned int>*seqDbr, unsigned int **elementLookupTable, unsigned short **elementScoreTable, int scoretype) {

    size_t dbSize = seqDbr->getSize();
#pragma omp parallel for schedule(static)
    for(size_t i = 0; i < dbSize; i++) {
        Log::printProgress(i);
        // seqDbr is descending sorted by length
        // the assumption is that clustering is B -> B (not A -> B)
        unsigned int clusterId = seqDbr->getDbKey(i);
        char *data = alnDbr->getDataByDBKey(clusterId);

        if (*data == '\0') { // check if file contains entry
            Debug(Debug::ERROR) << "ERROR: Sequence " << i
            << " does not contain any sequence for key " << clusterId
            << "!\n";

            continue;
        }
        while (*data != '\0' ) {
            char similarity[255+1];
            char dbKey[255 + 1];
            Util::parseKey(data, dbKey);
            unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);

            size_t curr_element = seqDbr->getId(key);
            if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                //column 1 = alignment score
                Util::parseByColumnNumber(data, similarity, 1);
                *elementScoreTable[i] = (unsigned short)(atof(similarity));
            }else {
                //column 2 = sequence identity
                Util::parseByColumnNumber(data, similarity, 2);
                *elementScoreTable[i] = (unsigned short)(atof(similarity)*1000);
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
}

size_t AlignmentSymmetry::findMissingLinks(unsigned int ** elementLookupTable, size_t * offsetTable, size_t dbSize, int threads) {
    // init memory for parallel merge
    unsigned short * tmpSize = new unsigned short[threads * dbSize];
    memset(tmpSize, 0, threads * dbSize * sizeof(unsigned short));
#pragma omp parallel for schedule(static)
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
                tmpSize[currElm * threads + thread_idx] += 1;
            }
        }
    }
    // merge size arrays
    size_t symmetricElementCount = 0;
    for(size_t setId = 0; setId < dbSize; setId++) {
        size_t elementSize = (offsetTable[setId + 1] - offsetTable[setId]);
        offsetTable[setId] = elementSize;
        for (int thread_idx = 0; thread_idx < threads; thread_idx++) {
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
