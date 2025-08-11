//
// Implemented by Martin Steinegger, Lars vdd
//
#include "AlignmentSymmetry.h"
#include <climits>
#include <new>
#include <algorithm>
#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "FastSort.h"
#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

#define LEN(x, y) (x[y+1] - x[y])

void AlignmentSymmetry::readInData(DBReader<unsigned int>*alnDbr, DBReader<unsigned int>*seqDbr,
                                   unsigned int **elementLookupTable, unsigned short **elementScoreTable,
                                   int scoretype, size_t *offsets) {
    const int alnType = alnDbr->getDbtype();
    const size_t dbSize = seqDbr->getSize();
    const size_t flushSize = 1000000;
    Debug::Progress progress(dbSize);
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbSize)/static_cast<double>(flushSize)));
    for(size_t it = 0; it < iterations; it++) {
        size_t start = it * flushSize;
        size_t bucketSize = std::min(dbSize - (it * flushSize), flushSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 100)
            for (size_t i = start; i < (start + bucketSize); i++) {
                progress.updateProgress();
                // seqDbr is descending sorted by length
                // the assumption is that clustering is B -> B (not A -> B)
                const unsigned int clusterId = seqDbr->getDbKey(i);
                char *data = alnDbr->getDataByDBKey(clusterId, thread_idx);

                if (*data == '\0') { // check if file contains entry
                    elementLookupTable[i][0] = seqDbr->getId(clusterId);
                    if (elementScoreTable != NULL) {
                        if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_ALIGNMENT_RES)) {
                            if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                                //column 1 = alignment score
                                elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                            } else {
                                //column 2 = sequence identity [0-1]
                                elementScoreTable[i][0] = (unsigned short) (1.0 * 1000.0f);
                            }
                        } else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_RES) ||
                                   Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_REV_RES)) {
                            //column 1 = alignment score or sequence identity [0-100]
                            elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                        } else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_CLUSTER_RES)) {
                            elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                        }
                    }
                    continue;
                }
                size_t setSize = LEN(offsets, i);
                size_t writePos = 0;
                while (*data != '\0') {
                    if (writePos >= setSize) {
                        Debug(Debug::ERROR) << "Set " << i
                                            << " has more elements than allocated (" << setSize
                                            << ")!\n";
                        continue;
                    }
                    char similarity[255 + 1];
                    char dbKey[255 + 1];
                    Util::parseKey(data, dbKey);
                    const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                    const size_t currElement = seqDbr->getId(key);
                    if (elementScoreTable != NULL) {
                        if (Parameters::isEqualDbtype(alnType,Parameters::DBTYPE_ALIGNMENT_RES)) {
                            if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                                //column 1 = alignment score
                                Util::parseByColumnNumber(data, similarity, 1);
                                elementScoreTable[i][writePos] = (unsigned short) (atof(similarity));
                            } else {
                                //column 2 = sequence identity [0-1]
                                Util::parseByColumnNumber(data, similarity, 2);
                                elementScoreTable[i][writePos] = (unsigned short) (atof(similarity) * 1000.0f);
                            }
                        }
                        else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_RES) ||
                                 Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_REV_RES)) {
                            //column 1 = alignment score or sequence identity [0-100]
                            Util::parseByColumnNumber(data, similarity, 1);
                            short sim = atoi(similarity);
                            elementScoreTable[i][writePos] = (unsigned short) (sim >0 ? sim : -sim);
                        }
                        else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_CLUSTER_RES)) {
                            elementScoreTable[i][writePos] = (unsigned short) (USHRT_MAX);
                        }
                        else {
                            Debug(Debug::ERROR) << "Alignment format is not supported!\n";
                            EXIT(EXIT_FAILURE);
                        }
                    }
                    if (currElement == UINT_MAX || currElement > seqDbr->getSize()) {
                        Debug(Debug::ERROR) << "Element " << dbKey
                                            << " contained in some alignment list, but not contained in the sequence database!\n";
                        EXIT(EXIT_FAILURE);
                    }
                    elementLookupTable[i][writePos] = currElement;
                    writePos++;
                    data = Util::skipLine(data);
                }
            }
        }
        alnDbr->remapData();
    }
}

void AlignmentSymmetry::readInDataSet(DBReader<unsigned int>*alnDbr, DBReader<unsigned int>*seqDbr,
                                   unsigned int **elementLookupTable, unsigned short **elementScoreTable,
                                   int scoretype, size_t *offsets, size_t *sourceOffsets, unsigned int **sourceLookupTable,  unsigned int *keyToSet, bool isfirst) {
    const int alnType = alnDbr->getDbtype();
    const size_t dbSize = seqDbr->getSize();
    const size_t flushSize = 1000000;
    Debug::Progress progress(dbSize);
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbSize)/static_cast<double>(flushSize)));
    for(size_t it = 0; it < iterations; it++) {
        size_t start = it * flushSize;
        size_t bucketSize = std::min(dbSize - (it * flushSize), flushSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 100)

            for (size_t i = start; i < (start + bucketSize); i++) {
                progress.updateProgress();
                // seqDbr is descending sorted by length
                // the assumption is that clustering is B -> B (not A -> B)
                const unsigned int clusterId = seqDbr->getDbKey(i);
                size_t start1 = sourceOffsets[clusterId];
                size_t end1 = sourceOffsets[clusterId+1];
                size_t len = end1 - start1;
                size_t isnull = 0;

                size_t writePos = 0;
                std::vector<bool> bitFlags(dbSize, false);
                for (size_t j = 0; j < len; ++j) {
                    unsigned int value = sourceLookupTable[clusterId][j];
                    if (value != UINT_MAX) {
                        const size_t alnId = alnDbr->getId(value);
                        char *data = alnDbr->getData(alnId, thread_idx);
                        if (*data == '\0') { // check if file contains entry
                            isnull++;
                            continue;
                        }
                        while (*data != '\0') {
                            char similarity[255 + 1];
                            char dbKey[255 + 1];
                            Util::parseKey(data, dbKey);
                            const size_t currElement = seqDbr->getId(keyToSet[(unsigned int) strtoul(dbKey, NULL, 10)]);
                            if(bitFlags[currElement]==0){
                                if (elementScoreTable != NULL) {
                                    if (Parameters::isEqualDbtype(alnType,Parameters::DBTYPE_ALIGNMENT_RES)) {
                                        if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                                            //column 1 = alignment score
                                            Util::parseByColumnNumber(data, similarity, 1);
                                            elementScoreTable[i][writePos] = (unsigned short) (atof(similarity));
                                        } else {
                                            //column 2 = sequence identity [0-1]
                                            Util::parseByColumnNumber(data, similarity, 2);
                                            elementScoreTable[i][writePos] = (unsigned short) (atof(similarity) * 1000.0f);
                                        }
                                    }
                                    else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_RES) ||
                                            Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_REV_RES)) {
                                        //column 1 = alignment score or sequence identity [0-100]
                                        Util::parseByColumnNumber(data, similarity, 1);
                                        short sim = atoi(similarity);
                                        elementScoreTable[i][writePos] = (unsigned short) (sim >0 ? sim : -sim);
                                    }
                                    else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_CLUSTER_RES)) {
                                        elementScoreTable[i][writePos] = (unsigned short) (USHRT_MAX);
                                    }
                                    else {
                                        Debug(Debug::ERROR) << "Alignment format is not supported!\n";
                                        EXIT(EXIT_FAILURE);
                                    }
                                }
                                if (currElement == UINT_MAX || currElement > seqDbr->getSize()) {
                                    Debug(Debug::ERROR) << "Element " << dbKey
                                                        << " contained in some alignment list, but not contained in the sequence database!\n";
                                    EXIT(EXIT_FAILURE);
                                }
                                elementLookupTable[i][writePos] = currElement;
                                bitFlags[currElement] = 1;
                                writePos++;
                            }
                            data = Util::skipLine(data);
                        }
                    }
                }
                if (isfirst) {
                    offsets[i] = writePos;
                }
                if (isnull == len) {
                    elementLookupTable[i][0] = seqDbr->getId(clusterId);
                    if (elementScoreTable != NULL) {
                        if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_ALIGNMENT_RES)) {
                            if (scoretype == Parameters::APC_ALIGNMENTSCORE) {
                                //column 1 = alignment score
                                elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                            } else {
                                //column 2 = sequence identity [0-1]
                                elementScoreTable[i][0] = (unsigned short) (1.0 * 1000.0f);
                            }
                        } else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_RES) ||
                                Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_PREFILTER_REV_RES)) {
                            //column 1 = alignment score or sequence identity [0-100]
                            elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                        } else if (Parameters::isEqualDbtype(alnType, Parameters::DBTYPE_CLUSTER_RES)) {
                            elementScoreTable[i][0] = (unsigned short) (USHRT_MAX);
                        }
                    }
                    isnull = 0;
                    if (isfirst) {
                        offsets[i] = 1;
                    }
                    continue;
                }
            }
        }
        alnDbr->remapData();
    }
}

size_t AlignmentSymmetry::findMissingLinks(unsigned int ** elementLookupTable, size_t * offsetTable, size_t dbSize, int threads) {
    // init memory for parallel merge
    unsigned int * tmpSize = new(std::nothrow) unsigned int[threads * dbSize];
    Util::checkAllocation(tmpSize, "Can not allocate memory in findMissingLinks");
    memset(tmpSize, 0, static_cast<size_t>(threads) * dbSize * sizeof(unsigned int));
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 1000)
        for (size_t setId = 0; setId < dbSize; setId++) {
            const size_t elementSize = LEN(offsetTable, setId);
            for (size_t elementId = 0; elementId < elementSize; elementId++) {
                const unsigned int currElm = elementLookupTable[setId][elementId];
                const unsigned int currElementSize = LEN(offsetTable, currElm);
                const bool elementFound = std::binary_search(elementLookupTable[currElm],
                                                            elementLookupTable[currElm] + currElementSize, setId);
                // this is a new connection since setId is not contained in currentElementSet
                if (elementFound == false) {                        
                    tmpSize[static_cast<size_t>(currElm) * static_cast<size_t>(threads) +
                            static_cast<size_t>(thread_idx)] += 1;
                }
            }
        }
    }
    // merge size arrays
    size_t symmetricElementCount = 0;
    for(size_t setId = 0; setId < dbSize; setId++) {
        offsetTable[setId] = LEN(offsetTable, setId);
        for (int thread_idx = 0; thread_idx < threads; thread_idx++) {
            offsetTable[setId] += tmpSize[static_cast<size_t>(setId) * static_cast<size_t>(threads) + static_cast<size_t>(thread_idx)];
        }
        symmetricElementCount += offsetTable[setId];
    }
    computeOffsetFromCounts(offsetTable, dbSize);

    // clear memory
    delete [] tmpSize;

    return symmetricElementCount;
}

void AlignmentSymmetry::addMissingLinks(unsigned int **elementLookupTable,
                                        size_t * offsetTableWithOutNewLinks, size_t * offsetTableWithNewLinks, size_t dbSize, unsigned short **elementScoreTable) {

    // iterate over all connections and check if it exists in the corresponding set
    // if not add it
    Debug::Progress progress(dbSize);

    for(size_t setId = 0; setId < dbSize; setId++) {
        progress.updateProgress();
        const size_t oldElementSize = LEN(offsetTableWithOutNewLinks, setId);
        const size_t newElementSize = LEN(offsetTableWithNewLinks, setId);
        if(oldElementSize > newElementSize){
            Debug(Debug::ERROR) << "SetId="<< setId <<
                                   " NewElementSize("<< newElementSize <<") <"
                                   " OldElementSize(" << oldElementSize <<") in addMissingLinks";
            EXIT(EXIT_FAILURE);
        }
        for(size_t elementId = 0; elementId < oldElementSize; elementId++) {
            const unsigned int currElm = elementLookupTable[setId][elementId];
            if(currElm == UINT_MAX || currElm > dbSize){
                Debug(Debug::ERROR) << "currElm > dbSize in element list (addMissingLinks). This should not happen.\n";
                EXIT(EXIT_FAILURE);
            }
            const unsigned int oldCurrElementSize = LEN(offsetTableWithOutNewLinks, currElm);
            const unsigned int newCurrElementSize = LEN(offsetTableWithNewLinks, currElm);
            bool found = false;
            // check if setId is already in set of currElm
            for(size_t pos = 0; pos < oldCurrElementSize && found == false; pos++){
                found = (elementLookupTable[currElm][pos] == setId);
            }
            // this is a new connection
            if(found == false){ // add connection if it could not be found
                // find pos to write
                size_t pos = oldCurrElementSize;
                while( pos < newCurrElementSize && elementLookupTable[currElm][pos] != UINT_MAX ){
                    pos++;
                }
                if(pos >= newCurrElementSize){
                    Debug(Debug::ERROR) << "pos(" << pos << ") > newCurrElementSize(" << newCurrElementSize << "). This should not happen.\n";
                    EXIT(EXIT_FAILURE);
                }
                elementLookupTable[currElm][pos] = setId;
                elementScoreTable[currElm][pos]=elementScoreTable[setId][elementId];
            }
        }
    }
}

// sort each element vector for bsearch
void AlignmentSymmetry::sortElements(unsigned int **elementLookupTable, size_t *elementOffsets, size_t dbSize) {
#pragma omp parallel for schedule(dynamic, 1000)
    for (size_t i = 0; i < dbSize; i++) {
        SORT_SERIAL(elementLookupTable[i], elementLookupTable[i] + LEN(elementOffsets, i));
    }
}
#undef LEN
