#include "DBReader.h"
#include "FastSort.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <cstddef>
#include <random>
#include <utility>

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
