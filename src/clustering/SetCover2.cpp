//
// Created by mad on 5/17/15.
//
#include "SetCover2.h"
#include "Debug.h"
#include "Util.h"

SetCover2::SetCover2(size_t elementSize, const unsigned int *setSizes,
                     unsigned int **inputSetsData,
                     unsigned short const **weights)
        : setSizes(setSizes), setData(inputSetsData), elementSize(elementSize) {
    weightRange = 0;

    offsetElements = new size_t[elementSize+1];

    // count occurrence of each element
    memset(offsetElements, 0, (elementSize + 1) * sizeof(size_t));
    countElementsInGraph(elementSize,(const unsigned int **) inputSetsData, setSizes, offsetElements);
    createOffsetTable(elementSize, offsetElements);

    // element size == set size
    elementCount = offsetElements[elementSize] + 1;
    elementLookupTable = new ElementLookup*[elementSize + 1];
    elementLookupData = new ElementLookup[offsetElements[elementSize] + 1];
    sets = new Set[elementSize];
    for(size_t i = 0; i < elementSize; i++) {
        sets[i].setId = i;
        sets[i].next = NULL;
        sets[i].last = NULL;
        sets[i].elements = NULL;
        sets[i].weight = 0;
    }

    // set element edge pointers by using the offset table
    for(size_t i = 0; i <= elementSize; i++) {
        elementLookupTable[i] = elementLookupData + offsetElements[i];
    }
//    Debug(Debug::INFO) << "Adding sets and edges...\n";

    for(size_t setId = 0; setId < elementSize; setId++) {
        sets[setId].elements = setData[setId];
        sets[setId].weight = setSizes[setId];
        weightRange = std::max(weightRange, sets[setId].weight);
        for(size_t elemPos = 0; elemPos < setSizes[setId]; elemPos++){
            const unsigned int elementId = setData[setId][elemPos];
            ElementLookup * elem = elementLookupTable[elementId];
            elem->setId = setId;
            elem->arrayPos = elemPos;
            elementLookupTable[elementId]++;
        }
    }

    // set element edge pointers by using the offset table
    for(size_t i = 0; i <= elementSize; i++) {
        elementLookupTable[i] = elementLookupData + offsetElements[i];
    }

    // init array by weight range
    this->orderedSetsByWeigh = new Set *[weightRange+1]; // score range
    memset(this->orderedSetsByWeigh, 0, sizeof(Set *) * (weightRange+1)); // init with 0 (no element is set)
    // add sets at weight position
    for(size_t i = 0; i < elementSize; i++) {
        addSetToWeightPosition(sets[i].weight, &sets[i]);
    }
}

SetCover2::~SetCover2(){
    delete [] sets;
    delete [] elementLookupTable;
    delete [] elementLookupData;
    delete [] offsetElements;
    delete [] orderedSetsByWeigh;
}

void SetCover2::addSetToWeightPosition(size_t weight, Set *set_to_add){
    if(weight > weightRange) {
        Debug(Debug::ERROR) << "ERROR: Weight is wrong in addSetToWeightPosition. Weight is " <<
        weight << " element_size:  " << weightRange << ".\n";
        EXIT(EXIT_FAILURE);
    }
    Set * weighted_position_start_set = this->orderedSetsByWeigh[weight];

    if(weighted_position_start_set == NULL) { // first element is not yet set
        set_to_add->next = NULL;
        set_to_add->last = NULL;
    }else{  // first element is already set
        weighted_position_start_set->last = set_to_add;
        set_to_add->next = weighted_position_start_set;
    }
    this->orderedSetsByWeigh[weight] = set_to_add;
}


void SetCover2::removeSet(Set * s) {
    const unsigned int setId = s->setId;

    unplug_set(s);
    // each element is a set
    //std::cout << setId << " " << sets[setId].weight << std::endl;
    for (size_t i = 0; i < setSizes[setId]; i++) {
        const unsigned int elementId = s->elements[i];
        if(elementId == UINT_MAX)
            continue;
        size_t elementSize = offsetElements[elementId + 1] - offsetElements[elementId];
        ElementLookup *allElements = elementLookupTable[elementId];
        for (size_t j = 0; j < elementSize; j++) {
            ElementLookup elementLookup = allElements[j];
            Set *parentSet = &sets[elementLookup.setId];

            if (parentSet->setId != setId && parentSet->elements[elementLookup.arrayPos] != UINT_MAX) {
                unplug_set(parentSet);
                //A => A, B, C
                //B => B, D
                //D => D
                //Resulting set:
                //A => A,B,C
                //B => D
                if (parentSet->setId == elementId) {
                    parentSet->weight = 0;
                    for (size_t x = 0; x < setSizes[parentSet->setId]; x++) {
                        parentSet->elements[x] = UINT_MAX;
                    }

                } else {
                    parentSet->weight -= 1;
                }

                parentSet->elements[elementLookup.arrayPos] = UINT_MAX;
                addSetToWeightPosition(parentSet->weight, parentSet);
            }
        }

    }
}

void SetCover2::unplug_set(Set * set_to_unplug){
    Set * last_set = set_to_unplug->last;
    Set * next_set = set_to_unplug->next;
    set_to_unplug->next = NULL;
    set_to_unplug->last = NULL;
    if(last_set == NULL && next_set == NULL){ // only one element left
        orderedSetsByWeigh[set_to_unplug->weight] = NULL;
    }else if(last_set == NULL){ // first set
        next_set->last = NULL;
        orderedSetsByWeigh[set_to_unplug->weight] = next_set;
    } else if (next_set == NULL) { // last set
        last_set->next = NULL;
    } else { // set in the middle
        last_set->next = next_set;
        next_set->last = last_set;
    }
}

void SetCover2::printGraph(){
    for(size_t setId = 0; setId < elementSize; setId++) {
        std::cout << setId << " ( ";
        for (size_t i = 0; i < setSizes[setId]; i++) {
            std::cout << sets[setId].elements[i] << " ";
        }
        std::cout << ") " << std::endl;
    }
}

std::list<SetCover2::Set *> SetCover2::execute_set_cover(){
    //covered_vertex is a bit set which holds the covered states for all vertexes
    std::list<SetCover2::Set *> result;
    for(size_t i = weightRange; i > 0;) {
        while(this->orderedSetsByWeigh[i] == NULL && i > 0){
            i--;
        }
        if(i == 0){ // if i == 0 than score = 0 -> no element in set anymore
            break;
        }
        Set *maxSet = this->orderedSetsByWeigh[i];

        removeSet(maxSet);
        result.push_back(maxSet); // O(1)
    }
    return result;
}

void SetCover2::countElementsInGraph(size_t elementSize, unsigned int const **graphData, unsigned int const * setSizes, size_t *countElements) {
    for(size_t i = 0; i < elementSize; i++){
        for(size_t j = 0; j < setSizes[i]; j++){
            size_t currId = graphData[i][j];
            countElements[currId]++;
        }
    }
}

void SetCover2::createOffsetTable(size_t elementSize, size_t *countElements) {
    size_t elementLenght = countElements[0];
    countElements[0] = 0;
    for(size_t i = 0; i < elementSize; i++) {
        size_t tmp = countElements[i+1];
        countElements[i+1] = countElements[i] + elementLenght;
        elementLenght = tmp;
    }
}
