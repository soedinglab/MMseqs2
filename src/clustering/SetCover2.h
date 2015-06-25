//
// Created by mad on 5/17/15.
//

#ifndef MMSEQS_SETCOVER2_H
#define MMSEQS_SETCOVER2_H


#include <stddef.h>
#include <list>
#include <limits.h>



class SetCover2 {

// edges point from an element to a sets
// this structure must es memory efficient as possible
// since we have billions of edges
    struct __attribute__((__packed__)) ElementLookup {
        unsigned int setId;
        unsigned short arrayPos;
        ElementLookup(unsigned int id, unsigned short arrayPos): setId(id), arrayPos(arrayPos) {};
        ElementLookup(): setId(UINT_MAX), arrayPos(USHRT_MAX) {};
    };


public:
    SetCover2(size_t elementSize, const unsigned int *setSizes,
              unsigned int **inputSetsData, unsigned short const **weights);
    ~SetCover2();


    struct Set {
        size_t weight;
        Set * last;
        Set * next;
        unsigned int * elements;
        unsigned int setId;
    };

    std::list<Set *> execute_set_cover();

    void printGraph();

private:
    size_t elementCount;
    // assumption set size == element size since each element is a set (in our case)
    unsigned int elementSize;
    // the maximal weight range
    size_t weightRange;
    // all sets
    Set * sets;
    unsigned int ** setData;
    const unsigned int *setSizes;
    // element lookup
    ElementLookup ** elementLookupTable;
    ElementLookup * elementLookupData;

    size_t * offsetElements;

    // keeps the elements ordered by its Weight
    Set **orderedSetsByWeigh;
    // methodes
    void removeSet(Set * s);
    void unplug_set(Set * set_to_remove);
    void addSetToWeightPosition(size_t weight,
                                Set *set_to_add);

    void countElementsInGraph(size_t elementSize, unsigned int const **graphData, unsigned int const *setSizes,
                              size_t *countElements);

    void createOffsetTable(size_t size, size_t *pInt);


};


#endif //MMSEQS_SETCOVER2_H
