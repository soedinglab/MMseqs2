//
// Created by lars on 10.06.15.
//

#ifndef MMSEQS_ALIGNMENTSYMMETRY_H
#define MMSEQS_ALIGNMENTSYMMETRY_H
#include <set>
#include <list>

#include "DBReader.h"

class AlignmentSymmetry {
public:
    static void readInData(DBReader<unsigned int>*pReader, DBReader<unsigned int>*pDBReader, unsigned int **pInt);
    static void readInData(DBReader<unsigned int>*pReader, DBReader<unsigned int>*pDBReader, unsigned int **pInt,unsigned short**elementScoreTable, int scoretype);
    static void computeOffsetTable(size_t *elementSizes, size_t dbSize);
    static size_t findMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize, int threads);
    static void addMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize,unsigned short**elementScoreTable);
    static void sortElements(unsigned int **elementLookupTable, size_t *offsets, size_t dbSize);

    template <typename T>
    static void setupLookupPointer(T *elements, T **elementLookupTable, size_t *elementOffset, size_t dbSize) {
        for(size_t i = 0; i < dbSize; i++) {
            elementLookupTable[i] = elements + elementOffset[i];
        }
    }
};
#endif //MMSEQS_ALIGNMENTSYMMETRY_H
