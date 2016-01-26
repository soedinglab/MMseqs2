//
// Created by lars on 10.06.15.
//

#ifndef MMSEQS_ALIGNMENTSYMMETRY_H
#define MMSEQS_ALIGNMENTSYMMETRY_H
#include <DBReader.h>
#include <DBWriter.h>
#include <set>
#include <list>
#include "SetElement.h"

class AlignmentSymmetry {
public:
    AlignmentSymmetry();
    static void readInData(DBReader<unsigned int>*pReader, DBReader<unsigned int>*pDBReader, unsigned int **pInt);
    static void readInData(DBReader<unsigned int>*pReader, DBReader<unsigned int>*pDBReader, unsigned int **pInt,unsigned short**elementScoreTable, int scoretype);
    static void computeOffsetTable(size_t *elementSizes, size_t dbSize);
    static void setupElementLookupPointer(unsigned int *elements, unsigned int **elementLookupTable, size_t *elementOffset,
                                          size_t dbSize);
    static void setupElementLookupPointerShort(unsigned short * elements, unsigned short ** elementLookupTable, size_t * elementOffset, size_t dbSize);
    static size_t findMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize, int threads);
    static void addMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize,unsigned short**elementScoreTable);

private:

};
#endif //MMSEQS_ALIGNMENTSYMMETRY_H
