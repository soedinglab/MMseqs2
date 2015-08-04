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
    AlignmentSymmetry(DBReader *seqDbr, DBReader *alnDbr, DBWriter *alnWr, int threads);

    void  execute();
    static void readInData(DBReader *pReader, DBReader *pDBReader, unsigned int **pInt);

    static void computeOffsetTable(size_t *elementSizes, size_t dbSize);

    static void setupElementLookupPointer(unsigned int *elements, unsigned int **elementLookupTable, size_t *elementOffset,
                                          size_t dbSize);
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    DBWriter* alnWr;

//datastructures

    size_t dbSize;

    int threads;

    size_t findMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize);

    void addMissingLinks(unsigned int **elementLookupTable, size_t *offsetTable, size_t dbSize);

    void reconstructSet(DBReader *alnDbr, DBReader *seqDbr, DBWriter *alnWr, const size_t *elementLookupTable,
                        const size_t *pInt, unsigned int **pInt1);
};
#endif //MMSEQS_ALIGNMENTSYMMETRY_H
