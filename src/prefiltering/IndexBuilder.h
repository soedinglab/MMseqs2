#ifndef MMSEQS_INDEXBUILDER_H
#define MMSEQS_INDEXBUILDER_H

#include "IndexTable.h"
#include "ExtendedSubstitutionMatrix.h"

class IndexBuilder {
public:
    static void fillDatabase(IndexTable *indexTable, SequenceLookup **externalLookup, BaseMatrix &subMat,
                             ScoreMatrix & three,  ScoreMatrix & two, Sequence *seq,
                             DBReader<unsigned int> *dbr, size_t dbFrom, size_t dbTo, int kmerThr,
                             bool mask, bool maskLowerCaseMode, float maskProb, int maskNrepeats, int targetSearchMode);
};

#endif
