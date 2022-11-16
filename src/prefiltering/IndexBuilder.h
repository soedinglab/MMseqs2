#ifndef MMSEQS_INDEXBUILDER_H
#define MMSEQS_INDEXBUILDER_H

#include "IndexTable.h"

class IndexBuilder {
public:
    static void fillDatabase(IndexTable *indexTable, SequenceLookup **maskedLookup, SequenceLookup **unmaskedLookup,
                             BaseMatrix &subMat, Sequence *seq,
                             DBReader<unsigned int> *dbr, size_t dbFrom, size_t dbTo, int kmerThr, bool mask, bool maskLowerCaseMode, float maskProb);
};

#endif
