#ifndef PREFILTERINGINDEXREADER_H
#define PREFILTERINGINDEXREADER_H

#include "IndexTable.h"
#include "DBReader.h"
#include <string>


struct PrefilteringIndexData {
    int kmerSize;
    int alphabetSize;
    int skip;
    int split;
    int local;
    int spacedKmer;
};


class PrefilteringIndexReader {
public:
    static const char*  CURRENT_VERSION;
    static unsigned int VERSION;
    static unsigned int ENTRIES;
    static unsigned int ENTRIESIZES;
    static unsigned int SEQINDEXDATA;
    static unsigned int SEQINDEXDATASIZE;
    static unsigned int SEQINDEXSEQSIZE;
    static unsigned int ENTRIESNUM;
    static unsigned int SEQCOUNT;
    static unsigned int META;

    static bool checkIfIndexFile(DBReader<unsigned int> *reader);

    static void createIndexFile(std::string outDb, DBReader<unsigned int> *dbr,
                                BaseMatrix * subMat, int maxSeqLen, bool spacedKmer, bool compBiasCorrection, int split, int alphabetSize, int kmerSize,
                                bool diagonalScoring, int threads);

    static DBReader<unsigned int> *openNewReader(DBReader<unsigned int> *dbr);

    static IndexTable *generateIndexTable(DBReader<unsigned int> *dbr, int split, bool diagonalScoring);

    static PrefilteringIndexData getMetadata(DBReader<unsigned int> *dbr);
};

#endif
