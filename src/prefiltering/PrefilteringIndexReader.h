#ifndef PREFILTERINGINDEXREADER_H
#define PREFILTERINGINDEXREADER_H
#include "IndexTable.h"
#include "../commons/DBReader.h"
#include <string>


struct PrefilteringIndexData{
    int kmerSize;
    int alphabetSize;
    int skip;
    uint64_t entrieSize;
    uint64_t samplingEntrieSize;
};


class PrefilteringIndexReader {
public:
    static const char * VERSION;
    static const char * ENTRIES;
    static const char * ENTRIESIZES;
    static const char * META;
    
    static bool checkIfIndexFile(DBReader * reader);

    static void createIndexFile(std::string outDb, DBReader * dbr, Sequence * seq,
                                int splitt, int alphabetSize, int kmerSize, int skip );
    
    static DBReader * openNewReader(DBReader * dbr);
    
    static IndexTable * generateIndexTable(DBReader * dbr);
    
    static PrefilteringIndexData getMetadata(DBReader * dbr);
};
#endif