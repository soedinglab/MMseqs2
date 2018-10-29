#ifndef PREFILTERINGINDEXREADER_H
#define PREFILTERINGINDEXREADER_H

#include "BaseMatrix.h"
#include "IndexTable.h"
#include "DBReader.h"
#include <string>

struct PrefilteringIndexData {
    int kmerSize;
    int alphabetSize;
    int mask;
    int spacedKmer;
    int kmerThr;
    int seqType;
    int headers;
};


class PrefilteringIndexReader {
public:
    static const char*  CURRENT_VERSION;
    static unsigned int VERSION;
    static unsigned int ENTRIES;
    static unsigned int ENTRIESOFFSETS;
    static unsigned int MASKEDSEQINDEXDATA;
    static unsigned int UNMASKEDSEQINDEXDATA;
    static unsigned int SEQINDEXDATASIZE;
    static unsigned int SEQINDEXSEQOFFSET;
    static unsigned int ENTRIESNUM;
    static unsigned int SEQCOUNT;
    static unsigned int META;
    static unsigned int SCOREMATRIXNAME;
    static unsigned int SCOREMATRIX2MER;
    static unsigned int SCOREMATRIX3MER;
    static unsigned int DBRINDEX;
    static unsigned int HDR1INDEX;
    static unsigned int HDR1DATA;
    static unsigned int HDR2INDEX;
    static unsigned int HDR2DATA;
    static unsigned int GENERATOR;

    static bool checkIfIndexFile(DBReader<unsigned int> *reader);

    static void createIndexFile(const std::string &outDb, DBReader<unsigned int> *dbr, DBReader<unsigned int> *hdbr1,
                                DBReader<unsigned int> *hdbr2, BaseMatrix *subMat, int maxSeqLen, bool spacedKmer, bool compBiasCorrection,
                                int alphabetSize, int kmerSize, int maskMode, int kmerThr);

    static DBReader<unsigned int> *openNewHeaderReader(DBReader<unsigned int>*dbr, unsigned int headerIdx, unsigned int dataIdx, bool touch);

    static DBReader<unsigned int> *openNewReader(DBReader<unsigned int> *dbr, bool touch);

    static SequenceLookup *getMaskedSequenceLookup(DBReader<unsigned int> *dbr, bool touch);

    static SequenceLookup *getUnmaskedSequenceLookup(DBReader<unsigned int> *dbr, bool touch);

    static IndexTable *generateIndexTable(DBReader<unsigned int> *dbr, bool touch);

    static void printSummary(DBReader<unsigned int> *dbr);

    static PrefilteringIndexData getMetadata(DBReader<unsigned int> *dbr);

    static std::string getSubstitutionMatrixName(DBReader<unsigned int> *dbr);

    static ScoreMatrix *get2MerScoreMatrix(DBReader<unsigned int> *dbr, bool touch);

    static ScoreMatrix *get3MerScoreMatrix(DBReader<unsigned int> *dbr, bool touch);

    static std::string searchForIndex(const std::string &pathToDB);

private:
    static void printMeta(int *meta);
};

#endif
