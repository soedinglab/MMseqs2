#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>

#include "DBReader.h"
#include "Matcher.h"
#include "Parameters.h"
#include "BaseMatrix.h"
#include "Sequence.h"
#include "SequenceLookup.h"

class Alignment {

public:

    Alignment(const std::string &querySeqDB, const std::string &querySeqDBIndex,
              const std::string &targetSeqDB, const std::string &targetSeqDBIndex,
              const std::string &prefDB, const std::string &prefDBIndex,
              const std::string &outDB, const std::string &outDBIndex,
              const Parameters &par);

    ~Alignment();

    //None MPI
    void run(const unsigned int maxAlnNum, const unsigned int maxRejected);

    //MPI function
    void run(const unsigned int mpiRank, const unsigned int mpiNumProc,
             const unsigned int maxAlnNum, const unsigned int maxRejected);

    //Run parallel
    void run(const std::string &outDB, const std::string &outDBIndex,
             const size_t dbFrom, const size_t dbSize,
             const unsigned int maxAlnNum, const unsigned int maxRejected);

    static void mergeAndRemoveTmpDatabases(const std::string &out, const std::string &outIndex,
                                           const std::vector<std::pair<std::string, std::string >> &vector);

private:
    // query sequence coverage threshold
    const double covThr;

    // query sequence coverage threshold
    const double targetCovThr;

    // e value threshold
    const double evalThr;

    // sequence identity threshold
    const double seqIdThr;

    // include id
    const bool includeIdentity;

    // includes backtrace to alignment
    const bool addBacktrace;

    // realign with different score matrix
    const bool realign;

    bool sameQTDB;

    // keeps state of the SW alignment mode (ALIGNMENT_MODE_SCORE_ONLY, ALIGNMENT_MODE_SCORE_COV or ALIGNMENT_MODE_SCORE_COV_SEQID)
    unsigned int swMode;
    const unsigned int threads;

    const std::string outDB;
    const std::string outDBIndex;

    const size_t maxSeqLen;
    const int querySeqType;
    const int targetSeqType;
    const bool compBiasCorrection;

    BaseMatrix *m;

    // needed for realignment
    BaseMatrix *realign_m;

    DBReader<unsigned int> *qseqdbr;
    SequenceLookup *qSeqLookup;

    DBReader<unsigned int> *tseqdbr;
    DBReader<unsigned int> *tidxdbr;
    SequenceLookup *tSeqLookup;

    DBReader<unsigned int> *prefdbr;

    bool templateDBIsIndex;

    void closeReaders();

    void setQuerySequence(Sequence &seq, size_t id, unsigned int key);

    void setTargetSequence(Sequence &seq, unsigned int key);

    size_t getQueryDbSize() const;

    size_t getTargetDbSize() const;

    size_t getTargetDbEntries() const;

    size_t getQueryDbEntries() const;
};

#endif
