#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>

#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "Parameters.h"
#include "BlastScoreUtils.h"


class Alignment {

public:

    Alignment (std::string& querySeqDB, std::string& querySeqDBIndex,
               std::string& targetSeqDB, std::string& targetSeqDBIndex,
               std::string& prefDB, std::string& prefDBIndex,
               std::string& outDB, std::string& outDBIndex,
               Parameters &par);

    ~Alignment();
    //None MPI
    void run (const unsigned int maxAlnNum, const unsigned int maxRejected);
    //MPI function
    void run (const unsigned int mpiRank, const unsigned int mpiNumProc,
              const unsigned int maxAlnNum, const unsigned int maxRejected);
    //Run parallel
    void run (const char * outDB, const char * outDBIndex,
              const size_t dbFrom, const size_t dbSize,
              const unsigned int maxAlnNum, const unsigned int maxRejected);

    static void mergeAndRemoveTmpDatabases(const std::string& out, const std::string& outIndex,
                                           const  std::vector<std::pair<std::string, std::string >>& vector);

private:

    // keeps state of the SW alignment mode (ALIGNMENT_MODE_SCORE_ONLY, ALIGNMENT_MODE_SCORE_COV or ALIGNMENT_MODE_SCORE_COV_SEQID)
    unsigned swMode;
    int threads;

    // sequence identity threshold
    double seqIdThr;

    // query sequence coverage threshold
    double covThr;

    // e value threshold
    double evalThr;

    BaseMatrix* m;

    // 2 Sequence objects for each thread: one for the query, one for the DB sequence
    Sequence** qSeqs;

    Sequence** dbSeqs;

    Matcher** matchers;

    // needed for realignment
    BaseMatrix* realign_m;
    Matcher** realigner;

    DBReader<unsigned int>* qseqdbr;

    DBReader<unsigned int>* tseqdbr;

    DBReader<unsigned int>* prefdbr;

    // buffers for the database keys (needed during the processing of the prefilterings lists)
    unsigned int* dbKeys;

    void closeReader();

    std::string outDB;
    std::string outDBIndex;


    bool sameQTDB;
    // include id
    bool includeIdentity;
    // merge fragments
    bool fragmentMerge;

    // includes backtrace to alignment
    bool addBacktrace;

    // realign with different score matrix
    bool realign;
	
};

#endif
