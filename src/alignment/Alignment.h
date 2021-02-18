#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "IndexReader.h"
#include "DBReader.h"
#include "Parameters.h"
#include "BaseMatrix.h"
#include "Matcher.h"

class Alignment {
public:
    Alignment(const std::string &querySeqDB,
              const std::string &targetSeqDB,
              const std::string &prefDB, const std::string &prefDBIndex,
              const std::string &outDB, const std::string &outDBIndex,
              const Parameters &par, const bool lcaAlign);
    ~Alignment();

    //Non-MPI
    void run();

    //MPI function
    void run(const unsigned int mpiRank, const unsigned int mpiNumProc);

    //Run parallel
    void run(const std::string &outDB, const std::string &outDBIndex, const size_t dbFrom, const size_t dbSize, bool merge);

    static bool checkCriteria(Matcher::result_t &res, bool isIdentity, double evalThr, double seqIdThr, int alnLenThr, int covMode, float covThr);

    static unsigned int initSWMode(unsigned int alignmentMode, float covThr, float seqIdThr);

private:
    // sequence coverage threshold
    double covThr;

    // sequence coverage threshold for canCov function (needed for realignment)
    double canCovThr;

    // query or query+target coverage mode
    const int covMode;

    // sets the mode how seq. id will be normalized
    const int seqIdMode;

    // e value threshold
    const double evalThr;

    // sequence identity threshold
    const double seqIdThr;

    // alignment length threshold
    const int alnLenThr;

    // include id
    const bool includeIdentity;

    // includes backtrace to alignment
    bool addBacktrace;

    // realign with different score matrix
    bool realign;
    float realignCov;

    bool sameQTDB;

    //to increase/decrease the threshold for finishing the alignment 
    float scoreBias;
    float realignScoreBias;
    int realignMaxSeqs;

    // keeps state of the SW alignment mode (ALIGNMENT_MODE_SCORE_ONLY, ALIGNMENT_MODE_SCORE_COV or ALIGNMENT_MODE_SCORE_COV_SEQID)
    unsigned int swMode;
    unsigned int realignSwMode;

    unsigned int lcaSwMode;

    unsigned int threads;
    unsigned int compressed;

    const std::string outDB;
    const std::string outDBIndex;

    size_t maxSeqLen;
    int querySeqType;
    int targetSeqType;
    bool compBiasCorrection;

    int altAlignment;
    int alignmentOutputMode;

    const unsigned int maxAccept;
    const unsigned int maxReject;
    const bool wrappedScoring;

    BaseMatrix *m;
    // costs to open a gap
    int gapOpen;
    // costs to extend a gap
    int gapExtend;
    // score difference to break alignment
    int zdrop;

    bool lcaAlign;

    // needed for realignment
    BaseMatrix *realign_m;

    DBReader<unsigned int> *qdbr;
    IndexReader * qDbrIdx;

    DBReader<unsigned int> *tdbr;
    IndexReader * tDbrIdx;

    DBReader<unsigned int> *prefdbr;

    bool reversePrefilterResult;

    static size_t estimateHDDMemoryConsumption(int dbSize, int maxSeqs);

    void computeAlternativeAlignment(unsigned int queryDbKey, Sequence &dbSeq,
                                     std::vector<Matcher::result_t> &vector, Matcher &matcher,
                                     float covThr, float evalThr, int swMode, int thread_idx);
};

#endif
