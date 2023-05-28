#ifndef PREFILTERING_H
#define PREFILTERING_H

#include "Parameters.h"
#include "DBReader.h"
#include "IndexTable.h"
#include "BaseMatrix.h"
#include "ScoreMatrix.h"
#include "PrefilteringIndexReader.h"
#include "QueryMatcher.h"

#include <string>
#include <list>
#include <utility>

class QueryMatcherTaxonomyHook;

struct KmerThreshold{
    int sequenceType;
    int kmerSize;
    float base;
    float sensPerStep;
};

extern std::vector<KmerThreshold> externalThreshold;


class Prefiltering {
public:
    Prefiltering(
            const std::string &queryDB,
            const std::string &queryDBIndex,
            const std::string &targetDB,
            const std::string &targetDBIndex,
            int querySeqType, int targetSeqType,
            const Parameters &par);

    ~Prefiltering();

    void runAllSplits(const std::string &resultDB, const std::string &resultDBIndex);

#ifdef HAVE_MPI
    void runMpiSplits(const std::string &resultDB, const std::string &resultDBIndex, const std::string &localTmpPath, const int runRandomId);
#endif

    int runSplits(const std::string &resultDB, const std::string &resultDBIndex, size_t fromSplit, size_t splitProcessCount, bool merge);

    // merge file
    void mergePrefilterSplits(const std::string &outDb, const std::string &outDBIndex,
                    const std::vector<std::pair<std::string, std::string>> &splitFiles);

    // get substitution matrix
    static BaseMatrix *getSubstitutionMatrix(const MultiParam<NuclAA<std::string>> &scoringMatrixFile, MultiParam<NuclAA<int>> alphabetSize, float bitFactor, bool profileState, bool isNucl);

    static void setupSplit(DBReader<unsigned int>& dbr, const int alphabetSize, const unsigned int querySeqType, const int threads,
                           const bool templateDBIsIndex, const size_t memoryLimit, const size_t qDbSize,
                           size_t& maxResListLen, int& kmerSize, int& split, int& splitMode);

    static int getKmerThreshold(const float sensitivity, const bool isProfile, const bool hasContextPseudoCnts,
                                const SeqProf<int> kmerScore, const int kmerSize);

    static void mergeTargetSplits(const std::string &outDB, const std::string &outDBIndex,
                                  const std::vector<std::pair<std::string, std::string>> &fileNames, unsigned int threads);

private:
    const std::string queryDB;
    const std::string queryDBIndex;
    const std::string targetDB;
    const std::string targetDBIndex;
    DBReader<unsigned int> *qdbr;
    DBReader<unsigned int> *tdbr;
    DBReader<unsigned int> *tidxdbr;
    bool sameQTDB;

    BaseMatrix *kmerSubMat;
    BaseMatrix *ungappedSubMat;
    ScoreMatrix _2merSubMatrix;
    ScoreMatrix _3merSubMatrix;
    IndexTable *indexTable;
    SequenceLookup *sequenceLookup;

    // parameter
    int splits;
    int kmerSize;
    std::string spacedKmerPattern;
    std::string localTmp;
    bool spacedKmer;
    int alphabetSize;
    bool templateDBIsIndex;
    int maskMode;
    int maskLowerCaseMode;
    float maskProb;
    int splitMode;
    int kmerThr;
    MultiParam<NuclAA<std::string>> scoringMatrixFile;
    MultiParam<NuclAA<std::string>> seedScoringMatrixFile;
    int targetSeqType;
    int targetSearchMode;
    bool takeOnlyBestKmer;
    size_t maxResListLen;

    const float sensitivity;
    size_t maxSeqLen;
    int querySeqType;
    const unsigned int diagonalScoring;
    const unsigned int minDiagScoreThr;
    bool aaBiasCorrection;
    float aaBiasCorrectionScale;
    const float covThr;
    const int covMode;
    const bool includeIdentical;
    int preloadMode;
    const unsigned int threads;
    int compressed;
    QueryMatcherTaxonomyHook* taxonomyHook;

    bool runSplit(const std::string &resultDB, const std::string &resultDBIndex, size_t split, bool merge);

    // compute kmer size and split size for index table
    static std::pair<int, int> optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr, int alphabetSize, int kmerSize,
                                             unsigned int querySeqType, unsigned int threads);

    // estimates memory consumption while runtime
    static size_t estimateMemoryConsumption(int split, size_t dbSize, size_t resSize,
                                            size_t maxHitsPerQuery,
                                            int alphabetSize, int kmerSize, unsigned int querySeqType,
                                            int threads);

    static size_t estimateHDDMemoryConsumption(size_t dbSize, size_t maxResListLen);

    ScoreMatrix getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize);


    // needed for index lookup
    void getIndexTable(int split, size_t dbFrom, size_t dbSize);

    void printStatistics(const statistics_t &stats, std::list<int> **reslens,
                         unsigned int resLensSize, size_t empty, size_t maxResults);

    bool isSameQTDB();
};

#endif

