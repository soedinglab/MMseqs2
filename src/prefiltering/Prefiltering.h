#ifndef PREFILTERING_H
#define PREFILTERING_H

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexTable.h"
#include "BaseMatrix.h"
#include "ScoreMatrix.h"
#include "PrefilteringIndexReader.h"
#include "QueryMatcher.h"

#include <string>
#include <list>
#include <utility>


class Prefiltering {
public:
    Prefiltering(
            const std::string &targetDB,
            const std::string &targetDBIndex,
            int querySeqType, int targetSeqType,
            const Parameters &par);

    ~Prefiltering();

    void runAllSplits(const std::string &queryDB, const std::string &queryDBIndex,
                      const std::string &resultDB, const std::string &resultDBIndex);

#ifdef HAVE_MPI
    void runMpiSplits(const std::string &queryDB, const std::string &queryDBIndex,
                      const std::string &resultDB, const std::string &resultDBIndex,
                      const std::string &localTmpPath);
#endif

    bool runSplits(const std::string &queryDB, const std::string &queryDBIndex,
                   const std::string &resultDB, const std::string &resultDBIndex,
                   size_t fromSplit, size_t splitProcessCount, bool merge);

    // merge file
    void mergeFiles(const std::string &outDb, const std::string &outDBIndex,
                    const std::vector<std::pair<std::string, std::string>> &splitFiles);

    // get substitution matrix
    static BaseMatrix *getSubstitutionMatrix(const std::string &scoringMatrixFile, size_t alphabetSize, float bitFactor, bool profileState, bool isNucl);

    static void setupSplit(DBReader<unsigned int>& dbr, const int alphabetSize, const unsigned int querySeqType, const int threads,
                           const bool templateDBIsIndex, const size_t maxResListLen, const size_t memoryLimit,
                           int *kmerSize, int *split, int *splitMode);

    static int getKmerThreshold(const float sensitivity, const bool isProfile, const int kmerScore, const int kmerSize);

private:
    static const size_t BUFFER_SIZE = 1000000;

    const std::string targetDB;
    const std::string targetDBIndex;
    DBReader<unsigned int> *tdbr;
    DBReader<unsigned int> *tidxdbr;

    BaseMatrix *kmerSubMat;
    BaseMatrix *ungappedSubMat;
    ScoreMatrix *_2merSubMatrix;
    ScoreMatrix *_3merSubMatrix;
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
    int splitMode;
    int kmerThr;
    std::string scoringMatrixFile;
    std::string seedScoringMatrixFile;
    int targetSeqType;
    bool takeOnlyBestKmer;


    const size_t maxResListLen;
    const int kmerScore;
    const float sensitivity;
    const size_t resListOffset;
    size_t maxSeqLen;
    int querySeqType;
    const bool diagonalScoring;
    const unsigned int minDiagScoreThr;
    bool aaBiasCorrection;
    const float covThr;
    const int covMode;
    const bool includeIdentical;
    int preloadMode;
    const unsigned int threads;
    const int compressed;

    bool runSplit(DBReader<unsigned int> *qdbr, const std::string &resultDB, const std::string &resultDBIndex,
                  size_t split, size_t splitCount, bool sameQTDB, bool merge);

    // compute kmer size and split size for index table
    static std::pair<int, int> optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr, int alphabetSize, int kmerSize,
                                             unsigned int querySeqType, unsigned int threads);

    // estimates memory consumption while runtime
    static size_t estimateMemoryConsumption(int split, size_t dbSize, size_t resSize,
                                            size_t maxHitsPerQuery,
                                            int alphabetSize, int kmerSize, unsigned int querySeqType,
                                            int threads);

    static size_t estimateHDDMemoryConsumption(size_t dbSize, size_t maxResListLen);

    ScoreMatrix *getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize);


    // needed for index lookup
    void getIndexTable(int split, size_t dbFrom, size_t dbSize);

    /*
     * Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
     * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
     */
    double setKmerThreshold(DBReader<unsigned int> *qdb);

    // write prefiltering to ffindex database
    void writePrefilterOutput(DBReader<unsigned int> *qdbr, DBWriter *dbWriter, unsigned int thread_idx, size_t id,
                              const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                              size_t resultOffsetPos, size_t maxResults);

    void printStatistics(const statistics_t &stats, std::list<int> **reslens,
                         unsigned int resLensSize, size_t empty, size_t maxResults);

    void mergeOutput(const std::string &outDb, const std::string &outDBIndex,
                     const std::vector<std::pair<std::string, std::string>> &filenames);

    bool isSameQTDB(const std::string &queryDB);

    void reopenTargetDb();

};

#endif
