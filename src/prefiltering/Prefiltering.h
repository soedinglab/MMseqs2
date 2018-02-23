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
                      const std::string &resultDB, const std::string &resultDBIndex);
#endif

    bool runSplits(const std::string &queryDB, const std::string &queryDBIndex,
                   const std::string &resultDB, const std::string &resultDBIndex,
                   size_t fromSplit, size_t splitProcessCount);

    // merge file
    void mergeFiles(const std::string &outDb, const std::string &outDBIndex,
                    const std::vector<std::pair<std::string, std::string>> &splitFiles);

    // get substitution matrix
    static BaseMatrix *getSubstitutionMatrix(const std::string &scoringMatrixFile, size_t alphabetSize,
                                             float bitFactor, bool ignoreX, bool profileState);

    static void setupSplit(DBReader<unsigned int>& dbr, const int alphabetSize, const int threads,
                           const bool templateDBIsIndex, const size_t maxResListLen, const size_t memoryLimit,
                           int *kmerSize, int *split, int *splitMode);

    static int getKmerThreshold(const float sensitivity, const int querySeqType,
                                const int kmerScore, const int kmerSize);

private:
    static const size_t BUFFER_SIZE = 1000000;

    const std::string targetDB;
    const std::string targetDBIndex;
    DBReader<unsigned int> *tdbr;
    DBReader<unsigned int> *tidxdbr;

    BaseMatrix *subMat;
    ScoreMatrix *_2merSubMatrix;
    ScoreMatrix *_3merSubMatrix;
    IndexTable *indexTable;

    // parameter
    int splits;
    int kmerSize;
    bool spacedKmer;
    int alphabetSize;
    bool templateDBIsIndex;
    int maskMode;
    int splitMode;
    int kmerThr;
    std::string scoringMatrixFile;
    int targetSeqType;
    bool takeOnlyBestKmer;


    const size_t maxResListLen;
    const int kmerScore;
    const float sensitivity;
    const size_t resListOffset;
    const size_t maxSeqLen;
    int querySeqType;
    const bool diagonalScoring;
    const unsigned int minDiagScoreThr;
    const bool aaBiasCorrection;
    const float covThr;
    const int covMode;
    const bool includeIdentical;
    const bool earlyExit;
    const bool noPreload;
    const unsigned int threads;

    bool runSplit(DBReader<unsigned int> *qdbr, const std::string &resultDB, const std::string &resultDBIndex,
                  size_t split, size_t splitCount, bool sameQTDB);

    // compute kmer size and split size for index table
    static std::pair<int, int> optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr, int alphabetSize, int kmerSize, unsigned int threads);

    // estimates memory consumption while runtime
    static size_t estimateMemoryConsumption(int split, size_t dbSize, size_t resSize,
                                            size_t maxHitsPerQuery,
                                            int alphabetSize, int kmerSize,
                                            int threads);

    static size_t estimateHDDMemoryConsumption(size_t dbSize, size_t maxResListLen);

    ScoreMatrix *getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize);


    // needed for index lookup
    IndexTable *getIndexTable(int split, size_t dbFrom, size_t dbSize, unsigned int threads);

    /*
     * Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
     * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
     */
    double setKmerThreshold(IndexTable *indexTable, DBReader<unsigned int> *qdb);

    // write prefiltering to ffindex database
    void writePrefilterOutput(DBReader<unsigned int> *qdbr, DBWriter *dbWriter, unsigned int thread_idx, size_t id,
                              const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                              bool diagonalScoring, size_t resultOffsetPos, size_t maxResults);

    void printStatistics(const statistics_t &stats, std::list<int> **reslens,
                         unsigned int resLensSize, size_t empty, size_t maxResults);

    void mergeOutput(const std::string &outDb, const std::string &outDBIndex,
                     const std::vector<std::pair<std::string, std::string>> &filenames);

    bool isSameQTDB(const std::string &queryDB);

    void reopenTargetDb();

};

#endif
