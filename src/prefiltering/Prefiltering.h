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
            const Parameters &par);

    ~Prefiltering();

    void runAllSplits(const std::string &queryDB, const std::string &queryDBIndex,
                      const std::string &resultDB, const std::string &resultDBIndex);
    void runSplits(const std::string &queryDB, const std::string &queryDBIndex,
                   const std::string &resultDB, const std::string &resultDBIndex,
                   size_t fromSplit, size_t splits);

    // merge file
    void mergeFiles(const std::vector<std::pair<std::string, std::string>> &splitFiles,
                    const std::string &outDb, const std::string &outDBIndex);

    // get substitution matrix
    static BaseMatrix *getSubstitutionMatrix(const std::string &scoringMatrixFile, size_t alphabetSize, float bitFactor, bool ignoreX);

    static void fillDatabase(DBReader<unsigned int> *dbr, Sequence *seq, IndexTable *indexTable, BaseMatrix *subMat,
                             size_t dbFrom, size_t dbTo, bool diagonalScoring, unsigned int threads);

    static void setupSplit(DBReader<unsigned int>& dbr, const int alphabetSize, const int threads,
                           const bool templateDBIsIndex, int *kmerSize, int *split, int *splitMode);

#ifdef HAVE_MPI
public:
#else
private:
#endif
int split;
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
    int kmerSize;
    bool spacedKmer;
    int alphabetSize;
    bool templateDBIsIndex;
    short kmerThr;
    double kmerMatchProb;
    int splitMode;
    std::string scoringMatrixFile;

    size_t maxResListLen;
    const int kmerScore;
    const float sensitivity;
    const size_t resListOffset;
    const size_t maxSeqLen;
    const int querySeqType;
    const int targetSeqType;
    const bool diagonalScoring;
    const unsigned int minDiagScoreThr;
    const bool aaBiasCorrection;
    const float covThr;
    const bool includeIdentical;
    const unsigned int threads;

    void runSplit(DBReader<unsigned int> *qdbr,
                  const std::string &resultDB, const std::string &resultDBIndex,
                  size_t dbFrom, size_t dbSize, bool sameQTDB);

    static IndexTable *generateIndexTable(DBReader<unsigned int> *dbr, Sequence *seq, BaseMatrix *subMat,
                                          int alphabetSize, int kmerSize, size_t dbFrom, size_t dbTo,
                                          bool diagonalScoring, unsigned int threads);

    // compute kmer size and split size for index table
    static std::pair<int, int> optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr, int alphabetSize, int kmerSize, unsigned int threads);

    // estimates memory consumption while runtime
    static size_t estimateMemoryConsumption(int split, size_t dbSize, size_t resSize, int alphabetSize, int kmerSize, unsigned int threads);

    ScoreMatrix *getScoreMatrix(const BaseMatrix& matrix, const size_t kmerSize);


    IndexTable *getIndexTable(int split, size_t dbFrom, size_t dbSize, unsigned int threads); // needed for index lookup

    /* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
     * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
     */
    std::pair<short, double> setKmerThreshold(IndexTable *indexTable, DBReader<unsigned int> *qdbr, Sequence **qseq,
                                              DBReader<unsigned int> *tdbr, int kmerScore,
                                              unsigned int effectiveKmerSize);

    // write prefiltering to ffindex database
    void writePrefilterOutput(DBReader<unsigned int> *qdbr, DBWriter *dbWriter, unsigned int thread_idx, size_t id,
                              const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                              bool diagonalScoring, size_t resultOffsetPos);

    // init QueryTemplateMatcher
    QueryMatcher **createQueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable, Sequence** qseq,
                                              unsigned int *seqLens, short kmerThr,
                                              double kmerMatchProb, int kmerSize,
                                              unsigned int effectiveKmerSize, size_t dbSize,
                                              bool aaBiasCorrection, bool diagonalScoring,
                                              unsigned int maxSeqLen, size_t maxHitsPerQuery);


    void printStatistics(const statistics_t &stats, std::list<int> **reslens, size_t empty);

    statistics_t computeStatisticForKmerThreshold(DBReader<unsigned int> *qdbr, Sequence **qseq,
                                                  IndexTable *indexTable, size_t querySetSize,
                                                  unsigned int *querySeqsIds, bool reverseQuery,
                                                  const int kmerThrMid);

    int getKmerThreshold(const float sensitivity);

    std::string searchForIndex(const std::string &pathToDB);

    void mergeOutput(const std::string &outDb, const std::string &outDBIndex,
                     const std::vector<std::pair<std::string, std::string>> &filenames);

};

#endif
