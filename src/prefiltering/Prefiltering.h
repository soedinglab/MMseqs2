#ifndef PREFILTERING_H
#define PREFILTERING_H

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexTable.h"
#include "BaseMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "QueryMatcher.h"

#include <string>
#include <list>
#include <utility>


class Prefiltering {
public:
    Prefiltering(const std::string &queryDB,
                 const std::string &queryDBIndex,
                 const std::string &targetDB,
                 const std::string &targetDBIndex,
                 const std::string &outDB,
                 const std::string &outDBIndex,
                 Parameters &par);

    ~Prefiltering();

    void run(size_t dbFrom, size_t dbSize, int splitMode, const std::string &resultDB,
             const std::string &resultDBIndex);

    void run(size_t fromSplit, size_t splits);

    void closeReader();
    // merge file
    static void mergeFiles(const std::vector<std::pair<std::string, std::string>> &splitFiles, int mode,
                           std::string outDb, std::string outDBIndex);

    static void mergeOutput(const std::string &outDb, const std::string &outDBIndex,
                            const std::vector<std::pair<std::string, std::string>> &filenames);

    IndexTable *getIndexTable(int split, size_t dbFrom, size_t dbSize, int threads); // needed for index lookup

    static IndexTable *generateIndexTable(DBReader<unsigned int> *dbr, Sequence *seq, BaseMatrix *subMat,
                                          int alphabetSize, int kmerSize, size_t dbFrom, size_t dbTo,
                                          bool diagonalScoring, int threads);

    static void fillDatabase(DBReader<unsigned int> *dbr, Sequence *seq, IndexTable *indexTable, BaseMatrix *subMat,
                             size_t dbFrom, size_t dbTo, bool diagonalScoring, int threads);

    // get substitution matrix
    static BaseMatrix *getSubstitutionMatrix(const std::string &scoringMatrixFile, int alphabetSize, float bitFactor,
                                             bool ignoreX);

    // compute kmer size and split size for index table
    static std::pair<int, int> optimizeSplit(size_t totalMemoryInByte, DBReader<unsigned int> *tdbr, int alphabetSize, int kmerSize, int threads);

    size_t getSplitMode(){
        return splitMode;
    }

    size_t setSplit(size_t split){
        return this->split = split;
    }

    int getSplit(){
        return split;
    }
private:
    static const size_t BUFFER_SIZE = 1000000;

    int threads;

    DBReader<unsigned int> *qdbr;
    DBReader<unsigned int> *tdbr;
    DBReader<unsigned int> *tidxdbr;


    Sequence **qseq;
    char *notEmpty;

    std::list<int> **reslens;
    BaseMatrix *subMat;
    ExtendedSubstitutionMatrix *_2merSubMatrix;
    ExtendedSubstitutionMatrix *_3merSubMatrix;

    std::string outDB;
    std::string outDBIndex;
    // parameter
    int kmerSize;
    const int kmerScore;
    bool spacedKmer;
    const float sensitivity;
    size_t resListOffset;
    size_t maxResListLen;
    int alphabetSize;
    size_t maxSeqLen;
    const int querySeqType;
    const int targetSeqType;
    bool templateDBIsIndex;
    const bool diagonalScoring;
    const unsigned int minDiagScoreThr;
    bool aaBiasCorrection;
    short kmerThr;
    double kmerMatchProb;
    int split;
    int splitMode;
    bool sameQTDB;
    bool includeIdentical;

    /* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
     * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
     */
    std::pair<short, double> setKmerThreshold(IndexTable *indexTable, DBReader<unsigned int> *qdbr,
                                              DBReader<unsigned int> *tdbr, const int kmerScore);

    // write prefiltering to ffindex database
    int writePrefilterOutput(DBWriter *dbWriter, int thread_idx, size_t id,
                             const std::pair<hit_t *, size_t> &prefResults, size_t seqIdOffset,
                             bool diagonalScoring, size_t resultOffsetPos);

    // init QueryTemplateMatcher
    QueryMatcher **createQueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable,
                                              unsigned int *seqLens, short kmerThr,
                                              double kmerMatchProb, int kmerSize,
                                              size_t effectiveKmerSize, size_t dbSize,
                                              bool aaBiasCorrection, bool diagonalScoring,
                                              unsigned int maxSeqLen, size_t maxHitsPerQuery);


    void printStatistics(const statistics_t &stats, size_t empty);

    statistics_t computeStatisticForKmerThreshold(IndexTable *indexTable, size_t querySetSize,
                                                  unsigned int *querySeqsIds, bool reverseQuery,
                                                  const size_t kmerThrMid);

    int getKmerThreshold(const float sensitivity, const int score);

    static size_t estimateMemoryConsumption(int split, size_t dbSize, size_t resSize, int alphabetSize, int kmerSize,
                                            int threads);

    std::string searchForIndex(const std::string &pathToDB);

};

#endif
