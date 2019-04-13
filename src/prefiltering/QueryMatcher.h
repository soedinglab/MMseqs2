//
// Created by mad on 5/26/15.
//

#ifndef MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
#define MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H

#include <cstdlib>
#include "itoa.h"
#include "EvalueComputation.h"
#include "CacheFriendlyOperations.h"
#include "UngappedAlignment.h"
#include "KmerGenerator.h"


struct statistics_t{
    double kmersPerPos;
    size_t dbMatches;
    size_t doubleMatches;
    size_t querySeqLen;
    size_t diagonalOverflow;
    size_t resultsPassedPrefPerSeq;
    statistics_t() : kmersPerPos(0.0) , dbMatches(0) , doubleMatches(0), querySeqLen(0), diagonalOverflow(0), resultsPassedPrefPerSeq(0) {};
    statistics_t(double kmersPerPos, size_t dbMatches,
                 size_t doubleMatches, size_t querySeqLen, size_t diagonalOverflow, size_t resultsPassedPrefPerSeq) : kmersPerPos(kmersPerPos),
                                                                                                                      dbMatches(dbMatches),
                                                                                                                      doubleMatches(doubleMatches),
                                                                                                                      querySeqLen(querySeqLen),
                                                                                                                      diagonalOverflow(diagonalOverflow),
                                                                                                                      resultsPassedPrefPerSeq(resultsPassedPrefPerSeq){};
};

struct hit_t {
    unsigned int seqId;
    int prefScore;
    unsigned short diagonal;

    static bool compareHitsByScoreAndId(hit_t first, hit_t second){
        if(first.prefScore > second.prefScore )
            return true;
        if(second.prefScore > first.prefScore )
            return false;
        if(first.seqId < second.seqId )
            return true;
        if(second.seqId < first.seqId )
            return false;
        return false;
    }


};



class QueryMatcher {
public:
    QueryMatcher(IndexTable *indexTable, SequenceLookup *sequenceLookup,
                 BaseMatrix *kmerSubMat, BaseMatrix *ungappedAlignmentSubMat,
                 unsigned int *seqLens, short kmerThr,
                 double kmerMatchProb, int kmerSize, size_t dbSize,
                 unsigned int maxSeqLen, unsigned int effectiveKmerSize,
                 size_t maxHitsPerQuery, bool aaBiasCorrection, bool diagonalScoring,
                 unsigned int minDiagScoreThr, bool takeOnlyBestKmer,size_t resListOffset);
    ~QueryMatcher();

    // returns result for the sequence
    // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
    std::pair<hit_t *, size_t>  matchQuery(Sequence * querySeq, unsigned int identityId);

    // find duplicates in the diagonal bins
    size_t evaluateBins(IndexEntryLocal **hitsByIndex, CounterResult *output,
                        size_t outputSize, unsigned short indexFrom, unsigned short indexTo, bool computeTotalScore);

    void updateScoreBins(CounterResult *result, size_t elementCount);

    // set substituion matrix for KmerGenerator
    void setProfileMatrix(ScoreMatrix **matrix){
        this->kmerGenerator->setDivideStrategy(matrix );
    }

    // set substitution matrix
    void setSubstitutionMatrix(ScoreMatrix * three, ScoreMatrix * two) {
        this->kmerGenerator->setDivideStrategy(three, two );
    }

    // get statistics
    const statistics_t * getStatistics(){
        return stats;
    }

    const static size_t SCORE_RANGE = 256;

    static unsigned int computeScoreThreshold(unsigned int * scoreSizes, size_t maxHitsPerQuery) {
        size_t foundHits = 0;
        size_t scoreThr = 0;
        for(scoreThr = SCORE_RANGE - 1; scoreThr > 0 ; scoreThr--){
            foundHits += scoreSizes[scoreThr];
            if(foundHits >= maxHitsPerQuery)
                break;
        }
        return scoreThr;
    }

// compute -log(p)
    static inline double computeLogProbability(const unsigned short rawScore, const unsigned int dbSeqLen,
                                               const double kmerMatchProb, const double logMatchProb,
                                               const double logScoreFactorial) {
        const double score = static_cast<double>(rawScore);
        const double dbSeqLenDbl = static_cast<double>(dbSeqLen);
        const double mu = kmerMatchProb * dbSeqLenDbl;
        const double mid_term = score * (logMatchProb + log(dbSeqLenDbl));
        const double first_term = -(mu * score /(score + 1));
        return first_term + mid_term - logScoreFactorial;
    }

    static hit_t parsePrefilterHit(char* data)
    {
        hit_t result;
        const char *wordCnt[255];
        size_t cols = Util::getWordsOfLine(data, wordCnt, 254);
        if (cols == 3) {
            result.seqId = Util::fast_atoi<unsigned int>(wordCnt[0]);
            result.prefScore = Util::fast_atoi<short>(wordCnt[1]);
            result.diagonal = static_cast<unsigned short>(Util::fast_atoi<short>(wordCnt[2]));
        }
        return result;
    }

    static std::vector<hit_t> parsePrefilterHits(char *data) {
        std::vector<hit_t> ret;
        while (*data != '\0') {
            hit_t result = parsePrefilterHit(data);
            ret.push_back(result);
            data = Util::skipLine(data);
        }
        return ret;
    }

    static size_t prefilterHitToBuffer(char *buff1, hit_t &h)
    {
        char * basePos = buff1;
        char * tmpBuff = Itoa::u32toa_sse2((uint32_t) h.seqId, buff1);
        *(tmpBuff-1) = '\t';
        int score = static_cast<int>(h.prefScore);
        tmpBuff = Itoa::i32toa_sse2(score, tmpBuff);
        *(tmpBuff-1) = '\t';
        int32_t diagonal = static_cast<short>(h.diagonal);
        tmpBuff = Itoa::i32toa_sse2(diagonal, tmpBuff);
        *(tmpBuff-1) = '\n';
        *(tmpBuff) = '\0';
        return tmpBuff - basePos;
    }

protected:
    const static int KMER_SCORE = 0;
    const static int UNGAPPED_DIAGONAL_SCORE = 1;

    // keeps stats for run
    statistics_t * stats;
    // scoring matrix for local amino acid bias correction
    BaseMatrix * kmerSubMat;
    // scoring matrix for ungapped alignment
    BaseMatrix * ungappedAlignmentSubMat;
    /* generates kmer lists */
    KmerGenerator * kmerGenerator;
    /* contains the sequences for a kmer */
    IndexTable * indexTable;
    // k of the k-mer
    int kmerSize;
    // local amino acid bias correction
    bool aaBiasCorrection;
    // take only best kmer
    bool takeOnlyBestKmer;
    // kmer threshold for kmer generator
    short kmerThr;

    unsigned int maxDbMatches;
    unsigned int dbSize;

    // result hit buffer
    //CacheFriendlyOperations * diagonalMatcher;
    unsigned int activeCounter;
#define CacheFriendlyOperations(x)  CacheFriendlyOperations<x> * cachedOperation##x
    CacheFriendlyOperations(2);
    CacheFriendlyOperations(4);
    CacheFriendlyOperations(8);
    CacheFriendlyOperations(16);
    CacheFriendlyOperations(32);
    CacheFriendlyOperations(64);
    CacheFriendlyOperations(128);
    CacheFriendlyOperations(256);
    CacheFriendlyOperations(512);
    CacheFriendlyOperations(1024);
    CacheFriendlyOperations(2048);
#undef CacheFriendlyOperations

    // matcher for diagonal
    UngappedAlignment *ungappedAlignment;

    // score distribution of current query
    unsigned int *scoreSizes;

    // result hit buffer
    hit_t *resList;

    // i position to hits pointer
    IndexEntryLocal **indexPointer;

    // keeps data in inner loop
    IndexEntryLocal * __restrict databaseHits;

    // evaluated bins
    CounterResult * foundDiagonals;

    // last data pointer (for overflow check)
    IndexEntryLocal * lastSequenceHit;

    // the following variables are needed to calculate the Z-score computation
    double mu;

    //log match prob (mu) of poisson distribution
    double logMatchProb;

    //pre computed score factorials
    // S_fact = score!
    double *logScoreFactorial;

    // max seq. per query
    size_t maxHitsPerQuery;

    // offset in the result list
    size_t resListOffset;

    //pointer to seqLens
    float *seqLens;

    // match sequence against the IndexTable
    size_t match(Sequence *seq, float *pDouble);

    // extract result from databaseHits
    template <int TYPE>
    std::pair<hit_t *, size_t> getResult(CounterResult * results,
                                         size_t resultSize,
                                         size_t maxHitPerQuery,
                                         const unsigned int id,
                                         const unsigned short thr,
                                         UngappedAlignment *ungappedAlignment,
                                         const int rescale);
    // compute double hits
    size_t getDoubleDiagonalMatches();

    float *compositionBias;

    // diagonal scoring active
    bool diagonalScoring;
    unsigned int minDiagScoreThr;
    // size of max diagonalMatcher result objects
    size_t counterResultSize;

    void initDiagonalMatcher(size_t dbsize, unsigned int maxDbMatches);

    void deleteDiagonalMatcher(unsigned int activeCounter);

    size_t mergeElements(bool diagonalScoring, CounterResult *foundDiagonals, size_t hitCounter);

    size_t keepMaxScoreElementOnly(CounterResult *foundDiagonals, size_t resultSize);

    size_t radixSortByScoreSize(const unsigned int *scoreSizes,
                              CounterResult *writePos, const unsigned int scoreThreshold,
                              const CounterResult *results, const size_t resultSize);

    std::pair<size_t, unsigned int> rescoreHits(Sequence * querySeq, unsigned int *scoreSizes, CounterResult *results,
                                                size_t resultSize, UngappedAlignment *align, int lowerBoundScore);
};

#endif //MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
