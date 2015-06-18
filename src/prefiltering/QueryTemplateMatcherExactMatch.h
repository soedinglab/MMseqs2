//
// Created by mad on 5/26/15.
//

#ifndef MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
#define MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H

#include <sys/cdefs.h>
#include "CountInt32Array.h"
#include "QueryTemplateMatcher.h"


class QueryTemplateMatcherExactMatch : public virtual QueryTemplateMatcher {
public:
    QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                   unsigned int *seqLens, short kmerThr,
                                   double kmerMatchProb, int kmerSize, size_t dbSize,
                                   unsigned int maxSeqLen);
    ~QueryTemplateMatcherExactMatch();

    // returns result for the sequence
    // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
    std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId);

    // find duplicates in the diagonal bins
    size_t evaluateBins(CounterResult *output);

protected:
    // match sequence against the IndexTable
    void match(Sequence* seq);
    const static unsigned int MAX_DB_MATCHES = 16777216;

    CounterResult * counterOutput;

    std::pair<hit_t *, size_t> getResult(const int l, const unsigned int id,
                                         const unsigned short thr);
    // result hit buffer
    CountInt32Array * counter;

    static bool compareCounterResult(CounterResult first, CounterResult second){
        return (first.id > second.id) ? true : false;
    }

    // result hit buffer
    hit_t *resList;

    // pointer to position to write in bin
    const static unsigned int BIN_COUNT = 16;
    unsigned int * diagonalBins[BIN_COUNT];
    const static unsigned int BIN_SIZE = MAX_DB_MATCHES / (BIN_COUNT / 2);
    unsigned int * __restrict binData;

    void reallocBinMemory(unsigned int const binCount, size_t const binSize);

    int checkForOverflow();

    void setupBinPointer();
    // the following variables are needed to calculate the Z-score computation
    double mu;
    double sqrtMu;
    // array to pre buffer diagonals
    static const int ENTRIES_BUFFER_SIZE = 131072;
    IndexEntryLocal * entriesBuffer;

    void fillDiagonals(IndexEntryLocal *pLocal, size_t pos);
};

#endif //MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
