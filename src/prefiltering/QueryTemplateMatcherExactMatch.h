//
// Created by mad on 5/26/15.
//

#ifndef MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
#define MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H

#include <sys/cdefs.h>
#include "QueryTemplateMatcher.h"


class QueryTemplateMatcherExactMatch : public virtual  QueryTemplateMatcher {
public:
    QueryTemplateMatcherExactMatch(BaseMatrix *m, IndexTable *indexTable,
                                                                       unsigned int *seqLens, short kmerThr,
                                                                       double kmerMatchProb, int kmerSize, size_t dbSize,
                                                                       unsigned int maxSeqLen, size_t maxHitsPerQuery);
    ~QueryTemplateMatcherExactMatch();

    // returns result for the sequence
    // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
    std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId);

protected:

    // match sequence against the IndexTable
    void match(Sequence* seq);

    // max hits per query
    size_t maxHitsPerQuery;


    // index the exact kmer
    //Indexer *idxer;

    unsigned int * foundSequences;
    const unsigned int MAX_DB_MATCHES = 10000000;

    std::pair<hit_t *, size_t> getResult(const int l, const unsigned int id,
                                                                         const unsigned short thr);

    // result hit buffer
    hit_t *resList;
};

#endif //MMSEQS_QUERYTEMPLATEMATCHEREXACTMATCH_H
