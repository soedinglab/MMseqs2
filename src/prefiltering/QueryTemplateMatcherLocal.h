#ifndef QUERY_TEMPLATE_MATCHER_LOCAL_H
#define QUERY_TEMPLATE_MATCHER_LOCAL_H

#include "QueryTemplateMatcher.h"
#include "QueryScoreLocal.h"


class QueryTemplateMatcherLocal : public virtual  QueryTemplateMatcher {
    public:
        QueryTemplateMatcherLocal(BaseMatrix *m,
                IndexTable *indexTable,
                unsigned int *seqLens,
                short kmerThr,
                double kmerMatchProb,
                int kmerSize,
                size_t effectiveKmerSize,
                size_t dbSize,
                bool aaBiasCorrection,
                unsigned int maxSeqLen,
                size_t maxHitsPerQuery);
        ~QueryTemplateMatcherLocal ();

        // returns result for the sequence
        // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
        std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId);
    
    protected:
    
        // match sequence against the IndexTable
        void match(Sequence* seq);

        /* calculates the score */
        QueryScoreLocal * queryScore;

        // max hits per query
        size_t maxHitsPerQuery;

        // fast mode for extremly fast search (20.000 faster than BLAST)
        bool fastMode;
};

#endif
