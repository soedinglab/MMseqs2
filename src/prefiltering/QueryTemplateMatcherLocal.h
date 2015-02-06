#ifndef QUERY_TEMPLATE_MATCHER_LOCAL_H
#define QUERY_TEMPLATE_MATCHER_LOCAL_H

#include "QueryTemplateMatcher.h"


class QueryTemplateMatcherLocal : public virtual  QueryTemplateMatcher {
    public:
        QueryTemplateMatcherLocal(BaseMatrix *m,
                IndexTable *indexTable,
                unsigned int *seqLens,
                short kmerThr,
                double kmerMatchProb,
                int kmerSize,
                int dbSize,
                bool aaBiasCorrection,
                int maxSeqLen,
                float zscoreThr);
        ~QueryTemplateMatcherLocal ();

        // returns result for the sequence
        // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
        std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId);
    
    protected:
    
        // match sequence against the IndexTable
        void match(Sequence* seq);
    

};

#endif
