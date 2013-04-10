#ifndef QUERY_TEMPLATE_MATCHER_H
#define QUERY_TEMPLATE_MATCHER_H

#include <list>
#include <iostream>
#include <cstring>


#include "BaseMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Sequence.h"
#include "QueryScore.h"
#include "IndexTable.h"
#include "KmerGenerator.h"


class QueryTemplateMatcher {
    public:
        QueryTemplateMatcher (ExtendedSubstitutionMatrix* _2merSubMatrix,
                ExtendedSubstitutionMatrix* _3merSubMatrix,
                IndexTable * indexTable,
                short kmerThr,
                float prefThr,
                int kmerSize,
                int dbSize,
                int alphabetSize); 

        ~QueryTemplateMatcher();
        // returns result for the sequence
        std::list<hit_t>*  matchQuery (Sequence * seq);

        void printStats();

    private:
        /* generates kmer lists */
        KmerGenerator * kmerGenerator;
        /* contains the sequences for a kmer */
        IndexTable * indexTable;
        /* calculates the score */
        QueryScore * queryScore;
        // k of the k-mer
        int kmerSize;
        // amino acid alphabet size
        int alphabetSize;




};

#endif
