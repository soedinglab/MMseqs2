#ifndef QUERY_TEMPLATE_MATCHER_H
#define QUERY_TEMPLATE_MATCHER_H

#include <list>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "SubstitutionMatrix.h"
#include "ReduceMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Sequence.h"
#include "QueryScore.h"
#include "IndexTable.h"
#include "KmerGenerator.h"

class QueryTemplateMatcher {
public:
    
    QueryTemplateMatcher (short kmerThreshold,
                          short queryScoreThreshold, 
                          size_t kmer_size,
                          SubstitutionMatrix * substitutionMatrix,
                          IndexTable * indexTable );
    ~QueryTemplateMatcher();
    // returns result for the sequence
    std::list<hit_t>*  matchQuery (Sequence * seq);
    
private:
    /* calculates the score */
    QueryScore * queryScore;
    /* generates kmer lists */
    KmerGenerator * kmerGenerator;
    /* contains the sequences for a kmer */
    IndexTable * indexTable;
    /* size of kmer  */
    size_t kmer_size;

    
    
};

#endif
