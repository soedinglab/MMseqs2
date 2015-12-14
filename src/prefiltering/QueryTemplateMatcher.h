#ifndef QUERY_TEMPLATE_MATCHER_H
#define QUERY_TEMPLATE_MATCHER_H

#include <list>
#include <iostream>
#include <cstring>
#include <stddef.h>


#include "BaseMatrix.h"
#include "Sequence.h"
#include "ExtendedSubstitutionMatrix.h"
#include "QueryScore.h"
#include "IndexTable.h"
#include "KmerGenerator.h"
#include "Indexer.h"

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

class QueryTemplateMatcher {
    public:
        QueryTemplateMatcher(BaseMatrix *m, IndexTable *indexTable, unsigned int *seqLens, short kmerThr, double kmerMatchProb, int kmerSize,
                size_t dbSize, bool aaBiasCorrection, int maxSeqLen) {
            this->m = m;
            this->indexTable = indexTable;
            this->kmerSize = kmerSize;
            this->kmerThr = kmerThr;
            this->kmerGenerator = new KmerGenerator(kmerSize, m->alphabetSize, kmerThr);
            this->aaBiasCorrection = aaBiasCorrection;
            this->stats = new statistics_t();
        }

        virtual ~QueryTemplateMatcher () {
            delete stats;
            delete kmerGenerator;
        }

        const statistics_t * getStatistics(){
            return stats;
        }

        // set substituion matrix for KmerGenerator
        void setProfileMatrix(ScoreMatrix **matrix){
            this->kmerGenerator->setDivideStrategy(matrix );
        }
    
        // set substitution matrix
        void setSubstitutionMatrix(ScoreMatrix * three, ScoreMatrix * two) {
            this->kmerGenerator->setDivideStrategy(three, two );
        }

        // returns result for the sequence
        // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
        virtual std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId) = 0;
    
    protected:
        // keeps stats for run
        statistics_t * stats;

        // scoring matrix for local amino acid bias correction
        BaseMatrix * m;
        /* generates kmer lists */
        KmerGenerator * kmerGenerator;
        /* contains the sequences for a kmer */
        IndexTable * indexTable;
        // k of the k-mer
        int kmerSize;
        // local amino acid bias correction
        bool aaBiasCorrection;

        // kmer threshold for kmer generator
        short kmerThr;


};

#endif
