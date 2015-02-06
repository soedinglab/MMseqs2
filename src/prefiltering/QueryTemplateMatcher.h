#ifndef QUERY_TEMPLATE_MATCHER_H
#define QUERY_TEMPLATE_MATCHER_H

#include <list>
#include <iostream>
#include <cstring>


#include "BaseMatrix.h"
#include "Sequence.h"
#include "ExtendedSubstitutionMatrix.h"
#include "QueryScore.h"
#include "QueryScoreGlobal.h"
#include "IndexTable.h"
#include "KmerGenerator.h"
#include "Indexer.h"

struct statistics_t{
    float kmersPerPos;
    size_t dbMatches;
    size_t doubleMatches;
    statistics_t() : kmersPerPos(0.0) , dbMatches(0) , doubleMatches(0) {};
};

class QueryTemplateMatcher {
    public:
        QueryTemplateMatcher(BaseMatrix *m,
                IndexTable *indexTable,
                unsigned int *seqLens,
                short kmerThr,
                double kmerMatchProb,
                int kmerSize,
                int dbSize,
                bool aaBiasCorrection,
                int maxSeqLen,
                float zscoreThr){
            this->m = m;
            this->indexTable = indexTable;
            this->kmerSize = kmerSize;
            this->kmerThr = kmerThr;
            this->kmerGenerator = new KmerGenerator(kmerSize, m->alphabetSize, kmerThr);
            this->aaBiasCorrection = aaBiasCorrection;
            this->stats = new statistics_t();
            this->deltaS = new float[maxSeqLen];
            memset(this->deltaS, 0, maxSeqLen * sizeof(float));
        }

        virtual ~QueryTemplateMatcher () {
            delete stats;
            delete[] deltaS;
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

        // calculate local amino acid bias correction score for each position in the sequence
        void calcLocalAaBiasCorrection(Sequence* seq){
            int windowSize = 40;
            if (seq->L < windowSize + 1)
                memset(this->deltaS, 0, seq->L * sizeof(float));
            else{
                float deltaS_i;
                int minPos;
                int maxPos;
                int _2d;
                // calculate local amino acid bias
                for (int i = 0; i < seq->L; i++){
                    deltaS_i = 0.0;
                    minPos = std::max(0, (i - windowSize/2));
                    maxPos = std::min(seq->L, (i + windowSize/2));
                    _2d = maxPos - minPos;
                    
                    // negative score for the amino acids in the neighborhood of i
                    for (int j = minPos; j < maxPos; j++){
                        if (j != i){
                            deltaS_i += m->subMatrix[seq->int_sequence[i]][seq->int_sequence[j]];
                        }
                    }
                    deltaS_i /= -1.0 * _2d;
                    // positive score for the background score distribution for i
                    for (int a = 0; a < m->alphabetSize; a++)
                        deltaS_i += m->pBack[a] * m->subMatrix[seq->int_sequence[i]][a];
                    
                    deltaS[i] = deltaS_i;
                }
            }
        }
        
        
        // returns result for the sequence
        // identityId is the id of the identitical sequence in the target database if there is any, UINT_MAX otherwise
        virtual std::pair<hit_t *, size_t>  matchQuery (Sequence * seq, unsigned int identityId) = 0;
    
    protected:
        // keeps stats for run
        statistics_t * stats;
        // match sequence against the IndexTable
        virtual void match(Sequence* seq) = 0;
        // scoring matrix for local amino acid bias correction
        BaseMatrix * m;
        /* generates kmer lists */
        KmerGenerator * kmerGenerator;
        /* contains the sequences for a kmer */
        IndexTable * indexTable;
        /* calculates the score */
        QueryScore * queryScore;
        // k of the k-mer
        int kmerSize;
        // local amino acid bias correction
        bool aaBiasCorrection;
        // local score correction values
        float* deltaS;
        // kmer threshold for kmer generator
        unsigned short kmerThr;


};

#endif
