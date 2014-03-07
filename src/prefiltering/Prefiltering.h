#ifndef PREFILTERING_H
#define PREFILTERING_H

#include <iostream>
#include <unistd.h>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/SubstitutionMatrix.h"
#include "../commons/Sequence.h"
#include "../commons/NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

class Prefiltering {

    public:
        
        Prefiltering(std::string queryDB, 
                std::string queryDBIndex,
                std::string targetDB, 
                std::string targetDBIndex,
                std::string outDB, 
                std::string ouDBIndex,
                std::string scoringMatrixFile, 
                float sensitivity, 
                int kmerSize, 
                int alphabetSize, 
                float zscoreThr, 
                size_t maxSeqLen, 
                int seqType, 
                bool aaBiasCorrection, 
                int skip);

        ~Prefiltering();

        void run (size_t maxResListLen);

        static IndexTable* getIndexTable(DBReader* dbr, Sequence* seq, int alphabetSize, int kmerSize, int dbFrom, int dbTo, int skip = 0);

    private:

        size_t BUFFER_SIZE;

        int threads;

        DBReader* qdbr;
        DBReader* tdbr;
        DBWriter* dbw;

        Sequence** seqs;
        std::list<int>** reslens;
        BaseMatrix* subMat;
        ExtendedSubstitutionMatrix* _2merSubMatrix;
        ExtendedSubstitutionMatrix* _3merSubMatrix;
        char** outBuffers;
        QueryTemplateMatcher** matchers;
        IndexTable* indexTable;

        int kmerSize;
        int alphabetSize;
        size_t maxSeqLen;

        BaseMatrix* getSubstitutionMatrix(std::string scoringMatrixFile, float bitFactor);

        /* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
         * K-mer similarity threshold is set to meet a certain DB match probability.
         * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
         */
        std::pair<short,double> setKmerThreshold(DBReader* dbr, double targetKmerMatchProb, double toleratedDeviation);


};

#endif
