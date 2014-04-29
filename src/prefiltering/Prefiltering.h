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
#include "../commons/Debug.h"
#include "../commons/Log.h"
#include "../commons/Util.h"

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
                int maxResListLen,
                int alphabetSize, 
                float zscoreThr, 
                size_t maxSeqLen, 
                int querySeqType,
                int targetSeqType,
                bool aaBiasCorrection,
                int split,
                int skip);

        ~Prefiltering();
        void run (size_t dbFrom,size_t dbSize,
                  std::string resultDB, std::string resultDBIndex);
        void run(int mpi_rank, int mpi_num_procs);
        void run ();
        void closeReader();
        void mergeOutput(std::vector<std::pair<std::string, std::string> > filenames);
        void removeDatabaes(std::vector<std::pair<std::string, std::string> > filenames);
        IndexTable * getIndexTable(int split, int splitCount); // needed for index lookup
    
        static IndexTable* generateIndexTable(DBReader* dbr, Sequence* seq,
                                         int alphabetSize, int kmerSize,
                                         size_t dbFrom, size_t dbTo, int skip = 0);
        // get substituion matrix
        static BaseMatrix* getSubstitutionMatrix(std::string scoringMatrixFile, int alphabetSize, float bitFactor);


    private:
        static const size_t BUFFER_SIZE = 1000000;

        int threads;

        DBReader* qdbr;
        DBReader* tdbr;
        DBReader* tidxdbr;


        Sequence** qseq;
        int* notEmpty;

        std::list<int>** reslens;
        BaseMatrix* subMat;
        ExtendedSubstitutionMatrix* _2merSubMatrix;
        ExtendedSubstitutionMatrix* _3merSubMatrix;

        std::string outDB;
        std::string outDBIndex;

        int kmerSize;
        int alphabetSize;
        float zscoreThr;
        size_t maxSeqLen;
        int querySeqType;
        int targetSeqType;
        bool templateDBIsIndex;
    
        bool aaBiasCorrection;
        short kmerThr;
        double kmerMatchProb;
        int split;
        int skip;
        size_t maxResListLen;
        // statistics
        size_t kmersPerPos;
        size_t resSize;
        size_t dbMatches;

        /* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
         * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
         */
        std::pair<short,double> setKmerThreshold(DBReader* qdbr, DBReader* tdbr, double targetKmerMatchProb, double toleratedDeviation);
        std::pair<std::string, std::string> createTmpFileNames(std::string db, std::string dbindex, int numb);
        // write prefiltering to ffindex database
        int writePrefilterOutput(DBWriter * dbWriter, int thread_idx, size_t id, std::pair<hit_t *,size_t> prefResults);
        // init QueryTemplateMatcher
        QueryTemplateMatcher ** createQueryTemplateMatcher ( BaseMatrix* m, IndexTable * indexTable,
                                   unsigned short * seqLens,
                                   short kmerThr, double kmerMatchProb,
                                   int kmerSize, int dbSize,
                                   bool aaBiasCorrection, int maxSeqLen,
                                   float zscoreThr);
    
    
        void printStatistics();

};

#endif
