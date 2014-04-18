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
                int seqType, 
                bool aaBiasCorrection,
                int splitSize,
                int skip);

        ~Prefiltering();
        void run (size_t dbFrom,size_t dbSize,
                  std::string resultDB, std::string resultDBIndex);
        void run(int mpi_rank, int mpi_num_procs);
        void run ();
        void closeReader();
        void mergeOutput(std::vector<std::pair<std::string, std::string> > filenames);
        void removeDatabaes(std::vector<std::pair<std::string, std::string> > filenames);
        static IndexTable* getIndexTable(DBReader* dbr, Sequence* seq,
                                         int alphabetSize, int kmerSize,
                                         size_t dbFrom, size_t dbTo, int skip = 0);

    private:
        static const size_t BUFFER_SIZE = 1000000;

        int threads;

        DBReader* qdbr;
        DBReader* tdbr;

        Sequence** seqs;
        int* notEmpty;

        std::list<int>** reslens;
        BaseMatrix* subMat;
        ExtendedSubstitutionMatrix* _2merSubMatrix;
        ExtendedSubstitutionMatrix* _3merSubMatrix;
        char** outBuffers;
        QueryTemplateMatcher** matchers;
        IndexTable* indexTable;

        std::string outDB;
        std::string outDBIndex;

        int kmerSize;
        int alphabetSize;
        float zscoreThr;
        size_t maxSeqLen;
        int seqType;
        bool aaBiasCorrection;
        short kmerThr;
        double kmerMatchProb;
        int splitSize;
        int skip;
        int maxResListLen;
        // statistics
        size_t kmersPerPos;
        size_t resSize;
        size_t dbMatches;
        BaseMatrix* getSubstitutionMatrix(std::string scoringMatrixFile, float bitFactor);

        /* Set the k-mer similarity threshold that regulates the length of k-mer lists for each k-mer in the query sequence.
         * As a result, the prefilter always has roughly the same speed for different k-mer and alphabet sizes.
         */
        std::pair<short,double> setKmerThreshold(DBReader* dbr, double targetKmerMatchProb, double toleratedDeviation);
        std::pair<std::string, std::string> createTmpFileNames(std::string db, std::string dbindex, int numb);
        // write prefiltering to ffindex database
        int writePrefilterOutput(DBWriter * dbWriter, int thread_idx, size_t id, std::list<hit_t>* prefResults);
        // if query key exists in target db it will be added to the result
        void addIdIfExistInTargetDb(int thread_idx, std::list<hit_t>* prefResults);

        void printStatistics();

};

#endif
