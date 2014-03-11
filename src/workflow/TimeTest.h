#ifndef TIME_TEST_H
#define TIME_TEST_H

#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/SubstitutionMatrix.h"
#include "../commons/Sequence.h"
#include "../commons/NucleotideMatrix.h"
#include "../prefiltering/ExtendedSubstitutionMatrix.h"
#include "../prefiltering/ReducedMatrix.h"
#include "../prefiltering/KmerGenerator.h"
#include "../prefiltering/QueryTemplateMatcher.h"
#include "../prefiltering/Prefiltering.h"

#ifdef OPENMP
#include <omp.h>
#endif

class TimeTest {

    public:

        TimeTest(std::string queryDB,
                std::string queryDBIndex,
                std::string targetDB,
                std::string targetDBIndex,
                std::string scoringMatrixFile,
                size_t maxSeqLen,
                std::string logFile);

        ~TimeTest();

        void runTimeTest();

    private:

        size_t BUFFER_SIZE;

        std::string logFile;

        int threads;

        size_t maxSeqLen;

        DBReader* qdbr;
        DBReader* tdbr;
        DBWriter* dbw;

        Sequence** seqs;

        std::string scoringMatrixFile;
};

#endif
