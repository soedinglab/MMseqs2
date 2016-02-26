#ifndef TIME_TEST_H
#define TIME_TEST_H

#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "DBReader.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "Sequence.h"
#include "NucleotideMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"
#include "Prefiltering.h"

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

        Sequence** seqs;

        std::string scoringMatrixFile;

        std::string targetDBIndex;
        std::string targetDB;
};

#endif
