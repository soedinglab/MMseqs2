#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#ifdef OPENMP
#include <omp.h>
#endif


#include <string>
#include <list>
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>

#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"
#include "Log.h"

#include "Matcher.h"
#include "Parameters.h"


class Alignment {

    public:

        Alignment (std::string querySeqDB, std::string querySeqDBIndex,
                std::string targetSeqDB, std::string targetSeqDBIndex,
                std::string prefDB, std::string prefDBIndex,
                std::string outDB, std::string outDBIndex,
                Parameters par);

        ~Alignment();

        void run (const unsigned int maxAlnNum, const unsigned int maxRejected);
        void run (const unsigned int maxAlnNum, const unsigned int maxRejected, const unsigned int mpiRank, const unsigned int mpiNumProc);

    private:

        // keeps state of alignment mode (SCORE_ONLY, SCORE_COV or SCORE_COV_SEQID)
        unsigned mode;
        int threads;

        size_t BUFFER_SIZE;

        // sequence identity threshold
        double seqIdThr;

        // query sequence coverage threshold
        double covThr;

        // e value threshold
        double evalThr;

        BaseMatrix* m;

        // 2 Sequence objects for each thread: one for the query, one for the DB sequence
        Sequence** qSeqs;

        Sequence** dbSeqs;

        Matcher** matchers;

        DBReader* qseqdbr;

        DBReader* tseqdbr;

        DBReader* prefdbr;

        DBWriter* dbw;

        // buffers for the database keys (needed during the processing of the prefilterings lists)
        char** dbKeys;

        // output buffers
        char** outBuffers;

};

#endif
