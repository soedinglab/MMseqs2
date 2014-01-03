#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#ifdef OPENMP
#include <omp.h>
#endif

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <string>
#include <list>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/NucleotideMatrix.h"
#include "../commons/SubstitutionMatrix.h"

#include "Matcher.h"

class Alignment {

    public:

        Alignment (std::string seqDB, std::string prefDB, std::string outDB, std::string matrixFile, double evalThr, double covThr, int maxSeqLen, int seqType);

        ~Alignment();

        void run (int maxAlnNum);

    private:

        int threads;

        size_t BUFFER_SIZE;

        double covThr;

        double evalThr;

        // 2 Sequence objects for each thread: one for the query, one for the DB sequence
        Sequence** qSeqs;
        Sequence** dbSeqs;

        Matcher** matchers;

        DBReader* seqdbr;

        DBReader* prefdbr;

        DBWriter* dbw;

        // buffers for the database keys (needed during the processing of the prefilterings lists)
        char** dbKeys;

        // output buffers
        char** outBuffers;

};

#endif
