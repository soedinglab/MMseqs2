#ifndef ALIGNMENT_H
#define ALIGNMENT_H



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

class Alignment {

    public:

        Alignment (std::string querySeqDB, std::string querySeqDBIndex,
                std::string targetSeqDB, std::string targetSeqDBIndex,
                std::string prefDB, std::string prefDBIndex,
                std::string outDB, std::string outDBIndex,
                std::string matrixFile, double evalThr, double covThr, int maxSeqLen, int seqType);

        ~Alignment();

        void run (int maxAlnNum, int maxRejected);

    private:

        int threads;

        size_t BUFFER_SIZE;

        double covThr;

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
