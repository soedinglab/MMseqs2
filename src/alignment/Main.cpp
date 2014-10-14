#include "Alignment.h"
#include "Parameters.h"
#include "CommandDeclarations.h"
#include <iostream>
#include <string>
#include <sys/time.h>
#ifdef OPENMP
#include <omp.h>
#endif

bool compareHits (Matcher::result_t first, Matcher::result_t second){
    if (first.score > second.score)
        return true;
    return false;
}


int alignment(int argc, const char *argv[])
{

    std::string usage("\nCalculates Smith-Waterman alignment scores between all sequences in the query database and the sequences of the target database which passed the prefiltering.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: alignment <queryDB> <targetDB> <prefResultsDB> <outDB> [opts]\n");
    std::vector<MMseqsParameter> perfPar = {
            Parameters::PARAM_E,
            Parameters::PARAM_C,
            Parameters::PARAM_MAX_SEQ_LEN,
            Parameters::PARAM_MAX_SEQS,
            Parameters::PARAM_MAX_REJECTED,
            Parameters::PARAM_NUCL,
            Parameters::PARAM_SUB_MAT,
            Parameters::PARAM_THREADS,
            Parameters::PARAM_V};
    Parameters par;
    par.parseParameters(argc, (char**)argv, usage, perfPar, 4);


    Debug::setDebugLevel(Debug::INFO);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Debug::setDebugLevel(par.verbosity);

    Debug(Debug::WARNING) << "Init data structures...\n";
    Alignment* aln = new Alignment(par.db1,           par.db1Index,
                                   par.db2,           par.db2Index,
                                   par.db3,           par.db3Index,
                                   par.db4,           par.db4Index,
                                   par);

    Debug(Debug::WARNING) << "Calculation of Smith-Waterman alignments.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    aln->run(par.maxResListLen, par.maxRejected);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for alignments calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
 
    delete aln;

    return 0;
}


