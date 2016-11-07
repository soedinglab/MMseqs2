#include "Alignment.h"
#include "Parameters.h"
#include "Debug.h"
#include <string>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "MMseqsMPI.h"


int align(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true, false, MMseqsParameter::COMMAND_ALIGN);

    MMseqsMPI::init(argc, argv);

    Debug::setDebugLevel(Debug::INFO);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Debug(Debug::WARNING) << "Init data structures...\n";
    Alignment* aln = new Alignment(par.db1,           par.db1Index,
                                   par.db2,           par.db2Index,
                                   par.db3,           par.db3Index,
                                   par.db4,           par.db4Index,
                                   par);

    Debug(Debug::WARNING) << "Calculation of Smith-Waterman alignments.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

#ifdef HAVE_MPI
    aln->run(MMseqsMPI::rank, MMseqsMPI::numProc, par.maxResListLen, par.maxRejected);
#else
    aln->run(par.maxResListLen, par.maxRejected);
#endif

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for alignments calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
 
    delete aln;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}


