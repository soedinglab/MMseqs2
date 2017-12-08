#include "Alignment.h"
#include "Parameters.h"
#include "Debug.h"
#include "Util.h"
#include "MMseqsMPI.h"

#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif


int align(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true, 0, MMseqsParameter::COMMAND_ALIGN);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Debug(Debug::INFO) << "Init data structures...\n";
    Alignment aln(par.db1, par.db1Index, par.db2, par.db2Index,
                  par.db3, par.db3Index, par.db4, par.db4Index, par);

    Debug(Debug::INFO) << "Calculation of Smith-Waterman alignments.\n";

#ifdef HAVE_MPI
    aln.run(MMseqsMPI::rank, MMseqsMPI::numProc, par.maxAccept, par.maxRejected);
#else
    aln.run(par.maxAccept, par.maxRejected);
#endif

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for alignments calculation: " << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
 
    EXIT(EXIT_SUCCESS);
}


