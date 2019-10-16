#include "Alignment.h"
#include "Parameters.h"
#include "Debug.h"
#include "Util.h"
#include "MMseqsMPI.h"

#ifdef OPENMP
#include <omp.h>
#endif


int align(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    Alignment aln(par.db1, par.db2,
                  par.db3, par.db3Index,
                  par.db4, par.db4Index, par);

    Debug(Debug::INFO) << "Calculation of alignments\n";

#ifdef HAVE_MPI
    aln.run(MMseqsMPI::rank, MMseqsMPI::numProc, par.maxAccept, par.maxRejected, par.wrappedScoring);
#else
    aln.run(par.maxAccept, par.maxRejected, par.wrappedScoring);
#endif

    return EXIT_SUCCESS;
}


