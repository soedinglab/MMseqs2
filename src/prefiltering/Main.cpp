#include <iostream>
#include <unistd.h>
#include <string>
#include <signal.h>
#include <execinfo.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "Prefiltering.h"
#include "CommandDeclarations.h"
#include "Parameters.h"

#include "MMseqsMPI.h"


int prefilter(int argc, const char **argv)
{
    MMseqsMPI::init(argc, argv);

    std::string usage("\nCalculates k-mer similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n");
    usage.append("USAGE: prefilter <queryDB> <targetDB> <outDB> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.prefilter, 3);

    struct timeval start, end;
    gettimeofday(&start, NULL);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#ifdef HAVE_MPI
    if(MMseqsMPI::isMaster())
    {
        Debug::setDebugLevel(par.verbosity);
    }
#else
    Debug::setDebugLevel(par.verbosity);
#endif
    Debug(Debug::WARNING) << "Initialising data structures...\n";
    Prefiltering* pref = new Prefiltering(par.db1,par.db1Index,
            par.db2,par.db2Index,
            par.db3,par.db3Index,par);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    gettimeofday(&start, NULL);

#ifdef HAVE_MPI
    pref->run(MMseqsMPI::rank, MMseqsMPI::numProc, par.splitMode);
#else
    pref->run();
#endif
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nOverall time for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    delete pref;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
