#include <iostream>
#include <string>
#include <sys/time.h>

#include "Prefiltering.h"
#include "Util.h"
#include "Parameters.h"

#include "MMseqsMPI.h"

#ifdef OPENMP
#include <omp.h>
#endif

int prefilter(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, false, MMseqsParameter::COMMAND_PREFILTER);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Debug(Debug::INFO) << "Initialising data structures...\n";

    Prefiltering pref(par.db2, par.db2Index, par);
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";

#ifdef HAVE_MPI
    pref.runMpiSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#else
    pref.runAllSplits(par.db1, par.db1Index, par.db3, par.db3Index);
#endif

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nOverall time for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    EXIT(EXIT_SUCCESS);
}
