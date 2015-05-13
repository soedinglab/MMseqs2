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

void mmseqs_debug_catch_signal(int sig_num)
{
    if(sig_num == SIGILL)
    {
        fprintf(stderr, "Your CPU does not support all the latest features that this version of mmseqs makes use of!\n"
                "Please run on a newer machine.");
        EXIT(sig_num);
    }
    else
    {
        fprintf (stderr, "\n\n-------------------------8<-----------------------\nExiting on error!\n");
        fprintf (stderr, "Signal %d received\n",sig_num);
        perror("ERROR (can be bogus)");

        fprintf(stderr, "Backtrace:");
        void *buffer[30];
        int nptrs = backtrace(buffer, 30);
        backtrace_symbols_fd(buffer, nptrs, 2);
        fprintf (stderr, "------------------------->8-----------------------\n\n"
                "Send the binary program that caused this error and the coredump (ls core.*).\n"
                "Or send the backtrace:"
                "\n$ gdb -ex=bt --batch PROGRAMM_NAME CORE_FILE\n"
                "If there is no core file, enable coredumps in your shell and run again:\n"
                "$ ulimit -c unlimited\n\n");
    }

    EXIT(1);
}

void mmseqs_cuticle_init()
{
    struct sigaction handler;
    handler.sa_handler = mmseqs_debug_catch_signal;
    sigemptyset(&handler.sa_mask);
    handler.sa_flags = 0;

    sigaction(SIGFPE, &handler, NULL);
    sigaction(SIGSEGV, &handler, NULL);
    sigaction(SIGBUS, &handler, NULL);
    sigaction(SIGABRT, &handler, NULL);
}

int prefilter(int argc, const char **argv)
{
    MMseqsMPI::init(argc, argv);

    std::string usage("\nCalculates k-mer similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n");
    usage.append("USAGE: prefilter <queryDB> <targetDB> <outDB> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.prefilter, 3);

    mmseqs_cuticle_init();

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
    pref->run(MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    pref->run();
#endif
    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "\nTime for prefiltering run: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    delete pref;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
