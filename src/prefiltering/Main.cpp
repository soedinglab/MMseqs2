#include <iostream>
#include <unistd.h>
#include <string>
#include <signal.h>
#include <execinfo.h>
#include "Prefiltering.h"
#include "CommandDeclarations.h"

#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#define MPI_MASTER 0
#endif


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

    std::string usage("\nCalculates k-mer similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n");
    usage.append("USAGE: prefilter <queryDB> <targetDB> <outDB> [opts]\n");

    std::vector<MMseqsParameter> perfPar = {
            Parameters::PARAM_S,
            Parameters::PARAM_K,
            Parameters::PARAM_K_SCORE,
            Parameters::PARAM_ALPH_SIZE,
            Parameters::PARAM_MAX_SEQ_LEN,
            Parameters::PARAM_PROFILE,
            Parameters::PARAM_NUCL,
            Parameters::PARAM_Z_SCORE,
            Parameters::PARAM_SKIP,
            Parameters::PARAM_MAX_SEQS,
            Parameters::PARAM_SPLIT,
            Parameters::PARAM_SEARCH_MODE,
            Parameters::PARAM_NO_COMP_BIAS_CORR,
            Parameters::PARAM_FAST_MODE,
            Parameters::PARAM_SPACED_KMER_MODE,
            Parameters::PARAM_SUB_MAT,
            Parameters::PARAM_THREADS,
            Parameters::PARAM_V};
    Parameters par;
    par.parseParameters(argc, argv, usage, perfPar, 3);



#ifdef HAVE_MPI
    int mpi_error,mpi_rank,mpi_num_procs;
    mpi_error = MPI_Init(&argc, &argv);
    mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);

    if(mpi_rank != MPI_MASTER)
        Debug::setDebugLevel(Debug::NOTHING);

    Debug(Debug::WARNING) << "MPI Init...\n";
    Debug(Debug::WARNING) << "Rank: " << mpi_rank << " Size: " << mpi_num_procs << "\n";

#endif
    mmseqs_cuticle_init();

    struct timeval start, end;
    gettimeofday(&start, NULL);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Debug::setDebugLevel(par.verbosity);

    Debug(Debug::WARNING) << "Initialising data structures...\n";
    Prefiltering* pref = new Prefiltering(par.db1,par.db1Index,
            par.db2,par.db2Index,
            par.db3,par.db3Index,par);

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for init: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n\n";
    gettimeofday(&start, NULL);
#ifdef HAVE_MPI
    pref->run(mpi_rank, mpi_num_procs);
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
