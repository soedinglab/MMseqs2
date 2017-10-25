#ifndef MMSEQS_MPI_H
#define MMSEQS_MPI_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

class MMseqsMPI {
public:
    static const int MASTER = 0;

    static bool active;
    static int rank;
    static int numProc;

    static void init(int argc, const char **argv);
    static inline bool isMaster() {
#ifdef HAVE_MPI
        return rank == MASTER;
#else
        return true;
#endif
    };
};

// if we are in an error case, do not call MPI_Finalize, it might still be in a Barrier
#ifdef HAVE_MPI
#define EXIT(exitCode) do {                  \
    int __status = (exitCode);               \
    if(MMseqsMPI::active && __status == 0) { \
        MPI_Finalize();                      \
        MMseqsMPI::active = false;           \
    }                                        \
    std::cerr.flush();                       \
    std::cout.flush();                       \
    exit(__status);                          \
} while(0)
#endif

#endif
