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
    static inline bool isMaster() { return rank == MASTER; };
};

// if we are in an error case, do not call MPI_Finalize, it might still be in a Barrier
#ifdef HAVE_MPI
#define EXIT(exitCode) do {             \
    if ((exitCode) != 0)               \
        exit(exitCode);                 \
    if(MMseqsMPI::active == true) {     \
        MPI_Finalize();                 \
        MMseqsMPI::active = false;      \
    }                                   \
    exit(exitCode);                     \
} while(0)
#endif

#endif
