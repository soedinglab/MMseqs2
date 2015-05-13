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

#ifdef HAVE_MPI
#define EXIT(exitCode) do {             \
    if(MMseqsMPI::active == true) {     \
        MPI_Finalize();                 \
        MMseqsMPI::active = false;      \
    }                                   \
    exit(exitCode);                     \
} while(0)
#endif

#endif