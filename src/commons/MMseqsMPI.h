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
    static void finalize();
    static bool broadcast(bool value);

    static inline bool isMaster() {
#ifdef HAVE_MPI
        return rank == MASTER;
#else
        return true;
#endif
    };

    static inline void barrier() {
#ifdef HAVE_MPI
        if (active && numProc > 1) {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
    };
};

// Do not finalize a failed multi-rank job: another rank might still be blocked in a
// collective. Abort the communicator so the complete job fails instead of hanging.
#ifdef HAVE_MPI
#define EXIT(exitCode) do {                        \
    int __status = (exitCode);                     \
    if(MMseqsMPI::active) {                        \
        if (__status == 0) {                       \
            MMseqsMPI::finalize();                  \
        } else if (MMseqsMPI::numProc > 1) {       \
            std::cerr.flush();                     \
            std::cout.flush();                     \
            MPI_Abort(MPI_COMM_WORLD, __status);   \
        }                                          \
    }                                              \
    std::cerr.flush();                             \
    std::cout.flush();                             \
    exit(__status);                                \
} while(0)
#endif

#endif
