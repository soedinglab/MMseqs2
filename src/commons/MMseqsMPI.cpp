#include "MMseqsMPI.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

#include <cstdlib>
#include <cstring>

bool MMseqsMPI::active = false;
int MMseqsMPI::rank = -1;
int MMseqsMPI::numProc = -1;

#ifdef HAVE_MPI
void MMseqsMPI::init(int argc, const char **argv) {
    MPI_Init(&argc, const_cast<char ***>(&argv));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    active = true;

    if(!isMaster()) {
        Parameters& par = Parameters::getInstance();
        par.verbosity = Debug::ERROR;
        Debug::setDebugLevel(Debug::ERROR);
    }

    Debug(Debug::INFO) << "MPI Init\n";
    Debug(Debug::INFO) << "Rank: " << rank << " Size: " << numProc << "\n";
}
#else
// Best-effort fail-fast: this binary was built without MPI (-DHAVE_MPI=0/unset), so it
// can only ever act as a single rank/process with no awareness of any other copies of
// itself. Launching it under a multi-task runner (srun -n>1, mpirun -np>1, ...) would
// silently run every task as an independent "rank 0 of 1" that reads the full input and
// writes to the same output paths -- exactly the unsafe scenario multi-node support is
// meant to avoid. There is no fully portable way to detect this from inside a plain
// (non-MPI) process, so this only checks common launcher environment variables set by
// OpenMPI, MPICH-family MPI implementations, and native (non-MPI) Slurm process launch;
// it is a safety net, not a guarantee.
static int getEnvTaskCount(const char *name) {
    const char *value = getenv(name);
    if (value == NULL || *value == '\0') {
        return -1;
    }
    char *end = NULL;
    long parsed = strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed <= 0) {
        return -1;
    }
    return static_cast<int>(parsed);
}

void MMseqsMPI::init(int, const char **) {
    rank = 0;

    const char *taskCountVars[] = {
        "OMPI_COMM_WORLD_SIZE",  // OpenMPI
        "PMI_SIZE",              // MPICH / Intel MPI / Hydra
        "MPI_LOCALNRANKS",       // some MPICH-family launchers
        "SLURM_NTASKS",          // native srun (no MPI library involved)
        "SLURM_NPROCS"           // legacy alias for SLURM_NTASKS
    };
    for (size_t i = 0; i < sizeof(taskCountVars) / sizeof(taskCountVars[0]); ++i) {
        const int taskCount = getEnvTaskCount(taskCountVars[i]);
        if (taskCount > 1) {
            Debug(Debug::ERROR) << "This mmseqs binary was built without MPI support, but appears to have "
                                << "been launched with " << taskCountVars[i] << "=" << taskCount
                                << " (multiple tasks/ranks). Every task would independently process the "
                                << "full input and write to the same output paths, corrupting the result.\n"
                                << "Rebuild with -DHAVE_MPI=1 to use --mpi-runner/RUNNER with more than one "
                                << "task, or launch this binary with a single task.\n";
            EXIT(EXIT_FAILURE);
        }
    }
}
#endif

