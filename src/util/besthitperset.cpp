#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"

#ifdef OPENMP
#include <omp.h>
#endif

int besthitperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, true);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    BestHitBySetFilter aggregation(par.db2, par.db3, par.db4, par.simpleBestHit, (unsigned int) par.threads);
    return aggregation.run();
}
