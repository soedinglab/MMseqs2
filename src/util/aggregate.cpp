#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"

#ifdef OPENMP
#include <omp.h>
#endif

int aggregate(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true, true);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Aggregation *aggregation = NULL;
    if (par.aggregationMode == 0 || par.aggregationMode == 1) {
        Debug(Debug::INFO) << "Aggregation by Best Hit" << "\n";
        aggregation = new BestHitAggregator(par.db2, par.db3, par.db4, par.aggregationMode == 1, (unsigned int) par.threads);
    } else if (par.aggregationMode == 2) {
        Debug(Debug::INFO) << "Aggregation of p-values" << "\n";
        aggregation = new PvalueAggregator(par.db1, par.db2, par.db3, par.db4, par.alpha, (unsigned int) par.threads);
    } else if (par.aggregationMode == 3) {
        Debug(Debug::INFO) << "Calculation of median intergene length" << "\n";
        aggregation = new HitDistanceAggregator(par.db1, par.db2, par.db3, par.db4, par.shortOutput, par.alpha, (unsigned int) par.threads);
    } else {
        Debug(Debug::ERROR) << "Unknown aggregation mode " << par.aggregationMode << "\n";
        return EXIT_FAILURE;
    }

    return aggregation->run();
}
