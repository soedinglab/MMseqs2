#include <iostream>
#include "Clustering.h"
#include "Parameters.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif


int clust(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    Clustering* clu = new Clustering(par.db1, par.db1Index, par.db2, par.db2Index,
                                     par.db3, par.db3Index, par.maxIteration,
                                     par.similarityScoreType, par.threads, par.compressed);

    clu->run(par.clusteringMode);

    delete clu;
    return 0;    
}

