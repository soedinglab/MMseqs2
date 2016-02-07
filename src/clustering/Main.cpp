
#include <iostream>
#include "Clustering.h"
#include "Parameters.h"
#include "Debug.h"

int cluster(int argc, const char ** argv)
{

    std::string usage("\nCalculates a clustering of a sequence database based on Smith Waterman alignment scores of the sequence pairs.\n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: cluster <sequenceDB> <alnResultsDB> <outDB> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.clustering, 3);

    Debug::setDebugLevel(par.verbosity);

    Clustering* clu = new Clustering(par.db1, par.db1Index, par.db2, par.db2Index,
                                     par.db3, par.db3Index, par.validateClustering, par.maxIteration, par.similarityScoreType, par.threads);

    clu->run(par.clusteringMode);

    delete clu;
    return 0;    
}

