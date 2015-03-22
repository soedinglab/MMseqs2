
#include <iostream>
#include <time.h>
#include "Clustering.h"
#include "Parameters.h"
#include "CommandDeclarations.h"

int cluster(int argc, const char ** argv)
{

    std::string usage("\nCalculates a clustering of a sequence database based on Smith Waterman alignment scores of the sequence pairs.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: mmseqs_clu <sequenceDB> <alnResultsDB> <outDB> [opts]\n");
    std::vector<MMseqsParameter> perfPar = {
        Parameters::PARAM_G,
        Parameters::PARAM_A,
        Parameters::PARAM_MAX_SEQS,
        Parameters::PARAM_V,
        Parameters::PARAM_MAXITERATIONS,
        Parameters::PARAM_CONVERGENCEITERATIONS,
        Parameters::PARAM_DAMPING,
        Parameters::PARAM_SIMILARITYSCORE};
    Parameters par;
    par.parseParameters(argc, argv, usage, perfPar, 3);

    Debug::setDebugLevel(par.verbosity);

    Clustering* clu = new Clustering(par.db1, par.db1Index, par.db2, par.db2Index,
                                     par.db3, par.db3Index, par.validateClustering, par.maxResListLen,par.maxIteration,par.convergenceIterations,par.dampingFactor,par.similarityScoreType);

    clu->run(par.clusteringMode);

    delete clu;
    return 0;    
}

