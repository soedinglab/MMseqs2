#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"
#include "Util.h"


int aggregate(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2, true, true);
    Aggregation *aggregFunction = nullptr;

    if(par.mode=="bestHit"){
        Debug(Debug::INFO) << "Aggregation by Best Hit" << "\n";
        aggregFunction = new BestHitAggregator(par.db1, par.db2, par.db3, par.setColumn, (unsigned int)par.threads) ;
    }
    else if (par.mode=="pval"){
        Debug(Debug::INFO) << "Aggregation of p-values" << "\n";
        aggregFunction = new PvalAggregator(par.db1, par.db2, (unsigned int)par.threads, par.db3, par.setColumn) ;
    }
    else if (par.mode=="clustering-index") {
        Debug(Debug::INFO) << "Calculation of clusterig index" << "\n";
        aggregFunction = new ClusteringAggregator(par.db1, par.db2, (unsigned int)par.threads) ;
    }
    else {
        Debug(Debug::ERROR) << "Unknown aggregation mode " << par.mode << "\n";
        EXIT(EXIT_FAILURE);
    }

    aggregFunction->runAggregate();

    return 0 ;
}
