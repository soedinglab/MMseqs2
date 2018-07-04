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
        aggregFunction = new BestHitAggregator(par.db1, par.db2, par.db3, par.setColumn, (unsigned int)par.threads, 1, par.simpleBestHitMode) ; // 3 for pval ; 1 for score
    }
    else if (par.mode=="pval"){
        Debug(Debug::INFO) << "Aggregation of p-values" << "\n";
        aggregFunction = new PvalAggregator(par.db1, par.db2, (unsigned int)par.threads, par.db3, par.db4, par.setColumn, par.alpha) ;
    }
    else if (par.mode=="clustering-index") {
        Debug(Debug::INFO) << "Calculation of median intergene length" << "\n";
        aggregFunction = new GetHitDistance(par.db1, par.db2, par.db3, par.db4, par.db5, (unsigned int)par.threads,par.alpha, par.setColumn) ;
    }
    else {
        Debug(Debug::ERROR) << "Unknown aggregation mode " << par.mode << "\n";
        EXIT(EXIT_FAILURE);
    }

    aggregFunction->runAggregate();

    return 0 ;
}
