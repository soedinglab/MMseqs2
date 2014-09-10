#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>

#include "WorkflowFunctions.h"
#include "Parameters.h"
extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

void runSearch(Parameters par, std::string queryDB, std::string targetDB, std::string outDB,std::string tmpDir)
{
    std::string queryDBIndex = queryDB + ".index";
    std::string targetDBIndex = targetDB + ".index";
    if (par.zscoreThr == 0.0)
        par.zscoreThr = getZscoreForSensitivity(par.sensitivity);
    
    std::list<std::string>* tmpFiles = new std::list<std::string>();
    
    std::string alnDB = runStep(queryDB, queryDBIndex, targetDB, targetDBIndex, tmpDir,
                                par, 1, 0, true, tmpFiles);
    
    std::string alnDBIndex = alnDB + ".index";
    std::string outDBIndex = outDB + ".index";
    
    // copy the clustering databases to the right location
    copy(alnDBIndex, outDBIndex);
    copy(alnDB, outDB);
    deleteTmpFiles(tmpFiles);
    delete tmpFiles;
}

int search (int argc, const char * argv[]){
    
    std::string usage("\nCompares all sequences in the query database with all sequences in the target database.\n");
    usage.append("Written by Martin Steinegger (martin.steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: search <queryDB> <targetDB> <outDB> <tmpDir> [opts]\n");
    std::vector<MMseqsParameter> perfPar = {
        Parameters::PARAM_S,
        Parameters::PARAM_Z_SCORE,
        Parameters::PARAM_SEARCH_MODE,
        Parameters::PARAM_PROFILE,
        Parameters::PARAM_NUCL,
        Parameters::PARAM_NO_SPACED_KMER,
        Parameters::PARAM_NO_COMP_BIAS_CORR,
        Parameters::PARAM_MAX_SEQS,
        Parameters::PARAM_MAX_SEQ_LEN,
        Parameters::PARAM_V};
    Parameters par;
    par.parseParameters(argc, (char**)argv, usage, perfPar, 4);
    
    Debug::setDebugLevel(par.verbosity);
    runSearch(par, par.db1, par.db2, par.db3, par.db4);
    
    return 0;
}
