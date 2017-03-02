#include <iostream>

#include <sys/time.h>
#include <cassert>
#include <FileUtil.h>
#include <update_clustering.sh.h>

#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"


int clusterupdate(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 5);

    bool targetCov = false;
    bool cov = false;
    for (size_t i = 0; i < par.clusterUpdate.size(); i++) {
        if (par.clusterUpdate[i].uniqid == Parameters::PARAM_TARGET_COV_ID && par.clusterUpdate[i].wasSet) {
            targetCov = true;
        }
        if (par.clusterUpdate[i].uniqid == Parameters::PARAM_C_ID && par.clusterUpdate[i].wasSet) {
            cov = true;
        }
    }
    if (cov && targetCov) {
        Debug(Debug::ERROR) << "The paramter -c can not be combined with --target-cov.\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<MMseqsParameter> clustParams;
    for (std::vector<MMseqsParameter>::iterator it = par.clusterUpdateClust.begin(); it != par.clusterUpdateClust.end(); ++it) {
        if ((cov && (*it).uniqid == Parameters::PARAM_TARGET_COV_ID)
                || (targetCov && (*it).uniqid == Parameters::PARAM_C_ID)) {
            continue;
        }
        clustParams.push_back(*it);
    }


    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
	
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterUpdateClust).c_str());

    FileUtil::writeFile(par.db5 + "/update_clustering.sh", update_clustering_sh, update_clustering_sh_len);
    std::string program(par.db5 + "/update_clustering.sh");
	cmd.execProgram(program.c_str(), 5, argv);

    // Should never get here
    assert(false);
    return 0;
}
