#include <cassert>

#include "Debug.h"
#include "Util.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include "update_clustering.sh.h"

int clusterupdate(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 6);

    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }

    if(par.recoverDeleted) {
        cmd.addVariable("RECOVER_DELETED", "TRUE");
    }
	
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusteringWorkflow).c_str());

    std::string scriptPath(par.db6);
    if(!FileUtil::directoryExists(par.db6.c_str())) {
        Debug(Debug::ERROR) << "Tmp directory " << par.db6 << " not found!\n";
        return EXIT_FAILURE;
    }

    scriptPath.append("/update_clustering.sh");
    FileUtil::writeFile(scriptPath, update_clustering_sh, update_clustering_sh_len);
    std::string program(par.db6 + "/update_clustering.sh");
	cmd.execProgram(program.c_str(), 6, argv);

    // Should never get here
    assert(false);
    return 0;
}
