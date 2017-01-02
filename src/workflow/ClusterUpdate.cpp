#include <cassert>
#include <Debug.h>

#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include "update_clustering.sh.h"

int clusterupdate(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 5);

    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
	
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterUpdateClust).c_str());

    std::string scriptPath(par.db5);
    if(!FileUtil::directoryExists(par.db5.c_str())) {
        Debug(Debug::ERROR) << "Tmp directory " << par.db5 << " not found!\n";
        return EXIT_FAILURE;
    }

    scriptPath.append("/update_clustering.sh");
    if(!FileUtil::fileExists(scriptPath.c_str())) {
        FileUtil::writeFile(scriptPath, update_clustering_sh, update_clustering_sh_len);
    } else {
        Debug(Debug::INFO) << "Update script already exists. Skipped overwriting.\n";
    }
    std::string program(par.db5 + "/update_clustering.sh");
	cmd.execProgram(program.c_str(), 5, argv);

    // Should never get here
    assert(false);
    return 0;
}
