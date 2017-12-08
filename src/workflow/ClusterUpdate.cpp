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

    int maxAccept = par.maxAccept;
    par.maxAccept = 1;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
    par.maxAccept = maxAccept;

    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusteringWorkflow).c_str());

    std::string scriptPath(par.db6);
    if(FileUtil::directoryExists(par.db6.c_str())==false){
        Debug(Debug::WARNING) << "Tmp " << par.db6 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db6.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db6 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::WARNING) << "Created dir " << par.db6 << "\n";
        }
    }

    scriptPath.append("/update_clustering.sh");
    FileUtil::writeFile(scriptPath, update_clustering_sh, update_clustering_sh_len);

	cmd.execProgram(scriptPath.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
