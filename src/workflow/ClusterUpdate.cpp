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
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RECOVER_DELETED", par.recoverDeleted ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());

    int maxAccept = par.maxAccept;
    par.maxAccept = 1;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
    par.maxAccept = maxAccept;

    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterworkflow).c_str());

    std::string scriptPath(par.db6);
    if(FileUtil::directoryExists(par.db6.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db6 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db6.c_str()) == false){
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db6 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db6 << "\n";
        }
    }

    scriptPath.append("/update_clustering.sh");
    FileUtil::writeFile(scriptPath, update_clustering_sh, update_clustering_sh_len);

	cmd.execProgram(scriptPath.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
