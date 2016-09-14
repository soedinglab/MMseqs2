#include <iostream>

#include <sys/time.h>
#include <cassert>
#include <FileUtil.h>
#include <update_clustering.sh.h>

#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"


int clusterupdate (int argc, const char * argv[]){
    std::string usage("\nUpdates the existing clustering of the previous database version with new sequences from the current version of the same database.\n");
    usage.append("Written by Clovis Galiez (clovis.galiez@mpibpc.mpg.de)\n");
    usage.append("USAGE: clusterupdate <oldDB> <newDB> <oldDB_clustering> <newDB_clustering> <tmpDir> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.clusterUpdate, 5);

    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
	
	cmd.addVariable("RUNNER", par.runner.c_str());
	cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch).c_str());
	cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterUpdateClust).c_str());

    FileUtil::writeFile(par.db5 + "/update_clustering.sh", update_clustering_sh, update_clustering_sh_len);
    std::string program(par.db5 + "/update_clustering.sh");
	cmd.execProgram(program.c_str(), 5, argv);

    // Should never get here
    assert(false);
    return 0;
}
