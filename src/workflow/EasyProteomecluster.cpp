#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

#include "easyproteomecluster.sh.h"

void setEasyproteomeclusterDefaults(Parameters *p) {
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->weightClusterCount = 1.0;
    p->proteomeSimThr = 0.9;
    p->proteomeRelativeSimThr = 0.0;
    p->proteomeCascadedClustering = 0;
    p->includeAlignFiles = true;
    // p->proteomeScoreMode = 0;
}

void setEasyproteomeclusterMustPassAlong(Parameters *p){
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
}

int easyproteomecluster(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEasyproteomeclusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyproteomeclusterMustPassAlong(&par);
    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("PROTEOME_HIDDEN_REPORT", par.proteomeHiddenReport ? "TRUE" : NULL);
    cmd.addVariable("WRITE_ALIGN_PROTEOME", par.includeAlignFiles ? "TRUE" : NULL);
    cmd.addVariable("CASCADED_PROTEOME_CLUSTERING", par.proteomeCascadedClustering ? "TRUE" : NULL);
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clusterworkflow, true).c_str());
    if (par.clusterModule == 0) {
        cmd.addVariable("CLUSTER_MODULE", "linclust");
    } else if (par.clusterModule == 1) {
        cmd.addVariable("CLUSTER_MODULE", "cluster");
    } else {
        assert(false);
    }
    cmd.addVariable("PROTEOMECLUSTER_PAR", par.createParameterString(par.proteomecluster).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyproteomecluster.sh";
    FileUtil::writeFile(program, easyproteomecluster_sh, easyproteomecluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    
    // Should never get here
    assert(false);
    return 0;
}