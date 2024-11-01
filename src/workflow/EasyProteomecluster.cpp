#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

#include "easyproteomecluster.sh.h"

void setEasyproteomeclusterDefaults(Parameters *p) {
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->proteomeSimThr = 0.9;
    p->proteomeRelativeSimThr = 0.0;
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

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clusterworkflow, true).c_str());
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
