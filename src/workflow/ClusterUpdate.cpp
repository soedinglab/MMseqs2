#include "Debug.h"
#include "Util.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include <cassert>

#include "update_clustering.sh.h"

extern void setClusterAutomagicParameters(Parameters& par);

void setClusterUpdateDefaults(Parameters *p) {
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}
void setClusterUpdateMustPassAlong(Parameters *p) {
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
}
int clusterupdate(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setClusterUpdateDefaults(&par);
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ_SCALE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_START_SENS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SENS_STEPS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_CLUSTER_REASSIGN.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_NUM_ITERATIONS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    setClusterUpdateMustPassAlong(&par);
    setClusterAutomagicParameters(par);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RECOVER_DELETED", par.recoverDeleted ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("DIFF_PAR", par.createParameterString(par.diff).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("RESULT2REPSEQ_PAR", par.createParameterString(par.result2repseq).c_str());

    cmd.addVariable("CLUST_PAR", par.createParameterString(par.clusterworkflow, true).c_str());

    int maxAccept = par.maxAccept;
    par.maxAccept = 1;
    par.PARAM_MAX_ACCEPT.wasSet = true;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.clusterUpdateSearch, true).c_str());
    par.maxAccept = maxAccept;

    std::string tmpDir = par.db6;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.clusterUpdate));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    std::string program = tmpDir + "/update_clustering.sh";
    FileUtil::writeFile(program, update_clustering_sh, update_clustering_sh_len);
	cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
