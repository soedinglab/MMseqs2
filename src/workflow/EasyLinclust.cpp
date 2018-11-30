#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
namespace linclust {
#include "easycluster.sh.h"
}

void setEasyLinclusterDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.001;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}

int easylinclust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEasyLinclusterDefaults(&par);
    par.overrideParameterDescription((Command &)command, par.PARAM_ADD_BACKTRACE.uniqid, NULL, NULL, par.PARAM_ADD_BACKTRACE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_ALT_ALIGNMENT.uniqid, NULL, NULL, par.PARAM_ALT_ALIGNMENT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL, par.PARAM_RESCORE_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL, par.PARAM_MAX_REJECTED.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_ACCEPT.uniqid, NULL, NULL, par.PARAM_MAX_ACCEPT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_S.uniqid, "Sensitivity will be automatically determined but can be adjusted", NULL, par.PARAM_S.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_INCLUDE_ONLY_EXTENDABLE.uniqid, NULL, NULL, par.PARAM_INCLUDE_ONLY_EXTENDABLE.category | MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.overrideParameterDescription((Command &)command, par.createdb[i]->uniqid, NULL, NULL, par.createdb[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, 3);

    if (FileUtil::directoryExists(par.db3.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db3.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create tmp folder " << par.db3 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db3 << "\n";
        }
    }

    size_t hash = par.hashParameter(par.filenames, *command.params);
    std::string tmpDir = par.db3 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    std::string createdbParam = par.createParameterString(par.createdb);
    cmd.addVariable("CREATEDB_PAR", createdbParam.c_str());
    std::string clusterParam = par.createParameterString(par.linclustworkflow, true);
    cmd.addVariable("CLUSTER_PAR", clusterParam.c_str());
    cmd.addVariable("CLUSTER_MODULE", "linclust");
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/easycluster.sh", linclust::easycluster_sh, linclust::easycluster_sh_len);
    std::string program(tmpDir + "/easycluster.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
