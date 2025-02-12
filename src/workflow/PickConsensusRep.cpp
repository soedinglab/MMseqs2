#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include <cassert>
#include <string>
// Include the embedded shell script.
#include "pickconsensusrep.sh.h"

// Minimal workflow function that runs the pickcenterrep workflow.
int pickconsensusrep(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    CommandCaller cmd;
    par.allowDeletion = 1;
    par.PARAM_ALLOW_DELETION.wasSet = true;
    cmd.addVariable("RESULT2MSA_PAR", par.createParameterString(par.result2msa, true).c_str());
    par.matchMode = 1;
    par.PARAM_MATCH_MODE.wasSet = true;
    cmd.addVariable("MSA2PROFILE_PAR", par.createParameterString(par.msa2profile, true).c_str());
    cmd.addVariable("RENAMEDBKEYS_PAR", par.createParameterString(par.renamedbkeys).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    // The temporary directory is provided as the 4th argument.
    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.mapworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    // Write out the embedded shell script to a file in the temporary directory.
    std::string program = tmpDir + "/pickconsensusrep.sh";
    FileUtil::writeFile(program, pickconsensusrep_sh, pickconsensusrep_sh_len);

    // Execute the shell script.
    cmd.execProgram(program.c_str(), par.filenames);

    // The shell script should not return; if it does, abort.
    assert(false);
    return 0;
}

