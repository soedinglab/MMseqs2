#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include <cassert>
#include <string>
// Include the embedded shell script.
#include "pickconsensusrepfast.sh.h"

// Fast profile-guided representative selection: reuses the clustering alignments
// (${clusterDB}_aln, produced with --include-align-files) to pick the best observed
// member per cluster and rewrite the cluster DB, without profile-vs-member realignment.
int pickconsensusrepfast(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    CommandCaller cmd;
    cmd.addVariable("PICKREP_PAR", par.createParameterString(par.pickrepprofile).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    // The temporary directory is provided as the 4th argument.
    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.pickconsensusrepfast));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    // Write out the embedded shell script to a file in the temporary directory.
    std::string program = tmpDir + "/pickconsensusrepfast.sh";
    FileUtil::writeFile(program, pickconsensusrepfast_sh, pickconsensusrepfast_sh_len);

    // Execute the shell script.
    cmd.execProgram(program.c_str(), par.filenames);

    // The shell script should not return; if it does, abort.
    assert(false);
    return 0;
}
