#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include <cassert>

#include "tsv2exprofiledb.sh.h"

int tsv2exprofiledb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, tsv2exprofiledb_sh, tsv2exprofiledb_sh_len);

    CommandCaller cmd;
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
