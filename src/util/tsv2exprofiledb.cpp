#include "Parameters.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Debug.h"

#include <cassert>

#include "tsv2exprofiledb.sh.h"

void setTsv2ExProfileDbDefaults(Parameters *p) {
    p->compressed = true;
}

int tsv2exprofiledb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setTsv2ExProfileDbDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string program = par.db2 + ".sh";
    FileUtil::writeFile(program, tsv2exprofiledb_sh, tsv2exprofiledb_sh_len);

    if (par.gpu) {
        Debug(Debug::INFO) << "Disabling compression for GPU-databases\n";
        par.compressed = false;
    }

    CommandCaller cmd;
    cmd.addVariable("COMPRESSED", par.compressed ? "TRUE" : NULL);
    cmd.addVariable("GPU", par.gpu ? "TRUE" : NULL);
    cmd.addVariable("THREADS", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.execProgram(FileUtil::getRealPathFromSymLink(program).c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
