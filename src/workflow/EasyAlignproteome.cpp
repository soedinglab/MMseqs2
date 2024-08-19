#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

#include "easyalignproteome.sh.h"

void setEasyAlignproteomeDefaults(Parameters *p) {
    p->proteomeSimThr = 0.9;
}

// void setEasyAlignproteomeMustPassAlong(Parameters *p){

// }

int easyalignproteome(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEasyAlignproteomeDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    // setEasyAlignproteomeMustPassAlong(&par);
    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("ALIGNPROTEOME_PAR", par.createParameterString(par.alignproteome,true).c_str()); // what?
    std::string program = tmpDir + "/easyalignproteome.sh";
    FileUtil::writeFile(program, easyalignproteome_sh, easyalignproteome_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    
    // Should never get here
    assert(false);
    return 0;
}
