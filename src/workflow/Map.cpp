#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "map.sh.h"

#include <cassert>

void setMapWorkflowDefaults(Parameters *p) {
    p->compBiasCorrection = 0;
    p->maskMode = 0;
    p->covThr = 0.95;
    p->covMode = 2;
    p->seqIdThr = 0.9;
    p->sensitivity = 2;
    p->rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    p->sortResults = true;
    //p->orfLongest = true;
    p->orfStartMode = 1;
    p->orfMinLength = 10;
    p->orfMaxLength = 32734;
}

int map(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMapWorkflowDefaults(&par);

    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.mapworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("RUNNER", par.runner.c_str());

    par.mapworkflow.push_back(&(par.PARAM_ALIGNMENT_MODE));
    par.alignmentMode = 4;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.mapworkflow).c_str());

    std::string program = tmpDir + "/map.sh";
    FileUtil::writeFile(program, map_sh, map_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
