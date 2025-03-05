#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
#include "rbh.sh.h"

#include <cassert>

void setRbhDefaults(Parameters *p) {
    p->compBiasCorrection = 0;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->maskMode = 0;
    p->orfStartMode = 1;
    p->orfMinLength = 10;
    p->orfMaxLength = 32734;
}

int rbh(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setRbhDefaults(&par);

    // set a lot of possibly misleading comments to EXPERT mode
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    // restore threads and verbosity
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);


    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("SEARCH_A_B_PAR", par.createParameterString(par.searchworkflow).c_str());
    int originalCovMode = par.covMode;
    par.covMode = Util::swapCoverageMode(par.covMode);
    cmd.addVariable("SEARCH_B_A_PAR", par.createParameterString(par.searchworkflow).c_str());
    par.covMode = originalCovMode;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERB_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    std::string program = tmpDir + "/rbh.sh";
    FileUtil::writeFile(program, rbh_sh, rbh_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
