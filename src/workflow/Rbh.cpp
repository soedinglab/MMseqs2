#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
#include "rbh.sh.h"


void setRbhDefaults(Parameters *p) {
    p->compBiasCorrection = 0;
    p->maskMode = 0;
    p->orfStartMode = 1;
    p->orfMinLength = 10;
    p->orfMaxLength = 32734;
}

int rbh(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setRbhDefaults(&par);

    // set a lot of possibly misleading comments to EXPERT mode
    par.overrideParameterDescription((Command &) command, par.PARAM_OVERLAP.uniqid, NULL, NULL,
                                     par.PARAM_OVERLAP.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.extractorfs[i]->uniqid, NULL, NULL, par.extractorfs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.translatenucs[i]->uniqid, NULL, NULL, par.translatenucs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    // restore threads and verbosity
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, 4);

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
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

    std::string program = par.db4 + "/rbh.sh";
    FileUtil::writeFile(program, rbh_sh, rbh_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
