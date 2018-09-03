#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

#include "easysearch.sh.h"


void setEasySearchDefaults(Parameters *p) {
    p->sensitivity = 5.7;
    p->removeTmpFiles = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}


int easysearch(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEasySearchDefaults(&par);

    par.overrideParameterDescription((Command &) command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL,
                                     par.PARAM_MAX_REJECTED.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_OVERLAP.uniqid, NULL, NULL,
                                     par.PARAM_OVERLAP.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.extractorfs[i].uniqid, NULL, NULL, par.extractorfs[i].category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.translatenucs[i].uniqid, NULL, NULL, par.translatenucs[i].category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.overrideParameterDescription((Command &)command, par.result2profile[i].uniqid, NULL, NULL, par.result2profile[i].category | MMseqsParameter::COMMAND_EXPERT);
    }

    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, 4);

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not crate tmp folder " << par.db4 << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }

    size_t hash = par.hashParameter(par.filenames, par.easysearchworkflow);

    std::string tmpDir = par.db4 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult).c_str());

    FileUtil::writeFile(tmpDir + "/easysearch.sh", easysearch_sh, easysearch_sh_len);
    std::string program(tmpDir + "/easysearch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
