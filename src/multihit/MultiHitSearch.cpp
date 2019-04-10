#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "multihitsearch.sh.h"

void setMultiHitSearchWorkflowDefaults(Parameters *p) {
    p->sensitivity = 7;
    p->covThr = 0.7;
    p->evalThr = 100;

    // TODO: Needs to be more than the count of target sets (10x?)
    p->maxSequences = 1500;
    p->scoreBias = 0.3;

    p->simpleBestHit = false;
}

int multihitsearch(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMultiHitSearchWorkflowDefaults(&par);

    par.overrideParameterDescription((Command &) command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL,
                                     par.PARAM_MAX_REJECTED.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_OVERLAP.uniqid, NULL, NULL,
                                     par.PARAM_OVERLAP.category | MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.extractorfs[i]->uniqid, NULL, NULL, par.extractorfs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.translatenucs[i]->uniqid, NULL, NULL, par.translatenucs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.overrideParameterDescription((Command &)command, par.result2profile[i]->uniqid, NULL, NULL, par.result2profile[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, 4, true);

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.multihitsearch);
    std::string tmpDir = par.db4 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    par.splitSeqByLen = false;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    cmd.addVariable("BESTHITBYSET_PAR", par.createParameterString(par.besthitbyset).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(tmpDir + "/multihitsearch.sh", multihitsearch_sh, multihitsearch_sh_len);
    std::string program(tmpDir + "/multihitsearch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
