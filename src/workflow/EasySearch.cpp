#include <cassert>
#include "LinsearchIndexReader.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
#include "easysearch.sh.h"

void setEasyLinsearchDefaults(Parameters *p) {
    p->shuffleDatabase = false;
    p->removeTmpFiles = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}

void setEasySearchDefaults(Parameters *p) {
    p->sensitivity = 5.7;
    p->removeTmpFiles = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}

int doeasysearch(int argc, const char **argv, const Command &command, bool linsearch) {
    Parameters &par = Parameters::getInstance();
    if(linsearch){
        setEasyLinsearchDefaults(&par);
    }else{
        setEasySearchDefaults(&par);
    }
    par.overrideParameterDescription((Command &) command, par.PARAM_ADD_BACKTRACE.uniqid, NULL, NULL, par.PARAM_ADD_BACKTRACE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL,
                                     par.PARAM_MAX_REJECTED.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_OVERLAP.uniqid, NULL, NULL,
                                     par.PARAM_OVERLAP.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL,
                                     par.PARAM_RESCORE_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.overrideParameterDescription((Command &)command, par.createdb[i]->uniqid, NULL, NULL, par.createdb[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
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
    par.parseParameters(argc, argv, command, 4);


    bool needBacktrace = false;
    {
        bool needSequenceDB = false;
        bool needFullHeaders = false;
        Parameters::getOutputFormat(par.outfmt, needSequenceDB, needBacktrace, needFullHeaders);
    }
    if(par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM){
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }

    std::string hash = SSTR(par.hashParameter(par.filenames, *command.params));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    if (linsearch) {
        const bool isIndex = LinsearchIndexReader::searchForIndex(par.db2).empty() == false;
        cmd.addVariable("INDEXEXT", isIndex ? ".linidx" : NULL);
        cmd.addVariable("SEARCH_MODULE", "linsearch");
        cmd.addVariable("LINSEARCH", "TRUE");
        cmd.addVariable("CREATELININDEX_PAR", par.createParameterString(par.createlinindex).c_str());
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.linsearchworkflow, true).c_str());
    } else {
        const bool isIndex = PrefilteringIndexReader::searchForIndex(par.db2).empty() == false;
        cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);
        cmd.addVariable("SEARCH_MODULE", "search");
        cmd.addVariable("LINSEARCH", NULL);
        cmd.addVariable("CREATELININDEX_PAR", NULL);
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());

    }
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult).c_str());

    std::string program = tmpDir + "/easysearch.sh";
    FileUtil::writeFile(program, easysearch_sh, easysearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}

int easysearch(int argc, const char **argv, const Command &command) {
    return doeasysearch(argc, argv, command, false);
}

int easylinsearch(int argc, const char **argv, const Command &command) {
    return doeasysearch(argc, argv, command, true);
}
