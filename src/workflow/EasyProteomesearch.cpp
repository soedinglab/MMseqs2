#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

#include "easyproteomesearch.sh.h"

void setEasyproteomesearchDefaults(Parameters *p) {
    p->removeTmpFiles = true; 
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->covThr = 0.8;
}

void setEasyproteomesearchMustPassAlong(Parameters *p){
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_C.wasSet = true;
}

int easyproteomesearch(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEasyproteomesearchDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyproteomesearchMustPassAlong(&par);
    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string target = par.filenames.back().c_str();
    cmd.addVariable("TARGET", target.c_str());
    par.filenames.pop_back();

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);
    
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb,true).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());
    cmd.addVariable("PARSEPROTEOMEALN_PAR", par.createParameterString(par.parseproteomealignments).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyproteomesearch.sh";
    FileUtil::writeFile(program, easyproteomesearch_sh, easyproteomesearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    
    // Should never get here
    assert(false);
    return 0;
}