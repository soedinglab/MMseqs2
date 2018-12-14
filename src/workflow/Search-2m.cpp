#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "search2m.sh.h"


void setSearch2MDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    p->orfStartMode = 0;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
}

int search2m(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setSearch2MDefaults(&par);
    par.parseParameters(argc, argv, command, 4);

    if(FileUtil::directoryExists(par.db4.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::ERROR) << "Could not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }


    std::string hash = SSTR(par.hashParameter(par.filenames, par.searchworkflow));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("SEARCH1_PAR", par.createParameterString(par.searchworkflow).c_str());
    std::vector<MMseqsParameter*> searchNoIterativeBest;
    for (size_t i = 0; i < par.searchworkflow.size(); i++){
        if (par.searchworkflow[i]->uniqid != par.PARAM_START_SENS.uniqid
            || par.searchworkflow[i]->uniqid != par.PARAM_SENS_STEPS.uniqid) {
            searchNoIterativeBest.push_back(par.searchworkflow[i]);
        }
    }
    cmd.addVariable("SEARCH2_PAR", par.createParameterString(searchNoIterativeBest).c_str());


    FileUtil::writeFile(tmpDir + "/search2m.sh", search2m_sh, search2m_sh_len);
    std::string program(tmpDir + "/search2m.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
