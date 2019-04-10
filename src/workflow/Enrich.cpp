#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "enrich.sh.h"

void setEnrichWorkflowDefaults(Parameters *p) {
    p->numIterations = 3;
    p->expansionMode = 1;
}

int enrich(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setEnrichWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 6);

    if (FileUtil::directoryExists(par.db6.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db6 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db6.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not crate tmp folder " << par.db6 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db6 << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(par.filenames, par.enrichworkflow));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db6+"/latest");
    }
    std::string tmpDir = par.db6+"/"+hash;
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
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    par.addBacktrace = true;

    int originalNumIterations = par.numIterations;
    par.numIterations = 1;
    par.sliceSearch = true;
    cmd.addVariable("PROF_SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    par.sliceSearch = false;
    par.numIterations = originalNumIterations;


    cmd.addVariable("PROF_PROF_PAR", par.createParameterString(par.result2profile).c_str());
    cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    // change once rescorediagonal supports profiles
    const bool isUngappedMode = false;
    cmd.addVariable("ALIGN_MODULE", "align");

    float originalEval = par.evalThr;
    par.evalThr = par.evalProfile;
    par.realign = false;
    for (int i = 0; i < par.numIterations; i++) {
        if (i == (par.numIterations - 1)) {
            par.evalThr = originalEval;
        }

        cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.prefilter).c_str());
        if (isUngappedMode) {
            int originalRescoreMode = par.rescoreMode;
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.align).c_str());
        }

        cmd.addVariable(std::string("EXPAND_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.expandaln).c_str());

        par.pca = 0.0;
        cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.result2profile).c_str());
        par.pca = 1.0;
    }

    std::string program(tmpDir + "/enrich.sh");
    FileUtil::writeFile(program, enrich_sh, enrich_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
