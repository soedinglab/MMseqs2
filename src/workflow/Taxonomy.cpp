#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "taxonomy.sh.h"

int taxonomy(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 6);

    if(FileUtil::directoryExists(par.db6.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db6 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db6.c_str()) == false){
            Debug(Debug::ERROR) << "Could not crate tmp folder " << par.db6 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db6 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.taxonomy);
    std::string tmpDir = par.db6+"/"+SSTR(hash);
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

    int alignmentMode = par.alignmentMode;
    if (par.lcaMode == Parameters::TAXONOMY_2BLCA) {
        // at least cov must be set for extractalignedregion
        int targetMode = (int)Parameters::ALIGNMENT_MODE_SCORE_COV;
        par.alignmentMode = std::max(par.alignmentMode, targetMode);
    }
    cmd.addVariable("SEARCH1_PAR", par.createParameterString(par.searchworkflow).c_str());
    par.alignmentMode = alignmentMode;

    if (par.lcaMode == Parameters::TAXONOMY_2BLCA) {
        std::vector<MMseqsParameter> searchNoIterativeBest;
        for (size_t i = 0; i < par.searchworkflow.size(); i++){
            if (par.searchworkflow[i].uniqid != par.PARAM_START_SENS.uniqid
             || par.searchworkflow[i].uniqid != par.PARAM_SENS_STEPS.uniqid) {
                searchNoIterativeBest.push_back(par.searchworkflow[i]);
            }
        }
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(searchNoIterativeBest).c_str());
    }

    if (par.lcaMode != Parameters::TAXONOMY_NO_LCA) {
        cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    }

    FileUtil::writeFile(tmpDir + "/taxonomy.sh", taxonomy_sh, taxonomy_sh_len);
    std::string program(tmpDir + "/taxonomy.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
