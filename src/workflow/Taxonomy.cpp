#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "taxonomy.sh.h"


void setTaxonomyDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 1;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
}

int taxonomy(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setTaxonomyDefaults(&par);
    par.parseParameters(argc, argv, command, 4);

    if(FileUtil::directoryExists(par.db4.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(par.filenames, par.taxonomy));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
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

    if (par.lcaMode == Parameters::TAXONOMY_2BLCA){
        par.sensSteps = 1;
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(par.searchworkflow).c_str());
    }else if(par.lcaMode == Parameters::TAXONOMY_2BLCA_APPROX){
        cmd.addVariable("APPROX_2BLCA", "1");
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(par.align).c_str());
    }

    if (par.lcaMode != Parameters::TAXONOMY_NO_LCA) {
        cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    }

    FileUtil::writeFile(tmpDir + "/taxonomy.sh", taxonomy_sh, taxonomy_sh_len);
    std::string program(tmpDir + "/taxonomy.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
