#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "taxonomy.sh.h"

void setTaxonomyWorkflowDefaults(Parameters *p) {
}

int taxonomy(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setTaxonomyWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 6);

    if(FileUtil::directoryExists(par.db6.c_str()) == false){
        Debug(Debug::ERROR) << "Temporary folder " << par.db6 << " does not exist or is not a directory!\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
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
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(par.searchworkflow).c_str());
    }

    if (par.lcaMode != Parameters::TAXONOMY_NO_LCA) {
        cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    }

    FileUtil::writeFile(par.db6 + "/taxonomy.sh", taxonomy_sh, taxonomy_sh_len);
    std::string program(par.db6 + "/taxonomy.sh");
    cmd.execProgram(program.c_str(), 6, argv);

    return EXIT_SUCCESS;
}
