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
void setTaxonomyMustPassAlong(Parameters *p) {
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_ORF_START_MODE.wasSet = true;
    p->PARAM_ORF_MIN_LENGTH.wasSet = true;
    p->PARAM_ORF_MAX_LENGTH.wasSet = true;
}

int taxonomy(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();

    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_NUM_ITERATIONS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_PICK_ID_FROM.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setTaxonomyDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);
    setTaxonomyMustPassAlong(&par);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.taxonomy));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    int alignmentMode = par.alignmentMode;
    if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA || par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA_APPROX) {
        // at least cov must be set for extractalignedregion
        int targetMode = (int)Parameters::ALIGNMENT_MODE_SCORE_COV;
        par.alignmentMode = std::max(par.alignmentMode, targetMode);
    }
    cmd.addVariable("SEARCH1_PAR", par.createParameterString(par.searchworkflow, true).c_str());
    par.alignmentMode = alignmentMode;

    if (par.taxonomySearchMode == Parameters::TAXONOMY_TOP_HIT) {
        cmd.addVariable("TOP_HIT", "1");
    }else if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA){
        par.sensSteps = 1;
        par.PARAM_SENS_STEPS.wasSet = true;
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(par.searchworkflow, true).c_str());
    }else if(par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA_APPROX){
        cmd.addVariable("APPROX_2BLCA", "1");
        cmd.addVariable("SEARCH2_PAR", par.createParameterString(par.align).c_str());
    }

    if (par.taxonomyOutpuMode == Parameters::TAXONOMY_OUTPUT_LCA) {
        cmd.addVariable("TAX_OUTPUT", "0" );
        cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    } else if (par.taxonomyOutpuMode == Parameters::TAXONOMY_OUTPUT_BOTH) {
        cmd.addVariable("TAX_OUTPUT", "2" );
        cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    } else {
        cmd.addVariable("TAX_OUTPUT", "1" );
    }

    FileUtil::writeFile(tmpDir + "/taxonomy.sh", taxonomy_sh, taxonomy_sh_len);
    std::string program(tmpDir + "/taxonomy.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
