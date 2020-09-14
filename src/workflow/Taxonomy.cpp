#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "lcaapprox2blca.sh.h"
#include "lcatophit.sh.h"

void setTaxonomyDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->sensitivity = 5.7;
    p->evalThr = 1;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
}
void setTaxonomyMustPassAlong(Parameters *p) {
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
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
    for (size_t i = 0; i < par.createdb.size(); i++) {
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++) {
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setTaxonomyDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);
    setTaxonomyMustPassAlong(&par);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.taxonomy));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA) {
        Debug(Debug::WARNING) << "2bLCA was replaced by Approx 2bLCA\n";
        par.taxonomySearchMode = Parameters::TAXONOMY_2BLCA_APPROX;
    }

    std::string program(tmpDir);
    int origAlignmentMode = par.alignmentMode;
    if (par.taxonomySearchMode == Parameters::TAXONOMY_TOP_HIT) {
        par.alignmentMode = std::max((int)Parameters::ALIGNMENT_MODE_SCORE_ONLY, par.alignmentMode);
        par.PARAM_ALIGNMENT_MODE.wasSet = true;
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());
        par.alignmentMode = origAlignmentMode;
        program.append("/lcatophit.sh");
        FileUtil::writeFile(program.c_str(), lcatophit_sh, lcatophit_sh_len);
    } else if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA_APPROX) {
        float origEvalue = par.evalThr;
        par.alignmentMode = Parameters::ALIGNMENT_MODE_UNGAPPED;
        par.PARAM_ALIGNMENT_MODE.wasSet = true;
        par.sortResults = true;
        par.PARAM_SORT_RESULTS.wasSet = true;
        par.evalThr = FLT_MAX;
        par.PARAM_E.wasSet = true;
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());
        par.evalThr = origEvalue;

        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        int origMaxAccept = par.maxAccept;
        par.maxAccept = 30;
        cmd.addVariable("ALN1_PAR", par.createParameterString(par.align).c_str());
        par.maxAccept = origMaxAccept;

        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
        cmd.addVariable("ALNTOP_PAR", par.createParameterString(par.align).c_str());

        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        cmd.addVariable("ALN2_PAR", par.createParameterString(par.align).c_str());
        par.alignmentMode = origAlignmentMode;

        program.append("/lcaapprox2blca.sh");
        FileUtil::writeFile(program.c_str(), lcaapprox2blca_sh, lcaapprox2blca_sh_len);
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

    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
