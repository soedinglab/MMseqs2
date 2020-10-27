#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "taxonomy.sh.h"

void setTaxonomyDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->sensitivity = 2;
    p->evalThr = 1;
    p->maxAccept = 30;
    p->maxRejected = 5;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
}
void setTaxonomyMustPassAlong(Parameters *p) {
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_MAX_ACCEPT.wasSet = true;
    p->PARAM_MAX_REJECTED.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
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
    if (par.taxonomySearchMode == Parameters::TAXONOMY_TOP_HIT) {
        cmd.addVariable("TOPHIT_MODE", "1");
    } else if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA_APPROX) {
        par.lcaSearch = true;
        par.PARAM_LCA_SEARCH.wasSet = true;
        cmd.addVariable("TOPHIT_MODE", NULL);
    }
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());

    program.append("/taxonomy.sh");
    FileUtil::writeFile(program.c_str(), taxonomy_sh, taxonomy_sh_len);

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
