#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "easytaxonomy.sh.h"

void setEasyTaxonomyDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->removeTmpFiles = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->createdbMode = Parameters::SEQUENCE_SPLIT_MODE_SOFT;
    p->writeLookup = false;
    p->sensitivity = 5.7;
    p->evalThr = 1;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
    p->orfFilter = 0;
}
void setEasyTaxonomyMustPassAlong(Parameters *p) {
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_ORF_START_MODE.wasSet = true;
    p->PARAM_ORF_MIN_LENGTH.wasSet = true;
    p->PARAM_ORF_MAX_LENGTH.wasSet = true;
    p->PARAM_ORF_FILTER.wasSet = true;
}

int easytaxonomy(int argc, const char **argv, const Command& command) {
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
    for (size_t i = 0; i < par.convertalignments.size(); i++){
        par.convertalignments[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.createtsv.size(); i++){
        par.createtsv[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyTaxonomyDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyTaxonomyMustPassAlong(&par);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("TARGET", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    par.taxonomyOutpuMode = Parameters::TAXONOMY_OUTPUT_ALIGNMENT;
    par.PARAM_TAX_OUTPUT_MODE.wasSet = true;
    cmd.addVariable("TAXONOMY_PAR", par.createParameterString(par.taxonomy, true).c_str());
    cmd.addVariable("CREATEDB_QUERY_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("CREATETSV_PAR", par.createParameterString(par.createtsv).c_str());
    par.evalThr = 100000000;
    cmd.addVariable("SWAPRESULT_PAR", par.createParameterString(par.swapresult).c_str());
    FileUtil::writeFile(tmpDir + "/easy-taxonomy.sh", easytaxonomy_sh, easytaxonomy_sh_len);
    std::string program(tmpDir + "/easy-taxonomy.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
