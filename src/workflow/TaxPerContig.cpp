#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "PrefilteringIndexReader.h"
#include "taxpercontig.sh.h"

void setTaxPerContigDefaults(Parameters *p) {
    p->orfStartMode = 1;
    p->showTaxLineage = 0;
    p->taxonomySearchMode = Parameters::TAXONOMY_2BLCA_APPROX;
    p->orfFilter = true;
}

void setTaxPerContigMustPassAlong(Parameters *p) {
    p->PARAM_ORF_START_MODE.wasSet = true;
    p->PARAM_TAXON_ADD_LINEAGE.wasSet = true;
    p->PARAM_LCA_MODE.wasSet = true;
}

int taxpercontig(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALIGNMENT_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_NUM_ITERATIONS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEARCH_TYPE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
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
    setTaxPerContigDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);
    setTaxPerContigMustPassAlong(&par);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);
    cmd.addVariable("TARGETDB_IDX", (indexStr == "") ? par.db2.c_str() : indexStr.c_str());
    cmd.addVariable("RUNNER", par.runner.c_str());
    par.translate = 1;
    cmd.addVariable("EXTRACT_ORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    int showTaxLineageOrig = par.showTaxLineage;
    // never show lineage for the orfs
    par.showTaxLineage = 0;
    par.taxonomyOutpuMode = 2;
    par.PARAM_TAX_OUTPUT_MODE.wasSet = true;
    cmd.addVariable("TAXONOMY_PAR", par.createParameterString(par.taxonomy, true).c_str());
    par.showTaxLineage = showTaxLineageOrig;
    cmd.addVariable("AGGREGATETAX_PAR", par.createParameterString(par.aggregatetax).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());

    cmd.addVariable("ORF_FILTER", par.orfFilter && par.orfFilterSens <= par.sensitivity ? "TRUE" : NULL);
    // --min-ungapped-score 3 -s 2 --diag-score 0 --max-seqs 1
    par.minDiagScoreThr = 3;
    par.sensitivity = par.orfFilterSens;
    par.diagonalScoring = false;
    par.maxResListLen = 1;
    cmd.addVariable("ORF_FILTER_PREFILTER", par.createParameterString(par.prefilter).c_str());

    par.evalThr = par.orfFilterEval;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("ORF_FILTER_RESCOREDIAGONAL", par.createParameterString(par.rescorediagonal).c_str());

    par.subDbMode = Parameters::SUBDB_MODE_SOFT;
    cmd.addVariable("CREATESUBDB_PAR", par.createParameterString(par.createsubdb).c_str());

    std::string program(tmpDir + "/taxpercontig.sh");
    FileUtil::writeFile(program, taxpercontig_sh, taxpercontig_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
