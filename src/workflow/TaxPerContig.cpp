#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"
#include "taxpercontig.sh.h"

void setTaxPerContigDefaults(Parameters *p) {
    p->orfStartMode = 1;
    p->translate = 1;
    p->taxonomyOutpuMode = 2;
    p->showTaxLineage = 0;
}

void setTaxPerContigMustPassAlong(Parameters *p) {
    p->PARAM_ORF_START_MODE.wasSet = true;
    p->PARAM_TRANSLATE.wasSet = true;
    p->PARAM_TAX_OUTPUT_MODE.wasSet = true;
    p->PARAM_TAXON_ADD_LINEAGE.wasSet = true;
}

int taxpercontig(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setTaxPerContigDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);
    setTaxPerContigMustPassAlong(&par);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    par.translate = 1;
    cmd.addVariable("EXTRACT_ORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    int showTaxLineageOrig = par.showTaxLineage;
    // never show lineage for the orfs
    par.showTaxLineage = 0;
    par.translate = 0;
    cmd.addVariable("TAXONOMY_PAR", par.createParameterString(par.taxonomy, true).c_str());
    par.showTaxLineage = showTaxLineageOrig;
    cmd.addVariable("AGGREGATETAX_PAR", par.createParameterString(par.aggregatetax).c_str());
    cmd.addVariable("VERBOSITY_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());

    std::string program(tmpDir + "/taxpercontig.sh");
    FileUtil::writeFile(program, taxpercontig_sh, taxpercontig_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // should never get here
    return EXIT_FAILURE;
}
