#include "Parameters.h"
#include <string>
#include <cassert>
#include <Util.h>
#include "assembler.sh.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"


void setAssemblerWorkflowDefaults(Parameters *p) {
    p->spacedKmer = false;
    p->maskMode = 1;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.95;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int assembler(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 3);
    if(FileUtil::directoryExists(par.db3.c_str())==false){
        Debug(Debug::ERROR) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        EXIT(EXIT_FAILURE);
    }
    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    // save some values to restore them later
    size_t alphabetSize = par.alphabetSize;
    size_t kmerSize = par.kmerSize;
    // # 1. Finding exact $k$-mer matches.
    int baseKmerSize = 14;
    par.kmerSize = baseKmerSize;
    par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    for(size_t i = 0; i < par.numIterations; i++){
        std::string key = "KMERMATCHER"+SSTR(i)+"_PAR";
    //    par.kmerSize = baseKmerSize - i;
        cmd.addVariable(key.c_str(), par.createParameterString(par.kmermatcher).c_str());
    }
    par.alphabetSize = alphabetSize;
    par.kmerSize = kmerSize;
    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    FileUtil::writeFile(par.db3 + "/assembler.sh", assembler_sh, assembler_sh_len);
    std::string program(par.db3 + "/assembler.sh");
    cmd.execProgram(program.c_str(), 3, argv);
    return 0;
}
