#include "Parameters.h"
#include <string>
#include <cassert>
#include <Util.h>
#include "linclust.sh.h"

#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"


void setLinclustWorkflowDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->evalThr = 0.001;
    p->seqIdThr = 0.9;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int linclust(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setLinclustWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 3);
    if(FileUtil::directoryExists(par.db3.c_str())==false){
        Debug(Debug::ERROR) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        EXIT(EXIT_FAILURE);
    }
    bool targetCov = false;
    bool cov = false;
    for (size_t i = 0; i < par.linclustworkflow.size(); i++) {
        if (par.linclustworkflow[i].uniqid == par.PARAM_TARGET_COV.uniqid && par.linclustworkflow[i].wasSet) {
            if(par.targetCovThr > 0.0 ){
                targetCov = true;
                par.covThr = 0.0;
            }
        }
        if (par.linclustworkflow[i].uniqid == par.PARAM_C.uniqid && par.linclustworkflow[i].wasSet) {
            cov = true;
        }
    }
    if(cov && targetCov){
        Debug(Debug::ERROR) << "The paramter -c can not be combined with --target-cov.\n";
        EXIT(EXIT_FAILURE);
    }
    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());

    // save some values to restore them later
    size_t alphabetSize = par.alphabetSize;
    size_t kmerSize = par.kmerSize;
    // # 1. Finding exact $k$-mer matches.
    par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    par.alphabetSize = alphabetSize;
    par.kmerSize = kmerSize;
    // # 2. Hamming distance pre-clustering
    par.rescoreMode = Parameters::RESCORE_MODE_HAMMING;
    par.filterHits = false;
    float prevSeqId = par.seqIdThr;
    // hamming distance does not work well with seq. id < 0.5 since it does not have an e-value criteria
    par.seqIdThr = std::max(0.5f, par.seqIdThr);
    // also coverage should not be under 0.5
    float prevCov = par.covThr;
    float prevTargetCov = par.targetCovThr;
    if(cov){
        par.covThr = std::max(0.5f, par.covThr);
    }
    if(targetCov){
        par.targetCovThr = std::max(0.5f, par.targetCovThr);
    }
    cmd.addVariable("HAMMING_PAR", par.createParameterString(par.rescorediagonal).c_str());
    // set it back to old value
    par.covThr = prevCov;
    par.targetCovThr = prevTargetCov;
    par.seqIdThr = prevSeqId;
    par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
    // # 3. Ungapped alignment filtering
    par.filterHits = true;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    // # 4. Local gapped sequence alignment.
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    // # 5. Clustering using greedy set cover.
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    FileUtil::writeFile(par.db3 + "/linclust.sh", linclust_sh, linclust_sh_len);
    std::string program(par.db3 + "/linclust.sh");
    cmd.execProgram(program.c_str(), 3, argv);
    return 0;
}
