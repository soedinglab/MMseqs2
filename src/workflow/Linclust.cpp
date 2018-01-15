#include "Parameters.h"
#include <string>
#include <cassert>
#include <Util.h>
#include <Sequence.h>
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
    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL, par.PARAM_RESCORE_MODE.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL, par.PARAM_MAX_REJECTED.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_ACCEPT.uniqid, NULL, NULL, par.PARAM_MAX_ACCEPT.category |MMseqsParameter::COMMAND_EXPERT );

    par.parseParameters(argc, argv, command, 3);

    int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());

    if(FileUtil::directoryExists(par.db3.c_str())==false){
        Debug(Debug::WARNING) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db3.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db3 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::WARNING) << "Created dir " << par.db3 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.linclustworkflow);
    std::string tmpDir = par.db3+"/"+SSTR(hash);
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    if(FileUtil::symlinkCreateOrRepleace(par.db3+"/latest", tmpDir) == false){
        Debug(Debug::WARNING) << "Could not link latest folder in tmp." << tmpDir << ".\n";
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
    bool kmerSizeWasSet = false;
    bool alphabetSizeWasSet = false;
    bool clusterModeSet = false;
    for (size_t i = 0; i < par.linclustworkflow.size(); i++) {
        if (par.linclustworkflow[i].uniqid == par.PARAM_K.uniqid && par.linclustworkflow[i].wasSet) {
            kmerSizeWasSet = true;
        }
        if (par.linclustworkflow[i].uniqid == par.PARAM_ALPH_SIZE.uniqid && par.linclustworkflow[i].wasSet) {
            alphabetSizeWasSet = true;
        }
        if (par.linclustworkflow[i].uniqid == par.PARAM_CLUSTER_MODE.uniqid && par.linclustworkflow[i].wasSet) {
            clusterModeSet = true;
        }
    }

    bool noneSymetric = (par.covMode == Parameters::COV_MODE_TARGET ||par.covMode == Parameters::COV_MODE_QUERY);
    if(clusterModeSet == false){
        if(noneSymetric){
            par.clusteringMode = Parameters::GREEDY;
        }else{
            par.clusteringMode = Parameters::SET_COVER;
        }
        std::string cluMode = (par.clusteringMode==Parameters::GREEDY) ? "GREEDY" : "SET COVER";
        Debug(Debug::WARNING) << "Set cluster mode " << cluMode << ".\n";
    }

    if(kmerSizeWasSet==false){
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    }
    if(alphabetSizeWasSet == false){
        par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    }

    // filter by diagonal in case of AA (do not filter for nucl, profiles, ...)
    if(queryDbType == Sequence::AMINO_ACIDS){
        cmd.addVariable("FILTER","1");
    }
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
    par.covThr = std::max(0.5f, par.covThr);
    cmd.addVariable("HAMMING_PAR", par.createParameterString(par.rescorediagonal).c_str());
    // set it back to old value
    par.covThr = prevCov;
    par.seqIdThr = prevSeqId;
    par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
    // # 3. Ungapped alignment filtering
    par.filterHits = true;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    // # 4. Local gapped sequence alignment.
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    // # 5. Clustering using greedy set cover.
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    FileUtil::writeFile(tmpDir + "/linclust.sh", linclust_sh, linclust_sh_len);
    std::string program(tmpDir + "/linclust.sh");
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}
