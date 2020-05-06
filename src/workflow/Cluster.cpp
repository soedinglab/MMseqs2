#include "Parameters.h"
#include "Util.h"
#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

#include "cascaded_clustering.sh.h"
#include "nucleotide_clustering.sh.h"
#include "clustering.sh.h"

#include <cassert>

void setWorkflowDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.001;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->maxResListLen = 20;
}

float setAutomaticThreshold(float seqId) {
    float sens;
    if (seqId <= 0.3) {
        sens = 6;
    } else if (seqId > 0.8) {
        sens = 1.0;
    } else {
        sens = 1.0 + (1.0 * (0.7 - seqId) * 10);
    }
    return sens;
}

int setAutomaticIterations(float sens){
    if (sens <= 2.0) {
        return 1;
    } else {
        return 3;
    }
}


void setNuclClusterDefaults(Parameters *p) {
    // leave ungapped alignment untouched
    if(p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED){
        p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    //p->orfLongest = true;
    p->exactKmerMatching = true;
    if ( p->PARAM_DIAGONAL_SCORING.wasSet == false) {
        p->diagonalScoring = 0;
    }
    if ( p->PARAM_STRAND.wasSet == false) {
        p->strand = 2;
    }
    if ( p->PARAM_K.wasSet == false) {
        p->kmerSize = 15;
    }
    if (  p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
}



int clusteringworkflow(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setWorkflowDefaults(&par);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_S.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    const int dbType = FileUtil::parseDbType(par.db1.c_str());
    bool isNucleotideDb = (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES));
    if(isNucleotideDb){
        setNuclClusterDefaults(&par);
    }
    bool sensitivitySet = par.PARAM_S.wasSet;
    bool compositionBiasSet = par.PARAM_NO_COMP_BIAS_CORR.wasSet;
    bool clusterModeSet = par.PARAM_CLUSTER_MODE.wasSet;
    bool clusterStepsSet = par.PARAM_CLUSTER_STEPS.wasSet;
    bool minDiagonalScoreSet = par.PARAM_MIN_DIAG_SCORE.wasSet;

    if (compositionBiasSet == false) {
        if(par.seqIdThr >= 0.7){
            par.compBiasCorrection = 0;
        }
    }

    if (minDiagonalScoreSet == false) {
        if(par.seqIdThr >= 0.7){
            par.minDiagScoreThr = 60;
        }
    }

    if (sensitivitySet == false && isNucleotideDb == false) {
        par.sensitivity = setAutomaticThreshold(par.seqIdThr);
        Debug(Debug::INFO) << "Set cluster sensitivity to -s " << par.sensitivity << "\n";
    }

    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE)) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases.\n";
        EXIT(EXIT_FAILURE);
    }


    const bool nonSymetric = (par.covMode == Parameters::COV_MODE_TARGET ||par.covMode == Parameters::COV_MODE_QUERY);
    if (clusterModeSet == false) {
        if (nonSymetric) {
            par.clusteringMode = Parameters::GREEDY_MEM;
        } else {
            par.clusteringMode = Parameters::SET_COVER;
        }
        std::string cluMode = (par.clusteringMode == Parameters::GREEDY_MEM) ? "GREEDY MEM" : "SET COVER";
        Debug(Debug::INFO) << "Set cluster mode " << cluMode << ".\n";
    }
    if (nonSymetric && par.clusteringMode != Parameters::GREEDY && par.clusteringMode != Parameters::GREEDY_MEM) {
        Debug(Debug::WARNING) << "Combining cluster mode " << par.clusteringMode
                              << " in combination with coverage mode " << par.covMode << " can produce wrong results.\n"
                              << "Please use --cov-mode 2\n";
    }
    if (par.singleStepClustering == false && par.clusteringMode == Parameters::CONNECTED_COMPONENT) {
        Debug(Debug::WARNING) << "Connected component clustering produces less clusters in a single step clustering.\n"
                              << "Please use --single-step-cluster";
    }
    if (clusterStepsSet == false) {
        par.clusterSteps = setAutomaticIterations(par.sensitivity);
        Debug(Debug::INFO) << "Set cluster iterations to " << par.clusterSteps << "\n";
    }

    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.clusterworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    const int originalRescoreMode = par.rescoreMode;

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
    par.rescoreMode = originalRescoreMode;
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    if(isNucleotideDb){
        par.forwardFrames= "1";
        par.reverseFrames= "1";
        par.searchType = 3;
        cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());
        int oldKmer = par.kmerSize;
        par.kmerSize = 0;
        cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
        par.kmerSize = oldKmer;
        if (par.PARAM_MAX_SEQS.wasSet == false) {
            par.maxResListLen = 300;
        }

        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable("ALIGNMENT_MODE_NOT_SET","TRUE");
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("RESCORE_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
            cmd.addVariable("THREADSANDCOMPRESS_PAR", par.createParameterString(par.threadsandcompression).c_str());
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        }
        cmd.addVariable("CLUSTER_PAR",   par.createParameterString(par.clust).c_str());
        cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment).c_str());
        std::string program = tmpDir + "/nucleotide_clustering.sh";
        FileUtil::writeFile(program, nucleotide_clustering_sh, nucleotide_clustering_sh_len);
        cmd.execProgram(program.c_str(), par.filenames);
    } else if (par.singleStepClustering == false) {
        // save some values to restore them later
        float targetSensitivity = par.sensitivity;
        MultiParam<int> alphabetSize = par.alphabetSize;
        par.alphabetSize = MultiParam<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5);
        int kmerSize = par.kmerSize;
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
        int maskMode = par.maskMode;
        par.maskMode = 0;
        cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
        par.alphabetSize = alphabetSize;
        par.kmerSize = kmerSize;
        par.maskMode = maskMode;
        // 1 is lowest sens
        par.sensitivity = ((par.clusterSteps - 1) == 0 ) ? par.sensitivity  : 1;
        int minDiagScoreThr = par.minDiagScoreThr;
        par.minDiagScoreThr = 0;
        par.diagonalScoring = 0;
        par.compBiasCorrection = 0;
        cmd.addVariable("PREFILTER0_PAR", par.createParameterString(par.prefilter).c_str());
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("ALIGNMENT0_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable("ALIGNMENT0_PAR", par.createParameterString(par.align).c_str());
        }
        cmd.addVariable("CLUSTER0_PAR",   par.createParameterString(par.clust).c_str());
        par.diagonalScoring = 1;
        par.compBiasCorrection = 1;
        par.minDiagScoreThr = minDiagScoreThr;
        float sensStepSize = (targetSensitivity - 1)/ (static_cast<float>(par.clusterSteps)-1);
        for(int step = 1; step < par.clusterSteps; step++){
            par.sensitivity =  1.0 + sensStepSize * step;

            cmd.addVariable(std::string("PREFILTER"+SSTR(step)+"_PAR").c_str(), par.createParameterString(par.prefilter).c_str());
            if (isUngappedMode) {
                par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
                cmd.addVariable(std::string("ALIGNMENT"+SSTR(step)+"_PAR").c_str(), par.createParameterString(par.rescorediagonal).c_str());
                par.rescoreMode = originalRescoreMode;
            } else {
                cmd.addVariable(std::string("ALIGNMENT"+SSTR(step)+"_PAR").c_str(), par.createParameterString(par.align).c_str());
            }
            cmd.addVariable(std::string("CLUSTER"  +SSTR(step)+"_PAR").c_str(), par.createParameterString(par.clust).c_str());
        }
        cmd.addVariable("STEPS", SSTR(par.clusterSteps).c_str());
        // correct for cascading clustering errors
        if(par.clusterReassignment){
            cmd.addVariable("REASSIGN","TRUE");
        }
        cmd.addVariable("THREADSANDCOMPRESS", par.createParameterString(par.threadsandcompression).c_str());
        cmd.addVariable("VERBCOMPRESS", par.createParameterString(par.verbandcompression).c_str());
        cmd.addVariable("ALIGNMENT_REASSIGN_PAR", par.createParameterString(par.align).c_str());

        std::string program = tmpDir + "/cascaded_clustering.sh";
        FileUtil::writeFile(program, cascaded_clustering_sh, cascaded_clustering_sh_len);
        cmd.execProgram(program.c_str(), par.filenames);
    } else {
        // same as above, clusthash needs a smaller alphabetsize
        MultiParam<int> alphabetSize = par.alphabetSize;
        par.alphabetSize = MultiParam<int> (Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5);
        float seqIdThr = par.seqIdThr;
        par.seqIdThr = (float)Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID/100.0f;
        cmd.addVariable("DETECTREDUNDANCY_PAR", par.createParameterString(par.clusthash).c_str());
        par.alphabetSize = alphabetSize;
        par.seqIdThr = seqIdThr;
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        if (isUngappedMode) {
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
        } else {
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        }
        cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
        FileUtil::writeFile(tmpDir + "/clustering.sh", clustering_sh, clustering_sh_len);
        std::string program(tmpDir+ "/clustering.sh");
        cmd.execProgram(program.c_str(), par.filenames);
    }

    // Unreachable
    assert(false);
    return 0;
}
