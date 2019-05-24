#include "Parameters.h"
#include "Util.h"
#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

#include "cascaded_clustering.sh.h"
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

int clusteringworkflow(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setWorkflowDefaults(&par);
    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL, par.PARAM_RESCORE_MODE.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_REJECTED.uniqid, NULL, NULL, par.PARAM_MAX_REJECTED.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_MAX_ACCEPT.uniqid, NULL, NULL, par.PARAM_MAX_ACCEPT.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_KMER_PER_SEQ.uniqid, NULL, NULL, par.PARAM_KMER_PER_SEQ.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_S.uniqid, "sensitivity will be automatically determined but can be adjusted", NULL,  par.PARAM_S.category |MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_INCLUDE_ONLY_EXTENDABLE.uniqid, NULL, NULL, par.PARAM_INCLUDE_ONLY_EXTENDABLE.category |MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, 3);
    if (FileUtil::directoryExists(par.db3.c_str())==false) {
        Debug(Debug::INFO) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db3.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db3 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db3 << "\n";
        }
    }

    bool sensitivitySet = false;
    bool compositionBiasSet = false;
    bool clusterModeSet = false;
    bool clusterStepsSet = false;
    bool minDiagonalScoreSet = false;

    for (size_t i = 0; i < par.clusterworkflow.size(); i++) {
        if (par.clusterworkflow[i]->uniqid == par.PARAM_S.uniqid && par.clusterworkflow[i]->wasSet) {
            sensitivitySet = true;
        }
        if (par.clusterworkflow[i]->uniqid == par.PARAM_CLUSTER_MODE.uniqid && par.clusterworkflow[i]->wasSet) {
            clusterModeSet = true;
        }
        if (par.clusterworkflow[i]->uniqid == par.PARAM_CLUSTER_STEPS.uniqid && par.clusterworkflow[i]->wasSet) {
            clusterStepsSet = true;
        }
        if (par.clusterworkflow[i]->uniqid == par.PARAM_NO_COMP_BIAS_CORR.uniqid && par.clusterworkflow[i]->wasSet) {
            compositionBiasSet = true;
        }
        if (par.clusterworkflow[i]->uniqid == par.PARAM_MIN_DIAG_SCORE.uniqid && par.clusterworkflow[i]->wasSet) {
            minDiagonalScoreSet = true;
        }
    }

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

    if (sensitivitySet == false) {
        par.sensitivity = setAutomaticThreshold(par.seqIdThr);
        Debug(Debug::INFO) << "Set cluster settings automatically to s=" << par.sensitivity << "\n";
    }

    const int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
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
        Debug(Debug::WARNING) << "combining cluster mode " << par.clusteringMode
                              << " in combination with coverage mode " << par.covMode << " can produce wrong results.\n"
                              << "Please use --cov-mode 2\n";
    }
    if (par.cascaded == true && par.clusteringMode == Parameters::CONNECTED_COMPONENT) {
        Debug(Debug::WARNING) << "connected component clustering produces less clusters in a single step clustering.\n"
                              << "Please use --single-step-cluster";
    }
    if (clusterStepsSet == false) {
        par.clusterSteps = setAutomaticIterations(par.sensitivity);
        Debug(Debug::INFO) << "Set cluster iterations to " << par.clusterSteps << "\n";
    }

    std::string hash = SSTR(par.hashParameter(par.filenames, par.clusterworkflow));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db3+"/latest");
    }
    std::string tmpDir = par.db3+"/"+hash;

    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    const int originalRescoreMode = par.rescoreMode;

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
    par.rescoreMode = originalRescoreMode;
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());

    if (par.cascaded) {
        // save some values to restore them later
        float targetSensitivity = par.sensitivity;
        int alphabetSize = par.alphabetSize;
        par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
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
        //cmd.addVariable("REASSIGN","TRUE");
        size_t olfMaxResSize = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        cmd.addVariable("ALIGNMENT_REASSIGN_PAR", par.createParameterString(par.align).c_str());
        par.maxResListLen = olfMaxResSize;
        // set parameter for first step
        FileUtil::writeFile(tmpDir + "/cascaded_clustering.sh", cascaded_clustering_sh, cascaded_clustering_sh_len);
        std::string program(tmpDir + "/cascaded_clustering.sh");
        cmd.execProgram(program.c_str(), par.filenames);
    } else {
        // same as above, clusthash needs a smaller alphabetsize
        size_t alphabetSize = par.alphabetSize;
        par.alphabetSize = Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE;
        cmd.addVariable("DETECTREDUNDANCY_PAR", par.createParameterString(par.clusthash).c_str());
        par.alphabetSize = alphabetSize;

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
