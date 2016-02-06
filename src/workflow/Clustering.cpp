#include "Parameters.h"

#include <string>
#include <cassert>

#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "WorkflowFunctions.h"
#include "FileUtil.h"


void setWorkflowDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.001;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

std::pair<float, bool> setAutomaticThreshold(float seqId) {
    float sens;
    bool cascaded = true;
    if (seqId <= 0.3) {
        sens = 7.0;
        cascaded = true;
    } else if (seqId > 0.7) {
        sens = 1.0;
    } else {
        sens = 1.0 + (1.5 * (0.7 - seqId) * 10);
    }
    if (sens <= 2.0) {
        cascaded = false;
    }
    return std::make_pair(sens, cascaded);
}

int clusteringworkflow(int argc, const char *argv[]) {
    std::string usage("\nCalculates the clustering of the sequences in the input database.\n");
    usage.append("Written by Martin Steinegger(martin.steinegger@mpibpc.mpg.de) Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs clusteringworkflow <sequenceDB> <outDB> <tmpDir> [opts]\n");


    Parameters par;
    setWorkflowDefaults(&par);
    par.parseParameters(argc, argv, usage, par.clusteringWorkflow, 3);

    Debug::setDebugLevel(par.verbosity);

    bool parameterSet = false;
    for (size_t i = 0; i < par.clusteringWorkflow.size(); i++) {
        if (par.clusteringWorkflow[i].uniqid == par.PARAM_S.uniqid && par.clusteringWorkflow[i].wasSet) {
            parameterSet = true;
        }
        if (par.clusteringWorkflow[i].uniqid == par.PARAM_CASCADED.uniqid && par.clusteringWorkflow[i].wasSet) {
            parameterSet = true;
        }
    }

    if (!par.noAutomaticThreshold && parameterSet == false) {
        std::pair<float, bool> settings = setAutomaticThreshold(par.seqIdThr);
        par.sensitivity = settings.first;
        par.cascaded = settings.second;
        Debug(Debug::WARNING) << "Set cluster settings automatic to s=" << par.sensitivity << " cascaded=" <<
        par.cascaded << "\n";
    }

    FileUtil::errorIfFileExist(par.db2.c_str());
    FileUtil::errorIfFileExist(par.db2Index.c_str());

    CommandCaller cmd;

    if(par.keepTempFiles) {
        cmd.addVariable("KEEP_TEMP", "TRUE");
    }

    if (par.cascaded) {
        float targetSensitivity = par.sensitivity;
        size_t maxResListLen = par.maxResListLen;

        // 1 is lowest sens
        par.clusteringMode = Parameters::GREEDY;
        par.sensitivity = 1;
        par.maxResListLen = 20;
        par.fragmentMerge = true;
        par.diagonalScoring = 0;
        par.compBiasCorrection = 0;
        cmd.addVariable("PREFILTER0_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT0_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("CLUSTER0_PAR", par.createParameterString(par.clustering).c_str());

        // set parameter for first step
        par.clusteringMode = Parameters::SET_COVER;
        par.sensitivity = targetSensitivity / 3.0;
        par.maxResListLen = 100;
        par.fragmentMerge = false;
        par.diagonalScoring = 1;
        par.compBiasCorrection = 1;
        cmd.addVariable("PREFILTER1_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT1_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("CLUSTER1_PAR", par.createParameterString(par.clustering).c_str());

        // set parameter for second step
        par.sensitivity = targetSensitivity * (2.0 / 3.0);
        par.maxResListLen = 200;
        cmd.addVariable("PREFILTER2_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT2_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("CLUSTER2_PAR", par.createParameterString(par.clustering).c_str());

        // set parameter for last step
        par.sensitivity = targetSensitivity;
        par.maxResListLen = maxResListLen;
        cmd.addVariable("PREFILTER3_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT3_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("CLUSTER3_PAR", par.createParameterString(par.clustering).c_str());

        std::string program(par.mmdir);
        program.append("/bin/cascaded_clustering.sh");
        cmd.execProgram(program.c_str(), 3, argv);

    } else {
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clustering).c_str());

        std::string program(par.mmdir);
        program.append("/bin/clustering.sh");
        cmd.execProgram(program.c_str(), 3, argv);
    }

    // Unreachable
    assert(false);

    return 0;
}
