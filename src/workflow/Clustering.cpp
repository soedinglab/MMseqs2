#include "Parameters.h"

#include <string>
#include <iostream>
#include <list>
#include <CommandCaller.h>

#include "WorkflowFunctions.h"


void setWorkflowDefaults(Parameters* p) {
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.001;
}

std::pair<float, bool> setAutomaticThreshold(float seqId){
    float sens;
    bool cascaded = true;
    if(seqId <= 0.3){
        sens = 7.0;
        cascaded = true;
    } else if (seqId > 0.7){
        sens = 1.0;
    } else {
        sens = 1.0 + (1.5 * (0.7 - seqId) * 10);
    }
    if(sens <= 2.0){
        cascaded = false;
    }
    return std::make_pair(sens, cascaded);
}

int clusteringworkflow (int argc, const char * argv[]) {

    std::string usage("\nCalculates the clustering of the sequences in the input database.\n");
    usage.append("Written by Martin Steinegger(martin.steinegger@mpibpc.mpg.de) Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_cluster <sequenceDB> <outDB> <tmpDir> [opts]\n");
    //            "--restart          \t[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
    //            "                \t     \tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering.\n"
    //            "CASCADED CLUSTERING OPTIONS:\n"
    /*            "--step          \t[int]\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging.\n"
     "\nRESTART OPTIONS NOTE:\n"
     "                \t     \tIt is assumed that all valid intermediate results exist.\n"
     "                \t     \tValid intermediate results are taken from the tmp directory specified by the user.\n"
     "                \t     \tFor the cascaded clustering, --restart and --step options can be combined.\n"*/

    Parameters par;
    setWorkflowDefaults(&par);
    par.parseParameters(argc, argv, usage, par.clusteringWorkflow, 3);
    Debug::setDebugLevel(par.verbosity);
    bool parameterSet = false;
    for(size_t i = 0; i < par.clusteringWorkflow.size(); i++) {
        if (par.clusteringWorkflow[i].uniqid == par.PARAM_S.uniqid && par.clusteringWorkflow[i].wasSet) {
            parameterSet = true;
        }
        if (par.clusteringWorkflow[i].uniqid == par.PARAM_CASCADED.uniqid && par.clusteringWorkflow[i].wasSet){
            parameterSet = true;
        }
    }
    if(parameterSet == false){
        std::pair<float, bool> settings = setAutomaticThreshold(par.seqIdThr);
        par.sensitivity = settings.first;
        par.cascaded = settings.second;
        Debug(Debug::WARNING) << "Set cluster settings automatic to s=" << par.sensitivity << " cascaded="<<par.cascaded << "\n";
    }
    DBWriter::errorIfFileExist(par.db2.c_str());

    if (par.cascaded) {
        CommandCaller cmd;

        float targetSensitivity = par.sensitivity;
        size_t maxResListLen = par.maxResListLen;

// set parameter for first step

        par.sensitivity   = 1; // 1 is lowest sens

        par.zscoreThr     = getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = 100;
        cmd.addVariable("PREFILTER1_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT1_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER1_PAR", par.createParameterString(par.clustering));

// set parameter for second step
        par.sensitivity = targetSensitivity / 2.0;
        par.zscoreThr =  getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = 200;
        cmd.addVariable("PREFILTER2_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT2_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER2_PAR",   par.createParameterString(par.clustering));

// set parameter for last step
        par.sensitivity =  targetSensitivity;
        par.zscoreThr = getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = maxResListLen;
        cmd.addVariable("PREFILTER3_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT3_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER3_PAR", par.createParameterString(par.clustering));

        cmd.callProgram(par.mmdir + "/bin/cascaded_clustering.sh",argv, 3);

    }   else{
        CommandCaller cmd;
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER_PAR",   par.createParameterString(par.clustering));
        cmd.callProgram(par.mmdir + "/bin/clustering.sh",argv, 3);
    }

    return 0;
}
