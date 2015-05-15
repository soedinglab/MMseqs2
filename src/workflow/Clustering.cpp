#include "Parameters.h"

#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>
#include <CommandCaller.h>

#include "WorkflowFunctions.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "Util.h"
#include "mergecluster.h"


void setWorkflowDefaults(Parameters* p) {
    p->searchMode = 1;
    p->spacedKmer = true;
    p->covThr = 0.8;
    p->evalThr = 0.001;
}

int clusteringworkflow (int argc, const char * argv[]) {

    std::string usage("\nCalculates the clustering of the sequences in the input database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
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
    std::vector<MMseqsParameter> params = par.combineList(par.prefilter, par.alignment);
    params.push_back(par.PARAM_CASCADED);
    params = par.combineList(params, par.clustering);

    par.parseParameters(argc, argv, usage, params, 3);
    Debug::setDebugLevel(par.verbosity);
    
    int restart = par.restart;
    int step = par.step;
    DBWriter::errorIfFileExist(par.db2.c_str());
    if (par.cascaded) {
        CommandCaller cmd;

        float targetSensitivity = par.sensitivity;
        size_t maxResListLen = par.maxResListLen;

        // set parameter for first step

        par.sensitivity   = 1; // 1 is lowest sens
        par.fastMode      = true; // be as fast a possible
        par.zscoreThr     = getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = 100;
        cmd.addVariable("PREFILTER1_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT1_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER1_PAR", par.createParameterString(par.clustering));
        par.fastMode      = false; // be as fast a possible

        // set parameter for second step
        par.sensitivity = targetSensitivity / 2.0;
        par.zscoreThr =  getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = 200;
        cmd.addVariable("PREFILTER2_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT2_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER2_PAR",   par.createParameterString(par.clustering));

        // set parameter for last step
        par.sensitivity = targetSensitivity;
        par.zscoreThr = getZscoreForSensitivity( par.sensitivity );
        par.maxResListLen = maxResListLen;
        cmd.addVariable("PREFILTER3_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT3_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER3_PAR", par.createParameterString(par.clustering));

        cmd.callProgram(par.mmdir + "/bin/cascaded_clustering.sh",(argv+1), 3);

    }   else{
        CommandCaller cmd;
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("CLUSTER_PAR",   par.createParameterString(par.clustering));
        cmd.callProgram(par.mmdir + "/bin/clustering.sh",(argv+1), 3);
    }

    return 0;
}