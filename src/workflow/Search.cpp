#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>
#include <CommandCaller.h>
#include "WorkflowFunctions.h"
#include "Parameters.h"
#include "cluster2profile.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

int search (int argc, const char * argv[]) {

    std::string usage("\nCompares all sequences in the query database with all sequences in the target database.\n");
    usage.append(
            "Written by Martin Steinegger (martin.steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: search <queryDB> <targetDB> <outDB> <tmpDir> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.searchworkflow, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Debug::setDebugLevel(par.verbosity);
    if (par.numIterations > 1) {
        CommandCaller cmd;
        cmd.addVariable("NUM_IT", SSTR(par.numIterations));
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("PROFILE_PAR",   par.createParameterString(par.createprofiledb));
        cmd.callProgram(par.mmdir + "/bin/blastpgp.sh", argv , 4);
    } else {
        CommandCaller cmd;
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.callProgram(par.mmdir + "/bin/blastp.sh", argv, 4);
    }
    return 0;
}
