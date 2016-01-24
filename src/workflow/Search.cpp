#include <string>
#include <cassert>

#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"

int search(int argc, const char *argv[]) {
    std::string usage("\nCompares all sequences in the query database with all sequences in the target database.\n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de) & Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs search <queryDB> <targetDB> <outDB> <tmpDir> [opts]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.searchworkflow, 4);

    Debug::setDebugLevel(par.verbosity);

    CommandCaller cmd;
    if (par.numIterations > 1) {
        cmd.addVariable("NUM_IT", SSTR(par.numIterations));
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.addVariable("PROFILE_PAR", par.createParameterString(par.createprofiledb));
        cmd.execProgram(par.mmdir + "/bin/blastpgp.sh", argc, argv);
    } else {
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter));
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment));
        cmd.execProgram(par.mmdir + "/bin/blastp.sh", argc, argv);
    }

    // Should never get here
    assert(false);
    return 0;
}
