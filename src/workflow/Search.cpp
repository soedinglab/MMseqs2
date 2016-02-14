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

    CommandCaller cmd;

    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());

    if (par.numIterations > 1) {
        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment).c_str());
        cmd.addVariable("PROFILE_PAR",   par.createParameterString(par.createprofiledb).c_str());
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.substractresult).c_str());
        std::string program(par.mmdir);
        program.append("/bin/blastpgp.sh");
        cmd.execProgram(program.c_str(), 4, argv);
    } else {
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.alignment).c_str());

        std::string program(par.mmdir);
        program.append("/bin/blastp.sh");
        cmd.execProgram(program.c_str(), 4, argv);
    }

    // Should never get here
    assert(false);
    return 0;
}
