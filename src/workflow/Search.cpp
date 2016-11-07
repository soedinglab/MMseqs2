#include <string>
#include <cassert>
#include <FileUtil.h>
#include <blastpgp.sh.h>
#include <blastp.sh.h>
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
int search(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true, false, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);

    CommandCaller cmd;

    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    std::string templateDB(par.db2);
    cmd.addVariable("TARGET_DB_PREF", templateDB.c_str());

    if (par.numIterations > 1) {
        for (size_t i = 0; i < par.searchworkflow.size(); i++) {
            if (par.searchworkflow[i].uniqid == par.PARAM_E_PROFILE.uniqid && par.searchworkflow[i].wasSet== false) {
                par.evalProfile = 0.001;
            }
        }
        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        cmd.addVariable("PROFILE", SSTR((par.profile) ? 1 : 0).c_str());
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        cmd.addVariable("PROFILE_PAR",   par.createParameterString(par.result2profile).c_str());
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
        FileUtil::writeFile(par.db4 + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
        std::string program(par.db4 + "/blastpgp.sh");
        cmd.execProgram(program.c_str(), 4, argv);
    } else {

        bool startSens  = false;
        for (size_t i = 0; i < par.searchworkflow.size(); i++) {
            if (par.searchworkflow[i].uniqid == par.PARAM_START_SENS.uniqid && par.searchworkflow[i].wasSet) {
                startSens = true;
            }
        }

        if(startSens){
            cmd.addVariable("START_SENS", SSTR(par.startSens).c_str());
            cmd.addVariable("TARGET_SENS", SSTR((int)par.sensitivity).c_str());
            cmd.addVariable("SENS_STEP_SIZE", SSTR(par.sensStepSize).c_str());
        } else {
            cmd.addVariable("START_SENS", SSTR((int)par.sensitivity).c_str());
            cmd.addVariable("TARGET_SENS", SSTR((int)par.sensitivity).c_str());
            cmd.addVariable("SENS_STEP_SIZE", SSTR(1).c_str());
        }

        std::vector<MMseqsParameter> prefilterWithoutS;
        for(size_t i = 0; i < par.prefilter.size(); i++){
            if(par.prefilter[i].uniqid != par.PARAM_S.uniqid ){
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(prefilterWithoutS).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());


        FileUtil::writeFile(par.db4 + "/blastp.sh", blastp_sh, blastp_sh_len);
        std::string program(par.db4 + "/blastp.sh");
        cmd.execProgram(program.c_str(), 4, argv);
    }

    // Should never get here
    assert(false);
    return 0;
}
