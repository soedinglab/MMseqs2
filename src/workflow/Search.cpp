#include <string>
#include <cassert>
#include <FileUtil.h>
#include <searchtargetprofile.sh.h>
#include <blastpgp.sh.h>
#include <blastp.sh.h>
#include <iomanip>
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"


void setSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
}


int search(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setSearchDefaults(&par);
    par.overrideParameterDescription((Command &)command, par.PARAM_COV_MODE.uniqid, NULL, NULL, par.PARAM_COV_MODE.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_C.uniqid, NULL, NULL, par.PARAM_C.category |MMseqsParameter::COMMAND_EXPERT );
    par.overrideParameterDescription((Command &)command, par.PARAM_MIN_SEQ_ID.uniqid, NULL, NULL, par.PARAM_MIN_SEQ_ID.category |MMseqsParameter::COMMAND_EXPERT );

    par.parseParameters(argc, argv, command, 4, false, false, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        if (par.targetProfile == true) {
            par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Iterative target-profile searches are not supported.\n";
            EXIT(EXIT_FAILURE);
        }

        par.addBacktrace = true;
        if (par.queryProfile == true) {
            for (size_t i = 0; i < par.searchworkflow.size(); i++) {
                if (par.searchworkflow[i].uniqid == par.PARAM_REALIGN.uniqid && par.searchworkflow[i].wasSet) {
                    par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
                    Debug(Debug::ERROR) << "Cannot realign query profiles.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            par.realign = false;
        }
    }
    par.printParameters(argc, argv, par.searchworkflow);
    if(FileUtil::directoryExists(par.db4.c_str())==false){
        Debug(Debug::WARNING) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::WARNING) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.searchworkflow);
    std::string tmpDir = par.db4+"/"+SSTR(hash);
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    if(FileUtil::symlinkCreateOrRepleace(par.db4+"/latest", tmpDir) == false){
        Debug(Debug::WARNING) << "Could not link latest folder in tmp." << tmpDir << ".\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    std::string templateDB(par.db2);
    cmd.addVariable("TARGET_DB_PREF", templateDB.c_str());
    if (par.targetProfile == true){
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        par.targetProfile = false;
        par.queryProfile = true;
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
        FileUtil::writeFile(tmpDir + "/searchtargetprofile.sh", searchtargetprofile_sh, searchtargetprofile_sh_len);
        std::string program(tmpDir + "/searchtargetprofile.sh");
        cmd.execProgram(program.c_str(), par.filenames);
    } else if (par.numIterations > 1) {
        for (size_t i = 0; i < par.searchworkflow.size(); i++) {
            if (par.searchworkflow[i].uniqid == par.PARAM_E_PROFILE.uniqid && par.searchworkflow[i].wasSet== false) {
                par.evalProfile = 0.1;
            }
        }
        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        cmd.addVariable("PROFILE", SSTR((par.queryProfile) ? 1 : 0).c_str());
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());

        float originalEval = par.evalThr;
        par.evalThr = par.evalProfile;
        for (int i = 0; i < par.numIterations; i++){
            if (i == 0 && par.queryProfile == false) {
                par.realign = true;
            }

            if (i > 0) {
                par.queryProfile = true;
                par.realign = false;
            }

            if (i == (par.numIterations - 1)) {
                par.evalThr = originalEval;
            }

            cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.prefilter).c_str());
            cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.align).c_str());
            cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),   par.createParameterString(par.result2profile).c_str());
        }

        FileUtil::writeFile(tmpDir + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
        std::string program(tmpDir + "/blastpgp.sh");
        cmd.execProgram(program.c_str(), par.filenames);
    } else {
        if(par.sensSteps > 1){
            if (par.startSens > par.sensitivity) {
                Debug(Debug::ERROR) << "--start-sens should not be greater -s.\n";
                EXIT(EXIT_FAILURE);
            }
            cmd.addVariable("SENSE_0", SSTR(par.startSens).c_str());
            float sensStepSize = (par.sensitivity - par.startSens)/ (static_cast<float>(par.sensSteps)-1);
            for(size_t step = 1; step < par.sensSteps; step++){
                std::string stepKey = "SENSE_" + SSTR(step);
                float stepSense =  par.startSens + sensStepSize * step;
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << stepSense;
                std::string value = stream.str();
                cmd.addVariable(stepKey.c_str(), value.c_str());
            }
            cmd.addVariable("STEPS", SSTR((int)par.sensSteps).c_str());
        } else {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << par.sensitivity;
            std::string sens = stream.str();
            cmd.addVariable("SENSE_0", sens.c_str());
            cmd.addVariable("STEPS", SSTR(1).c_str());
            //cmd.addVariable("SENS_STEP_SIZE", SSTR(1).c_str());
        }

        std::vector<MMseqsParameter> prefilterWithoutS;
        for(size_t i = 0; i < par.prefilter.size(); i++){
            if(par.prefilter[i].uniqid != par.PARAM_S.uniqid ){
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(prefilterWithoutS).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        FileUtil::writeFile(tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
        std::string program(tmpDir + "/blastp.sh");
        cmd.execProgram(program.c_str(), par.filenames);
    }

    // Should never get here
    assert(false);
    return 0;
}
