#include <string>
#include <cassert>
#include <FileUtil.h>
#include <searchtargetprofile.sh.h>
#include <blastpgp.sh.h>
#include <translated_search.sh.h>
#include <blastp.sh.h>
#include <iomanip>
#include <DBReader.h>
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"


void setSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    p->includeHeader = true;
    p->orfLongest = true;
    p->orfMinLength = 30;
    p->orfMaxLength = 98202; // 32734 AA (just to be sure)
}


int search(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setSearchDefaults(&par);
    par.overrideParameterDescription((Command &)command, par.PARAM_COV_MODE.uniqid, NULL, NULL, par.PARAM_COV_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_C.uniqid, NULL, NULL, par.PARAM_C.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_MIN_SEQ_ID.uniqid, NULL, NULL, par.PARAM_MIN_SEQ_ID.category | MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.extractorfs[i].uniqid, NULL, NULL, par.extractorfs[i].category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.translatenucs[i].uniqid, NULL, NULL, par.translatenucs[i].category | MMseqsParameter::COMMAND_EXPERT);
    }

    par.parseParameters(argc, argv, command, 4, false, 0, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);

    const int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    const int targetDbType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if (queryDbType == -1 || targetDbType == -1){
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    if (queryDbType == DBReader<unsigned int>::DBTYPE_PROFILE && targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE ){
        Debug(Debug::ERROR) << "Profile-Profile searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }

    if (queryDbType == DBReader<unsigned int>::DBTYPE_NUC && targetDbType == DBReader<unsigned int>::DBTYPE_NUC) {
        Debug(Debug::ERROR) << "Nucleotide-Nucleotide searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }

    const bool isTranslatedNuclSearch =
               (queryDbType == DBReader<unsigned int>::DBTYPE_NUC || targetDbType == DBReader<unsigned int>::DBTYPE_NUC);

    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        if (targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE) {
            par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Iterative target-profile searches are not supported.\n";
            EXIT(EXIT_FAILURE);
        }

        par.addBacktrace = true;
        if (queryDbType == DBReader<unsigned int>::DBTYPE_PROFILE) {
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
    if (FileUtil::directoryExists(par.db4.c_str())==false){
        Debug(Debug::WARNING) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::WARNING) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.searchworkflow);
    std::string tmpDir = par.db4+"/"+SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    if (FileUtil::symlinkCreateOrRepleace(par.db4+"/latest", tmpDir) == false){
        Debug(Debug::WARNING) << "Could not link latest folder in tmp." << tmpDir << ".\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    std::string program;
    cmd.addVariable("RUNNER", par.runner.c_str());
    if (targetDbType == DBReader<unsigned int>::DBTYPE_PROFILE){
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
        FileUtil::writeFile(tmpDir + "/searchtargetprofile.sh", searchtargetprofile_sh, searchtargetprofile_sh_len);
        program=std::string(tmpDir + "/searchtargetprofile.sh");
    } else if (par.numIterations > 1) {
        for (size_t i = 0; i < par.searchworkflow.size(); i++) {
            if (par.searchworkflow[i].uniqid == par.PARAM_E_PROFILE.uniqid && par.searchworkflow[i].wasSet== false) {
                par.evalProfile = 0.1;
            }
        }
        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        cmd.addVariable("PROFILE", SSTR((queryDbType == DBReader<unsigned int>::DBTYPE_PROFILE) ? 1 : 0).c_str());
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());

        float originalEval = par.evalThr;
        par.evalThr = par.evalProfile;
        for (int i = 0; i < par.numIterations; i++){
            if (i == 0 && queryDbType != DBReader<unsigned int>::DBTYPE_PROFILE) {
                par.realign = true;
            }

            if (i > 0) {
//                par.queryProfile = true;
                par.realign = false;
            }

            if (i == (par.numIterations - 1)) {
                par.evalThr = originalEval;
            }

            cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.prefilter).c_str());
            cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(), par.createParameterString(par.align).c_str());
            par.pca = 0.0;
            cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),   par.createParameterString(par.result2profile).c_str());
            par.pca = 1.0;

        }

        FileUtil::writeFile(tmpDir + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
        program = std::string(tmpDir + "/blastpgp.sh");
    } else {
        if (par.sensSteps > 1){
            if (par.startSens > par.sensitivity) {
                Debug(Debug::ERROR) << "--start-sens should not be greater -s.\n";
                EXIT(EXIT_FAILURE);
            }
            cmd.addVariable("SENSE_0", SSTR(par.startSens).c_str());
            float sensStepSize = (par.sensitivity - par.startSens)/ (static_cast<float>(par.sensSteps)-1);
            for (int step = 1; step < par.sensSteps; step++){
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
        }

        std::vector<MMseqsParameter> prefilterWithoutS;
        for (size_t i = 0; i < par.prefilter.size(); i++){
            if (par.prefilter[i].uniqid != par.PARAM_S.uniqid ){
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(prefilterWithoutS).c_str());
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        FileUtil::writeFile(tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
        program = std::string(tmpDir + "/blastp.sh");
    }

    if (isTranslatedNuclSearch==true){
        FileUtil::writeFile(tmpDir + "/translated_search.sh", translated_search_sh, translated_search_sh_len);
        if (queryDbType == DBReader<unsigned int>::DBTYPE_NUC) {
            cmd.addVariable("QUERY_NUCL", "TRUE");
        }
        if (targetDbType == DBReader<unsigned int>::DBTYPE_NUC) {
            cmd.addVariable("TARGET_NUCL", "TRUE");
        }
        cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
        cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
        cmd.addVariable("SEARCH", program.c_str());
        program = std::string(tmpDir + "/translated_search.sh");
    }
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
