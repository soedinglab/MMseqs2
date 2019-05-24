#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "map.sh.h"

void setMapWorkflowDefaults(Parameters *p) {
    p->compBiasCorrection = 0;
    p->maskMode = 0;
    p->covThr = 0.95;
    p->covMode = 2;
    p->seqIdThr = 0.9;
    p->sensitivity = 2;
    p->rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    p->sortResults = true;
    //p->orfLongest = true;
    p->orfStartMode = 1;
    p->orfMinLength = 10;
    p->orfMaxLength = 32734;
}

int map(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMapWorkflowDefaults(&par);

    par.overrideParameterDescription((Command &) command, par.PARAM_OVERLAP.uniqid, NULL, NULL,
                                     par.PARAM_OVERLAP.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_DB_OUTPUT.uniqid, NULL, NULL,
                                     par.PARAM_DB_OUTPUT.category | MMseqsParameter::COMMAND_EXPERT);

    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.extractorfs[i]->uniqid, NULL, NULL, par.extractorfs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.overrideParameterDescription((Command &)command, par.translatenucs[i]->uniqid, NULL, NULL, par.translatenucs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, 4);

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(par.filenames, par.mapworkflow));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("RUNNER", par.runner.c_str());

    par.mapworkflow.push_back(&(par.PARAM_ALIGNMENT_MODE));
    par.alignmentMode = 4;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.mapworkflow).c_str());

    std::string program = tmpDir + "/map.sh";
    FileUtil::writeFile(program, map_sh, map_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
