#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "multihitsearch.sh.h"

void setMultiHitSearchWorkflowDefaults(Parameters *p) {
    p->sensitivity = 7;
    p->cov = 0.7;
    p->evalThr = 100;

    // TODO: Needs to be count of target genomes
    p->maxSequences = 1500;
    p->scoreBias = 0.3;
}

int summerizeresultsbyset(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setMultiHitSearchWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 4);


    if(FileUtil::directoryExists(par.db4.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db4.c_str()) == false){
            Debug(Debug::ERROR) << "Could not crate tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.multihitsearch);
    std::string tmpDir = par.db4+"/"+SSTR(hash);
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    par.splitSeqByLen = false;
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    par.simpleBestHit = 1;
    cmd.addVariable("AGGREGATE_PAR", par.createParameterString(par.besthitbyset).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(tmpDir + "/multihitsearch.sh", multihitsearch_sh, multihitsearch_sh_len);
    std::string program(tmpDir + "/multihitsearch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
