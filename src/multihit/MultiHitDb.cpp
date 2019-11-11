#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "DBReader.h"
#include "CommandCaller.h"

#include "multihitdb.sh.h"

void setMultiHitDbWorkflowDefaults(Parameters *p) {
    p->orfMinLength = 30;
}

int multihitdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setMultiHitDbWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.filenames.back();
    par.filenames.pop_back();

    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << tmpDir << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << tmpDir << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(par.filenames, par.multihitdb));
    if(par.reuseLatest == true){
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest" );
    }
    tmpDir = tmpDir + "/" + hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    FileUtil::symlinkAlias(tmpDir, "latest");

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("OUTDB", outDb.c_str());
    cmd.addVariable("TMP_PATH", tmpDir.c_str());

    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());
    par.stat = "linecount";
    cmd.addVariable("RESULT2STATS_PAR", par.createParameterString(par.result2stats).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(tmpDir + "/multihitdb.sh", multihitdb_sh, multihitdb_sh_len);
    std::string program(tmpDir + "/multihitdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
