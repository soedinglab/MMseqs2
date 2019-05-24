#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "summarizeresultsbyset.sh.h"

int summerizeresultsbyset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 5);

    if (FileUtil::directoryExists(par.db5.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db5 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db5.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db5 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db5 << "\n";
        }
    }

    size_t hash = par.hashParameter(par.filenames, par.summerizeresultsbyset);
    std::string tmpDir = par.db5 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
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
    cmd.addVariable("RESULTSBYSET_PAR", par.createParameterString(par.summerizeresultsbyset).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(tmpDir + "/summerizeresultsbyset.sh", summarizeresultsbyset_sh, summarizeresultsbyset_sh_len);
    std::string program(tmpDir + "/summerizeresultsbyset.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
