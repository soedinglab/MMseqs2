#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "multihitdb.sh.h"

void setMultiHitDbWorkflowDefaults(Parameters *p) {
    p->orfMinLength = 30;
}

int multihitdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setMultiHitDbWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 3, true, Parameters::PARSE_VARIADIC);

    if (par.filenames.size() < 3) {
        Debug(Debug::ERROR) << "Not enough databases passed.\n";
        return EXIT_FAILURE;
    }

    std::string tmpDir(par.filenames.back());
    par.filenames.pop_back();
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    par.splitSeqByLen = false;
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());

    FileUtil::writeFile(tmpDir + "/multihitdb.sh", multihitdb_sh, multihitdb_sh_len);
    std::string program(tmpDir + "/multihitdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
