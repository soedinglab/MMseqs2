#include <CommandCaller.h>
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "createtaxdb.sh.h"

int createtaxdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1, true, Parameters::PARSE_VARIADIC);

    std::string tmp = par.filenames.back();
    if (FileUtil::directoryExists(tmp.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << tmp << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmp.c_str()) == false){
            Debug(Debug::ERROR) << "Can not create tmp folder " << tmp << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << tmp << "\n";
        }
    }

    CommandCaller cmd;

    cmd.addVariable("TMP_PATH", tmp.c_str());

    if(par.filenames.size() == 4) {
        cmd.addVariable("DOWNLOAD_DATA", "0");

    }else if(par.filenames.size() == 2) {
        cmd.addVariable("DOWNLOAD_DATA", "1");
    }

    FileUtil::writeFile(tmp + "/createindex.sh", createtaxdb_sh, createtaxdb_sh_len);
    std::string program(tmp + "/createindex.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
