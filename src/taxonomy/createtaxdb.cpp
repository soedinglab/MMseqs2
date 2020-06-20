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
    par.parseParameters(argc, argv, command, true, 0, 0);

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
    //[<i:taxMappingFile> <i:ncbi-taxdump-folder>]
    CommandCaller cmd;

    cmd.addVariable("TMP_PATH", tmp.c_str());
    if (par.taxMappingFile.empty()) {
        cmd.addVariable("DOWNLOAD_MAPPING", "1");
    }else{
        cmd.addVariable("DOWNLOAD_MAPPING", "0");
        cmd.addVariable("MAPPINGFILE", par.taxMappingFile.c_str());

    }
    if (par.ncbiTaxDump.empty()) {
        cmd.addVariable("DOWNLOAD_NCBITAXDUMP", "1");
    }else{
        cmd.addVariable("DOWNLOAD_NCBITAXDUMP", "0");
        cmd.addVariable("NCBITAXINFO", par.ncbiTaxDump.c_str());
    }
    cmd.addVariable("ARIA_NUM_CONN", SSTR(std::min(16, par.threads)).c_str());
    FileUtil::writeFile(tmp + "/createindex.sh", createtaxdb_sh, createtaxdb_sh_len);
    std::string program(tmp + "/createindex.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
