#include "Parameters.h"
#include <string>
#include <cassert>
#include <climits>
#include <Util.h>
#include <DBReader.h>
#include "createindex.sh.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

int createindex(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.orfStartMode = 0;
    par.orfMinLength = 30;
    par.orfMaxLength = 98202; // 32734 AA (just to be sure)
    par.kmerScore = 0; // extract all k-mers
    par.sensitivity = 7.5;
    par.parseParameters(argc, argv, command, 2);
    bool sensitivity = false;
    // only set kmerScore  to INT_MAX if -s was used
    for (size_t i = 0; i < par.createindex.size(); i++) {
        if (par.createindex[i].uniqid == par.PARAM_S.uniqid && par.createindex[i].wasSet) {
            par.kmerScore = INT_MAX;
            sensitivity=true;
            break;
        }
    }

    int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    if(dbType==-1){
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    if (dbType == Sequence::HMM_PROFILE && sensitivity == false){
        Debug(Debug::ERROR) << "Please adjust the sensitivity of your target profile index with -s.\n"
                               "Be aware that this searches can take huge amount of memory. \n";
        EXIT(EXIT_FAILURE);
    }

    size_t hash = par.hashParameter(par.filenames, par.createindex);
    std::string tmpDir = par.db2 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    if (dbType == Sequence::NUCLEOTIDES) {
        cmd.addVariable("NUCL", "TRUE");
    }
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("INDEX_PAR", par.createParameterString(par.indexdb).c_str());

    FileUtil::writeFile(par.db2 + "/createindex.sh", createindex_sh, createindex_sh_len);
    std::string program(par.db2 + "/createindex.sh");
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}
