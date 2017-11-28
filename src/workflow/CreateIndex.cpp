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
    par.includeHeader = true;
    par.orfLongest = true;
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

    if(dbType==DBReader<unsigned int>::DBTYPE_PROFILE && sensitivity == false){
        Debug(Debug::ERROR) << "Please adjust the sensitivity of your target profile index with -s.\n"
                               "Be aware that this searches can take huge amount of memory. \n";
        EXIT(EXIT_FAILURE);
    }

    if(FileUtil::directoryExists(par.db2.c_str())==false){
        Debug(Debug::WARNING) << "Tmp " << par.db2 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db2.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db2 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::WARNING) << "Created dir " << par.db2 << "\n";
        }
    }

    CommandCaller cmd;

    if(dbType==DBReader<unsigned int>::DBTYPE_NUC){
        cmd.addVariable("NUCL", "TRUE");
    }
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("INDEX_PAR", par.createParameterString(par.indexdb).c_str());

    FileUtil::writeFile(par.db2 + "/indexdb.sh", createindex_sh, createindex_sh_len);
    std::string program(par.db2 + "/indexdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}
