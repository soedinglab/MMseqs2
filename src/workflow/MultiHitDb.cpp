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

int multihitdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setMultiHitDbWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 2, true);


    if(FileUtil::directoryExists(par.db2.c_str())==false){
        Debug(Debug::INFO) << "Tmp " << par.db2 << " folder does not exist or is not a directory.\n";
        if(FileUtil::makeDir(par.db2.c_str()) == false){
            Debug(Debug::ERROR) << "Could not crate tmp folder " << par.db2 << ".\n";
            EXIT(EXIT_FAILURE);
        }else{
            Debug(Debug::INFO) << "Created dir " << par.db2 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.multihitdb);
    std::string tmpDir = par.db2+"/"+SSTR(hash);
    if(FileUtil::directoryExists(tmpDir.c_str())==false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    const int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }

    if (dbType == Sequence::NUCLEOTIDES) {
        cmd.addVariable("NUCL", "TRUE");
    }

    par.splitSeqByLen = false;
    par.clusterDB = true;
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());
    par.stat = "linecount";
    cmd.addVariable("RESULT2STATS_PAR", par.createParameterString(par.result2stats).c_str());

    FileUtil::writeFile(tmpDir + "/multihitdb.sh", multihitdb_sh, multihitdb_sh_len);
    std::string program(tmpDir + "/multihitdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
