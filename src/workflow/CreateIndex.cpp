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


int createindex(Parameters &par, std::string indexerModule, std::string flag) {
    bool sensitivity = false;
    // only set kmerScore  to INT_MAX if -s was used
    for (size_t i = 0; i < par.createindex.size(); i++) {
        if (par.createindex[i]->uniqid == par.PARAM_S.uniqid && par.createindex[i]->wasSet) {
            par.kmerScore = INT_MAX;
            sensitivity=true;
            break;
        }
    }

    int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    if (dbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }

    if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE) && sensitivity == false) {
        Debug(Debug::ERROR) << "Please adjust the sensitivity of your target profile index with -s.\n"
                               "Be aware that this searches can take huge amount of memory. \n";
        return EXIT_FAILURE;
    }

    if (FileUtil::directoryExists(par.db2.c_str())==false) {
        Debug(Debug::INFO) << "Tmp " << par.db2 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db2.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db2 << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db2 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.createindex);
    std::string tmpDir = par.db2 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    cmd.addVariable("INDEXER", indexerModule.c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());
    cmd.addVariable("SPLIT_SEQ_PAR", par.createParameterString(par.splitsequence).c_str());
    cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
    if(indexerModule == "kmerindexdb"){
        cmd.addVariable("INDEX_PAR", par.createParameterString(par.kmerindexdb).c_str());
    }else{
        cmd.addVariable("INDEX_PAR", par.createParameterString(par.indexdb).c_str());
    }
    if(flag.size() > 0){
        cmd.addVariable(flag.c_str(), "1");
    }

    FileUtil::writeFile(tmpDir + "/createindex.sh", createindex_sh, createindex_sh_len);
    std::string program(tmpDir + "/createindex.sh");
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}


int createlinindex(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.orfStartMode = 1;
    par.orfMinLength = 30;
    par.orfMaxLength = 98202; // 32734 AA (just to be sure)
    par.kmerScore = 0; // extract all k-mers
    par.maskMode = 0;

    par.parseParameters(argc, argv, command, 2, false);
    int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    bool isNucl = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);
    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES && par.PARAM_MAX_SEQ_LEN.wasSet == false){
        if(par.PARAM_MAX_SEQ_LEN.wasSet == false){
            par.maxSeqLen = 10000;
        }
    }
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);

    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO){
        Debug(Debug::ERROR) << "Database " << par.db1 << " is a nucleotide database. \n"
                            << "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
        return EXIT_FAILURE;
    }
    return createindex(par, "kmerindexdb", (isNucl == false) ? "" : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED) ? "TRANSLATED" : "LIN_NUCL");
}

int createindex(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.orfStartMode = 1;
    par.orfMinLength = 30;
    par.orfMaxLength = 98202; // 32734 AA (just to be sure)
    par.kmerScore = 0; // extract all k-mers
    par.sensitivity = 7.5;
    par.maskMode = 1;
    par.strand = 1;
    par.overrideParameterDescription((Command &) command, par.PARAM_MASK_RESIDUES.uniqid, "0: w/o low complexity masking, 1: with low complexity masking, 2: add both masked and unmasked sequences to index", "^[0-2]{1}", par.PARAM_MASK_RESIDUES.category);
    par.parseParameters(argc, argv, command, 2);
    //  0: reverse, 1: forward, 2: both
    switch (par.strand){
        case 0:
            par.forwardFrames= "";
            par.reverseFrames= "1";
            break;
        case 1:
            par.forwardFrames= "1";
            par.reverseFrames= "";
            break;
        case 2:
            par.forwardFrames= "1";
            par.reverseFrames= "1";
            break;
    }
    int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    bool isNucl = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);

    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES ){
        if ( par.PARAM_K.wasSet == false) {
            par.kmerSize = 15;
        }
        if ( par.PARAM_MAX_SEQ_LEN.wasSet == false) {
            par.maxSeqLen = 10000;
        }
    }

    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO){
        Debug(Debug::ERROR) << "Database " << par.db1 << " is a nucleotide database. \n"
                            << "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
        return EXIT_FAILURE;
    }
    return createindex(par, "indexdb",  (isNucl == false) ? "" : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED) ? "TRANSLATED" : "NUCL");
}
