#include "Parameters.h"
#include "Util.h"
#include "DBReader.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

#include "createindex.sh.h"

#include <string>
#include <cassert>
#include <climits>

int createindex(Parameters &par, const Command &command, const std::string &indexerModule, const std::string &flag) {
    bool sensitivity = false;
    // only set kmerScore  to INT_MAX if -s was used
    for (size_t i = 0; i < par.createindex.size(); i++) {
        if (par.createindex[i]->uniqid == par.PARAM_S.uniqid && par.createindex[i]->wasSet) {
            par.kmerScore = INT_MAX;
            sensitivity=true;
            break;
        }
    }

    int dbType = FileUtil::parseDbType(par.db1.c_str());
    if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE) && sensitivity == false) {
        Debug(Debug::ERROR) << "Please adjust the sensitivity of your target profile index with -s.\n"
                               "Be aware that this searches can take huge amount of memory. \n";
        return EXIT_FAILURE;
    }

    std::string tmpDir = par.db2;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.createindex));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("INDEXER", indexerModule.c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    par.translate = 1;
    cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());
    cmd.addVariable("SPLIT_SEQ_PAR", par.createParameterString(par.splitsequence).c_str());
    if(indexerModule == "kmerindexdb"){
        cmd.addVariable("INDEX_PAR", par.createParameterString(par.kmerindexdb).c_str());
    }else{
        cmd.addVariable("INDEX_PAR", par.createParameterString(par.indexdb).c_str());
    }
    if(flag.size() > 0){
        cmd.addVariable(flag.c_str(), "1");
    }

    std::string program(tmpDir + "/createindex.sh");
    FileUtil::writeFile(program, createindex_sh, createindex_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}


int createlinindex(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.orfStartMode = 1;
    par.orfMinLength = 30;
    par.orfMaxLength = 32734;
    par.kmerScore = 0; // extract all k-mers
    par.maskMode = 0;
    par.spacedKmer = false;
    // VTML has a slightly lower sensitivity in the regression test
    par.seedScoringMatrixFile = MultiParam<char*>("blosum62.out", "nucleotide.out");

    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    int dbType = FileUtil::parseDbType(par.db1.c_str());
    bool isNucl = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);
    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES && par.PARAM_MAX_SEQ_LEN.wasSet == false){
        if(par.PARAM_MAX_SEQ_LEN.wasSet == false){
            par.maxSeqLen = 10000;
        }
    }
    par.printParameters(command.cmd, argc, argv, *command.params);

    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO){
        Debug(Debug::WARNING) << "Database " << par.db1 << " is a nucleotide database. \n"
                            << "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
        return EXIT_FAILURE;
    }
    return createindex(par, command, "kmerindexdb", (isNucl == false) ? "" : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED||
                                                                                           par.searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN) ? "TRANSLATED" : "LIN_NUCL");
}

int createindex(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.orfStartMode = 1;
    par.orfMinLength = 30;
    par.orfMaxLength = 32734;
    par.kmerScore = 0; // extract all k-mers
    par.sensitivity = 7.5;
    par.maskMode = 1;

    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SPLIT.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);

    int dbType = FileUtil::parseDbType(par.db1.c_str());
    bool isNucl = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);

    if (par.PARAM_STRAND.wasSet == false) {
        par.strand = 1;
    }
    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES ){
        if ( par.PARAM_K.wasSet == false) {
            par.kmerSize = 15;
        }
        if ( par.PARAM_MAX_SEQ_LEN.wasSet == false) {
            par.maxSeqLen = 10000;
        }

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
    }
    par.printParameters(command.cmd, argc, argv, *command.params);
    if(isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO){
        Debug(Debug::WARNING) << "Database " << par.db1 << " is a nucleotide database. \n"
                            << "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
        return EXIT_FAILURE;
    }
    return createindex(par, command, "indexdb",  (isNucl == false) ? "" : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED||
                                                                  par.searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN) ? "TRANSLATED" : "NUCL");
}
