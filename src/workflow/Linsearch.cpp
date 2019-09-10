#include "DBReader.h"
#include "CommandCaller.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include "LinsearchIndexReader.h"

#include "linsearch.sh.h"

namespace Linsearch {
#include "translated_search.sh.h"
}

#include <climits>
#include <cassert>

void setLinsearchDefaults(Parameters *p) {
    p->spacedKmer = false;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    p->maskMode = 0;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
    p->evalProfile = 0.1;

    // VTML has a slightly lower sensitivity in the regression test
    p->seedScoringMatrixFile = MultiParam<char*>("blosum62.out", "nucleotide.out");
}


int linsearch(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinsearchDefaults(&par);
    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);

    const int queryDbType = FileUtil::parseDbType(par.db1.c_str());
    std::string indexStr = LinsearchIndexReader::searchForIndex(par.db2);
    if (indexStr.empty()) {
        Debug(Debug::ERROR) << par.db2 << " needs to be index.\n";
        Debug(Debug::ERROR) << "createlinindex " << par.db2 << ".\n";
        EXIT(EXIT_FAILURE);
    }
    int targetDbType = 0;
    if(indexStr != ""){
        DBReader<unsigned int> dbr(indexStr.c_str(), (indexStr+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
        targetDbType = data.seqType;
        dbr.close();
    }

    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR)
                << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-Profile searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }

    const bool isNuclSearch = (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)
                               && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));
//    if(isNuclSearch == true){
//        setNuclSearchDefaults(&par);
//    }else{
//        par.overrideParameterDescription((Command &) command, par.PARAM_STRAND.uniqid, NULL, NULL,
//                                         par.PARAM_STRAND.category | MMseqsParameter::COMMAND_EXPERT);
//    }


    par.filenames[1] = indexStr;
    const bool isTranslatedNuclSearch =
            isNuclSearch == false && (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) ||
                                      Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));

    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) ||
                           Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE))) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases.\n";
        EXIT(EXIT_FAILURE);
    }

    par.printParameters(command.cmd, argc, argv, par.searchworkflow);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.linsearchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("FILTER", "1");
    int oldCovMode = par.covMode;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    if (par.PARAM_COV_MODE.wasSet == false) {
        par.covMode = Parameters::COV_MODE_TARGET;
    }
    float oldCov = par.covThr;
    par.covThr = std::max(par.covThr, 0.9f);
    cmd.addVariable("RESCORE_FILTER_PAR", par.createParameterString(par.rescorediagonal).c_str());
    par.covMode = oldCovMode;
    par.covThr = oldCov;

    cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
    cmd.addVariable("KMERSEARCH_PAR", par.createParameterString(par.kmersearch).c_str());
    float oldEval = par.evalThr;
    par.evalThr = 100000;
    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
    par.evalThr = oldEval;
    cmd.addVariable("SWAPRESULT_PAR", par.createParameterString(par.swapresult).c_str());
    cmd.addVariable("NUCL", isNuclSearch ? "1" : NULL);

    std::string program = tmpDir + "/linsearch.sh";
    FileUtil::writeFile(program, linsearch_sh, linsearch_sh_len);

    if (isTranslatedNuclSearch == true) {
        cmd.addVariable("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : NULL);
        cmd.addVariable("QUERY_NUCL", Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) ? "TRUE" : NULL);
        cmd.addVariable("TARGET_NUCL", Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES) ? "TRUE" : NULL);
        par.translate = 1;
        cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
        cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment).c_str());
        cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
        cmd.addVariable("SEARCH", program.c_str());
        program = std::string(tmpDir + "/translated_search.sh");
        FileUtil::writeFile(program, Linsearch::translated_search_sh, Linsearch::translated_search_sh_len);
    }
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
