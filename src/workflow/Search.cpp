#include "DBReader.h"
#include "CommandCaller.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include "searchtargetprofile.sh.h"
#include "searchslicedtargetprofile.sh.h"
#include "blastpgp.sh.h"
#include "translated_search.sh.h"
#include "blastp.sh.h"
#include "blastn.sh.h"
#include "Parameters.h"

#include <iomanip>
#include <climits>
#include <cassert>


void setSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
    p->evalProfile = 0.1;
}


int computeSearchMode(int queryDbType, int targetDbType, int targetSrcDbType, int searchType) {
    // reject unvalid search
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-Profile searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }
    // index was used
    if(targetSrcDbType == -1) {
        if(Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
           Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES))
        {
            if(searchType == Parameters::SEARCH_TYPE_AUTO){
                // WARNING because its not really an error, just a req. parameter
                Debug(Debug::WARNING) << "It is unclear from the input if a translated or nucleotide search should be performed\n"
                                         "Please provide the parameter --search-type 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)\n";
                EXIT(EXIT_FAILURE);
            }
            // nucl/nucl
            // nucl/nucl translated
            if(searchType == Parameters::SEARCH_TYPE_TRANSLATED||searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN){
                return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED| Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
            }else if (searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES ){
                return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE| Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
            } else {
                Debug(Debug::ERROR) << "--search-type 1 (amino acid) can not used in combination with a nucleotide database\n "
                                       "The only possible options --search-types 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)\n";
                EXIT(EXIT_FAILURE);
            }
        }
        // protein/protein
        if (Parameters::isEqualDbtype(queryDbType,  Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // profile/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // profile/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // nucleotide/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // nucleotide/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }
    } else{

        // protein/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // profile/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_HMM_PROFILE)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // profile/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

        // nucleotide/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_HMM_PROFILE)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // nucleotide/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            (Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS))) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            (Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES))) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

        // nucl/nucl
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE | Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
        }

        // nucl/nucl translated
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

    }
    Debug(Debug::ERROR) << "Invalid input database and --search-type combination\n"
                        << "queryDbType: " << Parameters::getDbTypeName(queryDbType) << "\n"
                        << "targetDbType: " <<  Parameters::getDbTypeName(targetDbType) << "\n"
                        << "targetSrcDbType: " <<  Parameters::getDbTypeName(targetSrcDbType) << "\n"
                        << "searchMode: " << searchType << "\n";
    EXIT(EXIT_FAILURE);
}



void setNuclSearchDefaults(Parameters *p) {
    // leave ungapped alignment untouched
    if(p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED){
        p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    //p->orfLongest = true;
    p->exactKmerMatching = true;
//    if ( p->PARAM_DIAGONAL_SCORING.wasSet == false) {
//        p->diagonalScoring = 0;
//    }
    if ( p->PARAM_STRAND.wasSet == false) {
        p->strand = 2;
    }
    if ( p->PARAM_K.wasSet == false) {
        p->kmerSize = 15;
    }
    if (  p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
}


int search(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    setSearchDefaults(&par);
    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);

    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);

    int targetDbType = FileUtil::parseDbType(par.db2.c_str());
    std::string targetDB =  (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
    int targetSrcDbType = -1;
    if(indexStr != "" || Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB)){
        indexStr = par.db2;
        DBReader<unsigned int> dbr(targetDB.c_str(), (targetDB+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
        targetSrcDbType = data.srcSeqType;
        targetDbType = data.seqType;
    }
    const int queryDbType = FileUtil::parseDbType(par.db1.c_str());
    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR)
                << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    int searchMode = computeSearchMode(queryDbType, targetDbType, targetSrcDbType, par.searchType);
    if ((searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE) && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE)) {
        setNuclSearchDefaults(&par);
    } else{
        par.PARAM_STRAND.addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    // FIXME: use larger default k-mer size in target-profile case if memory is available
    // overwrite default kmerSize for target-profile searches and parse parameters again
    if (par.exhaustiveSearch == false && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) && par.PARAM_K.wasSet == false) {
        par.kmerSize = 5;
    }

    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE ))) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases\n";
        EXIT(EXIT_FAILURE);
    }

    if (isUngappedMode && par.lcaSearch) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with lca search\n";
        EXIT(EXIT_FAILURE);
    }

    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
            par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Iterative target-profile searches are not supported.\n";
            EXIT(EXIT_FAILURE);
        }

        par.addBacktrace = true;
        if (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE) {
            for (size_t i = 0; i < par.searchworkflow.size(); i++) {
                if (par.searchworkflow[i]->uniqid == par.PARAM_REALIGN.uniqid && par.searchworkflow[i]->wasSet) {
                    par.printUsageMessage(command,
                                          MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
                    Debug(Debug::ERROR) << "Cannot realign query profiles.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            par.realign = false;
        }
    }
    par.printParameters(command.cmd, argc, argv, par.searchworkflow);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    const int originalRescoreMode = par.rescoreMode;
    CommandCaller cmd;
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("VERB_COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    if (isUngappedMode) {
        cmd.addVariable("ALIGN_MODULE", "rescorediagonal");
    } else if (par.lcaSearch) {
        cmd.addVariable("ALIGN_MODULE", "lcaalign");
    } else {
        cmd.addVariable("ALIGN_MODULE", "align");
    }
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    std::string program;
    cmd.addVariable("RUNNER", par.runner.c_str());
//    cmd.addVariable("ALIGNMENT_DB_EXT", Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ) ? ".255" : "");
    par.filenames[1] = targetDB;
    if (par.exhaustiveSearch == true) {
        // By default (0), diskSpaceLimit (in bytes) will be set in the workflow to use as much as possible
        cmd.addVariable("AVAIL_DISK", SSTR(static_cast<size_t>(par.diskSpaceLimit)).c_str());

        // correct Eval threshold for inverted search
        const size_t queryDbSize = FileUtil::countLines(par.db1Index.c_str());
        const size_t targetDbSize = FileUtil::countLines(par.db2Index.c_str());
        par.evalThr *= ((float) queryDbSize)/targetDbSize;

        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = std::max((size_t)300, queryDbSize);
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        par.maxResListLen = maxResListLen;
        double originalEvalThr = par.evalThr;
        par.evalThr = std::numeric_limits<double>::max();
        cmd.addVariable("SWAPRES_PAR", par.createParameterString(par.swapresult).c_str());
        par.evalThr = originalEvalThr;
        cmd.addVariable("FILTER_PAR", par.createParameterString(par.filterresult).c_str());
        if(par.exhaustiveFilterMsa == 1){
            cmd.addVariable("FILTER_RESULT", "1");
        }
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
            par.alignmentOutputMode = Parameters::ALIGNMENT_OUTPUT_CLUSTER;
            cmd.addVariable("ALIGNMENT_IT_PAR", par.createParameterString(par.align).c_str());
        }

        par.covMode = originalCovMode;

        program = tmpDir + "/searchslicedtargetprofile.sh";
        FileUtil::writeFile(program, searchslicedtargetprofile_sh, searchslicedtargetprofile_sh_len);
    } else if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
        // we need to align all hits in case of target Profile hits
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        }
        par.covMode = originalCovMode;
        par.maxResListLen = maxResListLen;
        cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
        FileUtil::writeFile(tmpDir + "/searchtargetprofile.sh", searchtargetprofile_sh, searchtargetprofile_sh_len);
        program = std::string(tmpDir + "/searchtargetprofile.sh");
    } else if (par.numIterations > 1) {
        cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
        cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
        cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

        double originalEval = par.evalThr;
        par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr  : par.evalProfile;
        for (int i = 0; i < par.numIterations; i++) {
            if (i == 0 && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) == false) {
                par.realign = true;
            }

            if (i > 0) {
//                par.queryProfile = true;
                par.realign = false;
            }

            if (i == (par.numIterations - 1)) {
                par.evalThr = originalEval;
            }

            cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.prefilter).c_str());
            if (isUngappedMode) {
                par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
                cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                                par.createParameterString(par.rescorediagonal).c_str());
                par.rescoreMode = originalRescoreMode;
            } else {
                cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                                par.createParameterString(par.align).c_str());
            }
            par.pca = 0.0;
            cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),
                            par.createParameterString(par.result2profile).c_str());
            par.pca = 1.0;
        }

        FileUtil::writeFile(tmpDir + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
        program = std::string(tmpDir + "/blastpgp.sh");
    } else {
        if (par.sensSteps > 1) {
            if (par.startSens > par.sensitivity) {
                Debug(Debug::ERROR) << "--start-sens should not be greater -s.\n";
                EXIT(EXIT_FAILURE);
            }
            cmd.addVariable("SENSE_0", SSTR(par.startSens).c_str());
            float sensStepSize = (par.sensitivity - par.startSens) / (static_cast<float>(par.sensSteps) - 1);
            for (int step = 1; step < par.sensSteps; step++) {
                std::string stepKey = "SENSE_" + SSTR(step);
                float stepSense = par.startSens + sensStepSize * step;
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << stepSense;
                std::string value = stream.str();
                cmd.addVariable(stepKey.c_str(), value.c_str());
            }
            cmd.addVariable("STEPS", SSTR((int) par.sensSteps).c_str());
        } else {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << par.sensitivity;
            std::string sens = stream.str();
            cmd.addVariable("SENSE_0", sens.c_str());
            cmd.addVariable("STEPS", SSTR(1).c_str());
        }

        std::vector<MMseqsParameter*> prefilterWithoutS;
        for (size_t i = 0; i < par.prefilter.size(); i++) {
            if (par.prefilter[i]->uniqid != par.PARAM_S.uniqid) {
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(prefilterWithoutS).c_str());
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        }
        FileUtil::writeFile(tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
        program = std::string(tmpDir + "/blastp.sh");
    }

    if (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED|Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)) {
        cmd.addVariable("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : NULL);
        FileUtil::writeFile(tmpDir + "/translated_search.sh", translated_search_sh, translated_search_sh_len);
        cmd.addVariable("QUERY_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED) ? "TRUE" : NULL);
        cmd.addVariable("TARGET_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)  ? "TRUE" : NULL);
        cmd.addVariable("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
        par.subDbMode = 1;
        cmd.addVariable("CREATESUBDB_PAR", par.createParameterString(par.createsubdb).c_str());
        par.translate = 1;
        cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
        cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment).c_str());
        cmd.addVariable("SEARCH", program.c_str());
        program = std::string(tmpDir + "/translated_search.sh");
    }else if(searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE &&
            searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE){
        FileUtil::writeFile(tmpDir + "/blastn.sh", blastn_sh, blastn_sh_len);
        //  0: reverse, 1: forward, 2: both
        switch (par.strand){
            case 0:
                par.forwardFrames= "";
                par.reverseFrames= "1";
                cmd.addVariable("EXTRACTFRAMES","TRUE");
                break;
            case 1:
                par.forwardFrames= "1";
                par.reverseFrames= "";
                break;
            case 2:
                par.forwardFrames= "1";
                par.reverseFrames= "1";
                cmd.addVariable("EXTRACTFRAMES","TRUE");
                break;
        }
        cmd.addVariable("SPLITSEQUENCE_PAR", par.createParameterString(par.splitsequence).c_str());
        if(indexStr=="") {
            cmd.addVariable("NEEDTARGETSPLIT", "TRUE");
        }
        cmd.addVariable("NEEDQUERYSPLIT","TRUE");
        cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());
        cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment).c_str());
        cmd.addVariable("SEARCH", program.c_str());
        program = std::string(tmpDir + "/blastn.sh");

    }
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);

    return 0;
}
