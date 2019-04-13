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

#include <iomanip>
#include <climits>
#include <cassert>


void setSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    //p->orfLongest = true;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
    p->evalProfile = 0.1;
}



void setNuclSearchDefaults(Parameters *p) {
    // leave ungapped alignment untouched
    if(p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED){
        p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    //p->orfLongest = true;
    p->exactKmerMatching = true;

    if ( p->PARAM_STRAND.wasSet == false) {
        p->strand = 2;
    }
    if ( p->PARAM_K.wasSet == false) {
        p->kmerSize = 15;
    }
    if (  p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
    if( p->PARAM_GAP_OPEN.wasSet == false){
        p->gapOpen = 5;
    }
    if( p->PARAM_GAP_EXTEND.wasSet  == false){
        p->gapExtend = 2;
    }
}


int search(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    setSearchDefaults(&par);
    par.overrideParameterDescription((Command &) command, par.PARAM_COV_MODE.uniqid, NULL, NULL,
                                     par.PARAM_COV_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_C.uniqid, NULL, NULL,
                                     par.PARAM_C.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_MIN_SEQ_ID.uniqid, NULL, NULL,
                                     par.PARAM_MIN_SEQ_ID.category | MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++) {
        par.overrideParameterDescription((Command &) command, par.extractorfs[i]->uniqid, NULL, NULL,
                                         par.extractorfs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++) {
        par.overrideParameterDescription((Command &) command, par.translatenucs[i]->uniqid, NULL, NULL,
                                         par.translatenucs[i]->category | MMseqsParameter::COMMAND_EXPERT);
    }
    par.overrideParameterDescription((Command &) command, par.PARAM_THREADS.uniqid, NULL, NULL,
                                     par.PARAM_THREADS.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &) command, par.PARAM_V.uniqid, NULL, NULL,
                                     par.PARAM_V.category & ~MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, 4, false, 0,
                        MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);

    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);
    std::string targetDB =  (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
    int targetDbType;
    if(indexStr != ""){
        DBReader<unsigned int> dbr(targetDB.c_str(), (targetDB+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
        targetDbType = data.seqType;
    }else{
        targetDbType = DBReader<unsigned int>::parseDbType(targetDB.c_str());
    }

    const int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR)
                << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType,Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-Profile searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }

    bool isNuclNuclSearch = (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)
                            && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));
    if(isNuclNuclSearch == true && par.searchType == Parameters::SEARCH_TYPE_AUTO) {
        setNuclSearchDefaults(&par);
    } else if(isNuclNuclSearch == true && par.searchType == Parameters::SEARCH_TYPE_TRANSLATED) {
        isNuclNuclSearch = false;
    } else{
        par.overrideParameterDescription((Command &) command, par.PARAM_STRAND.uniqid, NULL, NULL,
                                         par.PARAM_STRAND.category | MMseqsParameter::COMMAND_EXPERT);
    }
    // FIXME: use larger default k-mer size in target-profile case if memory is available
    // overwrite default kmerSize for target-profile searches and parse parameters again
    if (par.sliceSearch == false && Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) && par.PARAM_K.wasSet == false) {
        par.kmerSize = 5;
    }

    const bool isTranslatedNuclSearch =
            isNuclNuclSearch==false && (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) ||
                                    Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));

    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(targetDbType,Parameters::DBTYPE_HMM_PROFILE))) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases.\n";
        EXIT(EXIT_FAILURE);
    }

    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
            par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Iterative target-profile searches are not supported.\n";
            EXIT(EXIT_FAILURE);
        }

        par.addBacktrace = true;
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE)) {
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

    if(isNuclNuclSearch == true && par.searchType == Parameters::SEARCH_TYPE_AUTO) {
        Debug(Debug::WARNING) << "Perform nucleotide search. \n";
        Debug(Debug::WARNING) << "To perform a translated search use --search-type 2 \n";
    }
        if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }

    std::string hash = SSTR(par.hashParameter(par.filenames, par.searchworkflow));
    if(par.reuseLatest){
        hash = FileUtil::getHashFromSymLink(par.db4+"/latest");
    }
    std::string tmpDir = par.db4+"/"+hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    const int originalRescoreMode = par.rescoreMode;

    CommandCaller cmd;
    cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    std::string program;
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("ALIGNMENT_DB_EXT", Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ) ? ".255" : "");
    par.filenames[1] = targetDB;
    if (par.sliceSearch == true) {
        if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) == false) {
            par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Sliced search only works with profiles as targets.\n";
            EXIT(EXIT_FAILURE);
        }
        cmd.addVariable("MAX_STEPS", SSTR(30).c_str());
        cmd.addVariable("MAX_RESULTS_PER_QUERY", SSTR(par.maxResListLen).c_str());
        size_t diskLimit;

        // By default (0), diskSpaceLimit will be set in the workflow to use as much as possible
        diskLimit = static_cast<size_t>(par.diskSpaceLimit) * 1024; // in kb
        cmd.addVariable("AVAIL_DISK", SSTR(static_cast<size_t>(diskLimit)).c_str());

        // --max-seqs and --offset-results are set inside the workflow
        std::vector<MMseqsParameter*> prefilter;
        for (size_t i = 0; i < par.prefilter.size(); i++){
            if (par.prefilter[i]->uniqid != par.PARAM_MAX_SEQS.uniqid && par.prefilter[i]->uniqid != par.PARAM_RES_LIST_OFFSET.uniqid){
                prefilter.push_back(par.prefilter[i]);
            }
        }

        // correct Eval threshold for inverted search
        const size_t queryDbSize = FileUtil::countLines(par.db1Index.c_str());
        const size_t targetDbSize = FileUtil::countLines(par.db2Index.c_str());
        par.evalThr *= ((float) queryDbSize)/targetDbSize;

        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        cmd.addVariable("PREFILTER_PAR", par.createParameterString(prefilter).c_str());
        float originalEvalThr = par.evalThr;
        par.evalThr = std::numeric_limits<float>::max();
        cmd.addVariable("SWAP_PAR", par.createParameterString(par.swapresult).c_str());
        par.evalThr = originalEvalThr;
        cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        cmd.addVariable("SORTRESULT_PAR", par.createParameterString(par.sortresult).c_str());
        cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
        cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
        par.covMode = originalCovMode;

        program = tmpDir + "/searchslicedtargetprofile.sh";
        FileUtil::writeFile(program, searchslicedtargetprofile_sh, searchslicedtargetprofile_sh_len);
    } else if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
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

        float originalEval = par.evalThr;
        par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr  : par.evalProfile;
        for (int i = 0; i < par.numIterations; i++) {
            if (i == 0 && Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) == false) {
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



    if (isTranslatedNuclSearch == true) {
        cmd.addVariable("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : NULL);
        FileUtil::writeFile(tmpDir + "/translated_search.sh", translated_search_sh, translated_search_sh_len);
        cmd.addVariable("QUERY_NUCL", Parameters::isEqualDbtype(queryDbType,Parameters::DBTYPE_NUCLEOTIDES) ? "TRUE" : NULL);
        cmd.addVariable("TARGET_NUCL", Parameters::isEqualDbtype(targetDbType,Parameters::DBTYPE_NUCLEOTIDES) ? "TRUE" : NULL);
        cmd.addVariable("ORF_PAR", par.createParameterString(par.extractorfs).c_str());
        cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment).c_str());
        cmd.addVariable("TRANSLATE_PAR", par.createParameterString(par.translatenucs).c_str());
        cmd.addVariable("SEARCH", program.c_str());
        program = std::string(tmpDir + "/translated_search.sh");
    }else if(isNuclNuclSearch== true){
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
