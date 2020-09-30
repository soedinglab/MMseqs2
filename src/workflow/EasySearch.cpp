#include <cassert>
#include "LinsearchIndexReader.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "Parameters.h"
#include "easysearch.sh.h"

void setEasySearchDefaults(Parameters *p, bool linsearch) {
    if (linsearch) {
        p->shuffleDatabase = false;
    }
    p->sensitivity = 5.7;
    p->removeTmpFiles = true;
    p->writeLookup = false;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->orfFilter = 0;
}

void setEasySearchMustPassAlong(Parameters *p, bool linsearch) {
    if (linsearch) {
        p->PARAM_SHUFFLE.wasSet = true;
    }
    p->PARAM_S.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_ORF_FILTER.wasSet = true;
}

int doeasysearch(int argc, const char **argv, const Command &command, bool linsearch) {
    Parameters &par = Parameters::getInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasySearchDefaults(&par, linsearch);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasySearchMustPassAlong(&par, linsearch);

    bool needBacktrace = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    bool needLookup = false;

    {
        bool needSequenceDB = false;
        bool needFullHeaders = false;
        bool needSource = false;
        Parameters::getOutputFormat(par.formatAlignmentMode, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                needLookup, needSource, needTaxonomyMapping, needTaxonomy);
    }

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM || par.greedyBestHits) {
        needBacktrace = true;
    }
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }
    if(needLookup){
        par.writeLookup = true;
    }

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string target = par.filenames.back().c_str();
    cmd.addVariable("TARGET", target.c_str());
    par.filenames.pop_back();

    if(needTaxonomy || needTaxonomyMapping){
        Parameters::checkIfTaxDbIsComplete(target);
    }

    if (linsearch) {
        const bool isIndex = LinsearchIndexReader::searchForIndex(target).empty() == false;
        cmd.addVariable("INDEXEXT", isIndex ? ".linidx" : NULL);
        cmd.addVariable("SEARCH_MODULE", "linsearch");
        cmd.addVariable("LINSEARCH", "TRUE");
        cmd.addVariable("CREATELININDEX_PAR", par.createParameterString(par.createlinindex).c_str());
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.linsearchworkflow, true).c_str());
    } else {
        const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
        cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);
        cmd.addVariable("SEARCH_MODULE", "search");
        cmd.addVariable("LINSEARCH", NULL);
        cmd.addVariable("CREATELININDEX_PAR", NULL);
        cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow, true).c_str());

    }
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : NULL);
    cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("CREATEDB_QUERY_PAR", par.createParameterString(par.createdb).c_str());
    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("CONVERT_PAR", par.createParameterString(par.convertalignments).c_str());
    cmd.addVariable("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult).c_str());

    std::string program = tmpDir + "/easysearch.sh";
    FileUtil::writeFile(program, easysearch_sh, easysearch_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}

int easysearch(int argc, const char **argv, const Command &command) {
    return doeasysearch(argc, argv, command, false);
}

int easylinsearch(int argc, const char **argv, const Command &command) {
    return doeasysearch(argc, argv, command, true);
}
