#include "Parameters.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <vector>
#include <algorithm>

#include "NodePlacement.h"

#include "lin8clust.sh.h"

void setLinclustOneshotWorkflowDefaults(Parameters *p) {
    p->covThr = 0.8;
    p->seqIdThr = 0.9;
    p->covMode = Parameters::COV_MODE_BIDIRECTIONAL;
    p->maskMode = 0;
    p->spacedKmer = false;
    p->clusteringMode = Parameters::GREEDY;
    p->clustHash = true;
}

static bool scriptOwns(const Parameters &par, const MMseqsParameter *p) {
    const int owned[] = {par.PARAM_LINCLUSTERDB_NODE_LIST.uniqid,
                         par.PARAM_LINCLUSTERDB_NODE_ID.uniqid,
                         par.PARAM_LINCLUSTERDB_NODE_COUNT.uniqid,
                         par.PARAM_LIN8_REP_RANK_BLOCK.uniqid,
                         par.PARAM_LIN8_REP_RANK_BLOCK_COUNT.uniqid,
                         par.PARAM_LIN8_REP_RANK_BLOCK_LOOKAHEAD.uniqid,
                         par.PARAM_LIN8_REP_RANK_BLOCKS.uniqid,
                         par.PARAM_LIN8_MONITOR_PID.uniqid,
                         par.PARAM_LIN8_MONITOR_INTERVAL.uniqid,
                         par.PARAM_THREADS.uniqid,
                         par.PARAM_MIN_SEQ_ID.uniqid,
                         par.PARAM_C.uniqid,
                         par.PARAM_COV_MODE.uniqid,
                         par.PARAM_INCLUDE_ALIGN_FILES.uniqid};
    for (size_t i = 0; i < sizeof(owned) / sizeof(owned[0]); i++) {
        if (p->uniqid == owned[i]) {
            return true;
        }
    }
    return false;
}

static std::string passOn(Parameters &par, const std::vector<MMseqsParameter *> &all) {
    std::vector<MMseqsParameter *> mine;
    for (size_t i = 0; i < all.size(); i++) {
        if (scriptOwns(par, all[i]) == false) {
            mine.push_back(all[i]);
        }
    }
    return par.createParameterString(mine, true);
}

int lin8clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinclustOneshotWorkflowDefaults(&par);
    par.overrideParameterDescription(par.PARAM_LINCLUSTERDB_NODE_LIST,
                                     "Every machine this runs on, in order. A machine finds itself "
                                     "in it by host name and that position is its number",
                                     NULL, par.PARAM_LINCLUSTERDB_NODE_LIST.category);
    par.overrideParameterDescription(par.PARAM_LINCLUSTERDB_NODE_ID,
                                     "Take this position in --node-list rather than the one the host "
                                     "name finds, for a scheduler that renames machines or a test "
                                     "putting several on one", NULL, par.PARAM_LINCLUSTERDB_NODE_ID.category);
    par.overrideParameterDescription(par.PARAM_CLUST_HASH,
                                     "Reduce exact duplicates away before the k-mer passes and put "
                                     "them back at the end", NULL, par.PARAM_CLUST_HASH.category);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    if (par.clusteringMode != Parameters::GREEDY && par.clusteringMode != Parameters::GREEDY_MEM) {
        Debug(Debug::ERROR) << "lin8clust clusters greedily by length, so --cluster-mode must "
                            << "be " << Parameters::GREEDY << " or " << Parameters::GREEDY_MEM << "\n";
        EXIT(EXIT_FAILURE);
    }

    const std::string tmpDir = par.filenames.back();
    par.filenames.pop_back();
    const std::string out = par.filenames.back();
    par.filenames.pop_back();
    if (FileUtil::directoryExists(tmpDir.c_str()) == false
        && FileUtil::makeDir(tmpDir.c_str()) == false) {
        Debug(Debug::ERROR) << "Cannot create tmp directory " << tmpDir << "\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    const NodePlacement node = NodePlacement::resolve(par);
    cmd.addVariable("NODES", SSTR(node.count).c_str());
    cmd.addVariable("NODE", SSTR(node.index).c_str());
    cmd.addVariable("REP_RANK_BLOCKS", SSTR(par.lin8RepRankBlocks).c_str());
    cmd.addVariable("PAIR_SPLIT_COUNT", SSTR(par.lin8RepRankBlockCount).c_str());
    cmd.addVariable("PAIR_SPLIT_LOOKAHEAD", par.PARAM_LIN8_REP_RANK_BLOCK_LOOKAHEAD.wasSet
                                                ? SSTR(par.lin8RepRankBlockLookahead).c_str()
                                                : NULL);
    cmd.addVariable("THREADS", SSTR(par.threads).c_str());
    cmd.addVariable("MONITOR_INTERVAL", SSTR(par.lin8MonitorInterval).c_str());
    cmd.addVariable("SEQID", SSTR(par.seqIdThr).c_str());
    cmd.addVariable("HASHSEQID", SSTR(std::max(0.9f, par.seqIdThr)).c_str());
    cmd.addVariable("COV", SSTR(par.covThr).c_str());
    cmd.addVariable("COVMODE", SSTR(par.covMode).c_str());
    cmd.addVariable("CLUSTHASH", par.clustHash ? "TRUE" : NULL);
    cmd.addVariable("INCLUDE_ALIGN_FILES", par.includeAlignFiles ? "TRUE" : NULL);
    cmd.addVariable("SWITCH_CONSENSUS_REP", par.switchConsensusRep ? "TRUE" : NULL);
    cmd.addVariable("REPSEQ", par.fastaSplits > 0 ? "1" : "0");
    cmd.addVariable("CREATEDB_PAR", passOn(par, par.lin8createdb).c_str());
    cmd.addVariable("HASH_PAR", passOn(par, par.lin8clusthash).c_str());
    cmd.addVariable("EXTRACT_PAR", passOn(par, par.lin8extractkmers).c_str());
    cmd.addVariable("GROUP_PAR", passOn(par, par.lin8assignedpairs).c_str());
    cmd.addVariable("PREF_PAR", passOn(par, par.lin8pref).c_str());
    cmd.addVariable("ALIGN_PAR", passOn(par, par.lin8align2clust).c_str());
    cmd.addVariable("ASSIGN_PAR", passOn(par, par.lin8align2clustmulti).c_str());
    cmd.addVariable("EXPAND_PAR", passOn(par, par.lin8mergehashredundancy).c_str());
    cmd.addVariable("CLUSTERDB_PAR", passOn(par, par.lin8createclusterdb).c_str());
    cmd.addVariable("PICKREP_PAR", passOn(par, par.lin8pickrepprofile).c_str());
    cmd.addVariable("TSV_PAR", passOn(par, par.lin8createtsv).c_str());
    cmd.addVariable("REPSEQ_PAR", passOn(par, par.lin8createrepseqfasta).c_str());

    std::string program = tmpDir + "/lin8clust.sh";
    const std::string programTmp = program + "." + SSTR(node.index) + ".tmp";
    FileUtil::writeFile(programTmp, lin8clust_sh, lin8clust_sh_len);
    FileUtil::publishAtomically(programTmp, program);
    par.filenames.push_back(out);
    par.filenames.push_back(tmpDir);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
