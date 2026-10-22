#include "Parameters.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "DBReader.h"

#include "lin8clust.sh.h"

#include <cassert>

extern void setLinclustWorkflowDefaults(Parameters *p);

int lin8clust(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setLinclustWorkflowDefaults(&par);
    par.clustHash = true;
    par.linclustIterations = 1;
    par.diagonalMode = Parameters::DIAGONAL_MODE_UNGAPPED;
    par.clusteringMode = Parameters::GREEDY;
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription(par.PARAM_S, "Sensitivity will be automatically determined but can be adjusted", NULL, par.PARAM_S.category | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.lin8clustworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    const int dbType = FileUtil::parseDbType(par.db1.c_str());
    if ((DBReader<DBKeyType>::getExtendedDbtype(dbType) & Parameters::DBTYPE_EXTENDED_LINFINITY) == 0) {
        Debug(Debug::ERROR) << "lin8clust clusters --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases only. "
                            << "Use createdb --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " or easy-lin8clust.\n";
        EXIT(EXIT_FAILURE);
    }

    const bool nonSymetric = (par.covMode == Parameters::COV_MODE_TARGET || par.covMode == Parameters::COV_MODE_QUERY);
    par.includeCountTable = false;
    par.countTableIteration = 0;
    if (par.PARAM_INCLUDE_ADJACENCY.wasSet == false && par.PARAM_NUM_ADJACENCY.wasSet == false) {
        par.adjIteration = nonSymetric ? Parameters::CLUST_LINEAR_DEFAULT_NUM_ADJACENCY
                                       : Parameters::CLUST_LINEAR_SYMMETRIC_NUM_ADJACENCY;
    }

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.nodeRunner.c_str());
    cmd.addVariable("CLUSTHASH", par.clustHash ? "TRUE" : NULL);
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    cmd.addVariable("ALIGN2CLUST_PAR", par.createParameterString(par.align2clust).c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.mergeclusters).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    // # 0. clust hash
    par.seqIdThr = std::max(0.9f, par.seqIdThr);
    par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    cmd.addVariable("CLUSTHASHFAST_PAR", par.createParameterString(par.clusthashfast).c_str());

    std::string program = tmpDir + "/lin8clust.sh";
    FileUtil::writeFile(program, lin8clust_sh, lin8clust_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Unreachable
    assert(false);
    return 0;
}
