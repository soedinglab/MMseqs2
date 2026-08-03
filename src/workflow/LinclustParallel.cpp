#include "ByteParser.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"

#include "linclustparallel.sh.h"

#include <cassert>

// Shared-filesystem parallel linclust.
//
// Same result as `linclust`, computed by many independent worker processes that
// coordinate only through files. The stages are ordinary MMseqs2 commands run
// with byte-identical arguments on every worker, so a stage maps onto a Slurm
// array job and workers may join late, die, or be restarted:
//
//   mmseqs linclustparallel seqDB clusters.tsv tmp --runner "srun -n 64"
//
// The parameter surface is deliberately narrower than `linclust`'s. See the
// comment on linclustparallelworkflow in Parameters.cpp for why, and the header
// of data/workflow/linclustparallel.sh for what the stages are.
// A ByteParser value the stages will parse back. `0` means "unset" to every one
// of them and has no unit; anything else carries an explicit B so it is not read
// as megabytes.
static std::string byteValue(size_t bytes) {
    return bytes == 0 ? std::string("0") : ByteParser::format(bytes, 'B', 'l');
}

void setLinclustParallelWorkflowDefaults(Parameters *p) {
    p->covThr = 0.8;
    p->covMode = Parameters::COV_MODE_TARGET;
    p->evalThr = 0.001;
    p->seqIdThr = 0.9;
}

int linclustparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinclustParallelWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db3;
    std::string hash =
        SSTR(par.hashParameter(command.databases, par.filenames, par.linclustparallelworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    // --runner, not --mpi-runner: this only has to start one worker per node, and
    // the workers coordinate through files rather than communicating.
    cmd.addVariable("RUNNER", par.workerRunner.c_str());
    cmd.addVariable("THREADS", SSTR(par.threads).c_str());
    cmd.addVariable("MIN_SEQ_ID", SSTR(par.seqIdThr).c_str());
    cmd.addVariable("COV", SSTR(par.covThr).c_str());
    cmd.addVariable("COV_MODE", SSTR(par.covMode).c_str());
    cmd.addVariable("EVAL", SSTR(par.evalThr).c_str());
    // Through ByteParser, not as a bare number: the stages parse these with it,
    // and it reads an unsuffixed value as *megabytes*, so passing the byte count
    // straight through would inflate both limits a millionfold.
    cmd.addVariable("SPLIT_MEMORY_LIMIT", byteValue(par.splitMemoryLimit).c_str());
    cmd.addVariable("SCRATCH_BUDGET", byteValue(par.scratchBudget).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

    std::string program = tmpDir + "/linclustparallel.sh";
    FileUtil::writeFile(program, linclustparallel_sh, linclustparallel_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Unreachable
    assert(false);
    return 0;
}
