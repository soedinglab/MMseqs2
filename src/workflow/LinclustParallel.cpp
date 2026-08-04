#include "ByteParser.h"
#include "CommandCaller.h"
#include "KmerPartition.h"
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
    // As linclust does (Linclust.cpp:18). Without it the alignment downgrades to
    // SCORE_COV at --min-seq-id 0, where stock stays SCORE_COV_SEQID.
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}

int linclustparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinclustParallelWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    // Only the coverage modes whose clustering this pipeline actually implements.
    //
    // For a symmetric --cov-mode, linclust switches to SET_COVER *and* enables the
    // count table's extra grouping rounds (Linclust.cpp:91-109). greedycluster
    // implements only the length-ordered greedy and the reduce runs no count-table
    // rounds, so accepting mode 0 would apply this pipeline's thresholds to a
    // different algorithm and return a silently different clustering -- measured
    // at 101 clusters against stock's 117 on a 389-sequence input, exit 0.
    if (par.covMode != Parameters::COV_MODE_TARGET && par.covMode != Parameters::COV_MODE_QUERY) {
        Debug(Debug::ERROR) << "--cov-mode " << par.covMode << " is not supported.\n"
                            << "linclust selects SET_COVER clustering and the count-table grouping "
                            << "rounds for symmetric coverage modes; this pipeline implements "
                            << "neither, so the result would differ from linclust without saying "
                            << "so. Use --cov-mode 1 (target) or 2 (query).\n";
        EXIT(EXIT_FAILURE);
    }

    // Said once, up front, rather than discovered several stages in.
    //
    // DBKeyType is a uint32_t unless the build sets MMSEQS_INT64_IDS, whose CMake
    // default is off. createdbparallel does refuse an input past the ceiling, but
    // only after scanning and planning the whole thing -- hours at the scales this
    // command exists for. Naming the limit before anything is launched turns that
    // into a one-line answer.
    if (static_cast<uint64_t>(DB_KEY_INVALID) - 1 < KmerRecord::MAX_ID) {
        Debug(Debug::WARNING)
            << "This build uses 32-bit database keys, so it can cluster at most "
            << static_cast<uint64_t>(DB_KEY_INVALID) << " sequences. The distributed pipeline is "
            << "designed for 1e11-1e12; build with -DMMSEQS_INT64_IDS=1 for those.\n";
    }

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
    // WORKER_RUNNER, not RUNNER: ten existing workflows set RUNNER from the *MPI*
    // --mpi-runner, and Parameters.cpp seeds par.runner from getenv("RUNNER"), so
    // reusing the name would hand every stage sub-process an MPI runner it never
    // asked for. --runner itself is separate from --mpi-runner because it only has
    // to start one worker per node; the workers coordinate through files.
    cmd.addVariable("WORKER_RUNNER", par.workerRunner.c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
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
