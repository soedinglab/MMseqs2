#include "Parameters.h"
#include "Util.h"
#include "DBWriter.h"
#include "CommandCaller.h"
#include "Command.h"
#include "Debug.h"
#include "FileUtil.h"
#include "MMseqsMPI.h"

#include "linclust.sh.h"

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include <dirent.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

void setLinclustWorkflowDefaults(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->evalThr = 0.001;
    p->seqIdThr = 0.9;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID; // set alignmentmode 3 as a default in linclust2
    p->linclustVersion = Parameters::LINCLUST_VERSION2;
    p->clustHash = false;
}

// ---------------------------------------------------------------------------------
// Linclust v2 in-process coordinator.
//
// A single outer `mpirun`/`srun` launches one "mmseqs linclust" process per rank and
// owns the MPI communicator for the whole run. Rank 0 drives every serial, filesystem-
// mutating stage (clusthash, createsubdb, mergeclusters, mergedbs/mvdb, pickconsensus-
// repfast, cleanup) while every other rank simply skips those blocks; kmermatcher and
// align2clust are the only collective stages and are called in-process on every rank so
// they can participate directly in the already-initialized communicator. There is no
// nested RUNNER/--mpi-runner spawn anywhere in this path.
// ---------------------------------------------------------------------------------
namespace {

// Splits a whitespace-separated parameter string, as produced by
// Parameters::createParameterString, into argv tokens. This mirrors the word-splitting
// a POSIX shell performs on an unquoted "$FOO_PAR" expansion, which is how these exact
// strings were consumed before Linclust v2 became an in-process coordinator;
// createParameterString base64-encodes any value containing whitespace specifically so
// this simple split remains lossless.
std::vector<std::string> tokenize(const std::string &s) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) {
        tokens.push_back(tok);
    }
    return tokens;
}

bool dbExists(const std::string &path) {
    return FileUtil::fileExists((path + ".dbtype").c_str());
}

bool shouldRunCollectiveStage(const std::string &outputDb) {
    const bool shouldRun = MMseqsMPI::isMaster() && dbExists(outputDb) == false;
    return MMseqsMPI::broadcast(shouldRun);
}

void appendTokens(std::vector<std::string> &args, const std::vector<std::string> &tokens) {
    args.insert(args.end(), tokens.begin(), tokens.end());
}

// Runs an MMseqs leaf command in-process (no fork/exec), reusing the current process'
// already-initialized MPI communicator (if any).
//
// Every MMseqsParameter is shared by uniqid across all commands that expose it (e.g.
// --threads is the same object for kmermatcher, align2clust, createsubdb, ...), and
// Parameters::parseParameters() never clears MMseqsParameter::wasSet itself -- it is a
// "seen once" latch meant to catch a flag repeated twice on one command line. Replaying
// parseParameters for a second in-process command whose argv repeats a flag already
// consumed by an earlier in-process command would otherwise misreport a "duplicate
// parameter" and abort. Parameters::resetWasSet() clears exactly the parameter list the
// upcoming command owns, right before it parses, which keeps every in-process invocation
// as clean as a fresh process without disturbing any other command's bookkeeping.
int runInProcess(const char *name, std::vector<std::string> &args) {
    const Command *cmd = getCommandByName(name);
    if (cmd == NULL) {
        Debug(Debug::ERROR) << "linclust: internal command \"" << name << "\" was not found.\n";
        EXIT(EXIT_FAILURE);
    }
    if (cmd->params != NULL) {
        Parameters::resetWasSet(*cmd->params);
    }
    std::vector<const char*> argv;
    argv.reserve(args.size());
    for (size_t i = 0; i < args.size(); ++i) {
        argv.push_back(args[i].c_str());
    }
    const char **argvData = argv.empty() ? NULL : argv.data();
    int status = cmd->commandFunction(static_cast<int>(argv.size()), argvData, *cmd);
    if (status != EXIT_SUCCESS) {
        Debug(Debug::ERROR) << "linclust: \"" << name << "\" failed with status " << status << ".\n";
        EXIT(EXIT_FAILURE);
    }
    return status;
}

void runSerialCommand(const char *name, const std::vector<std::string> &args) {
    const char *mmseqsBin = getenv("MMSEQS");
    if (mmseqsBin == NULL) {
        mmseqsBin = "mmseqs";
    }

    std::vector<char *> argv;
    argv.reserve(args.size() + 3);
    argv.push_back(const_cast<char *>(mmseqsBin));
    argv.push_back(const_cast<char *>(name));
    for (size_t i = 0; i < args.size(); ++i) {
        argv.push_back(const_cast<char *>(args[i].c_str()));
    }
    argv.push_back(NULL);

    pid_t pid = -1;
    const int spawnStatus = posix_spawnp(&pid, mmseqsBin, NULL, NULL, argv.data(), environ);
    if (spawnStatus != 0) {
        Debug(Debug::ERROR) << "linclust: could not start \"" << name << "\": "
                            << strerror(spawnStatus) << "\n";
        EXIT(EXIT_FAILURE);
    }

    int waitStatus = 0;
    while (waitpid(pid, &waitStatus, 0) == -1) {
        if (errno == EINTR) {
            continue;
        }
        Debug(Debug::ERROR) << "linclust: waiting for \"" << name << "\" failed: "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (WIFEXITED(waitStatus) == 0 || WEXITSTATUS(waitStatus) != EXIT_SUCCESS) {
        Debug(Debug::ERROR) << "linclust: \"" << name << "\" failed";
        if (WIFEXITED(waitStatus)) {
            Debug(Debug::ERROR) << " with status " << WEXITSTATUS(waitStatus);
        } else if (WIFSIGNALED(waitStatus)) {
            Debug(Debug::ERROR) << " after signal " << WTERMSIG(waitStatus);
        }
        Debug(Debug::ERROR) << ".\n";
        EXIT(EXIT_FAILURE);
    }
}

void removeDbIfExists(const std::string &path) {
    if (dbExists(path)) {
        std::vector<std::string> args(1, path);
        runSerialCommand("rmdb", args);
    }
}

// Non-conditional counterpart of removeDbIfExists(), for intermediates that this
// coordinator always produces by the time cleanup runs (matches the previous shell
// workflow, which removed these unconditionally rather than guarding every rmdb with a
// notExists check).
void removeDb(const std::string &path) {
    std::vector<std::string> args(1, path);
    runSerialCommand("rmdb", args);
}

void moveDb(const std::string &from, const std::string &to) {
    std::vector<std::string> args;
    args.push_back(from);
    args.push_back(to);
    runSerialCommand("mvdb", args);
}

// In-process equivalent of `awk '{print $1}' indexFile > outFile`: writes the first
// whitespace/tab-delimited column of a plain-text .index file, one key per line.
void writeFirstIndexColumn(const std::string &indexFile, const std::string &outFile) {
    FILE *in = fopen(indexFile.c_str(), "r");
    if (in == NULL) {
        Debug(Debug::ERROR) << "linclust: could not open " << indexFile << " for reading.\n";
        EXIT(EXIT_FAILURE);
    }
    FILE *out = fopen(outFile.c_str(), "w");
    if (out == NULL) {
        fclose(in);
        Debug(Debug::ERROR) << "linclust: could not open " << outFile << " for writing.\n";
        EXIT(EXIT_FAILURE);
    }
    char *line = NULL;
    size_t lineCap = 0;
    ssize_t len;
    while ((len = getline(&line, &lineCap, in)) != -1) {
        size_t fieldLen = 0;
        while (fieldLen < static_cast<size_t>(len)
               && line[fieldLen] != '\t' && line[fieldLen] != ' '
               && line[fieldLen] != '\n' && line[fieldLen] != '\r') {
            fieldLen++;
        }
        if (fieldLen > 0) {
            fwrite(line, 1, fieldLen, out);
            fputc('\n', out);
        }
    }
    free(line);
    fclose(out);
    fclose(in);
}

// Removes a flat scratch directory (no nested subdirectories are ever created inside
// it), mirroring `rm -rf` for the pickconsensusrepfast tmp directory. Same pattern as
// Align2ClustCheckpoint::cleanup().
void removeFlatDirectory(const std::string &dir) {
    DIR *d = opendir(dir.c_str());
    if (d != NULL) {
        struct dirent *entry;
        while ((entry = readdir(d)) != NULL) {
            const std::string name = entry->d_name;
            if (name == "." || name == "..") {
                continue;
            }
            FileUtil::remove((dir + "/" + name).c_str());
        }
        closedir(d);
    }
    rmdir(dir.c_str());
}

// Baked, tokenized argv fragments for every Linclust v2 stage, computed once up-front
// (identically on every rank, since every rank parses identical CLI argv) before any
// stage runs. This mirrors the original shell workflow, where every "*_PAR" string was
// built before the generated script was even written to disk.
struct LinclustV2Options {
    std::vector<std::string> kmermatcherPar;
    std::vector<std::string> kmermatcherPar2;
    std::vector<std::string> align2clustPar;
    std::vector<std::string> clusthashPar;
    std::vector<std::string> clusthashClustPar;
    std::vector<std::string> pickrepPar;
    std::vector<std::string> verbosityPar;
    bool clustHashEnabled = false;
    bool switchConsensusRep = false;
    bool keepSwitchAln = false;
    bool removeTmpFiles = false;
};

void runPickConsensusRepFastSubprocess(std::vector<std::string> &args) {
    runSerialCommand("pickconsensusrepfast", args);
}

// sourceDb/outputDb/tmpDir are taken by value on purpose: the caller passes par.db1,
// par.db2 and a path derived from par.db3, but Parameters is a process-wide singleton
// that every in-process leaf command (runInProcess -> parseParameters) overwrites with
// its own database arguments. Binding const references to those members would make the
// coordinator's own paths silently change to the last stage's argv.
int runLinclustV2(const std::string sourceDb, const std::string outputDb, const std::string tmpDir,
                   const LinclustV2Options &opts) {
    std::string input = sourceDb;

    // 0. clusthash: rank 0 only, deterministic given identical par on every rank.
    if (opts.clustHashEnabled) {
        const std::string clusthashDb = tmpDir + "/input_clusthash";
        const std::string clusthashClustDb = tmpDir + "/input_clusthash_clust";
        const std::string orderFile = tmpDir + "/order_clusthash_redundancy";
        const std::string redundancyDb = tmpDir + "/input_clusthash_redundancy";
        if (MMseqsMPI::isMaster()) {
            if (dbExists(clusthashDb) == false) {
                std::vector<std::string> args;
                args.push_back(input);
                args.push_back(clusthashDb);
                appendTokens(args, opts.clusthashPar);
                runSerialCommand("clusthash", args);
            }
            if (dbExists(clusthashClustDb) == false) {
                std::vector<std::string> args;
                args.push_back(input);
                args.push_back(clusthashDb);
                args.push_back(clusthashClustDb);
                appendTokens(args, opts.clusthashClustPar);
                runSerialCommand("clust", args);
            }
            writeFirstIndexColumn(clusthashClustDb + ".index", orderFile);
            if (dbExists(redundancyDb) == false) {
                std::vector<std::string> args;
                args.push_back(orderFile);
                args.push_back(sourceDb);
                args.push_back(redundancyDb);
                appendTokens(args, opts.verbosityPar);
                args.push_back("--subdb-mode");
                args.push_back("1");
                runSerialCommand("createsubdb", args);
            }
        }
        // Every rank needs the same substituted input for the (collective)
        // kmermatcher/align2clust calls below, not just rank 0.
        input = redundancyDb;
    }
    MMseqsMPI::barrier();

    // 1. k-mer matching: every rank, in-process.
    const std::string prefDb = tmpDir + "/pref";
    if (shouldRunCollectiveStage(prefDb)) {
        std::vector<std::string> args;
        args.push_back(input);
        args.push_back(prefDb);
        appendTokens(args, opts.kmermatcherPar);
        runInProcess("kmermatcher", args);
    }
    MMseqsMPI::barrier();

    // 2. Alignment + clustering (ungapped/gapped banded-block aligner): every rank.
    const std::string cluDb = tmpDir + "/clu";
    if (shouldRunCollectiveStage(cluDb)) {
        std::vector<std::string> args;
        args.push_back(input);
        args.push_back(prefDb);
        args.push_back(cluDb);
        appendTokens(args, opts.align2clustPar);
        runInProcess("align2clust", args);
    }
    MMseqsMPI::barrier();

    // 2a. Merge clusthash pre-clusters with alignment clusters: rank 0 only.
    std::string cluDbFinal = cluDb;
    if (opts.clustHashEnabled) {
        cluDbFinal = tmpDir + "/clu_merged";
        if (MMseqsMPI::isMaster() && dbExists(cluDbFinal) == false) {
            std::vector<std::string> args;
            args.push_back(sourceDb);
            args.push_back(cluDbFinal);
            args.push_back(tmpDir + "/input_clusthash_clust");
            args.push_back(cluDb);
            runSerialCommand("mergeclusters", args);
        }
    }
    MMseqsMPI::barrier();

    // 3. Refinement pass: re-cluster representative sequences. rank 0 only.
    const std::string inputRepDb = tmpDir + "/input_rep";
    if (MMseqsMPI::isMaster() && dbExists(inputRepDb) == false) {
        std::vector<std::string> args;
        args.push_back(cluDbFinal);
        args.push_back(input);
        args.push_back(inputRepDb);
        appendTokens(args, opts.verbosityPar);
        args.push_back("--subdb-mode");
        args.push_back("1");
        runSerialCommand("createsubdb", args);
    }
    MMseqsMPI::barrier();

    const std::string prefRepDb = tmpDir + "/pref_rep";
    if (shouldRunCollectiveStage(prefRepDb)) {
        std::vector<std::string> args;
        args.push_back(inputRepDb);
        args.push_back(prefRepDb);
        appendTokens(args, opts.kmermatcherPar2);
        runInProcess("kmermatcher", args);
    }
    MMseqsMPI::barrier();

    const std::string cluRepDb = tmpDir + "/clu_rep";
    if (shouldRunCollectiveStage(cluRepDb)) {
        std::vector<std::string> args;
        args.push_back(inputRepDb);
        args.push_back(prefRepDb);
        args.push_back(cluRepDb);
        appendTokens(args, opts.align2clustPar);
        args.push_back("--filter-cludb-file");
        args.push_back(cluDbFinal);
        args.push_back("--filter-seqdb-file");
        args.push_back(sourceDb);
        runInProcess("align2clust", args);
    }
    MMseqsMPI::barrier();

    // Final merge: rank 0 only. Matches the pre-existing v2 behavior of not passing any
    // extra mergeclusters flags (MERGECLU_PAR was never populated for Linclust v2).
    if (MMseqsMPI::isMaster() && dbExists(outputDb) == false) {
        std::vector<std::string> args;
        args.push_back(sourceDb);
        args.push_back(outputDb);
        args.push_back(cluDbFinal);
        args.push_back(cluRepDb);
        runSerialCommand("mergeclusters", args);
    }
    MMseqsMPI::barrier();

    if (MMseqsMPI::isMaster()) {
        // Expose alignment results (only produced when --include-align-files is set).
        // The two align2clust passes each emit their own alignments; union them keyed
        // by the final representatives (outputDb) and keep only lines whose target is
        // a member of that final cluster (--merge-filter-target), so the result has
        // exactly one entry per cluster containing exactly that cluster's
        // rep->member alignments.
        const std::string cluAlnDb = tmpDir + "/clu_aln";
        const std::string cluRepAlnDb = tmpDir + "/clu_rep_aln";
        if (dbExists(cluAlnDb)) {
            if (dbExists(cluRepAlnDb)) {
                std::vector<std::string> args;
                args.push_back(outputDb);
                args.push_back(outputDb + "_aln");
                args.push_back(cluAlnDb);
                args.push_back(cluRepAlnDb);
                args.push_back("--merge-filter-target");
                args.push_back("1");
                appendTokens(args, opts.verbosityPar);
                runSerialCommand("mergedbs", args);
            } else {
                moveDb(cluAlnDb, outputDb + "_aln");
            }
        }

        // Optionally replace representatives by the most profile-consistent observed
        // member, reusing the alignments in outputDb_aln (no profile-vs-member
        // realignment).
        if (opts.switchConsensusRep) {
            const std::string switchedDb = tmpDir + "/clu_switched";
            const std::string switchTmp = tmpDir + "/switch_tmp";
            std::vector<std::string> args;
            args.push_back(sourceDb);
            args.push_back(outputDb);
            args.push_back(switchedDb);
            args.push_back(switchTmp);
            appendTokens(args, opts.pickrepPar);
            runPickConsensusRepFastSubprocess(args);

            removeDb(outputDb);
            moveDb(switchedDb, outputDb);
            if (opts.keepSwitchAln == false) {
                removeDbIfExists(outputDb + "_aln");
            }
            removeFlatDirectory(switchTmp);
        }

        if (opts.removeTmpFiles) {
            removeDb(prefDb);
            removeDb(cluDb);
            removeDb(inputRepDb);
            removeDbIfExists(inputRepDb + "_h");
            removeDb(prefRepDb);
            removeDb(cluRepDb);
            removeDbIfExists(tmpDir + "/clu_aln");
            removeDbIfExists(tmpDir + "/clu_rep_aln");
            if (opts.clustHashEnabled) {
                removeDb(tmpDir + "/input_clusthash");
                removeDb(tmpDir + "/input_clusthash_clust");
                removeDb(tmpDir + "/input_clusthash_redundancy");
                removeDb(cluDbFinal);
                FileUtil::remove((tmpDir + "/order_clusthash_redundancy").c_str());
            }
        }
    }
    MMseqsMPI::barrier();

    return EXIT_SUCCESS;
}

} // namespace

int linclust(int argc, const char **argv, const Command& command) {
    // A single outer mpirun/srun launches this "linclust" process itself for Linclust
    // v2 (see runLinclustV2 above); establish the persistent communicator here, before
    // any argument parsing, exactly like every other MPI-aware command does.
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    setLinclustWorkflowDefaults(&par);
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription(par.PARAM_S, "Sensitivity will be automatically determined but can be adjusted", NULL, par.PARAM_S.category | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.linclustVersion == 1 && MMseqsMPI::numProc > 1) {
        Debug(Debug::ERROR) << "Linclust v1 (--linclust-version 1) does not support multi-rank execution "
                             << "(this process is one of " << MMseqsMPI::numProc << " ranks under an active "
                             << "mpirun/srun). Only Linclust v2 (the default) can run across multiple ranks; "
                             << "either drop --linclust-version 1 or launch this command with a single rank.\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.linclustVersion == 1) {
        // The legacy v1 implementation replaces this process with its shell workflow.
        // It does not reuse the outer communicator, so release a single-rank MPI
        // runtime before exec. Multi-rank v1 was rejected above.
        MMseqsMPI::finalize();
    } else if (par.runner.empty() == false) {
        Debug(Debug::ERROR) << "Linclust v2 no longer accepts --mpi-runner or RUNNER. Launch the complete "
                             << "command instead, for example: mpirun -np N mmseqs linclust <inputDB> "
                             << "<outputDB> <tmp>.\n";
        EXIT(EXIT_FAILURE);
    }

    std::string tmpBase = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.linclustworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpBase + "/latest");
    }
    // Every rank parses identical argv, so tmpBase/hash are identical everywhere; only
    // rank 0 may touch the filesystem here (concurrent mkdir from every rank would
    // race), everyone else just needs the resulting path once rank 0 is done creating
    // it, hence the barrier right after.
    std::string tmpDir;
    if (MMseqsMPI::isMaster()) {
        tmpDir = FileUtil::createTemporaryDirectory(tmpBase, hash);
    } else {
        tmpDir = tmpBase + "/" + hash;
    }
    MMseqsMPI::barrier();
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller *legacyCmd = NULL;
    if (par.linclustVersion == 1) {
        legacyCmd = new CommandCaller();
        legacyCmd->addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
        legacyCmd->addVariable("RUNNER", par.runner.c_str());
        legacyCmd->addVariable("LINCLUST_MODULE", "linclust1");
    }

    // Optionally switch to profile-consensus representatives after clustering.
    // This reuses the representative-to-member alignments, so force their creation.
    const bool writeAlnFiles = par.PARAM_INCLUDE_ALIGN_FILES.wasSet;
    if (par.switchConsensusRep && par.linclustVersion != Parameters::LINCLUST_VERSION2) {
        Debug(Debug::WARNING) << "--switch-consensus-rep requires --linclust-version 2; ignoring.\n";
        par.switchConsensusRep = false;
    }
    if (par.switchConsensusRep) {
        par.includeAlignFiles = true;
        par.addBacktrace = true;
    }
    // save some values to restore them later
    MultiParam<NuclAA<int>>alphabetSize = par.alphabetSize;
    size_t kmerSize = par.kmerSize;
    bool kmerSizeWasSet = false;
    bool alphabetSizeWasSet = false;
    bool clusterModeSet = false;
    bool includeCountTableSet = false;
    for (size_t i = 0; i < par.linclustworkflow.size(); i++) {
        if (par.linclustworkflow[i]->uniqid == par.PARAM_K.uniqid && par.linclustworkflow[i]->wasSet) {
            kmerSizeWasSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid && par.linclustworkflow[i]->wasSet) {
            alphabetSizeWasSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_CLUSTER_MODE.uniqid && par.linclustworkflow[i]->wasSet) {
            clusterModeSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_INCLUDE_COUNTTABLE.uniqid && par.linclustworkflow[i]->wasSet) {
            includeCountTableSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_NUM_COUNTS.uniqid && par.linclustworkflow[i]->wasSet) {
            includeCountTableSet = true;
        }

    }

    const bool nonSymetric = (par.covMode == Parameters::COV_MODE_TARGET || par.covMode == Parameters::COV_MODE_QUERY);
    if (clusterModeSet == false){
        if (nonSymetric) {
            par.clusteringMode = Parameters::GREEDY_MEM;
        } else {
            par.clusteringMode = Parameters::SET_COVER;
        }
        std::string cluMode = (par.clusteringMode==Parameters::GREEDY_MEM) ? "GREEDY MEM" : "SET COVER";
        Debug(Debug::INFO) << "Set cluster mode " << cluMode << ".\n";
    }

    if (includeCountTableSet == false) {
        if (nonSymetric) {
            par.includeCountTable = false;
            par.countTableIteration = 0;
        } else {
            par.includeCountTable = true;
        }
    }

    if (kmerSizeWasSet == false) {
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    }
    if (alphabetSizeWasSet == false) {
        par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5));
    }

    const int dbType = FileUtil::parseDbType(par.db1.c_str());
    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE)) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases.\n";
        EXIT(EXIT_FAILURE);
    }

    LinclustV2Options v2Opts;
    if (par.linclustVersion == 1) {
        legacyCmd->addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
        // filter by diagonal in case of AA (do not filter for nucl, profiles, ...)
        legacyCmd->addVariable("FILTER", Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_AMINO_ACIDS) ? "1" : NULL);
        legacyCmd->addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
        legacyCmd->addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
        legacyCmd->addVariable("VERBOSITYANDCOMPRESS", par.createParameterString(par.threadsandcompression).c_str());

        par.alphabetSize = alphabetSize;
        par.kmerSize = kmerSize;
        // # 2. Hamming distance pre-clustering
        par.rescoreMode = Parameters::RESCORE_MODE_HAMMING;
        par.filterHits = false;
        float prevSeqId = par.seqIdThr;
        // hamming distance does not work well with seq. id < 0.5 since it does not have an e-value criteria
        par.seqIdThr = std::max(0.5f, par.seqIdThr);
        // also coverage should not be under 0.5
        float prevCov = par.covThr;
        par.covThr = std::max(0.5f, par.covThr);
        legacyCmd->addVariable("HAMMING_PAR", par.createParameterString(par.rescorediagonal).c_str());
        // set it back to old value
        par.covThr = prevCov;
        par.seqIdThr = prevSeqId;
        par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;

        // # 3. Ungapped alignment filtering
        par.filterHits = true;
        legacyCmd->addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());

        // # 4. Local gapped sequence alignment.
        if (isUngappedMode) {
            const int originalRescoreMode = par.rescoreMode;
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            legacyCmd->addVariable("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal).c_str());
            par.rescoreMode = originalRescoreMode;
        } else {
            legacyCmd->addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
        }
        // # 5. Clustering using greedy set cover.
        legacyCmd->addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
        legacyCmd->addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());

    } else if (par.linclustVersion == 2) {
        par.alphabetSize = alphabetSize;
        par.kmerSize = kmerSize;
        bool prevspacedKmer = par.spacedKmer;
        bool prevmaskMode = par.maskMode;
        par.spacedKmer = false;
        par.maskMode = false;
        v2Opts.kmermatcherPar = tokenize(par.createParameterString(par.kmermatcher));

        v2Opts.verbosityPar = tokenize(par.createParameterString(par.onlyverbosity));
        v2Opts.align2clustPar = tokenize(par.createParameterString(par.align2clust));

        par.spacedKmer = true;
        par.kmersPerSequenceScale = 0.1;
        v2Opts.kmermatcherPar2 = tokenize(par.createParameterString(par.kmermatcher));

        par.spacedKmer = prevspacedKmer;
        par.maskMode = prevmaskMode;
    }
    float prevSeqId = par.seqIdThr;
    // # 0. clust hash
    par.seqIdThr = std::max(0.9f, par.seqIdThr);
    par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    const std::string clusthashParameterString = par.createParameterString(par.clusthash);
    std::vector<std::string> clusthashPar = tokenize(clusthashParameterString);
    par.seqIdThr = prevSeqId;
    par.alphabetSize = alphabetSize;
    const std::string clusthashClustParameterString = par.createParameterString(par.clust);
    std::vector<std::string> clusthashClustPar = tokenize(clusthashClustParameterString);

    const std::string pickrepParameterString = par.createParameterString(par.pickrepprofile);
    std::vector<std::string> pickrepPar = tokenize(pickrepParameterString);

    if (par.linclustVersion == 2) {
        v2Opts.clusthashPar = clusthashPar;
        v2Opts.clusthashClustPar = clusthashClustPar;
        v2Opts.pickrepPar = pickrepPar;
        v2Opts.clustHashEnabled = par.clustHash;
        v2Opts.switchConsensusRep = par.switchConsensusRep;
        v2Opts.keepSwitchAln = par.switchConsensusRep && writeAlnFiles;
        v2Opts.removeTmpFiles = par.removeTmpFiles;
        return runLinclustV2(par.db1, par.db2, tmpDir, v2Opts);
    }

    legacyCmd->addVariable("CLUSTHASH", par.clustHash ? "TRUE" : NULL);
    legacyCmd->addVariable("CLUSTHASH_PAR", clusthashParameterString.c_str());
    legacyCmd->addVariable("CLUSTHASH_CLUST_PAR", clusthashClustParameterString.c_str());
    legacyCmd->addVariable("SWITCH_CONSENSUS_REP", par.switchConsensusRep ? "TRUE" : NULL);
    legacyCmd->addVariable("KEEP_SWITCH_ALN", (par.switchConsensusRep && writeAlnFiles) ? "TRUE" : NULL);
    legacyCmd->addVariable("PICKREP_PAR", pickrepParameterString.c_str());

    std::string program = tmpDir + "/linclust.sh";
    FileUtil::writeFile(program, linclust_sh, linclust_sh_len);
    legacyCmd->execProgram(program.c_str(), par.filenames);

    // Unreachable
    assert(false);
    return 0;
}
