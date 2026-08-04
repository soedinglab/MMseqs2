/*
 * alignparallel -- the alignment stage of the distributed linclust.
 *
 * Claims one representative-key-range bucket of candidate edges at a time,
 * merges the duplicate copies the k-mer partitions produced, aligns each
 * surviving pair once, and writes the survivors.
 *
 * Why this is a separate stage rather than fused into the reduce, which is the
 * obvious design and the one tried first: alignment needs both sequences of every
 * pair, and a k-mer partition's pairs are scattered over the whole key space.
 * Measured on UniRef100, a single partition needed 1.61 GB of sequences but read
 * 83.9 GB to get them -- 52x amplification, effectively the whole database, once
 * per partition. Bucketing by representative key instead gives each worker a
 * contiguous slice of sequences, which is what makes the reads sequential.
 *
 * Two properties fall out of bucketing by representative, and both matter:
 *
 *   - Every copy of a pair, from however many k-mer partitions found it, lands in
 *     the same bucket. So the cross-partition duplicates (measured factor 1.35)
 *     are removed *before* aligning rather than paid for as duplicate alignments.
 *   - Having all copies together means stock's global accumulation of score per
 *     (pair, diagonal) can be reproduced exactly -- pick the diagonal the most
 *     k-mers agree on across all partitions -- rather than approximated within one
 *     partition.
 *
 * Neither the k-mer buckets nor a resident sequence index is needed here; the
 * sequences come through the dense companion index, addressed by key.
 */
#include "Alignment.h"
#include "Matcher.h"
#include "BlockAligner.h"
#include "CandidateEdge.h"
#include "Command.h"
#include "Debug.h"
#include "DenseIndex.h"
#include "EvalueComputation.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "Matcher.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "PartitionSequences.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

bool compareEdgeByPairAndDiagonal(const CandidateEdge &a, const CandidateEdge &b) {
    const uint64_t ra = a.getRep(), rb = b.getRep();
    if (ra != rb) return ra < rb;
    const uint64_t ma = a.getMember(), mb = b.getMember();
    if (ma != mb) return ma < mb;
    if (a.diagonal != b.diagonal) return a.diagonal < b.diagonal;
    return a.reverseStrand < b.reverseStrand;
}

// Collapses the copies of each pair into one edge, exactly as stock's merge does.
//
// Stock accumulates a score per (pair, diagonal) across every split and iteration
// and keeps the diagonal with the highest total (kmermatcher.cpp:1913-1925). Here
// the copies come from different k-mer partitions instead, but they are all
// present in this bucket, so the same accumulation is possible and the result is
// the same. Doing it per partition -- which is all the reduce can do -- gives a
// different diagonal for 0.86% of pairs; doing it here does not.
size_t mergePairCopies(std::vector<CandidateEdge> &edges) {
    if (edges.empty()) {
        return 0;
    }
    SORT_PARALLEL(edges.begin(), edges.end(), compareEdgeByPairAndDiagonal);

    size_t out = 0;
    size_t i = 0;
    while (i < edges.size()) {
        const uint64_t rep = edges[i].getRep();
        const uint64_t member = edges[i].getMember();

        int bestTotal = -1;
        int16_t bestDiagonal = edges[i].diagonal;
        uint8_t bestStrand = edges[i].reverseStrand;
        int runTotal = 0;
        int16_t runDiagonal = edges[i].diagonal;

        size_t j = i;
        while (j < edges.size() && edges[j].getRep() == rep && edges[j].getMember() == member) {
            if (edges[j].diagonal != runDiagonal) {
                runDiagonal = edges[j].diagonal;
                runTotal = 0;
            }
            runTotal += edges[j].score;
            // >= not >, so the largest diagonal wins a tie -- the same rule stock
            // applies in its merge (kmermatcher.cpp:1913), reached because the
            // edges are sorted by ascending diagonal.
            if (runTotal >= bestTotal) {
                bestTotal = runTotal;
                bestDiagonal = runDiagonal;
                bestStrand = edges[j].reverseStrand;
            }
            j++;
        }

        edges[out] = edges[i];
        edges[out].diagonal = bestDiagonal;
        edges[out].reverseStrand = bestStrand;
        edges[out].score = static_cast<uint16_t>(std::min(bestTotal, 65535));
        out++;
        i = j;
    }
    edges.resize(out);
    return out;
}

// Which worker's edges count, per k-mer partition.
//
// A partition redone after a crash leaves the dead worker's edges on disk as well
// as the redo's, and mergePairCopies sums matching (pair, diagonal) records, so
// keeping both inflates a diagonal's support. The queue named exactly one
// producer per item -- complete() keeps the first worker to record it -- so
// reading it back says which copy is authoritative.
//
// Waves are concatenated in order because each has its own queue over its own
// contiguous slice of partition space, so wave w's items are partitions
// [w * P / W, (w + 1) * P / W) and appending lands each at its own index.
std::vector<int64_t> readReduceAuthority(const std::string &edgeDir, unsigned int expectedWaves,
                                        unsigned int expectedPartitions) {
    std::vector<int64_t> authority;
    // Every wave, by count, not "until one is missing".
    //
    // This used to stop at the first absent reduce.<w>.queue and accept whatever
    // it had. Starting the align before the last wave had run therefore produced a
    // short authority vector that looked complete: every edge block naming a
    // partition past its end was kept (see EdgeBucketReader::readShard), the
    // missing partitions' edges simply were not there, and the result was a
    // successful but incomplete clustering that a restart would then preserve,
    // because the align output now looked finished too.
    for (unsigned int wave = 0; wave < expectedWaves; wave++) {
        std::vector<int64_t> workers;
        const std::string path = edgeDir + "/coord/reduce." + SSTR(wave) + ".queue";
        if (WorkQueue::readCompletedWorkers(path, workers) == false) {
            Debug(Debug::ERROR) << "The reduce recorded " << expectedWaves << " waves in "
                                << edgeDir << "/coord/edge.info, but " << path
                                << " does not exist. Run kmerreduceparallel for every wave "
                                << "before aligning.\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t i = 0; i < workers.size(); i++) {
            if (workers[i] < 0) {
                Debug(Debug::ERROR) << "Partition " << (authority.size() + i) << " of " << path
                                    << " was never completed. The reduce is unfinished; re-run it "
                                    << "before aligning.\n";
                EXIT(EXIT_FAILURE);
            }
        }
        authority.insert(authority.end(), workers.begin(), workers.end());
    }
    if (authority.size() != expectedPartitions) {
        Debug(Debug::ERROR) << "The reduce queues under " << edgeDir << "/coord cover "
                            << authority.size() << " partitions, but edge.info records "
                            << expectedPartitions << ". The edge directory mixes output from two "
                            << "different runs.\n";
        EXIT(EXIT_FAILURE);
    }
    return authority;
}

size_t readBucket(const std::string &edgeDir, unsigned int bucket,
                  const std::vector<int64_t> &authority, std::vector<CandidateEdge> &out) {
    const std::vector<std::string> shards = EdgeBucketWriter::shardFiles(edgeDir, bucket);
    for (size_t i = 0; i < shards.size(); i++) {
        EdgeBucketReader::readShard(shards[i], authority, out);
    }
    return out.size();
}

// The pass-2 acceptance gate stock applies with --filter-cludb-file.
//
// Before a representative q may take a member t, stock additionally requires that
// **every sequence of t's pass-1 cluster** also aligns to q (Align2clust.cpp:657-730,
// the `allpass` loop). It is strictly more conservative than the plain alignment,
// so leaving it out merges too much: measured, 902,641 clusters against stock's
// 902,795.
//
// t is a pass-1 representative, so it is dense in *sub-key* space -- which makes the
// lookup a CSR index over sub-keys: an offset per representative plus one key per
// sequence.
//
// That CSR is built once by createrepdb and *paged*, not loaded. Every align
// worker used to build it for itself -- the whole key map (8 B x R), an offset per
// representative (8 B x R) and every pass-1 member key (8 B x N) -- by streaming
// the pass-1 TSV twice with a binary search per line. At 1e11 with R/N = 0.4 that
// is ~1.44 TB resident per worker on a 2 TB node, before the bucket's own edges
// and sequence arena, and ~14.4 TB at 1e12; the two TSV passes were themselves
// 2 x 1e11 cache-missing probes into a 320 GB array, in every worker, over a
// ~2.8 TB file.
//
// Now a worker reads only the slice its current bucket refers to. Resident cost is
// O(bucket), and the whole thing is one sequential pread per contiguous run of
// sub-keys because the bucket's member sub-keys are already sorted.
struct FilterGate {
    // The bucket's slice, rebuilt per bucket by loadSlice().
    std::vector<uint64_t> subs;      // sorted, distinct member sub-keys of this bucket
    std::vector<uint64_t> selfKey;   // keymap[sub], parallel to subs
    std::vector<uint64_t> start;     // index into `members`, parallel to subs
    std::vector<uint64_t> count;     // cluster size, parallel to subs
    std::vector<uint64_t> members;   // the original keys of those clusters, concatenated
    PartitionSequences fullSeqs;

    uint64_t repCount;      // sub-keys the gate covers
    uint64_t memberCount;   // rows in the members file
    int keymapFd;
    int offsetsFd;
    int membersFd;

    FilterGate(const std::string &fullDb)
        : fullSeqs(fullDb), repCount(0), memberCount(0), keymapFd(-1), offsetsFd(-1),
          membersFd(-1) {}

    ~FilterGate() {
        if (keymapFd >= 0) close(keymapFd);
        if (offsetsFd >= 0) close(offsetsFd);
        if (membersFd >= 0) close(membersFd);
    }

    static int openReadOnly(const std::string &path) {
        const int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            Debug(Debug::ERROR) << "Cannot open " << path << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        return fd;
    }

    static void readAt(int fd, void *dst, size_t bytes, size_t offset, const std::string &what) {
        char *p = static_cast<char *>(dst);
        size_t done = 0;
        while (done < bytes) {
            const ssize_t got = pread(fd, p + done, bytes - done, static_cast<off_t>(offset + done));
            if (got <= 0) {
                if (got < 0 && errno == EINTR) continue;
                Debug(Debug::ERROR) << "Cannot read " << what << ": "
                                    << (got < 0 ? strerror(errno) : "unexpected end of file")
                                    << "\n";
                EXIT(EXIT_FAILURE);
            }
            done += static_cast<size_t>(got);
        }
    }

    // Opens the CSR and validates it against the key map, without reading either.
    void open(const std::string &keymapFile, const std::string &repDb) {
        const std::string offsetsPath = repDb + ".gate.offsets";
        const std::string membersPath = repDb + ".gate.members";
        if (FileUtil::fileExists(offsetsPath.c_str()) == false) {
            Debug(Debug::ERROR) << "No filter gate next to " << repDb << " (" << offsetsPath
                                << " is missing). It is built by createrepdb; re-run that stage "
                                << "with this build.\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t keymapBytes = FileUtil::getFileSize(keymapFile);
        repCount = keymapBytes / sizeof(uint64_t);
        const size_t offsetBytes = FileUtil::getFileSize(offsetsPath);
        if (offsetBytes != (repCount + 1) * sizeof(uint64_t)) {
            Debug(Debug::ERROR) << offsetsPath << " holds " << (offsetBytes / sizeof(uint64_t))
                                << " offsets but " << keymapFile << " holds " << repCount
                                << " representatives. They are from different runs of "
                                << "createrepdb.\n";
            EXIT(EXIT_FAILURE);
        }
        memberCount = FileUtil::getFileSize(membersPath) / sizeof(uint64_t);
        keymapFd = openReadOnly(keymapFile);
        offsetsFd = openReadOnly(offsetsPath);
        membersFd = openReadOnly(membersPath);
    }

    // Fetches the clusters of `wanted` (sorted, distinct sub-keys) and reports the
    // original keys they contain, so the caller can prefetch those sequences.
    void loadSlice(const std::vector<uint64_t> &wanted, std::vector<uint64_t> &memberKeysOut) {
        subs = wanted;
        selfKey.assign(subs.size(), 0);
        start.assign(subs.size(), 0);
        count.assign(subs.size(), 0);
        members.clear();
        memberKeysOut.clear();
        if (subs.empty()) {
            return;
        }

        // Offsets and key map, in coalesced ascending runs. Both are fixed width
        // and indexed by sub-key, so a run of nearby sub-keys is one pread.
        const size_t coalesce = 64 * 1024 / sizeof(uint64_t);
        std::vector<uint64_t> block;
        size_t i = 0;
        while (i < subs.size()) {
            size_t j = i;
            while (j + 1 < subs.size() && subs[j + 1] - subs[i] < coalesce) {
                j++;
            }
            // One past the last, because a cluster needs offsets[s] and offsets[s+1].
            const uint64_t from = subs[i];
            const uint64_t to = subs[j] + 2;
            block.resize(static_cast<size_t>(to - from));
            readAt(offsetsFd, block.data(), block.size() * sizeof(uint64_t),
                   static_cast<size_t>(from) * sizeof(uint64_t), "the filter gate offsets");
            for (size_t k = i; k <= j; k++) {
                const size_t at = static_cast<size_t>(subs[k] - from);
                start[k] = block[at];
                count[k] = block[at + 1] - block[at];
            }
            block.resize(static_cast<size_t>(subs[j] - subs[i] + 1));
            readAt(keymapFd, block.data(), block.size() * sizeof(uint64_t),
                   static_cast<size_t>(subs[i]) * sizeof(uint64_t), "the sub-key map");
            for (size_t k = i; k <= j; k++) {
                selfKey[k] = block[static_cast<size_t>(subs[k] - subs[i])];
            }
            i = j + 1;
        }

        // The member keys themselves. Clusters of ascending sub-keys occupy
        // ascending, mostly contiguous ranges of the members file, so this is a
        // forward scan; runs that are actually adjacent become one read.
        uint64_t total = 0;
        for (size_t k = 0; k < subs.size(); k++) {
            total += count[k];
        }
        members.resize(static_cast<size_t>(total));
        uint64_t at = 0;
        i = 0;
        while (i < subs.size()) {
            size_t j = i;
            while (j + 1 < subs.size() && start[j + 1] == start[j] + count[j]) {
                j++;
            }
            uint64_t run = 0;
            for (size_t k = i; k <= j; k++) {
                run += count[k];
            }
            if (run > 0) {
                readAt(membersFd, members.data() + at, static_cast<size_t>(run) * sizeof(uint64_t),
                       static_cast<size_t>(start[i]) * sizeof(uint64_t), "the filter gate members");
            }
            // start[] is rewritten as an index into the local arena.
            uint64_t local = at;
            for (size_t k = i; k <= j; k++) {
                start[k] = local;
                local += count[k];
            }
            at += run;
            i = j + 1;
        }

        memberKeysOut.assign(members.begin(), members.end());
    }

    // Index into the loaded slice for a sub-key, or -1 when the bucket did not ask
    // for it.
    int64_t find(uint64_t sub) const {
        const std::vector<uint64_t>::const_iterator it =
            std::lower_bound(subs.begin(), subs.end(), sub);
        if (it == subs.end() || *it != sub) {
            return -1;
        }
        return static_cast<int64_t>(it - subs.begin());
    }
};


// Runs stock's allpass loop for one candidate member.
//
// The member's whole pass-1 cluster must align to the representative -- ungapped
// first, then a full Smith-Waterman if that fails, exactly as
// Align2clust.cpp:669-726. A singleton pass-1 cluster passes trivially.
//
// `elementDiagonal` is the caller's, because stock's two call sites disagree: the
// ungapped-accept gate seeds every element with the *candidate pair's* diagonal
// (`:680`), while the gapped-accept gate uses 0 (`:816`). Using the pair's
// diagonal for a different sequence reads like an oversight upstream, but parity
// is the goal, so both are reproduced as they are. Passing 0 in both places is
// what left 9 of 1,000,000 sequences under a different representative than stock.
bool passesFilterGate(const FilterGate *gate, uint64_t subMember, Sequence &query,
                      Sequence &element, BlockAligner &aligner, Matcher &matcher,
                      Parameters &par, unsigned int swMode, short elementDiagonal) {
    if (gate == NULL || subMember >= gate->repCount) {
        return true;
    }
    const int64_t slot = gate->find(subMember);
    if (slot < 0) {
        // Not in this bucket's slice. The slice is built from exactly the member
        // sub-keys of the bucket's edges, so this is only reachable for a sub-key
        // no edge named -- for which the gate is never consulted.
        return true;
    }
    if (gate->count[static_cast<size_t>(slot)] <= 1) {
        return true;  // stock only runs the loop when numClu > 1
    }
    const uint64_t targetKey = gate->selfKey[static_cast<size_t>(slot)];
    const uint64_t from = gate->start[static_cast<size_t>(slot)];
    const uint64_t to = from + gate->count[static_cast<size_t>(slot)];
    for (uint64_t j = from; j < to; j++) {
        const uint64_t elementKey = gate->members[j];
        if (elementKey == targetKey) {
            continue;
        }
        unsigned int elementLen = 0;
        const char *elementSeq = gate->fullSeqs.get(elementKey, &elementLen);
        if (elementSeq == NULL || elementLen == 0) {
            return false;
        }
        element.mapSequence(0, static_cast<DBKeyType>(elementKey), elementSeq, elementLen);
        if (Util::canBeCovered(par.covThr, par.covMode, query.L, element.L) == false) {
            return false;
        }
        BlockAligner::UngappedAln_res ua = aligner.ungappedAlign(&element, elementDiagonal);
        const bool hasEvalue = (ua.eval <= par.evalThr);
        const bool hasAlnLen = (ua.alnLen >= par.alnLenThr);
        const bool hasCoverage = Util::hasCoverage(par.covThr, par.covMode, ua.qcov, ua.tcov);
        int identical = 0;
        for (int q = ua.qStart; q <= ua.qEnd; q++) {
            const char a = query.getSeqData()[q] & static_cast<unsigned char>(~0x20);
            const char b = elementSeq[ua.tStart + (q - ua.qStart)] & static_cast<unsigned char>(~0x20);
            identical += (a == b) ? 1 : 0;
        }
        const float elementSeqId =
            Util::computeSeqId(par.seqIdMode, identical, query.L, elementLen, ua.alnLen);
        const bool hasSeqId =
            elementSeqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());

        if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
            continue;
        }
        Matcher::result_t res = matcher.getSWResult(&element, static_cast<int>(elementDiagonal),
                                                    false, par.covMode, par.covThr, par.evalThr,
                                                    swMode, par.seqIdMode, false);
        if (Alignment::checkCriteria(res, false, par.evalThr, par.seqIdThr, par.alnLenThr,
                                     par.covMode, par.covThr) == false) {
            return false;
        }
    }
    return true;
}

}  // namespace

int alignparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const std::string seqDb = par.db1;
    const std::string edgeDir = par.db2;
    const std::string alnDir = par.db3;

    const int dbType = FileUtil::parseDbType(seqDb.c_str());
    if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        // Not a gap in the port: stock's align2clust, which this stage is the
        // parallel form of, builds a protein matrix and protein Sequences
        // unconditionally (Align2clust.cpp:406-520). Nucleotide linclust runs the
        // v1 path (rescorediagonal) instead, which has no v2 counterpart to port.
        // Erroring is therefore stricter than stock, which would silently score
        // DNA with BLOSUM.
        Debug(Debug::ERROR) << "alignparallel is amino acid only, as stock align2clust is. "
                            << "Nucleotide linclust uses the v1 alignment path.\n";
        EXIT(EXIT_FAILURE);
    }
    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);

    const std::string manifestPath = edgeDir + "/coord/edge.info";
    if (FileUtil::fileExists(manifestPath.c_str()) == false) {
        Debug(Debug::ERROR) << "No edge buckets at " << edgeDir << " (" << manifestPath
                            << " is missing). Run kmerreduceparallel first.\n";
        EXIT(EXIT_FAILURE);
    }
    unsigned int bucketCount = 0;
    uint64_t bucketSpan = 0;
    unsigned int partitionCount = 0;
    unsigned int waveCount = 0;
    {
        FILE *file = FileUtil::openFileOrDie(manifestPath.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(file, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "bucketCount") bucketCount = static_cast<unsigned int>(value);
            else if (key == "bucketSpan") bucketSpan = value;
            else if (key == "partitionCount") partitionCount = static_cast<unsigned int>(value);
            else if (key == "waveCount") waveCount = static_cast<unsigned int>(value);
        }
        fclose(file);
    }
    if (bucketCount == 0) {
        Debug(Debug::ERROR) << "Edge manifest " << manifestPath << " has no bucket count\n";
        EXIT(EXIT_FAILURE);
    }
    if (partitionCount == 0 || waveCount == 0) {
        Debug(Debug::ERROR) << "Edge manifest " << manifestPath << " does not record how many "
                            << "partitions and waves the reduce covered, so this stage cannot "
                            << "tell a finished reduce from a partial one. It was written by an "
                            << "older build; re-run kmerreduceparallel with this one.\n";
        EXIT(EXIT_FAILURE);
    }
    par.printParameters(command.cmd, argc, argv, *command.params);

    if (FileUtil::directoryExists(alnDir.c_str()) == false) {
        FileUtil::makeDir(alnDir.c_str());
    }
    const std::string coordDir = alnDir + "/coord";
    if (FileUtil::directoryExists(coordDir.c_str()) == false) {
        FileUtil::makeDir(coordDir.c_str());
    }

    BaseMatrix *subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);
    const uint64_t residues = info.dataSize - 2 * info.entryCount;
    EvalueComputation evaluer(residues, subMat);
    const std::string library = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? getCovSeqidQscPercMinDiag()
                                    : getCovSeqidQscPercMinDiagTargetCov();
    const float scorePerColThreshold = parsePrecisionLib(library, par.seqIdThr, par.covThr, 0.99);
    unsigned int maxSeqLen = info.maxSeqLen + 1;
    if (par.filterSeqDBFile.empty() == false) {
        // Stock sizes its aligner buffers with the larger of the two databases
        // (Align2clust.cpp:502), because the filter gate compares against
        // sequences drawn from the full database, not the pass-2 one.
        const DenseIndex::Info full = DenseIndex::readInfo(par.filterSeqDBFile);
        maxSeqLen = std::max(maxSeqLen, full.maxSeqLen + 1);
    }
    // Same x-drop as stock (Align2clust.cpp:419, MIN_SIZE 32).
    const int32_t xDrop = 32 * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid();
    const unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);

    // Pass 2 only: stock's --filter-cludb-file gate.
    FilterGate *gate = NULL;
    if (par.filterCluDBFile.empty() == false) {
        if (par.filterSeqDBFile.empty() || par.keyMapFile.empty()) {
            Debug(Debug::ERROR) << "--filter-cludb-file needs --filter-seqdb-file and --key-map\n";
            EXIT(EXIT_FAILURE);
        }
        gate = new FilterGate(par.filterSeqDBFile);
        // seqDb is the re-keyed representative database this pass runs on, which
        // is where createrepdb put the CSR. --filter-cludb-file is no longer read
        // here at all: it is the input the CSR was built from, kept as a parameter
        // so the stage still states which clustering it is gating against.
        gate->open(par.keyMapFile, seqDb);
        Debug(Debug::INFO) << "Filter gate: " << gate->repCount << " representatives, "
                           << gate->memberCount << " pass-1 members, paged from "
                           << seqDb << ".gate.*\n";
    }
    Debug(Debug::INFO) << "Aligning " << bucketCount << " edge buckets; score-per-column cutoff "
                       << scorePerColThreshold << "\n";

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined\n";

    PartitionSequences sequences(seqDb);
    uint64_t survivorCount = 0;
    const std::vector<int64_t> reduceAuthority =
        readReduceAuthority(edgeDir, waveCount, partitionCount);

    // The layout the greedy sweep must consume, recorded next to the output rather
    // than rediscovered from it. greedycluster used to infer the bucket count by
    // counting consecutive p<n>.edges files, so an align that had only produced a
    // prefix read as a complete, smaller layout -- and it recomputed bucketSpan
    // from that smaller count, silently sweeping the wrong key ranges.
    {
        FileLock layoutLock(coordDir + "/align.lock");
        layoutLock.lock();
        const std::string layoutPath = coordDir + "/align.info";
        if (FileUtil::fileExists(layoutPath.c_str()) == false) {
            const std::string tmp = layoutPath + ".tmp." + SSTR(getpid());
            FILE *f = FileUtil::openAndDelete(tmp.c_str(), "w");
            fprintf(f, "bucketCount\t%zu\n", (size_t)bucketCount);
            fprintf(f, "bucketSpan\t%zu\n", (size_t)bucketSpan);
            fprintf(f, "entryCount\t%zu\n", (size_t)info.entryCount);
            if (fclose(f) != 0 || rename(tmp.c_str(), layoutPath.c_str()) != 0) {
                Debug(Debug::ERROR) << "Cannot publish " << layoutPath << ": " << strerror(errno)
                                    << "\n";
                layoutLock.unlock();
                EXIT(EXIT_FAILURE);
            }
        }
        layoutLock.unlock();
    }

    {
        WorkQueue queue(coordDir + "/align.queue", static_cast<int64_t>(bucketCount));
        const bool finished = queue.drain(workerId, [&](size_t bucket) {
            std::vector<CandidateEdge> edges;
            readBucket(edgeDir, static_cast<unsigned int>(bucket), reduceAuthority, edges);
            const size_t raw = edges.size();
            if (raw == 0) {
                EdgeWriter empty(EdgeWriter::partitionPath(alnDir, static_cast<unsigned int>(bucket)), workerId);
                empty.close();
                return;
            }

            const size_t merged = mergePairCopies(edges);

            // Both endpoints, in ascending key order so the fetch is a forward
            // scan. Representatives are contiguous within a bucket by
            // construction; members are close to them because keys are
            // length-ranked and homologues have similar lengths.
            std::vector<uint64_t> needed;
            needed.reserve(edges.size() * 2);
            for (size_t i = 0; i < edges.size(); i++) {
                needed.push_back(edges[i].getRep());
                needed.push_back(edges[i].getMember());
            }
            SORT_PARALLEL(needed.begin(), needed.end());
            needed.erase(std::unique(needed.begin(), needed.end()), needed.end());
            sequences.load(needed);

            // The gate compares against the pass-1 cluster members of each target,
            // which live in the *full* database; fetch exactly those, in key order.
            if (gate != NULL) {
                // The bucket's member sub-keys, deduplicated *before* their
                // clusters are expanded. Expanding first and deduplicating after
                // -- which is what this used to do -- materialised one key per
                // (edge, cluster member) pair, so a bucket referring to a few
                // large pass-1 clusters many times over held far more than the
                // clusters themselves.
                std::vector<uint64_t> wanted;
                wanted.reserve(edges.size());
                for (size_t i = 0; i < edges.size(); i++) {
                    const uint64_t sub = edges[i].getMember();
                    if (sub < gate->repCount) {
                        wanted.push_back(sub);
                    }
                }
                SORT_PARALLEL(wanted.begin(), wanted.end());
                wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());
                std::vector<uint64_t> gateKeys;
                gate->loadSlice(wanted, gateKeys);
                SORT_PARALLEL(gateKeys.begin(), gateKeys.end());
                gateKeys.erase(std::unique(gateKeys.begin(), gateKeys.end()), gateKeys.end());
                gate->fullSeqs.load(gateKeys);
            }

            // Edges are sorted by (rep, member), so a representative's edges are
            // contiguous and its query profile is built once.
            std::vector<size_t> repStarts;
            for (size_t i = 0; i < edges.size(); i++) {
                if (i == 0 || edges[i].getRep() != edges[i - 1].getRep()) {
                    repStarts.push_back(i);
                }
            }
            repStarts.push_back(edges.size());
            std::vector<unsigned char> survives(edges.size(), 0);

#pragma omp parallel num_threads(par.threads)
            {
                Sequence query(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false,
                               par.compBiasCorrection);
                Sequence target(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false,
                                par.compBiasCorrection);
                BlockAligner aligner(Parameters::DBTYPE_AMINO_ACIDS, maxSeqLen, subMat, &fastMatrix,
                                     &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale,
                                     -par.gapOpen.values.aminoacid(),
                                     -par.gapExtend.values.aminoacid());
                Matcher matcher(Parameters::DBTYPE_AMINO_ACIDS, maxSeqLen, subMat, &evaluer,
                                par.compBiasCorrection, par.compBiasCorrectionScale,
                                par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(),
                                0.0, par.zdrop);
                Sequence element(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, subMat, 0, false,
                                 par.compBiasCorrection);

#pragma omp for schedule(dynamic, 16)
                for (size_t g = 0; g < repStarts.size() - 1; g++) {
                    const size_t from = repStarts[g];
                    const size_t to = repStarts[g + 1];
                    unsigned int repLen = 0;
                    const char *repSeq = sequences.get(edges[from].getRep(), &repLen);
                    if (repSeq == NULL || repLen == 0) {
                        continue;
                    }
                    query.mapSequence(0, static_cast<DBKeyType>(edges[from].getRep()), repSeq, repLen);
                    aligner.initQuery(&query);
                    matcher.initQuery(&query);

                    for (size_t e = from; e < to; e++) {
                        unsigned int memberLen = 0;
                        const char *memberSeq = sequences.get(edges[e].getMember(), &memberLen);
                        if (memberSeq == NULL || memberLen == 0) {
                            continue;
                        }
                        if (Util::canBeCovered(par.covThr, par.covMode, query.L,
                                               static_cast<int>(memberLen)) == false) {
                            continue;
                        }
                        target.mapSequence(0, static_cast<DBKeyType>(edges[e].getMember()),
                                           memberSeq, memberLen);

                        // Stock's two-stage acceptance (Align2clust.cpp:660-792): an
                        // ungapped alignment on the k-mer diagonal, and if that
                        // fails, a gapped banded alignment seeded from a
                        // three-residue exact match. Leaving the second stage out
                        // silently drops every pair that needs gaps to align.
                        BlockAligner::UngappedAln_res aln =
                            aligner.ungappedAlign(&target, edges[e].diagonal);

                        const bool hasEvalue = (aln.eval <= par.evalThr);
                        const bool hasAlnLen = (aln.alnLen >= par.alnLenThr);
                        const bool hasCoverage =
                            Util::hasCoverage(par.covThr, par.covMode, aln.qcov, aln.tcov);
                        float seqId = 0;
                        if (hasEvalue) {
                            int identical = 0;
                            for (int q = aln.qStart; q <= aln.qEnd; q++) {
                                const char a = repSeq[q] & static_cast<unsigned char>(~0x20);
                                const char b = memberSeq[aln.tStart + (q - aln.qStart)] &
                                               static_cast<unsigned char>(~0x20);
                                identical += (a == b) ? 1 : 0;
                            }
                            seqId = Util::computeSeqId(par.seqIdMode, identical, query.L, target.L,
                                                       aln.alnLen);
                        }
                        const bool hasSeqId =
                            seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());

                        if (hasAlnLen && hasCoverage && hasSeqId && hasEvalue) {
                            // Stock seeds this gate with the pair's diagonal
                            // (Align2clust.cpp:680), not with 0 as the gapped
                            // branch below does.
                            if (passesFilterGate(gate, edges[e].getMember(), query, element,
                                                 aligner, matcher, par, swMode,
                                                 edges[e].diagonal) == false) {
                                continue;
                            }
                            // The greedy ranks by alignment score, not k-mer count.
                            edges[e].score =
                                static_cast<uint16_t>(std::min(aln.bitScore, 65535));
                            survives[e] = 1;
                            continue;
                        }

                        // Ungapped failed. score-per-column is the gate on paying
                        // for a gapped alignment -- it is *not* a filter on hits the
                        // ungapped stage already accepted.
                        if (aln.diagonalLen <= 0 ||
                            static_cast<float>(aln.score) / static_cast<float>(aln.diagonalLen) <
                                scorePerColThreshold) {
                            continue;
                        }
                        if (aln.qStart == -1 || aln.tStart == -1 || aln.alnLen < 3) {
                            continue;
                        }
                        // Seed the band on the first three consecutive identities,
                        // starting one past it, exactly as stock does.
                        int seedQuery = static_cast<int>(aln.qStart);
                        int seedTarget = static_cast<int>(aln.tStart);
                        bool foundSeed = false;
                        for (int b = 0; b <= aln.alnLen - 3; b++) {
                            const int qp = static_cast<int>(aln.qStart) + b;
                            const int tp = static_cast<int>(aln.tStart) + b;
                            if (repSeq[qp] == memberSeq[tp] && repSeq[qp + 1] == memberSeq[tp + 1] &&
                                repSeq[qp + 2] == memberSeq[tp + 2]) {
                                seedQuery = qp + 1;
                                seedTarget = tp + 1;
                                foundSeed = true;
                                break;
                            }
                        }
                        if (foundSeed == false) {
                            continue;
                        }

                        std::string backtrace;
                        s_align gapped = aligner.bandedalign(&target, seedQuery, seedTarget,
                                                             backtrace, xDrop, par.covThr,
                                                             par.covMode);
                        const unsigned int gappedLen = backtrace.size();
                        const double gappedSeqId =
                            Util::computeSeqId(par.seqIdMode, gapped.identicalAACnt, query.L,
                                               memberLen, gappedLen);
                        Matcher::result_t result(
                            static_cast<DBKeyType>(edges[e].getMember()), gapped.score1,
                            gapped.qCov, gapped.tCov, gappedSeqId, gapped.evalue, gappedLen,
                            gapped.qStartPos1, gapped.qEndPos1, query.L, gapped.dbStartPos1,
                            gapped.dbEndPos1, memberLen, backtrace);
                        // isIdentity is always false here: self-edges are dropped in
                        // the reduce, so a representative is never its own member.
                        if (Alignment::checkCriteria(result, false, par.evalThr, par.seqIdThr,
                                                     par.alnLenThr, par.covMode, par.covThr)) {
                            // 0 here, matching Align2clust.cpp:816.
                            if (passesFilterGate(gate, edges[e].getMember(), query, element,
                                                 aligner, matcher, par, swMode, 0) == false) {
                                continue;
                            }
                            edges[e].score =
                                static_cast<uint16_t>(std::min<uint32_t>(gapped.score1, 65535));
                            survives[e] = 1;
                        }
                    }
                }
            }

            EdgeWriter writer(EdgeWriter::partitionPath(alnDir, static_cast<unsigned int>(bucket)), workerId);
            size_t kept = 0;
            for (size_t i = 0; i < edges.size(); i++) {
                if (survives[i]) {
                    writer.append(edges[i]);
                    kept++;
                }
            }
            writer.close();
            __sync_fetch_and_add(&survivorCount, static_cast<uint64_t>(kept));
            Debug(Debug::INFO) << "Bucket " << bucket << ": " << raw << " copies -> " << merged
                               << " pairs -> " << kept << " surviving, " << needed.size()
                               << " sequences (" << sequences.getBytes() / (1024 * 1024)
                               << " MB arena, " << sequences.getBytesRead() / (1024 * 1024)
                               << " MB read)\n";
        });
        if (finished == false) {
            Debug(Debug::ERROR) << "Align stage stalled: work remains but no bucket is claimable\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO) << "Worker " << workerId << " wrote " << survivorCount << " surviving edges\n";
    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
    delete subMat;
    return EXIT_SUCCESS;
}
