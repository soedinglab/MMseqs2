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
#include <cstring>
#include <limits>
#include <string>
#include <vector>

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

// Which worker's edges count for each k-mer partition.
//
// The reduce's work queue is the authority: complete() keeps the first worker to
// record an item done, so exactly one is named per partition, and any block from
// another worker is a copy left by one that died before recording it. One queue
// per extraction wave, and wave w holds the partitions after all earlier waves',
// so appending them in wave order indexes by partition directly.
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
std::vector<int64_t> readReduceAuthority(const std::string &edgeDir) {
    std::vector<int64_t> authority;
    for (unsigned int wave = 0;; wave++) {
        std::vector<int64_t> workers;
        const std::string path = edgeDir + "/coord/reduce." + SSTR(wave) + ".queue";
        if (WorkQueue::readCompletedWorkers(path, workers) == false) {
            break;
        }
        authority.insert(authority.end(), workers.begin(), workers.end());
    }
    if (authority.empty()) {
        Debug(Debug::WARNING) << "No reduce work queue under " << edgeDir
                              << "/coord; cannot tell a crashed worker's superseded edges from a "
                              << "second partition's. Duplicate support would be summed.\n";
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
// lookup a CSR index over sub-keys rather than anything indexed by the full key
// space: an offset per representative plus one key per sequence.
struct FilterGate {
    std::vector<uint64_t> keymap;        // sub-key -> original key
    std::vector<uint64_t> clusterStart;  // sub-key -> offset into members
    std::vector<uint64_t> members;       // original keys, grouped by pass-1 cluster
    PartitionSequences fullSeqs;

    FilterGate(const std::string &fullDb) : fullSeqs(fullDb) {}

    size_t size(uint64_t sub) const { return clusterStart[sub + 1] - clusterStart[sub]; }

    void load(const std::string &keymapFile, const std::string &pass1Tsv) {
        const size_t bytes = FileUtil::getFileSize(keymapFile);
        keymap.resize(bytes / sizeof(uint64_t));
        FILE *m = FileUtil::openFileOrDie(keymapFile.c_str(), "rb", true);
        if (fread(keymap.data(), sizeof(uint64_t), keymap.size(), m) != keymap.size()) {
            Debug(Debug::ERROR) << "Cannot read " << keymapFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(m);

        // Two streaming passes: count, then fill. The representative of a pass-1
        // cluster is looked up in the ascending key map by binary search.
        std::vector<uint64_t> counts(keymap.size() + 1, 0);
        for (int pass = 0; pass < 2; pass++) {
            FILE *f = FileUtil::openFileOrDie(pass1Tsv.c_str(), "r", true);
            char *line = NULL;
            size_t cap = 0;
            while (getline(&line, &cap, f) > 0) {
                char *tab = strchr(line, '\t');
                if (tab == NULL) continue;
                const uint64_t rep = strtoull(line, NULL, 10);
                const uint64_t member = strtoull(tab + 1, NULL, 10);
                const std::vector<uint64_t>::const_iterator it =
                    std::lower_bound(keymap.begin(), keymap.end(), rep);
                if (it == keymap.end() || *it != rep) continue;
                const size_t sub = static_cast<size_t>(it - keymap.begin());
                if (pass == 0) {
                    counts[sub]++;
                } else {
                    members[clusterStart[sub] + counts[sub]] = member;
                    counts[sub]++;
                }
            }
            free(line);
            fclose(f);
            if (pass == 0) {
                clusterStart.assign(keymap.size() + 1, 0);
                uint64_t total = 0;
                for (size_t i = 0; i < keymap.size(); i++) {
                    clusterStart[i] = total;
                    total += counts[i];
                }
                clusterStart[keymap.size()] = total;
                members.assign(total, 0);
                counts.assign(keymap.size() + 1, 0);
            }
        }
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
    if (gate == NULL || subMember >= gate->keymap.size()) {
        return true;
    }
    if (gate->size(subMember) <= 1) {
        return true;  // stock only runs the loop when numClu > 1
    }
    const uint64_t targetKey = gate->keymap[subMember];
    for (uint64_t j = gate->clusterStart[subMember]; j < gate->clusterStart[subMember + 1]; j++) {
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
    {
        FILE *file = FileUtil::openFileOrDie(manifestPath.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(file, "%63s\t%zu\n", name, &value) == 2) {
            if (std::string(name) == "bucketCount") {
                bucketCount = static_cast<unsigned int>(value);
            }
        }
        fclose(file);
    }
    if (bucketCount == 0) {
        Debug(Debug::ERROR) << "Edge manifest " << manifestPath << " has no bucket count\n";
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
        gate->load(par.keyMapFile, par.filterCluDBFile);
        Debug(Debug::INFO) << "Filter gate: " << gate->keymap.size() << " representatives, "
                           << gate->members.size() << " pass-1 members\n";
    }
    Debug(Debug::INFO) << "Aligning " << bucketCount << " edge buckets; score-per-column cutoff "
                       << scorePerColThreshold << "\n";

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined\n";

    PartitionSequences sequences(seqDb);
    uint64_t survivorCount = 0;
    const std::vector<int64_t> reduceAuthority = readReduceAuthority(edgeDir);

    {
        WorkQueue queue(coordDir + "/align.queue", static_cast<int64_t>(bucketCount));
        const bool finished = queue.drain(workerId, [&](size_t bucket) {
            std::vector<CandidateEdge> edges;
            readBucket(edgeDir, static_cast<unsigned int>(bucket), reduceAuthority, edges);
            const size_t raw = edges.size();
            if (raw == 0) {
                EdgeWriter empty(EdgeWriter::partitionPath(alnDir, static_cast<unsigned int>(bucket)));
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
                std::vector<uint64_t> gateKeys;
                for (size_t i = 0; i < edges.size(); i++) {
                    const uint64_t sub = edges[i].getMember();
                    if (sub >= gate->keymap.size()) continue;
                    for (uint64_t j = gate->clusterStart[sub]; j < gate->clusterStart[sub + 1]; j++) {
                        gateKeys.push_back(gate->members[j]);
                    }
                }
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

            EdgeWriter writer(EdgeWriter::partitionPath(alnDir, static_cast<unsigned int>(bucket)));
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
