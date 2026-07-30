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
        edges[out].score = static_cast<uint8_t>(std::min(bestTotal, 255));
        out++;
        i = j;
    }
    edges.resize(out);
    return out;
}

size_t readBucket(const std::string &edgeDir, unsigned int bucket,
                  std::vector<CandidateEdge> &out) {
    const std::vector<std::string> shards = EdgeBucketWriter::shardFiles(edgeDir, bucket);
    for (size_t i = 0; i < shards.size(); i++) {
        const size_t bytes = FileUtil::getFileSize(shards[i]);
        if (bytes % sizeof(CandidateEdge) != 0) {
            Debug(Debug::ERROR) << "Edge shard " << shards[i] << " is " << bytes
                                << " bytes, not a whole number of edges. It was probably written "
                                << "by an interrupted worker.\n";
            EXIT(EXIT_FAILURE);
        }
        const size_t count = bytes / sizeof(CandidateEdge);
        if (count == 0) {
            continue;
        }
        FILE *file = FileUtil::openFileOrDie(shards[i].c_str(), "rb", true);
        const size_t offset = out.size();
        out.resize(offset + count);
        if (fread(out.data() + offset, sizeof(CandidateEdge), count, file) != count) {
            Debug(Debug::ERROR) << "Cannot read edge shard " << shards[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(file);
    }
    return out.size();
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
        Debug(Debug::ERROR) << "alignparallel is implemented for amino acid databases only\n";
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
    const unsigned int maxSeqLen = info.maxSeqLen + 1;
    // Same x-drop as stock (Align2clust.cpp:419, MIN_SIZE 32).
    const int32_t xDrop = 32 * par.gapExtend.values.aminoacid() + par.gapOpen.values.aminoacid();
    Debug(Debug::INFO) << "Aligning " << bucketCount << " edge buckets; score-per-column cutoff "
                       << scorePerColThreshold << "\n";

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined\n";

    PartitionSequences sequences(seqDb);
    uint64_t survivorCount = 0;

    {
        WorkQueue queue(coordDir + "/align.queue", static_cast<int64_t>(bucketCount));
        const bool finished = queue.drain(workerId, [&](size_t bucket) {
            std::vector<CandidateEdge> edges;
            readBucket(edgeDir, static_cast<unsigned int>(bucket), edges);
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
                            // The greedy ranks by alignment score, not k-mer count.
                            edges[e].score =
                                static_cast<uint8_t>(std::min(aln.bitScore, 255));
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
                            edges[e].score =
                                static_cast<uint8_t>(std::min<uint32_t>(gapped.score1, 255));
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
