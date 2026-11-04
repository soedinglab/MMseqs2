#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "FastSort.h"
#include "EMCoverage.h"
#include "EMTargetGroups.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif

namespace {

struct HitData {
    // per query: the reads with >= 1 hit.  The filter rounds compact deleted reads away, so after
    // filtering these are the surviving reads only.
    std::vector<unsigned int> queryKey;
    std::vector<size_t> queryOffset;          // queryCount + 1 entries

    // per hit, hot
    std::vector<unsigned int> hitTarget;      // dense index into targetKey
    std::vector<double> hitPosterior;         // EM posterior gamma_qt, from alignment column 3

    // per hit, cold: read by the coverage pass, which needs the CIGAR and the target span
    std::vector<Matcher::result_t> hitResult;

    // per target
    std::vector<unsigned int> targetKey;      // sorted, unique
    std::vector<unsigned int> targetLength;   // dbLen
    std::vector<double> abundance;            // theta_t
    std::vector<size_t> topHitReads;          // reads whose top hit (highest posterior) is the target
    std::vector<double> expectedReads;        // sum of the posteriors on the target
    std::vector<char> dropped;                // removed by --min-reads
    std::vector<int> droppedRound;            // filter round that removed it (0 = kept)

    size_t queryCount() const { return queryOffset.empty() ? 0 : (queryOffset.size() - 1); }
    size_t hitCount() const { return hitTarget.size(); }
    size_t targetCount() const { return targetKey.size(); }
};

struct TargetStats {
    unsigned int key;
    unsigned int targetLength;
    double abundance;
    size_t topHitReads;
    double expectedReads;
    bool dropped;
    int droppedRound;
    // Posterior-weighted coverage, using the EM assignment probabilities gamma_qt.
    //   coverageDepth     = mean expected read depth over the target = (1/L) * sum_p depth_p
    //   coveredLength     = expected number of covered residues = sum_p [1 - prod_hits(1 - gamma)]
    //   breadth           = coveredLength / L  (derived at write time)
    //   coverageEvenness  = how uniformly that depth is spread over the target (EMCoverage::evenness)
    double coverageDepth;
    double coveredLength;
    double coverageEvenness;
};


struct ReclassTaxContext {
    HitData hits;
    size_t hitCount;            // hits loaded from the input (before filtering)
    size_t backtraceHitCount;   // hits whose coverage came from the CIGAR
    size_t reverseStrandHits;   // hits with dbStartPos > dbEndPos (reverse-strand target)
    bool hasBacktrace;          // input records carry a backtrace column
    bool hasOrfPosition;        // input records carry the ORF position columns
    std::vector<unsigned int> emptyQueryKey;    // queries of the input DB without any hit
    std::vector<unsigned int> deletedQueryKey;  // reads deleted by the filter (their top hit was removed)

    ReclassTaxContext() : hitCount(0), backtraceHitCount(0), reverseStrandHits(0), hasBacktrace(false),
                          hasOrfPosition(false) {}
};

static double clamp01(double value);

static const size_t SCALAR_REDUCTION_BLOCKS = 64;

static int usableThreads(int threads) {
#ifdef OPENMP
    return std::max(1, threads);
#else
    (void)threads;
    return 1;
#endif
}


struct TargetHitIndex {
    std::vector<size_t> offset;
    std::vector<unsigned int> hit;
};

static TargetHitIndex groupHitsByTarget(const HitData &hits) {
    const size_t targetCount = hits.targetCount();
    const size_t hitCount = hits.hitCount();
    if (hitCount > static_cast<size_t>(std::numeric_limits<unsigned int>::max())) {
        Debug(Debug::ERROR) << "Cannot index " << hitCount
                            << " hits for the coverage pass; split the alignment DB.\n";
        EXIT(EXIT_FAILURE);
    }

    TargetHitIndex byTarget;
    byTarget.offset.assign(targetCount + 1, 0);
    for (size_t h = 0; h < hitCount; ++h) {
        ++byTarget.offset[hits.hitTarget[h] + 1];
    }
    for (size_t t = 0; t < targetCount; ++t) {
        byTarget.offset[t + 1] += byTarget.offset[t];
    }
    byTarget.hit.assign(hitCount, 0);
    std::vector<size_t> cursor(byTarget.offset.begin(), byTarget.offset.end() - 1);
    for (size_t h = 0; h < hitCount; ++h) {
        byTarget.hit[cursor[hits.hitTarget[h]]++] = static_cast<unsigned int>(h);
    }
    return byTarget;
}

static size_t countRecords(const char *data) {
    size_t records = 0;
    char *cursor = const_cast<char *>(data);
    while (*cursor != '\0') {
        ++records;
        cursor = Util::skipLine(cursor);
    }
    return records;
}

static void loadAlignmentDb(DBReader<unsigned int> &reader, ReclassTaxContext &ctx, int threads) {
    HitData &hits = ctx.hits;
    const int nThreads = usableThreads(threads);

    std::vector<size_t> entryIds;
    entryIds.reserve(reader.getSize());
    for (size_t i = 0; i < reader.getSize(); ++i) {
        if (reader.getEntryLen(i) <= 1) {
            ctx.emptyQueryKey.push_back(reader.getDbKey(i));
            continue;
        }
        entryIds.push_back(i);
    }
    std::vector<size_t> recordsPerQuery(entryIds.size(), 0);

    Debug::Progress countProgress(entryIds.size());
#pragma omp parallel num_threads(nThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 64)
        for (size_t q = 0; q < entryIds.size(); ++q) {
            countProgress.updateProgress();
            recordsPerQuery[q] = countRecords(reader.getData(entryIds[q], thread_idx));
        }
    }

    hits.queryKey.resize(entryIds.size());
    hits.queryOffset.resize(entryIds.size() + 1, 0);
    for (size_t q = 0; q < entryIds.size(); ++q) {
        hits.queryKey[q] = reader.getDbKey(entryIds[q]);
        hits.queryOffset[q + 1] = hits.queryOffset[q] + recordsPerQuery[q];
    }
    const size_t hitCount = hits.queryOffset.back();
    hits.hitTarget.resize(hitCount);
    hits.hitPosterior.resize(hitCount);
    hits.hitResult.resize(hitCount);

    std::vector<size_t> backtraceHits(nThreads, 0);
    std::vector<size_t> reverseHits(nThreads, 0);
    std::vector<char> sawBacktrace(nThreads, 0);
    std::vector<char> sawOrfPosition(nThreads, 0);
    // dbKey -> dbLen over every target seen. The key set is the target set.
    std::vector<std::unordered_map<unsigned int, unsigned int> > targetLengthByThread(nThreads);

    Debug::Progress progress(entryIds.size());
#pragma omp parallel num_threads(nThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const char *entry[255];
        std::unordered_map<unsigned int, unsigned int> &targetLength = targetLengthByThread[thread_idx];

#pragma omp for schedule(dynamic, 64)
        for (size_t q = 0; q < entryIds.size(); ++q) {
            progress.updateProgress();
            char *data = reader.getData(entryIds[q], thread_idx);
            size_t h = hits.queryOffset[q];
            while (*data != '\0') {
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Debug(Debug::ERROR) << "Invalid alignment result record in query "
                                        << reader.getDbKey(entryIds[q]) << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT
                    || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                    sawBacktrace[thread_idx] = 1;
                }
                if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT
                    || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                    sawOrfPosition[thread_idx] = 1;
                }

                Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);

                // Always use the 3rd alignment column (seqId) as posterior.
                // `mmseqs em` writes that column with resultToBuffer(..., highPrecisionSeqId=true),
                hits.hitPosterior[h] = clamp01(static_cast<double>(result.seqId));
                hits.hitTarget[h] = result.dbKey;  // raw dbKey until the remap below
                targetLength[result.dbKey] = result.dbLen;

                if (result.backtrace.empty() == false) {
                    ++backtraceHits[thread_idx];
                }
                if (result.dbStartPos > result.dbEndPos) {
                    ++reverseHits[thread_idx];
                }

                hits.hitResult[h] = std::move(result);
                ++h;
                data = Util::skipLine(data);
            }
        }
    }

    ctx.hitCount = hitCount;
    for (int t = 0; t < nThreads; ++t) {
        ctx.backtraceHitCount += backtraceHits[t];
        ctx.reverseStrandHits += reverseHits[t];
        ctx.hasBacktrace = ctx.hasBacktrace || (sawBacktrace[t] != 0);
        ctx.hasOrfPosition = ctx.hasOrfPosition || (sawOrfPosition[t] != 0);
    }

    std::unordered_map<unsigned int, unsigned int> targetLengthByKey;
    for (int t = 0; t < nThreads; ++t) {
        targetLengthByKey.insert(targetLengthByThread[t].begin(), targetLengthByThread[t].end());
        std::unordered_map<unsigned int, unsigned int>().swap(targetLengthByThread[t]);
    }


    hits.targetKey.reserve(targetLengthByKey.size());
    for (std::unordered_map<unsigned int, unsigned int>::const_iterator it = targetLengthByKey.begin();
         it != targetLengthByKey.end(); ++it) {
        hits.targetKey.push_back(it->first);
    }
    SORT_SERIAL(hits.targetKey.begin(), hits.targetKey.end());
    hits.targetLength.resize(hits.targetKey.size());
    for (size_t t = 0; t < hits.targetKey.size(); ++t) {
        hits.targetLength[t] = targetLengthByKey[hits.targetKey[t]];
    }
    const size_t targetCount = hits.targetKey.size();
    hits.abundance.assign(targetCount, 0.0);
    hits.topHitReads.assign(targetCount, 0);
    hits.expectedReads.assign(targetCount, 0.0);
    hits.dropped.assign(targetCount, 0);
    hits.droppedRound.assign(targetCount, 0);

#pragma omp parallel for schedule(static) num_threads(nThreads)
    for (size_t h = 0; h < hitCount; ++h) {
        const std::vector<unsigned int>::const_iterator found =
                std::lower_bound(hits.targetKey.begin(), hits.targetKey.end(), hits.hitTarget[h]);
        hits.hitTarget[h] = static_cast<unsigned int>(found - hits.targetKey.begin());
    }
}


static const char *WARNING_RULE =
        "========================================================================\n";


static void warningBanner(const std::string &headline, const std::string &detail) {
    Debug(Debug::WARNING) << "\n" << WARNING_RULE
                          << "  WARNING: " << headline << "\n"
                          << detail
                          << WARNING_RULE << "\n";
}

static void reportInputDiagnostics(const ReclassTaxContext &ctx) {
    const HitData &hits = ctx.hits;
    Debug(Debug::INFO) << ctx.hitCount << " hits; coverage source: "
                       << ctx.backtraceHitCount << " from backtrace, "
                       << (ctx.hitCount - ctx.backtraceHitCount) << " from aligned span";
    if (ctx.reverseStrandHits > 0) {
        Debug(Debug::INFO) << "; " << ctx.reverseStrandHits << " reverse-strand";
    }
    Debug(Debug::INFO) << ".\n";
    if (ctx.backtraceHitCount < ctx.hitCount) {
        warningBanner("Some hits carry no backtrace.",
                      "           You should run 'mmseqs search' with the '-a' parameter.\n"
                      "\n"
                      "           " + SSTR(ctx.hitCount - ctx.backtraceHitCount) + " of "
                      + SSTR(ctx.hitCount) + " hits have no backtrace, so their coverage is\n"
                      "           measured over the whole aligned span instead of the CIGAR.\n"
                      "           Indels are ignored and coverage is over-counted by the gap\n"
                      "           length, so coverage_depth and breadth_of_coverage are too high.\n");
    }

    const size_t queryCount = hits.queryCount();
    if (queryCount == 0) {
        return;
    }
    double deviationSum = 0.0;
    double maxDeviation = 0.0;
    size_t offQueries = 0;
    for (size_t q = 0; q < queryCount; ++q) {
        double posteriorSum = 0.0;
        for (size_t h = hits.queryOffset[q]; h < hits.queryOffset[q + 1]; ++h) {
            posteriorSum += hits.hitPosterior[h];
        }
        const double deviation = std::fabs(posteriorSum - 1.0);
        deviationSum += deviation;
        maxDeviation = std::max(maxDeviation, deviation);
        if (deviation > 1e-3) {
            ++offQueries;
        }
    }
    const double meanDeviation = deviationSum / static_cast<double>(queryCount);
    Debug(Debug::INFO) << "Posterior check: mean |sum(gamma) - 1| = " << meanDeviation
                       << ", max = " << maxDeviation << ", "
                       << offQueries << "/" << queryCount << " queries off by > 1e-3.\n";
    if (meanDeviation > 0.01) {
        warningBanner("You should run 'mmseqs em' first.",
                      "\n"
                      "           Alignment column 3 does not sum to 1 per query, so it does not\n"
                      "           hold EM posteriors -- it still holds sequence identities.\n"
                      "           mean |sum(gamma) - 1| = " + SSTR(meanDeviation)
                      + " over " + SSTR(queryCount) + " queries.\n"
                      "\n"
                      "           Every number below is meaningless until this DB has been\n"
                      "           through 'mmseqs em'.\n");
    }
}


// ------------------------------------------------------------------ target filter
// The top hit of a read: highest posterior, then highest score, then Matcher::compareHits.
static size_t topHitOf(const HitData &hits, size_t q) {
    const size_t begin = hits.queryOffset[q];
    const size_t end = hits.queryOffset[q + 1];
    size_t best = begin;
    for (size_t h = begin + 1; h < end; ++h) {
        if (hits.hitPosterior[h] > hits.hitPosterior[best]) {
            best = h;
        } else if (hits.hitPosterior[h] == hits.hitPosterior[best]) {
            if (hits.hitResult[h].score > hits.hitResult[best].score) {
                best = h;
            } else if (hits.hitResult[h].score == hits.hitResult[best].score
                       && Matcher::compareHits(hits.hitResult[h], hits.hitResult[best])) {
                best = h;
            }
        }
    }
    return best;
}

// Per target: the number of reads whose top hit it is, and the expected read count (sum of gamma).
static void computeSupport(HitData &hits, std::vector<size_t> &top, int threads) {
    (void)threads;
    const size_t queryCount = hits.queryCount();
    const size_t targetCount = hits.targetCount();
    top.resize(queryCount);
#pragma omp parallel for schedule(static) num_threads(usableThreads(threads))
    for (size_t q = 0; q < queryCount; ++q) {
        top[q] = topHitOf(hits, q);
    }
    hits.topHitReads.assign(targetCount, 0);
    hits.expectedReads.assign(targetCount, 0.0);
    for (size_t q = 0; q < queryCount; ++q) {
        hits.topHitReads[hits.hitTarget[top[q]]] += 1;
        for (size_t h = hits.queryOffset[q]; h < hits.queryOffset[q + 1]; ++h) {
            hits.expectedReads[hits.hitTarget[h]] += hits.hitPosterior[h];
        }
    }
}

// --min-reads.  A target that is the top hit of fewer than minReads reads is removed together with
// those reads (deleted, never re-assigned to their second-best target); the hits of the surviving
// reads on removed targets are dropped and their posteriors renormalised.  Removing reads changes
// the support of the other targets (their secondary mass is gone and the read pool shrank), so the
// pass is repeated until every remaining target is the top hit of at least minReads reads.
static void filterTargets(ReclassTaxContext &ctx, int minReads, int threads) {
    HitData &hits = ctx.hits;
    const size_t targetCount = hits.targetCount();
    const size_t initialQueries = hits.queryCount();
    const size_t initialHits = hits.hitCount();
    std::vector<size_t> top;
    std::vector<size_t> savedTop(targetCount, 0);      // support of a removed target at removal time
    std::vector<double> savedExp(targetCount, 0.0);
    size_t removedTargets = 0;
    int round = 0;
    const int maxRounds = 1000;
    while (true) {
        computeSupport(hits, top, threads);
        const size_t queryCount = hits.queryCount();
        std::vector<char> dropNow(targetCount, 0);
        size_t byReads = 0;
        for (size_t t = 0; t < targetCount; ++t) {
            if (hits.dropped[t] || minReads <= 0) {
                continue;
            }
            if (hits.topHitReads[t] < static_cast<size_t>(minReads)) {
                dropNow[t] = 1;
                hits.dropped[t] = 1;
                hits.droppedRound[t] = round + 1;
                savedTop[t] = hits.topHitReads[t];
                savedExp[t] = hits.expectedReads[t];
                ++byReads;
            }
        }
        if (byReads == 0 || round >= maxRounds) {
            break;
        }
        ++round;
        removedTargets += byReads;

        // in-place compaction: reads whose top hit was removed go, hits on removed targets go
        size_t writeHit = 0;
        size_t writeQuery = 0;
        size_t deletedReads = 0;
        std::vector<size_t> newOffset;
        newOffset.reserve(queryCount + 1);
        newOffset.push_back(0);
        for (size_t q = 0; q < queryCount; ++q) {
            const size_t begin = hits.queryOffset[q];
            const size_t end = hits.queryOffset[q + 1];
            if (dropNow[hits.hitTarget[top[q]]]) {
                ctx.deletedQueryKey.push_back(hits.queryKey[q]);
                ++deletedReads;
                continue;
            }
            const size_t first = writeHit;
            double mass = 0.0;
            for (size_t h = begin; h < end; ++h) {
                if (hits.dropped[hits.hitTarget[h]]) {
                    continue;
                }
                if (writeHit != h) {
                    hits.hitTarget[writeHit] = hits.hitTarget[h];
                    hits.hitPosterior[writeHit] = hits.hitPosterior[h];
                    hits.hitResult[writeHit] = std::move(hits.hitResult[h]);
                }
                mass += hits.hitPosterior[writeHit];
                ++writeHit;
            }
            if (mass > 0.0) {
                for (size_t h = first; h < writeHit; ++h) {
                    hits.hitPosterior[h] /= mass;
                }
            }
            hits.queryKey[writeQuery] = hits.queryKey[q];
            ++writeQuery;
            newOffset.push_back(writeHit);
        }
        hits.queryKey.resize(writeQuery);
        hits.queryOffset.swap(newOffset);
        hits.hitTarget.resize(writeHit);
        hits.hitPosterior.resize(writeHit);
        hits.hitResult.resize(writeHit);

        Debug(Debug::INFO) << "Filter round " << round << ": removed " << byReads
                           << " targets with < " << minReads << " top-hit reads and the "
                           << deletedReads << " reads whose top hit they were; "
                           << hits.queryCount() << " reads and " << hits.hitCount() << " hits remain.\n";
    }
    if (round >= maxRounds) {
        Debug(Debug::WARNING) << "Target filter stopped after " << maxRounds << " rounds.\n";
    }
    for (size_t t = 0; t < targetCount; ++t) {
        if (hits.dropped[t]) {
            hits.topHitReads[t] = savedTop[t];
            hits.expectedReads[t] = savedExp[t];
        }
    }
    const size_t keptTargets = targetCount - removedTargets;
    Debug(Debug::INFO) << "--min-reads " << minReads << ": " << removedTargets << " of " << targetCount
                       << " targets removed in " << round << " round(s), " << keptTargets << " kept; "
                       << ctx.deletedQueryKey.size() << " of " << initialQueries << " reads deleted, "
                       << (initialHits - hits.hitCount()) << " of " << initialHits << " hits dropped.\n";
    if (keptTargets == 0 && targetCount > 0) {
        warningBanner("The filter removed every target.",
                      "           --min-reads " + SSTR(minReads) + " is above the top-hit read count of\n"
                      "           every target, so nothing is left to report.\n");
    }
}


static void computeAbundanceFromPosterior(HitData &hits, const TargetHitIndex &byTarget, int threads) {
    (void)threads;  // only read by the num_threads() clause, which vanishes without OPENMP
    const size_t targetCount = hits.targetCount();
    const size_t queryCount = hits.queryCount();
    const double denom = (queryCount > 0) ? static_cast<double>(queryCount) : 1.0;
    hits.abundance.assign(targetCount, 0.0);

#pragma omp parallel for schedule(dynamic, 256) num_threads(usableThreads(threads))
    for (size_t t = 0; t < targetCount; ++t) {
        double mass = 0.0;
        for (size_t g = byTarget.offset[t]; g < byTarget.offset[t + 1]; ++g) {
            mass += hits.hitPosterior[byTarget.hit[g]];
        }
        hits.abundance[t] = mass / denom;
    }
}

static std::vector<TargetStats> collectTargetStats(const ReclassTaxContext &ctx,
                                                   const TargetHitIndex &byTarget,
                                                   int threads) {
    const HitData &hits = ctx.hits;
    const size_t targetCount = hits.targetCount();
    const int nThreads = usableThreads(threads);
    (void)nThreads;

    // Posterior-weighted coverage per target, on a map sized to the target length. Every position is
    // confined to [0, targetLength), so coveredLength <= targetLength by construction (no clamping).
    std::vector<TargetStats> out(targetCount);
    Debug::Progress covProgress(targetCount);
#pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (size_t t = 0; t < targetCount; ++t) {
        covProgress.updateProgress();
        TargetStats &stats = out[t];
        stats.key = hits.targetKey[t];
        stats.targetLength = hits.targetLength[t];
        stats.abundance = hits.abundance[t];
        stats.topHitReads = hits.topHitReads[t];
        stats.expectedReads = hits.expectedReads[t];
        stats.dropped = (hits.dropped[t] != 0);
        stats.droppedRound = hits.droppedRound[t];
        stats.coverageDepth = 0.0;
        stats.coveredLength = 0.0;
        stats.coverageEvenness = 0.0;
        if (stats.targetLength == 0 || stats.dropped) {
            continue;
        }
        // notCov[p] = prod over hits covering p of (1 - gamma); starts at 1 (uncovered).
        // depth[p]  = sum over hits covering p of gamma; the expected read depth at p.
        std::vector<double> notCov(stats.targetLength, 1.0);
        std::vector<double> depth(stats.targetLength, 0.0);
        for (size_t g = byTarget.offset[t]; g < byTarget.offset[t + 1]; ++g) {
            const size_t h = byTarget.hit[g];
            const Matcher::result_t &result = hits.hitResult[h];
            if (result.backtrace.empty() == false) {
                EMCoverage::paintBacktrace(notCov, depth, result.backtrace,
                                           result.dbStartPos, result.dbEndPos,
                                           hits.hitPosterior[h]);
            } else {
                EMCoverage::paintSpan(notCov, depth,
                                      result.dbStartPos, result.dbEndPos,
                                      hits.hitPosterior[h]);
            }
        }
        stats.coveredLength = EMCoverage::expectedCoveredLength(notCov);
        stats.coverageDepth = EMCoverage::meanDepth(depth);
        stats.coverageEvenness = EMCoverage::evenness(depth);
    }

    // kept targets by abundance, then the removed ones by their support at removal time
    std::sort(out.begin(), out.end(), [](const TargetStats &lhs, const TargetStats &rhs) {
        if (lhs.dropped != rhs.dropped) {
            return lhs.dropped < rhs.dropped;
        }
        if (lhs.abundance != rhs.abundance) {
            return lhs.abundance > rhs.abundance;
        }
        if (lhs.topHitReads != rhs.topHitReads) {
            return lhs.topHitReads > rhs.topHitReads;
        }
        return lhs.key < rhs.key;
    });
    return out;
}

static double clamp01(double value) {
    return std::max(0.0, std::min(1.0, value));
}

static void convertAbundanceToPercent(std::vector<TargetStats> &stats) {
    double total = 0.0;
    for (size_t i = 0; i < stats.size(); ++i) {
        if (stats[i].dropped == false) {
            total += stats[i].abundance;
        }
    }

    for (size_t i = 0; i < stats.size(); ++i) {
        if (stats[i].dropped || total <= 0.0) {
            stats[i].abundance = 0.0;
        } else {
            stats[i].abundance = 100.0 * (stats[i].abundance / total);
        }
    }
}

static const char *headerForKey(DBReader<unsigned int> &headerReader, unsigned int key, unsigned int threadIdx) {
    size_t id = headerReader.getId(key);
    if (id == UINT_MAX) {
        return NULL;
    }
    return headerReader.getData(id, threadIdx);
}

static std::string identifierForKey(DBReader<unsigned int> &headerReader, unsigned int key, unsigned int threadIdx) {
    const char *header = headerForKey(headerReader, key, threadIdx);
    if (header == NULL) {
        return SSTR(key);
    }
    std::string parsed = Util::parseFastaHeader(header);
    return parsed.empty() ? SSTR(key) : parsed;
}

// abundance_pct of the kept targets sums to 100 over the surviving reads.  A removed target is
// listed with Drop=y, abundance_pct 0 and the support (reads_top_hit, reads_expected) it had in the
// round that removed it, so the reason for the removal stays visible.
static void writeProteinStats(const std::vector<TargetStats> &stats,
                              DBReader<unsigned int> &targetHeaderReader,
                              const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("target_key\ttarget_id\tabundance_pct\treads_top_hit\treads_expected\tcoverage_depth\tbreadth_of_coverage\t"
          "coverage_evenness\tDrop(y/n)\tdrop_round\tcovered_length\ttarget_length\n", handle);

    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int key = stats[i].key;
        const std::string targetId = identifierForKey(targetHeaderReader, key, 0);

        const double coveredLength = stats[i].coveredLength;
        const double breadthOfCoverage = (stats[i].targetLength > 0)
                ? coveredLength / static_cast<double>(stats[i].targetLength)
                : 0.0;

        fprintf(handle, "%u\t%s\t%.12g\t%zu\t%.12g\t%.12g\t%.12g\t%.12g\t%s\t%d\t%.12g\t%u\n",
                key,
                targetId.c_str(),
                stats[i].abundance,
                stats[i].topHitReads,
                stats[i].expectedReads,
                stats[i].coverageDepth,
                breadthOfCoverage,
                stats[i].coverageEvenness,
                stats[i].dropped ? "y" : "n",
                stats[i].droppedRound,
                coveredLength,
                stats[i].targetLength);
    }

    fclose(handle);
}

// --ref-group 1: the per-target rows summed per genome (--target-groups, as in `mmseqs em`); a target outside every
// genome is a row of its own.  length = every target of the genome in targetDB (hit or not), the L_g of the EM.
// abundance_pct = share of the reads (as per target); coverage_pct = share of reads_expected / length, i.e. the
// genomes' relative coverage depth (cell abundance rather than read share).
static void writeGroupStats(const std::vector<TargetStats> &stats,
                            const std::string &targetDb,
                            const std::string &spec,
                            DBReader<unsigned int> &targetHeaderReader,
                            const std::string &path) {
    DBReader<unsigned int> tdb(targetDb.c_str(), (targetDb + ".index").c_str(), 1, DBReader<unsigned int>::USE_INDEX);
    tdb.open(DBReader<unsigned int>::NOSORT);
    unsigned int maxKey = 0;
    for (size_t i = 0; i < tdb.getSize(); ++i) {
        maxKey = std::max(maxKey, tdb.getDbKey(i));
    }
    EMTargetGroups::Groups groups = EMTargetGroups::load(targetDb, spec, maxKey);
    struct Row {
        std::string name;
        double abundance = 0.0, readsExpected = 0.0, length = 0.0;
        size_t readsTopHit = 0, targets = 0, targetsWithReads = 0;
    };
    std::vector<Row> rows(groups.names.size());
    for (size_t g = 0; g < rows.size(); ++g) {
        rows[g].name = groups.names[g];
    }
    for (size_t i = 0; i < tdb.getSize(); ++i) {
        const unsigned int g = groups.groupOfKey[tdb.getDbKey(i)];
        if (g != EMTargetGroups::NONE) {
            rows[g].length += static_cast<double>(tdb.getSeqLen(i));
            rows[g].targets += 1;
        }
    }
    tdb.close();
    Debug(Debug::INFO) << "Genomes from " << groups.how << ": " << groups.mapped << " targets in "
                       << groups.names.size() << " genomes.\n";
    for (size_t i = 0; i < stats.size(); ++i) {
        if (stats[i].dropped) {
            continue;
        }
        const unsigned int key = stats[i].key;
        unsigned int g = key <= maxKey ? groups.groupOfKey[key] : EMTargetGroups::NONE;
        if (g == EMTargetGroups::NONE) {
            g = static_cast<unsigned int>(rows.size());
            rows.push_back(Row());
            rows[g].name = identifierForKey(targetHeaderReader, key, 0);
            rows[g].length = stats[i].targetLength;
            rows[g].targets = 1;
        }
        rows[g].abundance += stats[i].abundance;
        rows[g].readsExpected += stats[i].expectedReads;
        rows[g].readsTopHit += stats[i].topHitReads;
        rows[g].targetsWithReads += 1;
    }
    double densitySum = 0.0;
    for (size_t g = 0; g < rows.size(); ++g) {
        densitySum += rows[g].length > 0.0 ? rows[g].readsExpected / rows[g].length : 0.0;
    }
    std::vector<size_t> order;
    for (size_t g = 0; g < rows.size(); ++g) {
        if (rows[g].targetsWithReads > 0) {
            order.push_back(g);
        }
    }
    std::sort(order.begin(), order.end(), [&rows](size_t a, size_t b) {
        if (rows[a].abundance != rows[b].abundance) {
            return rows[a].abundance > rows[b].abundance;
        }
        return rows[a].name < rows[b].name;
    });
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("genome\tabundance_pct\tcoverage_pct\treads_expected\treads_top_hit\tdensity\tlength\ttargets_with_reads\ttargets\n", handle);
    for (size_t k = 0; k < order.size(); ++k) {
        const Row &r = rows[order[k]];
        const double density = r.length > 0.0 ? r.readsExpected / r.length : 0.0;
        fprintf(handle, "%s\t%.12g\t%.12g\t%.12g\t%zu\t%.12g\t%.0f\t%zu\t%zu\n", r.name.c_str(), r.abundance,
                densitySum > 0.0 ? 100.0 * density / densitySum : 0.0, r.readsExpected, r.readsTopHit, density,
                r.length, r.targetsWithReads, r.targets);
    }
    fclose(handle);
    Debug(Debug::INFO) << order.size() << " genomes (or lone targets) with reads written to " << path << ".\n";
}

// The surviving reads with their surviving hits, column 3 = renormalised posterior, best first.
// Deleted reads and reads that had no hit are written as empty entries, so the output holds every
// key of the input DB and `mmseqs convertalis` / `createsubdb` work on it unchanged.
static void writeFilteredDb(ReclassTaxContext &ctx,
                            int dbType,
                            const std::string &outDb,
                            const std::string &outIndex,
                            int threads,
                            bool compress) {
    HitData &hits = ctx.hits;
    const size_t queryCount = hits.queryCount();
    const size_t hitCount = hits.hitCount();
    DBWriter writer(outDb.c_str(), outIndex.c_str(), threads, compress, dbType);
    writer.open();

#pragma omp parallel for schedule(static) num_threads(usableThreads(threads))
    for (size_t h = 0; h < hitCount; ++h) {
        hits.hitResult[h].seqId = static_cast<float>(hits.hitPosterior[h]);
    }

    Debug::Progress progress(queryCount);
#pragma omp parallel num_threads(usableThreads(threads))
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024 + 32768 * 4];
        std::vector<size_t> order;

#pragma omp for schedule(dynamic, 5)
        for (size_t q = 0; q < queryCount; ++q) {
            progress.updateProgress();
            const size_t begin = hits.queryOffset[q];
            const size_t end = hits.queryOffset[q + 1];
            order.clear();
            for (size_t h = begin; h < end; ++h) {
                order.push_back(h);
            }
            SORT_SERIAL(order.begin(), order.end(), [&hits](size_t lhs, size_t rhs) {
                if (hits.hitPosterior[lhs] != hits.hitPosterior[rhs]) {
                    return hits.hitPosterior[lhs] > hits.hitPosterior[rhs];
                }
                if (hits.hitResult[lhs].score != hits.hitResult[rhs].score) {
                    return hits.hitResult[lhs].score > hits.hitResult[rhs].score;
                }
                return Matcher::compareHits(hits.hitResult[lhs], hits.hitResult[rhs]);
            });
            writer.writeStart(thread_idx);
            for (size_t j = 0; j < order.size(); ++j) {
                // The backtrace was parsed with readCompressed=true.
                size_t len = Matcher::resultToBuffer(buffer, hits.hitResult[order[j]], ctx.hasBacktrace,
                                                     false, ctx.hasOrfPosition, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(hits.queryKey[q], thread_idx);
        }
    }
    for (size_t i = 0; i < ctx.deletedQueryKey.size(); ++i) {
        writer.writeStart(0);
        writer.writeEnd(ctx.deletedQueryKey[i], 0);
    }
    for (size_t i = 0; i < ctx.emptyQueryKey.size(); ++i) {
        writer.writeStart(0);
        writer.writeEnd(ctx.emptyQueryKey[i], 0);
    }
    writer.close();
    Debug(Debug::INFO) << "Filtered alignment DB: " << queryCount << " reads with hits, "
                       << ctx.deletedQueryKey.size() << " deleted reads and " << ctx.emptyQueryKey.size()
                       << " reads without a hit written as empty entries.\n";
}

// Expected reads assigned directly to each taxon, kept fractional.

static std::unordered_map<TaxID, double> directCountsFromStats(const std::vector<TargetStats> &stats,
                                                               MappingReader &mapping,
                                                               size_t assignedReads) {
    std::unordered_map<TaxID, double> directCounts;
    directCounts.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        const TaxID taxId = mapping.lookup(stats[i].key);
        if (taxId == 0) {
            continue;
        }
        directCounts[taxId] += stats[i].abundance * static_cast<double>(assignedReads) / 100.0;
    }
    return directCounts;
}


static std::unordered_map<TaxID, double> directCountsByLca(const HitData &hits,
                                                          MappingReader &mapping,
                                                          NcbiTaxonomy *taxonomy,
                                                          double threshold,
                                                          int threads) {
    const size_t queryCount = hits.queryCount();
    const size_t blocks = std::min<size_t>(SCALAR_REDUCTION_BLOCKS, std::max<size_t>(1, queryCount));
    std::vector<std::unordered_map<TaxID, double> > partial(blocks);
    std::vector<size_t> bounds(blocks + 1, queryCount);
    bounds[0] = 0;
    for (size_t b = 1; b < blocks; ++b) {
        bounds[b] = queryCount * b / blocks;
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(usableThreads(threads))
    for (size_t b = 0; b < blocks; ++b) {
        std::unordered_map<TaxID, double> &local = partial[b];
        std::vector<std::pair<double, TaxID> > candidates;
        std::vector<TaxID> chosen;
        for (size_t q = bounds[b]; q < bounds[b + 1]; ++q) {
            const size_t begin = hits.queryOffset[q];
            const size_t end = hits.queryOffset[q + 1];
            candidates.clear();
            double total = 0.0;
            for (size_t h = begin; h < end; ++h) {
                const unsigned int t = hits.hitTarget[h];
                const TaxID taxId = mapping.lookup(hits.targetKey[t]);
                if (taxId == 0) {
                    continue;
                }
                const double p = hits.hitPosterior[h];
                if (p <= 0.0) {
                    continue;
                }
                candidates.push_back(std::make_pair(p, taxId));
                total += p;
            }
            if (candidates.empty() || total <= 0.0) {
                continue;                           // falls through to the `U` row
            }
            // Descending posterior; ties broken by taxid so the output is reproducible.
            SORT_SERIAL(candidates.begin(), candidates.end(),
                        [](const std::pair<double, TaxID> &a, const std::pair<double, TaxID> &b) {
                            return (a.first != b.first) ? (a.first > b.first) : (a.second < b.second);
                        });
            chosen.clear();
            double cumulative = 0.0;
            for (size_t i = 0; i < candidates.size(); ++i) {
                chosen.push_back(candidates[i].second);
                cumulative += candidates[i].first / total;
                if (cumulative >= threshold) {
                    break;
                }
            }
            TaxID assigned = chosen[0];
            if (chosen.size() > 1) {
                const TaxonNode *node = taxonomy->LCA(chosen);
                if (node != NULL) {
                    assigned = node->taxId;
                }
            }
            local[assigned] += 1.0;
        }
    }

    std::unordered_map<TaxID, double> directCounts;
    for (size_t b = 0; b < blocks; ++b) {
        for (std::unordered_map<TaxID, double>::const_iterator it = partial[b].begin(); it != partial[b].end(); ++it) {
            directCounts[it->first] += it->second;
        }
    }
    return directCounts;
}

static std::unordered_map<TaxID, double> cladeCountsFromDirect(const std::unordered_map<TaxID, double> &directCounts,
                                                              NcbiTaxonomy *taxonomy) {
    std::unordered_map<TaxID, double> cladeCounts;
    for (std::unordered_map<TaxID, double>::const_iterator it = directCounts.begin(); it != directCounts.end(); ++it) {
        cladeCounts[it->first] += it->second;
        if (taxonomy->nodeExists(it->first) == false) {
            continue;
        }
        const TaxonNode *node = taxonomy->taxonNode(it->first);
        while (node->parentTaxId != node->taxId && taxonomy->nodeExists(node->parentTaxId)) {
            node = taxonomy->taxonNode(node->parentTaxId);
            cladeCounts[node->taxId] += it->second;
        }
    }
    return cladeCounts;
}

static char rankCodeOf(NcbiTaxonomy *taxonomy, const TaxonNode *node) {
    if (node == NULL) {
        return '-';
    }
    const char *rankStr = taxonomy->getString(node->rankIdx);
    if (rankStr == NULL) {
        return '-';
    }
    std::map<std::string, char>::const_iterator rankIt = NcbiShortRanks.find(std::string(rankStr));
    return (rankIt != NcbiShortRanks.end()) ? rankIt->second : '-';
}

static double lookupCount(const std::unordered_map<TaxID, double> &counts, TaxID taxId) {
    std::unordered_map<TaxID, double>::const_iterator it = counts.find(taxId);
    return (it != counts.end()) ? it->second : 0.0;
}

static void writeKrakenReport(const std::unordered_map<TaxID, double> &directCounts,
                             const std::unordered_map<TaxID, double> &cladeCounts,
                             const std::unordered_map<TaxID, std::vector<TaxID> > &parentToChildren,
                             NcbiTaxonomy *taxonomy,
                             size_t totalReads,
                             const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    const double denom = (totalReads > 0) ? static_cast<double>(totalReads) : 1.0;

    // Unclassified row, so the percentage column sums to 100.
    const double classified = lookupCount(cladeCounts, 1);
    const double unclassified = std::max(0.0, denom - classified);
    fprintf(handle, "%.4f\t%.2f\t%.2f\tU\t0\tunclassified\n",
            100.0 * unclassified / denom, unclassified, unclassified);

    std::vector<TaxID> stack;
    std::vector<int> depthStack;
    stack.push_back(1);
    depthStack.push_back(0);

    while (stack.empty() == false) {
        const TaxID taxId = stack.back();
        stack.pop_back();
        const int depth = depthStack.back();
        depthStack.pop_back();

        const double cladeCount = lookupCount(cladeCounts, taxId);
        if (cladeCount <= 0.0) {
            continue;
        }
        const double directCount = lookupCount(directCounts, taxId);
        const TaxonNode *node = taxonomy->taxonNode(taxId, false);
        const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";

        for (int i = 0; i < depth; ++i) {
            fputs("  ", handle);
        }
        fprintf(handle, "%.4f\t%.2f\t%.2f\t%c\t%u\t%s\n",
                100.0 * cladeCount / denom, cladeCount, directCount,
                rankCodeOf(taxonomy, node), static_cast<unsigned int>(taxId), name);

        std::unordered_map<TaxID, std::vector<TaxID> >::const_iterator childIt = parentToChildren.find(taxId);
        if (childIt == parentToChildren.end()) {
            continue;
        }
        std::vector<std::pair<double, TaxID> > children;
        children.reserve(childIt->second.size());
        for (size_t i = 0; i < childIt->second.size(); ++i) {
            const TaxID child = childIt->second[i];
            const double count = lookupCount(cladeCounts, child);
            if (count > 0.0) {
                children.push_back(std::make_pair(count, child));
            }
        }
        std::sort(children.begin(), children.end(), [](const std::pair<double, TaxID> &lhs,
                                                      const std::pair<double, TaxID> &rhs) {
            if (lhs.first != rhs.first) {
                return lhs.first < rhs.first;
            }
            return lhs.second > rhs.second;
        });
        for (size_t i = 0; i < children.size(); ++i) {
            stack.push_back(children[i].second);
            depthStack.push_back(depth + 1);
        }
    }

    fclose(handle);
}

static void writeBrackenReport(const std::unordered_map<TaxID, double> &directCounts,
                              const std::unordered_map<TaxID, double> &cladeCounts,
                              NcbiTaxonomy *taxonomy,
                              size_t totalReads,
                              const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n", handle);

    const double denom = (totalReads > 0) ? static_cast<double>(totalReads) : 1.0;

    // Sorted output: iterating the hash map directly made the row order differ between
    // runs on identical input, which makes results impossible to diff.
    std::vector<std::pair<double, TaxID> > rows;
    rows.reserve(cladeCounts.size());
    for (std::unordered_map<TaxID, double>::const_iterator it = cladeCounts.begin(); it != cladeCounts.end(); ++it) {
        if (it->second > 0.0) {
            rows.push_back(std::make_pair(it->second, it->first));
        }
    }
    std::sort(rows.begin(), rows.end(), [](const std::pair<double, TaxID> &lhs,
                                          const std::pair<double, TaxID> &rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first > rhs.first;
        }
        return lhs.second < rhs.second;
    });

    for (size_t i = 0; i < rows.size(); ++i) {
        const TaxID taxId = rows[i].second;
        const double cladeCount = rows[i].first;
        const double directCount = lookupCount(directCounts, taxId);
        const TaxonNode *node = taxonomy->taxonNode(taxId, false);
        const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";
        fprintf(handle, "%s\t%u\t%c\t%.2f\t%.2f\t%.2f\t%.12g\n",
                name,
                static_cast<unsigned int>(taxId),
                rankCodeOf(taxonomy, node),
                directCount,
                std::max(0.0, cladeCount - directCount),
                cladeCount,
                cladeCount / denom);
    }

    fclose(handle);
}
}

int emabundance(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const bool withTaxonomy = (par.reclassifyTaxonomy == 1);
    if (withTaxonomy) {
        // fail before loading the alignments: without these files openTaxonomy / MappingReader stop with a bare
        // "file not found" that does not say how to make them
        std::vector<std::string> missing;
        if (!FileUtil::fileExists((par.db2 + "_mapping").c_str())) {
            missing.push_back(par.db2 + "_mapping");
        }
        if (!FileUtil::fileExists((par.db2 + "_taxonomy").c_str())
            && !(FileUtil::fileExists((par.db2 + "_nodes.dmp").c_str())
                 && FileUtil::fileExists((par.db2 + "_names.dmp").c_str())
                 && FileUtil::fileExists((par.db2 + "_merged.dmp").c_str()))) {
            missing.push_back(par.db2 + "_taxonomy");
        }
        if (!missing.empty()) {
            Debug(Debug::ERROR) << "--taxonomy 1 needs the taxonomy of the target DB, but " << missing[0]
                                << (missing.size() > 1 ? " and " + missing[1] + " are" : " is") << " missing.\n"
                                << "You need to run 'mmseqs createtaxdb " << par.db2 << " tmp' before you run em-abundance "
                                << "(or use --taxonomy 0 for the per-target table, which needs no taxonomy).\n";
            EXIT(EXIT_FAILURE);
        }
    }
    NcbiTaxonomy *taxonomy = NULL;
    MappingReader *mapping = NULL;
    if (withTaxonomy) {
        taxonomy = NcbiTaxonomy::openTaxonomy(par.db2);
        mapping = new MappingReader(par.db2);
    }

    ReclassTaxContext ctx;
    loadAlignmentDb(reader, ctx, par.threads);
    Debug(Debug::INFO) << "Loaded " << ctx.hits.queryCount() << " queries with hits and "
                       << ctx.hits.targetCount() << " unique targets.\n";
    reportInputDiagnostics(ctx);

    const int minReads = std::max(0, par.abundanceMinReads);
    if (minReads > 0) {
        Debug(Debug::INFO) << "Removing targets that are the top hit of < " << minReads
                           << " reads, together with those reads...\n";
    }
    filterTargets(ctx, minReads, par.threads);
    const size_t queryCount = ctx.hits.queryCount();

    Debug(Debug::INFO) << "Computing abundance from posterior...\n";
    const TargetHitIndex byTarget = groupHitsByTarget(ctx.hits);
    computeAbundanceFromPosterior(ctx.hits, byTarget, par.threads);

    Debug(Debug::INFO) << "Collecting per-target statistics...\n";
    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx, byTarget, par.threads);
    convertAbundanceToPercent(allTargetStats);

    Debug(Debug::INFO) << "Writing the filtered alignment DB...\n";
    writeFilteredDb(ctx, reader.getDbtype(), par.db4, par.db4Index, par.threads, par.compressed);

    Debug(Debug::INFO) << "Writing output...\n";
    if (withTaxonomy) {
        // The percentage basis is every read in the query DB, not just the reads that got a
        // hit (queryCount). Reads with no hit at all (and the reads deleted by the filter) have
        // to be counted somewhere for the report to be comparable with Kraken2/Bracken, and
        // they land in the `U` row.
        DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                           DBReader<unsigned int>::USE_INDEX);
        queryReader.open(DBReader<unsigned int>::NOSORT);
        const size_t totalReads = queryReader.getSize();
        queryReader.close();
        Debug(Debug::INFO) << totalReads << " reads in the query DB, " << queryCount
                           << " with at least one hit after the filter.\n";

        std::vector<TargetStats> targetStats = allTargetStats;
        targetStats.erase(std::remove_if(targetStats.begin(), targetStats.end(), [](const TargetStats &entry) {
            return entry.dropped;
        }), targetStats.end());

        const double lcaThreshold = std::max(0.0, std::min(1.0, par.abundanceLcaThreshold));
        if (lcaThreshold != par.abundanceLcaThreshold) {
            Debug(Debug::WARNING) << "--lca-threshold " << par.abundanceLcaThreshold
                                  << " is outside [0, 1] and was clamped to " << lcaThreshold
                                  << ". It is a cumulative posterior, not a percentage.\n";
        }
        std::unordered_map<TaxID, double> directCounts;
        if (lcaThreshold > 0.0) {
            Debug(Debug::INFO) << "LCA fallback: --lca-threshold " << lcaThreshold
                               << ", assigning each read to the LCA of its top hits.\n";
            directCounts = directCountsByLca(ctx.hits, *mapping, taxonomy, lcaThreshold, par.threads);
        } else {
            directCounts = directCountsFromStats(targetStats, *mapping, queryCount);
        }
        const std::unordered_map<TaxID, double> cladeCounts =
                cladeCountsFromDirect(directCounts, taxonomy);
        const std::unordered_map<TaxID, std::vector<TaxID> > parentToChildren = taxonomy->getParentToChildren();

        writeKrakenReport(directCounts, cladeCounts, parentToChildren, taxonomy, totalReads, par.db5);
        writeBrackenReport(directCounts, cladeCounts, taxonomy, totalReads, par.db5 + ".bracken");
    } else {
        DBReader<unsigned int> targetHeaderReader((par.db2 + "_h").c_str(), (par.db2 + "_h.index").c_str(), par.threads,
                                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        targetHeaderReader.open(DBReader<unsigned int>::NOSORT);
        writeProteinStats(allTargetStats, targetHeaderReader, par.db5);
        targetHeaderReader.close();
    }
    if (par.emRefGroup) {
        DBReader<unsigned int> targetHeaderReader((par.db2 + "_h").c_str(), (par.db2 + "_h.index").c_str(), par.threads,
                                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        targetHeaderReader.open(DBReader<unsigned int>::NOSORT);
        writeGroupStats(allTargetStats, par.db2, par.emTargetGroups, targetHeaderReader, par.db5 + ".groups.tsv");
        targetHeaderReader.close();
    }

    delete mapping;
    delete taxonomy;
    reader.close();
    return EXIT_SUCCESS;
}
